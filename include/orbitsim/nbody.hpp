#pragma once

#include "orbitsim/integrators.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace orbitsim
{

    class NBodySimulation
    {
    public:
        struct Config
        {
            double gravitational_constant{orbitsim::kGravitationalConstant_SI};
            double softening_length_m{0.0};
            DOPRI5Options spacecraft_integrator{};
        };

        NBodySimulation() = default;
        explicit NBodySimulation(Config cfg) : cfg_(std::move(cfg)) {}

        double time_s() const { return time_s_; }

        std::vector<MassiveBody> &massive_bodies() { return massive_; }
        const std::vector<MassiveBody> &massive_bodies() const { return massive_; }

        std::vector<Spacecraft> &spacecraft() { return spacecraft_; }
        const std::vector<Spacecraft> &spacecraft() const { return spacecraft_; }

        void step(const double dt_s)
        {
            if (!(dt_s != 0.0) || !std::isfinite(dt_s))
            {
                return;
            }

            // Snapshot massive body endpoints for interpolation.
            std::vector<State> start_states;
            start_states.reserve(massive_.size());
            for (const auto &b: massive_)
            {
                start_states.push_back(b.state);
            }

            // Advance massive bodies using symplectic 4th order integrator.
            symplectic4_step(massive_, dt_s, cfg_.gravitational_constant, cfg_.softening_length_m);

            std::vector<State> end_states;
            end_states.reserve(massive_.size());
            for (const auto &b: massive_)
            {
                end_states.push_back(b.state);
            }

            // Spacecraft: integrate in the time-varying gravity field using DOPRI5 with interpolation.
            for (auto &sc: spacecraft_)
            {
                const SpacecraftKinematics y0{sc.state.position_m, sc.state.velocity_mps};

                auto accel = [&](double t_eval_s, const Vec3 &pos_m, const Vec3 & /*vel_mps*/) -> Vec3 {
                    const double tau = (t_eval_s - time_s_) / dt_s;
                    Vec3 a{0.0, 0.0, 0.0};
                    const double eps2 = cfg_.softening_length_m * cfg_.softening_length_m;
                    for (std::size_t i = 0; i < massive_.size(); ++i)
                    {
                        const Vec3 p =
                                hermite_position(start_states[i].position_m, start_states[i].velocity_mps,
                                                 end_states[i].position_m, end_states[i].velocity_mps, dt_s, tau);
                        const Vec3 dr = p - pos_m;
                        const double r2 = glm::dot(dr, dr) + eps2;
                        if (!(r2 > 0.0) || !std::isfinite(r2))
                        {
                            continue;
                        }
                        const double inv_r = 1.0 / std::sqrt(r2);
                        const double inv_r3 = inv_r * inv_r * inv_r;
                        a += (cfg_.gravitational_constant * massive_[i].mass_kg) * inv_r3 * dr;
                    }
                    return a;
                };

                DOPRI5Stats stats{};
                const SpacecraftKinematics y1 =
                        dopri5_integrate_interval(time_s_, y0, dt_s, accel, cfg_.spacecraft_integrator, &stats);
                (void) stats;

                sc.state.position_m = y1.position_m;
                sc.state.velocity_mps = y1.velocity_mps;
                advance_spin(sc.state.spin, dt_s);
            }

            time_s_ += dt_s;
        }

    private:
        Config cfg_{};
        double time_s_{0.0};
        std::vector<MassiveBody> massive_{};
        std::vector<Spacecraft> spacecraft_{};
    };

} // namespace orbitsim
