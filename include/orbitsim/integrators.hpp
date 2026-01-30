#pragma once

#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

namespace orbitsim
{

    inline void advance_spin(SpinState &spin, const double dt_s)
    {
        spin.axis = normalized_or(spin.axis, Vec3{0.0, 1.0, 0.0});
        if (std::isfinite(spin.rate_rad_per_s) && std::isfinite(dt_s) && std::isfinite(spin.angle_rad))
        {
            spin.angle_rad += spin.rate_rad_per_s * dt_s;
        }
    }

    inline void compute_nbody_accelerations(const std::vector<MassiveBody> &bodies, const double G,
                                            const double softening_length_m, std::vector<Vec3> &out_acc_mps2)
    {
        const std::size_t n = bodies.size();
        out_acc_mps2.assign(n, Vec3{0.0, 0.0, 0.0});
        const double eps2 = softening_length_m * softening_length_m;

        for (std::size_t i = 0; i < n; ++i)
        {
            for (std::size_t j = i + 1; j < n; ++j)
            {
                const Vec3 dr = bodies[j].state.position_m - bodies[i].state.position_m;
                const double r2 = glm::dot(dr, dr) + eps2;
                if (!(r2 > 0.0) || !std::isfinite(r2))
                {
                    continue;
                }
                const double inv_r = 1.0 / std::sqrt(r2);
                const double inv_r3 = inv_r * inv_r * inv_r;
                const Vec3 a_dir = (G * inv_r3) * dr;

                out_acc_mps2[i] += a_dir * bodies[j].mass_kg;
                out_acc_mps2[j] -= a_dir * bodies[i].mass_kg;
            }
        }
    }

    namespace detail
    {
        inline void velocity_verlet_step(std::vector<MassiveBody> &bodies, const double dt_s, const double G,
                                         const double softening_length_m)
        {
            std::vector<Vec3> acc;
            compute_nbody_accelerations(bodies, G, softening_length_m, acc);

            for (std::size_t i = 0; i < bodies.size(); ++i)
            {
                bodies[i].state.velocity_mps += 0.5 * dt_s * acc[i];
            }
            for (auto &b: bodies)
            {
                b.state.position_m += dt_s * b.state.velocity_mps;
            }

            compute_nbody_accelerations(bodies, G, softening_length_m, acc);
            for (std::size_t i = 0; i < bodies.size(); ++i)
            {
                bodies[i].state.velocity_mps += 0.5 * dt_s * acc[i];
            }

            for (auto &b: bodies)
            {
                advance_spin(b.state.spin, dt_s);
            }
        }
    } // namespace detail

    inline void symplectic4_step(std::vector<MassiveBody> &bodies, const double dt_s, const double G,
                                 const double softening_length_m)
    {
        // Yoshida 4th order composition of velocity-Verlet:
        // step(dt) = VV(w1*dt) ∘ VV(w0*dt) ∘ VV(w1*dt)
        const double cbrt2 = std::cbrt(2.0);
        const double w1 = 1.0 / (2.0 - cbrt2);
        const double w0 = -cbrt2 * w1;

        detail::velocity_verlet_step(bodies, w1 * dt_s, G, softening_length_m);
        detail::velocity_verlet_step(bodies, w0 * dt_s, G, softening_length_m);
        detail::velocity_verlet_step(bodies, w1 * dt_s, G, softening_length_m);
    }

    struct DOPRI5Options
    {
        bool adaptive{true};
        double abs_tol{1e-3}; // meters or (m/s) depending on component (scaled internally)
        double rel_tol{1e-9};
        int max_substeps{64};
        // If > 0, clamps the internal step size |h| (helps prevent "big jumps" over long coast intervals).
        // Note: if max_step_s is too small relative to (dt_s / max_substeps), the integrator will fall back
        // to a final single step to finish the interval.
        double max_step_s{0.0};
        double min_step_s{1e-6};
    };

    struct SpacecraftKinematics
    {
        Vec3 position_m{0.0, 0.0, 0.0};
        Vec3 velocity_mps{0.0, 0.0, 0.0};
    };

    inline double dopri5_error_norm(const SpacecraftKinematics &y, const SpacecraftKinematics &y_next,
                                    const SpacecraftKinematics &err, const double abs_tol, const double rel_tol)
    {
        auto scale = [&](double a, double b) { return abs_tol + rel_tol * std::max(std::abs(a), std::abs(b)); };

        const std::array<double, 6> yv{y.position_m.x,   y.position_m.y,   y.position_m.z,
                                       y.velocity_mps.x, y.velocity_mps.y, y.velocity_mps.z};
        const std::array<double, 6> yn{y_next.position_m.x,   y_next.position_m.y,   y_next.position_m.z,
                                       y_next.velocity_mps.x, y_next.velocity_mps.y, y_next.velocity_mps.z};
        const std::array<double, 6> ev{err.position_m.x,   err.position_m.y,   err.position_m.z,
                                       err.velocity_mps.x, err.velocity_mps.y, err.velocity_mps.z};

        double sum = 0.0;
        for (std::size_t i = 0; i < 6; ++i)
        {
            const double s = scale(yv[i], yn[i]);
            if (!(s > 0.0) || !std::isfinite(s))
            {
                continue;
            }
            const double r = ev[i] / s;
            sum += r * r;
        }
        return std::sqrt(sum / 6.0);
    }

    inline SpacecraftKinematics dopri5_single_step(
            const double t_s, const SpacecraftKinematics &y, const double h_s,
            const std::function<Vec3(double /*t_s*/, const Vec3 & /*pos_m*/, const Vec3 & /*vel_mps*/)> &accel_mps2,
            SpacecraftKinematics *out_err)
    {
        // Dormand–Prince 5(4) tableau (ode45).
        auto f = [&](double t, const SpacecraftKinematics &s) {
            SpacecraftKinematics dy;
            dy.position_m = s.velocity_mps;
            dy.velocity_mps = accel_mps2(t, s.position_m, s.velocity_mps);
            return dy;
        };

        const SpacecraftKinematics k1 = f(t_s, y);
        const SpacecraftKinematics k2 =
                f(t_s + h_s * (1.0 / 5.0), {y.position_m + h_s * (1.0 / 5.0) * k1.position_m,
                                            y.velocity_mps + h_s * (1.0 / 5.0) * k1.velocity_mps});
        const SpacecraftKinematics k3 =
                f(t_s + h_s * (3.0 / 10.0),
                  {y.position_m + h_s * ((3.0 / 40.0) * k1.position_m + (9.0 / 40.0) * k2.position_m),
                   y.velocity_mps + h_s * ((3.0 / 40.0) * k1.velocity_mps + (9.0 / 40.0) * k2.velocity_mps)});
        const SpacecraftKinematics k4 =
                f(t_s + h_s * (4.0 / 5.0),
                  {y.position_m + h_s * ((44.0 / 45.0) * k1.position_m + (-56.0 / 15.0) * k2.position_m +
                                         (32.0 / 9.0) * k3.position_m),
                   y.velocity_mps + h_s * ((44.0 / 45.0) * k1.velocity_mps + (-56.0 / 15.0) * k2.velocity_mps +
                                           (32.0 / 9.0) * k3.velocity_mps)});
        const SpacecraftKinematics k5 = f(
                t_s + h_s * (8.0 / 9.0),
                {y.position_m + h_s * ((19372.0 / 6561.0) * k1.position_m + (-25360.0 / 2187.0) * k2.position_m +
                                       (64448.0 / 6561.0) * k3.position_m + (-212.0 / 729.0) * k4.position_m),
                 y.velocity_mps + h_s * ((19372.0 / 6561.0) * k1.velocity_mps + (-25360.0 / 2187.0) * k2.velocity_mps +
                                         (64448.0 / 6561.0) * k3.velocity_mps + (-212.0 / 729.0) * k4.velocity_mps)});
        const SpacecraftKinematics k6 =
                f(t_s + h_s,
                  {y.position_m + h_s * ((9017.0 / 3168.0) * k1.position_m + (-355.0 / 33.0) * k2.position_m +
                                         (46732.0 / 5247.0) * k3.position_m + (49.0 / 176.0) * k4.position_m +
                                         (-5103.0 / 18656.0) * k5.position_m),
                   y.velocity_mps + h_s * ((9017.0 / 3168.0) * k1.velocity_mps + (-355.0 / 33.0) * k2.velocity_mps +
                                           (46732.0 / 5247.0) * k3.velocity_mps + (49.0 / 176.0) * k4.velocity_mps +
                                           (-5103.0 / 18656.0) * k5.velocity_mps)});
        const SpacecraftKinematics k7 =
                f(t_s + h_s,
                  {y.position_m + h_s * ((35.0 / 384.0) * k1.position_m + (500.0 / 1113.0) * k3.position_m +
                                         (125.0 / 192.0) * k4.position_m + (-2187.0 / 6784.0) * k5.position_m +
                                         (11.0 / 84.0) * k6.position_m),
                   y.velocity_mps + h_s * ((35.0 / 384.0) * k1.velocity_mps + (500.0 / 1113.0) * k3.velocity_mps +
                                           (125.0 / 192.0) * k4.velocity_mps + (-2187.0 / 6784.0) * k5.velocity_mps +
                                           (11.0 / 84.0) * k6.velocity_mps)});

        // 5th order solution.
        const SpacecraftKinematics y5{
                y.position_m + h_s * ((35.0 / 384.0) * k1.position_m + (500.0 / 1113.0) * k3.position_m +
                                      (125.0 / 192.0) * k4.position_m + (-2187.0 / 6784.0) * k5.position_m +
                                      (11.0 / 84.0) * k6.position_m),
                y.velocity_mps + h_s * ((35.0 / 384.0) * k1.velocity_mps + (500.0 / 1113.0) * k3.velocity_mps +
                                        (125.0 / 192.0) * k4.velocity_mps + (-2187.0 / 6784.0) * k5.velocity_mps +
                                        (11.0 / 84.0) * k6.velocity_mps),
        };

        // 4th order solution (embedded).
        const SpacecraftKinematics y4{
                y.position_m + h_s * ((5179.0 / 57600.0) * k1.position_m + (7571.0 / 16695.0) * k3.position_m +
                                      (393.0 / 640.0) * k4.position_m + (-92097.0 / 339200.0) * k5.position_m +
                                      (187.0 / 2100.0) * k6.position_m + (1.0 / 40.0) * k7.position_m),
                y.velocity_mps + h_s * ((5179.0 / 57600.0) * k1.velocity_mps + (7571.0 / 16695.0) * k3.velocity_mps +
                                        (393.0 / 640.0) * k4.velocity_mps + (-92097.0 / 339200.0) * k5.velocity_mps +
                                        (187.0 / 2100.0) * k6.velocity_mps + (1.0 / 40.0) * k7.velocity_mps),
        };

        if (out_err != nullptr)
        {
            out_err->position_m = y5.position_m - y4.position_m;
            out_err->velocity_mps = y5.velocity_mps - y4.velocity_mps;
        }
        return y5;
    }

    inline SpacecraftKinematics dopri5_integrate_interval(
            const double t0_s, const SpacecraftKinematics &y0, const double dt_s,
            const std::function<Vec3(double /*t_s*/, const Vec3 & /*pos_m*/, const Vec3 & /*vel_mps*/)> &accel_mps2,
            const DOPRI5Options &opt)
    {
        if (!(dt_s != 0.0) || !std::isfinite(dt_s))
        {
            return y0;
        }

        const double dir = (dt_s > 0.0) ? 1.0 : -1.0;
        const double t_end = t0_s + dt_s;
        double t = t0_s;
        SpacecraftKinematics y = y0;

        const double eff_max_step_s =
                (opt.max_step_s > 0.0) ? std::max(opt.max_step_s, opt.min_step_s) : 0.0;
        const auto clamp_h = [&](double h_s) -> double {
            if (!(eff_max_step_s > 0.0))
            {
                return h_s;
            }
            if (std::abs(h_s) > eff_max_step_s)
            {
                return eff_max_step_s * ((h_s >= 0.0) ? 1.0 : -1.0);
            }
            return h_s;
        };

        double h = dt_s;
        if (opt.adaptive)
        {
            const double abs_dt = std::abs(dt_s);
            if (abs_dt <= opt.min_step_s)
            {
                return dopri5_single_step(t0_s, y0, dt_s, accel_mps2, nullptr);
            }
            h = abs_dt * dir;
        }
        h = clamp_h(h);

        const double safety = 0.9;
        const double min_fac = 0.2;
        const double max_fac = 5.0;

        int substeps = 0;
        while ((dir > 0.0) ? (t < t_end) : (t > t_end))
        {
            if (++substeps > std::max(1, opt.max_substeps))
            {
                return dopri5_single_step(t, y, t_end - t, accel_mps2, nullptr);
            }
            const double remaining = t_end - t;
            if (std::abs(h) > std::abs(remaining))
            {
                h = remaining;
            }
            h = clamp_h(h);

            SpacecraftKinematics err{};
            const SpacecraftKinematics y_next = dopri5_single_step(t, y, h, accel_mps2, &err);

            if (!opt.adaptive)
            {
                y = y_next;
                t += h;
                continue;
            }

            const double err_norm = dopri5_error_norm(y, y_next, err, opt.abs_tol, opt.rel_tol);
            if (!(err_norm >= 0.0) || !std::isfinite(err_norm))
            {
                y = y_next;
                t += h;
                continue;
            }

            if (err_norm <= 1.0)
            {
                y = y_next;
                t += h;
                const double fac =
                        (err_norm == 0.0) ? max_fac : std::clamp(safety * std::pow(err_norm, -0.2), min_fac, max_fac);
                h *= fac;
                h = clamp_h(h);
            }
            else
            {
                const double fac = std::clamp(safety * std::pow(err_norm, -0.2), min_fac, 1.0);
                h *= fac;
                h = clamp_h(h);
                if (std::abs(h) < opt.min_step_s)
                {
                    h = opt.min_step_s * dir;
                }
            }
        }

        return y;
    }

} // namespace orbitsim
