#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/events.hpp"
#include "orbitsim/integrators.hpp"
#include "orbitsim/maneuvers.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace orbitsim
{

    class GameSimulation
    {
    public:
        struct Config
        {
            double gravitational_constant{orbitsim::kGravitationalConstant_SI};
            double softening_length_m{0.0};
            DOPRI5Options spacecraft_integrator{};
            EventOptions events{};
            bool enable_events{true};
        };

        GameSimulation() = default;
        explicit GameSimulation(Config cfg) : cfg_(std::move(cfg)) {}

        double time_s() const { return time_s_; }

        std::vector<MassiveBody> &massive_bodies() { return massive_; }
        const std::vector<MassiveBody> &massive_bodies() const { return massive_; }

        std::vector<Spacecraft> &spacecraft() { return spacecraft_; }
        const std::vector<Spacecraft> &spacecraft() const { return spacecraft_; }

        ManeuverPlan &maneuver_plan() { return plan_; }
        const ManeuverPlan &maneuver_plan() const { return plan_; }

        std::size_t select_primary_by_max_accel(const std::size_t spacecraft_index) const
        {
            if (spacecraft_index >= spacecraft_.size() || massive_.empty())
            {
                return 0;
            }
            const Vec3 p = spacecraft_[spacecraft_index].state.position_m;
            const double eps2 = cfg_.softening_length_m * cfg_.softening_length_m;

            std::size_t best = 0;
            double best_a = -1.0;
            for (std::size_t i = 0; i < massive_.size(); ++i)
            {
                const Vec3 dr = massive_[i].state.position_m - p;
                const double r2 = glm::dot(dr, dr) + eps2;
                if (!(r2 > 0.0) || !std::isfinite(r2))
                {
                    continue;
                }
                const double amag = (cfg_.gravitational_constant * massive_[i].mass_kg) / r2;
                if (amag > best_a)
                {
                    best_a = amag;
                    best = i;
                }
            }
            return best;
        }

        void step(const double dt_s)
        {
            if (!(dt_s != 0.0) || !std::isfinite(dt_s))
            {
                return;
            }
            if (!cfg_.enable_events || cfg_.events.max_event_splits_per_step <= 0)
            {
                do_step_no_events_(dt_s);
                return;
            }
            if (dt_s < 0.0)
            {
                // Backwards integration with event splitting is not supported.
                do_step_no_events_(dt_s);
                return;
            }

            sort_segments_by_start(plan_);

            double remaining = dt_s;
            int splits = 0;
            while (remaining > 0.0 && splits++ < std::max(1, cfg_.events.max_event_splits_per_step))
            {
                // Preview a full remaining step to build a single ephemeris segment for event search.
                std::vector<MassiveBody> massive_preview = massive_;
                std::vector<Spacecraft> spacecraft_preview = spacecraft_;
                double t_preview = time_s_;

                CelestialEphemerisSegment eph_preview{};
                preview_step_no_events_(massive_preview, spacecraft_preview, t_preview, remaining, &eph_preview);

                std::optional<Event> best;
                for (const auto &sc: spacecraft_)
                {
                    auto propagate_sc = [&](const Spacecraft &sc_start, const double t0_s,
                                            const double dt_sc_s) -> Spacecraft {
                        return propagate_spacecraft_(sc_start, eph_preview, t0_s, dt_sc_s);
                    };
                    const std::optional<Event> e = find_earliest_event_in_interval(
                            massive_, eph_preview, sc, time_s_, remaining, plan_, cfg_.events, propagate_sc);
                    if (e.has_value() && (!best.has_value() || e->t_event_s < best->t_event_s))
                    {
                        best = e;
                    }
                }

                if (!best.has_value())
                {
                    do_step_no_events_(remaining);
                    return;
                }

                double dt_event = best->t_event_s - time_s_;
                const double min_step = std::max(0.0, cfg_.events.time_tol_s);
                if (!(dt_event > min_step) || !std::isfinite(dt_event))
                {
                    dt_event = min_step;
                }
                if (dt_event > remaining)
                {
                    dt_event = remaining;
                }

                do_step_no_events_(dt_event);
                remaining -= dt_event;
            }

            if (remaining > 0.0)
            {
                do_step_no_events_(remaining);
            }
        }

    private:
        static inline void snapshot_states_(const std::vector<MassiveBody> &massive, std::vector<State> *out)
        {
            if (out == nullptr)
            {
                return;
            }
            out->clear();
            out->reserve(massive.size());
            for (const auto &b: massive)
            {
                out->push_back(b.state);
            }
        }

        inline CelestialEphemerisSegment make_segment_(const std::vector<State> &start, const std::vector<State> &end,
                                                       const double t0_s, const double dt_s) const
        {
            CelestialEphemerisSegment eph;
            eph.t0_s = t0_s;
            eph.dt_s = dt_s;
            eph.start = start;
            eph.end = end;
            return eph;
        }

        inline Vec3 gravity_accel_mps2_(const CelestialEphemerisSegment &eph, const double t_s, const Vec3 &pos_m) const
        {
            Vec3 a{0.0, 0.0, 0.0};
            const double eps2 = cfg_.softening_length_m * cfg_.softening_length_m;
            for (std::size_t i = 0; i < massive_.size(); ++i)
            {
                const Vec3 p = eph.body_position_at(i, t_s);
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
        }

        inline Spacecraft propagate_spacecraft_(const Spacecraft &sc0, const CelestialEphemerisSegment &eph,
                                                const double t0_s, const double dt_s) const
        {
            Spacecraft out = sc0;
            if (!(dt_s != 0.0) || !std::isfinite(dt_s))
            {
                return out;
            }
            if (dt_s < 0.0)
            {
                // Backwards integration currently supports gravity-only (no maneuvers, no prop changes).
                const SpacecraftKinematics y0{out.state.position_m, out.state.velocity_mps};
                auto accel = [&](double t_eval_s, const Vec3 &pos_m, const Vec3 & /*vel_mps*/) -> Vec3 {
                    return gravity_accel_mps2_(eph, t_eval_s, pos_m);
                };
                DOPRI5Stats stats{};
                const SpacecraftKinematics y1 =
                        dopri5_integrate_interval(t0_s, y0, dt_s, accel, cfg_.spacecraft_integrator, &stats);
                (void) stats;

                out.state.position_m = y1.position_m;
                out.state.velocity_mps = y1.velocity_mps;
                advance_spin(out.state.spin, dt_s);
                return out;
            }

            const double t_end_s = t0_s + dt_s;
            double t = t0_s;
            double remaining = dt_s;

            SpacecraftKinematics y{out.state.position_m, out.state.velocity_mps};

            auto primary_for = [&](const BurnSegment &seg, const std::size_t sc_index_guess) -> std::size_t {
                if (seg.primary_index < massive_.size())
                {
                    return seg.primary_index;
                }
                return select_primary_by_max_accel(sc_index_guess);
            };

            // Best-effort index for auto-primary selection (only used when segment primary is invalid).
            std::size_t sc_index_guess = 0;

            while (remaining > 0.0)
            {
                const double boundary = next_burn_boundary_after(plan_, t, t_end_s);
                double h = std::min(remaining, boundary - t);
                if (!(h > 0.0) || !std::isfinite(h))
                {
                    break;
                }

                const BurnSegment *seg = active_burn_at(plan_, t);
                const bool has_engine = (seg != nullptr) && (seg->engine_index < out.engines.size());
                const Engine engine = has_engine ? out.engines[seg->engine_index] : Engine{};

                const double throttle = (seg != nullptr) ? clamp01(seg->throttle_0_1) : 0.0;
                const double min_thr = std::max(0.0, engine.min_throttle_0_1);
                const bool thrust_on = has_engine && (engine.max_thrust_N > 0.0) && (engine.isp_s > 0.0) &&
                                       (throttle >= min_thr) && (out.prop_mass_kg > 0.0);

                const double thrust_N = thrust_on ? (engine.max_thrust_N * throttle) : 0.0;
                const double mdot_kgps = thrust_on ? (thrust_N / (engine.isp_s * kStandardGravity_mps2)) : 0.0;

                if (mdot_kgps > 0.0 && out.prop_mass_kg > 0.0)
                {
                    const double h_burn = out.prop_mass_kg / mdot_kgps;
                    if (h_burn > 0.0 && h_burn < h)
                    {
                        // Burn until prop exhausted; the loop will then coast for the remaining time.
                        h = h_burn;
                    }
                }

                const double t_segment_start = t;
                const double prop_start_kg = out.prop_mass_kg;
                const std::size_t primary = (seg != nullptr) ? primary_for(*seg, sc_index_guess) : 0;
                const Vec3 dir_rtn = (seg != nullptr) ? seg->dir_rtn_unit : Vec3{0.0, 0.0, 0.0};

                auto accel = [&](double t_eval_s, const Vec3 &pos_m, const Vec3 &vel_mps) -> Vec3 {
                    Vec3 a = gravity_accel_mps2_(eph, t_eval_s, pos_m);
                    if (!(thrust_N > 0.0))
                    {
                        return a;
                    }

                    const double dt_prop = t_eval_s - t_segment_start;
                    double prop_kg = prop_start_kg - mdot_kgps * dt_prop;
                    if (prop_kg <= 0.0)
                    {
                        return a;
                    }

                    const double mass_kg = out.dry_mass_kg + prop_kg;
                    if (!(mass_kg > 0.0) || !std::isfinite(mass_kg))
                    {
                        return a;
                    }

                    const Vec3 dir_i = burn_dir_inertial_unit(eph, primary, t_eval_s, pos_m, vel_mps, dir_rtn);
                    const double dir2 = glm::dot(dir_i, dir_i);
                    if (!(dir2 > 0.0) || !std::isfinite(dir2))
                    {
                        return a;
                    }

                    a += (thrust_N / mass_kg) * dir_i;
                    return a;
                };

                DOPRI5Stats stats{};
                y = dopri5_integrate_interval(t, y, h, accel, cfg_.spacecraft_integrator, &stats);
                (void) stats;

                if (mdot_kgps > 0.0)
                {
                    out.prop_mass_kg = std::max(0.0, out.prop_mass_kg - mdot_kgps * h);
                }

                out.state.position_m = y.position_m;
                out.state.velocity_mps = y.velocity_mps;
                advance_spin(out.state.spin, h);

                t += h;
                remaining -= h;
            }

            return out;
        }

        inline void preview_step_no_events_(std::vector<MassiveBody> &massive, std::vector<Spacecraft> &spacecraft,
                                            double &t_s, const double dt_s, CelestialEphemerisSegment *out_eph) const
        {
            if (out_eph != nullptr)
            {
                *out_eph = {};
            }
            if (!(dt_s != 0.0) || !std::isfinite(dt_s))
            {
                return;
            }

            std::vector<State> start_states;
            snapshot_states_(massive, &start_states);

            symplectic4_step(massive, dt_s, cfg_.gravitational_constant, cfg_.softening_length_m);

            std::vector<State> end_states;
            snapshot_states_(massive, &end_states);

            CelestialEphemerisSegment eph = make_segment_(start_states, end_states, t_s, dt_s);
            if (out_eph != nullptr)
            {
                *out_eph = eph;
            }

            // Note: spacecraft preview ignores the maneuver plan mass/prop changes in the output container; event
            // search re-propagates from the true start state via propagate_spacecraft_ anyway.
            (void) spacecraft;
            t_s += dt_s;
        }

        inline void do_step_no_events_(const double dt_s)
        {
            if (!(dt_s != 0.0) || !std::isfinite(dt_s))
            {
                return;
            }

            std::vector<State> start_states;
            snapshot_states_(massive_, &start_states);

            symplectic4_step(massive_, dt_s, cfg_.gravitational_constant, cfg_.softening_length_m);

            std::vector<State> end_states;
            snapshot_states_(massive_, &end_states);

            const CelestialEphemerisSegment eph = make_segment_(start_states, end_states, time_s_, dt_s);

            for (auto &sc: spacecraft_)
            {
                sc = propagate_spacecraft_(sc, eph, time_s_, dt_s);
            }

            time_s_ += dt_s;
        }

        Config cfg_{};
        double time_s_{0.0};
        std::vector<MassiveBody> massive_{};
        std::vector<Spacecraft> spacecraft_{};
        ManeuverPlan plan_{};
    };

} // namespace orbitsim
