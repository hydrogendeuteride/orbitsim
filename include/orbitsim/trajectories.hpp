#pragma once

#include "orbitsim/game_sim.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

namespace orbitsim
{

    struct TrajectorySample
    {
        double t_s{0.0};
        Vec3 position_m{0.0, 0.0, 0.0};
        Vec3 velocity_mps{0.0, 0.0, 0.0};
    };

    struct TrajectoryOptions
    {
        // Trajectory begins at the simulation's current time t0 = sim.time_s().
        double duration_s{3600.0};

        // If <= 0, it is derived from duration_s / (max_samples - 1).
        double sample_dt_s{10.0};

        // If > 0, uses this dt to step massive bodies when building a prediction ephemeris. Otherwise uses sample_dt_s
        // (or the derived sample dt).
        double celestial_dt_s{0.0};

        // Hard cap to protect against unbounded allocations.
        std::size_t max_samples{2048};

        bool include_start{true};
        bool include_end{true};

        // If set, returned state is relative to this body's state at the same time.
        std::optional<std::size_t> origin_body_index{};
    };

    namespace detail
    {
        inline double compute_sample_dt_(const TrajectoryOptions &opt)
        {
            if (opt.sample_dt_s > 0.0 && std::isfinite(opt.sample_dt_s))
            {
                return opt.sample_dt_s;
            }
            if (opt.max_samples <= 1)
            {
                return 0.0;
            }
            if (!(opt.duration_s > 0.0) || !std::isfinite(opt.duration_s))
            {
                return 0.0;
            }
            return opt.duration_s / static_cast<double>(opt.max_samples - 1);
        }

        inline double compute_celestial_dt_(const TrajectoryOptions &opt)
        {
            if (opt.celestial_dt_s > 0.0 && std::isfinite(opt.celestial_dt_s))
            {
                return opt.celestial_dt_s;
            }
            return compute_sample_dt_(opt);
        }

        inline void snapshot_states_(const std::vector<MassiveBody> &massive, std::vector<State> *out)
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

        inline CelestialEphemeris build_celestial_ephemeris_(const GameSimulation &sim, const TrajectoryOptions &opt)
        {
            CelestialEphemeris eph;
            if (!(opt.duration_s > 0.0) || !std::isfinite(opt.duration_s))
            {
                return eph;
            }
            if (opt.max_samples == 0)
            {
                return eph;
            }

            const double dt = compute_celestial_dt_(opt);
            if (!(dt > 0.0) || !std::isfinite(dt))
            {
                return eph;
            }

            std::vector<MassiveBody> massive = sim.massive_bodies();
            double t = sim.time_s();

            const double t_end = t + opt.duration_s;
            if (!std::isfinite(t_end) || !(t_end > t))
            {
                return eph;
            }

            std::vector<State> start_states;
            std::vector<State> end_states;

            while (t < t_end && eph.segments.size() < opt.max_samples)
            {
                double h = dt;
                if (t + h > t_end)
                {
                    h = t_end - t;
                }
                if (!(h > 0.0) || !std::isfinite(h))
                {
                    break;
                }

                snapshot_states_(massive, &start_states);
                symplectic4_step(massive, h, sim.config().gravitational_constant, sim.config().softening_length_m);
                snapshot_states_(massive, &end_states);

                eph.segments.push_back(CelestialEphemerisSegment{
                        .t0_s = t,
                        .dt_s = h,
                        .start = start_states,
                        .end = end_states,
                });
                t += h;
            }

            return eph;
        }

        inline std::optional<State> origin_state_at_(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                     const std::optional<std::size_t> origin_body_index,
                                                     const double t_s)
        {
            if (!origin_body_index.has_value())
            {
                return std::nullopt;
            }
            const std::size_t origin = *origin_body_index;
            if (origin >= sim.massive_bodies().size())
            {
                return std::nullopt;
            }
            if (!eph.empty())
            {
                return eph.body_state_at(origin, t_s);
            }
            return sim.massive_bodies()[origin].state;
        }

        inline TrajectorySample make_sample_(const double t_s, const State &s, const std::optional<State> &origin_state)
        {
            TrajectorySample out;
            out.t_s = t_s;
            out.position_m = s.position_m;
            out.velocity_mps = s.velocity_mps;
            if (origin_state.has_value())
            {
                out.position_m -= origin_state->position_m;
                out.velocity_mps -= origin_state->velocity_mps;
            }
            return out;
        }

        inline Vec3 gravity_accel_mps2_(const CelestialEphemeris &eph, const std::vector<MassiveBody> &bodies,
                                        const double G, const double softening_length_m, const double t_s,
                                        const Vec3 &pos_m)
        {
            Vec3 a{0.0, 0.0, 0.0};
            if (!(G > 0.0) || !std::isfinite(G))
            {
                return a;
            }
            const double eps2 = softening_length_m * softening_length_m;
            for (std::size_t i = 0; i < bodies.size(); ++i)
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
                a += (G * bodies[i].mass_kg) * inv_r3 * dr;
            }
            return a;
        }

        inline Vec3 burn_dir_inertial_unit_(const CelestialEphemeris &eph, const std::size_t primary_index,
                                            const double t_s, const Vec3 &sc_pos_m, const Vec3 &sc_vel_mps,
                                            const Vec3 &dir_rtn_unit)
        {
            const Vec3 r_primary_m = eph.body_position_at(primary_index, t_s);
            const Vec3 v_primary_mps = eph.body_velocity_at(primary_index, t_s);

            const Vec3 r_rel = sc_pos_m - r_primary_m;
            const Vec3 v_rel = sc_vel_mps - v_primary_mps;
            const RtnFrame f = compute_rtn_frame(r_rel, v_rel);

            const Vec3 dir_i = dir_rtn_unit.x * f.R + dir_rtn_unit.y * f.T + dir_rtn_unit.z * f.N;
            return normalized_or(dir_i, Vec3{0.0, 0.0, 0.0});
        }

        inline Spacecraft propagate_spacecraft_in_ephemeris_(const Spacecraft &sc0,
                                                             const std::vector<MassiveBody> &bodies,
                                                             const CelestialEphemeris &eph, const ManeuverPlan &plan,
                                                             const GameSimulation::Config &cfg, const double t0_s,
                                                             const double dt_s, const std::size_t spacecraft_index)
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
                    return gravity_accel_mps2_(eph, bodies, cfg.gravitational_constant, cfg.softening_length_m,
                                               t_eval_s, pos_m);
                };
                DOPRI5Stats stats{};
                const SpacecraftKinematics y1 =
                        dopri5_integrate_interval(t0_s, y0, dt_s, accel, cfg.spacecraft_integrator, &stats);
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

            auto primary_for = [&](const BurnSegment &seg, const double t_s, const Vec3 &pos_m) -> std::size_t {
                if (seg.primary_index < bodies.size())
                {
                    return seg.primary_index;
                }
                if (bodies.empty())
                {
                    return 0;
                }

                const double eps2 = cfg.softening_length_m * cfg.softening_length_m;
                std::size_t best = 0;
                double best_a = -1.0;
                for (std::size_t i = 0; i < bodies.size(); ++i)
                {
                    const Vec3 dr = eph.body_position_at(i, t_s) - pos_m;
                    const double r2 = glm::dot(dr, dr) + eps2;
                    if (!(r2 > 0.0) || !std::isfinite(r2))
                    {
                        continue;
                    }
                    const double amag = (cfg.gravitational_constant * bodies[i].mass_kg) / r2;
                    if (amag > best_a)
                    {
                        best_a = amag;
                        best = i;
                    }
                }
                return best;
            };

            while (remaining > 0.0)
            {
                const double boundary = next_burn_boundary_after(plan, spacecraft_index, t, t_end_s);
                double h = std::min(remaining, boundary - t);
                if (!(h > 0.0) || !std::isfinite(h))
                {
                    break;
                }

                const BurnSegment *seg = active_burn_at(plan, spacecraft_index, t);
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
                const std::size_t primary = (seg != nullptr) ? primary_for(*seg, t, y.position_m) : 0;
                const Vec3 dir_rtn = (seg != nullptr) ? seg->dir_rtn_unit : Vec3{0.0, 0.0, 0.0};

                auto accel = [&](double t_eval_s, const Vec3 &pos_m, const Vec3 &vel_mps) -> Vec3 {
                    Vec3 a = gravity_accel_mps2_(eph, bodies, cfg.gravitational_constant, cfg.softening_length_m,
                                                 t_eval_s, pos_m);
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

                    if (bodies.empty())
                    {
                        return a;
                    }

                    const Vec3 dir_i = burn_dir_inertial_unit_(eph, primary, t_eval_s, pos_m, vel_mps, dir_rtn);
                    const double dir2 = glm::dot(dir_i, dir_i);
                    if (!(dir2 > 0.0) || !std::isfinite(dir2))
                    {
                        return a;
                    }

                    a += (thrust_N / mass_kg) * dir_i;
                    return a;
                };

                DOPRI5Stats stats{};
                y = dopri5_integrate_interval(t, y, h, accel, cfg.spacecraft_integrator, &stats);
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

        inline std::vector<TrajectorySample> predict_body_trajectory_from_ephemeris_(const GameSimulation &sim,
                                                                                     const CelestialEphemeris &eph,
                                                                                     const std::size_t body_index,
                                                                                     const TrajectoryOptions &opt)
        {
            std::vector<TrajectorySample> out;
            if (opt.max_samples == 0)
            {
                return out;
            }
            if (body_index >= sim.massive_bodies().size())
            {
                return out;
            }
            if (eph.empty())
            {
                return out;
            }

            const double dt = compute_sample_dt_(opt);
            if (!(dt > 0.0) || !std::isfinite(dt))
            {
                return out;
            }

            const double t0 = sim.time_s();
            const double t_end = t0 + opt.duration_s;
            if (!std::isfinite(t_end) || !(t_end > t0))
            {
                return out;
            }

            out.reserve(std::min<std::size_t>(opt.max_samples, 1024));

            double t = t0;
            if (opt.include_start)
            {
                const State s = eph.body_state_at(body_index, t0);
                const std::optional<State> origin = origin_state_at_(sim, eph, opt.origin_body_index, t0);
                out.push_back(make_sample_(t0, s, origin));
            }

            while (t < t_end && out.size() < opt.max_samples)
            {
                double h = dt;
                if (t + h > t_end)
                {
                    h = t_end - t;
                }
                if (!(h > 0.0) || !std::isfinite(h))
                {
                    break;
                }
                t += h;

                if (!opt.include_end && !(t < t_end))
                {
                    break;
                }

                const State s = eph.body_state_at(body_index, t);
                const std::optional<State> origin = origin_state_at_(sim, eph, opt.origin_body_index, t);
                out.push_back(make_sample_(t, s, origin));
            }

            return out;
        }

        inline std::vector<TrajectorySample>
        predict_spacecraft_trajectory_from_ephemeris_(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                      const std::size_t spacecraft_index, const TrajectoryOptions &opt)
        {
            std::vector<TrajectorySample> out;
            if (opt.max_samples == 0)
            {
                return out;
            }
            if (spacecraft_index >= sim.spacecraft().size())
            {
                return out;
            }

            const double dt = compute_sample_dt_(opt);
            if (!(dt > 0.0) || !std::isfinite(dt))
            {
                return out;
            }

            const double t0 = sim.time_s();
            const double t_end = t0 + opt.duration_s;
            if (!std::isfinite(t_end) || !(t_end > t0))
            {
                return out;
            }

            out.reserve(std::min<std::size_t>(opt.max_samples, 1024));

            Spacecraft sc = sim.spacecraft()[spacecraft_index];
            double t = t0;

            if (opt.include_start)
            {
                const std::optional<State> origin = origin_state_at_(sim, eph, opt.origin_body_index, t);
                out.push_back(make_sample_(t, sc.state, origin));
            }

            while (t < t_end && out.size() < opt.max_samples)
            {
                double h = dt;
                if (t + h > t_end)
                {
                    h = t_end - t;
                }
                if (!(h > 0.0) || !std::isfinite(h))
                {
                    break;
                }

                sc = propagate_spacecraft_in_ephemeris_(sc, sim.massive_bodies(), eph, sim.maneuver_plan(),
                                                        sim.config(), t, h, spacecraft_index);
                t += h;

                if (!opt.include_end && !(t < t_end))
                {
                    break;
                }

                const std::optional<State> origin = origin_state_at_(sim, eph, opt.origin_body_index, t);
                out.push_back(make_sample_(t, sc.state, origin));
            }

            return out;
        }
    } // namespace detail

    inline CelestialEphemeris build_celestial_ephemeris(const GameSimulation &sim, const TrajectoryOptions &opt = {})
    {
        return detail::build_celestial_ephemeris_(sim, opt);
    }

    inline std::vector<TrajectorySample>
    predict_body_trajectory(const GameSimulation &sim, const std::size_t body_index, const TrajectoryOptions &opt = {})
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return detail::predict_body_trajectory_from_ephemeris_(sim, eph, body_index, opt);
    }

    inline std::vector<TrajectorySample> predict_spacecraft_trajectory(const GameSimulation &sim,
                                                                       const std::size_t spacecraft_index,
                                                                       const TrajectoryOptions &opt = {})
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return detail::predict_spacecraft_trajectory_from_ephemeris_(sim, eph, spacecraft_index, opt);
    }

    inline std::vector<TrajectorySample> predict_body_trajectory(const GameSimulation &sim,
                                                                 const CelestialEphemeris &eph,
                                                                 const std::size_t body_index,
                                                                 const TrajectoryOptions &opt = {})
    {
        return detail::predict_body_trajectory_from_ephemeris_(sim, eph, body_index, opt);
    }

    inline std::vector<TrajectorySample> predict_spacecraft_trajectory(const GameSimulation &sim,
                                                                       const CelestialEphemeris &eph,
                                                                       const std::size_t spacecraft_index,
                                                                       const TrajectoryOptions &opt = {})
    {
        return detail::predict_spacecraft_trajectory_from_ephemeris_(sim, eph, spacecraft_index, opt);
    }

} // namespace orbitsim
