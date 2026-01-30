#pragma once

#include "orbitsim/detail/spacecraft_propagation.hpp"
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

                snapshot_states(massive, &start_states);
                symplectic4_step(massive, h, sim.config().gravitational_constant, sim.config().softening_length_m);
                snapshot_states(massive, &end_states);

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

                sc = propagate_spacecraft_in_ephemeris(
                        sc, sim.massive_bodies(), eph, sim.maneuver_plan(), sim.config().gravitational_constant,
                        sim.config().softening_length_m, sim.config().spacecraft_integrator, t, h, spacecraft_index);
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
