#pragma once

#include "orbitsim/detail/spacecraft_propagation.hpp"
#include "orbitsim/synodic.hpp"
#include "orbitsim/game_sim.hpp"
#include "orbitsim/time_utils.hpp"
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

        // If > 0, spacecraft trajectories use this dt instead of sample_dt_s.
        // This is useful for smooth orbit-line rendering of fast local orbits (e.g. LEO)
        // while keeping body sampling coarse for long-horizon predictions.
        double spacecraft_sample_dt_s{0.0};

        // If > 0, uses this dt to step massive bodies when building a prediction ephemeris. Otherwise uses sample_dt_s
        // (or the derived sample dt).
        double celestial_dt_s{0.0};

        // Hard cap to protect against unbounded allocations.
        std::size_t max_samples{2048};

        bool include_start{true};
        bool include_end{true};

        // If set, returned state is relative to this body's state at the same time.
        std::optional<BodyId> origin_body_id{};
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

        inline double compute_spacecraft_sample_dt_(const TrajectoryOptions &opt)
        {
            if (opt.spacecraft_sample_dt_s > 0.0 && std::isfinite(opt.spacecraft_sample_dt_s))
            {
                return opt.spacecraft_sample_dt_s;
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

            std::vector<BodyId> body_ids;
            body_ids.reserve(sim.massive_bodies().size());
            for (const auto &b: sim.massive_bodies())
            {
                body_ids.push_back(b.id);
            }
            eph.set_body_ids(std::move(body_ids));

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
                                                     const std::optional<BodyId> origin_body_id,
                                                     const double t_s)
        {
            if (!origin_body_id.has_value())
            {
                return std::nullopt;
            }
            const BodyId origin = *origin_body_id;
            const MassiveBody *body = sim.body_by_id(origin);
            if (body == nullptr)
            {
                return std::nullopt;
            }
            if (!eph.empty())
            {
                return eph.body_state_at_by_id(origin, t_s);
            }
            return body->state;
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
                                                                                     const BodyId body_id,
                                                                                     const TrajectoryOptions &opt)
        {
            std::vector<TrajectorySample> out;
            if (opt.max_samples == 0)
            {
                return out;
            }
            if (eph.empty())
            {
                return out;
            }
            std::size_t body_index = 0;
            if (!eph.body_index_for_id(body_id, &body_index))
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
                const std::optional<State> origin = origin_state_at_(sim, eph, opt.origin_body_id, t0);
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
                const std::optional<State> origin = origin_state_at_(sim, eph, opt.origin_body_id, t);
                out.push_back(make_sample_(t, s, origin));
            }

            return out;
        }

        inline std::vector<TrajectorySample>
        predict_spacecraft_trajectory_from_ephemeris_by_id_(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                           const SpacecraftId spacecraft_id, const TrajectoryOptions &opt)
        {
            std::vector<TrajectorySample> out;
            if (opt.max_samples == 0)
            {
                return out;
            }
            const Spacecraft *sc_ptr = sim.spacecraft_by_id(spacecraft_id);
            if (sc_ptr == nullptr)
            {
                return out;
            }

            const double dt = compute_spacecraft_sample_dt_(opt);
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

            Spacecraft sc = *sc_ptr;
            double t = t0;

            if (opt.include_start)
            {
                const std::optional<State> origin = origin_state_at_(sim, eph, opt.origin_body_id, t);
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
                        sim.config().softening_length_m, sim.config().spacecraft_integrator, t, h);
                t += h;

                if (!opt.include_end && !(t < t_end))
                {
                    break;
                }

                const std::optional<State> origin = origin_state_at_(sim, eph, opt.origin_body_id, t);
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
    predict_body_trajectory(const GameSimulation &sim, const BodyId body_id, const TrajectoryOptions &opt = {})
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return detail::predict_body_trajectory_from_ephemeris_(sim, eph, body_id, opt);
    }

    inline std::vector<TrajectorySample> predict_spacecraft_trajectory(const GameSimulation &sim,
                                                                       const SpacecraftId spacecraft_id,
                                                                       const TrajectoryOptions &opt = {})
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return detail::predict_spacecraft_trajectory_from_ephemeris_by_id_(sim, eph, spacecraft_id, opt);
    }

    inline std::vector<TrajectorySample> predict_body_trajectory(const GameSimulation &sim,
                                                                 const CelestialEphemeris &eph,
                                                                 const BodyId body_id,
                                                                 const TrajectoryOptions &opt = {})
    {
        return detail::predict_body_trajectory_from_ephemeris_(sim, eph, body_id, opt);
    }

    inline std::vector<TrajectorySample> predict_spacecraft_trajectory(const GameSimulation &sim,
                                                                 const CelestialEphemeris &eph,
                                                                 const SpacecraftId spacecraft_id,
                                                                 const TrajectoryOptions &opt = {})
    {
        return detail::predict_spacecraft_trajectory_from_ephemeris_by_id_(sim, eph, spacecraft_id, opt);
    }

    // -------------------------------------------------------------------------
    // TrajectoryOptionsBuilder: fluent interface for creating TrajectoryOptions
    // -------------------------------------------------------------------------

    /// @brief Fluent builder for creating TrajectoryOptions objects.
    /// @example
    ///   auto opts = trajectory_options()
    ///       .duration(days(20.0))
    ///       .sample_dt(minutes(10.0))
    ///       .celestial_dt(minutes(5.0))
    ///       .max_samples(100'000);
    class TrajectoryOptionsBuilder
    {
    public:
        TrajectoryOptionsBuilder() = default;

        /// @brief Set the trajectory duration [s].
        TrajectoryOptionsBuilder &duration(const double duration_s)
        {
            opt_.duration_s = duration_s;
            return *this;
        }

        /// @brief Set the sample time step [s].
        TrajectoryOptionsBuilder &sample_dt(const double sample_dt_s)
        {
            opt_.sample_dt_s = sample_dt_s;
            return *this;
        }

        /// @brief Set the spacecraft-specific sample time step [s].
        TrajectoryOptionsBuilder &spacecraft_sample_dt(const double spacecraft_sample_dt_s)
        {
            opt_.spacecraft_sample_dt_s = spacecraft_sample_dt_s;
            return *this;
        }

        /// @brief Set the celestial body integration time step [s].
        TrajectoryOptionsBuilder &celestial_dt(const double celestial_dt_s)
        {
            opt_.celestial_dt_s = celestial_dt_s;
            return *this;
        }

        /// @brief Set the maximum number of samples.
        TrajectoryOptionsBuilder &max_samples(const std::size_t max_samples)
        {
            opt_.max_samples = max_samples;
            return *this;
        }

        /// @brief Set whether to include the start point.
        TrajectoryOptionsBuilder &include_start(const bool include)
        {
            opt_.include_start = include;
            return *this;
        }

        /// @brief Set whether to include the end point.
        TrajectoryOptionsBuilder &include_end(const bool include)
        {
            opt_.include_end = include;
            return *this;
        }

        /// @brief Set the origin body for relative coordinates.
        TrajectoryOptionsBuilder &origin(const BodyId origin_body_id)
        {
            opt_.origin_body_id = origin_body_id;
            return *this;
        }

        /// @brief Clear the origin body (use absolute coordinates).
        TrajectoryOptionsBuilder &no_origin()
        {
            opt_.origin_body_id = std::nullopt;
            return *this;
        }

        /// @brief Implicit conversion to TrajectoryOptions.
        operator TrajectoryOptions() const { return opt_; }

        /// @brief Explicit conversion to TrajectoryOptions.
        TrajectoryOptions build() const { return opt_; }

    private:
        TrajectoryOptions opt_{};
    };

    /// @brief Start building TrajectoryOptions with the fluent interface.
    inline TrajectoryOptionsBuilder trajectory_options() { return TrajectoryOptionsBuilder{}; }

    /// @brief Convert a single inertial trajectory sample into a synodic rotating frame.
    /// @note The input sample must be in the same inertial frame as the body states used to construct the SynodicFrame.
    inline TrajectorySample inertial_sample_to_synodic(const TrajectorySample &sample_in, const SynodicFrame &frame)
    {
        const State s_in = make_state(sample_in.position_m, sample_in.velocity_mps);
        const State s_rot = inertial_state_to_frame(s_in, frame);
        return TrajectorySample{.t_s = sample_in.t_s, .position_m = s_rot.position_m, .velocity_mps = s_rot.velocity_mps};
    }

    /// @brief Convert inertial trajectory samples into the time-varying synodic frame of bodies (A,B).
    /// @note The input samples must be inertial (i.e., not already offset by TrajectoryOptions::origin_body_id).
    inline std::vector<TrajectorySample> trajectory_to_synodic(const std::vector<TrajectorySample> &samples_in,
                                                               const CelestialEphemeris &eph,
                                                               const MassiveBody &body_a,
                                                               const MassiveBody &body_b)
    {
        std::vector<TrajectorySample> out;
        out.reserve(samples_in.size());

        for (const auto &s: samples_in)
        {
            const std::optional<SynodicFrame> frame = make_synodic_frame_at(eph, body_a, body_b, s.t_s);
            if (!frame.has_value())
            {
                return {};
            }
            out.push_back(inertial_sample_to_synodic(s, *frame));
        }

        return out;
    }

} // namespace orbitsim
