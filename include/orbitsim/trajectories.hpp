#pragma once

#include "orbitsim/detail/spacecraft_propagation.hpp"
#include "orbitsim/events.hpp"
#include "orbitsim/game_sim.hpp"
#include "orbitsim/nodes.hpp"
#include "orbitsim/time_utils.hpp"
#include "orbitsim/trajectory_types.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

namespace orbitsim
{

    /**
     * @brief Configuration for trajectory prediction and sampling.
     *
     * Controls duration and sample intervals for trajectory prediction functions.
     */

    struct TrajectoryOptions
    {
        double duration_s{3600.0};  ///< Prediction duration from sim.time_s()

        double sample_dt_s{10.0};   ///< Sample interval; if <= 0, derived from duration/(max_samples-1)

        /// Spacecraft sample interval; if > 0, overrides sample_dt_s for spacecraft.
        /// Useful for smooth orbit-line rendering of fast local orbits (e.g. LEO).
        double spacecraft_sample_dt_s{0.0};

        /// Massive body integration step; if <= 0, uses sample_dt_s.
        double celestial_dt_s{0.0};

        std::size_t max_samples{2048};  ///< Hard cap to prevent unbounded allocations

        bool include_start{true};   ///< Include sample at t0
        bool include_end{true};     ///< Include sample at t0 + duration
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

        inline TrajectorySample make_sample_(const double t_s, const State &s)
        {
            return TrajectorySample{.t_s = t_s, .position_m = s.position_m, .velocity_mps = s.velocity_mps};
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
                out.push_back(make_sample_(t0, s));
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
                out.push_back(make_sample_(t, s));
            }

            return out;
        }

        inline std::vector<TrajectorySample>
        predict_spacecraft_trajectory_from_ephemeris_by_id_(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                            const SpacecraftId spacecraft_id,
                                                            const TrajectoryOptions &opt)
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
                out.push_back(make_sample_(t, sc.state));
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

                out.push_back(make_sample_(t, sc.state));
            }

            return out;
        }

        inline bool is_terminal_event_(const Event &e)
        {
            return (e.type == EventType::Impact) && (e.crossing == Crossing::Enter);
        }

        inline std::vector<Event> predict_spacecraft_events_from_ephemeris_by_id_(const GameSimulation &sim,
                                                                                  const CelestialEphemeris &eph,
                                                                                  const SpacecraftId spacecraft_id,
                                                                                  const TrajectoryOptions &traj_opt,
                                                                                  const EventOptions &event_opt)
        {
            std::vector<Event> out;
            if (traj_opt.duration_s <= 0.0 || !std::isfinite(traj_opt.duration_s))
            {
                return out;
            }
            if (eph.empty())
            {
                return out;
            }
            const Spacecraft *sc_ptr = sim.spacecraft_by_id(spacecraft_id);
            if (sc_ptr == nullptr)
            {
                return out;
            }

            const double t0 = sim.time_s();
            const double t_end = t0 + traj_opt.duration_s;
            if (!std::isfinite(t_end) || !(t_end > t0))
            {
                return out;
            }

            Spacecraft sc = *sc_ptr;
            double t = t0;

            for (const auto &seg: eph.segments)
            {
                const double seg_start = seg.t0_s;
                const double seg_end = seg.t0_s + seg.dt_s;
                if (!(seg_end > seg_start) || !std::isfinite(seg_start) || !std::isfinite(seg_end))
                {
                    continue;
                }
                if (!(seg_end > t))
                {
                    continue;
                }
                if (!(seg_start < t_end))
                {
                    break;
                }

                // Align to the segment time if needed.
                if (t < seg_start)
                {
                    t = seg_start;
                }

                double remaining = std::min(seg_end, t_end) - t;
                if (!(remaining > 0.0) || !std::isfinite(remaining))
                {
                    continue;
                }

                auto propagate_sc = [&](const Spacecraft &sc_start, const double t0_s,
                                        const double dt_sc_s) -> Spacecraft {
                    return propagate_spacecraft_in_ephemeris(sc_start, sim.massive_bodies(), seg, sim.maneuver_plan(),
                                                             sim.config().gravitational_constant,
                                                             sim.config().softening_length_m,
                                                             sim.config().spacecraft_integrator, t0_s, dt_sc_s);
                };

                while (remaining > 0.0)
                {
                    const std::optional<Event> e = find_earliest_event_in_interval(
                            sim.massive_bodies(), seg, sc, t, remaining, sim.maneuver_plan(), event_opt, propagate_sc);
                    if (!e.has_value())
                    {
                        sc = propagate_sc(sc, t, remaining);
                        t += remaining;
                        remaining = 0.0;
                        break;
                    }

                    double dt_event = e->t_event_s - t;
                    if (!std::isfinite(dt_event) || dt_event < 0.0)
                    {
                        dt_event = 0.0;
                    }
                    if (dt_event > remaining)
                    {
                        dt_event = remaining;
                    }

                    sc = propagate_sc(sc, t, dt_event);
                    t += dt_event;

                    Event logged = *e;
                    logged.t_event_s = t;
                    out.push_back(logged);

                    if (is_terminal_event_(logged))
                    {
                        return out;
                    }

                    // Nudge forward to avoid re-detecting the same boundary at the exact same time.
                    double dt_eps = std::max(0.0, event_opt.time_tol_s);
                    if (!(dt_eps > 0.0) || !std::isfinite(dt_eps))
                    {
                        dt_eps = 1e-6;
                    }
                    dt_eps = std::min(dt_eps, std::min(seg_end, t_end) - t);
                    if (!(dt_eps > 0.0))
                    {
                        remaining = 0.0;
                        break;
                    }
                    sc = propagate_sc(sc, t, dt_eps);
                    t += dt_eps;
                    remaining = std::min(seg_end, t_end) - t;
                    if (!(remaining > 0.0) || !std::isfinite(remaining))
                    {
                        remaining = 0.0;
                        break;
                    }
                }

                if (!(t < t_end))
                {
                    break;
                }
            }

            return out;
        }

        inline std::optional<Vec3> spin_axis_unit_(const MassiveBody &body)
        {
            const Vec3 axis = normalized_or(body.state.spin.axis, Vec3{0.0, 0.0, 0.0});
            const double a2 = glm::dot(axis, axis);
            if (!(a2 > 0.0) || !std::isfinite(a2))
            {
                return std::nullopt;
            }
            return axis;
        }

        inline std::optional<Vec3> orbit_plane_normal_unit_(const Vec3 &r_rel_m, const Vec3 &v_rel_mps)
        {
            const Vec3 h = glm::cross(r_rel_m, v_rel_mps);
            const double h2 = glm::dot(h, h);
            if (!(h2 > 1e-24) || !std::isfinite(h2))
            {
                return std::nullopt;
            }
            return h / std::sqrt(h2);
        }

        template<class EphemerisLike>
        inline std::optional<std::size_t> body_index_for_id_(const EphemerisLike &eph, const BodyId body_id)
        {
            std::size_t idx = 0;
            if constexpr (requires { eph.body_index_for_id(body_id, &idx); })
            {
                if (eph.body_index_for_id(body_id, &idx))
                {
                    return idx;
                }
            }
            return std::nullopt;
        }

        inline std::vector<NodeEvent> predict_plane_nodes_from_ephemeris_by_id_(
                const GameSimulation &sim, const CelestialEphemeris &eph, const SpacecraftId spacecraft_id,
                const BodyId primary_body_id, const Vec3 &plane_normal_unit_i, const SpacecraftId target_spacecraft_id,
                const TrajectoryOptions &traj_opt, const EventOptions &opt)
        {
            std::vector<NodeEvent> out;
            if (!(traj_opt.duration_s > 0.0) || !std::isfinite(traj_opt.duration_s))
            {
                return out;
            }
            if (eph.empty())
            {
                return out;
            }
            const Spacecraft *sc_ptr = sim.spacecraft_by_id(spacecraft_id);
            if (sc_ptr == nullptr)
            {
                return out;
            }
            const MassiveBody *primary_ptr = sim.body_by_id(primary_body_id);
            if (primary_ptr == nullptr)
            {
                return out;
            }

            const double t0 = sim.time_s();
            const double t_end = t0 + traj_opt.duration_s;
            if (!std::isfinite(t_end) || !(t_end > t0))
            {
                return out;
            }

            Spacecraft sc = *sc_ptr;
            double t = t0;

            for (const auto &seg: eph.segments)
            {
                const double seg_start = seg.t0_s;
                const double seg_end = seg.t0_s + seg.dt_s;
                if (!(seg_end > seg_start) || !std::isfinite(seg_start) || !std::isfinite(seg_end))
                {
                    continue;
                }
                if (!(seg_end > t))
                {
                    continue;
                }
                if (!(seg_start < t_end))
                {
                    break;
                }

                if (t < seg_start)
                {
                    t = seg_start;
                }

                double remaining = std::min(seg_end, t_end) - t;
                if (!(remaining > 0.0) || !std::isfinite(remaining))
                {
                    continue;
                }

                auto propagate_sc = [&](const Spacecraft &sc_start, const double t0_s,
                                        const double dt_sc_s) -> Spacecraft {
                    return propagate_spacecraft_in_ephemeris(sc_start, sim.massive_bodies(), seg, sim.maneuver_plan(),
                                                             sim.config().gravitational_constant,
                                                             sim.config().softening_length_m,
                                                             sim.config().spacecraft_integrator, t0_s, dt_sc_s);
                };

                while (remaining > 0.0)
                {
                    const std::optional<NodeEvent> e = find_earliest_plane_node_in_interval_(
                            sim.massive_bodies(), seg, sc, t, remaining, primary_body_id, plane_normal_unit_i,
                            target_spacecraft_id, opt, propagate_sc);

                    if (!e.has_value())
                    {
                        sc = propagate_sc(sc, t, remaining);
                        t += remaining;
                        remaining = 0.0;
                        break;
                    }

                    double dt_event = e->t_event_s - t;
                    if (!std::isfinite(dt_event) || dt_event < 0.0)
                    {
                        dt_event = 0.0;
                    }
                    if (dt_event > remaining)
                    {
                        dt_event = remaining;
                    }

                    sc = propagate_sc(sc, t, dt_event);
                    t += dt_event;

                    NodeEvent logged = *e;
                    logged.t_event_s = t;
                    out.push_back(logged);

                    double dt_eps = std::max(0.0, opt.time_tol_s);
                    if (!(dt_eps > 0.0) || !std::isfinite(dt_eps))
                    {
                        dt_eps = 1e-6;
                    }
                    dt_eps = std::min(dt_eps, std::min(seg_end, t_end) - t);
                    if (!(dt_eps > 0.0))
                    {
                        remaining = 0.0;
                        break;
                    }

                    sc = propagate_sc(sc, t, dt_eps);
                    t += dt_eps;
                    remaining = std::min(seg_end, t_end) - t;
                    if (!(remaining > 0.0) || !std::isfinite(remaining))
                    {
                        remaining = 0.0;
                        break;
                    }
                }

                if (!(t < t_end))
                {
                    break;
                }
            }

            return out;
        }
    } // namespace detail

    /**
     * @brief Build precomputed ephemeris for all massive bodies.
     *
     * Steps the N-body simulation forward and stores piecewise-Hermite segments
     * for efficient interpolation during spacecraft trajectory prediction.
     */
    inline CelestialEphemeris build_celestial_ephemeris(const GameSimulation &sim, const TrajectoryOptions &opt = {})
    {
        return detail::build_celestial_ephemeris_(sim, opt);
    }

    /**
     * @brief Predict future trajectory of a massive body.
     *
     * Builds ephemeris internally and samples body positions at regular intervals.
     * Use the overload with pre-built ephemeris when predicting multiple bodies.
     */
    inline std::vector<TrajectorySample> predict_body_trajectory(const GameSimulation &sim, const BodyId body_id,
                                                                 const TrajectoryOptions &opt = {})
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return detail::predict_body_trajectory_from_ephemeris_(sim, eph, body_id, opt);
    }

    /**
     * @brief Predict future trajectory of a spacecraft.
     *
     * Propagates spacecraft through precomputed massive body ephemeris,
     * applying scheduled maneuvers from the simulation's maneuver plan.
     */
    inline std::vector<TrajectorySample> predict_spacecraft_trajectory(const GameSimulation &sim,
                                                                       const SpacecraftId spacecraft_id,
                                                                       const TrajectoryOptions &opt = {})
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return detail::predict_spacecraft_trajectory_from_ephemeris_by_id_(sim, eph, spacecraft_id, opt);
    }

    /** @brief Predict body trajectory using pre-built ephemeris. */
    inline std::vector<TrajectorySample> predict_body_trajectory(const GameSimulation &sim,
                                                                 const CelestialEphemeris &eph, const BodyId body_id,
                                                                 const TrajectoryOptions &opt = {})
    {
        return detail::predict_body_trajectory_from_ephemeris_(sim, eph, body_id, opt);
    }

    /** @brief Predict spacecraft trajectory using pre-built ephemeris. */
    inline std::vector<TrajectorySample> predict_spacecraft_trajectory(const GameSimulation &sim,
                                                                       const CelestialEphemeris &eph,
                                                                       const SpacecraftId spacecraft_id,
                                                                       const TrajectoryOptions &opt = {})
    {
        return detail::predict_spacecraft_trajectory_from_ephemeris_by_id_(sim, eph, spacecraft_id, opt);
    }

    /**
     * @brief Predict future events (SOI crossings, impacts, etc.) for a spacecraft.
     *
     * Propagates spacecraft and detects boundary crossings. Returns events
     * in chronological order. Stops early on terminal events (e.g. impact).
     */
    inline std::vector<Event> predict_spacecraft_events(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                        const SpacecraftId spacecraft_id, const TrajectoryOptions &opt,
                                                        const EventOptions &event_opt)
    {
        return detail::predict_spacecraft_events_from_ephemeris_by_id_(sim, eph, spacecraft_id, opt, event_opt);
    }

    inline std::vector<Event> predict_spacecraft_events(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                        const SpacecraftId spacecraft_id,
                                                        const TrajectoryOptions &opt = {})
    {
        return detail::predict_spacecraft_events_from_ephemeris_by_id_(sim, eph, spacecraft_id, opt,
                                                                       sim.config().events);
    }

    inline std::vector<Event> predict_spacecraft_events(const GameSimulation &sim, const SpacecraftId spacecraft_id,
                                                        const TrajectoryOptions &opt, const EventOptions &event_opt)
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return detail::predict_spacecraft_events_from_ephemeris_by_id_(sim, eph, spacecraft_id, opt, event_opt);
    }

    inline std::vector<Event> predict_spacecraft_events(const GameSimulation &sim, const SpacecraftId spacecraft_id,
                                                        const TrajectoryOptions &opt = {})
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return detail::predict_spacecraft_events_from_ephemeris_by_id_(sim, eph, spacecraft_id, opt,
                                                                       sim.config().events);
    }

    /**
     * @brief Predict ascending/descending node crossings through a body's equatorial plane.
     *
     * The equatorial plane is defined by the body's spin axis. Useful for
     * planning inclination change maneuvers.
     */
    inline std::vector<NodeEvent> predict_equatorial_nodes(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                           const SpacecraftId spacecraft_id,
                                                           const BodyId primary_body_id, const TrajectoryOptions &opt,
                                                           const EventOptions &node_opt)
    {
        const MassiveBody *primary = sim.body_by_id(primary_body_id);
        if (primary == nullptr)
        {
            return {};
        }
        const std::optional<Vec3> axis = detail::spin_axis_unit_(*primary);
        if (!axis.has_value())
        {
            return {};
        }
        return detail::predict_plane_nodes_from_ephemeris_by_id_(sim, eph, spacecraft_id, primary_body_id, *axis,
                                                                 kInvalidSpacecraftId, opt, node_opt);
    }

    inline std::vector<NodeEvent> predict_equatorial_nodes(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                           const SpacecraftId spacecraft_id,
                                                           const BodyId primary_body_id,
                                                           const TrajectoryOptions &opt = {})
    {
        return predict_equatorial_nodes(sim, eph, spacecraft_id, primary_body_id, opt, sim.config().events);
    }

    inline std::vector<NodeEvent> predict_equatorial_nodes(const GameSimulation &sim, const SpacecraftId spacecraft_id,
                                                           const BodyId primary_body_id, const TrajectoryOptions &opt,
                                                           const EventOptions &node_opt)
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return predict_equatorial_nodes(sim, eph, spacecraft_id, primary_body_id, opt, node_opt);
    }

    inline std::vector<NodeEvent> predict_equatorial_nodes(const GameSimulation &sim, const SpacecraftId spacecraft_id,
                                                           const BodyId primary_body_id,
                                                           const TrajectoryOptions &opt = {})
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return predict_equatorial_nodes(sim, eph, spacecraft_id, primary_body_id, opt, sim.config().events);
    }

    /**
     * @brief Predict node crossings through another spacecraft's orbital plane.
     *
     * The target plane is defined by the target spacecraft's angular momentum
     * vector at the current time. Useful for rendezvous plane-change planning.
     */
    inline std::vector<NodeEvent> predict_target_plane_nodes(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                             const SpacecraftId spacecraft_id,
                                                             const SpacecraftId target_spacecraft_id,
                                                             const BodyId primary_body_id, const TrajectoryOptions &opt,
                                                             const EventOptions &node_opt)
    {
        const Spacecraft *target_ptr = sim.spacecraft_by_id(target_spacecraft_id);
        const MassiveBody *primary_ptr = sim.body_by_id(primary_body_id);
        if (target_ptr == nullptr || primary_ptr == nullptr)
        {
            return {};
        }
        const double t0 = sim.time_s();
        const State primary_state = eph.empty() ? primary_ptr->state : eph.body_state_at_by_id(primary_body_id, t0);
        const Vec3 r_rel = target_ptr->state.position_m - primary_state.position_m;
        const Vec3 v_rel = target_ptr->state.velocity_mps - primary_state.velocity_mps;
        const std::optional<Vec3> n = detail::orbit_plane_normal_unit_(r_rel, v_rel);
        if (!n.has_value())
        {
            return {};
        }
        return detail::predict_plane_nodes_from_ephemeris_by_id_(sim, eph, spacecraft_id, primary_body_id, *n,
                                                                 target_spacecraft_id, opt, node_opt);
    }

    inline std::vector<NodeEvent> predict_target_plane_nodes(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                             const SpacecraftId spacecraft_id,
                                                             const SpacecraftId target_spacecraft_id,
                                                             const BodyId primary_body_id,
                                                             const TrajectoryOptions &opt = {})
    {
        return predict_target_plane_nodes(sim, eph, spacecraft_id, target_spacecraft_id, primary_body_id, opt,
                                          sim.config().events);
    }

    inline std::vector<NodeEvent> predict_target_plane_nodes(const GameSimulation &sim,
                                                             const SpacecraftId spacecraft_id,
                                                             const SpacecraftId target_spacecraft_id,
                                                             const BodyId primary_body_id, const TrajectoryOptions &opt,
                                                             const EventOptions &node_opt)
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return predict_target_plane_nodes(sim, eph, spacecraft_id, target_spacecraft_id, primary_body_id, opt,
                                          node_opt);
    }

    inline std::vector<NodeEvent> predict_target_plane_nodes(const GameSimulation &sim,
                                                             const SpacecraftId spacecraft_id,
                                                             const SpacecraftId target_spacecraft_id,
                                                             const BodyId primary_body_id,
                                                             const TrajectoryOptions &opt = {})
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return predict_target_plane_nodes(sim, eph, spacecraft_id, target_spacecraft_id, primary_body_id, opt,
                                          sim.config().events);
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

        /// @brief Implicit conversion to TrajectoryOptions.
        operator TrajectoryOptions() const { return opt_; }

        /// @brief Explicit conversion to TrajectoryOptions.
        TrajectoryOptions build() const { return opt_; }

    private:
        TrajectoryOptions opt_{};
    };

    /// @brief Start building TrajectoryOptions with the fluent interface.
    inline TrajectoryOptionsBuilder trajectory_options() { return TrajectoryOptionsBuilder{}; }

} // namespace orbitsim
