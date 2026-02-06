#pragma once

#include "orbitsim/detail/spacecraft_propagation.hpp"
#include "orbitsim/events.hpp"
#include "orbitsim/game_sim.hpp"
#include "orbitsim/nodes.hpp"
#include "orbitsim/spacecraft_state_cache.hpp"
#include "orbitsim/trajectories.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

namespace orbitsim
{

    namespace detail
    {
        inline bool is_terminal_event_(const Event &e)
        {
            return e.type == EventType::Impact && e.crossing == Crossing::Enter;
        }

        inline std::optional<Vec3> spin_axis_unit_(const MassiveBody &body)
        {
            const Vec3 axis = normalized_or(body.state.spin.axis, Vec3{0.0, 0.0, 0.0});
            if (!(glm::dot(axis, axis) > 0.0))
                return std::nullopt;
            return axis;
        }

        inline std::optional<Vec3> orbit_plane_normal_unit_(const Vec3 &r_rel_m, const Vec3 &v_rel_mps)
        {
            const Vec3 h = glm::cross(r_rel_m, v_rel_mps);
            const double h2 = glm::dot(h, h);
            if (!(h2 > 1e-24))
                return std::nullopt;
            return h / std::sqrt(h2);
        }

        /// Compute nudge timestep to avoid re-detecting the same event boundary.
        inline double event_nudge_dt_(const EventOptions &opt)
        {
            return (opt.time_tol_s > 0.0) ? opt.time_tol_s : 1e-6;
        }

        /// Build a SpacecraftStateCache from common simulation parameters.
        inline SpacecraftStateCache<CelestialEphemeris> make_sc_cache_(
            const GameSimulation &sim, const CelestialEphemeris &eph,
            const double t0, const double t_end, const double lookup_dt)
        {
            return SpacecraftStateCache<CelestialEphemeris>(
                sim.massive_bodies(), eph, sim.maneuver_plan(),
                sim.config().gravitational_constant, sim.config().softening_length_m,
                sim.config().spacecraft_integrator, t0, t_end,
                [&sim](const SpacecraftId id) { return sim.spacecraft_by_id(id); },
                SpacecraftStateCache<CelestialEphemeris>::Options{.lookup_dt_s = lookup_dt});
        }

        // ---- Common ephemeris event-scan loop ----

        /// Walks ephemeris segments, propagating a spacecraft and collecting
        /// events detected by `find_event`.  `is_terminal` can stop early.
        ///
        /// FindEventFn : (seg, sc, t, remaining, propagate_sc, *sc1) -> optional<EventT>
        /// IsTerminalFn: (EventT) -> bool
        template<class EventT, class FindEventFn, class IsTerminalFn>
        inline std::vector<EventT> scan_ephemeris_for_events_(
            const GameSimulation &sim, const CelestialEphemeris &eph,
            const SpacecraftId spacecraft_id, const TrajectoryOptions &traj_opt,
            const EventOptions &opt, FindEventFn find_event, IsTerminalFn is_terminal)
        {
            std::vector<EventT> out;

            const Spacecraft *sc_ptr = sim.spacecraft_by_id(spacecraft_id);
            if (!sc_ptr || eph.empty() || !(traj_opt.duration_s > 0.0))
                return out;

            const double t0 = sim.time_s();
            const double t_end = t0 + traj_opt.duration_s;
            Spacecraft sc = *sc_ptr;
            double t = t0;

            auto sc_cache = make_sc_cache_(sim, eph, t0, t_end, compute_spacecraft_lookup_dt_(traj_opt));
            const SpacecraftStateLookup sc_lookup = sc_cache.lookup();
            const double nudge = event_nudge_dt_(opt);

            for (const auto &seg : eph.segments)
            {
                const double seg_end = seg.t0_s + seg.dt_s;
                if (seg_end <= t) continue;
                if (seg.t0_s >= t_end) break;

                t = std::max(t, seg.t0_s);
                double remaining = std::min(seg_end, t_end) - t;
                if (remaining <= 0.0) continue;

                auto propagate_sc = [&](const Spacecraft &sc_start, const double t0_s,
                                        const double dt_sc_s) -> Spacecraft {
                    return propagate_spacecraft_in_ephemeris(
                        sc_start, sim.massive_bodies(), seg, sim.maneuver_plan(),
                        sim.config().gravitational_constant, sim.config().softening_length_m,
                        sim.config().spacecraft_integrator, t0_s, dt_sc_s, sc_lookup);
                };

                while (remaining > 0.0)
                {
                    Spacecraft sc1{};
                    const auto e = find_event(seg, sc, t, remaining, propagate_sc, &sc1);

                    if (!e)
                    {
                        sc = sc1;
                        t += remaining;
                        break;
                    }

                    double dt_event = e->t_event_s - t;
                    if (!(dt_event >= 0.0)) dt_event = 0.0;
                    if (dt_event > remaining) dt_event = remaining;

                    sc = propagate_sc(sc, t, dt_event);
                    t += dt_event;

                    EventT logged = *e;
                    logged.t_event_s = t;
                    out.push_back(logged);

                    if (is_terminal(logged))
                        return out;

                    const double step = std::min(nudge, std::min(seg_end, t_end) - t);
                    if (step <= 0.0) break;

                    sc = propagate_sc(sc, t, step);
                    t += step;
                    remaining = std::min(seg_end, t_end) - t;
                }

                if (t >= t_end) break;
            }

            return out;
        }

        // ---- Generic spacecraft event prediction (SOI, atmosphere, impact) ----

        inline std::vector<Event> predict_spacecraft_events_from_ephemeris_by_id_(
            const GameSimulation &sim, const CelestialEphemeris &eph,
            const SpacecraftId spacecraft_id, const TrajectoryOptions &traj_opt,
            const EventOptions &event_opt)
        {
            return scan_ephemeris_for_events_<Event>(
                sim, eph, spacecraft_id, traj_opt, event_opt,
                [&](const auto &seg, const Spacecraft &sc, double t, double remaining,
                    auto propagate_sc, Spacecraft *sc1) {
                    return find_earliest_event_in_interval(
                        sim.massive_bodies(), seg, sc, t, remaining, sim.maneuver_plan(),
                        event_opt, propagate_sc, sc1);
                },
                [](const Event &e) { return is_terminal_event_(e); });
        }

        // ---- Apsis event detection (periapsis / apoapsis) ----

        template<class EphemerisLike, class Propagator>
        inline std::optional<ApsisEvent> find_earliest_apsis_event_in_interval_(
            const std::vector<MassiveBody> &bodies, const EphemerisLike &eph,
            const Spacecraft &sc0, const double t0_s, const double dt_s,
            const BodyId primary_body_id, const EventOptions &opt,
            Propagator propagate_sc, Spacecraft *out_sc1 = nullptr)
        {
            if (!(dt_s > 0.0) || opt.max_bisect_iters <= 0)
                return std::nullopt;

            const auto primary_idx = body_index_for_id(bodies, primary_body_id);
            if (!primary_idx)
                return std::nullopt;

            const std::size_t pi = *primary_idx;
            const double t1_s = t0_s + dt_s;

            auto rv_dot = [&](const Spacecraft &sc, const double t_s) -> double {
                const Vec3 r_rel = sc.state.position_m - eph.body_position_at(pi, t_s);
                const Vec3 v_rel = sc.state.velocity_mps - eph.body_velocity_at(pi, t_s);
                return glm::dot(r_rel, v_rel);
            };

            const Spacecraft sc1 = propagate_sc(sc0, t0_s, dt_s);
            if (out_sc1)
                *out_sc1 = sc1;

            const double f0 = rv_dot(sc0, t0_s);
            const double f1 = rv_dot(sc1, t1_s);
            if (!std::isfinite(f0) || !std::isfinite(f1))
                return std::nullopt;

            if (!((f0 <= 0.0 && f1 >= 0.0) || (f0 >= 0.0 && f1 <= 0.0)))
                return std::nullopt;

            const double time_tol = std::max(0.0, opt.time_tol_s);
            double a = t0_s, b = t1_s, fa = f0;

            for (int it = 0; it < opt.max_bisect_iters && std::abs(b - a) > time_tol; ++it)
            {
                const double m = 0.5 * (a + b);
                const double fm = rv_dot(propagate_sc(sc0, t0_s, m - t0_s), m);
                if (!std::isfinite(fm))
                    break;

                if ((fa <= 0.0 && fm >= 0.0) || (fa >= 0.0 && fm <= 0.0))
                    b = m;
                else
                {
                    a = m;
                    fa = fm;
                }
            }

            const ApsisKind kind = (f0 <= 0.0 && f1 >= 0.0) ? ApsisKind::Periapsis : ApsisKind::Apoapsis;
            return ApsisEvent{
                .t_event_s = 0.5 * (a + b),
                .primary_body_id = primary_body_id,
                .spacecraft_id = sc0.id,
                .kind = kind,
            };
        }

        inline std::vector<ApsisEvent> predict_apsis_events_from_ephemeris_by_id_(
            const GameSimulation &sim, const CelestialEphemeris &eph,
            const SpacecraftId spacecraft_id, const BodyId primary_body_id,
            const TrajectoryOptions &traj_opt, const EventOptions &opt)
        {
            if (!sim.has_body(primary_body_id))
                return {};
            return scan_ephemeris_for_events_<ApsisEvent>(
                sim, eph, spacecraft_id, traj_opt, opt,
                [&](const auto &seg, const Spacecraft &sc, double t, double remaining,
                    auto propagate_sc, Spacecraft *sc1) {
                    return find_earliest_apsis_event_in_interval_(
                        sim.massive_bodies(), seg, sc, t, remaining, primary_body_id,
                        opt, propagate_sc, sc1);
                },
                [](const ApsisEvent &) { return false; });
        }

        // ---- Plane node prediction (ascending / descending nodes) ----

        inline std::vector<NodeEvent> predict_plane_nodes_from_ephemeris_by_id_(
            const GameSimulation &sim, const CelestialEphemeris &eph,
            const SpacecraftId spacecraft_id, const BodyId primary_body_id,
            const Vec3 &plane_normal_unit_i, const SpacecraftId target_spacecraft_id,
            const TrajectoryOptions &traj_opt, const EventOptions &opt)
        {
            if (!sim.body_by_id(primary_body_id))
                return {};
            return scan_ephemeris_for_events_<NodeEvent>(
                sim, eph, spacecraft_id, traj_opt, opt,
                [&](const auto &seg, const Spacecraft &sc, double t, double remaining,
                    auto propagate_sc, Spacecraft *sc1) {
                    return find_earliest_plane_node_in_interval_(
                        sim.massive_bodies(), seg, sc, t, remaining, primary_body_id,
                        plane_normal_unit_i, target_spacecraft_id, opt, propagate_sc, sc1);
                },
                [](const NodeEvent &) { return false; });
        }

    } // namespace detail

    // ---- Public API helpers ----

    namespace detail
    {
        inline const EventOptions &resolve_event_opt_(
            const GameSimulation &sim, const std::optional<EventOptions> &opt)
        {
            return opt ? *opt : sim.config().events;
        }
    } // namespace detail

    // ---- Public API: Spacecraft Events ----

    inline std::vector<Event> predict_spacecraft_events(
        const GameSimulation &sim, const CelestialEphemeris &eph,
        const SpacecraftId spacecraft_id,
        const TrajectoryOptions &opt = {},
        std::optional<EventOptions> event_opt = std::nullopt)
    {
        return detail::predict_spacecraft_events_from_ephemeris_by_id_(
            sim, eph, spacecraft_id, opt, detail::resolve_event_opt_(sim, event_opt));
    }

    inline std::vector<Event> predict_spacecraft_events(
        const GameSimulation &sim, const SpacecraftId spacecraft_id,
        const TrajectoryOptions &opt = {},
        std::optional<EventOptions> event_opt = std::nullopt)
    {
        return predict_spacecraft_events(
            sim, build_celestial_ephemeris(sim, opt), spacecraft_id, opt, event_opt);
    }

    // ---- Public API: Apsis Events ----

    inline std::vector<ApsisEvent> predict_apsis_events(
        const GameSimulation &sim, const CelestialEphemeris &eph,
        const SpacecraftId spacecraft_id, const BodyId primary_body_id,
        const TrajectoryOptions &opt = {},
        std::optional<EventOptions> event_opt = std::nullopt)
    {
        return detail::predict_apsis_events_from_ephemeris_by_id_(
            sim, eph, spacecraft_id, primary_body_id, opt,
            detail::resolve_event_opt_(sim, event_opt));
    }

    inline std::vector<ApsisEvent> predict_apsis_events(
        const GameSimulation &sim, const SpacecraftId spacecraft_id,
        const BodyId primary_body_id,
        const TrajectoryOptions &opt = {},
        std::optional<EventOptions> event_opt = std::nullopt)
    {
        return predict_apsis_events(
            sim, build_celestial_ephemeris(sim, opt), spacecraft_id, primary_body_id, opt, event_opt);
    }

    // ---- Public API: Equatorial Nodes ----

    inline std::vector<NodeEvent> predict_equatorial_nodes(
        const GameSimulation &sim, const CelestialEphemeris &eph,
        const SpacecraftId spacecraft_id, const BodyId primary_body_id,
        const TrajectoryOptions &opt = {},
        std::optional<EventOptions> node_opt = std::nullopt)
    {
        const MassiveBody *primary = sim.body_by_id(primary_body_id);
        if (!primary)
            return {};
        const auto axis = detail::spin_axis_unit_(*primary);
        if (!axis)
            return {};
        return detail::predict_plane_nodes_from_ephemeris_by_id_(
            sim, eph, spacecraft_id, primary_body_id, *axis, kInvalidSpacecraftId,
            opt, detail::resolve_event_opt_(sim, node_opt));
    }

    inline std::vector<NodeEvent> predict_equatorial_nodes(
        const GameSimulation &sim, const SpacecraftId spacecraft_id,
        const BodyId primary_body_id,
        const TrajectoryOptions &opt = {},
        std::optional<EventOptions> node_opt = std::nullopt)
    {
        return predict_equatorial_nodes(
            sim, build_celestial_ephemeris(sim, opt), spacecraft_id, primary_body_id, opt, node_opt);
    }

    // ---- Public API: Target Plane Nodes ----

    inline std::vector<NodeEvent> predict_target_plane_nodes(
        const GameSimulation &sim, const CelestialEphemeris &eph,
        const SpacecraftId spacecraft_id, const SpacecraftId target_spacecraft_id,
        const BodyId primary_body_id,
        const TrajectoryOptions &opt = {},
        std::optional<EventOptions> node_opt = std::nullopt)
    {
        const Spacecraft *target_ptr = sim.spacecraft_by_id(target_spacecraft_id);
        const MassiveBody *primary_ptr = sim.body_by_id(primary_body_id);
        if (!target_ptr || !primary_ptr)
            return {};

        const double t0 = sim.time_s();
        const auto &resolved_opt = detail::resolve_event_opt_(sim, node_opt);
        const State ps = eph.empty() ? primary_ptr->state : eph.body_state_at_by_id(primary_body_id, t0);
        const auto n = detail::orbit_plane_normal_unit_(
            target_ptr->state.position_m - ps.position_m,
            target_ptr->state.velocity_mps - ps.velocity_mps);
        if (!n)
            return {};

        return detail::predict_plane_nodes_from_ephemeris_by_id_(
            sim, eph, spacecraft_id, primary_body_id, *n, target_spacecraft_id, opt, resolved_opt);
    }

    inline std::vector<NodeEvent> predict_target_plane_nodes(
        const GameSimulation &sim, const SpacecraftId spacecraft_id,
        const SpacecraftId target_spacecraft_id, const BodyId primary_body_id,
        const TrajectoryOptions &opt = {},
        std::optional<EventOptions> node_opt = std::nullopt)
    {
        return predict_target_plane_nodes(
            sim, build_celestial_ephemeris(sim, opt), spacecraft_id, target_spacecraft_id,
            primary_body_id, opt, node_opt);
    }

    // ---- Public API: Maneuver Events ----

    inline std::vector<ManeuverEvent> predict_maneuver_events(
        const GameSimulation &sim, const SpacecraftId spacecraft_id, const TrajectoryOptions &opt = {})
    {
        std::vector<ManeuverEvent> out;
        if (opt.duration_s <= 0.0 || !std::isfinite(opt.duration_s))
            return out;

        const double t0 = sim.time_s();
        const double t_end = t0 + opt.duration_s;
        const ManeuverPlan &plan = sim.maneuver_plan();

        for (std::size_t i = 0; i < plan.segments.size(); ++i)
        {
            const BurnSegment &seg = plan.segments[i];
            if (!segment_applies_to_spacecraft(seg, spacecraft_id))
                continue;

            if (seg.t_start_s >= t0 && seg.t_start_s <= t_end)
                out.push_back({.t_event_s = seg.t_start_s, .type = ManeuverEventType::BurnStart,
                               .spacecraft_id = spacecraft_id, .index = i});
            if (seg.t_end_s >= t0 && seg.t_end_s <= t_end)
                out.push_back({.t_event_s = seg.t_end_s, .type = ManeuverEventType::BurnEnd,
                               .spacecraft_id = spacecraft_id, .index = i});
        }

        for (std::size_t i = 0; i < plan.impulses.size(); ++i)
        {
            const ImpulseSegment &imp = plan.impulses[i];
            if (!segment_applies_to_spacecraft(imp, spacecraft_id))
                continue;
            if (imp.t_s >= t0 && imp.t_s <= t_end)
                out.push_back({.t_event_s = imp.t_s, .type = ManeuverEventType::Impulse,
                               .spacecraft_id = spacecraft_id, .index = i});
        }

        auto type_rank = [](ManeuverEventType t) -> int {
            switch (t)
            {
            case ManeuverEventType::BurnStart: return 0;
            case ManeuverEventType::Impulse:   return 1;
            case ManeuverEventType::BurnEnd:   return 2;
            }
            return 3;
        };

        std::stable_sort(out.begin(), out.end(), [&](const ManeuverEvent &a, const ManeuverEvent &b) {
            if (a.t_event_s != b.t_event_s)
                return a.t_event_s < b.t_event_s;
            const int ao = type_rank(a.type), bo = type_rank(b.type);
            return ao != bo ? ao < bo : a.index < b.index;
        });

        return out;
    }

} // namespace orbitsim
