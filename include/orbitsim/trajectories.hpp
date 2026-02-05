#pragma once

#include "orbitsim/detail/spacecraft_propagation.hpp"
#include "orbitsim/events.hpp"
#include "orbitsim/game_sim.hpp"
#include "orbitsim/nodes.hpp"
#include "orbitsim/spacecraft_state_cache.hpp"
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

        /// Spacecraft state-lookup interval used for LVLH/relative-frame features.
        /// If > 0, the library may build/interpolate a target spacecraft trajectory at this interval to provide
        /// time-aware LVLH frames without re-integrating targets for every lookup.
        /// If <= 0, defaults to `spacecraft_sample_dt_s` (or `sample_dt_s`).
        double spacecraft_lookup_dt_s{0.0};

        /// Massive body integration step; if <= 0, uses sample_dt_s.
        double celestial_dt_s{0.0};

        std::size_t max_samples{2048};  ///< Hard cap to prevent unbounded allocations

        bool include_start{true};   ///< Include sample at t0
        bool include_end{true};     ///< Include sample at t0 + duration

        /// @brief If true, stop sampling at the first detected impact (Impact + Enter).
        ///
        /// This is useful for orbit-line rendering to avoid drawing trajectories after a surface collision.
        bool stop_on_impact{false};
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

        inline double compute_spacecraft_lookup_dt_(const TrajectoryOptions &opt)
        {
            if (opt.spacecraft_lookup_dt_s > 0.0 && std::isfinite(opt.spacecraft_lookup_dt_s))
            {
                return opt.spacecraft_lookup_dt_s;
            }
            return compute_spacecraft_sample_dt_(opt);
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

            const double lookup_dt = compute_spacecraft_lookup_dt_(opt);
            SpacecraftStateCache<CelestialEphemeris> sc_cache(
                    sim.massive_bodies(),
                    eph,
                    sim.maneuver_plan(),
                    sim.config().gravitational_constant,
                    sim.config().softening_length_m,
                    sim.config().spacecraft_integrator,
                    t0,
                    t_end,
                    [&sim](const SpacecraftId id) { return sim.spacecraft_by_id(id); },
                    SpacecraftStateCache<CelestialEphemeris>::Options{.lookup_dt_s = lookup_dt});
            const SpacecraftStateLookup sc_lookup = sc_cache.lookup();

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

                const Spacecraft sc_next = propagate_spacecraft_in_ephemeris(
                        sc, sim.massive_bodies(), eph, sim.maneuver_plan(), sim.config().gravitational_constant,
                        sim.config().softening_length_m, sim.config().spacecraft_integrator, t, h, sc_lookup);

                if (opt.stop_on_impact)
                {
                    const double t0_step = t;
                    const double t1_step = t + h;

                    std::optional<Event> best_impact;
                    for (const auto &body: sim.massive_bodies())
                    {
                        double threshold_m = 0.0;
                        if (!boundary_threshold_m(body, EventType::Impact, &threshold_m))
                        {
                            continue;
                        }

                        const Vec3 bp0 = eph.body_position_at_by_id(body.id, t0_step);
                        const Vec3 bp1 = eph.body_position_at_by_id(body.id, t1_step);
                        const double d0 = glm::length(sc.state.position_m - bp0);
                        const double d1 = glm::length(sc_next.state.position_m - bp1);
                        if (!std::isfinite(d0) || !std::isfinite(d1))
                        {
                            continue;
                        }
                        const double f0 = d0 - threshold_m;
                        const double f1 = d1 - threshold_m;
                        if (!((f0 <= 0.0 && f1 >= 0.0) || (f0 >= 0.0 && f1 <= 0.0)))
                        {
                            continue;
                        }

                        const Crossing crossing = (f0 > 0.0 && f1 <= 0.0) ? Crossing::Enter : Crossing::Exit;
                        if (crossing != Crossing::Enter)
                        {
                            continue;
                        }

                        const double t_event = bisect_crossing_time_s(
                                t0_step, t1_step, f0, sim.config().events, [&](const double t_s) -> double {
                                    const Spacecraft sc_at = propagate_spacecraft_in_ephemeris(
                                            sc, sim.massive_bodies(), eph, sim.maneuver_plan(),
                                            sim.config().gravitational_constant, sim.config().softening_length_m,
                                            sim.config().spacecraft_integrator, t0_step, t_s - t0_step, sc_lookup);
                                    const Vec3 bp = eph.body_position_at_by_id(body.id, t_s);
                                    return glm::length(sc_at.state.position_m - bp) - threshold_m;
                                });
                        if (!std::isfinite(t_event))
                        {
                            continue;
                        }

                        const Event e = Event{
                                .type = EventType::Impact,
                                .body_id = body.id,
                                .crossing = Crossing::Enter,
                                .t_event_s = t_event,
                                .spacecraft_id = spacecraft_id,
                        };
                        if (!best_impact.has_value() || e.t_event_s < best_impact->t_event_s)
                        {
                            best_impact = e;
                        }
                    }

                    if (best_impact.has_value())
                    {
                        const double dt_imp = best_impact->t_event_s - t0_step;
                        const Spacecraft sc_imp = propagate_spacecraft_in_ephemeris(
                                sc, sim.massive_bodies(), eph, sim.maneuver_plan(), sim.config().gravitational_constant,
                                sim.config().softening_length_m, sim.config().spacecraft_integrator, t0_step, dt_imp,
                                sc_lookup);

                        t = best_impact->t_event_s;
                        out.push_back(make_sample_(t, sc_imp.state));
                        break;
                    }
                }

                sc = sc_next;
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

            const double lookup_dt = compute_spacecraft_lookup_dt_(traj_opt);
            SpacecraftStateCache<CelestialEphemeris> sc_cache(
                    sim.massive_bodies(),
                    eph,
                    sim.maneuver_plan(),
                    sim.config().gravitational_constant,
                    sim.config().softening_length_m,
                    sim.config().spacecraft_integrator,
                    t0,
                    t_end,
                    [&sim](const SpacecraftId id) { return sim.spacecraft_by_id(id); },
                    SpacecraftStateCache<CelestialEphemeris>::Options{.lookup_dt_s = lookup_dt});
            const SpacecraftStateLookup sc_lookup = sc_cache.lookup();

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
                                                             sim.config().spacecraft_integrator, t0_s, dt_sc_s, sc_lookup);
                };

                while (remaining > 0.0)
                {
                    Spacecraft sc1{};
                    const std::optional<Event> e = find_earliest_event_in_interval(
                            sim.massive_bodies(), seg, sc, t, remaining, sim.maneuver_plan(), event_opt, propagate_sc,
                            &sc1);
                    if (!e.has_value())
                    {
                        sc = sc1;
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

        inline std::optional<std::size_t> body_index_for_id_vec_(const std::vector<MassiveBody> &bodies,
                                                                 const BodyId body_id)
        {
            if (body_id == kInvalidBodyId)
            {
                return std::nullopt;
            }
            for (std::size_t i = 0; i < bodies.size(); ++i)
            {
                if (bodies[i].id == body_id)
                {
                    return i;
                }
            }
            return std::nullopt;
        }

        template<class EphemerisLike, class Propagator>
        inline std::optional<ApsisEvent>
        find_earliest_apsis_event_in_interval_(const std::vector<MassiveBody> &bodies, const EphemerisLike &eph,
                                               const Spacecraft &sc0, const double t0_s, const double dt_s,
                                               const BodyId primary_body_id, const EventOptions &opt,
                                               Propagator propagate_sc, Spacecraft *out_sc1 = nullptr)
        {
            if (!(dt_s > 0.0) || !std::isfinite(dt_s) || !(opt.max_bisect_iters > 0))
            {
                return std::nullopt;
            }

            const std::optional<std::size_t> primary_index_opt = body_index_for_id_vec_(bodies, primary_body_id);
            if (!primary_index_opt.has_value())
            {
                return std::nullopt;
            }
            const std::size_t primary_index = *primary_index_opt;

            const double t1_s = t0_s + dt_s;

            auto rv_dot_m2ps = [&](const Spacecraft &sc, const double t_s) -> double {
                const Vec3 rp = eph.body_position_at(primary_index, t_s);
                const Vec3 vp = eph.body_velocity_at(primary_index, t_s);
                const Vec3 r_rel = sc.state.position_m - rp;
                const Vec3 v_rel = sc.state.velocity_mps - vp;
                return glm::dot(r_rel, v_rel);
            };

            const Spacecraft sc1 = propagate_sc(sc0, t0_s, dt_s);
            if (out_sc1 != nullptr)
            {
                *out_sc1 = sc1;
            }

            const double f0 = rv_dot_m2ps(sc0, t0_s);
            const double f1 = rv_dot_m2ps(sc1, t1_s);
            if (!std::isfinite(f0) || !std::isfinite(f1))
            {
                return std::nullopt;
            }

            const bool sign_change = (f0 <= 0.0 && f1 >= 0.0) || (f0 >= 0.0 && f1 <= 0.0);
            if (!sign_change)
            {
                return std::nullopt;
            }

            // Bisection refine on the sign-change root of dot(r_rel, v_rel). Converge based on time tolerance only.
            const double time_tol_s = std::max(0.0, opt.time_tol_s);
            double a = t0_s;
            double b = t1_s;
            double fa = f0;

            for (int it = 0; it < opt.max_bisect_iters; ++it)
            {
                const double m = 0.5 * (a + b);
                const Spacecraft scm = propagate_sc(sc0, t0_s, m - t0_s);
                const double fm = rv_dot_m2ps(scm, m);
                if (!std::isfinite(fm))
                {
                    break;
                }
                if (std::abs(b - a) <= time_tol_s)
                {
                    a = m;
                    b = m;
                    break;
                }

                const bool left = (fa <= 0.0 && fm >= 0.0) || (fa >= 0.0 && fm <= 0.0);
                if (left)
                {
                    b = m;
                }
                else
                {
                    a = m;
                    fa = fm;
                }
            }

            const double t_event_s = 0.5 * (a + b);
            if (!std::isfinite(t_event_s))
            {
                return std::nullopt;
            }

            const ApsisKind kind = (f0 <= 0.0 && f1 >= 0.0) ? ApsisKind::Periapsis : ApsisKind::Apoapsis;
            return ApsisEvent{.t_event_s = t_event_s,
                              .primary_body_id = primary_body_id,
                              .spacecraft_id = sc0.id,
                              .kind = kind};
        }

        inline std::vector<ApsisEvent> predict_apsis_events_from_ephemeris_by_id_(
                const GameSimulation &sim, const CelestialEphemeris &eph, const SpacecraftId spacecraft_id,
                const BodyId primary_body_id, const TrajectoryOptions &traj_opt, const EventOptions &opt)
        {
            std::vector<ApsisEvent> out;
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
            if (!sim.has_body(primary_body_id))
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

            const double lookup_dt = compute_spacecraft_lookup_dt_(traj_opt);
            SpacecraftStateCache<CelestialEphemeris> sc_cache(
                    sim.massive_bodies(),
                    eph,
                    sim.maneuver_plan(),
                    sim.config().gravitational_constant,
                    sim.config().softening_length_m,
                    sim.config().spacecraft_integrator,
                    t0,
                    t_end,
                    [&sim](const SpacecraftId id) { return sim.spacecraft_by_id(id); },
                    SpacecraftStateCache<CelestialEphemeris>::Options{.lookup_dt_s = lookup_dt});
            const SpacecraftStateLookup sc_lookup = sc_cache.lookup();

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
                                                             sim.config().spacecraft_integrator, t0_s, dt_sc_s, sc_lookup);
                };

                while (remaining > 0.0)
                {
                    Spacecraft sc1{};
                    const std::optional<ApsisEvent> e = find_earliest_apsis_event_in_interval_(
                            sim.massive_bodies(), seg, sc, t, remaining, primary_body_id, opt, propagate_sc, &sc1);
                    if (!e.has_value())
                    {
                        sc = sc1;
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

                    ApsisEvent logged = *e;
                    logged.t_event_s = t;
                    out.push_back(logged);

                    // Nudge forward to avoid repeated detection at the same time.
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

            const double lookup_dt = compute_spacecraft_lookup_dt_(traj_opt);
            SpacecraftStateCache<CelestialEphemeris> sc_cache(
                    sim.massive_bodies(),
                    eph,
                    sim.maneuver_plan(),
                    sim.config().gravitational_constant,
                    sim.config().softening_length_m,
                    sim.config().spacecraft_integrator,
                    t0,
                    t_end,
                    [&sim](const SpacecraftId id) { return sim.spacecraft_by_id(id); },
                    SpacecraftStateCache<CelestialEphemeris>::Options{.lookup_dt_s = lookup_dt});
            const SpacecraftStateLookup sc_lookup_node = sc_cache.lookup();

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
                                                             sim.config().spacecraft_integrator, t0_s, dt_sc_s, sc_lookup_node);
                };

                while (remaining > 0.0)
                {
                    Spacecraft sc1{};
                    const std::optional<NodeEvent> e = find_earliest_plane_node_in_interval_(
                            sim.massive_bodies(), seg, sc, t, remaining, primary_body_id, plane_normal_unit_i,
                            target_spacecraft_id, opt, propagate_sc, &sc1);

                    if (!e.has_value())
                    {
                        sc = sc1;
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
     * @brief Predict periapsis/apoapsis events for a spacecraft relative to a primary body.
     *
     * Uses the general N-body trajectory (with maneuvers) and detects local extrema of distance to the primary
     * by finding roots of dot(r_rel, v_rel).
     */
    inline std::vector<ApsisEvent> predict_apsis_events(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                        const SpacecraftId spacecraft_id, const BodyId primary_body_id,
                                                        const TrajectoryOptions &opt, const EventOptions &event_opt)
    {
        return detail::predict_apsis_events_from_ephemeris_by_id_(sim, eph, spacecraft_id, primary_body_id, opt,
                                                                 event_opt);
    }

    inline std::vector<ApsisEvent> predict_apsis_events(const GameSimulation &sim, const CelestialEphemeris &eph,
                                                        const SpacecraftId spacecraft_id, const BodyId primary_body_id,
                                                        const TrajectoryOptions &opt = {})
    {
        return predict_apsis_events(sim, eph, spacecraft_id, primary_body_id, opt, sim.config().events);
    }

    inline std::vector<ApsisEvent> predict_apsis_events(const GameSimulation &sim, const SpacecraftId spacecraft_id,
                                                        const BodyId primary_body_id, const TrajectoryOptions &opt,
                                                        const EventOptions &event_opt)
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return predict_apsis_events(sim, eph, spacecraft_id, primary_body_id, opt, event_opt);
    }

    inline std::vector<ApsisEvent> predict_apsis_events(const GameSimulation &sim, const SpacecraftId spacecraft_id,
                                                        const BodyId primary_body_id, const TrajectoryOptions &opt = {})
    {
        const CelestialEphemeris eph = detail::build_celestial_ephemeris_(sim, opt);
        return predict_apsis_events(sim, eph, spacecraft_id, primary_body_id, opt, sim.config().events);
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

    /**
     * @brief Predict maneuver timeline events (burn start/end, impulse time) for a spacecraft.
     *
     * This does not propagate the trajectory; it simply filters the simulation's ManeuverPlan to the interval
     * `[sim.time_s(), sim.time_s() + opt.duration_s]` and returns the resulting boundary events in time order.
     */
    inline std::vector<ManeuverEvent> predict_maneuver_events(const GameSimulation &sim, const SpacecraftId spacecraft_id,
                                                              const TrajectoryOptions &opt = {})
    {
        std::vector<ManeuverEvent> out;
        if (!(opt.duration_s > 0.0) || !std::isfinite(opt.duration_s))
        {
            return out;
        }

        const double t0 = sim.time_s();
        const double t_end = t0 + opt.duration_s;
        if (!std::isfinite(t_end) || !(t_end > t0))
        {
            return out;
        }

        const ManeuverPlan &plan = sim.maneuver_plan();

        for (std::size_t i = 0; i < plan.segments.size(); ++i)
        {
            const BurnSegment &seg = plan.segments[i];
            if (!segment_applies_to_spacecraft(seg, spacecraft_id))
            {
                continue;
            }

            if (seg.t_start_s >= t0 && seg.t_start_s <= t_end && std::isfinite(seg.t_start_s))
            {
                out.push_back(ManeuverEvent{.t_event_s = seg.t_start_s,
                                            .type = ManeuverEventType::BurnStart,
                                            .spacecraft_id = spacecraft_id,
                                            .index = i});
            }
            if (seg.t_end_s >= t0 && seg.t_end_s <= t_end && std::isfinite(seg.t_end_s))
            {
                out.push_back(ManeuverEvent{.t_event_s = seg.t_end_s,
                                            .type = ManeuverEventType::BurnEnd,
                                            .spacecraft_id = spacecraft_id,
                                            .index = i});
            }
        }

        for (std::size_t i = 0; i < plan.impulses.size(); ++i)
        {
            const ImpulseSegment &imp = plan.impulses[i];
            if (!segment_applies_to_spacecraft(imp, spacecraft_id))
            {
                continue;
            }
            if (imp.t_s >= t0 && imp.t_s <= t_end && std::isfinite(imp.t_s))
            {
                out.push_back(ManeuverEvent{.t_event_s = imp.t_s,
                                            .type = ManeuverEventType::Impulse,
                                            .spacecraft_id = spacecraft_id,
                                            .index = i});
            }
        }

        const auto type_order = [](const ManeuverEventType type) -> int {
            switch (type)
            {
            case ManeuverEventType::BurnStart:
                return 0;
            case ManeuverEventType::Impulse:
                return 1;
            case ManeuverEventType::BurnEnd:
                return 2;
            }
            return 3;
        };

        std::stable_sort(out.begin(), out.end(), [&](const ManeuverEvent &a, const ManeuverEvent &b) {
            if (a.t_event_s != b.t_event_s)
            {
                return a.t_event_s < b.t_event_s;
            }
            const int ao = type_order(a.type);
            const int bo = type_order(b.type);
            if (ao != bo)
            {
                return ao < bo;
            }
            return a.index < b.index;
        });

        return out;
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

        /// @brief Set the spacecraft lookup time step [s] used for LVLH/relative-frame features.
        TrajectoryOptionsBuilder &spacecraft_lookup_dt(const double spacecraft_lookup_dt_s)
        {
            opt_.spacecraft_lookup_dt_s = spacecraft_lookup_dt_s;
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

        /// @brief Set whether to stop sampling at the first impact event.
        TrajectoryOptionsBuilder &stop_on_impact(const bool stop)
        {
            opt_.stop_on_impact = stop;
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
