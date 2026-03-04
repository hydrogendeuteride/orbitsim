#pragma once

#include "orbitsim/detail/spacecraft_propagation.hpp"
#include "orbitsim/events.hpp"
#include "orbitsim/game_sim.hpp"
#include "orbitsim/spacecraft_state_cache.hpp"
#include "orbitsim/time_utils.hpp"
#include "orbitsim/trajectory_segments.hpp"
#include "orbitsim/trajectory_types.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <optional>
#include <vector>

namespace orbitsim
{

    struct TrajectoryOptions
    {
        double duration_s{3600.0};  ///< Prediction duration from sim.time_s()
        double sample_dt_s{10.0};   ///< Sample interval; if <= 0, derived from duration/(max_samples-1)

        /// Spacecraft sample interval; if > 0, overrides sample_dt_s for spacecraft.
        double spacecraft_sample_dt_s{0.0};

        /// Spacecraft state-lookup interval for LVLH/relative-frame features.
        /// If <= 0, defaults to spacecraft_sample_dt_s (or sample_dt_s).
        double spacecraft_lookup_dt_s{0.0};

        /// Massive body integration step; if <= 0, uses sample_dt_s.
        double celestial_dt_s{0.0};

        std::size_t max_samples{2048};  ///< Hard cap to prevent unbounded allocations
        bool include_start{true};       ///< Include sample at t0
        bool include_end{true};         ///< Include sample at t0 + duration
        bool stop_on_impact{false};     ///< Stop sampling at first detected impact
    };

    namespace detail
    {
        inline double compute_sample_dt_(const TrajectoryOptions &opt)
        {
            if (opt.sample_dt_s > 0.0)
                return opt.sample_dt_s;
            if (opt.max_samples <= 1 || !(opt.duration_s > 0.0))
                return 0.0;
            return opt.duration_s / static_cast<double>(opt.max_samples - 1);
        }

        inline double compute_celestial_dt_(const TrajectoryOptions &opt)
        {
            if (opt.celestial_dt_s > 0.0)
                return opt.celestial_dt_s;
            return compute_sample_dt_(opt);
        }

        inline double compute_spacecraft_sample_dt_(const TrajectoryOptions &opt)
        {
            if (opt.spacecraft_sample_dt_s > 0.0)
                return opt.spacecraft_sample_dt_s;
            return compute_sample_dt_(opt);
        }

        inline double compute_spacecraft_lookup_dt_(const TrajectoryOptions &opt)
        {
            if (opt.spacecraft_lookup_dt_s > 0.0)
                return opt.spacecraft_lookup_dt_s;
            return compute_spacecraft_sample_dt_(opt);
        }

        inline double compute_segment_dt_(const TrajectorySegmentOptions &opt)
        {
            if (opt.max_segments == 0)
                return 0.0;
            if (opt.lookup_dt_s > 0.0)
                return opt.lookup_dt_s;
            if (!(opt.duration_s > 0.0))
                return 0.0;
            return opt.duration_s / static_cast<double>(opt.max_segments);
        }

        inline CelestialEphemeris build_celestial_ephemeris_(const GameSimulation &sim, const TrajectoryOptions &opt)
        {
            CelestialEphemeris eph;

            const double dt = compute_celestial_dt_(opt);
            if (!(dt > 0.0) || !(opt.duration_s > 0.0) || opt.max_samples == 0)
                return eph;

            std::vector<BodyId> body_ids;
            body_ids.reserve(sim.massive_bodies().size());
            for (const auto &b : sim.massive_bodies())
                body_ids.push_back(b.id);
            eph.set_body_ids(std::move(body_ids));

            std::vector<MassiveBody> massive = sim.massive_bodies();
            double t = sim.time_s();
            const double t_end = t + opt.duration_s;

            std::vector<State> start_states;
            std::vector<State> end_states;

            while (t < t_end && eph.segments.size() < opt.max_samples)
            {
                const double h = std::min(dt, t_end - t);
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

        inline std::vector<TrajectorySample> predict_body_trajectory_from_ephemeris_(
            const GameSimulation &sim, const CelestialEphemeris &eph,
            const BodyId body_id, const TrajectoryOptions &opt)
        {
            std::vector<TrajectorySample> out;

            std::size_t body_index = 0;
            if (eph.empty() || !eph.body_index_for_id(body_id, &body_index))
                return out;

            const double dt = compute_sample_dt_(opt);
            if (!(dt > 0.0))
                return out;

            const double t0 = sim.time_s();
            const double t_end = t0 + opt.duration_s;
            out.reserve(std::min<std::size_t>(opt.max_samples, 1024));

            if (opt.include_start)
                out.push_back(make_sample_(t0, eph.body_state_at(body_index, t0)));

            double t = t0;
            while (t < t_end && out.size() < opt.max_samples)
            {
                t += std::min(dt, t_end - t);
                if (!opt.include_end && t >= t_end)
                    break;
                out.push_back(make_sample_(t, eph.body_state_at(body_index, t)));
            }

            return out;
        }

        /// Find the earliest impact (Enter crossing) time within a spacecraft propagation step.
        inline std::optional<double> find_impact_time_(
            const Spacecraft &sc0, const Spacecraft &sc1,
            const std::vector<MassiveBody> &bodies, const CelestialEphemeris &eph,
            const ManeuverPlan &plan, const double G, const double softening,
            const DOPRI5Options &integrator, const double t0, const double h,
            const SpacecraftStateLookup &sc_lookup, const EventOptions &event_opt)
        {
            std::optional<double> best_t;

            for (const auto &body : bodies)
            {
                double threshold_m = 0.0;
                if (!boundary_threshold_m(body, EventType::Impact, &threshold_m))
                    continue;

                const Vec3 bp0 = eph.body_position_at_by_id(body.id, t0);
                const Vec3 bp1 = eph.body_position_at_by_id(body.id, t0 + h);
                const double f0 = glm::length(sc0.state.position_m - bp0) - threshold_m;
                const double f1 = glm::length(sc1.state.position_m - bp1) - threshold_m;

                // Only detect Enter crossings (outside → inside surface)
                if (!(f0 > 0.0 && f1 <= 0.0))
                    continue;

                const double t_event = bisect_crossing_time_s(
                    t0, t0 + h, f0, event_opt, [&](const double t_s) -> double {
                        const Spacecraft sc_at = propagate_spacecraft_in_ephemeris(
                            sc0, bodies, eph, plan, G, softening, integrator, t0, t_s - t0, sc_lookup);
                        return glm::length(sc_at.state.position_m - eph.body_position_at_by_id(body.id, t_s))
                               - threshold_m;
                    });

                if (std::isfinite(t_event) && (!best_t || t_event < *best_t))
                    best_t = t_event;
            }

            return best_t;
        }

        inline std::vector<TrajectorySegment> predict_spacecraft_trajectory_segments_from_ephemeris_by_id_(
            const GameSimulation &sim, const CelestialEphemeris &eph,
            const SpacecraftId spacecraft_id, const TrajectorySegmentOptions &opt)
        {
            std::vector<TrajectorySegment> out;

            const Spacecraft *sc_ptr = sim.spacecraft_by_id(spacecraft_id);
            if (!sc_ptr)
                return out;

            const double dt = compute_segment_dt_(opt);
            if (!(dt > 0.0))
                return out;

            const double t0 = sim.time_s();
            const double t_end = t0 + opt.duration_s;
            out.reserve(std::min<std::size_t>(opt.max_segments, 1024));

            Spacecraft sc = *sc_ptr;
            double t = t0;
            bool first_segment = true;

            const double sc_lookup_dt = (opt.lookup_dt_s > 0.0) ? opt.lookup_dt_s : dt;
            SpacecraftStateCache<CelestialEphemeris> sc_cache(
                sim.massive_bodies(), eph, sim.maneuver_plan(),
                sim.config().gravitational_constant, sim.config().softening_length_m,
                sim.config().spacecraft_integrator, t0, t_end,
                [&sim](const SpacecraftId id) { return sim.spacecraft_by_id(id); },
                SpacecraftStateCache<CelestialEphemeris>::Options{
                    .lookup_dt_s = sc_lookup_dt,
                    .max_segments = std::max<std::size_t>(opt.max_segments, 1),
                });
            const SpacecraftStateLookup sc_lookup = sc_cache.lookup();

            while (t < t_end && out.size() < opt.max_segments)
            {
                const double h = std::min(dt, t_end - t);
                if (!(h > 0.0) || !std::isfinite(h))
                    break;

                const Spacecraft sc_next = propagate_spacecraft_in_ephemeris(
                    sc, sim.massive_bodies(), eph, sim.maneuver_plan(),
                    sim.config().gravitational_constant, sim.config().softening_length_m,
                    sim.config().spacecraft_integrator, t, h, sc_lookup);

                const bool reaches_end = (t + h >= t_end);
                const bool emit_segment =
                    !(first_segment && !opt.include_start) &&
                    !(reaches_end && !opt.include_end);

                if (opt.stop_on_impact)
                {
                    const auto impact_t = find_impact_time_(
                        sc, sc_next, sim.massive_bodies(), eph, sim.maneuver_plan(),
                        sim.config().gravitational_constant, sim.config().softening_length_m,
                        sim.config().spacecraft_integrator, t, h, sc_lookup, sim.config().events);

                    if (impact_t)
                    {
                        const double h_imp = *impact_t - t;
                        if (emit_segment && h_imp > 0.0 && std::isfinite(h_imp) && out.size() < opt.max_segments)
                        {
                            const Spacecraft sc_imp = propagate_spacecraft_in_ephemeris(
                                sc, sim.massive_bodies(), eph, sim.maneuver_plan(),
                                sim.config().gravitational_constant, sim.config().softening_length_m,
                                sim.config().spacecraft_integrator, t, h_imp, sc_lookup);

                            out.push_back(TrajectorySegment{
                                .t0_s = t,
                                .dt_s = h_imp,
                                .start = sc.state,
                                .end = sc_imp.state,
                                .flags = 0u,
                            });
                        }
                        break;
                    }
                }

                if (emit_segment)
                {
                    out.push_back(TrajectorySegment{
                        .t0_s = t,
                        .dt_s = h,
                        .start = sc.state,
                        .end = sc_next.state,
                        .flags = 0u,
                    });
                }

                sc = sc_next;
                t += h;
                first_segment = false;
            }

            return out;
        }

        inline State state_at_segment_time_(const TrajectorySegment &seg, const double t_s)
        {
            if (!(seg.dt_s > 0.0) || !std::isfinite(seg.dt_s))
                return seg.start;

            const double tau = clamp01((t_s - seg.t0_s) / seg.dt_s);
            const double t_eval_s = seg.t0_s + tau * seg.dt_s;

            State out{};
            out.position_m = hermite_position(seg.start.position_m, seg.start.velocity_mps,
                                              seg.end.position_m, seg.end.velocity_mps, seg.dt_s, tau);
            out.velocity_mps = hermite_velocity_mps(seg.start.position_m, seg.start.velocity_mps,
                                                    seg.end.position_m, seg.end.velocity_mps, seg.dt_s, tau);
            out.spin.axis = normalized_or(seg.start.spin.axis, Vec3{0.0, 1.0, 0.0});
            out.spin.rate_rad_per_s = seg.start.spin.rate_rad_per_s;
            out.spin.angle_rad = seg.start.spin.angle_rad + seg.start.spin.rate_rad_per_s * (t_eval_s - seg.t0_s);
            return out;
        }

        inline State state_at_segments_time_(const std::vector<TrajectorySegment> &segments, const double t_s)
        {
            if (segments.empty())
                return {};

            if (!(t_s > segments.front().t0_s))
                return segments.front().start;

            const auto it = std::upper_bound(
                segments.begin(),
                segments.end(),
                t_s,
                [](const double t, const TrajectorySegment &seg) { return t < seg.t0_s; });

            if (it == segments.begin())
                return state_at_segment_time_(segments.front(), t_s);
            if (it == segments.end())
                return state_at_segment_time_(segments.back(), t_s);
            return state_at_segment_time_(*std::prev(it), t_s);
        }

        inline std::vector<TrajectorySample> sample_trajectory_segments_uniform_dt_(
            const std::vector<TrajectorySegment> &segments, const double sample_dt_s,
            const std::size_t max_samples, const bool include_start, const bool include_end)
        {
            std::vector<TrajectorySample> out;
            if (segments.empty() || !(sample_dt_s > 0.0) || max_samples == 0)
                return out;

            const double t0 = segments.front().t0_s;
            const TrajectorySegment &last = segments.back();
            const double t_end = last.t0_s + last.dt_s;
            if (!(t_end >= t0))
                return out;

            out.reserve(std::min<std::size_t>(max_samples, 1024));

            if (include_start && out.size() < max_samples)
                out.push_back(make_sample_(t0, segments.front().start));

            double t = t0;
            while (t < t_end && out.size() < max_samples)
            {
                t += std::min(sample_dt_s, t_end - t);
                if (!include_end && t >= t_end)
                    break;
                out.push_back(make_sample_(t, state_at_segments_time_(segments, t)));
            }

            return out;
        }

        inline std::vector<TrajectorySample> predict_spacecraft_trajectory_from_ephemeris_by_id_(
            const GameSimulation &sim, const CelestialEphemeris &eph,
            const SpacecraftId spacecraft_id, const TrajectoryOptions &opt)
        {
            std::vector<TrajectorySample> out;

            const Spacecraft *sc_ptr = sim.spacecraft_by_id(spacecraft_id);
            if (!sc_ptr)
                return out;

            const double dt = compute_spacecraft_sample_dt_(opt);
            if (!(dt > 0.0))
                return out;

            const double t0 = sim.time_s();
            const double t_end = t0 + opt.duration_s;
            out.reserve(std::min<std::size_t>(opt.max_samples, 1024));

            Spacecraft sc = *sc_ptr;
            double t = t0;

            SpacecraftStateCache<CelestialEphemeris> sc_cache(
                sim.massive_bodies(), eph, sim.maneuver_plan(),
                sim.config().gravitational_constant, sim.config().softening_length_m,
                sim.config().spacecraft_integrator, t0, t_end,
                [&sim](const SpacecraftId id) { return sim.spacecraft_by_id(id); },
                SpacecraftStateCache<CelestialEphemeris>::Options{.lookup_dt_s = compute_spacecraft_lookup_dt_(opt)});
            const SpacecraftStateLookup sc_lookup = sc_cache.lookup();

            if (opt.include_start)
                out.push_back(make_sample_(t, sc.state));

            while (t < t_end && out.size() < opt.max_samples)
            {
                const double h = std::min(dt, t_end - t);
                const Spacecraft sc_next = propagate_spacecraft_in_ephemeris(
                    sc, sim.massive_bodies(), eph, sim.maneuver_plan(),
                    sim.config().gravitational_constant, sim.config().softening_length_m,
                    sim.config().spacecraft_integrator, t, h, sc_lookup);

                if (opt.stop_on_impact)
                {
                    const auto impact_t = find_impact_time_(
                        sc, sc_next, sim.massive_bodies(), eph, sim.maneuver_plan(),
                        sim.config().gravitational_constant, sim.config().softening_length_m,
                        sim.config().spacecraft_integrator, t, h, sc_lookup, sim.config().events);

                    if (impact_t)
                    {
                        const Spacecraft sc_imp = propagate_spacecraft_in_ephemeris(
                            sc, sim.massive_bodies(), eph, sim.maneuver_plan(),
                            sim.config().gravitational_constant, sim.config().softening_length_m,
                            sim.config().spacecraft_integrator, t, *impact_t - t, sc_lookup);
                        out.push_back(make_sample_(*impact_t, sc_imp.state));
                        break;
                    }
                }

                sc = sc_next;
                t += h;

                if (!opt.include_end && t >= t_end)
                    break;

                out.push_back(make_sample_(t, sc.state));
            }

            return out;
        }

    } // namespace detail

    // ---- Public API ----

    inline CelestialEphemeris build_celestial_ephemeris(const GameSimulation &sim, const TrajectoryOptions &opt = {})
    {
        return detail::build_celestial_ephemeris_(sim, opt);
    }

    inline std::vector<TrajectorySample> predict_body_trajectory(
        const GameSimulation &sim, const BodyId body_id, const TrajectoryOptions &opt = {})
    {
        return detail::predict_body_trajectory_from_ephemeris_(sim, detail::build_celestial_ephemeris_(sim, opt), body_id, opt);
    }

    inline std::vector<TrajectorySample> predict_body_trajectory(
        const GameSimulation &sim, const CelestialEphemeris &eph,
        const BodyId body_id, const TrajectoryOptions &opt = {})
    {
        return detail::predict_body_trajectory_from_ephemeris_(sim, eph, body_id, opt);
    }

    inline std::vector<TrajectorySample> predict_spacecraft_trajectory(
        const GameSimulation &sim, const SpacecraftId spacecraft_id, const TrajectoryOptions &opt = {})
    {
        return detail::predict_spacecraft_trajectory_from_ephemeris_by_id_(
            sim, detail::build_celestial_ephemeris_(sim, opt), spacecraft_id, opt);
    }

    inline std::vector<TrajectorySample> predict_spacecraft_trajectory(
        const GameSimulation &sim, const CelestialEphemeris &eph,
        const SpacecraftId spacecraft_id, const TrajectoryOptions &opt = {})
    {
        return detail::predict_spacecraft_trajectory_from_ephemeris_by_id_(sim, eph, spacecraft_id, opt);
    }

    inline std::vector<TrajectorySegment> predict_spacecraft_trajectory_segments(
        const GameSimulation &sim, const SpacecraftId spacecraft_id,
        const TrajectorySegmentOptions &opt = {})
    {
        const double dt = detail::compute_segment_dt_(opt);

        TrajectoryOptions eph_opt{};
        eph_opt.duration_s = opt.duration_s;
        eph_opt.sample_dt_s = dt;
        eph_opt.spacecraft_sample_dt_s = dt;
        eph_opt.spacecraft_lookup_dt_s = (opt.lookup_dt_s > 0.0) ? opt.lookup_dt_s : dt;
        eph_opt.celestial_dt_s = dt;
        eph_opt.max_samples = std::max<std::size_t>(opt.max_segments, 1);
        eph_opt.include_start = true;
        eph_opt.include_end = true;
        eph_opt.stop_on_impact = false;

        return detail::predict_spacecraft_trajectory_segments_from_ephemeris_by_id_(
            sim, detail::build_celestial_ephemeris_(sim, eph_opt), spacecraft_id, opt);
    }

    inline std::vector<TrajectorySegment> predict_spacecraft_trajectory_segments(
        const GameSimulation &sim, const CelestialEphemeris &eph,
        const SpacecraftId spacecraft_id, const TrajectorySegmentOptions &opt = {})
    {
        return detail::predict_spacecraft_trajectory_segments_from_ephemeris_by_id_(sim, eph, spacecraft_id, opt);
    }

    inline std::vector<TrajectorySample> sample_trajectory_segments_uniform_dt(
        const std::vector<TrajectorySegment> &segments, const double sample_dt_s,
        const std::size_t max_samples = 2048, const bool include_start = true, const bool include_end = true)
    {
        return detail::sample_trajectory_segments_uniform_dt_(
            segments, sample_dt_s, max_samples, include_start, include_end);
    }

    inline std::vector<TrajectorySample> sample_trajectory_segments_uniform_dt(
        const std::vector<TrajectorySegment> &segments, const TrajectoryOptions &opt)
    {
        return detail::sample_trajectory_segments_uniform_dt_(
            segments, detail::compute_spacecraft_sample_dt_(opt), opt.max_samples, opt.include_start, opt.include_end);
    }

    // ---- Fluent Builder ----

    class TrajectorySegmentOptionsBuilder
    {
    public:
        TrajectorySegmentOptionsBuilder() = default;

        TrajectorySegmentOptionsBuilder &duration(const double v) { opt_.duration_s = v; return *this; }
        TrajectorySegmentOptionsBuilder &max_segments(const std::size_t v) { opt_.max_segments = v; return *this; }
        TrajectorySegmentOptionsBuilder &include_start(const bool v) { opt_.include_start = v; return *this; }
        TrajectorySegmentOptionsBuilder &include_end(const bool v) { opt_.include_end = v; return *this; }
        TrajectorySegmentOptionsBuilder &stop_on_impact(const bool v) { opt_.stop_on_impact = v; return *this; }
        TrajectorySegmentOptionsBuilder &lookup_dt(const double v) { opt_.lookup_dt_s = v; return *this; }

        operator TrajectorySegmentOptions() const { return opt_; }
        TrajectorySegmentOptions build() const { return opt_; }

    private:
        TrajectorySegmentOptions opt_{};
    };

    class TrajectoryOptionsBuilder
    {
    public:
        TrajectoryOptionsBuilder() = default;

        TrajectoryOptionsBuilder &duration(const double v) { opt_.duration_s = v; return *this; }
        TrajectoryOptionsBuilder &sample_dt(const double v) { opt_.sample_dt_s = v; return *this; }
        TrajectoryOptionsBuilder &spacecraft_sample_dt(const double v) { opt_.spacecraft_sample_dt_s = v; return *this; }
        TrajectoryOptionsBuilder &spacecraft_lookup_dt(const double v) { opt_.spacecraft_lookup_dt_s = v; return *this; }
        TrajectoryOptionsBuilder &celestial_dt(const double v) { opt_.celestial_dt_s = v; return *this; }
        TrajectoryOptionsBuilder &max_samples(const std::size_t v) { opt_.max_samples = v; return *this; }
        TrajectoryOptionsBuilder &include_start(const bool v) { opt_.include_start = v; return *this; }
        TrajectoryOptionsBuilder &include_end(const bool v) { opt_.include_end = v; return *this; }
        TrajectoryOptionsBuilder &stop_on_impact(const bool v) { opt_.stop_on_impact = v; return *this; }

        operator TrajectoryOptions() const { return opt_; }
        TrajectoryOptions build() const { return opt_; }

    private:
        TrajectoryOptions opt_{};
    };

    inline TrajectorySegmentOptionsBuilder trajectory_segment_options() { return {}; }
    inline TrajectoryOptionsBuilder trajectory_options() { return {}; }

} // namespace orbitsim
