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
#include <limits>
#include <optional>
#include <utility>
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

        inline bool finite_state_(const State &state)
        {
            return detail::finite3_(state.position_m) && detail::finite3_(state.velocity_mps);
        }

        inline double log_lerp_(const double a, const double b, const double u)
        {
            const double a_safe = std::max(a, 1.0e-12);
            const double b_safe = std::max(b, 1.0e-12);
            const double alpha = clamp01(u);
            return std::exp(std::log(a_safe) + (std::log(b_safe) - std::log(a_safe)) * alpha);
        }

        template<typename Diagnostics>
        inline void reset_adaptive_diagnostics_(Diagnostics *diag)
        {
            if (diag == nullptr)
            {
                return;
            }
            *diag = Diagnostics{};
        }

        template<typename Diagnostics>
        inline void record_adaptive_accept_(Diagnostics *diag, const double dt_s)
        {
            if (diag == nullptr)
            {
                return;
            }

            diag->covered_duration_s += std::max(0.0, dt_s);
            ++diag->accepted_segments;
            if (diag->accepted_segments == 1)
            {
                diag->min_dt_s = dt_s;
                diag->max_dt_s = dt_s;
            }
            else
            {
                diag->min_dt_s = std::min(diag->min_dt_s, dt_s);
                diag->max_dt_s = std::max(diag->max_dt_s, dt_s);
            }
            if (diag->accepted_segments > 0)
            {
                diag->avg_dt_s = diag->covered_duration_s / static_cast<double>(diag->accepted_segments);
            }
        }

        template<typename Diagnostics>
        inline void record_adaptive_reject_(Diagnostics *diag)
        {
            if (diag != nullptr)
            {
                ++diag->rejected_splits;
            }
        }

        template<typename Diagnostics>
        inline void record_adaptive_forced_boundary_(Diagnostics *diag)
        {
            if (diag != nullptr)
            {
                ++diag->forced_boundary_splits;
            }
        }

        template<typename Diagnostics>
        inline void mark_adaptive_hard_cap_(Diagnostics *diag)
        {
            if (diag != nullptr)
            {
                diag->hard_cap_hit = true;
            }
        }

        template<typename Diagnostics>
        inline void mark_adaptive_cancelled_(Diagnostics *diag)
        {
            if (diag != nullptr)
            {
                diag->cancelled = true;
            }
        }

        inline bool cancel_requested_(const std::function<bool()> &cancel_requested)
        {
            return static_cast<bool>(cancel_requested) && cancel_requested();
        }

        struct LocalDynamicsScale
        {
            double r_m{1.0};
            double v_mps{1.0};
            double tau_dyn_s{1.0};
        };

        template<typename BodyStateAt>
        inline LocalDynamicsScale local_dynamics_scale_(const std::vector<MassiveBody> &bodies,
                                                        const State &query_state,
                                                        const double gravitational_constant,
                                                        const double softening_length_m,
                                                        BodyStateAt &&body_state_at,
                                                        const std::optional<std::size_t> skip_index = std::nullopt)
        {
            const double query_r = glm::length(query_state.position_m);
            const double query_v = glm::length(query_state.velocity_mps);

            LocalDynamicsScale fallback{};
            fallback.r_m = (std::isfinite(query_r) && query_r > 1.0) ? query_r : 1.0;
            fallback.v_mps = (std::isfinite(query_v) && query_v > 1.0e-6) ? query_v : 1.0;
            fallback.tau_dyn_s = std::max(1.0, fallback.r_m / fallback.v_mps);

            if (bodies.empty())
            {
                return fallback;
            }

            const double eps2 = softening_length_m * softening_length_m;
            std::size_t best_index = bodies.size();
            double best_metric = -1.0;
            Vec3 best_dr{0.0, 0.0, 0.0};
            Vec3 best_dv{0.0, 0.0, 0.0};

            for (std::size_t i = 0; i < bodies.size(); ++i)
            {
                if (skip_index.has_value() && *skip_index == i)
                {
                    continue;
                }

                const double mass_kg = std::isfinite(bodies[i].mass_kg) ? bodies[i].mass_kg : 0.0;
                if (!(mass_kg > 0.0))
                {
                    continue;
                }

                const State body_state = body_state_at(i);
                if (!finite_state_(body_state))
                {
                    continue;
                }

                const Vec3 dr = query_state.position_m - body_state.position_m;
                const double r2 = glm::dot(dr, dr) + eps2;
                if (!(r2 > 0.0) || !std::isfinite(r2))
                {
                    continue;
                }

                const double metric = mass_kg / r2;
                if (metric > best_metric)
                {
                    best_metric = metric;
                    best_index = i;
                    best_dr = dr;
                    best_dv = query_state.velocity_mps - body_state.velocity_mps;
                }
            }

            if (best_index >= bodies.size())
            {
                return fallback;
            }

            const double r_local = glm::length(best_dr);
            const double v_local = glm::length(best_dv);
            if (!(r_local > 0.0) || !std::isfinite(r_local))
            {
                return fallback;
            }

            const double mu_local = gravitational_constant * bodies[best_index].mass_kg;
            const double tau_adv_s =
                    (std::isfinite(v_local) && v_local > 1.0e-6) ? (r_local / v_local) : std::numeric_limits<double>::infinity();
            const double tau_orbit_s =
                    (std::isfinite(mu_local) && mu_local > 0.0) ? std::sqrt((r_local * r_local * r_local) / mu_local)
                                                                 : std::numeric_limits<double>::infinity();

            LocalDynamicsScale out{};
            out.r_m = r_local;
            out.v_mps = (std::isfinite(v_local) && v_local > 1.0e-6) ? v_local : fallback.v_mps;
            out.tau_dyn_s = std::min(tau_adv_s, tau_orbit_s);
            if (!(out.tau_dyn_s > 0.0) || !std::isfinite(out.tau_dyn_s))
            {
                out.tau_dyn_s = fallback.tau_dyn_s;
            }
            return out;
        }

        inline std::pair<double, double> adaptive_tolerance_at_(const AdaptiveToleranceRamp &tolerance,
                                                                const double time_alpha,
                                                                const LocalDynamicsScale &scale)
        {
            const double pos_tol_abs = log_lerp_(tolerance.pos_near_m, tolerance.pos_far_m, time_alpha);
            const double vel_tol_abs = log_lerp_(tolerance.vel_near_mps, tolerance.vel_far_mps, time_alpha);
            const double pos_tol = std::max(pos_tol_abs, tolerance.rel_pos_floor * std::max(scale.r_m, 1.0));
            const double vel_floor = pos_tol / std::max(scale.tau_dyn_s, 1.0);
            const double vel_scale = std::max(scale.v_mps, vel_floor);
            const double vel_tol = std::max(vel_tol_abs, tolerance.rel_vel_floor * vel_scale);
            return {pos_tol, vel_tol};
        }

        struct MidpointError
        {
            bool finite{false};
            double pos_error_m{std::numeric_limits<double>::infinity()};
            double vel_error_mps{std::numeric_limits<double>::infinity()};
        };

        inline MidpointError midpoint_error_(const State &start,
                                             const State &end,
                                             const State &midpoint_ref,
                                             const double dt_s)
        {
            MidpointError out{};
            if (!(dt_s > 0.0) || !finite_state_(start) || !finite_state_(end) || !finite_state_(midpoint_ref))
            {
                return out;
            }

            const Vec3 pos_mid = hermite_position(
                    start.position_m, start.velocity_mps, end.position_m, end.velocity_mps, dt_s, 0.5);
            const Vec3 vel_mid = hermite_velocity_mps(
                    start.position_m, start.velocity_mps, end.position_m, end.velocity_mps, dt_s, 0.5);

            if (!detail::finite3_(pos_mid) || !detail::finite3_(vel_mid))
            {
                return out;
            }

            out.finite = true;
            out.pos_error_m = glm::length(pos_mid - midpoint_ref.position_m);
            out.vel_error_mps = glm::length(vel_mid - midpoint_ref.velocity_mps);
            return out;
        }

        inline double error_ratio_(const MidpointError &error, const double pos_tol, const double vel_tol)
        {
            if (!error.finite)
            {
                return std::numeric_limits<double>::infinity();
            }

            const double pos_ratio = (pos_tol > 0.0) ? (error.pos_error_m / pos_tol) : std::numeric_limits<double>::infinity();
            const double vel_ratio =
                    (vel_tol > 0.0) ? (error.vel_error_mps / vel_tol) : std::numeric_limits<double>::infinity();
            return std::max(pos_ratio, vel_ratio);
        }

        inline bool adaptive_accept_interval_(const bool within_tolerance,
                                              const double error_ratio,
                                              const double dt_s,
                                              const double min_dt_s,
                                              const std::size_t accepted_segments,
                                              const std::size_t soft_max_segments,
                                              const std::size_t hard_max_segments,
                                              bool *out_hard_cap)
        {
            if (out_hard_cap != nullptr)
            {
                *out_hard_cap = false;
            }

            if (within_tolerance)
            {
                return true;
            }

            if (dt_s <= (min_dt_s * (1.0 + 1.0e-9)))
            {
                return true;
            }

            if (hard_max_segments > 0 && accepted_segments + 1 >= hard_max_segments)
            {
                if (out_hard_cap != nullptr)
                {
                    *out_hard_cap = true;
                }
                return true;
            }

            if (soft_max_segments > 0 && accepted_segments + 1 >= soft_max_segments && error_ratio <= 4.0)
            {
                return true;
            }

            return false;
        }

        inline State state_after_impulses_at_(const Spacecraft &sc_before,
                                              const std::vector<MassiveBody> &bodies,
                                              const CelestialEphemeris &eph,
                                              const ManeuverPlan &plan,
                                              const double gravitational_constant,
                                              const double softening_length_m,
                                              const double t_s,
                                              const SpacecraftStateLookup &sc_lookup)
        {
            State out = sc_before.state;
            if (plan.impulses.empty() || bodies.empty())
            {
                return out;
            }

            const auto auto_primary_at = [&](const double query_time_s, const Vec3 &pos_m) -> std::size_t {
                return auto_select_primary_index(
                        bodies, pos_m, [&](const std::size_t i) -> Vec3 { return eph.body_position_at(i, query_time_s); },
                        softening_length_m);
            };

            const double eps_t = 1.0e-9;
            for (const auto &impulse : plan.impulses)
            {
                if (!segment_applies_to_spacecraft(impulse, sc_before.id))
                {
                    continue;
                }
                if (!(std::abs(impulse.t_s - t_s) <= eps_t))
                {
                    continue;
                }

                const std::size_t primary = resolve_primary_index(
                        bodies,
                        impulse,
                        t_s,
                        out.position_m,
                        sc_lookup,
                        [&](const double query_time_s, const Vec3 &pos_m) -> std::size_t {
                            return auto_primary_at(query_time_s, pos_m);
                        });
                out.velocity_mps += rtn_vector_to_inertial(
                        eph,
                        bodies,
                        primary,
                        t_s,
                        out.position_m,
                        out.velocity_mps,
                        impulse.dv_rtn_mps,
                        impulse.rtn_frame,
                        sc_lookup);
            }

            return out;
        }

        struct BoundaryTime
        {
            double t_s{0.0};
            std::uint32_t flags{kTrajectorySegmentFlagNone};
        };

        inline std::vector<BoundaryTime> collect_adaptive_boundaries_(const GameSimulation &sim,
                                                                      const SpacecraftId spacecraft_id,
                                                                      const AdaptiveSegmentOptions &opt)
        {
            std::vector<BoundaryTime> out;
            if (!(opt.duration_s > 0.0))
            {
                return out;
            }

            const double t0 = sim.time_s();
            const double t_end = t0 + opt.duration_s;
            out.reserve(2 + opt.forced_split_times_s.size() + sim.maneuver_plan().impulses.size() +
                        sim.maneuver_plan().segments.size() * 2);
            out.push_back(BoundaryTime{.t_s = t0, .flags = kTrajectorySegmentFlagNone});

            for (const double t_split : opt.forced_split_times_s)
            {
                if (std::isfinite(t_split) && t_split > t0 && t_split < t_end)
                {
                    out.push_back(BoundaryTime{.t_s = t_split, .flags = kTrajectorySegmentFlagForcedBoundary});
                }
            }

            if (opt.split_at_impulses)
            {
                for (const auto &impulse : sim.maneuver_plan().impulses)
                {
                    if (!segment_applies_to_spacecraft(impulse, spacecraft_id))
                    {
                        continue;
                    }
                    if (std::isfinite(impulse.t_s) && impulse.t_s > t0 && impulse.t_s < t_end)
                    {
                        out.push_back(BoundaryTime{.t_s = impulse.t_s, .flags = kTrajectorySegmentFlagImpulseBoundary});
                    }
                }
            }

            if (opt.split_at_burn_boundaries)
            {
                for (const auto &segment : sim.maneuver_plan().segments)
                {
                    if (!segment_applies_to_spacecraft(segment, spacecraft_id))
                    {
                        continue;
                    }
                    if (std::isfinite(segment.t_start_s) && segment.t_start_s > t0 && segment.t_start_s < t_end)
                    {
                        out.push_back(BoundaryTime{.t_s = segment.t_start_s, .flags = kTrajectorySegmentFlagBurnBoundary});
                    }
                    if (std::isfinite(segment.t_end_s) && segment.t_end_s > t0 && segment.t_end_s < t_end)
                    {
                        out.push_back(BoundaryTime{.t_s = segment.t_end_s, .flags = kTrajectorySegmentFlagBurnBoundary});
                    }
                }
            }

            out.push_back(BoundaryTime{.t_s = t_end, .flags = kTrajectorySegmentFlagNone});

            std::sort(out.begin(),
                      out.end(),
                      [](const BoundaryTime &a, const BoundaryTime &b) { return a.t_s < b.t_s; });

            std::vector<BoundaryTime> merged;
            merged.reserve(out.size());
            for (const BoundaryTime &boundary : out)
            {
                if (!merged.empty() && std::abs(merged.back().t_s - boundary.t_s) <= 1.0e-9)
                {
                    merged.back().flags |= boundary.flags;
                    continue;
                }
                merged.push_back(boundary);
            }

            return merged;
        }

        inline AdaptiveEphemerisOptions sanitize_adaptive_ephemeris_options_(AdaptiveEphemerisOptions opt)
        {
            opt.min_dt_s = std::max(1.0e-6, opt.min_dt_s);
            if (!(opt.duration_s > 0.0))
            {
                return opt;
            }
            if (!(opt.max_dt_s > 0.0))
            {
                opt.max_dt_s = opt.duration_s;
            }
            opt.max_dt_s = std::max(opt.max_dt_s, opt.min_dt_s);
            if (opt.soft_max_segments == 0)
            {
                opt.soft_max_segments = std::max<std::size_t>(1, opt.hard_max_segments);
            }
            if (opt.hard_max_segments == 0)
            {
                opt.hard_max_segments = std::max<std::size_t>(opt.soft_max_segments, 1);
            }
            opt.hard_max_segments = std::max(opt.hard_max_segments, opt.soft_max_segments);
            return opt;
        }

        inline AdaptiveSegmentOptions sanitize_adaptive_segment_options_(AdaptiveSegmentOptions opt)
        {
            opt.min_dt_s = std::max(1.0e-6, opt.min_dt_s);
            if (!(opt.duration_s > 0.0))
            {
                return opt;
            }
            if (!(opt.max_dt_s > 0.0))
            {
                opt.max_dt_s = opt.duration_s;
            }
            opt.max_dt_s = std::max(opt.max_dt_s, opt.min_dt_s);
            if (!(opt.lookup_max_dt_s > 0.0))
            {
                opt.lookup_max_dt_s = opt.max_dt_s;
            }
            opt.lookup_max_dt_s = std::max(opt.lookup_max_dt_s, opt.min_dt_s);
            if (opt.soft_max_segments == 0)
            {
                opt.soft_max_segments = std::max<std::size_t>(1, opt.hard_max_segments);
            }
            if (opt.hard_max_segments == 0)
            {
                opt.hard_max_segments = std::max<std::size_t>(opt.soft_max_segments, 1);
            }
            opt.hard_max_segments = std::max(opt.hard_max_segments, opt.soft_max_segments);
            return opt;
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

        inline std::optional<double> find_impact_time_(
            const Spacecraft &sc0, const Spacecraft &sc1,
            const std::vector<MassiveBody> &bodies, const CelestialEphemeris &eph,
            const ManeuverPlan &plan, const double G, const double softening,
            const DOPRI5Options &integrator, const double t0, const double h,
            const SpacecraftStateLookup &sc_lookup, const EventOptions &event_opt);

        struct AdaptiveEphemerisStepResult
        {
            std::vector<MassiveBody> midpoint_massive{};
            std::vector<MassiveBody> end_massive{};
            std::vector<State> start_states{};
            std::vector<State> midpoint_states{};
            std::vector<State> end_states{};
        };

        inline bool simulate_adaptive_ephemeris_step_(const GameSimulation &sim,
                                                      const std::vector<MassiveBody> &start_massive,
                                                      const double dt_s,
                                                      AdaptiveEphemerisStepResult *out)
        {
            if (out == nullptr || !(dt_s > 0.0))
            {
                return false;
            }

            out->midpoint_massive = start_massive;
            out->end_massive = start_massive;
            snapshot_states(start_massive, &out->start_states);

            symplectic4_step(
                    out->midpoint_massive,
                    0.5 * dt_s,
                    sim.config().gravitational_constant,
                    sim.config().softening_length_m);
            snapshot_states(out->midpoint_massive, &out->midpoint_states);

            out->end_massive = out->midpoint_massive;
            symplectic4_step(
                    out->end_massive,
                    0.5 * dt_s,
                    sim.config().gravitational_constant,
                    sim.config().softening_length_m);
            snapshot_states(out->end_massive, &out->end_states);

            return out->start_states.size() == out->midpoint_states.size() &&
                   out->start_states.size() == out->end_states.size();
        }

        inline bool append_adaptive_ephemeris_interval_(const GameSimulation &sim,
                                                        const AdaptiveEphemerisOptions &opt,
                                                        const double base_t0_s,
                                                        const std::vector<MassiveBody> &start_massive,
                                                        const double t0_s,
                                                        const double dt_s,
                                                        std::vector<CelestialEphemerisSegment> *segments_out,
                                                        AdaptiveEphemerisDiagnostics *diag)
        {
            if (segments_out == nullptr || !(dt_s > 0.0))
            {
                return false;
            }
            if (cancel_requested_(opt.cancel_requested))
            {
                mark_adaptive_cancelled_(diag);
                return false;
            }

            AdaptiveEphemerisStepResult step{};
            if (!simulate_adaptive_ephemeris_step_(sim, start_massive, dt_s, &step))
            {
                return false;
            }

            const double t_mid_s = t0_s + (0.5 * dt_s);
            const double time_alpha =
                    (opt.duration_s > 0.0) ? ((t_mid_s - base_t0_s) / opt.duration_s) : 0.0;

            bool within_tolerance = true;
            double error_ratio = 0.0;
            for (std::size_t i = 0; i < step.start_states.size(); ++i)
            {
                const MidpointError error = midpoint_error_(
                        step.start_states[i], step.end_states[i], step.midpoint_states[i], dt_s);
                const LocalDynamicsScale scale = local_dynamics_scale_(
                        start_massive,
                        step.midpoint_states[i],
                        sim.config().gravitational_constant,
                        sim.config().softening_length_m,
                        [&](const std::size_t j) -> State { return step.midpoint_states[j]; },
                        i);
                const auto [pos_tol, vel_tol] = adaptive_tolerance_at_(opt.tolerance, time_alpha, scale);
                const double ratio = error_ratio_(error, pos_tol, vel_tol);
                error_ratio = std::max(error_ratio, ratio);
                within_tolerance =
                        within_tolerance && error.finite && error.pos_error_m <= pos_tol && error.vel_error_mps <= vel_tol;
            }

            const bool dt_exceeds_cap = opt.max_dt_s > 0.0 && dt_s > (opt.max_dt_s * (1.0 + 1.0e-9));
            bool hard_cap_accept = false;
            const bool accept = !dt_exceeds_cap &&
                                (within_tolerance ||
                                 adaptive_accept_interval_(
                                         within_tolerance,
                                         error_ratio,
                                         dt_s,
                                         std::max(1.0e-6, opt.min_dt_s),
                                         segments_out->size(),
                                         opt.soft_max_segments,
                                         opt.hard_max_segments,
                                         &hard_cap_accept));

            if (accept || !(0.5 * dt_s > 0.0))
            {
                if (hard_cap_accept)
                {
                    mark_adaptive_hard_cap_(diag);
                }
                segments_out->push_back(CelestialEphemerisSegment{
                        .t0_s = t0_s,
                        .dt_s = dt_s,
                        .start = step.start_states,
                        .end = step.end_states,
                });
                record_adaptive_accept_(diag, dt_s);
                return true;
            }

            record_adaptive_reject_(diag);
            const std::size_t left_segment_count = segments_out->size();
            if (!append_adaptive_ephemeris_interval_(
                        sim,
                        opt,
                        base_t0_s,
                        start_massive,
                        t0_s,
                        0.5 * dt_s,
                        segments_out,
                        diag))
            {
                return false;
            }
            std::vector<MassiveBody> refined_midpoint_massive = step.midpoint_massive;
            if (segments_out->size() > left_segment_count &&
                segments_out->back().end.size() == refined_midpoint_massive.size())
            {
                const std::vector<State> &refined_midpoint_states = segments_out->back().end;
                for (std::size_t i = 0; i < refined_midpoint_massive.size(); ++i)
                {
                    refined_midpoint_massive[i].state = refined_midpoint_states[i];
                }
            }
            return append_adaptive_ephemeris_interval_(
                    sim,
                    opt,
                    base_t0_s,
                    refined_midpoint_massive,
                    t_mid_s,
                    0.5 * dt_s,
                    segments_out,
                    diag);
        }

        inline TrajectorySample make_sample_(const double t_s, const State &s)
        {
            return TrajectorySample{.t_s = t_s, .position_m = s.position_m, .velocity_mps = s.velocity_mps};
        }

        inline bool push_adaptive_spacecraft_segment_(const State &segment_start,
                                                      const State &segment_end,
                                                      const double t0_s,
                                                      const double dt_s,
                                                      const std::uint32_t flags,
                                                      std::vector<TrajectorySegment> *segments_out,
                                                      AdaptiveSegmentDiagnostics *diag)
        {
            if (segments_out == nullptr || !(dt_s > 0.0) || !finite_state_(segment_start) || !finite_state_(segment_end))
            {
                return false;
            }

            segments_out->push_back(TrajectorySegment{
                    .t0_s = t0_s,
                    .dt_s = dt_s,
                    .start = segment_start,
                    .end = segment_end,
                    .flags = flags,
            });
            record_adaptive_accept_(diag, dt_s);
            return true;
        }

        inline bool append_adaptive_spacecraft_interval_(const GameSimulation &sim,
                                                         const CelestialEphemeris &eph,
                                                         const AdaptiveSegmentOptions &opt,
                                                         const Spacecraft &sc0,
                                                         const SpacecraftStateLookup &sc_lookup,
                                                         const double base_t0_s,
                                                         const double t0_s,
                                                         const double dt_s,
                                                         const std::uint32_t start_flags,
                                                         std::vector<TrajectorySegment> *segments_out,
                                                         AdaptiveSegmentDiagnostics *diag,
                                                         Spacecraft *out_end_spacecraft,
                                                         bool *out_terminal_impact)
        {
            if (segments_out == nullptr || out_end_spacecraft == nullptr || !(dt_s > 0.0))
            {
                return false;
            }
            if (out_terminal_impact != nullptr)
            {
                *out_terminal_impact = false;
            }
            if (cancel_requested_(opt.cancel_requested))
            {
                mark_adaptive_cancelled_(diag);
                return false;
            }

            const double min_dt_s = std::max(1.0e-6, opt.min_dt_s);
            const double lookup_dt_s = (opt.lookup_max_dt_s > 0.0) ? opt.lookup_max_dt_s : std::max(opt.max_dt_s, min_dt_s);

            const State segment_start_state = state_after_impulses_at_(
                    sc0,
                    sim.massive_bodies(),
                    eph,
                    sim.maneuver_plan(),
                    sim.config().gravitational_constant,
                    sim.config().softening_length_m,
                    t0_s,
                    sc_lookup);

            const Spacecraft sc_mid = propagate_spacecraft_in_ephemeris(
                    sc0,
                    sim.massive_bodies(),
                    eph,
                    sim.maneuver_plan(),
                    sim.config().gravitational_constant,
                    sim.config().softening_length_m,
                    sim.config().spacecraft_integrator,
                    t0_s,
                    0.5 * dt_s,
                    sc_lookup);
            const Spacecraft sc_end = propagate_spacecraft_in_ephemeris(
                    sc_mid,
                    sim.massive_bodies(),
                    eph,
                    sim.maneuver_plan(),
                    sim.config().gravitational_constant,
                    sim.config().softening_length_m,
                    sim.config().spacecraft_integrator,
                    t0_s + (0.5 * dt_s),
                    0.5 * dt_s,
                    sc_lookup);

            if (opt.stop_on_impact)
            {
                const auto impact_t = find_impact_time_(
                        sc0,
                        sc_end,
                        sim.massive_bodies(),
                        eph,
                        sim.maneuver_plan(),
                        sim.config().gravitational_constant,
                        sim.config().softening_length_m,
                        sim.config().spacecraft_integrator,
                        t0_s,
                        dt_s,
                        sc_lookup,
                        sim.config().events);
                if (impact_t.has_value() && *impact_t >= t0_s && *impact_t <= (t0_s + dt_s))
                {
                    const double impact_dt_s = *impact_t - t0_s;
                    if (impact_dt_s > 0.0)
                    {
                        const Spacecraft sc_imp = propagate_spacecraft_in_ephemeris(
                                sc0,
                                sim.massive_bodies(),
                                eph,
                                sim.maneuver_plan(),
                                sim.config().gravitational_constant,
                                sim.config().softening_length_m,
                                sim.config().spacecraft_integrator,
                                t0_s,
                                impact_dt_s,
                                sc_lookup);
                        if (!push_adaptive_spacecraft_segment_(
                                    segment_start_state,
                                    sc_imp.state,
                                    t0_s,
                                    impact_dt_s,
                                    start_flags | kTrajectorySegmentFlagImpactTerminal,
                                    segments_out,
                                    diag))
                        {
                            return false;
                        }
                        *out_end_spacecraft = sc_imp;
                        if (out_terminal_impact != nullptr)
                        {
                            *out_terminal_impact = true;
                        }
                        return true;
                    }
                }
            }

            const double t_mid_s = t0_s + (0.5 * dt_s);
            const double time_alpha =
                    (opt.duration_s > 0.0) ? ((t_mid_s - base_t0_s) / opt.duration_s) : 0.0;

            MidpointError error = midpoint_error_(segment_start_state, sc_end.state, sc_mid.state, dt_s);
            const LocalDynamicsScale scale = local_dynamics_scale_(
                    sim.massive_bodies(),
                    sc_mid.state,
                    sim.config().gravitational_constant,
                    sim.config().softening_length_m,
                    [&](const std::size_t i) -> State {
                        std::size_t eph_index = 0;
                        if (eph.body_index_for_id(sim.massive_bodies()[i].id, &eph_index))
                        {
                            return eph.body_state_at(eph_index, t_mid_s);
                        }
                        return sim.massive_bodies()[i].state;
                    });
            const auto [pos_tol, vel_tol] = adaptive_tolerance_at_(opt.tolerance, time_alpha, scale);
            const double ratio = error_ratio_(error, pos_tol, vel_tol);
            const bool within_tolerance = error.finite && error.pos_error_m <= pos_tol && error.vel_error_mps <= vel_tol;
            const bool dt_exceeds_cap = opt.max_dt_s > 0.0 && dt_s > (opt.max_dt_s * (1.0 + 1.0e-9));

            bool hard_cap_accept = false;
            const bool accept = !dt_exceeds_cap &&
                                (within_tolerance ||
                                 adaptive_accept_interval_(
                                         within_tolerance,
                                         ratio,
                                         dt_s,
                                         min_dt_s,
                                         segments_out->size(),
                                         opt.soft_max_segments,
                                         opt.hard_max_segments,
                                         &hard_cap_accept));
            if (accept || !(0.5 * dt_s > 0.0) || lookup_dt_s <= 0.0)
            {
                if (hard_cap_accept)
                {
                    mark_adaptive_hard_cap_(diag);
                }
                if (!push_adaptive_spacecraft_segment_(
                            segment_start_state,
                            sc_end.state,
                            t0_s,
                            dt_s,
                            start_flags,
                            segments_out,
                            diag))
                {
                    return false;
                }
                *out_end_spacecraft = sc_end;
                return true;
            }

            record_adaptive_reject_(diag);

            Spacecraft left_end{};
            bool left_terminal = false;
            if (!append_adaptive_spacecraft_interval_(
                        sim,
                        eph,
                        opt,
                        sc0,
                        sc_lookup,
                        base_t0_s,
                        t0_s,
                        0.5 * dt_s,
                        start_flags,
                        segments_out,
                        diag,
                        &left_end,
                        &left_terminal))
            {
                return false;
            }
            if (left_terminal)
            {
                *out_end_spacecraft = left_end;
                if (out_terminal_impact != nullptr)
                {
                    *out_terminal_impact = true;
                }
                return true;
            }

            Spacecraft right_end{};
            bool right_terminal = false;
            if (!append_adaptive_spacecraft_interval_(
                        sim,
                        eph,
                        opt,
                        left_end,
                        sc_lookup,
                        base_t0_s,
                        t_mid_s,
                        0.5 * dt_s,
                        kTrajectorySegmentFlagNone,
                        segments_out,
                        diag,
                        &right_end,
                        &right_terminal))
            {
                return false;
            }

            *out_end_spacecraft = right_end;
            if (out_terminal_impact != nullptr)
            {
                *out_terminal_impact = right_terminal;
            }
            return true;
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
                const auto eval_f = [&](const double t_s) -> double {
                    const Spacecraft sc_at = propagate_spacecraft_in_ephemeris(
                        sc0, bodies, eph, plan, G, softening, integrator, t0, t_s - t0, sc_lookup);
                    return glm::length(sc_at.state.position_m - eph.body_position_at_by_id(body.id, t_s))
                           - threshold_m;
                };

                const std::optional<detail::CrossingBracket> bracket =
                        detail::first_crossing_bracket_s(t0, t0 + h, f0, f1, event_opt, eval_f);
                if (!bracket.has_value())
                {
                    continue;
                }

                if (detail::crossing_from_bracket_(*bracket, event_opt, eval_f) != Crossing::Enter)
                {
                    continue;
                }

                const double t_event = bisect_crossing_time_s(
                    bracket->t0_s, bracket->t1_s, bracket->f0, event_opt, eval_f);

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
            return trajectory_segment_state_at(seg, t_s);
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

    inline CelestialEphemeris build_celestial_ephemeris_adaptive(
            const GameSimulation &sim,
            const AdaptiveEphemerisOptions &opt,
            AdaptiveEphemerisDiagnostics *out_diag = nullptr)
    {
        detail::reset_adaptive_diagnostics_(out_diag);

        AdaptiveEphemerisOptions clean_opt = detail::sanitize_adaptive_ephemeris_options_(opt);
        CelestialEphemeris eph;
        if (!(clean_opt.duration_s > 0.0) || sim.massive_bodies().empty())
        {
            return eph;
        }

        std::vector<BodyId> body_ids;
        body_ids.reserve(sim.massive_bodies().size());
        for (const auto &body : sim.massive_bodies())
        {
            body_ids.push_back(body.id);
        }
        eph.set_body_ids(std::move(body_ids));

        const std::vector<MassiveBody> start_massive = sim.massive_bodies();
        detail::append_adaptive_ephemeris_interval_(
                sim,
                clean_opt,
                sim.time_s(),
                start_massive,
                sim.time_s(),
                clean_opt.duration_s,
                &eph.segments,
                out_diag);
        return eph;
    }

    inline std::vector<TrajectorySegment> predict_spacecraft_trajectory_segments_adaptive(
            const GameSimulation &sim,
            const CelestialEphemeris &eph,
            const SpacecraftId spacecraft_id,
            const AdaptiveSegmentOptions &opt,
            AdaptiveSegmentDiagnostics *out_diag = nullptr)
    {
        detail::reset_adaptive_diagnostics_(out_diag);

        std::vector<TrajectorySegment> out;
        const Spacecraft *sc_ptr = sim.spacecraft_by_id(spacecraft_id);
        if (sc_ptr == nullptr)
        {
            return out;
        }

        AdaptiveSegmentOptions clean_opt = detail::sanitize_adaptive_segment_options_(opt);
        if (!(clean_opt.duration_s > 0.0))
        {
            return out;
        }

        const double t0_s = sim.time_s();
        const double t_end_s = t0_s + clean_opt.duration_s;
        const double lookup_dt_s =
                (clean_opt.lookup_max_dt_s > 0.0) ? clean_opt.lookup_max_dt_s : std::max(clean_opt.max_dt_s, clean_opt.min_dt_s);

        SpacecraftStateCache<CelestialEphemeris> sc_cache(
                sim.massive_bodies(),
                eph,
                sim.maneuver_plan(),
                sim.config().gravitational_constant,
                sim.config().softening_length_m,
                sim.config().spacecraft_integrator,
                t0_s,
                t_end_s,
                [&sim](const SpacecraftId id) { return sim.spacecraft_by_id(id); },
                typename SpacecraftStateCache<CelestialEphemeris>::Options{
                        .lookup_dt_s = lookup_dt_s,
                        .max_segments = std::max<std::size_t>(clean_opt.hard_max_segments, 1),
                });
        const SpacecraftStateLookup sc_lookup = sc_cache.lookup();

        const std::vector<detail::BoundaryTime> boundaries =
                detail::collect_adaptive_boundaries_(sim, spacecraft_id, clean_opt);
        if (boundaries.size() < 2)
        {
            return out;
        }

        Spacecraft current_sc = *sc_ptr;
        bool terminal_impact = false;
        for (std::size_t i = 1; i < boundaries.size(); ++i)
        {
            if (detail::cancel_requested_(clean_opt.cancel_requested))
            {
                detail::mark_adaptive_cancelled_(out_diag);
                break;
            }

            const double interval_t0_s = boundaries[i - 1].t_s;
            const double interval_dt_s = boundaries[i].t_s - boundaries[i - 1].t_s;
            if (!(interval_dt_s > 0.0))
            {
                continue;
            }

            if (boundaries[i - 1].flags != kTrajectorySegmentFlagNone)
            {
                detail::record_adaptive_forced_boundary_(out_diag);
            }

            Spacecraft interval_end{};
            if (!detail::append_adaptive_spacecraft_interval_(
                        sim,
                        eph,
                        clean_opt,
                        current_sc,
                        sc_lookup,
                        t0_s,
                        interval_t0_s,
                        interval_dt_s,
                        boundaries[i - 1].flags,
                        &out,
                        out_diag,
                        &interval_end,
                        &terminal_impact))
            {
                break;
            }

            current_sc = interval_end;
            if (terminal_impact)
            {
                break;
            }
        }

        return out;
    }

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
