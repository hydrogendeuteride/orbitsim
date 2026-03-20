#pragma once

#include "orbitsim/coordinate_frames.hpp"
#include "orbitsim/frame_spec.hpp"
#include "orbitsim/spacecraft_lookup.hpp"
#include "orbitsim/synodic.hpp"
#include "orbitsim/trajectories.hpp"
#include "orbitsim/trajectory_segments.hpp"
#include "orbitsim/trajectory_types.hpp"

#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace orbitsim
{
    namespace detail
    {
        inline FrameSegmentTransformOptions sanitize_frame_segment_transform_options_(FrameSegmentTransformOptions opt)
        {
            opt.min_dt_s = std::max(1.0e-6, opt.min_dt_s);
            if (!(opt.max_dt_s > 0.0))
            {
                opt.max_dt_s = opt.min_dt_s;
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

        inline const MassiveBody *body_by_id_(const std::vector<MassiveBody> &bodies, const BodyId id)
        {
            if (id == kInvalidBodyId)
            {
                return nullptr;
            }
            for (const auto &body : bodies)
            {
                if (body.id == id)
                {
                    return &body;
                }
            }
            return nullptr;
        }

        inline std::optional<State> inertial_state_to_frame_spec_state_(
                const State &state_in,
                const double t_s,
                const CelestialEphemeris &eph,
                const std::vector<MassiveBody> &bodies,
                const TrajectoryFrameSpec &frame,
                const SpacecraftStateLookup &sc_lookup)
        {
            switch (frame.type)
            {
            case TrajectoryFrameType::Inertial:
                return state_in;
            case TrajectoryFrameType::BodyCenteredInertial:
            {
                const MassiveBody *body = body_by_id_(bodies, frame.primary_body_id);
                if (body == nullptr)
                {
                    return std::nullopt;
                }
                return inertial_state_to_frame(state_in, make_body_centered_inertial_frame_at(eph, *body, t_s));
            }
            case TrajectoryFrameType::BodyFixed:
            {
                const MassiveBody *body = body_by_id_(bodies, frame.primary_body_id);
                if (body == nullptr)
                {
                    return std::nullopt;
                }
                const std::optional<RotatingFrame> body_fixed = make_body_fixed_frame_at(eph, *body, t_s);
                if (!body_fixed.has_value())
                {
                    return std::nullopt;
                }
                return inertial_state_to_frame(state_in, *body_fixed);
            }
            case TrajectoryFrameType::Synodic:
            {
                const MassiveBody *body_a = body_by_id_(bodies, frame.primary_body_id);
                const MassiveBody *body_b = body_by_id_(bodies, frame.secondary_body_id);
                if (body_a == nullptr || body_b == nullptr)
                {
                    return std::nullopt;
                }
                const std::optional<SynodicFrame> synodic = make_synodic_frame_at(eph, *body_a, *body_b, t_s);
                if (!synodic.has_value())
                {
                    return std::nullopt;
                }
                return inertial_state_to_frame(state_in, *synodic);
            }
            case TrajectoryFrameType::LVLH:
            {
                if (!sc_lookup || frame.target_spacecraft_id == kInvalidSpacecraftId || bodies.empty())
                {
                    return std::nullopt;
                }

                const std::optional<State> target_state = sc_lookup(frame.target_spacecraft_id, t_s);
                if (!target_state.has_value())
                {
                    return std::nullopt;
                }

                const MassiveBody *primary = body_by_id_(bodies, frame.primary_body_id);
                if (primary == nullptr)
                {
                    const std::size_t primary_index = auto_select_primary_index(
                            bodies,
                            target_state->position_m,
                            [&](const std::size_t i) -> Vec3 {
                                std::size_t eph_index = 0;
                                if (eph.body_index_for_id(bodies[i].id, &eph_index))
                                {
                                    return eph.body_position_at(eph_index, t_s);
                                }
                                return bodies[i].state.position_m;
                            });
                    primary = &bodies[primary_index];
                }

                State primary_state = primary->state;
                if (!eph.empty())
                {
                    primary_state = eph.body_state_at_by_id(primary->id, t_s);
                }

                const std::optional<RotatingFrame> lvlh = make_lvlh_frame(primary_state, *target_state);
                if (!lvlh.has_value())
                {
                    return std::nullopt;
                }
                return inertial_state_to_frame(state_in, *lvlh);
            }
            }

            return std::nullopt;
        }

        inline LocalDynamicsScale frame_state_scale_(const State &state)
        {
            const double r = glm::length(state.position_m);
            const double v = glm::length(state.velocity_mps);
            LocalDynamicsScale out{};
            out.r_m = (std::isfinite(r) && r > 1.0) ? r : 1.0;
            out.v_mps = (std::isfinite(v) && v > 1.0e-6) ? v : 1.0;
            out.tau_dyn_s = std::max(1.0, out.r_m / out.v_mps);
            return out;
        }

        inline TrajectorySegment make_inertial_subsegment_(const TrajectorySegment &segment,
                                                           const double sub_t0_s,
                                                           const double sub_t1_s)
        {
            return TrajectorySegment{
                    .t0_s = sub_t0_s,
                    .dt_s = sub_t1_s - sub_t0_s,
                    .start = state_at_segment_time_(segment, sub_t0_s),
                    .end = state_at_segment_time_(segment, sub_t1_s),
                    .flags = kTrajectorySegmentFlagNone,
            };
        }

        inline bool append_transformed_segment_(const TrajectorySegment &segment_inertial,
                                                const CelestialEphemeris &eph,
                                                const std::vector<MassiveBody> &bodies,
                                                const TrajectoryFrameSpec &frame,
                                                const FrameSegmentTransformOptions &opt,
                                                const double base_t0_s,
                                                const double base_duration_s,
                                                const SpacecraftStateLookup &sc_lookup,
                                                const std::uint32_t start_flags,
                                                const std::uint32_t extra_flags,
                                                std::vector<TrajectorySegment> *segments_out,
                                                FrameSegmentTransformDiagnostics *diag)
        {
            if (segments_out == nullptr || !(segment_inertial.dt_s > 0.0))
            {
                return false;
            }
            if (cancel_requested_(opt.cancel_requested))
            {
                mark_adaptive_cancelled_(diag);
                return false;
            }

            const double segment_t1_s = segment_inertial.t0_s + segment_inertial.dt_s;
            for (const double forced_t_s : opt.forced_split_times_s)
            {
                if (!std::isfinite(forced_t_s) || !(forced_t_s > segment_inertial.t0_s) || !(forced_t_s < segment_t1_s))
                {
                    continue;
                }

                const TrajectorySegment left = make_inertial_subsegment_(segment_inertial, segment_inertial.t0_s, forced_t_s);
                const TrajectorySegment right = make_inertial_subsegment_(segment_inertial, forced_t_s, segment_t1_s);
                record_adaptive_forced_boundary_(diag);
                if (!append_transformed_segment_(
                            left,
                            eph,
                            bodies,
                            frame,
                            opt,
                            base_t0_s,
                            base_duration_s,
                            sc_lookup,
                            start_flags,
                            extra_flags,
                            segments_out,
                            diag))
                {
                    return false;
                }
                return append_transformed_segment_(
                        right,
                        eph,
                        bodies,
                        frame,
                        opt,
                        base_t0_s,
                        base_duration_s,
                        sc_lookup,
                        kTrajectorySegmentFlagForcedBoundary,
                        extra_flags,
                        segments_out,
                        diag);
            }

            const std::optional<State> start_state = inertial_state_to_frame_spec_state_(
                    segment_inertial.start, segment_inertial.t0_s, eph, bodies, frame, sc_lookup);
            const std::optional<State> end_state = inertial_state_to_frame_spec_state_(
                    segment_inertial.end, segment_t1_s, eph, bodies, frame, sc_lookup);
            const State midpoint_inertial = state_at_segment_time_(segment_inertial, segment_inertial.t0_s + 0.5 * segment_inertial.dt_s);
            const std::optional<State> midpoint_state = inertial_state_to_frame_spec_state_(
                    midpoint_inertial,
                    segment_inertial.t0_s + 0.5 * segment_inertial.dt_s,
                    eph,
                    bodies,
                    frame,
                    sc_lookup);
            if (!start_state.has_value() || !end_state.has_value() || !midpoint_state.has_value())
            {
                return false;
            }

            const double time_alpha =
                    (base_duration_s > 0.0) ? (((segment_inertial.t0_s + 0.5 * segment_inertial.dt_s) - base_t0_s) / base_duration_s)
                                            : 0.0;
            const MidpointError error = midpoint_error_(
                    *start_state,
                    *end_state,
                    *midpoint_state,
                    segment_inertial.dt_s);
            const LocalDynamicsScale scale = frame_state_scale_(*midpoint_state);
            const auto [pos_tol, vel_tol] = adaptive_tolerance_at_(opt.tolerance, time_alpha, scale);
            const double ratio = error_ratio_(error, pos_tol, vel_tol);
            const bool within_tolerance = error.finite && error.pos_error_m <= pos_tol && error.vel_error_mps <= vel_tol;
            const bool dt_exceeds_cap = opt.max_dt_s > 0.0 && segment_inertial.dt_s > (opt.max_dt_s * (1.0 + 1.0e-9));

            bool hard_cap_accept = false;
            const bool accept = (!dt_exceeds_cap && within_tolerance) ||
                                adaptive_accept_interval_(
                                        within_tolerance && !dt_exceeds_cap,
                                        ratio,
                                        segment_inertial.dt_s,
                                        opt.min_dt_s,
                                        segments_out->size(),
                                        opt.soft_max_segments,
                                        opt.hard_max_segments,
                                        &hard_cap_accept);
            if (accept || !(0.5 * segment_inertial.dt_s > 0.0))
            {
                if (hard_cap_accept)
                {
                    mark_adaptive_hard_cap_(diag);
                }
                segments_out->push_back(TrajectorySegment{
                        .t0_s = segment_inertial.t0_s,
                        .dt_s = segment_inertial.dt_s,
                        .start = *start_state,
                        .end = *end_state,
                        .flags = start_flags | extra_flags,
                });
                record_adaptive_accept_(diag, segment_inertial.dt_s);
                return true;
            }

            record_adaptive_reject_(diag);
            if (diag != nullptr)
            {
                ++diag->frame_resegmentation_count;
            }

            const double split_t_s = segment_inertial.t0_s + 0.5 * segment_inertial.dt_s;
            const TrajectorySegment left = make_inertial_subsegment_(segment_inertial, segment_inertial.t0_s, split_t_s);
            const TrajectorySegment right = make_inertial_subsegment_(segment_inertial, split_t_s, segment_t1_s);
            const std::uint32_t split_flags = extra_flags | kTrajectorySegmentFlagFrameResegmented;

            if (!append_transformed_segment_(
                        left,
                        eph,
                        bodies,
                        frame,
                        opt,
                        base_t0_s,
                        base_duration_s,
                        sc_lookup,
                        start_flags,
                        split_flags,
                        segments_out,
                        diag))
            {
                return false;
            }
            return append_transformed_segment_(
                    right,
                    eph,
                    bodies,
                    frame,
                    opt,
                    base_t0_s,
                    base_duration_s,
                    sc_lookup,
                    kTrajectorySegmentFlagNone,
                    split_flags,
                    segments_out,
                    diag);
        }
    } // namespace detail

    // -------------------------------------------------------------------------
    // Trajectory <-> frame transforms (display/analysis only; does not affect dynamics).
    // -------------------------------------------------------------------------

    /// @brief Convert a single inertial trajectory sample into a (possibly rotating) frame.
    inline TrajectorySample inertial_sample_to_frame(const TrajectorySample &sample_in, const RotatingFrame &frame)
    {
        const State s_in = make_state(sample_in.position_m, sample_in.velocity_mps);
        const State s_out = inertial_state_to_frame(s_in, frame);
        return TrajectorySample{.t_s = sample_in.t_s, .position_m = s_out.position_m, .velocity_mps = s_out.velocity_mps};
    }

    /// @brief Convert a single trajectory sample expressed in a frame back to inertial coordinates.
    inline TrajectorySample frame_sample_to_inertial(const TrajectorySample &sample_frame, const RotatingFrame &frame)
    {
        const State s_frame = make_state(sample_frame.position_m, sample_frame.velocity_mps);
        const State s_out = frame_state_to_inertial(s_frame, frame);
        return TrajectorySample{.t_s = sample_frame.t_s, .position_m = s_out.position_m, .velocity_mps = s_out.velocity_mps};
    }

    /// @brief Convert a single trajectory sample from one frame to another (via inertial).
    inline TrajectorySample frame_sample_to_frame(const TrajectorySample &sample_in, const RotatingFrame &frame_from,
                                                  const RotatingFrame &frame_to)
    {
        const State s_from = make_state(sample_in.position_m, sample_in.velocity_mps);
        const State s_to = frame_state_to_frame(s_from, frame_from, frame_to);
        return TrajectorySample{.t_s = sample_in.t_s, .position_m = s_to.position_m, .velocity_mps = s_to.velocity_mps};
    }

    /// @brief Convert inertial trajectory samples into a constant frame.
    inline std::vector<TrajectorySample> trajectory_to_frame(const std::vector<TrajectorySample> &samples_in,
                                                             const RotatingFrame &frame)
    {
        std::vector<TrajectorySample> out;
        out.reserve(samples_in.size());
        for (const auto &s: samples_in)
        {
            out.push_back(inertial_sample_to_frame(s, frame));
        }
        return out;
    }

    /// @brief Convert trajectory samples expressed in a constant frame back to inertial coordinates.
    inline std::vector<TrajectorySample> trajectory_from_frame(const std::vector<TrajectorySample> &samples_in_frame,
                                                               const RotatingFrame &frame)
    {
        std::vector<TrajectorySample> out;
        out.reserve(samples_in_frame.size());
        for (const auto &s: samples_in_frame)
        {
            out.push_back(frame_sample_to_inertial(s, frame));
        }
        return out;
    }

    /// @brief Convert inertial trajectory samples into the time-varying body-centered inertial frame of a body (ECI/BCI).
    /// @note The input samples must be inertial (i.e., not already expressed in another frame).
    inline std::vector<TrajectorySample> trajectory_to_body_centered_inertial(const std::vector<TrajectorySample> &samples_in,
                                                                              const CelestialEphemeris &eph,
                                                                              const MassiveBody &body)
    {
        std::vector<TrajectorySample> out;
        out.reserve(samples_in.size());

        for (const auto &s: samples_in)
        {
            const RotatingFrame frame = make_body_centered_inertial_frame_at(eph, body, s.t_s);
            out.push_back(inertial_sample_to_frame(s, frame));
        }

        return out;
    }

    /// @brief Convert inertial trajectory samples into the time-varying body-fixed rotating frame of a body (ECEF).
    /// @note The input samples must be inertial (i.e., not already expressed in another frame).
    inline std::vector<TrajectorySample> trajectory_to_body_fixed(const std::vector<TrajectorySample> &samples_in,
                                                                  const CelestialEphemeris &eph,
                                                                  const MassiveBody &body)
    {
        std::vector<TrajectorySample> out;
        out.reserve(samples_in.size());

        for (const auto &s: samples_in)
        {
            const std::optional<RotatingFrame> frame = make_body_fixed_frame_at(eph, body, s.t_s);
            if (!frame.has_value())
            {
                return {};
            }
            out.push_back(inertial_sample_to_frame(s, *frame));
        }

        return out;
    }

    /// @brief Alias for trajectory_to_body_centered_inertial (commonly called ECI when the body is Earth).
    inline std::vector<TrajectorySample> trajectory_to_eci(const std::vector<TrajectorySample> &samples_in,
                                                           const CelestialEphemeris &eph,
                                                           const MassiveBody &body)
    {
        return trajectory_to_body_centered_inertial(samples_in, eph, body);
    }

    /// @brief Alias for trajectory_to_body_fixed (commonly called ECEF when the body is Earth).
    inline std::vector<TrajectorySample> trajectory_to_ecef(const std::vector<TrajectorySample> &samples_in,
                                                            const CelestialEphemeris &eph,
                                                            const MassiveBody &body)
    {
        return trajectory_to_body_fixed(samples_in, eph, body);
    }

    /// @brief Convert a single inertial trajectory sample into a synodic rotating frame.
    /// @note The input sample must be in the same inertial frame as the body states used to construct the SynodicFrame.
    inline TrajectorySample inertial_sample_to_synodic(const TrajectorySample &sample_in, const SynodicFrame &frame)
    {
        return inertial_sample_to_frame(sample_in, frame);
    }

    /// @brief Convert inertial trajectory samples into the time-varying synodic frame of bodies (A,B).
    /// @note The input samples must be inertial (i.e., not already expressed in another frame).
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

    /// @brief Convert inertial trajectory samples into a time-varying frame specified by TrajectoryFrameSpec.
    ///
    /// @note The input samples must be inertial (i.e., not already expressed in another frame).
    inline std::vector<TrajectorySample> trajectory_to_frame_spec(const std::vector<TrajectorySample> &samples_in,
                                                                  const CelestialEphemeris &eph,
                                                                  const std::vector<MassiveBody> &bodies,
                                                                  const TrajectoryFrameSpec &frame,
                                                                  const SpacecraftStateLookup &sc_lookup = nullptr)
    {
        if (samples_in.empty())
        {
            return samples_in;
        }

        const auto body_by_id = [&](const BodyId id) -> const MassiveBody * {
            if (id == kInvalidBodyId)
            {
                return nullptr;
            }
            for (const auto &b: bodies)
            {
                if (b.id == id)
                {
                    return &b;
                }
            }
            return nullptr;
        };

        switch (frame.type)
        {
        case TrajectoryFrameType::Inertial:
            return samples_in;
        case TrajectoryFrameType::BodyCenteredInertial:
        {
            const MassiveBody *body = body_by_id(frame.primary_body_id);
            if (body == nullptr)
            {
                return {};
            }
            return trajectory_to_body_centered_inertial(samples_in, eph, *body);
        }
        case TrajectoryFrameType::BodyFixed:
        {
            const MassiveBody *body = body_by_id(frame.primary_body_id);
            if (body == nullptr)
            {
                return {};
            }
            return trajectory_to_body_fixed(samples_in, eph, *body);
        }
        case TrajectoryFrameType::Synodic:
        {
            const MassiveBody *body_a = body_by_id(frame.primary_body_id);
            const MassiveBody *body_b = body_by_id(frame.secondary_body_id);
            if (body_a == nullptr || body_b == nullptr)
            {
                return {};
            }
            return trajectory_to_synodic(samples_in, eph, *body_a, *body_b);
        }
        case TrajectoryFrameType::LVLH:
        {
            if (!sc_lookup)
            {
                return {};
            }
            if (frame.target_spacecraft_id == kInvalidSpacecraftId)
            {
                return {};
            }
            if (bodies.empty())
            {
                return {};
            }

            const MassiveBody *primary_fixed = body_by_id(frame.primary_body_id);

            std::vector<std::optional<std::size_t>> eph_indices;
            if (!eph.empty())
            {
                eph_indices.reserve(bodies.size());
                for (const auto &b: bodies)
                {
                    std::size_t idx = 0;
                    if (eph.body_index_for_id(b.id, &idx))
                    {
                        eph_indices.push_back(idx);
                    }
                    else
                    {
                        eph_indices.push_back(std::nullopt);
                    }
                }
            }

            std::vector<TrajectorySample> out;
            out.reserve(samples_in.size());

            for (const auto &s: samples_in)
            {
                const std::optional<State> target_state = sc_lookup(frame.target_spacecraft_id, s.t_s);
                if (!target_state.has_value())
                {
                    return {};
                }

                const MassiveBody *primary = primary_fixed;
                std::optional<std::size_t> primary_eph_index;

                if (primary == nullptr)
                {
                    const std::size_t best = auto_select_primary_index(
                            bodies, target_state->position_m,
                            [&](std::size_t i) -> Vec3 {
                                if (!eph.empty() && i < eph_indices.size() && eph_indices[i].has_value())
                                    return eph.body_position_at(*eph_indices[i], s.t_s);
                                return bodies[i].state.position_m;
                            });

                    primary = &bodies[best];
                    if (!eph.empty() && best < eph_indices.size())
                    {
                        primary_eph_index = eph_indices[best];
                    }
                }
                else if (!eph.empty())
                {
                    for (std::size_t i = 0; i < bodies.size(); ++i)
                    {
                        if (bodies[i].id == primary->id)
                        {
                            if (i < eph_indices.size())
                            {
                                primary_eph_index = eph_indices[i];
                            }
                            break;
                        }
                    }
                }

                const State primary_state = [&]() -> State {
                    if (eph.empty() || !primary_eph_index.has_value())
                    {
                        return primary->state;
                    }
                    return eph.body_state_at(*primary_eph_index, s.t_s);
                }();
                const std::optional<RotatingFrame> lvlh = make_lvlh_frame(primary_state, *target_state);
                if (!lvlh.has_value())
                {
                    return {};
                }
                out.push_back(inertial_sample_to_frame(s, *lvlh));
            }

            return out;
        }
        }

        return samples_in;
    }

    inline std::vector<TrajectorySegment> transform_trajectory_segments_to_frame_spec(
            const std::vector<TrajectorySegment> &inertial_segments,
            const CelestialEphemeris &eph,
            const std::vector<MassiveBody> &bodies,
            const TrajectoryFrameSpec &frame,
            const FrameSegmentTransformOptions &opt,
            const SpacecraftStateLookup &sc_lookup = nullptr,
            FrameSegmentTransformDiagnostics *out_diag = nullptr)
    {
        detail::reset_adaptive_diagnostics_(out_diag);

        std::vector<TrajectorySegment> out;
        if (inertial_segments.empty())
        {
            return out;
        }

        if (frame.type == TrajectoryFrameType::Inertial)
        {
            out = inertial_segments;
            for (const auto &segment : out)
            {
                detail::record_adaptive_accept_(out_diag, segment.dt_s);
            }
            return out;
        }

        const FrameSegmentTransformOptions clean_opt = detail::sanitize_frame_segment_transform_options_(opt);
        const double base_t0_s = inertial_segments.front().t0_s;
        const double base_t1_s = inertial_segments.back().t0_s + inertial_segments.back().dt_s;
        const double base_duration_s = std::max(0.0, base_t1_s - base_t0_s);

        out.reserve(inertial_segments.size());
        for (const TrajectorySegment &segment : inertial_segments)
        {
            if (!detail::append_transformed_segment_(
                        segment,
                        eph,
                        bodies,
                        frame,
                        clean_opt,
                        base_t0_s,
                        base_duration_s,
                        sc_lookup,
                        segment.flags,
                        kTrajectorySegmentFlagNone,
                        &out,
                        out_diag))
            {
                break;
            }
        }

        return out;
    }

} // namespace orbitsim
