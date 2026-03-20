#pragma once

#include "orbitsim/types.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace orbitsim
{

    enum TrajectorySegmentFlagBits : std::uint32_t
    {
        kTrajectorySegmentFlagNone = 0u,
        kTrajectorySegmentFlagForcedBoundary = 1u << 0,
        kTrajectorySegmentFlagImpulseBoundary = 1u << 1,
        kTrajectorySegmentFlagBurnBoundary = 1u << 2,
        kTrajectorySegmentFlagImpactTerminal = 1u << 3,
        kTrajectorySegmentFlagFrameResegmented = 1u << 4,
    };

    struct TrajectorySegment
    {
        double t0_s{0.0};
        double dt_s{0.0};
        State start{};
        State end{};
        std::uint32_t flags{0};
    };

    struct TrajectorySegmentOptions
    {
        double duration_s{3600.0};
        std::size_t max_segments{2048};
        bool include_start{true};
        bool include_end{true};
        bool stop_on_impact{false};
        double lookup_dt_s{0.0};
    };

    struct AdaptiveToleranceRamp
    {
        double pos_near_m{1.0};
        double pos_far_m{1000.0};
        double vel_near_mps{1.0e-3};
        double vel_far_mps{1.0};
        double rel_pos_floor{1.0e-9};
        double rel_vel_floor{1.0e-9};
    };

    struct AdaptiveSegmentOptions
    {
        double duration_s{3600.0};
        double min_dt_s{1.0};
        double max_dt_s{60.0};
        double lookup_max_dt_s{0.0};
        std::size_t soft_max_segments{2048};
        std::size_t hard_max_segments{8192};
        AdaptiveToleranceRamp tolerance{};
        bool stop_on_impact{false};
        bool split_at_impulses{true};
        bool split_at_burn_boundaries{true};
        std::vector<double> forced_split_times_s{};
        std::function<bool()> cancel_requested{};
    };

    struct AdaptiveEphemerisOptions
    {
        double duration_s{3600.0};
        double min_dt_s{1.0};
        double max_dt_s{60.0};
        std::size_t soft_max_segments{2048};
        std::size_t hard_max_segments{8192};
        AdaptiveToleranceRamp tolerance{};
        std::function<bool()> cancel_requested{};
    };

    struct FrameSegmentTransformOptions
    {
        double min_dt_s{1.0};
        double max_dt_s{60.0};
        std::size_t soft_max_segments{2048};
        std::size_t hard_max_segments{8192};
        AdaptiveToleranceRamp tolerance{};
        std::vector<double> forced_split_times_s{};
        std::function<bool()> cancel_requested{};
    };

    struct AdaptiveSegmentDiagnostics
    {
        double covered_duration_s{0.0};
        std::size_t accepted_segments{0};
        std::size_t rejected_splits{0};
        std::size_t forced_boundary_splits{0};
        double min_dt_s{0.0};
        double max_dt_s{0.0};
        double avg_dt_s{0.0};
        bool hard_cap_hit{false};
        bool cancelled{false};
    };

    struct AdaptiveEphemerisDiagnostics
    {
        double covered_duration_s{0.0};
        std::size_t accepted_segments{0};
        std::size_t rejected_splits{0};
        std::size_t forced_boundary_splits{0};
        double min_dt_s{0.0};
        double max_dt_s{0.0};
        double avg_dt_s{0.0};
        bool hard_cap_hit{false};
        bool cancelled{false};
    };

    struct FrameSegmentTransformDiagnostics
    {
        double covered_duration_s{0.0};
        std::size_t accepted_segments{0};
        std::size_t rejected_splits{0};
        std::size_t forced_boundary_splits{0};
        std::size_t frame_resegmentation_count{0};
        double min_dt_s{0.0};
        double max_dt_s{0.0};
        double avg_dt_s{0.0};
        bool hard_cap_hit{false};
        bool cancelled{false};
    };

} // namespace orbitsim
