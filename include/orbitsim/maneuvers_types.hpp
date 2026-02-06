#pragma once

#include "orbitsim/frame_spec.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

namespace orbitsim
{

    inline constexpr SpacecraftId kAllSpacecraft = std::numeric_limits<SpacecraftId>::max();

    // -------------------------------------------------------------------------
    // RTN direction constants for common burn orientations
    // -------------------------------------------------------------------------

    inline constexpr Vec3 kPrograde{0.0, 1.0, 0.0}; // +T: along velocity
    inline constexpr Vec3 kRetrograde{0.0, -1.0, 0.0}; // -T: against velocity
    inline constexpr Vec3 kRadialOut{1.0, 0.0, 0.0}; // +R: away from primary
    inline constexpr Vec3 kRadialIn{-1.0, 0.0, 0.0}; // -R: toward primary
    inline constexpr Vec3 kNormal{0.0, 0.0, 1.0}; // +N: normal to orbital plane
    inline constexpr Vec3 kAntiNormal{0.0, 0.0, -1.0}; // -N: opposite normal

    /**
     * @brief Continuous thrust segment over a time interval.
     *
     * Thrust direction is specified in RTN (Radial-Tangential-Normal) frame
     * relative to a primary body, with the RTN basis computed in `rtn_frame`.
     *
     * If primary_body_id is invalid, the body with highest gravitational
     * acceleration is auto-selected.
     */
    struct BurnSegment
    {
        double t_start_s{0.0};
        double t_end_s{0.0};
        BodyId primary_body_id{kInvalidBodyId}; ///< RTN frame primary; invalid = auto-select
        TrajectoryFrameSpec rtn_frame{}; ///< Frame used to compute RTN basis (default: inertial)
        Vec3 dir_rtn_unit{0.0, 0.0, 0.0}; ///< Thrust direction in (R, T, N)
        double throttle_0_1{0.0}; ///< Throttle level [0, 1]
        std::size_t engine_index{0}; ///< Index into Spacecraft::engines
        SpacecraftId spacecraft_id{kAllSpacecraft}; ///< Target spacecraft; kAllSpacecraft = all
    };

    /**
     * @brief Instantaneous delta-v at a single time (ideal maneuver node).
     *
     * Useful for maneuver planning and debugging. The delta-v is expressed
     * in RTN components relative to the chosen primary body, with the RTN
     * basis computed in `rtn_frame`.
     */
    struct ImpulseSegment
    {
        double t_s{0.0};
        BodyId primary_body_id{kInvalidBodyId}; ///< RTN frame primary; invalid = auto-select
        TrajectoryFrameSpec rtn_frame{}; ///< Frame used to compute RTN basis (default: inertial)
        Vec3 dv_rtn_mps{0.0, 0.0, 0.0}; ///< Delta-v in (R, T, N) [m/s]
        SpacecraftId spacecraft_id{kAllSpacecraft};
    };

    /** @brief Collection of scheduled burns and impulses for spacecraft. */
    struct ManeuverPlan
    {
        std::vector<BurnSegment> segments{}; ///< Continuous thrust burns
        std::vector<ImpulseSegment> impulses{}; ///< Instantaneous delta-v maneuvers
    };

    inline bool segment_applies_to_spacecraft(const BurnSegment &seg, const SpacecraftId spacecraft_id)
    {
        return seg.spacecraft_id == kAllSpacecraft || seg.spacecraft_id == spacecraft_id;
    }

    inline bool segment_applies_to_spacecraft(const ImpulseSegment &seg, const SpacecraftId spacecraft_id)
    {
        return seg.spacecraft_id == kAllSpacecraft || seg.spacecraft_id == spacecraft_id;
    }

    inline void sort_segments_by_start(ManeuverPlan &plan)
    {
        std::sort(plan.segments.begin(), plan.segments.end(),
                  [](const BurnSegment &a, const BurnSegment &b) { return a.t_start_s < b.t_start_s; });
    }

    inline void sort_impulses_by_time(ManeuverPlan &plan)
    {
        std::sort(plan.impulses.begin(), plan.impulses.end(),
                  [](const ImpulseSegment &a, const ImpulseSegment &b) { return a.t_s < b.t_s; });
    }

    inline void sort_plan(ManeuverPlan &plan)
    {
        sort_segments_by_start(plan);
        sort_impulses_by_time(plan);
    }

} // namespace orbitsim
