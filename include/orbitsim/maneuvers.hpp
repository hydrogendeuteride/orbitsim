#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace orbitsim
{

    struct BurnSegment
    {
        double t_start_s{0.0};
        double t_end_s{0.0};
        std::size_t primary_index{0};
        Vec3 dir_rtn_unit{0.0, 0.0, 0.0}; // Components in (R, T, N).
        double throttle_0_1{0.0};
        std::size_t engine_index{0};
    };

    struct ManeuverPlan
    {
        std::vector<BurnSegment> segments{};
    };

    inline void sort_segments_by_start(ManeuverPlan &plan)
    {
        std::sort(plan.segments.begin(), plan.segments.end(),
                  [](const BurnSegment &a, const BurnSegment &b) { return a.t_start_s < b.t_start_s; });
    }

    inline const BurnSegment *active_burn_at(const ManeuverPlan &plan, const double t_s)
    {
        for (const auto &seg: plan.segments)
        {
            if (seg.t_start_s <= t_s && t_s < seg.t_end_s)
            {
                return &seg;
            }
        }
        return nullptr;
    }

    inline double next_burn_boundary_after(const ManeuverPlan &plan, const double t_s, const double t_end_s)
    {
        double bound = t_end_s;
        for (const auto &seg: plan.segments)
        {
            if (seg.t_start_s > t_s && seg.t_start_s < bound)
            {
                bound = seg.t_start_s;
            }
            if (seg.t_end_s > t_s && seg.t_end_s < bound && seg.t_start_s <= t_s && t_s < seg.t_end_s)
            {
                bound = seg.t_end_s;
            }
        }
        return bound;
    }

    inline Vec3 burn_dir_inertial_unit(const CelestialEphemerisSegment &eph, const std::size_t primary_index,
                                       const double t_s, const Vec3 &sc_pos_m, const Vec3 &sc_vel_mps,
                                       const Vec3 &dir_rtn_unit)
    {
        const Vec3 r_primary_m = eph.body_position_at(primary_index, t_s);
        const Vec3 v_primary_mps = eph.body_velocity_at(primary_index, t_s);

        const Vec3 r_rel = sc_pos_m - r_primary_m;
        const Vec3 v_rel = sc_vel_mps - v_primary_mps;
        const RtnFrame f = compute_rtn_frame(r_rel, v_rel);

        const Vec3 dir_i = dir_rtn_unit.x * f.R + dir_rtn_unit.y * f.T + dir_rtn_unit.z * f.N;
        return normalized_or(dir_i, Vec3{0.0, 0.0, 0.0});
    }

} // namespace orbitsim
