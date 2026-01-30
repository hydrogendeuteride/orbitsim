#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

namespace orbitsim
{

    inline constexpr std::size_t kAllSpacecraft = std::numeric_limits<std::size_t>::max();

    struct BurnSegment
    {
        double t_start_s{0.0};
        double t_end_s{0.0};
        std::size_t primary_index{0};
        Vec3 dir_rtn_unit{0.0, 0.0, 0.0}; // Components in (R, T, N).
        double throttle_0_1{0.0};
        std::size_t engine_index{0};
        std::size_t spacecraft_index{kAllSpacecraft}; // Which spacecraft this burn targets; kAllSpacecraft = all.
    };

    struct ManeuverPlan
    {
        std::vector<BurnSegment> segments{};
    };

    inline bool segment_applies_to_spacecraft(const BurnSegment &seg, const std::size_t spacecraft_index)
    {
        return seg.spacecraft_index == kAllSpacecraft || seg.spacecraft_index == spacecraft_index;
    }

    inline void sort_segments_by_start(ManeuverPlan &plan)
    {
        std::sort(plan.segments.begin(), plan.segments.end(),
                  [](const BurnSegment &a, const BurnSegment &b) { return a.t_start_s < b.t_start_s; });
    }

    namespace detail
    {
        template<class Pred>
        inline const BurnSegment *active_burn_at_if_(const ManeuverPlan &plan, const double t_s, Pred pred)
        {
            for (const auto &seg: plan.segments)
            {
                if (!pred(seg))
                {
                    continue;
                }
                if (seg.t_start_s <= t_s && t_s < seg.t_end_s)
                {
                    return &seg;
                }
            }
            return nullptr;
        }

        template<class Pred>
        inline double next_burn_boundary_after_if_(const ManeuverPlan &plan, const double t_s, const double t_end_s,
                                                  Pred pred)
        {
            double bound = t_end_s;
            for (const auto &seg: plan.segments)
            {
                if (!pred(seg))
                {
                    continue;
                }
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
    } // namespace detail

    inline const BurnSegment *active_burn_at(const ManeuverPlan &plan, const double t_s)
    {
        return detail::active_burn_at_if_(plan, t_s, [](const BurnSegment &) { return true; });
    }

    inline const BurnSegment *active_burn_at(const ManeuverPlan &plan, const std::size_t spacecraft_index, const double t_s)
    {
        return detail::active_burn_at_if_(plan, t_s, [&](const BurnSegment &seg) {
            return segment_applies_to_spacecraft(seg, spacecraft_index);
        });
    }

    inline double next_burn_boundary_after(const ManeuverPlan &plan, const double t_s, const double t_end_s)
    {
        return detail::next_burn_boundary_after_if_(plan, t_s, t_end_s, [](const BurnSegment &) { return true; });
    }

    inline double next_burn_boundary_after(const ManeuverPlan &plan, const std::size_t spacecraft_index, const double t_s,
                                          const double t_end_s)
    {
        return detail::next_burn_boundary_after_if_(plan, t_s, t_end_s, [&](const BurnSegment &seg) {
            return segment_applies_to_spacecraft(seg, spacecraft_index);
        });
    }

    namespace detail
    {
        template<class EphemerisLike>
        inline Vec3 burn_dir_inertial_unit_(const EphemerisLike &eph, const std::size_t primary_index, const double t_s,
                                            const Vec3 &sc_pos_m, const Vec3 &sc_vel_mps, const Vec3 &dir_rtn_unit)
        {
            const Vec3 r_primary_m = eph.body_position_at(primary_index, t_s);
            const Vec3 v_primary_mps = eph.body_velocity_at(primary_index, t_s);

            const Vec3 r_rel = sc_pos_m - r_primary_m;
            const Vec3 v_rel = sc_vel_mps - v_primary_mps;
            const RtnFrame f = compute_rtn_frame(r_rel, v_rel);

            const Vec3 dir_i = dir_rtn_unit.x * f.R + dir_rtn_unit.y * f.T + dir_rtn_unit.z * f.N;
            return normalized_or(dir_i, Vec3{0.0, 0.0, 0.0});
        }
    } // namespace detail

    inline Vec3 burn_dir_inertial_unit(const CelestialEphemerisSegment &eph, const std::size_t primary_index,
                                       const double t_s, const Vec3 &sc_pos_m, const Vec3 &sc_vel_mps,
                                       const Vec3 &dir_rtn_unit)
    {
        return detail::burn_dir_inertial_unit_(eph, primary_index, t_s, sc_pos_m, sc_vel_mps, dir_rtn_unit);
    }

    inline Vec3 burn_dir_inertial_unit(const CelestialEphemeris &eph, const std::size_t primary_index, const double t_s,
                                       const Vec3 &sc_pos_m, const Vec3 &sc_vel_mps, const Vec3 &dir_rtn_unit)
    {
        return detail::burn_dir_inertial_unit_(eph, primary_index, t_s, sc_pos_m, sc_vel_mps, dir_rtn_unit);
    }

} // namespace orbitsim
