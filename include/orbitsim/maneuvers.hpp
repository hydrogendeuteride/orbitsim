#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>
#include <vector>

namespace orbitsim
{

    inline constexpr SpacecraftId kAllSpacecraft = std::numeric_limits<SpacecraftId>::max();

    // -------------------------------------------------------------------------
    // RTN direction constants for common burn orientations
    // -------------------------------------------------------------------------

    inline constexpr Vec3 kPrograde{0.0, 1.0, 0.0};   // +T: along velocity
    inline constexpr Vec3 kRetrograde{0.0, -1.0, 0.0}; // -T: against velocity
    inline constexpr Vec3 kRadialOut{1.0, 0.0, 0.0};  // +R: away from primary
    inline constexpr Vec3 kRadialIn{-1.0, 0.0, 0.0};  // -R: toward primary
    inline constexpr Vec3 kNormal{0.0, 0.0, 1.0};     // +N: normal to orbital plane
    inline constexpr Vec3 kAntiNormal{0.0, 0.0, -1.0}; // -N: opposite normal

    struct BurnSegment
    {
        double t_start_s{0.0};
        double t_end_s{0.0};
        BodyId primary_body_id{kInvalidBodyId}; // Burn RTN frame primary; invalid = auto-select.
        Vec3 dir_rtn_unit{0.0, 0.0, 0.0}; // Components in (R, T, N).
        double throttle_0_1{0.0};
        std::size_t engine_index{0};
        SpacecraftId spacecraft_id{kAllSpacecraft}; // Which spacecraft this burn targets; kAllSpacecraft = all.
    };

    // Instantaneous delta-v applied at a single time (useful for maneuver nodes / debugging).
    // dv_rtn_mps is expressed in RTN components relative to the chosen primary at t_s.
    struct ImpulseSegment
    {
        double t_s{0.0};
        BodyId primary_body_id{kInvalidBodyId}; // RTN frame primary; invalid = auto-select.
        Vec3 dv_rtn_mps{0.0, 0.0, 0.0}; // Components in (R, T, N) [m/s].
        SpacecraftId spacecraft_id{kAllSpacecraft};
    };

    struct ManeuverPlan
    {
        std::vector<BurnSegment> segments{};
        std::vector<ImpulseSegment> impulses{};
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

    inline const BurnSegment *active_burn_at(const ManeuverPlan &plan, const SpacecraftId spacecraft_id, const double t_s)
    {
        return detail::active_burn_at_if_(plan, t_s, [&](const BurnSegment &seg) {
            return segment_applies_to_spacecraft(seg, spacecraft_id);
        });
    }

    inline double next_burn_boundary_after(const ManeuverPlan &plan, const SpacecraftId spacecraft_id, const double t_s,
                                          const double t_end_s)
    {
        return detail::next_burn_boundary_after_if_(plan, t_s, t_end_s, [&](const BurnSegment &seg) {
            return segment_applies_to_spacecraft(seg, spacecraft_id);
        });
    }

    namespace detail
    {
        template<class Pred>
        inline double next_impulse_time_after_if_(const ManeuverPlan &plan, const double t_s, const double t_end_s, Pred pred)
        {
            double bound = t_end_s;
            for (const auto &seg: plan.impulses)
            {
                if (!pred(seg))
                {
                    continue;
                }
                if (seg.t_s > t_s && seg.t_s < bound)
                {
                    bound = seg.t_s;
                }
            }
            return bound;
        }

        template<class EphemerisLike>
        inline Vec3 rtn_vector_to_inertial_(const EphemerisLike &eph, const std::size_t primary_index, const double t_s,
                                            const Vec3 &sc_pos_m, const Vec3 &sc_vel_mps, const Vec3 &v_rtn)
        {
            const Vec3 r_primary_m = eph.body_position_at(primary_index, t_s);
            const Vec3 v_primary_mps = eph.body_velocity_at(primary_index, t_s);

            const Vec3 r_rel = sc_pos_m - r_primary_m;
            const Vec3 v_rel = sc_vel_mps - v_primary_mps;
            const RtnFrame f = compute_rtn_frame(r_rel, v_rel);

            return v_rtn.x * f.R + v_rtn.y * f.T + v_rtn.z * f.N;
        }

        template<class EphemerisLike>
        inline Vec3 burn_dir_inertial_unit_(const EphemerisLike &eph, const std::size_t primary_index, const double t_s,
                                            const Vec3 &sc_pos_m, const Vec3 &sc_vel_mps, const Vec3 &dir_rtn_unit)
        {
            const Vec3 dir_i =
                    rtn_vector_to_inertial_(eph, primary_index, t_s, sc_pos_m, sc_vel_mps, dir_rtn_unit);
            return normalized_or(dir_i, Vec3{0.0, 0.0, 0.0});
        }
    } // namespace detail

    inline Vec3 rtn_vector_to_inertial(const CelestialEphemerisSegment &eph, const std::size_t primary_index,
                                       const double t_s, const Vec3 &sc_pos_m, const Vec3 &sc_vel_mps,
                                       const Vec3 &v_rtn)
    {
        return detail::rtn_vector_to_inertial_(eph, primary_index, t_s, sc_pos_m, sc_vel_mps, v_rtn);
    }

    inline Vec3 rtn_vector_to_inertial(const CelestialEphemeris &eph, const std::size_t primary_index, const double t_s,
                                       const Vec3 &sc_pos_m, const Vec3 &sc_vel_mps, const Vec3 &v_rtn)
    {
        return detail::rtn_vector_to_inertial_(eph, primary_index, t_s, sc_pos_m, sc_vel_mps, v_rtn);
    }

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

    inline double next_impulse_time_after(const ManeuverPlan &plan, const SpacecraftId spacecraft_id, const double t_s,
                                          const double t_end_s)
    {
        return detail::next_impulse_time_after_if_(plan, t_s, t_end_s, [&](const ImpulseSegment &seg) {
            return segment_applies_to_spacecraft(seg, spacecraft_id);
        });
    }

    // -------------------------------------------------------------------------
    // BurnBuilder: fluent interface for creating BurnSegment
    // -------------------------------------------------------------------------

    /// @brief Fluent builder for creating BurnSegment objects.
    /// @example
    ///   auto seg = burn().start(hours(2.0)).duration(minutes(20.0)).prograde().spacecraft(sc_id);
    class BurnBuilder
    {
    public:
        BurnBuilder() = default;

        /// @brief Set burn start time [s].
        BurnBuilder &start(const double t_start_s)
        {
            seg_.t_start_s = t_start_s;
            return *this;
        }

        /// @brief Set burn end time [s].
        BurnBuilder &end(const double t_end_s)
        {
            seg_.t_end_s = t_end_s;
            return *this;
        }

        /// @brief Set burn duration [s]. Computes end time from start + duration.
        BurnBuilder &duration(const double duration_s)
        {
            seg_.t_end_s = seg_.t_start_s + duration_s;
            return *this;
        }

        /// @brief Set the RTN direction unit vector.
        BurnBuilder &direction(const Vec3 &dir_rtn_unit)
        {
            seg_.dir_rtn_unit = dir_rtn_unit;
            return *this;
        }

        /// @brief Set direction to prograde (+T).
        BurnBuilder &prograde()
        {
            seg_.dir_rtn_unit = kPrograde;
            return *this;
        }

        /// @brief Set direction to retrograde (-T).
        BurnBuilder &retrograde()
        {
            seg_.dir_rtn_unit = kRetrograde;
            return *this;
        }

        /// @brief Set direction to radial out (+R).
        BurnBuilder &radial_out()
        {
            seg_.dir_rtn_unit = kRadialOut;
            return *this;
        }

        /// @brief Set direction to radial in (-R).
        BurnBuilder &radial_in()
        {
            seg_.dir_rtn_unit = kRadialIn;
            return *this;
        }

        /// @brief Set direction to normal (+N).
        BurnBuilder &normal()
        {
            seg_.dir_rtn_unit = kNormal;
            return *this;
        }

        /// @brief Set direction to anti-normal (-N).
        BurnBuilder &anti_normal()
        {
            seg_.dir_rtn_unit = kAntiNormal;
            return *this;
        }

        /// @brief Set throttle level [0, 1].
        BurnBuilder &throttle(const double throttle_0_1)
        {
            seg_.throttle_0_1 = throttle_0_1;
            return *this;
        }

        /// @brief Set throttle to full (1.0).
        BurnBuilder &full_throttle()
        {
            seg_.throttle_0_1 = 1.0;
            return *this;
        }

        /// @brief Set the primary body for RTN frame computation.
        BurnBuilder &primary(const BodyId primary_body_id)
        {
            seg_.primary_body_id = primary_body_id;
            return *this;
        }

        /// @brief Set the target spacecraft.
        BurnBuilder &spacecraft(const SpacecraftId spacecraft_id)
        {
            seg_.spacecraft_id = spacecraft_id;
            return *this;
        }

        /// @brief Set the engine index.
        BurnBuilder &engine(const std::size_t engine_index)
        {
            seg_.engine_index = engine_index;
            return *this;
        }

        /// @brief Implicit conversion to BurnSegment.
        operator BurnSegment() const { return seg_; }

        /// @brief Explicit conversion to BurnSegment.
        BurnSegment build() const { return seg_; }

    private:
        BurnSegment seg_{.throttle_0_1 = 1.0}; // Default full throttle.
    };

    /// @brief Start building a BurnSegment with the fluent interface.
    inline BurnBuilder burn() { return BurnBuilder{}; }

    // -------------------------------------------------------------------------
    // ImpulseBuilder: fluent interface for creating ImpulseSegment
    // -------------------------------------------------------------------------

    class ImpulseBuilder
    {
    public:
        ImpulseBuilder() = default;

        ImpulseBuilder &time(const double t_s)
        {
            seg_.t_s = t_s;
            return *this;
        }

        ImpulseBuilder &dv_rtn(const Vec3 &dv_rtn_mps)
        {
            seg_.dv_rtn_mps = dv_rtn_mps;
            return *this;
        }

        ImpulseBuilder &primary(const BodyId primary_body_id)
        {
            seg_.primary_body_id = primary_body_id;
            return *this;
        }

        ImpulseBuilder &spacecraft(const SpacecraftId spacecraft_id)
        {
            seg_.spacecraft_id = spacecraft_id;
            return *this;
        }

        ImpulseBuilder &prograde(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{0.0, dv_mps, 0.0};
            return *this;
        }

        ImpulseBuilder &retrograde(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{0.0, -dv_mps, 0.0};
            return *this;
        }

        ImpulseBuilder &radial_out(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{dv_mps, 0.0, 0.0};
            return *this;
        }

        ImpulseBuilder &radial_in(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{-dv_mps, 0.0, 0.0};
            return *this;
        }

        ImpulseBuilder &normal(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{0.0, 0.0, dv_mps};
            return *this;
        }

        ImpulseBuilder &anti_normal(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{0.0, 0.0, -dv_mps};
            return *this;
        }

        operator ImpulseSegment() const { return seg_; }
        ImpulseSegment build() const { return seg_; }

    private:
        ImpulseSegment seg_{};
    };

    inline ImpulseBuilder impulse() { return ImpulseBuilder{}; }

    // -------------------------------------------------------------------------
    // Convenience factory functions for common burn types
    // -------------------------------------------------------------------------

    /// @brief Create a prograde burn (velocity-increasing).
    /// @param t_start_s Start time [s].
    /// @param duration_s Burn duration [s].
    /// @param spacecraft_id Target spacecraft (default: all spacecraft).
    /// @param primary_body_id Primary body for RTN frame (default: auto-select).
    /// @param throttle_0_1 Throttle level [0, 1] (default: 1.0).
    /// @param engine_index Engine index (default: 0).
    inline BurnSegment prograde_burn(const double t_start_s, const double duration_s,
                                     const SpacecraftId spacecraft_id = kAllSpacecraft,
                                     const BodyId primary_body_id = kInvalidBodyId,
                                     const double throttle_0_1 = 1.0, const std::size_t engine_index = 0)
    {
        return BurnSegment{.t_start_s = t_start_s,
                           .t_end_s = t_start_s + duration_s,
                           .primary_body_id = primary_body_id,
                           .dir_rtn_unit = kPrograde,
                           .throttle_0_1 = throttle_0_1,
                           .engine_index = engine_index,
                           .spacecraft_id = spacecraft_id};
    }

    /// @brief Create a retrograde burn (velocity-decreasing).
    /// @param t_start_s Start time [s].
    /// @param duration_s Burn duration [s].
    /// @param spacecraft_id Target spacecraft (default: all spacecraft).
    /// @param primary_body_id Primary body for RTN frame (default: auto-select).
    /// @param throttle_0_1 Throttle level [0, 1] (default: 1.0).
    /// @param engine_index Engine index (default: 0).
    inline BurnSegment retrograde_burn(const double t_start_s, const double duration_s,
                                       const SpacecraftId spacecraft_id = kAllSpacecraft,
                                       const BodyId primary_body_id = kInvalidBodyId,
                                       const double throttle_0_1 = 1.0, const std::size_t engine_index = 0)
    {
        return BurnSegment{.t_start_s = t_start_s,
                           .t_end_s = t_start_s + duration_s,
                           .primary_body_id = primary_body_id,
                           .dir_rtn_unit = kRetrograde,
                           .throttle_0_1 = throttle_0_1,
                           .engine_index = engine_index,
                           .spacecraft_id = spacecraft_id};
    }

    /// @brief Create a normal burn (plane change, +N direction).
    /// @param t_start_s Start time [s].
    /// @param duration_s Burn duration [s].
    /// @param spacecraft_id Target spacecraft (default: all spacecraft).
    /// @param primary_body_id Primary body for RTN frame (default: auto-select).
    /// @param throttle_0_1 Throttle level [0, 1] (default: 1.0).
    /// @param engine_index Engine index (default: 0).
    inline BurnSegment normal_burn(const double t_start_s, const double duration_s,
                                   const SpacecraftId spacecraft_id = kAllSpacecraft,
                                   const BodyId primary_body_id = kInvalidBodyId, const double throttle_0_1 = 1.0,
                                   const std::size_t engine_index = 0)
    {
        return BurnSegment{.t_start_s = t_start_s,
                           .t_end_s = t_start_s + duration_s,
                           .primary_body_id = primary_body_id,
                           .dir_rtn_unit = kNormal,
                           .throttle_0_1 = throttle_0_1,
                           .engine_index = engine_index,
                           .spacecraft_id = spacecraft_id};
    }

} // namespace orbitsim
