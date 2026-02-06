#pragma once

#include "orbitsim/maneuvers_types.hpp"
#include "orbitsim/math.hpp"

#include <cstddef>

namespace orbitsim
{

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

        /// @brief Set the reference frame used to compute the RTN basis.
        BurnBuilder &rtn_frame(const TrajectoryFrameSpec &frame)
        {
            seg_.rtn_frame = frame;
            return *this;
        }

        /// @brief Compute RTN basis in the simulation inertial frame.
        BurnBuilder &rtn_inertial()
        {
            seg_.rtn_frame = TrajectoryFrameSpec::inertial();
            return *this;
        }

        /// @brief Compute RTN basis in a body-centered inertial frame (translation only; no rotation).
        BurnBuilder &rtn_body_centered_inertial(const BodyId body_id)
        {
            seg_.rtn_frame = TrajectoryFrameSpec::body_centered_inertial(body_id);
            return *this;
        }

        /// @brief Compute RTN basis in a body-fixed rotating frame (ECEF-style).
        BurnBuilder &rtn_body_fixed(const BodyId body_id)
        {
            seg_.rtn_frame = TrajectoryFrameSpec::body_fixed(body_id);
            return *this;
        }

        /// @brief Compute RTN basis in a two-body synodic rotating frame derived from (A,B).
        BurnBuilder &rtn_synodic(const BodyId body_a_id, const BodyId body_b_id)
        {
            seg_.rtn_frame = TrajectoryFrameSpec::synodic(body_a_id, body_b_id);
            return *this;
        }

        /// @brief Compute RTN basis in LVLH frame of a target spacecraft.
        /// @param target_sc_id Target spacecraft whose LVLH frame is used.
        /// @param primary_body_id Primary body for RTN computation (default: auto-select).
        BurnBuilder &rtn_lvlh(const SpacecraftId target_sc_id, const BodyId primary_body_id = kInvalidBodyId)
        {
            seg_.rtn_frame = TrajectoryFrameSpec::lvlh(target_sc_id, primary_body_id);
            if (primary_body_id != kInvalidBodyId)
            {
                seg_.primary_body_id = primary_body_id;
            }
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

    /// @brief Fluent builder for creating ImpulseSegment objects.
    /// @example
    ///   auto imp = impulse().time(hours(1.0)).prograde(100.0).spacecraft(sc_id);
    class ImpulseBuilder
    {
    public:
        ImpulseBuilder() = default;

        /// @brief Set impulse time [s].
        ImpulseBuilder &time(const double t_s)
        {
            seg_.t_s = t_s;
            return *this;
        }

        /// @brief Set delta-v in RTN components [m/s].
        ImpulseBuilder &dv_rtn(const Vec3 &dv_rtn_mps)
        {
            seg_.dv_rtn_mps = dv_rtn_mps;
            return *this;
        }

        /// @brief Set delta-v direction (unit vector in RTN) and magnitude [m/s].
        ImpulseBuilder &direction(const Vec3 &dir_rtn_unit, const double magnitude_mps)
        {
            seg_.dv_rtn_mps = normalized_or(dir_rtn_unit, Vec3{0.0, 1.0, 0.0}) * magnitude_mps;
            return *this;
        }

        /// @brief Set the primary body for RTN frame computation.
        ImpulseBuilder &primary(const BodyId primary_body_id)
        {
            seg_.primary_body_id = primary_body_id;
            return *this;
        }

        /// @brief Set the reference frame used to compute the RTN basis.
        ImpulseBuilder &rtn_frame(const TrajectoryFrameSpec &frame)
        {
            seg_.rtn_frame = frame;
            return *this;
        }

        /// @brief Compute RTN basis in the simulation inertial frame.
        ImpulseBuilder &rtn_inertial()
        {
            seg_.rtn_frame = TrajectoryFrameSpec::inertial();
            return *this;
        }

        /// @brief Compute RTN basis in a body-centered inertial frame (translation only; no rotation).
        ImpulseBuilder &rtn_body_centered_inertial(const BodyId body_id)
        {
            seg_.rtn_frame = TrajectoryFrameSpec::body_centered_inertial(body_id);
            return *this;
        }

        /// @brief Compute RTN basis in a body-fixed rotating frame (ECEF-style).
        ImpulseBuilder &rtn_body_fixed(const BodyId body_id)
        {
            seg_.rtn_frame = TrajectoryFrameSpec::body_fixed(body_id);
            return *this;
        }

        /// @brief Compute RTN basis in a two-body synodic rotating frame derived from (A,B).
        ImpulseBuilder &rtn_synodic(const BodyId body_a_id, const BodyId body_b_id)
        {
            seg_.rtn_frame = TrajectoryFrameSpec::synodic(body_a_id, body_b_id);
            return *this;
        }

        /// @brief Compute RTN basis in LVLH frame of a target spacecraft.
        /// @param target_sc_id Target spacecraft whose LVLH frame is used.
        /// @param primary_body_id Primary body for RTN computation (default: auto-select).
        ImpulseBuilder &rtn_lvlh(const SpacecraftId target_sc_id, const BodyId primary_body_id = kInvalidBodyId)
        {
            seg_.rtn_frame = TrajectoryFrameSpec::lvlh(target_sc_id, primary_body_id);
            if (primary_body_id != kInvalidBodyId)
            {
                seg_.primary_body_id = primary_body_id;
            }
            return *this;
        }

        /// @brief Set the target spacecraft.
        ImpulseBuilder &spacecraft(const SpacecraftId spacecraft_id)
        {
            seg_.spacecraft_id = spacecraft_id;
            return *this;
        }

        /// @brief Set prograde delta-v (+T direction).
        ImpulseBuilder &prograde(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{0.0, dv_mps, 0.0};
            return *this;
        }

        /// @brief Set retrograde delta-v (-T direction).
        ImpulseBuilder &retrograde(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{0.0, -dv_mps, 0.0};
            return *this;
        }

        /// @brief Set radial-out delta-v (+R direction).
        ImpulseBuilder &radial_out(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{dv_mps, 0.0, 0.0};
            return *this;
        }

        /// @brief Set radial-in delta-v (-R direction).
        ImpulseBuilder &radial_in(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{-dv_mps, 0.0, 0.0};
            return *this;
        }

        /// @brief Set normal delta-v (+N direction).
        ImpulseBuilder &normal(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{0.0, 0.0, dv_mps};
            return *this;
        }

        /// @brief Set anti-normal delta-v (-N direction).
        ImpulseBuilder &anti_normal(const double dv_mps)
        {
            seg_.dv_rtn_mps = Vec3{0.0, 0.0, -dv_mps};
            return *this;
        }

        /// @brief Implicit conversion to ImpulseSegment.
        operator ImpulseSegment() const { return seg_; }

        /// @brief Explicit conversion to ImpulseSegment.
        ImpulseSegment build() const { return seg_; }

    private:
        ImpulseSegment seg_{};
    };

    /// @brief Start building an ImpulseSegment with the fluent interface.
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
                                     const BodyId primary_body_id = kInvalidBodyId, const double throttle_0_1 = 1.0,
                                     const std::size_t engine_index = 0)
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
                                       const BodyId primary_body_id = kInvalidBodyId, const double throttle_0_1 = 1.0,
                                       const std::size_t engine_index = 0)
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
