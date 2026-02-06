#pragma once

#include "orbitsim/coordinate_frames.hpp"
#include "orbitsim/ephemeris.hpp"
#include "orbitsim/maneuvers_builders.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/spacecraft_lookup.hpp"
#include "orbitsim/synodic.hpp"

#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

namespace orbitsim
{

    namespace detail
    {
        inline bool body_index_for_id_(const std::vector<MassiveBody> &bodies, const BodyId id, std::size_t *out_index)
        {
            if (out_index == nullptr || id == kInvalidBodyId)
            {
                return false;
            }
            for (std::size_t i = 0; i < bodies.size(); ++i)
            {
                if (bodies[i].id == id)
                {
                    *out_index = i;
                    return true;
                }
            }
            return false;
        }

        template<class EphemerisLike>
        inline std::optional<RotatingFrame>
        make_rtn_reference_frame_at_(const EphemerisLike &eph, const std::vector<MassiveBody> &bodies,
                                     const TrajectoryFrameSpec &spec, const double t_s)
        {
            switch (spec.type)
            {
                case TrajectoryFrameType::Inertial:
                    return RotatingFrame{};
                case TrajectoryFrameType::BodyCenteredInertial:
                {
                    std::size_t idx = 0;
                    if (!body_index_for_id_(bodies, spec.primary_body_id, &idx))
                    {
                        return std::nullopt;
                    }
                    return make_body_centered_inertial_frame(eph.body_state_at(idx, t_s));
                }
                case TrajectoryFrameType::BodyFixed:
                {
                    std::size_t idx = 0;
                    if (!body_index_for_id_(bodies, spec.primary_body_id, &idx))
                    {
                        return std::nullopt;
                    }
                    return make_body_fixed_frame(eph.body_state_at(idx, t_s));
                }
                case TrajectoryFrameType::Synodic:
                {
                    std::size_t ia = 0;
                    std::size_t ib = 0;
                    if (!body_index_for_id_(bodies, spec.primary_body_id, &ia) ||
                        !body_index_for_id_(bodies, spec.secondary_body_id, &ib))
                    {
                        return std::nullopt;
                    }

                    const State a_state = eph.body_state_at(ia, t_s);
                    const State b_state = eph.body_state_at(ib, t_s);
                    const std::optional<SynodicFrame> syn =
                            make_synodic_frame(a_state, bodies[ia].mass_kg, b_state, bodies[ib].mass_kg);
                    if (!syn.has_value())
                    {
                        return std::nullopt;
                    }
                    return RotatingFrame{*syn};
                }
                case TrajectoryFrameType::LVLH:
                    return RotatingFrame{};
            }

            return std::nullopt;
        }

        template<class EphemerisLike>
        inline Vec3 rtn_vector_to_inertial_in_frame_(const EphemerisLike &eph, const std::vector<MassiveBody> &bodies,
                                                     const std::size_t primary_index, const double t_s,
                                                     const Vec3 &sc_pos_m, const Vec3 &sc_vel_mps, const Vec3 &v_rtn,
                                                     const TrajectoryFrameSpec &rtn_frame,
                                                     const SpacecraftStateLookup &sc_lookup = nullptr)
        {
            // For LVLH frame, use target spacecraft's state to compute RTN basis.
            Vec3 rtn_basis_pos_m = sc_pos_m;
            Vec3 rtn_basis_vel_mps = sc_vel_mps;
            if (rtn_frame.type == TrajectoryFrameType::LVLH)
            {
                if (!sc_lookup)
                {
                    return Vec3{0.0, 0.0, 0.0};
                }
                const std::optional<State> target_state = sc_lookup(rtn_frame.target_spacecraft_id, t_s);
                if (!target_state.has_value())
                {
                    return Vec3{0.0, 0.0, 0.0};
                }
                rtn_basis_pos_m = target_state->position_m;
                rtn_basis_vel_mps = target_state->velocity_mps;
            }

            RotatingFrame ref_frame{};
            if (rtn_frame.type != TrajectoryFrameType::LVLH)
            {
                const std::optional<RotatingFrame> ref_frame_opt =
                        make_rtn_reference_frame_at_(eph, bodies, rtn_frame, t_s);
                if (!ref_frame_opt.has_value())
                {
                    return Vec3{0.0, 0.0, 0.0};
                }
                ref_frame = *ref_frame_opt;
            }

            const State sc_f = inertial_state_to_frame(make_state(rtn_basis_pos_m, rtn_basis_vel_mps), ref_frame);
            const State primary_f = inertial_state_to_frame(eph.body_state_at(primary_index, t_s), ref_frame);

            const Vec3 r_rel_f = sc_f.position_m - primary_f.position_m;
            const Vec3 v_rel_f = sc_f.velocity_mps - primary_f.velocity_mps;
            const RtnFrame f = compute_rtn_frame(r_rel_f, v_rel_f);

            const Vec3 v_frame = v_rtn.x * f.R + v_rtn.y * f.T + v_rtn.z * f.N;
            return frame_vector_to_inertial(ref_frame, v_frame);
        }

        template<class EphemerisLike>
        inline Vec3 burn_dir_inertial_unit_in_frame_(const EphemerisLike &eph, const std::vector<MassiveBody> &bodies,
                                                     const std::size_t primary_index, const double t_s,
                                                     const Vec3 &sc_pos_m, const Vec3 &sc_vel_mps,
                                                     const Vec3 &dir_rtn_unit, const TrajectoryFrameSpec &rtn_frame,
                                                     const SpacecraftStateLookup &sc_lookup = nullptr)
        {
            const Vec3 dir_i = rtn_vector_to_inertial_in_frame_(eph, bodies, primary_index, t_s, sc_pos_m, sc_vel_mps,
                                                                dir_rtn_unit, rtn_frame, sc_lookup);
            return normalized_or(dir_i, Vec3{0.0, 0.0, 0.0});
        }

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

                const bool seg_active = (seg.t_start_s <= t_s && t_s < seg.t_end_s);
                if (seg_active && seg.t_end_s > t_s && seg.t_end_s < bound)
                {
                    bound = seg.t_end_s;
                }
            }
            return bound;
        }
    } // namespace detail

    /** @brief Find active burn segment at time t_s for a spacecraft, or nullptr if none. */
    inline const BurnSegment *active_burn_at(const ManeuverPlan &plan, const SpacecraftId spacecraft_id,
                                             const double t_s)
    {
        return detail::active_burn_at_if_(
                plan, t_s, [&](const BurnSegment &seg) { return segment_applies_to_spacecraft(seg, spacecraft_id); });
    }

    /** @brief Find next burn start/end time after t_s, or t_end_s if none. */
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
        inline double next_impulse_time_after_if_(const ManeuverPlan &plan, const double t_s, const double t_end_s,
                                                  Pred pred)
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
    } // namespace detail

    /** @brief Convert vector from RTN (computed in a chosen reference frame) to inertial coordinates.
     *  @param sc_lookup Optional callback for LVLH frame to look up target spacecraft state.
     */
    inline Vec3 rtn_vector_to_inertial(const CelestialEphemerisSegment &eph, const std::vector<MassiveBody> &bodies,
                                       const std::size_t primary_index, const double t_s, const Vec3 &sc_pos_m,
                                       const Vec3 &sc_vel_mps, const Vec3 &v_rtn,
                                       const TrajectoryFrameSpec &rtn_frame = {},
                                       const SpacecraftStateLookup &sc_lookup = nullptr)
    {
        return detail::rtn_vector_to_inertial_in_frame_(eph, bodies, primary_index, t_s, sc_pos_m, sc_vel_mps, v_rtn,
                                                        rtn_frame, sc_lookup);
    }

    inline Vec3 rtn_vector_to_inertial(const CelestialEphemeris &eph, const std::vector<MassiveBody> &bodies,
                                       const std::size_t primary_index, const double t_s, const Vec3 &sc_pos_m,
                                       const Vec3 &sc_vel_mps, const Vec3 &v_rtn,
                                       const TrajectoryFrameSpec &rtn_frame = {},
                                       const SpacecraftStateLookup &sc_lookup = nullptr)
    {
        return detail::rtn_vector_to_inertial_in_frame_(eph, bodies, primary_index, t_s, sc_pos_m, sc_vel_mps, v_rtn,
                                                        rtn_frame, sc_lookup);
    }

    /** @brief Convert RTN direction (computed in a chosen reference frame) to normalized inertial direction.
     *  @param sc_lookup Optional callback for LVLH frame to look up target spacecraft state.
     */
    inline Vec3 burn_dir_inertial_unit(const CelestialEphemerisSegment &eph, const std::vector<MassiveBody> &bodies,
                                       const std::size_t primary_index, const double t_s, const Vec3 &sc_pos_m,
                                       const Vec3 &sc_vel_mps, const Vec3 &dir_rtn_unit,
                                       const TrajectoryFrameSpec &rtn_frame = {},
                                       const SpacecraftStateLookup &sc_lookup = nullptr)
    {
        return detail::burn_dir_inertial_unit_in_frame_(eph, bodies, primary_index, t_s, sc_pos_m, sc_vel_mps,
                                                        dir_rtn_unit, rtn_frame, sc_lookup);
    }

    inline Vec3 burn_dir_inertial_unit(const CelestialEphemeris &eph, const std::vector<MassiveBody> &bodies,
                                       const std::size_t primary_index, const double t_s, const Vec3 &sc_pos_m,
                                       const Vec3 &sc_vel_mps, const Vec3 &dir_rtn_unit,
                                       const TrajectoryFrameSpec &rtn_frame = {},
                                       const SpacecraftStateLookup &sc_lookup = nullptr)
    {
        return detail::burn_dir_inertial_unit_in_frame_(eph, bodies, primary_index, t_s, sc_pos_m, sc_vel_mps,
                                                        dir_rtn_unit, rtn_frame, sc_lookup);
    }

    /** @brief Find next impulse time after t_s, or t_end_s if none. */
    inline double next_impulse_time_after(const ManeuverPlan &plan, const SpacecraftId spacecraft_id, const double t_s,
                                          const double t_end_s)
    {
        return detail::next_impulse_time_after_if_(plan, t_s, t_end_s, [&](const ImpulseSegment &seg) {
            return segment_applies_to_spacecraft(seg, spacecraft_id);
        });
    }

} // namespace orbitsim
