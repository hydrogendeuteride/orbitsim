#pragma once

#include "orbitsim/detail/spacecraft_propagation.hpp"
#include "orbitsim/spacecraft_lookup.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace orbitsim
{

    /**
     * @brief On-demand time-indexed spacecraft state cache for LVLH/relative-frame features.
     *
     * Builds piecewise trajectory segments for queried spacecraft IDs by propagating them forward through the same
     * gravity + maneuver model as the main simulation/predictors, then answers state queries via Hermite interpolation.
     *
     * This is primarily intended to provide time-aware target spacecraft states for LVLH maneuvers and plotting
     * without re-integrating the target trajectory from scratch on every lookup.
     *
     * @tparam EphemerisLike Type providing `body_state_at(index, t_s)` and `body_position_at(index, t_s)`.
     */
    template<class EphemerisLike>
    class SpacecraftStateCache
    {
    public:
        struct Options
        {
            /// @brief Internal propagation step for cached segments. Smaller = more accurate, larger = faster.
            double lookup_dt_s{0.0};
            /// @brief Soft cap on the number of cached segments per spacecraft (0 = unlimited).
            std::size_t max_segments{4096};
        };

        SpacecraftStateCache(const std::vector<MassiveBody> &bodies, const EphemerisLike &eph, const ManeuverPlan &plan,
                             const double gravitational_constant, const double softening_length_m,
                             const DOPRI5Options &spacecraft_integrator, const double t0_s, const double t_end_s,
                             std::function<const Spacecraft *(SpacecraftId)> spacecraft_by_id, Options opt = {})
            : bodies_(bodies),
              eph_(eph),
              plan_(plan),
              gravitational_constant_(gravitational_constant),
              softening_length_m_(softening_length_m),
              spacecraft_integrator_(spacecraft_integrator),
              t0_s_(t0_s),
              t_end_s_(t_end_s),
              spacecraft_by_id_(std::move(spacecraft_by_id)),
              opt_(opt)
        {
        }

        /// @brief Time-aware lookup function suitable for passing into maneuver/propagation utilities.
        SpacecraftStateLookup lookup() const
        {
            return [this](const SpacecraftId id, const double t_s) { return this->state_at(id, t_s); };
        }

        /// @brief Get inertial spacecraft state at time `t_s` (clamped to cache interval).
        std::optional<State> state_at(const SpacecraftId id, const double t_s) const
        {
            if (id == kInvalidSpacecraftId)
            {
                return std::nullopt;
            }
            if (!std::isfinite(t_s))
            {
                return std::nullopt;
            }

            const double t_clamped = std::clamp(t_s, t0_s_, t_end_s_);

            Track *track = ensure_track_(id);
            if (track == nullptr)
            {
                return std::nullopt;
            }

            // If this track is already being extended (cycle), do not attempt to extend further.
            if (building_.contains(id))
            {
                const double t_safe = std::min(t_clamped, track->t_last_s);
                return state_from_track_(*track, t_safe);
            }

            if (t_clamped > track->t_last_s && track->t_last_s < t_end_s_)
            {
                extend_track_to_(*track, id, t_clamped);
            }

            if (track->segments.empty())
            {
                // Align with propagation semantics: apply any instantaneous impulses at the query time.
                building_.insert(id);
                const SpacecraftStateLookup sc_lookup = lookup();
                const State s = state_after_impulses_at_(id, t_clamped, track->sc_last.state, sc_lookup);
                building_.erase(id);
                return s;
            }

            return state_from_track_(*track, t_clamped);
        }

    private:
        struct Segment
        {
            double t0_s{0.0};
            double dt_s{0.0};
            State start{};
            State end{};
        };

        struct Track
        {
            Spacecraft sc_last{};
            double t_last_s{0.0};
            std::vector<Segment> segments{};
        };

        std::size_t auto_primary_at_(const double t_s, const Vec3 &pos_m) const
        {
            return auto_select_primary_index(
                    bodies_, pos_m, [&](std::size_t i) { return eph_.body_position_at(i, t_s); },
                    softening_length_m_);
        }

        std::size_t primary_for_impulse_(const ImpulseSegment &seg, const double t_s, const Vec3 &pos_m,
                                         const SpacecraftStateLookup &sc_lookup) const
        {
            return detail::resolve_primary_index(bodies_, seg, t_s, pos_m, sc_lookup,
                    [this](const double t, const Vec3 &p) { return auto_primary_at_(t, p); });
        }

        Vec3 impulse_dv_inertial_at_(const SpacecraftId id, const double t_s, const Vec3 &pos_m, const Vec3 &vel_mps,
                                     const SpacecraftStateLookup &sc_lookup) const
        {
            Vec3 dv_i{0.0, 0.0, 0.0};
            if (plan_.impulses.empty() || bodies_.empty())
            {
                return dv_i;
            }

            const double eps_t = 1e-9;
            for (const auto &imp: plan_.impulses)
            {
                if (!segment_applies_to_spacecraft(imp, id))
                {
                    continue;
                }
                if (!(std::abs(imp.t_s - t_s) <= eps_t))
                {
                    continue;
                }
                const std::size_t primary = primary_for_impulse_(imp, t_s, pos_m, sc_lookup);
                dv_i += rtn_vector_to_inertial(eph_, bodies_, primary, t_s, pos_m, vel_mps, imp.dv_rtn_mps,
                                               imp.rtn_frame, sc_lookup);
            }
            return dv_i;
        }

        State state_after_impulses_at_(const SpacecraftId id, const double t_s, const State &state_before,
                                       const SpacecraftStateLookup &sc_lookup) const
        {
            State out = state_before;
            const Vec3 dv_i = impulse_dv_inertial_at_(id, t_s, state_before.position_m, state_before.velocity_mps,
                                                      sc_lookup);
            if (glm::dot(dv_i, dv_i) > 0.0 && std::isfinite(dv_i.x) && std::isfinite(dv_i.y) && std::isfinite(dv_i.z))
            {
                out.velocity_mps += dv_i;
            }
            return out;
        }

        Track *ensure_track_(const SpacecraftId id) const
        {
            auto it = tracks_.find(id);
            if (it != tracks_.end())
            {
                return &it->second;
            }
            if (!spacecraft_by_id_)
            {
                return nullptr;
            }
            const Spacecraft *sc_ptr = spacecraft_by_id_(id);
            if (sc_ptr == nullptr)
            {
                return nullptr;
            }

            Track t;
            t.sc_last = *sc_ptr;
            t.sc_last.id = id;
            t.t_last_s = t0_s_;
            tracks_.insert({id, std::move(t)});
            return &tracks_.find(id)->second;
        }

        void extend_track_to_(Track &track, const SpacecraftId id, const double t_target_s) const
        {
            const double t_goal_s = std::min(t_target_s, t_end_s_);
            if (!(t_goal_s > track.t_last_s))
            {
                return;
            }

            building_.insert(id);
            const SpacecraftStateLookup sc_lookup = lookup();

            while (track.t_last_s < t_goal_s)
            {
                double remaining_s = t_goal_s - track.t_last_s;
                if (!(remaining_s > 0.0) || !std::isfinite(remaining_s))
                {
                    break;
                }

                const double burn_boundary_s = next_burn_boundary_after(plan_, id, track.t_last_s, t_goal_s);
                const double impulse_boundary_s = next_impulse_time_after(plan_, id, track.t_last_s, t_goal_s);
                const double boundary_s = std::min(t_goal_s, std::min(burn_boundary_s, impulse_boundary_s));
                const double remaining_to_boundary_s = boundary_s - track.t_last_s;
                if (!(remaining_to_boundary_s > 0.0) || !std::isfinite(remaining_to_boundary_s))
                {
                    break;
                }

                double dt_s = remaining_s;
                if (opt_.max_segments > 0 && track.segments.size() >= opt_.max_segments)
                {
                    // Coarsen: take as large a step as possible without crossing burn/impulse boundaries.
                    dt_s = std::min(dt_s, remaining_to_boundary_s);
                }
                else
                {
                    if (opt_.lookup_dt_s > 0.0 && std::isfinite(opt_.lookup_dt_s))
                    {
                        dt_s = std::min(dt_s, opt_.lookup_dt_s);
                    }
                    dt_s = std::min(dt_s, remaining_to_boundary_s);
                }

                if (!(dt_s > 0.0) || !std::isfinite(dt_s))
                {
                    break;
                }

                const Spacecraft sc0 = track.sc_last;
                const State start_state = state_after_impulses_at_(id, track.t_last_s, sc0.state, sc_lookup);
                const Spacecraft sc1 = detail::propagate_spacecraft_in_ephemeris(
                        sc0, bodies_, eph_, plan_, gravitational_constant_, softening_length_m_, spacecraft_integrator_,
                        track.t_last_s, dt_s, sc_lookup);

                track.segments.push_back(Segment{
                        .t0_s = track.t_last_s,
                        .dt_s = dt_s,
                        .start = start_state,
                        .end = sc1.state,
                });

                track.sc_last = sc1;
                track.t_last_s += dt_s;
            }

            building_.erase(id);
        }

        static State interpolate_state_(const Segment &seg, const double t_s)
        {
            State out{};
            if (!(seg.dt_s > 0.0) || !std::isfinite(seg.dt_s))
            {
                return seg.start;
            }
            const double tau = clamp01((t_s - seg.t0_s) / seg.dt_s);
            const double t_eval_s = seg.t0_s + tau * seg.dt_s;

            out.position_m = hermite_position(seg.start.position_m, seg.start.velocity_mps, seg.end.position_m,
                                              seg.end.velocity_mps, seg.dt_s, tau);
            out.velocity_mps = hermite_velocity_mps(seg.start.position_m, seg.start.velocity_mps, seg.end.position_m,
                                                    seg.end.velocity_mps, seg.dt_s, tau);

            out.spin.axis = normalized_or(seg.start.spin.axis, Vec3{0.0, 1.0, 0.0});
            out.spin.rate_rad_per_s = seg.start.spin.rate_rad_per_s;
            out.spin.angle_rad = seg.start.spin.angle_rad + seg.start.spin.rate_rad_per_s * (t_eval_s - seg.t0_s);
            return out;
        }

        static State state_from_track_(const Track &track, const double t_s)
        {
            if (track.segments.empty())
            {
                // Unextended track: sc_last is the initial state at t0_s_.
                return track.sc_last.state;
            }

            if (!(t_s > track.t_last_s))
            {
                if (!(t_s > track.segments.front().t0_s))
                {
                    return track.segments.front().start;
                }

                const auto it = std::upper_bound(
                        track.segments.begin(),
                        track.segments.end(),
                        t_s,
                        [](const double t, const Segment &seg) { return t < seg.t0_s; });

                if (it == track.segments.begin())
                {
                    return interpolate_state_(track.segments.front(), t_s);
                }
                if (it == track.segments.end())
                {
                    return interpolate_state_(track.segments.back(), t_s);
                }
                return interpolate_state_(*std::prev(it), t_s);
            }

            // If requested beyond built range, return the last known state.
            return track.sc_last.state;
        }

        const std::vector<MassiveBody> &bodies_;
        const EphemerisLike &eph_;
        const ManeuverPlan &plan_;
        double gravitational_constant_{kGravitationalConstant_SI};
        double softening_length_m_{0.0};
        DOPRI5Options spacecraft_integrator_{};
        double t0_s_{0.0};
        double t_end_s_{0.0};
        std::function<const Spacecraft *(SpacecraftId)> spacecraft_by_id_{};
        Options opt_{};

        mutable std::unordered_map<SpacecraftId, Track> tracks_{};
        mutable std::unordered_set<SpacecraftId> building_{};
    };

} // namespace orbitsim
