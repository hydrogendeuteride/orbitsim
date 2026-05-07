#pragma once

#include "orbitsim/coordinate_frames.hpp"
#include "orbitsim/kepler_trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace orbitsim
{
    struct KeplerImpulse
    {
        double t_s{0.0};
        Vec3 dv_rtn_mps{0.0, 0.0, 0.0};
    };

    struct KeplerManeuverDiagnostics
    {
        std::size_t candidate_impulses{0};
        std::size_t impulses_applied{0};
        std::size_t arcs_built{0};
        std::size_t failed_impulse_index{0};
        KeplerStatus first_failure{KeplerStatus::Ok};
    };

    inline Vec3 rtn_delta_v_to_inertial(const Vec3 &r_rel_m, const Vec3 &v_rel_mps, const Vec3 &dv_rtn_mps)
    {
        if (!detail::kepler_finite3_(r_rel_m) || !detail::kepler_finite3_(v_rel_mps) ||
            !detail::kepler_finite3_(dv_rtn_mps))
        {
            return Vec3{0.0, 0.0, 0.0};
        }

        const RtnFrame f = compute_rtn_frame(r_rel_m, v_rel_mps);
        return dv_rtn_mps.x * f.R + dv_rtn_mps.y * f.T + dv_rtn_mps.z * f.N;
    }

    inline State apply_impulse_rtn(const State &pre_impulse_relative, const Vec3 &dv_rtn_mps)
    {
        State out = pre_impulse_relative;
        out.velocity_mps +=
                rtn_delta_v_to_inertial(pre_impulse_relative.position_m, pre_impulse_relative.velocity_mps, dv_rtn_mps);
        return out;
    }

    namespace detail
    {
        inline bool kepler_times_equal_(const double a, const double b)
        {
            return std::abs(a - b) <= 1e-12;
        }

        inline bool kepler_time_in_arc_half_open_(const double t_s, const double t0_s, const double t1_s)
        {
            if (t1_s >= t0_s)
            {
                return t_s >= t0_s && t_s < t1_s;
            }
            return t_s <= t0_s && t_s > t1_s;
        }
    } // namespace detail

    inline std::vector<KeplerArc> build_maneuvered_kepler_arcs(const KeplerArc &base_arc,
                                                               std::span<const KeplerImpulse> impulses,
                                                               const KeplerPropagationOptions &propagation_options = {},
                                                               KeplerManeuverDiagnostics *diagnostics = nullptr)
    {
        if (diagnostics)
        {
            *diagnostics = {};
            diagnostics->candidate_impulses = impulses.size();
        }

        std::vector<KeplerArc> out;
        if (!kepler_arc_valid(base_arc))
        {
            if (diagnostics)
            {
                diagnostics->first_failure = KeplerStatus::InvalidInitialState;
            }
            return out;
        }

        std::vector<KeplerImpulse> active_impulses;
        active_impulses.reserve(impulses.size());
        for (const KeplerImpulse &impulse: impulses)
        {
            if (std::isfinite(impulse.t_s) && detail::kepler_finite3_(impulse.dv_rtn_mps) &&
                detail::kepler_time_in_arc_half_open_(impulse.t_s, base_arc.t0_s, base_arc.t1_s))
            {
                active_impulses.push_back(impulse);
            }
        }

        const bool forward = base_arc.t1_s >= base_arc.t0_s;
        std::stable_sort(active_impulses.begin(), active_impulses.end(), [forward](const KeplerImpulse &a,
                                                                                   const KeplerImpulse &b) {
            return forward ? (a.t_s < b.t_s) : (a.t_s > b.t_s);
        });

        if (active_impulses.empty())
        {
            out.push_back(base_arc);
            if (diagnostics)
            {
                diagnostics->arcs_built = out.size();
            }
            return out;
        }

        double cursor_t_s = base_arc.t0_s;
        State cursor_state = base_arc.state0_relative;
        for (std::size_t i = 0; i < active_impulses.size(); ++i)
        {
            const KeplerImpulse &impulse = active_impulses[i];
            if (!detail::kepler_times_equal_(impulse.t_s, cursor_t_s))
            {
                KeplerArc pre_arc = base_arc;
                pre_arc.t0_s = cursor_t_s;
                pre_arc.t1_s = impulse.t_s;
                pre_arc.state0_relative = cursor_state;
                out.push_back(pre_arc);

                const KeplerArcSample pre_impulse =
                        sample_kepler_arc_state(pre_arc, impulse.t_s, propagation_options);
                if (!pre_impulse.ok())
                {
                    if (diagnostics)
                    {
                        diagnostics->first_failure = pre_impulse.diagnostics.status;
                        diagnostics->failed_impulse_index = i;
                        diagnostics->arcs_built = out.size();
                    }
                    return out;
                }
                cursor_state = pre_impulse.state_relative;
                cursor_t_s = impulse.t_s;
            }

            cursor_state = apply_impulse_rtn(cursor_state, impulse.dv_rtn_mps);
            if (diagnostics)
            {
                ++diagnostics->impulses_applied;
            }
        }

        if (!detail::kepler_times_equal_(cursor_t_s, base_arc.t1_s))
        {
            KeplerArc final_arc = base_arc;
            final_arc.t0_s = cursor_t_s;
            final_arc.t1_s = base_arc.t1_s;
            final_arc.state0_relative = cursor_state;
            out.push_back(final_arc);
        }

        if (diagnostics)
        {
            diagnostics->arcs_built = out.size();
        }
        return out;
    }
} // namespace orbitsim
