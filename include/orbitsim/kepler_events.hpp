#pragma once

#include "orbitsim/kepler_trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <vector>

namespace orbitsim
{
    enum class KeplerApsisKind : uint8_t
    {
        Periapsis,
        Apoapsis,
    };

    struct KeplerApsisEvent
    {
        KeplerApsisKind kind{KeplerApsisKind::Periapsis};
        double t_s{0.0};
        State state_relative{};
        double radius_m{0.0};
    };

    struct KeplerApsisEventOptions
    {
        double min_eccentricity{1.0e-6};
        double near_parabolic_epsilon{1.0e-8};
        double time_epsilon_s{1.0e-6};
        KeplerPropagationOptions propagation{};
    };

    enum class KeplerNodeKind : uint8_t
    {
        Ascending,
        Descending,
    };

    struct KeplerNodeEvent
    {
        KeplerNodeKind kind{KeplerNodeKind::Ascending};
        double t_s{0.0};
        State state_relative{};
    };

    struct KeplerNodeEventOptions
    {
        double min_inclination_rad{1.0e-6};
        double near_parabolic_epsilon{1.0e-8};
        double time_epsilon_s{1.0e-6};
        KeplerPropagationOptions propagation{};
    };

    enum class KeplerRadiusCrossingKind : uint8_t
    {
        Enter,
        Exit,
    };

    struct KeplerRadiusCrossingEvent
    {
        KeplerRadiusCrossingKind kind{KeplerRadiusCrossingKind::Enter};
        double t_s{0.0};
        State state_relative{};
        double radius_m{0.0};
    };

    struct KeplerRadiusCrossingOptions
    {
        double max_step_s{300.0};
        double refine_tolerance_s{0.25};
        double distance_tolerance_m{0.01};
        KeplerPropagationOptions propagation{};
    };

    struct KeplerClosestApproachEvent
    {
        double t_s{0.0};
        State subject_state_relative{};
        State subject_state_inertial{};
        State target_state_inertial{};
        BodyId target_body_id{kInvalidBodyId};
        double separation_m{0.0};
    };

    struct KeplerClosestApproachOptions
    {
        double max_step_s{300.0};
        double refine_tolerance_s{0.25};
        KeplerPropagationOptions propagation{};
    };

    namespace detail
    {
        inline bool time_in_kepler_arc_(const KeplerArc &arc,
                                        const double t_s,
                                        const double eps_s)
        {
            if (!std::isfinite(t_s) || !std::isfinite(eps_s))
            {
                return false;
            }
            if (arc.t1_s >= arc.t0_s)
            {
                return t_s + eps_s >= arc.t0_s && t_s <= arc.t1_s + eps_s;
            }
            return t_s - eps_s <= arc.t0_s && t_s >= arc.t1_s - eps_s;
        }

        inline bool push_apsis_event_(std::vector<KeplerApsisEvent> &out,
                                      const KeplerArc &arc,
                                      const KeplerApsisKind kind,
                                      const double t_s,
                                      const KeplerPropagationOptions &propagation)
        {
            const KeplerArcSample sample = sample_kepler_arc_state(arc, t_s, propagation);
            if (!sample.ok() || !kepler_finite3_(sample.state_relative.position_m))
            {
                return false;
            }

            const double radius_m = glm::length(sample.state_relative.position_m);
            if (!(radius_m > 0.0) || !std::isfinite(radius_m))
            {
                return false;
            }

            out.push_back(KeplerApsisEvent{
                    .kind = kind,
                    .t_s = sample.t_s,
                    .state_relative = sample.state_relative,
                    .radius_m = radius_m,
            });
            return true;
        }

        inline bool push_node_event_(std::vector<KeplerNodeEvent> &out,
                                     const KeplerArc &arc,
                                     const KeplerNodeKind kind,
                                     const double t_s,
                                     const KeplerPropagationOptions &propagation)
        {
            const KeplerArcSample sample = sample_kepler_arc_state(arc, t_s, propagation);
            if (!sample.ok() || !kepler_finite3_(sample.state_relative.position_m))
            {
                return false;
            }

            out.push_back(KeplerNodeEvent{
                    .kind = kind,
                    .t_s = sample.t_s,
                    .state_relative = sample.state_relative,
            });
            return true;
        }

        inline double radius_delta_m_(const KeplerArc &arc,
                                      const double t_s,
                                      const double radius_m,
                                      const KeplerPropagationOptions &propagation)
        {
            const KeplerArcSample sample = sample_kepler_arc_state(arc, t_s, propagation);
            if (!sample.ok() || !kepler_finite3_(sample.state_relative.position_m))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return glm::length(sample.state_relative.position_m) - radius_m;
        }

        inline bool radius_crossed_(const double a,
                                    const double b,
                                    const double tolerance_m)
        {
            if (!std::isfinite(a) || !std::isfinite(b))
            {
                return false;
            }
            const double eps = std::max(0.0, tolerance_m);
            return (a > eps && b <= eps) || (a < -eps && b >= -eps);
        }
    } // namespace detail

    inline std::vector<KeplerApsisEvent> find_kepler_apsis_events(
            const KeplerArc &arc,
            const KeplerApsisEventOptions &options = {})
    {
        std::vector<KeplerApsisEvent> out{};
        if (!kepler_arc_valid(arc))
        {
            return out;
        }

        const OrbitalElements elements =
                orbital_elements_from_relative_state(arc.mu_m3_s2,
                                                     arc.state0_relative.position_m,
                                                     arc.state0_relative.velocity_mps);
        const double e = elements.eccentricity;
        if (!(e > options.min_eccentricity) ||
            !std::isfinite(e) ||
            !std::isfinite(elements.semi_major_axis_m) ||
            !std::isfinite(elements.mean_anomaly_rad) ||
            std::abs(e - 1.0) <= options.near_parabolic_epsilon)
        {
            return out;
        }

        const double safe_eps_s = std::max(0.0, options.time_epsilon_s);
        if (e < 1.0 - options.near_parabolic_epsilon)
        {
            const OrbitScalars scalars = orbit_scalars_from_elements(arc.mu_m3_s2, elements);
            if (!(scalars.mean_motion_radps > 0.0) ||
                !std::isfinite(scalars.mean_motion_radps))
            {
                return out;
            }

            const auto event_time = [&](const double event_true_anomaly_rad) {
                const KeplerAnomalyResult anomaly =
                        kepler_anomalies_from_true_anomaly(e, event_true_anomaly_rad);
                if (!anomaly.valid || !std::isfinite(anomaly.mean_anomaly_rad))
                {
                    return std::numeric_limits<double>::quiet_NaN();
                }

                double dM = wrap_angle_0_2pi(anomaly.mean_anomaly_rad -
                                             elements.mean_anomaly_rad);
                if (std::abs(dM) <= 1.0e-12)
                {
                    dM = 0.0;
                }
                return arc.t0_s + dM / scalars.mean_motion_radps;
            };

            const double periapsis_t_s = event_time(0.0);
            if (detail::time_in_kepler_arc_(arc, periapsis_t_s, safe_eps_s))
            {
                detail::push_apsis_event_(out,
                                          arc,
                                          KeplerApsisKind::Periapsis,
                                          std::clamp(periapsis_t_s,
                                                     std::min(arc.t0_s, arc.t1_s),
                                                     std::max(arc.t0_s, arc.t1_s)),
                                          options.propagation);
            }

            const double apoapsis_t_s = event_time(std::acos(-1.0));
            if (detail::time_in_kepler_arc_(arc, apoapsis_t_s, safe_eps_s))
            {
                detail::push_apsis_event_(out,
                                          arc,
                                          KeplerApsisKind::Apoapsis,
                                          std::clamp(apoapsis_t_s,
                                                     std::min(arc.t0_s, arc.t1_s),
                                                     std::max(arc.t0_s, arc.t1_s)),
                                          options.propagation);
            }
        }
        else if (e > 1.0 + options.near_parabolic_epsilon)
        {
            const double abs_a = -elements.semi_major_axis_m;
            if (!(abs_a > 0.0) || !std::isfinite(abs_a))
            {
                return out;
            }

            const double mean_motion_radps = std::sqrt(arc.mu_m3_s2 /
                                                       (abs_a * abs_a * abs_a));
            if (!(mean_motion_radps > 0.0) || !std::isfinite(mean_motion_radps))
            {
                return out;
            }

            const double periapsis_t_s =
                    arc.t0_s - elements.mean_anomaly_rad / mean_motion_radps;
            if (detail::time_in_kepler_arc_(arc, periapsis_t_s, safe_eps_s))
            {
                detail::push_apsis_event_(out,
                                          arc,
                                          KeplerApsisKind::Periapsis,
                                          std::clamp(periapsis_t_s,
                                                     std::min(arc.t0_s, arc.t1_s),
                                                     std::max(arc.t0_s, arc.t1_s)),
                                          options.propagation);
            }
        }

        std::sort(out.begin(),
                  out.end(),
                  [](const KeplerApsisEvent &a, const KeplerApsisEvent &b) {
                      if (a.t_s == b.t_s)
                      {
                          return static_cast<uint8_t>(a.kind) <
                                 static_cast<uint8_t>(b.kind);
                      }
                      return a.t_s < b.t_s;
                  });
        return out;
    }

    inline std::vector<KeplerNodeEvent> find_kepler_node_events(
            const KeplerArc &arc,
            const Vec3 &reference_axis_unit_i,
            const KeplerNodeEventOptions &options = {})
    {
        std::vector<KeplerNodeEvent> out{};
        if (!kepler_arc_valid(arc))
        {
            return out;
        }

        const OrbitalElements elements =
                orbital_elements_from_relative_state_about_axis(arc.mu_m3_s2,
                                                                arc.state0_relative.position_m,
                                                                arc.state0_relative.velocity_mps,
                                                                reference_axis_unit_i);
        const double e = elements.eccentricity;
        if (!std::isfinite(e) ||
            !std::isfinite(elements.semi_major_axis_m) ||
            !std::isfinite(elements.mean_anomaly_rad) ||
            !std::isfinite(elements.inclination_rad) ||
            !std::isfinite(elements.arg_periapsis_rad) ||
            elements.inclination_rad <= options.min_inclination_rad ||
            std::abs(e - 1.0) <= options.near_parabolic_epsilon ||
            !(e < 1.0 - options.near_parabolic_epsilon))
        {
            return out;
        }

        const OrbitScalars scalars = orbit_scalars_from_elements(arc.mu_m3_s2, elements);
        if (!(scalars.mean_motion_radps > 0.0) ||
            !std::isfinite(scalars.mean_motion_radps))
        {
            return out;
        }

        const double safe_eps_s = std::max(0.0, options.time_epsilon_s);
        const auto event_time = [&](const double event_true_anomaly_rad) {
            const KeplerAnomalyResult anomaly =
                    kepler_anomalies_from_true_anomaly(e, event_true_anomaly_rad);
            if (!anomaly.valid || !std::isfinite(anomaly.mean_anomaly_rad))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }

            double dM = wrap_angle_0_2pi(anomaly.mean_anomaly_rad -
                                         elements.mean_anomaly_rad);
            if (std::abs(dM) <= 1.0e-12)
            {
                dM = 0.0;
            }
            return arc.t0_s + dM / scalars.mean_motion_radps;
        };

        const double pi = std::acos(-1.0);
        const double node_true_anomaly[2]{
                wrap_angle_0_2pi(-elements.arg_periapsis_rad),
                wrap_angle_0_2pi(pi - elements.arg_periapsis_rad),
        };
        const KeplerNodeKind node_kind[2]{
                KeplerNodeKind::Ascending,
                KeplerNodeKind::Descending,
        };

        for (int i = 0; i < 2; ++i)
        {
            const double t_s = event_time(node_true_anomaly[i]);
            if (!detail::time_in_kepler_arc_(arc, t_s, safe_eps_s))
            {
                continue;
            }

            detail::push_node_event_(out,
                                     arc,
                                     node_kind[i],
                                     std::clamp(t_s,
                                                std::min(arc.t0_s, arc.t1_s),
                                                std::max(arc.t0_s, arc.t1_s)),
                                     options.propagation);
        }

        std::sort(out.begin(),
                  out.end(),
                  [](const KeplerNodeEvent &a, const KeplerNodeEvent &b) {
                      if (a.t_s == b.t_s)
                      {
                          return static_cast<uint8_t>(a.kind) <
                                 static_cast<uint8_t>(b.kind);
                      }
                      return a.t_s < b.t_s;
                  });
        return out;
    }

    inline std::optional<KeplerRadiusCrossingEvent> find_next_kepler_radius_crossing(
            const KeplerArc &arc,
            const double radius_m,
            const double t_limit_s,
            const KeplerRadiusCrossingOptions &options = {})
    {
        if (!kepler_arc_valid(arc) ||
            !(radius_m > 0.0) ||
            !std::isfinite(radius_m) ||
            !std::isfinite(t_limit_s) ||
            t_limit_s <= arc.t0_s)
        {
            return std::nullopt;
        }

        const double t_end_s = std::min(t_limit_s, arc.t1_s);
        if (!std::isfinite(t_end_s) || t_end_s <= arc.t0_s)
        {
            return std::nullopt;
        }

        const double max_step_s = std::max(1.0e-6, options.max_step_s);
        const double refine_tolerance_s = std::max(0.0, options.refine_tolerance_s);
        const double distance_tolerance_m = std::max(0.0, options.distance_tolerance_m);

        double t0_s = arc.t0_s;
        double f0 = detail::radius_delta_m_(arc, t0_s, radius_m, options.propagation);
        if (!std::isfinite(f0))
        {
            return std::nullopt;
        }

        while (t0_s < t_end_s)
        {
            const double t1_s = std::min(t_end_s, t0_s + max_step_s);
            const double f1 = detail::radius_delta_m_(arc, t1_s, radius_m, options.propagation);
            if (!std::isfinite(f1))
            {
                return std::nullopt;
            }

            if (detail::radius_crossed_(f0, f1, distance_tolerance_m))
            {
                const KeplerRadiusCrossingKind kind =
                        f0 > f1 ? KeplerRadiusCrossingKind::Enter : KeplerRadiusCrossingKind::Exit;
                double a_s = t0_s;
                double b_s = t1_s;
                double fa = f0;
                for (int i = 0; i < 64; ++i)
                {
                    const double mid_s = 0.5 * (a_s + b_s);
                    const double fm = detail::radius_delta_m_(arc, mid_s, radius_m, options.propagation);
                    if (!std::isfinite(fm))
                    {
                        break;
                    }
                    if (std::abs(b_s - a_s) <= refine_tolerance_s ||
                        std::abs(fm) <= distance_tolerance_m)
                    {
                        a_s = mid_s;
                        b_s = mid_s;
                        break;
                    }
                    if ((fa > 0.0 && fm <= 0.0) || (fa < 0.0 && fm >= 0.0))
                    {
                        b_s = mid_s;
                    }
                    else
                    {
                        a_s = mid_s;
                        fa = fm;
                    }
                }

                const double event_t_s = 0.5 * (a_s + b_s);
                const KeplerArcSample sample =
                        sample_kepler_arc_state(arc, event_t_s, options.propagation);
                if (!sample.ok() || !detail::kepler_finite3_(sample.state_relative.position_m))
                {
                    return std::nullopt;
                }
                return KeplerRadiusCrossingEvent{
                        .kind = kind,
                        .t_s = sample.t_s,
                        .state_relative = sample.state_relative,
                        .radius_m = glm::length(sample.state_relative.position_m),
                };
            }

            t0_s = t1_s;
            f0 = f1;
        }

        return std::nullopt;
    }

    inline std::vector<KeplerRadiusCrossingEvent> find_kepler_radius_crossings(
            const KeplerArc &arc,
            const double radius_m,
            const double t_limit_s,
            const KeplerRadiusCrossingOptions &options = {})
    {
        std::vector<KeplerRadiusCrossingEvent> out{};
        if (!kepler_arc_valid(arc) ||
            !(radius_m > 0.0) ||
            !std::isfinite(radius_m) ||
            !std::isfinite(t_limit_s) ||
            t_limit_s <= arc.t0_s)
        {
            return out;
        }

        const double t_end_s = std::min(t_limit_s, arc.t1_s);
        if (!std::isfinite(t_end_s) || t_end_s <= arc.t0_s)
        {
            return out;
        }

        const double max_step_s = std::max(1.0e-6, options.max_step_s);
        const double refine_tolerance_s = std::max(0.0, options.refine_tolerance_s);
        const double distance_tolerance_m = std::max(0.0, options.distance_tolerance_m);
        const double nudge_s = std::max(1.0e-6, refine_tolerance_s);

        double t0_s = arc.t0_s;
        double f0 = detail::radius_delta_m_(arc, t0_s, radius_m, options.propagation);
        if (!std::isfinite(f0))
        {
            return out;
        }

        while (t0_s < t_end_s)
        {
            const double t1_s = std::min(t_end_s, t0_s + max_step_s);
            const double f1 = detail::radius_delta_m_(arc, t1_s, radius_m, options.propagation);
            if (!std::isfinite(f1))
            {
                break;
            }

            if (detail::radius_crossed_(f0, f1, distance_tolerance_m))
            {
                const KeplerRadiusCrossingKind kind =
                        f0 > f1 ? KeplerRadiusCrossingKind::Enter : KeplerRadiusCrossingKind::Exit;
                double a_s = t0_s;
                double b_s = t1_s;
                double fa = f0;
                for (int i = 0; i < 64; ++i)
                {
                    const double mid_s = 0.5 * (a_s + b_s);
                    const double fm = detail::radius_delta_m_(arc, mid_s, radius_m, options.propagation);
                    if (!std::isfinite(fm))
                    {
                        break;
                    }
                    if (std::abs(b_s - a_s) <= refine_tolerance_s ||
                        std::abs(fm) <= distance_tolerance_m)
                    {
                        a_s = mid_s;
                        b_s = mid_s;
                        break;
                    }
                    if ((fa > 0.0 && fm <= 0.0) || (fa < 0.0 && fm >= 0.0))
                    {
                        b_s = mid_s;
                    }
                    else
                    {
                        a_s = mid_s;
                        fa = fm;
                    }
                }

                const double event_t_s = 0.5 * (a_s + b_s);
                const KeplerArcSample sample =
                        sample_kepler_arc_state(arc, event_t_s, options.propagation);
                if (!sample.ok() || !detail::kepler_finite3_(sample.state_relative.position_m))
                {
                    break;
                }

                out.push_back(KeplerRadiusCrossingEvent{
                        .kind = kind,
                        .t_s = sample.t_s,
                        .state_relative = sample.state_relative,
                        .radius_m = glm::length(sample.state_relative.position_m),
                });

                t0_s = std::min(t_end_s, event_t_s + nudge_s);
                f0 = detail::radius_delta_m_(arc, t0_s, radius_m, options.propagation);
                if (!std::isfinite(f0))
                {
                    break;
                }
                continue;
            }

            t0_s = t1_s;
            f0 = f1;
        }

        return out;
    }

    inline std::optional<KeplerClosestApproachEvent> find_kepler_closest_approach_to_state_on_arc(
            const KeplerArc &arc,
            const std::function<bool(double, State &)> &primary_state_at,
            const std::function<bool(double, State &)> &target_state_at,
            const double t_limit_s,
            const KeplerClosestApproachOptions &options = {})
    {
        if (!kepler_arc_valid(arc) ||
            !primary_state_at ||
            !target_state_at ||
            !std::isfinite(t_limit_s) ||
            t_limit_s <= arc.t0_s)
        {
            return std::nullopt;
        }

        const double t_end_s = std::min(t_limit_s, arc.t1_s);
        if (!std::isfinite(t_end_s) || t_end_s <= arc.t0_s)
        {
            return std::nullopt;
        }

        auto sample_sep2 = [&](const double t_s,
                               State *subject_relative,
                               State *subject_inertial,
                               State *target_inertial) {
            const KeplerArcSample sample =
                    sample_kepler_arc_state(arc, t_s, options.propagation);
            if (!sample.ok() ||
                !detail::kepler_finite3_(sample.state_relative.position_m) ||
                !detail::kepler_finite3_(sample.state_relative.velocity_mps))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }

            State primary_state{};
            State target_state{};
            if (!primary_state_at(t_s, primary_state) ||
                !target_state_at(t_s, target_state) ||
                !detail::kepler_finite3_(primary_state.position_m) ||
                !detail::kepler_finite3_(primary_state.velocity_mps) ||
                !detail::kepler_finite3_(target_state.position_m))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }

            State inertial = sample.state_relative;
            inertial.position_m += primary_state.position_m;
            inertial.velocity_mps += primary_state.velocity_mps;
            if (!detail::kepler_finite3_(inertial.position_m) ||
                !detail::kepler_finite3_(inertial.velocity_mps))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }

            const Vec3 delta = inertial.position_m - target_state.position_m;
            const double sep2 = glm::dot(delta, delta);
            if (!(sep2 >= 0.0) || !std::isfinite(sep2))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }

            if (subject_relative)
            {
                *subject_relative = sample.state_relative;
            }
            if (subject_inertial)
            {
                *subject_inertial = inertial;
            }
            if (target_inertial)
            {
                *target_inertial = target_state;
            }
            return sep2;
        };

        const double max_step_s = std::max(1.0e-6, options.max_step_s);
        double best_t_s = arc.t0_s;
        double best_sep2 = std::numeric_limits<double>::infinity();
        double best_prev_t_s = arc.t0_s;
        double best_next_t_s = std::min(t_end_s, arc.t0_s + max_step_s);

        double prev_t_s = arc.t0_s;
        double prev_sep2 = sample_sep2(prev_t_s, nullptr, nullptr, nullptr);
        if (!std::isfinite(prev_sep2))
        {
            return std::nullopt;
        }
        best_sep2 = prev_sep2;

        while (prev_t_s < t_end_s)
        {
            const double t_s = std::min(t_end_s, prev_t_s + max_step_s);
            const double sep2 = sample_sep2(t_s, nullptr, nullptr, nullptr);
            if (!std::isfinite(sep2))
            {
                return std::nullopt;
            }
            if (sep2 < best_sep2)
            {
                best_t_s = t_s;
                best_sep2 = sep2;
                best_prev_t_s = prev_t_s;
                best_next_t_s = std::min(t_end_s, t_s + max_step_s);
            }
            prev_t_s = t_s;
            prev_sep2 = sep2;
            (void) prev_sep2;
        }

        double a_s = std::max(arc.t0_s, best_prev_t_s);
        double b_s = std::min(t_end_s, best_next_t_s);
        if (b_s < a_s)
        {
            std::swap(a_s, b_s);
        }

        const double tolerance_s = std::max(0.0, options.refine_tolerance_s);
        constexpr double kInvPhi = 0.6180339887498948482;
        constexpr double kInvPhi2 = 0.3819660112501051518;
        double c_s = a_s + kInvPhi2 * (b_s - a_s);
        double d_s = a_s + kInvPhi * (b_s - a_s);
        double fc = sample_sep2(c_s, nullptr, nullptr, nullptr);
        double fd = sample_sep2(d_s, nullptr, nullptr, nullptr);
        if (!std::isfinite(fc) || !std::isfinite(fd))
        {
            return std::nullopt;
        }

        for (int i = 0; i < 64 && std::abs(b_s - a_s) > tolerance_s; ++i)
        {
            if (fc < fd)
            {
                b_s = d_s;
                d_s = c_s;
                fd = fc;
                c_s = a_s + kInvPhi2 * (b_s - a_s);
                fc = sample_sep2(c_s, nullptr, nullptr, nullptr);
            }
            else
            {
                a_s = c_s;
                c_s = d_s;
                fc = fd;
                d_s = a_s + kInvPhi * (b_s - a_s);
                fd = sample_sep2(d_s, nullptr, nullptr, nullptr);
            }
            if (!std::isfinite(fc) || !std::isfinite(fd))
            {
                return std::nullopt;
            }
        }

        const double candidates[5]{
                arc.t0_s,
                t_end_s,
                best_t_s,
                c_s,
                d_s,
        };
        for (const double t_s : candidates)
        {
            if (t_s < arc.t0_s || t_s > t_end_s)
            {
                continue;
            }
            const double sep2 = sample_sep2(t_s, nullptr, nullptr, nullptr);
            if (std::isfinite(sep2) && sep2 < best_sep2)
            {
                best_sep2 = sep2;
                best_t_s = t_s;
            }
        }

        State subject_relative{};
        State subject_inertial{};
        State target_inertial{};
        const double sep2 = sample_sep2(best_t_s,
                                        &subject_relative,
                                        &subject_inertial,
                                        &target_inertial);
        if (!std::isfinite(sep2))
        {
            return std::nullopt;
        }

        return KeplerClosestApproachEvent{
                .t_s = best_t_s,
                .subject_state_relative = subject_relative,
                .subject_state_inertial = subject_inertial,
                .target_state_inertial = target_inertial,
                .separation_m = std::sqrt(std::max(0.0, sep2)),
        };
    }

    inline std::optional<KeplerClosestApproachEvent> find_kepler_closest_approach_on_arc(
            const KeplerArc &arc,
            const BodyId target_body_id,
            const std::function<bool(BodyId, double, State &)> &state_at,
            const double t_limit_s,
            const KeplerClosestApproachOptions &options = {})
    {
        if (target_body_id == kInvalidBodyId || !state_at)
        {
            return std::nullopt;
        }

        std::optional<KeplerClosestApproachEvent> event =
                find_kepler_closest_approach_to_state_on_arc(
                        arc,
                        [&](const double t_s, State &out_state) {
                            return state_at(arc.primary_body_id, t_s, out_state);
                        },
                        [&](const double t_s, State &out_state) {
                            return state_at(target_body_id, t_s, out_state);
                        },
                        t_limit_s,
                        options);
        if (event.has_value())
        {
            event->target_body_id = target_body_id;
        }
        return event;
    }
} // namespace orbitsim
