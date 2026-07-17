#pragma once

#include "orbitsim/kepler_trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

namespace orbitsim
{
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

    namespace detail
    {
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

        inline KeplerRadiusCrossingKind crossing_kind_(const State &state)
        {
            const double radial_speed_mps =
                    glm::dot(state.position_m, state.velocity_mps);
            return radial_speed_mps < 0.0 ? KeplerRadiusCrossingKind::Enter
                                          : KeplerRadiusCrossingKind::Exit;
        }

        inline bool push_crossing_(std::vector<KeplerRadiusCrossingEvent> &out,
                                                const KeplerArc &arc,
                                                const double t_s,
                                                const KeplerPropagationOptions &propagation)
        {
            const KeplerArcSample sample = sample_kepler_arc_state(arc, t_s, propagation);
            if (!sample.ok() ||
                !kepler_finite3_(sample.state_relative.position_m) ||
                !kepler_finite3_(sample.state_relative.velocity_mps))
            {
                return false;
            }

            const double radius_m = glm::length(sample.state_relative.position_m);
            if (!(radius_m > 0.0) || !std::isfinite(radius_m))
            {
                return false;
            }

            out.push_back(KeplerRadiusCrossingEvent{
                    .kind = crossing_kind_(sample.state_relative),
                    .t_s = sample.t_s,
                    .state_relative = sample.state_relative,
                    .radius_m = radius_m,
            });
            return true;
        }

        inline double signed_true_anomaly_(const double true_anomaly_rad)
        {
            const double two_pi = 2.0 * std::acos(-1.0);
            double nu = wrap_angle_0_2pi(true_anomaly_rad);
            if (nu > std::acos(-1.0))
            {
                nu -= two_pi;
            }
            return nu;
        }

        inline double barker_value_(const double true_anomaly_rad)
        {
            const double d = std::tan(0.5 * signed_true_anomaly_(true_anomaly_rad));
            return d + (d * d * d) / 3.0;
        }

        inline double refine_crossing_(const KeplerArc &arc,
                                                   const double radius_m,
                                                   const double t_guess_s,
                                                   const double t_min_s,
                                                   const double t_max_s,
                                                   const KeplerRadiusCrossingOptions &options)
        {
            if (!std::isfinite(t_guess_s) ||
                !std::isfinite(t_min_s) ||
                !std::isfinite(t_max_s) ||
                !(t_max_s >= t_min_s))
            {
                return t_guess_s;
            }

            double t_s = std::clamp(t_guess_s, t_min_s, t_max_s);
            const double distance_tolerance_m = std::max(0.0, options.distance_tolerance_m);
            const double time_tolerance_s = std::max(0.0, options.refine_tolerance_s);
            const double span_s = std::max(1.0e-9, t_max_s - t_min_s);
            for (int i = 0; i < 8; ++i)
            {
                const KeplerArcSample sample =
                        sample_kepler_arc_state(arc, t_s, options.propagation);
                if (!sample.ok() ||
                    !kepler_finite3_(sample.state_relative.position_m) ||
                    !kepler_finite3_(sample.state_relative.velocity_mps))
                {
                    break;
                }

                const double r_m = glm::length(sample.state_relative.position_m);
                if (!(r_m > 0.0) || !std::isfinite(r_m))
                {
                    break;
                }

                const double delta_m = r_m - radius_m;
                if (std::abs(delta_m) <= distance_tolerance_m)
                {
                    break;
                }

                const double radial_speed_mps =
                        glm::dot(sample.state_relative.position_m,
                                 sample.state_relative.velocity_mps) /
                        r_m;
                if (!(std::abs(radial_speed_mps) > 0.0) ||
                    !std::isfinite(radial_speed_mps))
                {
                    break;
                }

                const double step_s = std::clamp(delta_m / radial_speed_mps,
                                                 -0.25 * span_s,
                                                 0.25 * span_s);
                const double next_t_s = std::clamp(t_s - step_s, t_min_s, t_max_s);
                if (std::abs(next_t_s - t_s) <= time_tolerance_s)
                {
                    t_s = next_t_s;
                    break;
                }
                t_s = next_t_s;
            }
            return t_s;
        }

        struct AnalyticCrossings
        {
            // applied == true means the analytic solve was valid; empty events
            // then proves there is no crossing, so callers may skip the
            // numeric fallback scan.
            std::vector<KeplerRadiusCrossingEvent> events{};
            bool applied{false};
        };

        inline void dedupe_crossings_(std::vector<KeplerRadiusCrossingEvent> &events)
        {
            std::sort(events.begin(),
                      events.end(),
                      [](const KeplerRadiusCrossingEvent &a, const KeplerRadiusCrossingEvent &b) {
                          if (a.t_s == b.t_s)
                          {
                              return static_cast<uint8_t>(a.kind) <
                                     static_cast<uint8_t>(b.kind);
                          }
                          return a.t_s < b.t_s;
                      });
            events.erase(std::unique(events.begin(),
                                     events.end(),
                                     [](const KeplerRadiusCrossingEvent &a,
                                        const KeplerRadiusCrossingEvent &b) {
                                         return std::abs(a.t_s - b.t_s) <= 1.0e-6;
                                     }),
                         events.end());
        }

        inline AnalyticCrossings parabolic_crossings_(
                const KeplerArc &arc,
                const OrbitalElements &elements,
                const double radius_m,
                const double t_end_s,
                const KeplerRadiusCrossingOptions &options)
        {
            AnalyticCrossings out{};
            const Vec3 h_vec = glm::cross(arc.state0_relative.position_m,
                                          arc.state0_relative.velocity_mps);
            const double h2 = glm::dot(h_vec, h_vec);
            if (!(h2 > 0.0) ||
                !std::isfinite(h2) ||
                !(arc.mu_m3_s2 > 0.0) ||
                !std::isfinite(arc.mu_m3_s2))
            {
                return out;
            }

            const double p_m = h2 / arc.mu_m3_s2;
            const double q_m = 0.5 * p_m;
            if (!(q_m > 0.0) || !std::isfinite(q_m))
            {
                return out;
            }

            const double e = elements.eccentricity;
            if (!(e > 0.0) || !std::isfinite(e))
            {
                return out;
            }
            out.applied = true;
            if (radius_m < q_m)
            {
                return out;
            }

            const double cos_nu = (p_m / radius_m - 1.0) / e;
            constexpr double kCosEps = 1.0e-12;
            if (cos_nu < -1.0 - kCosEps || cos_nu > 1.0 + kCosEps)
            {
                return out;
            }

            const double clamped_cos_nu = std::clamp(cos_nu, -1.0, 1.0);
            const double nu0 = std::acos(clamped_cos_nu);
            const double true_anomalies[2]{-nu0, nu0};
            const double b0 = barker_value_(elements.true_anomaly_rad);
            const double time_scale_s = std::sqrt((2.0 * q_m * q_m * q_m) / arc.mu_m3_s2);
            if (!std::isfinite(b0) ||
                !(time_scale_s > 0.0) ||
                !std::isfinite(time_scale_s))
            {
                return out;
            }

            const double min_event_t_s =
                    arc.t0_s + std::max(1.0e-9, std::max(0.0, options.refine_tolerance_s) * 1.0e-6);
            for (const double true_anomaly_rad : true_anomalies)
            {
                const double b = barker_value_(true_anomaly_rad);
                if (!std::isfinite(b))
                {
                    continue;
                }

                double t_s = arc.t0_s + (b - b0) * time_scale_s;
                if (t_s <= min_event_t_s || t_s > t_end_s + 1.0e-9)
                {
                    continue;
                }

                t_s = refine_crossing_(arc,
                                                   radius_m,
                                                   t_s,
                                                   min_event_t_s,
                                                   t_end_s,
                                                   options);
                (void) push_crossing_(out.events, arc, t_s, options.propagation);
            }

            dedupe_crossings_(out.events);
            return out;
        }

        inline AnalyticCrossings analytic_crossings_(
                const KeplerArc &arc,
                const double radius_m,
                const double t_end_s,
                const KeplerRadiusCrossingOptions &options,
                const bool first_only = false)
        {
            AnalyticCrossings out{};
            const OrbitalElements elements =
                    orbital_elements_from_relative_state(arc.mu_m3_s2,
                                                         arc.state0_relative.position_m,
                                                         arc.state0_relative.velocity_mps);
            const double e = elements.eccentricity;
            const double a = elements.semi_major_axis_m;
            if (!(e > 1.0e-12) || !std::isfinite(e))
            {
                return out;
            }
            if (std::abs(e - 1.0) <= 1.0e-8)
            {
                return parabolic_crossings_(arc,
                                                        elements,
                                                        radius_m,
                                                        t_end_s,
                                                        options);
            }
            if (!std::isfinite(a) || !std::isfinite(elements.mean_anomaly_rad))
            {
                return out;
            }

            const double p_m = a * (1.0 - e * e);
            if (!(p_m > 0.0) || !std::isfinite(p_m))
            {
                return out;
            }

            const double abs_a = std::abs(a);
            const double mean_motion_radps =
                    std::sqrt(arc.mu_m3_s2 / (abs_a * abs_a * abs_a));
            if (!(mean_motion_radps > 0.0) || !std::isfinite(mean_motion_radps))
            {
                return out;
            }
            out.applied = true;

            const double cos_nu = (p_m / radius_m - 1.0) / e;
            constexpr double kCosEps = 1.0e-12;
            if (cos_nu < -1.0 - kCosEps || cos_nu > 1.0 + kCosEps)
            {
                return out;
            }

            const double clamped_cos_nu = std::clamp(cos_nu, -1.0, 1.0);
            const double nu0 = std::acos(clamped_cos_nu);
            const double two_pi = 2.0 * std::acos(-1.0);
            const double true_anomalies[2]{nu0, two_pi - nu0};
            const double min_event_t_s =
                    arc.t0_s + std::max(1.0e-9, std::max(0.0, options.refine_tolerance_s) * 1.0e-6);

            for (const double true_anomaly_rad : true_anomalies)
            {
                const KeplerAnomalyResult anomaly =
                        kepler_anomalies_from_true_anomaly(e, true_anomaly_rad);
                if (!anomaly.valid || !std::isfinite(anomaly.mean_anomaly_rad))
                {
                    // Anomaly solve failed: report not-applied so callers fall
                    // back to the numeric scan.
                    return {};
                }

                if (e < 1.0)
                {
                    double dM = wrap_angle_0_2pi(anomaly.mean_anomaly_rad -
                                                 elements.mean_anomaly_rad);
                    if (std::abs(dM) <= 1.0e-12)
                    {
                        dM = 0.0;
                    }
                    const double period_s = two_pi / mean_motion_radps;
                    if (!(period_s > 0.0) || !std::isfinite(period_s))
                    {
                        continue;
                    }

                    double t_s = arc.t0_s + dM / mean_motion_radps;
                    while (t_s <= min_event_t_s)
                    {
                        t_s += period_s;
                    }
                    for (; t_s <= t_end_s + 1.0e-9; t_s += period_s)
                    {
                        if (!push_crossing_(out.events, arc, t_s,
                                                         options.propagation) ||
                            first_only)
                        {
                            break;
                        }
                    }
                }
                else
                {
                    const double dM = anomaly.mean_anomaly_rad - elements.mean_anomaly_rad;
                    const double t_s = arc.t0_s + dM / mean_motion_radps;
                    if (t_s > min_event_t_s && t_s <= t_end_s + 1.0e-9)
                    {
                        (void) push_crossing_(out.events, arc, t_s,
                                                           options.propagation);
                    }
                }
            }

            dedupe_crossings_(out.events);
            return out;
        }
    } // namespace detail

    inline std::optional<KeplerRadiusCrossingEvent> next_kepler_crossing(
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

        const detail::AnalyticCrossings analytic =
                detail::analytic_crossings_(arc, radius_m, t_end_s, options, true);
        if (!analytic.events.empty())
        {
            return analytic.events.front();
        }
        if (analytic.applied)
        {
            return std::nullopt;
        }

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

    inline std::vector<KeplerRadiusCrossingEvent> find_kepler_crossings(
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

        detail::AnalyticCrossings analytic =
                detail::analytic_crossings_(arc, radius_m, t_end_s, options);
        if (!analytic.events.empty() || analytic.applied)
        {
            return std::move(analytic.events);
        }

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
} // namespace orbitsim
