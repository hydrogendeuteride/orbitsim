#pragma once

#include "orbitsim/kepler_trajectory.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <vector>

namespace orbitsim
{
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
        // Long searches may widen max_step_s enough to stay within max_steps;
        // Hermite extrema keep interior CA candidates visible in each interval.
        double max_step_s{300.0};
        double refine_tolerance_s{0.25};
        std::size_t max_steps{4096};
        KeplerPropagationOptions propagation{};
    };

    namespace detail
    {
        struct CaSample
        {
            double t_s{0.0};
            double sep2_m2{0.0};
            double range_rate_m2_s{0.0};
            State subject_relative{};
            State subject_inertial{};
            State target_inertial{};
        };

        inline std::optional<CaSample> ca_sample_(
                const KeplerArc &arc,
                const std::function<bool(double, State &)> &primary_state_at,
                const std::function<bool(double, State &)> &target_state_at,
                const double t_s,
                const KeplerPropagationOptions &propagation)
        {
            const KeplerArcSample sample = sample_kepler_arc_state(arc, t_s, propagation);
            if (!sample.ok() ||
                !kepler_finite3_(sample.state_relative.position_m) ||
                !kepler_finite3_(sample.state_relative.velocity_mps))
            {
                return std::nullopt;
            }

            State primary{};
            State target{};
            if (!primary_state_at(t_s, primary) ||
                !target_state_at(t_s, target) ||
                !kepler_finite3_(primary.position_m) ||
                !kepler_finite3_(primary.velocity_mps) ||
                !kepler_finite3_(target.position_m) ||
                !kepler_finite3_(target.velocity_mps))
            {
                return std::nullopt;
            }

            State inertial = sample.state_relative;
            inertial.position_m += primary.position_m;
            inertial.velocity_mps += primary.velocity_mps;
            if (!kepler_finite3_(inertial.position_m) ||
                !kepler_finite3_(inertial.velocity_mps))
            {
                return std::nullopt;
            }

            const Vec3 dr_m = inertial.position_m - target.position_m;
            const Vec3 dv_mps = inertial.velocity_mps - target.velocity_mps;
            const double sep2_m2 = glm::dot(dr_m, dr_m);
            const double range_rate_m2_s = glm::dot(dr_m, dv_mps);
            if (!(sep2_m2 >= 0.0) ||
                !std::isfinite(sep2_m2) ||
                !std::isfinite(range_rate_m2_s))
            {
                return std::nullopt;
            }

            return CaSample{
                    .t_s = t_s,
                    .sep2_m2 = sep2_m2,
                    .range_rate_m2_s = range_rate_m2_s,
                    .subject_relative = sample.state_relative,
                    .subject_inertial = inertial,
                    .target_inertial = target,
            };
        }

        struct CaExtrema
        {
            std::array<double, 2> u{};
            std::size_t count{0};
        };

        // Cubic Hermite interpolation of D^2 from D^2 and
        // d(D^2)/dt = 2 * dot(dr, dv) at both endpoints. Its derivative is
        // quadratic, so every predicted interior extremum falls out without
        // another propagated sample.
        inline CaExtrema ca_hermite_extrema_(const CaSample &low, const CaSample &high)
        {
            CaExtrema out{};
            const double h_s = high.t_s - low.t_s;
            if (!(h_s > 0.0) || !std::isfinite(h_s))
            {
                return out;
            }

            const double dy = high.sep2_m2 - low.sep2_m2;
            const double m0 = 2.0 * low.range_rate_m2_s;
            const double m1 = 2.0 * high.range_rate_m2_s;
            const double qa = 3.0 * (-2.0 * dy + h_s * (m0 + m1));
            const double qb = 2.0 * (3.0 * dy - h_s * (2.0 * m0 + m1));
            const double qc = h_s * m0;
            if (!std::isfinite(qa) || !std::isfinite(qb) || !std::isfinite(qc))
            {
                return out;
            }

            const double coef_scale = std::max({1.0, std::abs(qa), std::abs(qb), std::abs(qc)});
            const double coef_eps =
                    32.0 * std::numeric_limits<double>::epsilon() * coef_scale;
            auto keep = [&out](const double u) {
                constexpr double kEndpointEps = 1.0e-12;
                if (!std::isfinite(u) || u <= kEndpointEps || u >= 1.0 - kEndpointEps)
                {
                    return;
                }
                if (out.count > 0u && std::abs(out.u[out.count - 1u] - u) <= kEndpointEps)
                {
                    return;
                }
                out.u[out.count++] = u;
            };

            if (std::abs(qa) <= coef_eps)
            {
                if (std::abs(qb) > coef_eps)
                {
                    keep(-qc / qb);
                }
                return out;
            }

            const double qb2 = qb * qb;
            const double four_ac = 4.0 * qa * qc;
            double disc = qb2 - four_ac;
            const double disc_eps = 64.0 * std::numeric_limits<double>::epsilon() *
                                    std::max({1.0, std::abs(qb2), std::abs(four_ac)});
            if (!std::isfinite(disc) || disc < -disc_eps)
            {
                return out;
            }
            disc = std::max(0.0, disc);
            const double q = -0.5 * (qb + std::copysign(std::sqrt(disc), qb));
            const double u0 = q / qa;
            const double u1 = std::abs(q) > coef_eps ? qc / q : u0;
            keep(std::min(u0, u1));
            keep(std::max(u0, u1));
            return out;
        }

        // Bisect a range-rate sign change (approaching -> receding) down to
        // tolerance_s and return the endpoint closest to the true minimum.
        template <typename SampleAt>
        std::optional<CaSample> ca_refine_root_(const SampleAt &sample_at,
                                                CaSample low,
                                                CaSample high,
                                                const double tolerance_s)
        {
            if (!(low.range_rate_m2_s < 0.0) || !(high.range_rate_m2_s >= 0.0))
            {
                return std::nullopt;
            }

            for (int i = 0; i < 96 && high.t_s - low.t_s > tolerance_s; ++i)
            {
                const std::optional<CaSample> mid = sample_at(0.5 * (low.t_s + high.t_s));
                if (!mid.has_value())
                {
                    return std::nullopt;
                }
                (mid->range_rate_m2_s >= 0.0 ? high : low) = *mid;
            }

            return std::abs(low.range_rate_m2_s) < std::abs(high.range_rate_m2_s) ? low
                                                                                  : high;
        }
    } // namespace detail

    inline std::optional<KeplerClosestApproachEvent> find_kepler_approach(
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

        const auto sample_at = [&](const double t_s) {
            return detail::ca_sample_(arc, primary_state_at, target_state_at, t_s,
                                      options.propagation);
        };

        const double configured_step_s = std::max(1.0e-6, std::abs(options.max_step_s));
        const std::size_t max_steps = std::max<std::size_t>(1u, options.max_steps);
        const double budget_step_s =
                (t_end_s - arc.t0_s) / static_cast<double>(max_steps);
        const double step_s = std::max(configured_step_s, budget_step_s);
        const double tolerance_s = std::max(1.0e-9, std::abs(options.refine_tolerance_s));

        std::optional<detail::CaSample> best{};
        const auto inspect_bracket = [&](const detail::CaSample &low,
                                         const detail::CaSample &high) -> bool {
            if (!(low.range_rate_m2_s < 0.0) || !(high.range_rate_m2_s >= 0.0))
            {
                return true;
            }
            const std::optional<detail::CaSample> root =
                    detail::ca_refine_root_(sample_at, low, high, tolerance_s);
            if (!root.has_value())
            {
                return false;
            }
            if (!best.has_value() || root->sep2_m2 < best->sep2_m2)
            {
                best = root;
            }
            return true;
        };

        std::optional<detail::CaSample> previous = sample_at(arc.t0_s);
        if (!previous.has_value())
        {
            return std::nullopt;
        }

        while (previous->t_s < t_end_s)
        {
            const double t_s = std::min(t_end_s, previous->t_s + step_s);
            if (!(t_s > previous->t_s))
            {
                break;
            }

            const std::optional<detail::CaSample> next = sample_at(t_s);
            if (!next.has_value())
            {
                return std::nullopt;
            }

            const detail::CaExtrema extrema = detail::ca_hermite_extrema_(*previous, *next);
            if (extrema.count == 0u)
            {
                if (!inspect_bracket(*previous, *next))
                {
                    return std::nullopt;
                }
                previous = next;
                continue;
            }

            // Sample each predicted extremum plus the midpoints between them
            // so every range-rate sign change lands in some bracket.
            std::vector<double> probe_us{};
            probe_us.reserve(2u * extrema.count + 1u);
            double low_u = 0.0;
            for (std::size_t i = 0u; i < extrema.count; ++i)
            {
                probe_us.push_back(0.5 * (low_u + extrema.u[i]));
                probe_us.push_back(extrema.u[i]);
                low_u = extrema.u[i];
            }
            probe_us.push_back(0.5 * (low_u + 1.0));

            std::vector<detail::CaSample> samples{};
            samples.reserve(probe_us.size() + 2u);
            samples.push_back(*previous);
            for (const double u : probe_us)
            {
                const std::optional<detail::CaSample> probe =
                        sample_at(previous->t_s + (t_s - previous->t_s) * u);
                if (!probe.has_value())
                {
                    return std::nullopt;
                }
                samples.push_back(*probe);
            }
            samples.push_back(*next);

            for (std::size_t i = 1u; i < samples.size(); ++i)
            {
                if (!inspect_bracket(samples[i - 1u], samples[i]))
                {
                    return std::nullopt;
                }
            }
            previous = next;
        }

        if (!best.has_value())
        {
            return std::nullopt;
        }

        return KeplerClosestApproachEvent{
                .t_s = best->t_s,
                .subject_state_relative = best->subject_relative,
                .subject_state_inertial = best->subject_inertial,
                .target_state_inertial = best->target_inertial,
                .separation_m = std::sqrt(std::max(0.0, best->sep2_m2)),
        };
    }

    inline std::optional<KeplerClosestApproachEvent> find_kepler_approach(
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
                find_kepler_approach(
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
