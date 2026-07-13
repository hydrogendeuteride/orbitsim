#pragma once

#include "orbitsim/kepler.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace orbitsim
{
    struct KeplerArc
    {
        double mu_m3_s2{0.0};
        BodyId primary_body_id{kInvalidBodyId};
        double t0_s{0.0};
        double t1_s{0.0};
        State state0_relative{};
    };

    struct KeplerArcSample
    {
        double t_s{0.0};
        State state_relative{};
        KeplerDiagnostics diagnostics{};

        [[nodiscard]] bool ok() const { return diagnostics.status == KeplerStatus::Ok; }
    };

    struct KeplerTrajectoryOptions
    {
        double sample_dt_s{60.0};
        std::size_t max_samples{4096};
        bool include_start{true};
        bool include_end{true};
        KeplerPropagationOptions propagation{};
    };

    struct KeplerTrajectoryDiagnostics
    {
        std::size_t requested_samples{0};
        std::size_t accepted_samples{0};
        std::size_t failed_sample_index{0};
        KeplerStatus first_failure{KeplerStatus::Ok};
    };

    inline bool kepler_arc_valid(const KeplerArc &arc)
    {
        return arc.mu_m3_s2 > 0.0 && std::isfinite(arc.mu_m3_s2) && std::isfinite(arc.t0_s) &&
               std::isfinite(arc.t1_s) && detail::kepler_finite3_(arc.state0_relative.position_m) &&
               detail::kepler_finite3_(arc.state0_relative.velocity_mps);
    }

    namespace detail
    {
        inline bool solve_hyperbolic_anomaly_(const double eccentricity,
                                              const double mean_anomaly_rad,
                                              double &out_hyperbolic_anomaly_rad)
        {
            if (!(eccentricity > 1.0) || !std::isfinite(eccentricity) || !std::isfinite(mean_anomaly_rad))
            {
                return false;
            }

            double h = std::asinh(mean_anomaly_rad / eccentricity);
            if (!std::isfinite(h))
            {
                return false;
            }

            constexpr int kMaxIterations = 32;
            constexpr double kTolerance = 1.0e-12;
            for (int i = 0; i < kMaxIterations; ++i)
            {
                const double sinh_h = std::sinh(h);
                const double cosh_h = std::cosh(h);
                const double f = eccentricity * sinh_h - h - mean_anomaly_rad;
                const double df = eccentricity * cosh_h - 1.0;
                if (!std::isfinite(f) || !(std::abs(df) > 0.0) || !std::isfinite(df))
                {
                    return false;
                }

                const double step = f / df;
                h -= step;
                if (!std::isfinite(h))
                {
                    return false;
                }
                if (std::abs(step) <= kTolerance * std::max(1.0, std::abs(h)))
                {
                    out_hyperbolic_anomaly_rad = h;
                    return true;
                }
            }

            out_hyperbolic_anomaly_rad = h;
            return std::isfinite(h);
        }

        inline bool sample_hyperbolic_arc_state_from_elements_(const KeplerArc &arc,
                                                               const double t_s,
                                                               State &out_state)
        {
            if (!kepler_arc_valid(arc) || !std::isfinite(t_s))
            {
                return false;
            }

            const OrbitalElements elements =
                    orbital_elements_from_relative_state(arc.mu_m3_s2,
                                                         arc.state0_relative.position_m,
                                                         arc.state0_relative.velocity_mps);
            if (!(elements.eccentricity > 1.0) ||
                !(elements.semi_major_axis_m < 0.0) ||
                !std::isfinite(elements.semi_major_axis_m) ||
                !std::isfinite(elements.mean_anomaly_rad))
            {
                return false;
            }

            const double abs_a = -elements.semi_major_axis_m;
            const double mean_motion_radps = std::sqrt(arc.mu_m3_s2 / (abs_a * abs_a * abs_a));
            if (!(mean_motion_radps > 0.0) || !std::isfinite(mean_motion_radps))
            {
                return false;
            }

            const double mean_anomaly_rad =
                    elements.mean_anomaly_rad + mean_motion_radps * (t_s - arc.t0_s);
            double hyperbolic_anomaly_rad = 0.0;
            if (!solve_hyperbolic_anomaly_(elements.eccentricity,
                                           mean_anomaly_rad,
                                           hyperbolic_anomaly_rad))
            {
                return false;
            }

            const double tanh_half_h = std::tanh(0.5 * hyperbolic_anomaly_rad);
            const double true_anomaly_rad =
                    2.0 * std::atan(std::sqrt((elements.eccentricity + 1.0) /
                                              (elements.eccentricity - 1.0)) *
                                     tanh_half_h);
            if (!std::isfinite(true_anomaly_rad))
            {
                return false;
            }

            OrbitalElements propagated = elements;
            propagated.true_anomaly_rad = true_anomaly_rad;
            out_state = relative_state_from_orbital_elements(arc.mu_m3_s2, propagated);
            out_state.spin = arc.state0_relative.spin;
            return kepler_finite3_(out_state.position_m) &&
                   kepler_finite3_(out_state.velocity_mps);
        }
    } // namespace detail

    inline KeplerArcSample sample_kepler_arc_state(const KeplerArc &arc, const double t_s,
                                                   const KeplerPropagationOptions &opt = {})
    {
        KeplerArcSample out;
        out.t_s = t_s;
        if (!kepler_arc_valid(arc) || !std::isfinite(t_s))
        {
            out.state_relative = arc.state0_relative;
            out.diagnostics.status = KeplerStatus::InvalidInitialState;
            return out;
        }

        // Closed orbits are periodic; fold large offsets near t0 so the
        // universal-variable solve stays well-conditioned on long horizons.
        double dt_s = t_s - arc.t0_s;
        const double r0 = glm::length(arc.state0_relative.position_m);
        const double v0 = glm::length(arc.state0_relative.velocity_mps);
        if (r0 > 0.0 && std::isfinite(r0) && std::isfinite(v0))
        {
            const double alpha = 2.0 / r0 - (v0 * v0) / arc.mu_m3_s2;
            if (alpha > 0.0 && std::isfinite(alpha))
            {
                const double two_pi = 2.0 * std::acos(-1.0);
                const double period_s =
                        two_pi / (std::sqrt(arc.mu_m3_s2) * alpha * std::sqrt(alpha));
                if (period_s > 0.0 && std::isfinite(period_s))
                {
                    dt_s = std::remainder(dt_s, period_s);
                }
            }
        }

        const KeplerPropagationResult step =
                propagate_kepler_universal_safe(arc.mu_m3_s2, arc.state0_relative, dt_s, opt);
        out.state_relative = step.state;
        out.diagnostics = step.diagnostics;
        if (!out.ok())
        {
            State fallback_state{};
            if (detail::sample_hyperbolic_arc_state_from_elements_(arc, t_s, fallback_state))
            {
                out.state_relative = fallback_state;
                out.diagnostics.status = KeplerStatus::Ok;
                out.diagnostics.regime = KeplerOrbitRegime::Hyperbolic;
                out.diagnostics.used_fallback = true;
            }
        }
        return out;
    }

    inline std::vector<KeplerArcSample> build_kepler_arc_samples(const KeplerArc &arc,
                                                                 const KeplerTrajectoryOptions &opt = {},
                                                                 KeplerTrajectoryDiagnostics *diagnostics = nullptr)
    {
        if (diagnostics)
        {
            *diagnostics = {};
        }

        std::vector<KeplerArcSample> out;
        if (!kepler_arc_valid(arc) || !(opt.sample_dt_s > 0.0) || !std::isfinite(opt.sample_dt_s) ||
            opt.max_samples == 0u)
        {
            if (diagnostics)
            {
                diagnostics->first_failure = KeplerStatus::InvalidInitialState;
            }
            return out;
        }

        const double duration_s = arc.t1_s - arc.t0_s;
        const double abs_duration_s = std::abs(duration_s);
        const double dir = (duration_s >= 0.0) ? 1.0 : -1.0;

        auto push_sample = [&](const double t_s) {
            if (out.size() >= opt.max_samples)
            {
                return false;
            }

            KeplerArcSample sample = sample_kepler_arc_state(arc, t_s, opt.propagation);
            if (diagnostics)
            {
                ++diagnostics->requested_samples;
            }

            out.push_back(sample);
            if (sample.ok())
            {
                if (diagnostics)
                {
                    ++diagnostics->accepted_samples;
                }
                return true;
            }

            if (diagnostics && diagnostics->first_failure == KeplerStatus::Ok)
            {
                diagnostics->first_failure = sample.diagnostics.status;
                diagnostics->failed_sample_index = out.size() - 1u;
            }
            return false;
        };

        if (opt.include_start)
        {
            if (!push_sample(arc.t0_s))
            {
                return out;
            }
        }

        if (abs_duration_s == 0.0)
        {
            return out;
        }

        const std::size_t reserve_for_end = opt.include_end ? 1u : 0u;
        double offset_s = opt.sample_dt_s;
        while (offset_s < abs_duration_s && out.size() + reserve_for_end < opt.max_samples)
        {
            if (!push_sample(arc.t0_s + dir * offset_s))
            {
                return out;
            }
            offset_s += opt.sample_dt_s;
        }

        if (opt.include_end && out.size() < opt.max_samples)
        {
            const bool already_at_end = !out.empty() && std::abs(out.back().t_s - arc.t1_s) <= 1e-12;
            if (!already_at_end)
            {
                (void) push_sample(arc.t1_s);
            }
        }

        return out;
    }

    inline std::vector<Vec3> build_kepler_arc_polyline(const KeplerArc &arc,
                                                       const KeplerTrajectoryOptions &opt = {},
                                                       KeplerTrajectoryDiagnostics *diagnostics = nullptr)
    {
        const std::vector<KeplerArcSample> samples = build_kepler_arc_samples(arc, opt, diagnostics);
        std::vector<Vec3> out;
        out.reserve(samples.size());
        for (const KeplerArcSample &sample: samples)
        {
            if (sample.ok())
            {
                out.push_back(sample.state_relative.position_m);
            }
        }
        return out;
    }
} // namespace orbitsim
