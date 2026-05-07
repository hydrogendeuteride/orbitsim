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

        const KeplerPropagationResult step =
                propagate_kepler_universal_safe(arc.mu_m3_s2, arc.state0_relative, t_s - arc.t0_s, opt);
        out.state_relative = step.state;
        out.diagnostics = step.diagnostics;
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
