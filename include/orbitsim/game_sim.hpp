#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/detail/spacecraft_propagation.hpp"
#include "orbitsim/events.hpp"
#include "orbitsim/integrators.hpp"
#include "orbitsim/maneuvers.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace orbitsim
{

    class GameSimulation
    {
    public:
        struct Config
        {
            double gravitational_constant{orbitsim::kGravitationalConstant_SI};
            double softening_length_m{0.0};
            DOPRI5Options spacecraft_integrator{};
            EventOptions events{};
            bool enable_events{true};
        };

        GameSimulation() = default;
        explicit GameSimulation(Config cfg) : cfg_(std::move(cfg)) {}

        double time_s() const { return time_s_; }
        const Config &config() const { return cfg_; }

        std::vector<MassiveBody> &massive_bodies() { return massive_; }
        const std::vector<MassiveBody> &massive_bodies() const { return massive_; }

        std::vector<Spacecraft> &spacecraft() { return spacecraft_; }
        const std::vector<Spacecraft> &spacecraft() const { return spacecraft_; }

        ManeuverPlan &maneuver_plan() { return plan_; }
        const ManeuverPlan &maneuver_plan() const { return plan_; }

        std::size_t select_primary_by_max_accel(const std::size_t spacecraft_index) const
        {
            if (spacecraft_index >= spacecraft_.size() || massive_.empty())
            {
                return 0;
            }
            const Vec3 p = spacecraft_[spacecraft_index].state.position_m;
            const double eps2 = cfg_.softening_length_m * cfg_.softening_length_m;

            std::size_t best = 0;
            double best_a = -1.0;
            for (std::size_t i = 0; i < massive_.size(); ++i)
            {
                const Vec3 dr = massive_[i].state.position_m - p;
                const double r2 = glm::dot(dr, dr) + eps2;
                if (!(r2 > 0.0) || !std::isfinite(r2))
                {
                    continue;
                }
                const double amag = (cfg_.gravitational_constant * massive_[i].mass_kg) / r2;
                if (amag > best_a)
                {
                    best_a = amag;
                    best = i;
                }
            }
            return best;
        }

        void step(const double dt_s)
        {
            if (!(dt_s != 0.0) || !std::isfinite(dt_s))
            {
                return;
            }
            if (!cfg_.enable_events || cfg_.events.max_event_splits_per_step <= 0)
            {
                do_step_no_events_(dt_s);
                return;
            }
            if (dt_s < 0.0)
            {
                // Backwards integration with event splitting is not supported.
                do_step_no_events_(dt_s);
                return;
            }

            sort_segments_by_start(plan_);

            double remaining = dt_s;
            int splits = 0;
            while (remaining > 0.0 && splits++ < std::max(1, cfg_.events.max_event_splits_per_step))
            {
                // Preview a full remaining step to build a single ephemeris segment for event search.
                std::vector<MassiveBody> massive_preview = massive_;
                std::vector<Spacecraft> spacecraft_preview = spacecraft_;
                double t_preview = time_s_;

                CelestialEphemerisSegment eph_preview{};
                preview_step_no_events_(massive_preview, spacecraft_preview, t_preview, remaining, &eph_preview);

                std::optional<Event> best;
                for (std::size_t sc_index = 0; sc_index < spacecraft_.size(); ++sc_index)
                {
                    const auto &sc = spacecraft_[sc_index];
                    auto propagate_sc = [&](const Spacecraft &sc_start, const double t0_s,
                                            const double dt_sc_s) -> Spacecraft {
                        return propagate_spacecraft_(sc_start, eph_preview, t0_s, dt_sc_s, sc_index);
                    };
                    const std::optional<Event> e = find_earliest_event_in_interval(
                            massive_, eph_preview, sc, time_s_, remaining, plan_, cfg_.events, propagate_sc, sc_index);
                    if (e.has_value() && (!best.has_value() || e->t_event_s < best->t_event_s))
                    {
                        best = e;
                    }
                }

                if (!best.has_value())
                {
                    // No events found: apply preview result for massive bodies to avoid recomputation,
                    // but still need to propagate spacecraft since preview_step_no_events_ skips them.
                    massive_ = std::move(massive_preview);

                    for (std::size_t sc_index = 0; sc_index < spacecraft_.size(); ++sc_index)
                    {
                        spacecraft_[sc_index] = propagate_spacecraft_(spacecraft_[sc_index], eph_preview, time_s_, remaining, sc_index);
                    }

                    time_s_ = t_preview;
                    return;
                }

                double dt_event = best->t_event_s - time_s_;
                const double min_step = std::max(0.0, cfg_.events.time_tol_s);
                if (!(dt_event > min_step) || !std::isfinite(dt_event))
                {
                    dt_event = min_step;
                }
                if (dt_event > remaining)
                {
                    dt_event = remaining;
                }

                do_step_no_events_(dt_event);
                remaining -= dt_event;
            }

            if (remaining > 0.0)
            {
                do_step_no_events_(remaining);
            }
        }

	    private:
	        inline CelestialEphemerisSegment make_segment_(const std::vector<State> &start, const std::vector<State> &end,
	                                                       const double t0_s, const double dt_s) const
	        {
            CelestialEphemerisSegment eph;
            eph.t0_s = t0_s;
            eph.dt_s = dt_s;
            eph.start = start;
            eph.end = end;
            return eph;
	        }

	        inline Spacecraft propagate_spacecraft_(const Spacecraft &sc0, const CelestialEphemerisSegment &eph,
	                                                const double t0_s, const double dt_s,
	                                                const std::size_t spacecraft_index) const
	        {
	            return detail::propagate_spacecraft_in_ephemeris(
	                    sc0, massive_, eph, plan_, cfg_.gravitational_constant, cfg_.softening_length_m,
	                    cfg_.spacecraft_integrator, t0_s, dt_s, spacecraft_index);
	        }

	        inline void preview_step_no_events_(std::vector<MassiveBody> &massive, std::vector<Spacecraft> &spacecraft,
	                                            double &t_s, const double dt_s, CelestialEphemerisSegment *out_eph) const
	        {
            if (out_eph != nullptr)
            {
                *out_eph = {};
            }
            if (!(dt_s != 0.0) || !std::isfinite(dt_s))
            {
                return;
            }

	            std::vector<State> start_states;
	            detail::snapshot_states(massive, &start_states);

            symplectic4_step(massive, dt_s, cfg_.gravitational_constant, cfg_.softening_length_m);

	            std::vector<State> end_states;
	            detail::snapshot_states(massive, &end_states);

            CelestialEphemerisSegment eph = make_segment_(start_states, end_states, t_s, dt_s);
            if (out_eph != nullptr)
            {
                *out_eph = eph;
            }

            // Note: spacecraft preview ignores the maneuver plan mass/prop changes in the output container; event
            // search re-propagates from the true start state via propagate_spacecraft_ anyway.
            (void) spacecraft;
            t_s += dt_s;
        }

	        inline void do_step_no_events_(const double dt_s)
	        {
            if (!(dt_s != 0.0) || !std::isfinite(dt_s))
            {
                return;
            }

	            std::vector<State> start_states;
	            detail::snapshot_states(massive_, &start_states);

            symplectic4_step(massive_, dt_s, cfg_.gravitational_constant, cfg_.softening_length_m);

	            std::vector<State> end_states;
	            detail::snapshot_states(massive_, &end_states);

            const CelestialEphemerisSegment eph = make_segment_(start_states, end_states, time_s_, dt_s);

            for (std::size_t sc_index = 0; sc_index < spacecraft_.size(); ++sc_index)
            {
                spacecraft_[sc_index] = propagate_spacecraft_(spacecraft_[sc_index], eph, time_s_, dt_s, sc_index);
            }

            time_s_ += dt_s;
        }

        Config cfg_{};
        double time_s_{0.0};
        std::vector<MassiveBody> massive_{};
        std::vector<Spacecraft> spacecraft_{};
        ManeuverPlan plan_{};
    };

} // namespace orbitsim
