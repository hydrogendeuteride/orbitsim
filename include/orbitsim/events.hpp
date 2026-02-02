#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/maneuvers.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>

namespace orbitsim
{

    enum class EventType
    {
        Impact,
        AtmosphereBoundary,
        SoiBoundary,
        Proximity,
    };

    enum class Crossing
    {
        Enter,
        Exit,
    };

    struct Event
    {
        EventType type{EventType::Impact};
        BodyId body_id{kInvalidBodyId};
        Crossing crossing{Crossing::Enter};
        double t_event_s{0.0};
        SpacecraftId spacecraft_id{kInvalidSpacecraftId};
        SpacecraftId other_spacecraft_id{kInvalidSpacecraftId}; // Used for Proximity events.
    };

    struct EventOptions
    {
        double time_tol_s{1e-3};
        double dist_tol_m{1e-2};
        int max_bisect_iters{64};
        int max_event_splits_per_step{8};
    };

    inline bool boundary_threshold_m(const MassiveBody &body, const EventType type, double *out_threshold_m)
    {
        if (out_threshold_m == nullptr)
        {
            return false;
        }
        *out_threshold_m = 0.0;
        if (!(body.radius_m > 0.0) || !std::isfinite(body.radius_m))
        {
            return false;
        }

        if (type == EventType::Impact)
        {
            const double terrain = (std::isfinite(body.terrain_max_height_m) ? body.terrain_max_height_m : 0.0);
            *out_threshold_m = body.radius_m + std::max(0.0, terrain);
            return (*out_threshold_m > 0.0);
        }
        if (type == EventType::AtmosphereBoundary)
        {
            if (!(body.atmosphere_top_height_m > 0.0) || !std::isfinite(body.atmosphere_top_height_m))
            {
                return false;
            }
            *out_threshold_m = body.radius_m + body.atmosphere_top_height_m;
            return (*out_threshold_m > 0.0);
        }
        if (type == EventType::SoiBoundary)
        {
            if (!(body.soi_radius_m > 0.0) || !std::isfinite(body.soi_radius_m))
            {
                return false;
            }
            *out_threshold_m = body.soi_radius_m;
            return (*out_threshold_m > 0.0);
        }
        return false;
    }

    template<class Propagator>
    inline std::optional<Event>
    find_earliest_event_in_interval(const std::vector<MassiveBody> &bodies, const CelestialEphemerisSegment &eph,
                                    const Spacecraft &sc0, const double t0_s, const double dt_s,
                                    const ManeuverPlan &plan, const EventOptions &opt, Propagator propagate_sc)
    {
        (void) plan;
        if (!(dt_s > 0.0) || !std::isfinite(dt_s) || !(opt.max_bisect_iters > 0))
        {
            return std::nullopt;
        }

        const double t1_s = t0_s + dt_s;

        auto dist_to_body = [&](const std::size_t body_index, const double t_s, const Vec3 &sc_pos_m) -> double {
            const Vec3 bp = eph.body_position_at(body_index, t_s);
            return glm::length(sc_pos_m - bp);
        };

        auto bisection_time = [&](const std::size_t body_index, const EventType type, const double threshold_m,
                                  const double fa0, const double fb1) -> double {
            double a = t0_s;
            double b = t1_s;
            double fa = fa0;
            double fb = fb1;

            for (int it = 0; it < opt.max_bisect_iters; ++it)
            {
                const double m = 0.5 * (a + b);
                const Spacecraft scm = propagate_sc(sc0, t0_s, m - t0_s);
                const double fm = dist_to_body(body_index, m, scm.state.position_m) - threshold_m;

                if (!std::isfinite(fm))
                {
                    break;
                }
                if (std::abs(b - a) <= std::max(0.0, opt.time_tol_s) || std::abs(fm) <= std::max(0.0, opt.dist_tol_m))
                {
                    return m;
                }
                const bool left = (fa <= 0.0 && fm >= 0.0) || (fa >= 0.0 && fm <= 0.0);
                if (left)
                {
                    b = m;
                    fb = fm;
                }
                else
                {
                    a = m;
                    fa = fm;
                }
            }
            return 0.5 * (a + b);
        };

        const Spacecraft sc1 = propagate_sc(sc0, t0_s, dt_s);

        std::optional<Event> best;
        for (std::size_t body_index = 0; body_index < bodies.size(); ++body_index)
        {
            for (EventType type: {EventType::Impact, EventType::AtmosphereBoundary, EventType::SoiBoundary})
            {
                double threshold_m = 0.0;
                if (!boundary_threshold_m(bodies[body_index], type, &threshold_m))
                {
                    continue;
                }

                const double d0 = dist_to_body(body_index, t0_s, sc0.state.position_m);
                const double d1 = dist_to_body(body_index, t1_s, sc1.state.position_m);
                if (!std::isfinite(d0) || !std::isfinite(d1))
                {
                    continue;
                }
                const double f0 = d0 - threshold_m;
                const double f1 = d1 - threshold_m;

                if (!((f0 <= 0.0 && f1 >= 0.0) || (f0 >= 0.0 && f1 <= 0.0)))
                {
                    continue;
                }

                const Crossing crossing = (f0 > 0.0 && f1 <= 0.0) ? Crossing::Enter : Crossing::Exit;
                const double t_event = bisection_time(body_index, type, threshold_m, f0, f1);
                if (!std::isfinite(t_event))
                {
                    continue;
                }

                if (!best.has_value() || t_event < best->t_event_s)
                {
                    best = Event{
                            .type = type,
                            .body_id = bodies[body_index].id,
                            .crossing = crossing,
                            .t_event_s = t_event,
                            .spacecraft_id = sc0.id,
                    };
                }
            }
        }

        return best;
    }

    template<class Propagator>
    inline std::optional<Event> find_earliest_proximity_event_in_interval(
            const CelestialEphemerisSegment &eph,
            const Spacecraft &center0,
            const Spacecraft &target0,
            const double t0_s,
            const double dt_s,
            const double threshold_m,
            const EventOptions &opt,
            Propagator propagate_sc)
    {
        if (!(dt_s > 0.0) || !std::isfinite(dt_s) || !(threshold_m > 0.0) || !std::isfinite(threshold_m) ||
            !(opt.max_bisect_iters > 0))
        {
            return std::nullopt;
        }

        const double t1_s = t0_s + dt_s;

        auto dist_m = [&](const double t_s, const Vec3 &pos_center_m, const Vec3 &pos_target_m) -> double {
            (void) eph;
            (void) t_s;
            return glm::length(pos_target_m - pos_center_m);
        };

        auto bisection_time = [&](const double fa0, const double fb1) -> double {
            double a = t0_s;
            double b = t1_s;
            double fa = fa0;
            double fb = fb1;

            for (int it = 0; it < opt.max_bisect_iters; ++it)
            {
                const double m = 0.5 * (a + b);
                const Spacecraft cm = propagate_sc(center0, t0_s, m - t0_s);
                const Spacecraft tm = propagate_sc(target0, t0_s, m - t0_s);
                const double fm = dist_m(m, cm.state.position_m, tm.state.position_m) - threshold_m;

                if (!std::isfinite(fm))
                {
                    break;
                }
                if (std::abs(b - a) <= std::max(0.0, opt.time_tol_s) || std::abs(fm) <= std::max(0.0, opt.dist_tol_m))
                {
                    return m;
                }
                const bool left = (fa <= 0.0 && fm >= 0.0) || (fa >= 0.0 && fm <= 0.0);
                if (left)
                {
                    b = m;
                    fb = fm;
                }
                else
                {
                    a = m;
                    fa = fm;
                }
            }
            return 0.5 * (a + b);
        };

        const Spacecraft center1 = propagate_sc(center0, t0_s, dt_s);
        const Spacecraft target1 = propagate_sc(target0, t0_s, dt_s);

        const double d0 = dist_m(t0_s, center0.state.position_m, target0.state.position_m);
        const double d1 = dist_m(t1_s, center1.state.position_m, target1.state.position_m);
        if (!std::isfinite(d0) || !std::isfinite(d1))
        {
            return std::nullopt;
        }
        const double f0 = d0 - threshold_m;
        const double f1 = d1 - threshold_m;

        if (!((f0 <= 0.0 && f1 >= 0.0) || (f0 >= 0.0 && f1 <= 0.0)))
        {
            return std::nullopt;
        }

        const Crossing crossing = (f0 > 0.0 && f1 <= 0.0) ? Crossing::Enter : Crossing::Exit;
        const double t_event = bisection_time(f0, f1);
        if (!std::isfinite(t_event))
        {
            return std::nullopt;
        }

        return Event{
                .type = EventType::Proximity,
                .body_id = kInvalidBodyId,
                .crossing = crossing,
                .t_event_s = t_event,
                .spacecraft_id = target0.id,
                .other_spacecraft_id = center0.id,
        };
    }

} // namespace orbitsim
