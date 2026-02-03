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

    /** @brief Types of trajectory events that can be detected. */
    enum class EventType
    {
        Impact, ///< Spacecraft enters body surface (radius + terrain)
        AtmosphereBoundary, ///< Spacecraft crosses atmosphere top boundary
        SoiBoundary, ///< Spacecraft crosses sphere of influence boundary
        Proximity, ///< Spacecraft-to-spacecraft distance threshold crossing
    };

    /** @brief Direction of boundary crossing. */
    enum class Crossing
    {
        Enter, ///< Crossing inward (distance decreasing through threshold)
        Exit, ///< Crossing outward (distance increasing through threshold)
    };

    /** @brief Describes a detected trajectory event. */
    struct Event
    {
        EventType type{EventType::Impact};
        BodyId body_id{kInvalidBodyId}; ///< Relevant body (invalid for Proximity)
        Crossing crossing{Crossing::Enter};
        double t_event_s{0.0}; ///< Time of event
        SpacecraftId spacecraft_id{kInvalidSpacecraftId}; ///< Spacecraft that triggered event
        SpacecraftId other_spacecraft_id{kInvalidSpacecraftId}; ///< Other spacecraft (Proximity only)
    };

    /** @brief Parameters for event detection algorithms. */
    struct EventOptions
    {
        double time_tol_s{1e-3}; ///< Bisection convergence tolerance in time
        double dist_tol_m{1e-2}; ///< Bisection convergence tolerance in distance
        int max_bisect_iters{64}; ///< Max iterations for root finding
        int max_event_splits_per_step{8}; ///< Max timestep subdivisions per step() call
    };

    namespace detail
    {
        template<class EvalF>
        inline double bisect_crossing_time_s(const double t0_s, const double t1_s, const double f0,
                                             const EventOptions &opt, EvalF eval_f)
        {
            double a = t0_s;
            double b = t1_s;
            double fa = f0;

            for (int it = 0; it < opt.max_bisect_iters; ++it)
            {
                const double m = 0.5 * (a + b);
                const double fm = eval_f(m);

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
                }
                else
                {
                    a = m;
                    fa = fm;
                }
            }
            return 0.5 * (a + b);
        }

    } // namespace detail

    /**
     * @brief Get the distance threshold for a given event type on a body.
     *
     * - Impact: radius + terrain_max_height
     * - AtmosphereBoundary: radius + atmosphere_top_height
     * - SoiBoundary: soi_radius
     *
     * @return true if threshold is valid for this body/event combination
     */
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

    /**
     * @brief Find the earliest boundary event for a spacecraft in a time interval.
     *
     * Checks all bodies for Impact, AtmosphereBoundary, and SoiBoundary crossings.
     * Uses bisection to refine the crossing time when a sign change in
     * (distance - threshold) is detected between interval endpoints.
     *
     * @tparam Propagator Callable: (sc0, t0_s, dt_s) -> Spacecraft
     * @param bodies All massive bodies to check against
     * @param eph Ephemeris segment for body positions
     * @param sc0 Spacecraft state at t0_s
     * @param t0_s Start time of interval
     * @param dt_s Interval duration (must be positive)
     * @param plan Maneuver plan (currently unused, reserved)
     * @param opt Event detection options
     * @param propagate_sc Propagator function to get spacecraft state at any time
     * @return Earliest event if any crossing detected, nullopt otherwise
     */
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
                const double t_event = detail::bisect_crossing_time_s(
                        t0_s, t1_s, f0, opt, [&](const double t_s) -> double {
                            const Spacecraft sc_at = propagate_sc(sc0, t0_s, t_s - t0_s);
                            return dist_to_body(body_index, t_s, sc_at.state.position_m) - threshold_m;
                        });
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

    /**
     * @brief Find proximity event between two spacecraft in a time interval.
     *
     * Detects when the distance between center and target spacecraft crosses
     * the given threshold. Uses bisection for precise crossing time.
     *
     * @tparam Propagator Callable: (sc0, t0_s, dt_s) -> Spacecraft
     * @param eph Ephemeris segment (for body positions if needed)
     * @param center0 Reference spacecraft at t0_s
     * @param target0 Target spacecraft at t0_s
     * @param t0_s Start time of interval
     * @param dt_s Interval duration (must be positive)
     * @param threshold_m Distance threshold for event
     * @param opt Event detection options
     * @param propagate_sc Propagator function
     * @return Proximity event if crossing detected, nullopt otherwise
     */
    template<class Propagator>
    inline std::optional<Event>
    find_earliest_proximity_event_in_interval(const CelestialEphemerisSegment &eph, const Spacecraft &center0,
                                              const Spacecraft &target0, const double t0_s, const double dt_s,
                                              const double threshold_m, const EventOptions &opt,
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
        const double t_event = detail::bisect_crossing_time_s(
                t0_s, t1_s, f0, opt, [&](const double t_s) -> double {
                    const Spacecraft cm = propagate_sc(center0, t0_s, t_s - t0_s);
                    const Spacecraft tm = propagate_sc(target0, t0_s, t_s - t0_s);
                    return dist_m(t_s, cm.state.position_m, tm.state.position_m) - threshold_m;
                });
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
