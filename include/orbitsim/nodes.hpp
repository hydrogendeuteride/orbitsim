#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/events.hpp"
#include "orbitsim/types.hpp"

#include <cmath>
#include <optional>
#include <vector>

namespace orbitsim
{

    /// @brief Node crossing direction relative to a reference plane.
    enum class NodeCrossing
    {
        Ascending,
        Descending,
    };

    /**
     * @brief Event describing a spacecraft crossing a plane ("node" event).
     *
     * The plane is typically defined relative to a primary body and a plane normal in inertial coordinates.
     * If `target_spacecraft_id` is set, the plane may represent another spacecraft's orbital plane.
     */
    struct NodeEvent
    {
        double t_event_s{0.0};
        BodyId primary_body_id{kInvalidBodyId};
        SpacecraftId spacecraft_id{kInvalidSpacecraftId};
        SpacecraftId target_spacecraft_id{kInvalidSpacecraftId}; // For "target plane" nodes
        NodeCrossing crossing{NodeCrossing::Ascending};
    };

    namespace detail
    {
        template<class EphemerisLike, class Propagator>
        /**
         * @brief Find the earliest plane-node crossing in `[t0_s, t0_s + dt_s]` using bisection on signed distance.
         *
         * A crossing is detected when the signed distance to the plane changes sign (or either endpoint is within
         * `opt.dist_tol_m`). The event time is refined with bisection, where the spacecraft state is re-propagated
         * from `sc0` at each midpoint via `propagate_sc`.
         *
         * The returned crossing direction is determined from the relative velocity projected onto the plane normal
         * at the event time.
         *
         * @tparam EphemerisLike Type providing `body_position_at(index,t)` and `body_velocity_at(index,t)`.
         * @tparam Propagator Callable with signature `(Spacecraft sc0, double t0_s, double dt_s) -> Spacecraft`.
         *
         * @param bodies Massive bodies (used for id -> index mapping).
         * @param eph Ephemeris-like provider for primary body state.
         * @param sc0 Spacecraft state at `t0_s`.
         * @param t0_s Interval start time [s].
         * @param dt_s Interval length [s] (must be > 0).
         * @param primary_body_id Primary body id defining the relative position for plane distance.
         * @param plane_normal_unit_i Plane normal (expected unit length) expressed in inertial coordinates.
         * @param target_spacecraft_id Optional target spacecraft id (metadata for "target plane" nodes).
         * @param opt Event detection tolerances and max bisection iterations.
         * @param propagate_sc Propagator used to evaluate the spacecraft state at candidate times.
         * @return NodeEvent if a crossing is detected/refined; `std::nullopt` otherwise.
         */
        inline std::optional<NodeEvent> find_earliest_plane_node_in_interval_(
                const std::vector<MassiveBody> &bodies,
                const EphemerisLike &eph,
                const Spacecraft &sc0,
                const double t0_s,
                const double dt_s,
                const BodyId primary_body_id,
                const Vec3 &plane_normal_unit_i,
                const SpacecraftId target_spacecraft_id,
                const EventOptions &opt,
                Propagator propagate_sc,
                Spacecraft *out_sc1 = nullptr)
        {
            if (!(dt_s > 0.0) || !std::isfinite(dt_s) || !(opt.max_bisect_iters > 0))
            {
                return std::nullopt;
            }
            const double n2 = glm::dot(plane_normal_unit_i, plane_normal_unit_i);
            if (!(n2 > 0.0) || !std::isfinite(n2) || !finite3_(plane_normal_unit_i))
            {
                return std::nullopt;
            }

            const std::optional<std::size_t> primary_index_opt = body_index_for_id(bodies, primary_body_id);
            if (!primary_index_opt.has_value())
            {
                return std::nullopt;
            }
            const std::size_t primary_index = *primary_index_opt;

            const double t1_s = t0_s + dt_s;

            auto signed_dist_m = [&](const Spacecraft &sc, const double t_s) -> double {
                const Vec3 rp = eph.body_position_at(primary_index, t_s);
                const Vec3 r_rel = sc.state.position_m - rp;
                return glm::dot(plane_normal_unit_i, r_rel);
            };

            const Spacecraft sc1 = propagate_sc(sc0, t0_s, dt_s);
            if (out_sc1 != nullptr)
            {
                *out_sc1 = sc1;
            }
            const double f0 = signed_dist_m(sc0, t0_s);
            const double f1 = signed_dist_m(sc1, t1_s);
            if (!std::isfinite(f0) || !std::isfinite(f1))
            {
                return std::nullopt;
            }

            const double dist_tol_m = std::max(0.0, opt.dist_tol_m);
            const double time_tol_s = std::max(0.0, opt.time_tol_s);

            // Require a sign change (or a near-zero endpoint) to declare a crossing in this interval.
            const bool near0 = (std::abs(f0) <= dist_tol_m) || (std::abs(f1) <= dist_tol_m);
            const bool sign_change = (f0 <= 0.0 && f1 >= 0.0) || (f0 >= 0.0 && f1 <= 0.0);
            if (!(near0 || sign_change))
            {
                return std::nullopt;
            }

            auto bisection_time = [&](double fa0, double fb1) -> double {
                double a = t0_s;
                double b = t1_s;
                double fa = fa0;
                double fb = fb1;

                for (int it = 0; it < opt.max_bisect_iters; ++it)
                {
                    const double m = 0.5 * (a + b);
                    const Spacecraft scm = propagate_sc(sc0, t0_s, m - t0_s);
                    const double fm = signed_dist_m(scm, m);

                    if (!std::isfinite(fm))
                    {
                        break;
                    }
                    if (std::abs(b - a) <= time_tol_s || std::abs(fm) <= dist_tol_m)
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

            const double t_event = bisection_time(f0, f1);
            if (!std::isfinite(t_event))
            {
                return std::nullopt;
            }

            const Spacecraft sce = propagate_sc(sc0, t0_s, t_event - t0_s);
            const Vec3 vp = eph.body_velocity_at(primary_index, t_event);
            const Vec3 v_rel = sce.state.velocity_mps - vp;
            if (!finite3_(v_rel))
            {
                return std::nullopt;
            }

            NodeEvent out;
            out.t_event_s = t_event;
            out.primary_body_id = primary_body_id;
            out.spacecraft_id = sc0.id;
            out.target_spacecraft_id = target_spacecraft_id;
            out.crossing = (glm::dot(plane_normal_unit_i, v_rel) >= 0.0) ? NodeCrossing::Ascending
                                                                         : NodeCrossing::Descending;
            return out;
        }

    } // namespace detail

} // namespace orbitsim
