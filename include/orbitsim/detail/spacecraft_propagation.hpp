#pragma once

#include "orbitsim/integrators.hpp"
#include "orbitsim/maneuvers.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

namespace orbitsim::detail
{

    /** @brief Extract State from each MassiveBody into a vector. */
    inline void snapshot_states(const std::vector<MassiveBody> &massive, std::vector<State> *out)
    {
        if (out == nullptr)
        {
            return;
        }
        out->clear();
        out->reserve(massive.size());
        for (const auto &b: massive)
        {
            out->push_back(b.state);
        }
    }

    /**
     * @brief Compute gravitational acceleration at a point from all massive bodies.
     *
     * Queries body positions from the ephemeris at time t_s, then sums
     * contributions using Newton's law with optional softening.
     *
     * @tparam EphemerisLike Type with body_position_at(index, t_s) method
     * @param eph Ephemeris providing body positions at arbitrary times
     * @param bodies Massive bodies (for mass values)
     * @param G Gravitational constant
     * @param softening_length_m Softening to avoid singularities
     * @param t_s Evaluation time
     * @param pos_m Position at which to compute acceleration
     * @return Total gravitational acceleration (m/s^2)
     */
    template<class EphemerisLike>
    inline Vec3 gravity_accel_mps2(const EphemerisLike &eph, const std::vector<MassiveBody> &bodies, const double G,
                                   const double softening_length_m, const double t_s, const Vec3 &pos_m)
    {
        Vec3 a{0.0, 0.0, 0.0};
        if (!(G > 0.0) || !std::isfinite(G))
        {
            return a;
        }

        const double eps2 = softening_length_m * softening_length_m;
        for (std::size_t i = 0; i < bodies.size(); ++i)
        {
            const Vec3 p = eph.body_position_at(i, t_s);
            const Vec3 dr = p - pos_m;
            const double r2 = glm::dot(dr, dr) + eps2;
            if (!(r2 > 0.0) || !std::isfinite(r2))
            {
                continue;
            }
            const double inv_r = 1.0 / std::sqrt(r2);
            const double inv_r3 = inv_r * inv_r * inv_r;
            a += (G * bodies[i].mass_kg) * inv_r3 * dr;
        }
        return a;
    }

    /**
     * @brief Resolve the primary body index for a maneuver segment (burn or impulse).
     *
     * Checks, in order: explicit primary_body_id, LVLH frame's primary_body_id,
     * LVLH target spacecraft's dominant body, then falls back to auto-selection.
     *
     * @tparam Segment BurnSegment or ImpulseSegment (must have primary_body_id, rtn_frame)
     * @tparam AutoPrimaryFn Callable (double t_s, Vec3 pos_m) -> std::size_t
     */
    template<class Segment, class AutoPrimaryFn>
    inline std::size_t resolve_primary_index(const std::vector<MassiveBody> &bodies, const Segment &seg,
                                             const double t_s, const Vec3 &pos_m,
                                             const SpacecraftStateLookup &sc_lookup, AutoPrimaryFn &&auto_primary_at)
    {
        const auto primary_opt = body_index_for_id(bodies, seg.primary_body_id);
        if (primary_opt)
        {
            return *primary_opt;
        }

        if (seg.rtn_frame.type == TrajectoryFrameType::LVLH)
        {
            const auto frame_primary_opt = body_index_for_id(bodies, seg.rtn_frame.primary_body_id);
            if (frame_primary_opt)
            {
                return *frame_primary_opt;
            }

            if (sc_lookup)
            {
                const std::optional<State> target_state = sc_lookup(seg.rtn_frame.target_spacecraft_id, t_s);
                if (target_state.has_value())
                {
                    return auto_primary_at(t_s, target_state->position_m);
                }
            }
        }

        return auto_primary_at(t_s, pos_m);
    }

    /**
     * @brief Propagate a spacecraft through a precomputed ephemeris segment.
     *
     * Integrates spacecraft motion under:
     * - Gravitational attraction from all massive bodies (via ephemeris lookup)
     * - Thrust from active burns in the maneuver plan (in RTN frame)
     * - Instantaneous impulse maneuvers at scheduled times
     *
     * Automatically subdivides at burn boundaries and handles propellant depletion.
     * Backward integration (dt_s < 0) is gravity-only (no maneuvers).
     *
     * @tparam EphemerisLike Type providing body_position_at(index, t_s)
     * @param sc0 Initial spacecraft state
     * @param bodies Massive bodies (for mass and ID lookup)
     * @param eph Ephemeris for body position interpolation
     * @param plan Scheduled burns and impulses
     * @param gravitational_constant G constant
     * @param softening_length_m Softening length
     * @param spacecraft_integrator DOPRI5 options
     * @param t0_s Start time
     * @param dt_s Integration interval (can be negative)
     * @param sc_lookup Optional callback to look up other spacecraft states for LVLH frame
     * @return Updated spacecraft with new state and reduced propellant
     */
    template<class EphemerisLike>
    inline Spacecraft propagate_spacecraft_in_ephemeris(
            const Spacecraft &sc0, const std::vector<MassiveBody> &bodies, const EphemerisLike &eph,
            const ManeuverPlan &plan, const double gravitational_constant, const double softening_length_m,
            const DOPRI5Options &spacecraft_integrator, const double t0_s, const double dt_s,
            const SpacecraftStateLookup &sc_lookup = nullptr)
    {
        Spacecraft out = sc0;
        if (!(dt_s != 0.0) || !std::isfinite(dt_s))
        {
            return out;
        }

        if (dt_s < 0.0)
        {
            // Backwards integration currently supports gravity-only (no maneuvers, no prop changes).
            const SpacecraftKinematics y0{out.state.position_m, out.state.velocity_mps};
            auto accel = [&](double t_eval_s, const Vec3 &pos_m, const Vec3 & /*vel_mps*/) -> Vec3 {
                return gravity_accel_mps2(eph, bodies, gravitational_constant, softening_length_m, t_eval_s, pos_m);
            };
            const SpacecraftKinematics y1 =
                    dopri5_integrate_interval(t0_s, y0, dt_s, accel, spacecraft_integrator);

            out.state.position_m = y1.position_m;
            out.state.velocity_mps = y1.velocity_mps;
            advance_spin(out.state.spin, dt_s);
            return out;
        }

        const double t_end_s = t0_s + dt_s;
        double t = t0_s;
        double remaining = dt_s;

        const SpacecraftId spacecraft_id = sc0.id;
        SpacecraftKinematics y{out.state.position_m, out.state.velocity_mps};

        auto auto_primary_at = [&](const double t_s, const Vec3 &pos_m) -> std::size_t {
            return auto_select_primary_index(
                    bodies, pos_m, [&](std::size_t i) { return eph.body_position_at(i, t_s); }, softening_length_m);
        };

        auto resolve_primary = [&](const auto &seg, const double t_s, const Vec3 &pos_m) -> std::size_t {
            return resolve_primary_index(bodies, seg, t_s, pos_m, sc_lookup, auto_primary_at);
        };

        auto apply_impulses_at = [&](const double t_s, const Vec3 &pos_m, const Vec3 &vel_mps) -> Vec3 {
            Vec3 dv_i{0.0, 0.0, 0.0};
            if (plan.impulses.empty())
            {
                return dv_i;
            }
            if (bodies.empty())
            {
                return dv_i;
            }

            // Since time steps are constructed from impulse boundaries, an absolute-equality check is usually OK,
            // but we allow a tiny tolerance to avoid missing due to floating-point addition.
            const double eps_t = 1e-9;
            for (const auto &imp: plan.impulses)
            {
                if (!segment_applies_to_spacecraft(imp, spacecraft_id))
                {
                    continue;
                }
                if (!(std::abs(imp.t_s - t_s) <= eps_t))
                {
                    continue;
                }
                const std::size_t primary = resolve_primary(imp, t_s, pos_m);
                dv_i += rtn_vector_to_inertial(eph, bodies, primary, t_s, pos_m, vel_mps, imp.dv_rtn_mps, imp.rtn_frame, sc_lookup);
            }
            return dv_i;
        };

        while (remaining > 0.0)
        {
            // Apply any instantaneous impulses at the current time.
            {
                const Vec3 dv_i = apply_impulses_at(t, y.position_m, y.velocity_mps);
                if (glm::dot(dv_i, dv_i) > 0.0 && std::isfinite(dv_i.x) && std::isfinite(dv_i.y) && std::isfinite(dv_i.z))
                {
                    y.velocity_mps += dv_i;
                }
            }

            const double burn_boundary = next_burn_boundary_after(plan, spacecraft_id, t, t_end_s);
            const double impulse_boundary = next_impulse_time_after(plan, spacecraft_id, t, t_end_s);
            const double boundary = std::min(burn_boundary, impulse_boundary);

            double h = std::min(remaining, boundary - t);
            if (!(h > 0.0) || !std::isfinite(h))
            {
                break;
            }

            const BurnSegment *seg = active_burn_at(plan, spacecraft_id, t);
            const bool has_engine = (seg != nullptr) && (seg->engine_index < out.engines.size());
            const Engine engine = has_engine ? out.engines[seg->engine_index] : Engine{};

            const double throttle = (seg != nullptr) ? clamp01(seg->throttle_0_1) : 0.0;
            const double min_thr = std::max(0.0, engine.min_throttle_0_1);
            const bool thrust_on = has_engine && (engine.max_thrust_N > 0.0) && (engine.isp_s > 0.0) &&
                                   (throttle >= min_thr) && (out.prop_mass_kg > 0.0) && !bodies.empty();

            const double thrust_N = thrust_on ? (engine.max_thrust_N * throttle) : 0.0;
            const double mdot_kgps = thrust_on ? (thrust_N / (engine.isp_s * kStandardGravity_mps2)) : 0.0;

            if (mdot_kgps > 0.0 && out.prop_mass_kg > 0.0)
            {
                const double h_burn = out.prop_mass_kg / mdot_kgps;
                if (h_burn > 0.0 && h_burn < h)
                {
                    // Burn until prop exhausted; the loop will then coast for the remaining time.
                    h = h_burn;
                }
            }

            const double t_segment_start = t;
            const double prop_start_kg = out.prop_mass_kg;
            const std::size_t primary = (seg != nullptr) ? resolve_primary(*seg, t, y.position_m) : 0;
            const Vec3 dir_rtn = (seg != nullptr) ? seg->dir_rtn_unit : Vec3{0.0, 0.0, 0.0};
            const TrajectoryFrameSpec rtn_frame = (seg != nullptr) ? seg->rtn_frame : TrajectoryFrameSpec{};

            auto accel = [&](double t_eval_s, const Vec3 &pos_m, const Vec3 &vel_mps) -> Vec3 {
                Vec3 a = gravity_accel_mps2(eph, bodies, gravitational_constant, softening_length_m, t_eval_s, pos_m);
                if (!(thrust_N > 0.0))
                {
                    return a;
                }

                const double dt_prop = t_eval_s - t_segment_start;
                double prop_kg = prop_start_kg - mdot_kgps * dt_prop;
                if (prop_kg <= 0.0)
                {
                    return a;
                }

                const double mass_kg = out.dry_mass_kg + prop_kg;
                if (!(mass_kg > 0.0) || !std::isfinite(mass_kg))
                {
                    return a;
                }

                const Vec3 dir_i = burn_dir_inertial_unit(eph, bodies, primary, t_eval_s, pos_m, vel_mps, dir_rtn, rtn_frame, sc_lookup);
                const double dir2 = glm::dot(dir_i, dir_i);
                if (!(dir2 > 0.0) || !std::isfinite(dir2))
                {
                    return a;
                }

                a += (thrust_N / mass_kg) * dir_i;
                return a;
            };

            y = dopri5_integrate_interval(t, y, h, accel, spacecraft_integrator);

            if (mdot_kgps > 0.0)
            {
                out.prop_mass_kg = std::max(0.0, out.prop_mass_kg - mdot_kgps * h);
            }

            out.state.position_m = y.position_m;
            out.state.velocity_mps = y.velocity_mps;
            advance_spin(out.state.spin, h);

            t += h;
            remaining -= h;
        }

        return out;
    }

} // namespace orbitsim::detail
