#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/game_sim.hpp"
#include "orbitsim/kepler_trajectory.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>

namespace orbitsim
{

    /**
     * @brief SOI-based primary selection options for rails/patched-conics mode.
     *
     * In rails mode you typically want:
     * - Switch into smaller/local SOIs immediately (e.g., planet -> moon).
     * - Switch out of the current SOI with hysteresis to avoid thrashing near the boundary.
     */
    struct SoiSwitchOptions
    {
        /// Entry test uses `distance <= enter_scale * soi_radius`.
        double enter_scale{1.0};
        /// Exit hysteresis uses `distance <= exit_scale * soi_radius` to keep current primary.
        /// Should be >= enter_scale.
        double exit_scale{1.02};

        /// If multiple SOIs contain the spacecraft, prefer the smallest SOI (most local body).
        bool prefer_smallest_soi{true};

        /// If no SOI contains the spacecraft (or SOIs are not configured), fall back to max-accel selection.
        bool fallback_to_max_accel{true};
    };

    struct SoiTransitionSearchOptions
    {
        double max_step_s{300.0};
        double refine_tolerance_s{0.25};
        std::size_t max_steps{10000};
        SoiSwitchOptions switch_options{};
        KeplerPropagationOptions propagation{};
    };

    struct SoiTransitionSearchResult
    {
        bool found{false};
        BodyId from_primary_body_id{kInvalidBodyId};
        BodyId to_primary_body_id{kInvalidBodyId};
        double t_s{std::numeric_limits<double>::quiet_NaN()};
        std::size_t tested_samples{0};
        KeplerStatus first_failure{KeplerStatus::Ok};
    };

    enum class SoiRadiusModel : std::uint8_t
    {
        Laplace = 0,
        Hill,
    };

    /**
     * @brief Compute the classical Laplace sphere-of-influence radius.
     *
     * The mass ratio is dimensionless, so callers may also pass gravitational
     * parameters (mu) instead of masses as long as both arguments use the same convention.
     *
     * Formula:
     *   r_soi = a * pow(m_child / m_parent, 2/5)
     */
    inline std::optional<double>
    compute_laplace_soi_radius(const double child_mass_kg, const double parent_mass_kg,
                               const double semi_major_axis_m)
    {
        if (!(child_mass_kg > 0.0) || !std::isfinite(child_mass_kg) || !(parent_mass_kg > 0.0) ||
            !std::isfinite(parent_mass_kg) || !(semi_major_axis_m > 0.0) || !std::isfinite(semi_major_axis_m))
        {
            return std::nullopt;
        }

        const double mass_ratio = child_mass_kg / parent_mass_kg;
        if (!(mass_ratio > 0.0) || !std::isfinite(mass_ratio))
        {
            return std::nullopt;
        }

        const double soi_radius_m = semi_major_axis_m * std::pow(mass_ratio, 2.0 / 5.0);
        if (!(soi_radius_m > 0.0) || !std::isfinite(soi_radius_m))
        {
            return std::nullopt;
        }
        return soi_radius_m;
    }

    /**
     * @brief Compute the Hill sphere radius using periapsis distance.
     *
     * Formula:
     *   r_hill = a * (1 - e) * cbrt(m_child / (3 * m_parent))
     */
    inline std::optional<double>
    compute_hill_soi_radius(const double child_mass_kg, const double parent_mass_kg, const double semi_major_axis_m,
                            const double eccentricity = 0.0)
    {
        if (!(child_mass_kg > 0.0) || !std::isfinite(child_mass_kg) || !(parent_mass_kg > 0.0) ||
            !std::isfinite(parent_mass_kg) || !(semi_major_axis_m > 0.0) || !std::isfinite(semi_major_axis_m) ||
            !(eccentricity >= 0.0) || !std::isfinite(eccentricity) || !(eccentricity < 1.0))
        {
            return std::nullopt;
        }

        const double mass_ratio = child_mass_kg / (3.0 * parent_mass_kg);
        if (!(mass_ratio > 0.0) || !std::isfinite(mass_ratio))
        {
            return std::nullopt;
        }

        const double periapsis_distance_m = semi_major_axis_m * (1.0 - eccentricity);
        if (!(periapsis_distance_m > 0.0) || !std::isfinite(periapsis_distance_m))
        {
            return std::nullopt;
        }

        const double soi_radius_m = periapsis_distance_m * std::cbrt(mass_ratio);
        if (!(soi_radius_m > 0.0) || !std::isfinite(soi_radius_m))
        {
            return std::nullopt;
        }
        return soi_radius_m;
    }

    namespace soi_detail
    {
        using detail::finite3_;

        inline bool near_equal_(const double a, const double b, const double rel_tol = 1e-12,
                                const double abs_tol = 1e-12)
        {
            const double scale = std::max(std::abs(a), std::abs(b));
            return std::abs(a - b) <= std::max(abs_tol, rel_tol * scale);
        }

        inline bool body_position_at_(const GameSimulation &sim, const CelestialEphemeris &eph, const BodyId id,
                                      const double t_s, Vec3 *out_pos_m)
        {
            if (out_pos_m == nullptr)
            {
                return false;
            }
            *out_pos_m = Vec3{0.0, 0.0, 0.0};
            const MassiveBody *b = sim.body_by_id(id);
            if (b == nullptr)
            {
                return false;
            }
            if (!eph.empty())
            {
                *out_pos_m = eph.body_position_at_by_id(id, t_s);
                return finite3_(*out_pos_m);
            }
            *out_pos_m = b->state.position_m;
            return finite3_(*out_pos_m);
        }

        inline std::optional<BodyId>
        select_body_in_soi_(const GameSimulation &sim, const CelestialEphemeris &eph, const Vec3 &sc_pos_m,
                            const double t_s, const SoiSwitchOptions &opt)
        {
            if (!finite3_(sc_pos_m))
            {
                return std::nullopt;
            }

            const double enter_scale = (std::isfinite(opt.enter_scale) ? std::max(0.0, opt.enter_scale) : 0.0);

            std::optional<BodyId> best_id;
            double best_soi = 0.0;
            double best_d = 0.0;
            double best_ratio = 0.0;

            for (const auto &b: sim.massive_bodies())
            {
                if (!(b.soi_radius_m > 0.0) || !std::isfinite(b.soi_radius_m))
                {
                    continue;
                }

                Vec3 bp;
                if (!body_position_at_(sim, eph, b.id, t_s, &bp))
                {
                    continue;
                }

                const Vec3 dr = sc_pos_m - bp;
                const double d2 = glm::dot(dr, dr);
                if (!(d2 >= 0.0) || !std::isfinite(d2))
                {
                    continue;
                }
                const double d = std::sqrt(d2);
                const double threshold = enter_scale * b.soi_radius_m;
                if (!(d <= threshold))
                {
                    continue;
                }

                const double ratio = (b.soi_radius_m > 0.0) ? (d / b.soi_radius_m) : std::numeric_limits<double>::infinity();

                if (!best_id.has_value())
                {
                    best_id = b.id;
                    best_soi = b.soi_radius_m;
                    best_d = d;
                    best_ratio = ratio;
                    continue;
                }

                if (opt.prefer_smallest_soi)
                {
                    if (b.soi_radius_m < best_soi || (near_equal_(b.soi_radius_m, best_soi) && d < best_d))
                    {
                        best_id = b.id;
                        best_soi = b.soi_radius_m;
                        best_d = d;
                        best_ratio = ratio;
                    }
                }
                else
                {
                    if (ratio < best_ratio || (near_equal_(ratio, best_ratio) && d < best_d))
                    {
                        best_id = b.id;
                        best_soi = b.soi_radius_m;
                        best_d = d;
                        best_ratio = ratio;
                    }
                }
            }

            return best_id;
        }

        inline std::optional<double>
        soi_radius_for_id_(const GameSimulation &sim, const BodyId id)
        {
            const MassiveBody *b = sim.body_by_id(id);
            if (b == nullptr)
            {
                return std::nullopt;
            }
            if (!(b->soi_radius_m > 0.0) || !std::isfinite(b->soi_radius_m))
            {
                return std::nullopt;
            }
            return b->soi_radius_m;
        }

        inline bool has_other_soi_candidate_(const GameSimulation &sim, const BodyId current_primary_body_id)
        {
            for (const MassiveBody &body: sim.massive_bodies())
            {
                if (body.id != current_primary_body_id &&
                    body.soi_radius_m > 0.0 &&
                    std::isfinite(body.soi_radius_m))
                {
                    return true;
                }
            }
            return false;
        }

        inline std::optional<double>
        distance_to_body_m_(const GameSimulation &sim, const CelestialEphemeris &eph, const Vec3 &sc_pos_m,
                            const double t_s, const BodyId body_id)
        {
            Vec3 bp;
            if (!body_position_at_(sim, eph, body_id, t_s, &bp))
            {
                return std::nullopt;
            }
            const Vec3 dr = sc_pos_m - bp;
            const double d2 = glm::dot(dr, dr);
            if (!(d2 >= 0.0) || !std::isfinite(d2))
            {
                return std::nullopt;
            }
            return std::sqrt(d2);
        }

        inline std::optional<BodyId>
        select_primary_by_max_accel_(const GameSimulation &sim, const CelestialEphemeris &eph, const Vec3 &sc_pos_m,
                                     const double t_s)
        {
            if (sim.massive_bodies().empty())
            {
                return std::nullopt;
            }
            const double G = sim.config().gravitational_constant;
            if (!(G > 0.0) || !std::isfinite(G))
            {
                return sim.massive_bodies().front().id;
            }

            const double eps2 = sim.config().softening_length_m * sim.config().softening_length_m;

            std::optional<BodyId> best_id;
            double best_a = -1.0;
            for (const auto &b: sim.massive_bodies())
            {
                if (!(b.mass_kg > 0.0) || !std::isfinite(b.mass_kg))
                {
                    continue;
                }
                Vec3 bp;
                if (!body_position_at_(sim, eph, b.id, t_s, &bp))
                {
                    continue;
                }

                const Vec3 dr = bp - sc_pos_m;
                const double r2 = glm::dot(dr, dr) + eps2;
                if (!(r2 > 0.0) || !std::isfinite(r2))
                {
                    continue;
                }
                const double amag = (G * b.mass_kg) / r2;
                if (amag > best_a)
                {
                    best_a = amag;
                    best_id = b.id;
                }
            }
            return best_id;
        }
    } // namespace soi_detail

    /**
     * @brief Select/maintain the current patched-conics primary body during rails timewarp.
     *
     * Rules:
     * - If within the SOI of any body, choose the most local one (smallest SOI by default).
     * - If already inside a small/local SOI, keep it until leaving by `exit_scale` (hysteresis).
     * - If no SOI is configured/contains the spacecraft, optionally fall back to "max accel" selection.
     *
     * @param sim Game simulation providing body definitions and configuration.
     * @param eph Ephemeris for body motion; if empty, uses current body states in sim.
     * @param sc_pos_m Spacecraft inertial position at time t_s.
     * @param t_s Time for ephemeris lookup.
     * @param current_primary_body_id Current primary body (may be invalid).
     * @param opt Switching options (hysteresis and tie-breaking).
     * @return Selected primary body id (or invalid if none).
     */
    inline BodyId select_primary_body_id_rails(const GameSimulation &sim, const CelestialEphemeris &eph,
                                               const Vec3 &sc_pos_m, const double t_s,
                                               const BodyId current_primary_body_id, const SoiSwitchOptions &opt = {})
    {
        const std::optional<BodyId> best_enter = soi_detail::select_body_in_soi_(sim, eph, sc_pos_m, t_s, opt);
        if (best_enter.has_value())
        {
            if (current_primary_body_id == *best_enter)
            {
                return *best_enter;
            }

            const std::optional<double> best_soi = soi_detail::soi_radius_for_id_(sim, *best_enter);
            const std::optional<double> cur_soi = soi_detail::soi_radius_for_id_(sim, current_primary_body_id);

            // If we can switch into a more local SOI, do it immediately (planet -> moon).
            if (best_soi.has_value() && cur_soi.has_value() && (*best_soi < *cur_soi))
            {
                return *best_enter;
            }
            if (!cur_soi.has_value())
            {
                return *best_enter;
            }

            // Otherwise keep the current primary while still within exit hysteresis, else move to best_enter.
            const double exit_scale = (std::isfinite(opt.exit_scale) ? std::max(0.0, opt.exit_scale) : 0.0);
            const std::optional<double> dcur =
                    soi_detail::distance_to_body_m_(sim, eph, sc_pos_m, t_s, current_primary_body_id);
            if (dcur.has_value() && (*dcur <= exit_scale * (*cur_soi)))
            {
                return current_primary_body_id;
            }
            return *best_enter;
        }

        // No SOI contains the spacecraft: keep current if still within exit hysteresis, else fall back.
        const std::optional<double> cur_soi = soi_detail::soi_radius_for_id_(sim, current_primary_body_id);
        if (cur_soi.has_value())
        {
            const double exit_scale = (std::isfinite(opt.exit_scale) ? std::max(0.0, opt.exit_scale) : 0.0);
            const std::optional<double> dcur =
                    soi_detail::distance_to_body_m_(sim, eph, sc_pos_m, t_s, current_primary_body_id);
            if (dcur.has_value() && (*dcur <= exit_scale * (*cur_soi)))
            {
                return current_primary_body_id;
            }
        }

        if (opt.fallback_to_max_accel)
        {
            const std::optional<BodyId> best = soi_detail::select_primary_by_max_accel_(sim, eph, sc_pos_m, t_s);
            return best.value_or(kInvalidBodyId);
        }
        return kInvalidBodyId;
    }

    inline BodyId select_primary_body_id_rails(const GameSimulation &sim, const CelestialEphemeris &eph,
                                               const SpacecraftId spacecraft_id, const double t_s,
                                               const BodyId current_primary_body_id, const SoiSwitchOptions &opt = {})
    {
        const Spacecraft *sc = sim.spacecraft_by_id(spacecraft_id);
        if (sc == nullptr)
        {
            return kInvalidBodyId;
        }
        return select_primary_body_id_rails(sim, eph, sc->state.position_m, t_s, current_primary_body_id, opt);
    }

    inline std::optional<Vec3> sample_kepler_arc_inertial_position_at(const GameSimulation &sim,
                                                                      const CelestialEphemeris &eph,
                                                                      const KeplerArc &arc,
                                                                      const double t_s,
                                                                      const KeplerPropagationOptions &propagation,
                                                                      KeplerStatus *out_status = nullptr)
    {
        if (out_status)
        {
            *out_status = KeplerStatus::Ok;
        }

        const KeplerArcSample sample = sample_kepler_arc_state(arc, t_s, propagation);
        if (!sample.ok())
        {
            if (out_status)
            {
                *out_status = sample.diagnostics.status;
            }
            return std::nullopt;
        }

        Vec3 primary_pos_m{0.0, 0.0, 0.0};
        if (!soi_detail::body_position_at_(sim, eph, arc.primary_body_id, t_s, &primary_pos_m))
        {
            if (out_status)
            {
                *out_status = KeplerStatus::InvalidFinalState;
            }
            return std::nullopt;
        }

        const Vec3 inertial_position_m = primary_pos_m + sample.state_relative.position_m;
        if (!detail::finite3_(inertial_position_m))
        {
            if (out_status)
            {
                *out_status = KeplerStatus::InvalidFinalState;
            }
            return std::nullopt;
        }
        return inertial_position_m;
    }

    inline SoiTransitionSearchResult find_next_soi_transition_on_kepler_arc(const GameSimulation &sim,
                                                                           const CelestialEphemeris &eph,
                                                                           const KeplerArc &arc,
                                                                           const BodyId current_primary_body_id,
                                                                           const double t_limit_s,
                                                                           const SoiTransitionSearchOptions &options = {})
    {
        SoiTransitionSearchResult out{};
        out.from_primary_body_id = current_primary_body_id;
        if (!kepler_arc_valid(arc) ||
            current_primary_body_id == kInvalidBodyId ||
            !std::isfinite(t_limit_s))
        {
            out.first_failure = KeplerStatus::InvalidInitialState;
            return out;
        }

        const double arc_span_s = arc.t1_s - arc.t0_s;
        const double requested_span_s = t_limit_s - arc.t0_s;
        if (arc_span_s == 0.0 ||
            requested_span_s == 0.0 ||
            ((arc_span_s > 0.0) != (requested_span_s > 0.0)))
        {
            return out;
        }
        if (!soi_detail::has_other_soi_candidate_(sim, current_primary_body_id))
        {
            return out;
        }

        const double dir = arc_span_s > 0.0 ? 1.0 : -1.0;
        const double max_abs_t_s = dir > 0.0 ? std::min(arc.t1_s, t_limit_s) : std::max(arc.t1_s, t_limit_s);
        const double abs_duration_s = std::abs(max_abs_t_s - arc.t0_s);
        if (!(abs_duration_s > 0.0) || !std::isfinite(abs_duration_s))
        {
            return out;
        }

        const double step_s = std::clamp(std::abs(options.max_step_s), 1.0e-6, abs_duration_s);
        const double tolerance_s = std::clamp(std::abs(options.refine_tolerance_s), 1.0e-9, step_s);
        const std::size_t max_steps = std::max<std::size_t>(1u, options.max_steps);

        SoiSwitchOptions transition_switch = options.switch_options;
        transition_switch.fallback_to_max_accel = false;

        auto selected_primary_at = [&](const double t_s, KeplerStatus &status) -> BodyId {
            const std::optional<Vec3> inertial_position_m =
                    sample_kepler_arc_inertial_position_at(sim,
                                                           eph,
                                                           arc,
                                                           t_s,
                                                           options.propagation,
                                                           &status);
            if (!inertial_position_m.has_value())
            {
                return kInvalidBodyId;
            }
            return select_primary_body_id_rails(sim,
                                                eph,
                                                *inertial_position_m,
                                                t_s,
                                                current_primary_body_id,
                                                transition_switch);
        };

        double prev_t_s = arc.t0_s;
        KeplerStatus prev_status = KeplerStatus::Ok;
        (void) selected_primary_at(prev_t_s, prev_status);
        if (prev_status != KeplerStatus::Ok)
        {
            out.first_failure = prev_status;
            return out;
        }

        for (std::size_t step_index = 1u; step_index <= max_steps; ++step_index)
        {
            const double offset_s = std::min(abs_duration_s, step_s * static_cast<double>(step_index));
            const double t_s = arc.t0_s + dir * offset_s;
            KeplerStatus sample_status = KeplerStatus::Ok;
            const BodyId selected = selected_primary_at(t_s, sample_status);
            ++out.tested_samples;
            if (sample_status != KeplerStatus::Ok)
            {
                out.first_failure = sample_status;
                return out;
            }

            if (selected != kInvalidBodyId && selected != current_primary_body_id)
            {
                double low_t_s = prev_t_s;
                double high_t_s = t_s;
                BodyId high_primary = selected;
                for (int i = 0; i < 96 && std::abs(high_t_s - low_t_s) > tolerance_s; ++i)
                {
                    const double mid_t_s = 0.5 * (low_t_s + high_t_s);
                    KeplerStatus mid_status = KeplerStatus::Ok;
                    const BodyId mid_primary = selected_primary_at(mid_t_s, mid_status);
                    ++out.tested_samples;
                    if (mid_status != KeplerStatus::Ok)
                    {
                        out.first_failure = mid_status;
                        return out;
                    }
                    if (mid_primary != kInvalidBodyId && mid_primary != current_primary_body_id)
                    {
                        high_t_s = mid_t_s;
                        high_primary = mid_primary;
                    }
                    else
                    {
                        low_t_s = mid_t_s;
                    }
                }

                out.found = true;
                out.to_primary_body_id = high_primary;
                out.t_s = high_t_s;
                return out;
            }

            if (offset_s >= abs_duration_s)
            {
                break;
            }
            prev_t_s = t_s;
        }

        return out;
    }

} // namespace orbitsim
