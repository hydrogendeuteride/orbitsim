#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/game_sim.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
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

    namespace soi_detail
    {
        inline bool finite3_(const Vec3 &v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

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

} // namespace orbitsim
