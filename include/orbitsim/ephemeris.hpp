#pragma once

#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <unordered_map>
#include <vector>

namespace orbitsim
{
    /**
     * @brief Single time-segment of precomputed massive body states.
     *
     * Stores start and end states for all bodies over interval [t0_s, t0_s + dt_s].
     * Interpolates using cubic Hermite to get smooth position/velocity at any time.
     */
    class CelestialEphemerisSegment
    {
    public:
        double t0_s{0.0};            ///< Start time of this segment
        double dt_s{0.0};            ///< Duration of this segment
        std::vector<State> start{};  ///< Body states at t0_s
        std::vector<State> end{};    ///< Body states at t0_s + dt_s

        /**
         * @brief Interpolate body state at arbitrary time within this segment.
         *
         * Uses Hermite interpolation for position/velocity, linear for spin angle.
         * Times outside [t0_s, t0_s + dt_s] are clamped to segment boundaries.
         */
        inline State body_state_at(const std::size_t body_index, const double t_s) const
        {
            State out{};
            if (body_index >= start.size() || body_index >= end.size())
            {
                return out;
            }
            if (!(dt_s != 0.0) || !std::isfinite(dt_s))
            {
                return start[body_index];
            }

            const double tau = clamp01((t_s - t0_s) / dt_s);
            const double t_eval_s = t0_s + tau * dt_s;
            const State &s0 = start[body_index];
            const State &s1 = end[body_index];

            out.position_m =
                    hermite_position(s0.position_m, s0.velocity_mps, s1.position_m, s1.velocity_mps, dt_s, tau);
            out.velocity_mps =
                    hermite_velocity_mps(s0.position_m, s0.velocity_mps, s1.position_m, s1.velocity_mps, dt_s, tau);

            out.spin.axis = normalized_or(s0.spin.axis, Vec3{0.0, 1.0, 0.0});
            out.spin.rate_rad_per_s = s0.spin.rate_rad_per_s;
            out.spin.angle_rad = s0.spin.angle_rad + s0.spin.rate_rad_per_s * (t_eval_s - t0_s);
            return out;
        }

        inline Vec3 body_position_at(const std::size_t body_index, const double t_s) const
        {
            return body_state_at(body_index, t_s).position_m;
        }

        inline Vec3 body_velocity_at(const std::size_t body_index, const double t_s) const
        {
            return body_state_at(body_index, t_s).velocity_mps;
        }
    };

    /**
     * @brief Piecewise ephemeris for all massive bodies over a time span.
     *
     * Contains multiple contiguous segments; binary search finds the correct
     * segment for a given query time. Used for trajectory prediction without
     * re-running full N-body simulation.
     */
    class CelestialEphemeris
    {
    public:
        std::vector<CelestialEphemerisSegment> segments{};            ///< Contiguous time segments
        std::vector<BodyId> body_ids{};                               ///< Body IDs in segment order
        std::unordered_map<BodyId, std::size_t> body_id_to_index{};   ///< Fast ID -> index lookup

        inline bool empty() const { return segments.empty(); }

        /** @brief Set body IDs and rebuild the ID-to-index lookup table. */
        inline void set_body_ids(std::vector<BodyId> ids)
        {
            body_ids = std::move(ids);
            body_id_to_index.clear();
            body_id_to_index.reserve(body_ids.size());
            for (std::size_t i = 0; i < body_ids.size(); ++i)
            {
                body_id_to_index[body_ids[i]] = i;
            }
        }

        /** @brief Look up array index for a body ID. Returns false if not found. */
        inline bool body_index_for_id(const BodyId id, std::size_t *out_index) const
        {
            if (out_index == nullptr)
            {
                return false;
            }
            *out_index = 0;
            auto it = body_id_to_index.find(id);
            if (it == body_id_to_index.end())
            {
                return false;
            }
            *out_index = it->second;
            return true;
        }

        inline double t0_s() const
        {
            if (segments.empty())
            {
                return 0.0;
            }
            return segments.front().t0_s;
        }

        inline double t_end_s() const
        {
            if (segments.empty())
            {
                return 0.0;
            }
            const auto &last = segments.back();
            return last.t0_s + last.dt_s;
        }

        /**
         * @brief Query body state at any time, finding correct segment via binary search.
         *
         * Times before first segment clamp to first; times after last clamp to last.
         */
        inline State body_state_at(const std::size_t body_index, const double t_s) const
        {
            State out{};
            if (segments.empty())
            {
                return out;
            }

            const auto it = std::upper_bound(
                    segments.begin(),
                    segments.end(),
                    t_s,
                    [](const double t, const CelestialEphemerisSegment &seg) { return t < seg.t0_s; });

            if (it == segments.begin())
            {
                return segments.front().body_state_at(body_index, t_s);
            }
            if (it == segments.end())
            {
                return segments.back().body_state_at(body_index, t_s);
            }
            return std::prev(it)->body_state_at(body_index, t_s);
        }

        inline Vec3 body_position_at(const std::size_t body_index, const double t_s) const
        {
            return body_state_at(body_index, t_s).position_m;
        }

        inline Vec3 body_velocity_at(const std::size_t body_index, const double t_s) const
        {
            return body_state_at(body_index, t_s).velocity_mps;
        }

        inline State body_state_at_by_id(const BodyId body_id, const double t_s) const
        {
            std::size_t index = 0;
            if (!body_index_for_id(body_id, &index))
            {
                return {};
            }
            return body_state_at(index, t_s);
        }

        inline Vec3 body_position_at_by_id(const BodyId body_id, const double t_s) const
        {
            return body_state_at_by_id(body_id, t_s).position_m;
        }

        inline Vec3 body_velocity_at_by_id(const BodyId body_id, const double t_s) const
        {
            return body_state_at_by_id(body_id, t_s).velocity_mps;
        }
    };
} // namespace orbitsim
