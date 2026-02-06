#pragma once

#include "orbitsim/types.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace orbitsim
{

    inline double wrap_angle_0_2pi(const double rad)
    {
        if (!std::isfinite(rad))
        {
            return 0.0;
        }
        const double two_pi = 2.0 * std::acos(-1.0);
        double x = std::fmod(rad, two_pi);
        if (x < 0.0)
        {
            x += two_pi;
        }
        return x;
    }

    /**
     * @brief Select the body index with highest gravitational acceleration at query_pos_m.
     *
     * Ranks bodies by mass_kg / r^2 (G is a constant positive multiplier and does not affect ranking).
     *
     * @tparam BodyPositionAt Callable: (std::size_t index) -> Vec3 returning body position.
     * @param bodies Massive bodies to consider
     * @param query_pos_m Query position (e.g. spacecraft position)
     * @param body_position_at Functor returning the position of body at given index
     * @param softening_length_m Optional softening length added in quadrature to distance
     * @return Index of the body with the strongest gravitational pull at query_pos_m
     */
    template<class BodyPositionAt>
    inline std::size_t auto_select_primary_index(const std::vector<MassiveBody> &bodies, const Vec3 &query_pos_m,
                                                  BodyPositionAt body_position_at,
                                                  const double softening_length_m = 0.0)
    {
        if (bodies.empty())
        {
            return 0;
        }

        const double eps2 = softening_length_m * softening_length_m;
        std::size_t best = 0;
        double best_metric = -1.0;
        for (std::size_t i = 0; i < bodies.size(); ++i)
        {
            const double mass_kg = std::isfinite(bodies[i].mass_kg) ? bodies[i].mass_kg : 0.0;
            if (!(mass_kg >= 0.0))
            {
                continue;
            }

            const Vec3 dr = body_position_at(i) - query_pos_m;
            const double r2 = glm::dot(dr, dr) + eps2;
            if (!(r2 > 0.0) || !std::isfinite(r2))
            {
                continue;
            }
            const double metric = mass_kg / r2;
            if (metric > best_metric)
            {
                best_metric = metric;
                best = i;
            }
        }
        return best;
    }

} // namespace orbitsim
