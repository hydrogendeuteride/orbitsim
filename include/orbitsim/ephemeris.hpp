#pragma once

#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <cstddef>
#include <vector>

namespace orbitsim
{
    class CelestialEphemerisSegment
    {
    public:
        double t0_s{0.0};
        double dt_s{0.0};
        std::vector<State> start{};
        std::vector<State> end{};

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

            const double tau = (t_s - t0_s) / dt_s;
            const State &s0 = start[body_index];
            const State &s1 = end[body_index];

            out.position_m =
                    hermite_position(s0.position_m, s0.velocity_mps, s1.position_m, s1.velocity_mps, dt_s, tau);
            out.velocity_mps =
                    hermite_velocity_mps(s0.position_m, s0.velocity_mps, s1.position_m, s1.velocity_mps, dt_s, tau);

            out.spin.axis = normalized_or(s0.spin.axis, Vec3{0.0, 1.0, 0.0});
            out.spin.rate_rad_per_s = s0.spin.rate_rad_per_s;
            out.spin.angle_rad = s0.spin.angle_rad + s0.spin.rate_rad_per_s * (t_s - t0_s);
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
} // namespace orbitsim
