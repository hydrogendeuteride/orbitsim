#pragma once

#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <cmath>
#include <utility>

namespace orbitsim
{

    /// @brief Result of computing a circular orbit state relative to a central body.
    struct RelativeOrbitState
    {
        Vec3 position_m{0.0, 0.0, 0.0};
        Vec3 velocity_mps{0.0, 0.0, 0.0};
    };

    /// @brief Compute relative position/velocity for a circular orbit around a central body.
    /// @param central_mass_kg Mass of the central body [kg].
    /// @param orbital_radius_m Orbital radius (distance from center of central body) [m].
    /// @param inclination_rad Orbital inclination [rad]. 0 = equatorial, positive tilts velocity into +Z.
    /// @param arg_latitude_rad Argument of latitude (angle from +X axis in orbital plane) [rad].
    /// @return Position and velocity vectors relative to the central body.
    inline RelativeOrbitState circular_orbit_relative_state(const double central_mass_kg, const double orbital_radius_m,
                                                            const double inclination_rad = 0.0,
                                                            const double arg_latitude_rad = 0.0)
    {
        const double mu = kGravitationalConstant_SI * central_mass_kg;
        const double v_circ = std::sqrt(mu / orbital_radius_m);

        // Position in orbital plane (before inclination rotation):
        // At arg_latitude=0: position is along +X, velocity is along +Y.
        const double cos_u = std::cos(arg_latitude_rad);
        const double sin_u = std::sin(arg_latitude_rad);

        // Orbital plane position: r = R * [cos(u), sin(u), 0]
        // Orbital plane velocity: v = V_circ * [-sin(u), cos(u), 0]
        // Then rotate the Y-component by inclination into Y-Z plane.
        const double cos_i = std::cos(inclination_rad);
        const double sin_i = std::sin(inclination_rad);

        RelativeOrbitState out;
        out.position_m = Vec3{orbital_radius_m * cos_u, orbital_radius_m * sin_u * cos_i,
                              orbital_radius_m * sin_u * sin_i};
        out.velocity_mps = Vec3{-v_circ * sin_u, v_circ * cos_u * cos_i, v_circ * cos_u * sin_i};
        return out;
    }

    /// @brief Result of computing barycentric states for a two-body circular orbit.
    struct TwoBodyBarycentricStates
    {
        State state_a; // State of body A in barycentric frame.
        State state_b; // State of body B in barycentric frame.
    };

    /// @brief Compute barycentric states for two bodies in a mutual circular orbit.
    /// @param mass_a_kg Mass of body A [kg].
    /// @param mass_b_kg Mass of body B [kg].
    /// @param separation_m Distance between the two bodies [m].
    /// @param inclination_rad Orbital inclination of the relative orbit [rad].
    /// @param arg_latitude_rad Initial argument of latitude (position angle) [rad].
    /// @return States of both bodies in the barycentric reference frame.
    inline TwoBodyBarycentricStates two_body_circular_barycentric(const double mass_a_kg, const double mass_b_kg,
                                                                  const double separation_m,
                                                                  const double inclination_rad = 0.0,
                                                                  const double arg_latitude_rad = 0.0)
    {
        const double m_tot = mass_a_kg + mass_b_kg;
        if (!(m_tot > 0.0) || !std::isfinite(m_tot))
        {
            return {};
        }

        // Relative orbit (A -> B direction).
        const RelativeOrbitState rel = circular_orbit_relative_state(m_tot, separation_m, inclination_rad, arg_latitude_rad);

        // Barycentric positions: r_a = -(m_b/m_tot)*r_rel, r_b = (m_a/m_tot)*r_rel
        const double frac_a = mass_b_kg / m_tot;
        const double frac_b = mass_a_kg / m_tot;

        TwoBodyBarycentricStates out;
        out.state_a.position_m = -frac_a * rel.position_m;
        out.state_a.velocity_mps = -frac_a * rel.velocity_mps;
        out.state_b.position_m = frac_b * rel.position_m;
        out.state_b.velocity_mps = frac_b * rel.velocity_mps;
        return out;
    }

} // namespace orbitsim
