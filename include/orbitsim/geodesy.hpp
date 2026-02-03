#pragma once

#include "orbitsim/coordinate_frames.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <cmath>
#include <optional>

namespace orbitsim
{

    /**
     * @brief Spherical geodetic coordinates on a rotating body.
     *
     * - `latitude_rad`: [-pi/2, +pi/2]
     * - `longitude_rad`: [-pi, +pi]
     * - `altitude_m`: above `body_radius_m`
     */
    struct GeodeticCoord
    {
        double latitude_rad{0.0};
        double longitude_rad{0.0};
        double altitude_m{0.0};
    };

    /**
     * @brief Convert an inertial position to spherical geodetic coordinates on a rotating body.
     *
     * Uses the provided body-fixed frame to compute latitude/longitude from the body-fixed position,
     * and altitude as `|r_bf| - body_radius_m`.
     *
     * @param body_fixed_frame Body-fixed rotating frame (ECEF-like).
     * @param pos_inertial_m Position in inertial coordinates [m].
     * @param body_radius_m Spherical body radius [m].
     * @return Geodetic coordinate if inputs are valid; `std::nullopt` otherwise.
     */
    inline std::optional<GeodeticCoord> geodetic_from_inertial(const RotatingFrame &body_fixed_frame,
                                                               const Vec3 &pos_inertial_m,
                                                               const double body_radius_m)
    {
        if (!(body_radius_m > 0.0) || !std::isfinite(body_radius_m))
        {
            return std::nullopt;
        }
        if (!body_fixed_frame.valid())
        {
            return std::nullopt;
        }

        const Vec3 r_bf_m = inertial_position_to_frame(body_fixed_frame, pos_inertial_m);
        const double r = glm::length(r_bf_m);
        if (!(r > 0.0) || !std::isfinite(r))
        {
            return std::nullopt;
        }

        GeodeticCoord out;
        out.latitude_rad = std::asin(std::clamp(r_bf_m.z / r, -1.0, 1.0));
        out.longitude_rad = std::atan2(r_bf_m.y, r_bf_m.x);
        out.altitude_m = r - body_radius_m;
        return out;
    }

    /**
     * @brief Convert spherical geodetic coordinates to an inertial position.
     *
     * Computes body-fixed Cartesian position from (lat, lon, radius+alt) and transforms it into inertial
     * using `body_fixed_frame`.
     *
     * @return Inertial position [m], or `{0,0,0}` if inputs are invalid.
     */
    inline Vec3 inertial_position_from_geodetic(const RotatingFrame &body_fixed_frame,
                                                const GeodeticCoord &coord,
                                                const double body_radius_m)
    {
        if (!(body_radius_m > 0.0) || !std::isfinite(body_radius_m) || !body_fixed_frame.valid())
        {
            return Vec3{0.0, 0.0, 0.0};
        }

        const double r = body_radius_m + coord.altitude_m;
        if (!(r > 0.0) || !std::isfinite(r))
        {
            return Vec3{0.0, 0.0, 0.0};
        }

        const double clat = std::cos(coord.latitude_rad);
        const double slat = std::sin(coord.latitude_rad);
        const double clon = std::cos(coord.longitude_rad);
        const double slon = std::sin(coord.longitude_rad);

        const Vec3 r_bf_m{r * clat * clon, r * clat * slon, r * slat};
        return frame_position_to_inertial(body_fixed_frame, r_bf_m);
    }

    /**
     * @brief Inertial velocity of a surface-fixed point on a rotating body.
     *
     * Treats the point as rigidly attached to the rotating frame:
     * `v = v_origin + ω × (r - r_origin)`.
     *
     * @return Inertial velocity [m/s], or `{0,0,0}` if `body_fixed_frame` is invalid.
     */
    inline Vec3 inertial_velocity_of_fixed_point(const RotatingFrame &body_fixed_frame, const Vec3 &pos_inertial_m)
    {
        if (!body_fixed_frame.valid())
        {
            return Vec3{0.0, 0.0, 0.0};
        }
        const Vec3 r_rel_m = pos_inertial_m - body_fixed_frame.origin_position_m;
        return body_fixed_frame.origin_velocity_mps + glm::cross(body_fixed_frame.omega_inertial_radps, r_rel_m);
    }

    /**
     * @brief Construct a local NED (north-east-down) frame at a geodetic location on the body's surface.
     *
     * Returned as a RotatingFrame:
     * - origin: inertial position/velocity of the surface point
     * - axes: `ex_i = North`, `ey_i = East`, `ez_i = Down` (all expressed in inertial coordinates)
     * - `omega_inertial_radps` matches the body's spin so `inertial_state_to_frame` accounts for non-inertial terms
     *
     * @return NED frame if inputs are valid; `std::nullopt` otherwise.
     */
    inline std::optional<RotatingFrame> make_ned_frame(const RotatingFrame &body_fixed_frame,
                                                       const GeodeticCoord &coord,
                                                       const double body_radius_m)
    {
        if (!body_fixed_frame.valid())
        {
            return std::nullopt;
        }
        if (!(body_radius_m > 0.0) || !std::isfinite(body_radius_m))
        {
            return std::nullopt;
        }

        // Basis in body-fixed coordinates.
        const double clat = std::cos(coord.latitude_rad);
        const double slat = std::sin(coord.latitude_rad);
        const double clon = std::cos(coord.longitude_rad);
        const double slon = std::sin(coord.longitude_rad);

        const Vec3 up_bf{clat * clon, clat * slon, slat};
        const Vec3 east_bf{-slon, clon, 0.0};
        const Vec3 north_bf{-slat * clon, -slat * slon, clat};
        const Vec3 down_bf = -up_bf;

        RotatingFrame out;
        out.origin_position_m = inertial_position_from_geodetic(body_fixed_frame, coord, body_radius_m);
        out.origin_velocity_mps = inertial_velocity_of_fixed_point(body_fixed_frame, out.origin_position_m);

        out.ex_i = normalized_or(frame_vector_to_inertial(body_fixed_frame, north_bf), Vec3{0.0, 1.0, 0.0});
        out.ey_i = normalized_or(frame_vector_to_inertial(body_fixed_frame, east_bf), Vec3{1.0, 0.0, 0.0});
        out.ez_i = normalized_or(frame_vector_to_inertial(body_fixed_frame, down_bf), Vec3{0.0, 0.0, 1.0});

        out.omega_inertial_radps = body_fixed_frame.omega_inertial_radps;
        return out;
    }

} // namespace orbitsim
