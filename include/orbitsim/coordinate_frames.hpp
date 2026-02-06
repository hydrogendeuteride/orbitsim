#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/frame_utils.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <cmath>
#include <optional>

namespace orbitsim
{

    /** @brief Radial-Tangential-Normal frame basis vectors. */
    struct RtnFrame
    {
        Vec3 R{1.0, 0.0, 0.0}; ///< Radial (away from primary)
        Vec3 T{0.0, 1.0, 0.0}; ///< Tangential (in orbital plane, roughly velocity direction)
        Vec3 N{0.0, 0.0, 1.0}; ///< Normal (angular momentum direction)
    };

    /**
     * @brief Compute RTN frame from relative position and velocity.
     *
     * R = r_hat, N = h_hat (angular momentum), T = N x R.
     * Handles degenerate cases (rectilinear, zero velocity) gracefully.
     */
    inline RtnFrame compute_rtn_frame(const Vec3 &r_rel_m, const Vec3 &v_rel_mps)
    {
        RtnFrame f;
        f.R = normalized_or(r_rel_m, Vec3{1.0, 0.0, 0.0});

        const Vec3 h = glm::cross(r_rel_m, v_rel_mps);
        const double h2 = glm::dot(h, h);
        if (h2 > 1e-24 && std::isfinite(h2))
        {
            f.N = h / std::sqrt(h2);
            f.T = normalized_or(glm::cross(f.N, f.R), Vec3{0.0, 1.0, 0.0});
            return f;
        }

        Vec3 vhat = normalized_or(v_rel_mps, Vec3{0.0, 1.0, 0.0});
        Vec3 t = vhat - glm::dot(vhat, f.R) * f.R;
        if (glm::dot(t, t) <= 1e-24 || !std::isfinite(glm::dot(t, t)))
        {
            const Vec3 axis = (std::abs(f.R.x) < 0.9) ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
            t = glm::cross(f.R, axis);
        }
        f.T = normalized_or(t, Vec3{0.0, 1.0, 0.0});
        f.N = normalized_or(glm::cross(f.R, f.T), Vec3{0.0, 0.0, 1.0});
        f.T = normalized_or(glm::cross(f.N, f.R), f.T);
        return f;
    }

    /**
     * @brief Rigid translating + rotating coordinate frame relative to the library inertial frame.
     *
     * Convention:
     * - Basis vectors `ex_i/ey_i/ez_i` are expressed in inertial coordinates (orthonormal).
     * - `omega_inertial_radps` is the frame angular velocity expressed in inertial coordinates.
     * - Inertial -> frame for velocities applies the non-inertial term: `v_f = R^T v_i - ω_f × r_f`.
     */
    struct RotatingFrame
    {
        /// Frame origin position expressed in inertial coordinates.
        Vec3 origin_position_m{0.0, 0.0, 0.0};
        /// Frame origin velocity expressed in inertial coordinates.
        Vec3 origin_velocity_mps{0.0, 0.0, 0.0};

        /// Frame X-axis expressed in inertial coordinates.
        Vec3 ex_i{1.0, 0.0, 0.0};
        /// Frame Y-axis expressed in inertial coordinates.
        Vec3 ey_i{0.0, 1.0, 0.0};
        /// Frame Z-axis expressed in inertial coordinates.
        Vec3 ez_i{0.0, 0.0, 1.0};

        /// Frame angular velocity expressed in inertial coordinates.
        Vec3 omega_inertial_radps{0.0, 0.0, 0.0};

        /// @brief Returns true if all frame fields are finite.
        inline bool valid() const
        {
            const auto finite3 = [](const Vec3 &v) {
                return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
            };

            return finite3(origin_position_m) && finite3(origin_velocity_mps) && finite3(ex_i) && finite3(ey_i) &&
                   finite3(ez_i) && finite3(omega_inertial_radps);
        }
    };

    inline Vec3 inertial_vector_to_frame(const RotatingFrame &frame, const Vec3 &v_in)
    {
        return Vec3{glm::dot(frame.ex_i, v_in), glm::dot(frame.ey_i, v_in), glm::dot(frame.ez_i, v_in)};
    }

    inline Vec3 frame_vector_to_inertial(const RotatingFrame &frame, const Vec3 &v_frame)
    {
        return frame.ex_i * v_frame.x + frame.ey_i * v_frame.y + frame.ez_i * v_frame.z;
    }

    inline Vec3 inertial_position_to_frame(const RotatingFrame &frame, const Vec3 &pos_in_m)
    {
        return inertial_vector_to_frame(frame, pos_in_m - frame.origin_position_m);
    }

    inline Vec3 frame_position_to_inertial(const RotatingFrame &frame, const Vec3 &pos_frame_m)
    {
        return frame.origin_position_m + frame_vector_to_inertial(frame, pos_frame_m);
    }

    /**
     * @brief Transform a full state from inertial to a rotating/translating frame.
     *
     * Applies translation/rotation for position, and applies the non-inertial velocity correction
     * `v_f = R^T (v_i - v_origin) - ω_f × r_f`.
     *
     * @param state_in State expressed in inertial coordinates.
     * @param frame Target frame definition (basis and angular velocity).
     * @return Equivalent state expressed in the target frame.
     */
    inline State inertial_state_to_frame(const State &state_in, const RotatingFrame &frame)
    {
        State out = state_in;

        const Vec3 r_in_rel_m = state_in.position_m - frame.origin_position_m;
        const Vec3 v_in_rel_mps = state_in.velocity_mps - frame.origin_velocity_mps;

        out.position_m = inertial_vector_to_frame(frame, r_in_rel_m);

        const Vec3 omega_frame_radps = inertial_vector_to_frame(frame, frame.omega_inertial_radps);
        out.velocity_mps =
                inertial_vector_to_frame(frame, v_in_rel_mps) - glm::cross(omega_frame_radps, out.position_m);

        return out;
    }

    /**
     * @brief Transform a full state from a rotating/translating frame to inertial.
     *
     * Applies translation/rotation for position, and inverts the non-inertial correction
     * `v_i = v_origin + R (v_f + ω_f × r_f)`.
     *
     * @param state_frame State expressed in the rotating frame.
     * @param frame Source frame definition (basis and angular velocity).
     * @return Equivalent state expressed in inertial coordinates.
     */
    inline State frame_state_to_inertial(const State &state_frame, const RotatingFrame &frame)
    {
        State out = state_frame;

        out.position_m = frame_position_to_inertial(frame, state_frame.position_m);

        const Vec3 omega_frame_radps = inertial_vector_to_frame(frame, frame.omega_inertial_radps);
        const Vec3 v_in_rel_mps = frame_vector_to_inertial(
                frame, state_frame.velocity_mps + glm::cross(omega_frame_radps, state_frame.position_m));
        out.velocity_mps = frame.origin_velocity_mps + v_in_rel_mps;

        return out;
    }

    // -------------------------------------------------------------------------
    // Frame-to-frame helpers (implemented via inertial as the intermediate).
    // -------------------------------------------------------------------------

    inline Vec3 frame_vector_to_frame(const RotatingFrame &frame_from, const RotatingFrame &frame_to,
                                      const Vec3 &v_from)
    {
        const Vec3 v_in = frame_vector_to_inertial(frame_from, v_from);
        return inertial_vector_to_frame(frame_to, v_in);
    }

    inline Vec3 frame_position_to_frame(const RotatingFrame &frame_from, const RotatingFrame &frame_to,
                                        const Vec3 &pos_from_m)
    {
        const Vec3 pos_in_m = frame_position_to_inertial(frame_from, pos_from_m);
        return inertial_position_to_frame(frame_to, pos_in_m);
    }

    inline State frame_state_to_frame(const State &state_from, const RotatingFrame &frame_from,
                                      const RotatingFrame &frame_to)
    {
        const State s_in = frame_state_to_inertial(state_from, frame_from);
        return inertial_state_to_frame(s_in, frame_to);
    }

    // -------------------------------------------------------------------------
    // ECI/BCI: Body-centered inertial frame (translation only; no rotation).
    // -------------------------------------------------------------------------

    inline RotatingFrame make_body_centered_inertial_frame(const State &body_inertial_state)
    {
        RotatingFrame f;
        f.origin_position_m = body_inertial_state.position_m;
        f.origin_velocity_mps = body_inertial_state.velocity_mps;
        f.ex_i = Vec3{1.0, 0.0, 0.0};
        f.ey_i = Vec3{0.0, 1.0, 0.0};
        f.ez_i = Vec3{0.0, 0.0, 1.0};
        f.omega_inertial_radps = Vec3{0.0, 0.0, 0.0};
        return f;
    }

    inline RotatingFrame make_body_centered_inertial_frame(const MassiveBody &body)
    {
        return make_body_centered_inertial_frame(body.state);
    }

    /// @brief Body-centered inertial frame at time `t_s` (uses ephemeris if available).
    inline RotatingFrame make_body_centered_inertial_frame_at(const CelestialEphemeris &eph, const MassiveBody &body,
                                                              const double t_s)
    {
        if (!eph.empty())
        {
            return make_body_centered_inertial_frame(eph.body_state_at_by_id(body.id, t_s));
        }
        return make_body_centered_inertial_frame(body);
    }

    // -------------------------------------------------------------------------
    // ECEF: Body-fixed rotating frame from SpinState (generic "planet-fixed").
    //
    // Since a full body orientation model isn't present yet, this defines a deterministic
    // body-fixed basis from (axis, angle):
    // - ez_i: spin axis (in inertial coordinates)
    // - ex_i: a "prime meridian" reference perpendicular to ez_i, rotated by angle about ez_i
    // - ey_i: completes right-handed basis
    // -------------------------------------------------------------------------

    /**
     * @brief Rotate vector about a unit axis by an angle (Rodrigues' rotation formula).
     *
     * @param v Vector to rotate.
     * @param axis_unit Rotation axis (must be unit length).
     * @param angle_rad Right-hand rotation angle about `axis_unit`.
     */
    inline Vec3 rotate_about_axis(const Vec3 &v, const Vec3 &axis_unit, const double angle_rad)
    {
        const double c = std::cos(angle_rad);
        const double s = std::sin(angle_rad);
        return v * c + glm::cross(axis_unit, v) * s + axis_unit * (glm::dot(axis_unit, v) * (1.0 - c));
    }

    /**
     * @brief Construct a deterministic body-fixed frame from the state's spin axis/angle/rate.
     *
     * The basis is constructed as:
     * - `ez_i`: spin axis (in inertial coords)
     * - `ex_i`: "prime meridian" reference perpendicular to `ez_i`, rotated by `angle_rad` about `ez_i`
     * - `ey_i`: completes a right-handed orthonormal basis
     *
     * Returns `std::nullopt` if the spin axis is not valid/finite or the basis construction becomes
     * degenerate (e.g., near-zero perpendicular component).
     */
    inline std::optional<RotatingFrame> make_body_fixed_frame(const State &body_inertial_state)
    {
        const Vec3 axis = normalized_or(body_inertial_state.spin.axis, Vec3{0.0, 0.0, 0.0});
        const double axis2 = glm::dot(axis, axis);
        if (!(axis2 > 0.0) || !std::isfinite(axis2))
        {
            return std::nullopt;
        }

        RotatingFrame f;
        f.origin_position_m = body_inertial_state.position_m;
        f.origin_velocity_mps = body_inertial_state.velocity_mps;

        const Vec3 ez = axis;

        // Choose a deterministic reference for the prime meridian at angle=0.
        const Vec3 a = (std::abs(ez.x) < 0.9) ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
        Vec3 ex0 = a - glm::dot(a, ez) * ez;
        ex0 = normalized_or(ex0, Vec3{0.0, 0.0, 0.0});
        const double ex02 = glm::dot(ex0, ex0);
        if (!(ex02 > 0.0) || !std::isfinite(ex02))
        {
            return std::nullopt;
        }

        const double theta = wrap_angle_0_2pi(body_inertial_state.spin.angle_rad);
        Vec3 ex = rotate_about_axis(ex0, ez, theta);
        ex = normalized_or(ex, ex0);

        Vec3 ey = glm::cross(ez, ex);
        ey = normalized_or(ey, Vec3{0.0, 0.0, 0.0});
        ex = normalized_or(glm::cross(ey, ez), ex);
        ey = normalized_or(glm::cross(ez, ex), ey);

        f.ex_i = ex;
        f.ey_i = ey;
        f.ez_i = ez;

        if (std::isfinite(body_inertial_state.spin.rate_rad_per_s))
        {
            f.omega_inertial_radps = ez * body_inertial_state.spin.rate_rad_per_s;
        }
        else
        {
            f.omega_inertial_radps = Vec3{0.0, 0.0, 0.0};
        }

        return f;
    }

    inline std::optional<RotatingFrame> make_body_fixed_frame(const MassiveBody &body)
    {
        return make_body_fixed_frame(body.state);
    }

    /// @brief Body-fixed frame at time `t_s` (uses ephemeris if available).
    inline std::optional<RotatingFrame> make_body_fixed_frame_at(const CelestialEphemeris &eph, const MassiveBody &body,
                                                                 const double t_s)
    {
        if (!eph.empty())
        {
            return make_body_fixed_frame(eph.body_state_at_by_id(body.id, t_s));
        }
        return make_body_fixed_frame(body);
    }

    // -------------------------------------------------------------------------
    // LVLH: Local-Vertical/Local-Horizontal frame centered on the spacecraft.
    //
    // - ex_i: radial (+R) from primary to spacecraft
    // - ez_i: orbit normal (+N) from (r×v)
    // - ey_i: completes right-handed (approximately +T prograde)
    // - omega: instantaneous frame angular velocity ω = (r×v)/|r|^2
    // -------------------------------------------------------------------------

    /**
     * @brief Construct a spacecraft-centered LVLH frame relative to a primary body.
     *
     * Uses relative position/velocity `r = r_sc - r_primary` and `v = v_sc - v_primary`:
     * - `ex_i`: +R (radial), along `r`
     * - `ez_i`: +N (orbit normal), along `r×v`
     * - `ey_i`: completes right-handed basis (approximately +T, prograde)
     * - `omega_inertial_radps`: instantaneous angular velocity `ω = (r×v)/|r|^2`
     *
     * Returns `std::nullopt` for degenerate relative geometry (near-zero `|r|` or `|r×v|`).
     */
    inline std::optional<RotatingFrame> make_lvlh_frame(const State &primary_state_in, const State &sc_state_in)
    {
        const Vec3 r_rel_m = sc_state_in.position_m - primary_state_in.position_m;
        const Vec3 v_rel_mps = sc_state_in.velocity_mps - primary_state_in.velocity_mps;

        const double r2 = glm::dot(r_rel_m, r_rel_m);
        if (!(r2 > 0.0) || !std::isfinite(r2))
        {
            return std::nullopt;
        }

        const Vec3 ex = normalized_or(r_rel_m, Vec3{0.0, 0.0, 0.0});
        const Vec3 h = glm::cross(r_rel_m, v_rel_mps);
        const double h2 = glm::dot(h, h);
        if (!(h2 > 1e-24) || !std::isfinite(h2))
        {
            return std::nullopt;
        }

        const Vec3 ez = h / std::sqrt(h2);
        Vec3 ey = normalized_or(glm::cross(ez, ex), Vec3{0.0, 0.0, 0.0});
        const Vec3 ex_ortho = normalized_or(glm::cross(ey, ez), ex);
        ey = normalized_or(glm::cross(ez, ex_ortho), ey);

        RotatingFrame f;
        f.origin_position_m = sc_state_in.position_m;
        f.origin_velocity_mps = sc_state_in.velocity_mps;
        f.ex_i = ex_ortho;
        f.ey_i = ey;
        f.ez_i = ez;
        f.omega_inertial_radps = h / r2;
        return f;
    }

} // namespace orbitsim
