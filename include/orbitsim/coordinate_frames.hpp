#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <cmath>
#include <optional>

namespace orbitsim
{

    // -------------------------------------------------------------------------
    // Generic rotating/translating frame relative to the library's inertial frame.
    //
    // Convention:
    // - Basis vectors ex_i/ey_i/ez_i are expressed in inertial coordinates (orthonormal).
    // - omega_inertial_radps is the frame angular velocity expressed in inertial coords.
    // - Transforming inertial -> frame applies the non-inertial velocity term: v' = R^T v - ω'×r'
    // -------------------------------------------------------------------------

    struct RotatingFrame
    {
        Vec3 origin_position_m{0.0, 0.0, 0.0};
        Vec3 origin_velocity_mps{0.0, 0.0, 0.0};

        Vec3 ex_i{1.0, 0.0, 0.0};
        Vec3 ey_i{0.0, 1.0, 0.0};
        Vec3 ez_i{0.0, 0.0, 1.0};

        Vec3 omega_inertial_radps{0.0, 0.0, 0.0};

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

    inline Vec3 rotate_about_axis(const Vec3 &v, const Vec3 &axis_unit, const double angle_rad)
    {
        const double c = std::cos(angle_rad);
        const double s = std::sin(angle_rad);
        return v * c + glm::cross(axis_unit, v) * s + axis_unit * (glm::dot(axis_unit, v) * (1.0 - c));
    }

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
