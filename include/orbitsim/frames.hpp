#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <cmath>
#include <limits>
#include <optional>

namespace orbitsim
{

    struct SynodicFrame
    {
        Vec3 origin_position_m{0.0, 0.0, 0.0};
        Vec3 origin_velocity_mps{0.0, 0.0, 0.0};

        // Rotating-frame axes expressed in inertial coordinates (orthonormal basis).
        // ex_i points from body A to body B.
        Vec3 ex_i{1.0, 0.0, 0.0};
        Vec3 ey_i{0.0, 1.0, 0.0};
        Vec3 ez_i{0.0, 0.0, 1.0};

        // Frame angular velocity expressed in inertial coordinates [rad/s].
        Vec3 omega_inertial_radps{0.0, 0.0, 0.0};

        // Convenience: |rB - rA| [m] and mu = mB / (mA + mB).
        double separation_m{0.0};
        double mu{0.0};

        inline bool valid() const
        {
            return (separation_m > 0.0) && std::isfinite(separation_m) && std::isfinite(mu) &&
                   std::isfinite(origin_position_m.x) && std::isfinite(origin_position_m.y) &&
                   std::isfinite(origin_position_m.z) && std::isfinite(origin_velocity_mps.x) &&
                   std::isfinite(origin_velocity_mps.y) && std::isfinite(origin_velocity_mps.z) &&
                   std::isfinite(ex_i.x) && std::isfinite(ex_i.y) && std::isfinite(ex_i.z) && std::isfinite(ey_i.x) &&
                   std::isfinite(ey_i.y) && std::isfinite(ey_i.z) && std::isfinite(ez_i.x) && std::isfinite(ez_i.y) &&
                   std::isfinite(ez_i.z) && std::isfinite(omega_inertial_radps.x) &&
                   std::isfinite(omega_inertial_radps.y) && std::isfinite(omega_inertial_radps.z);
        }
    };

    struct Cr3bpLagrangePoints
    {
        // Positions in the synodic rotating frame with origin at the barycenter [m].
        Vec3 primary_m{0.0, 0.0, 0.0};
        Vec3 secondary_m{0.0, 0.0, 0.0};
        Vec3 L1_m{0.0, 0.0, 0.0};
        Vec3 L2_m{0.0, 0.0, 0.0};
        Vec3 L3_m{0.0, 0.0, 0.0};
        Vec3 L4_m{0.0, 0.0, 0.0};
        Vec3 L5_m{0.0, 0.0, 0.0};
    };

    inline Vec3 inertial_vector_to_synodic(const SynodicFrame &frame, const Vec3 &v_in)
    {
        return Vec3{glm::dot(frame.ex_i, v_in), glm::dot(frame.ey_i, v_in), glm::dot(frame.ez_i, v_in)};
    }

    inline Vec3 synodic_vector_to_inertial(const SynodicFrame &frame, const Vec3 &v_rot)
    {
        return frame.ex_i * v_rot.x + frame.ey_i * v_rot.y + frame.ez_i * v_rot.z;
    }

    inline Vec3 inertial_position_to_synodic(const SynodicFrame &frame, const Vec3 &pos_in_m)
    {
        return inertial_vector_to_synodic(frame, pos_in_m - frame.origin_position_m);
    }

    inline Vec3 synodic_position_to_inertial(const SynodicFrame &frame, const Vec3 &pos_rot_m)
    {
        return frame.origin_position_m + synodic_vector_to_inertial(frame, pos_rot_m);
    }

    inline State inertial_state_to_synodic(const State &state_in, const SynodicFrame &frame)
    {
        State out = state_in;

        const Vec3 r_in_rel_m = state_in.position_m - frame.origin_position_m;
        const Vec3 v_in_rel_mps = state_in.velocity_mps - frame.origin_velocity_mps;

        out.position_m = inertial_vector_to_synodic(frame, r_in_rel_m);

        const Vec3 omega_rot_radps = inertial_vector_to_synodic(frame, frame.omega_inertial_radps);
        out.velocity_mps =
                inertial_vector_to_synodic(frame, v_in_rel_mps) - glm::cross(omega_rot_radps, out.position_m);

        return out;
    }

    inline State synodic_state_to_inertial(const State &state_rot, const SynodicFrame &frame)
    {
        State out = state_rot;

        out.position_m = synodic_position_to_inertial(frame, state_rot.position_m);

        const Vec3 omega_rot_radps = inertial_vector_to_synodic(frame, frame.omega_inertial_radps);
        const Vec3 v_in_rel_mps = synodic_vector_to_inertial(
                frame, state_rot.velocity_mps + glm::cross(omega_rot_radps, state_rot.position_m));
        out.velocity_mps = frame.origin_velocity_mps + v_in_rel_mps;

        return out;
    }

    inline std::optional<SynodicFrame> make_synodic_frame(const State &a_state, const double mass_a_kg,
                                                          const State &b_state, const double mass_b_kg)
    {
        const double m_tot = mass_a_kg + mass_b_kg;
        if (!(m_tot > 0.0) || !std::isfinite(m_tot))
        {
            return std::nullopt;
        }

        const Vec3 r_rel_m = b_state.position_m - a_state.position_m;
        const double r2 = glm::dot(r_rel_m, r_rel_m);
        if (!(r2 > 0.0) || !std::isfinite(r2))
        {
            return std::nullopt;
        }
        const double r = std::sqrt(r2);

        const Vec3 v_rel_mps = b_state.velocity_mps - a_state.velocity_mps;

        SynodicFrame out;
        out.separation_m = r;
        out.mu = mass_b_kg / m_tot;

        out.origin_position_m = (mass_a_kg * a_state.position_m + mass_b_kg * b_state.position_m) / m_tot;
        out.origin_velocity_mps = (mass_a_kg * a_state.velocity_mps + mass_b_kg * b_state.velocity_mps) / m_tot;

        out.ex_i = r_rel_m / r;

        Vec3 ez = glm::cross(r_rel_m, v_rel_mps);
        const double ez2 = glm::dot(ez, ez);
        if (!(ez2 > 1e-24) || !std::isfinite(ez2))
        {
            const Vec3 up = (std::abs(out.ex_i.z) < 0.9) ? Vec3{0.0, 0.0, 1.0} : Vec3{0.0, 1.0, 0.0};
            ez = glm::cross(out.ex_i, up);
        }
        out.ez_i = normalized_or(ez, Vec3{0.0, 0.0, 1.0});
        out.ey_i = normalized_or(glm::cross(out.ez_i, out.ex_i), Vec3{0.0, 1.0, 0.0});
        out.ez_i = normalized_or(glm::cross(out.ex_i, out.ey_i), out.ez_i);
        out.ey_i = normalized_or(glm::cross(out.ez_i, out.ex_i), out.ey_i);

        const Vec3 h = glm::cross(r_rel_m, v_rel_mps);
        out.omega_inertial_radps = h / r2;

        return out;
    }

    inline std::optional<SynodicFrame> make_synodic_frame(const MassiveBody &a, const MassiveBody &b)
    {
        return make_synodic_frame(a.state, a.mass_kg, b.state, b.mass_kg);
    }

    inline std::optional<SynodicFrame> make_synodic_frame_at(const CelestialEphemeris &eph, const MassiveBody &a,
                                                             const MassiveBody &b, const double t_s)
    {
        if (!eph.empty())
        {
            const State a_state = eph.body_state_at_by_id(a.id, t_s);
            const State b_state = eph.body_state_at_by_id(b.id, t_s);
            return make_synodic_frame(a_state, a.mass_kg, b_state, b.mass_kg);
        }
        return make_synodic_frame(a, b);
    }

    namespace detail
    {
        inline double cr3bp_collinear_equation_nd_(const double x, const double mu)
        {
            const double x1 = x + mu;
            const double x2 = x - (1.0 - mu);

            const double ax1 = std::abs(x1);
            const double ax2 = std::abs(x2);
            if (!(ax1 > 0.0) || !(ax2 > 0.0) || !std::isfinite(ax1) || !std::isfinite(ax2))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }

            const double inv1 = 1.0 / (ax1 * ax1 * ax1);
            const double inv2 = 1.0 / (ax2 * ax2 * ax2);

            return x - (1.0 - mu) * (x1 * inv1) - mu * (x2 * inv2);
        }

        inline std::optional<double> solve_cr3bp_collinear_root_nd_(const double mu, double x_min, double x_max)
        {
            if (!(x_min < x_max) || !std::isfinite(x_min) || !std::isfinite(x_max) || !std::isfinite(mu))
            {
                return std::nullopt;
            }

            double a = x_min;
            double b = x_max;
            double fa = cr3bp_collinear_equation_nd_(a, mu);
            double fb = cr3bp_collinear_equation_nd_(b, mu);
            if (!std::isfinite(fa) || !std::isfinite(fb))
            {
                return std::nullopt;
            }
            if (fa == 0.0)
            {
                return a;
            }
            if (fb == 0.0)
            {
                return b;
            }
            if (fa * fb > 0.0)
            {
                return std::nullopt;
            }

            for (int iter = 0; iter < 128; ++iter)
            {
                const double m = 0.5 * (a + b);
                const double fm = cr3bp_collinear_equation_nd_(m, mu);
                if (!std::isfinite(fm))
                {
                    return std::nullopt;
                }
                if (fm == 0.0)
                {
                    return m;
                }

                if (fa * fm < 0.0)
                {
                    b = m;
                    fb = fm;
                }
                else
                {
                    a = m;
                    fa = fm;
                }

                if (std::abs(b - a) <= 1e-14)
                {
                    break;
                }
            }

            return 0.5 * (a + b);
        }
    } // namespace detail

    inline std::optional<Cr3bpLagrangePoints> cr3bp_lagrange_points_m(const SynodicFrame &frame)
    {
        if (!frame.valid())
        {
            return std::nullopt;
        }
        if (!(frame.mu > 0.0) || !(frame.mu < 1.0))
        {
            return std::nullopt;
        }

        const double a_m = frame.separation_m;
        if (!(a_m > 0.0) || !std::isfinite(a_m))
        {
            return std::nullopt;
        }

        const double mu = frame.mu;
        const double x_primary_m = (-mu) * a_m;
        const double x_secondary_m = (1.0 - mu) * a_m;

        const double delta = 1e-12;
        const double x1 = -mu;
        const double x2 = 1.0 - mu;

        const std::optional<double> xL1 = detail::solve_cr3bp_collinear_root_nd_(mu, x1 + delta, x2 - delta);
        const std::optional<double> xL2 = detail::solve_cr3bp_collinear_root_nd_(mu, x2 + delta, x2 + 5.0);
        const std::optional<double> xL3 = detail::solve_cr3bp_collinear_root_nd_(mu, x1 - 5.0, x1 - delta);

        if (!xL1.has_value() || !xL2.has_value() || !xL3.has_value())
        {
            return std::nullopt;
        }

        Cr3bpLagrangePoints out;
        out.primary_m = Vec3{x_primary_m, 0.0, 0.0};
        out.secondary_m = Vec3{x_secondary_m, 0.0, 0.0};

        out.L1_m = Vec3{(*xL1) * a_m, 0.0, 0.0};
        out.L2_m = Vec3{(*xL2) * a_m, 0.0, 0.0};
        out.L3_m = Vec3{(*xL3) * a_m, 0.0, 0.0};

        const double s3 = std::sqrt(3.0);
        out.L4_m = Vec3{(0.5 - mu) * a_m, 0.5 * s3 * a_m, 0.0};
        out.L5_m = Vec3{(0.5 - mu) * a_m, -0.5 * s3 * a_m, 0.0};

        return out;
    }

    inline std::optional<Cr3bpLagrangePoints> cr3bp_lagrange_points_m_at(const CelestialEphemeris &eph,
                                                                         const MassiveBody &a, const MassiveBody &b,
                                                                         const double t_s)
    {
        const std::optional<SynodicFrame> frame = make_synodic_frame_at(eph, a, b, t_s);
        if (!frame.has_value())
        {
            return std::nullopt;
        }
        return cr3bp_lagrange_points_m(*frame);
    }

} // namespace orbitsim
