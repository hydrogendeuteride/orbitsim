#pragma once

#include "orbitsim/coordinate_frames.hpp"
#include "orbitsim/ephemeris.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <cmath>
#include <limits>
#include <optional>

namespace orbitsim
{

    /**
     * @brief Two-body synodic (co-rotating) frame derived from two inertial body states (A,B).
     *
     * This is a specialized RotatingFrame primarily for display/analysis. It intentionally does not
     * implement CR3BP dynamics; it only provides a consistent translating/rotating frame definition.
     */
    struct SynodicFrame : RotatingFrame
    {
        /// Separation distance `|r_B - r_A|` [m].
        double separation_m{0.0};
        /// Mass ratio `μ = m_B / (m_A + m_B)` (dimensionless).
        double mu{0.0};

        /// @brief Returns true if the frame and derived parameters are finite and well-formed.
        inline bool valid() const
        {
            return RotatingFrame::valid() && (separation_m > 0.0) && std::isfinite(separation_m) && std::isfinite(mu);
        }
    };

    /**
     * @brief CR3BP Lagrange points in the synodic rotating frame (origin at the barycenter).
     *
     * All positions are expressed in the synodic frame coordinates and are dimensional [m].
     */
    struct Cr3bpLagrangePoints
    {
        /// Primary body position [m].
        Vec3 primary_m{0.0, 0.0, 0.0};
        /// Secondary body position [m].
        Vec3 secondary_m{0.0, 0.0, 0.0};
        /// Collinear point between the primaries [m].
        Vec3 L1_m{0.0, 0.0, 0.0};
        /// Collinear point beyond the secondary [m].
        Vec3 L2_m{0.0, 0.0, 0.0};
        /// Collinear point beyond the primary [m].
        Vec3 L3_m{0.0, 0.0, 0.0};
        /// Triangular point leading the secondary by 60° [m].
        Vec3 L4_m{0.0, 0.0, 0.0};
        /// Triangular point trailing the secondary by 60° [m].
        Vec3 L5_m{0.0, 0.0, 0.0};
    };

    /**
     * @brief Construct a synodic (barycentric, co-rotating) frame from two inertial body states.
     *
     * The frame is defined as:
     * - Origin: barycenter position/velocity of (A,B)
     * - `ex_i`: along the relative position from A to B
     * - `ez_i`: along `r×v` (orbit normal), falling back to a deterministic perpendicular if degenerate
     * - `ey_i`: completes a right-handed orthonormal basis
     * - `omega_inertial_radps`: instantaneous angular velocity `ω = (r×v)/|r|^2` (in inertial coords)
     *
     * @return `std::nullopt` if masses are invalid/non-finite, bodies are coincident, or geometry is degenerate.
     */
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
        /// @brief CR3BP collinear equilibrium equation (non-dimensional, x-axis with y=z=0).
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

        /**
         * @brief Solve for a collinear CR3BP root in a bracket using bisection (non-dimensional x).
         *
         * @param mu Mass ratio `μ` (dimensionless).
         * @param x_min Lower bracket bound.
         * @param x_max Upper bracket bound.
         * @return Root if the bracket is valid and a sign change exists; `std::nullopt` otherwise.
         */
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

    /**
     * @brief Compute CR3BP L1-L5 positions in meters from a synodic frame.
     *
     * Uses the synodic frame's `separation_m` and `mu`:
     * - L4/L5 are placed at the classical equilateral locations
     * - L1/L2/L3 are solved along the x-axis using bisection on the collinear equilibrium equation
     *
     * @return `std::nullopt` if `frame` is invalid, `mu` is not in (0,1), or root finding fails.
     */
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

    /// @brief Compute CR3BP L1-L5 at time `t_s` (builds a synodic frame via ephemeris if available).
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
