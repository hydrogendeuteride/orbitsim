#pragma once

#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>

namespace orbitsim
{
    /** @brief Options for Kepler propagator Newton iteration. */
    struct KeplerOptions
    {
        int max_iterations{64};       ///< Maximum Newton iterations
        double abs_tolerance{1e-12};  ///< Convergence tolerance
    };

    /** @brief Result of Kepler propagation. */
    struct KeplerStepResult
    {
        Vec3 position_m{0.0, 0.0, 0.0};
        Vec3 velocity_mps{0.0, 0.0, 0.0};
        bool converged{false};  ///< True if Newton iteration converged
        int iterations{0};      ///< Number of iterations used
    };

    /**
     * @brief Propagate two-body orbit using universal variable formulation.
     *
     * Analytically propagates position and velocity under point-mass gravity.
     * Works for all orbit types (elliptic, parabolic, hyperbolic) without
     * special-casing. Uses Newton iteration to solve Kepler's equation.
     *
     * @param mu_m3_s2 Gravitational parameter (G*M) of the central body
     * @param r0_m Initial position relative to central body
     * @param v0_mps Initial velocity relative to central body
     * @param dt_s Time to propagate (can be negative for backward propagation)
     * @param opt Solver options
     * @return Propagated state with convergence info
     */
    inline KeplerStepResult propagate_kepler_universal(const double mu_m3_s2, const Vec3 &r0_m, const Vec3 &v0_mps,
                                                       const double dt_s, const KeplerOptions &opt = {})
    {
        KeplerStepResult out;
        out.position_m = r0_m;
        out.velocity_mps = v0_mps;

        if (!(mu_m3_s2 > 0.0) || !std::isfinite(mu_m3_s2) || !std::isfinite(dt_s))
        {
            return out;
        }

        const double r0 = glm::length(r0_m);
        const double v0 = glm::length(v0_mps);
        if (!(r0 > 0.0) || !std::isfinite(r0) || !std::isfinite(v0))
        {
            return out;
        }

        const double sqrt_mu = std::sqrt(mu_m3_s2);
        const double vr0 = glm::dot(r0_m, v0_mps) / r0;
        const double alpha = (2.0 / r0) - ((v0 * v0) / mu_m3_s2);

        double chi = 0.0;
        if (std::abs(alpha) > 1e-12)
        {
            chi = sqrt_mu * std::abs(alpha) * dt_s;
        }
        else
        {
            chi = sqrt_mu * dt_s / r0;
        }

        // Newton iterations for universal anomaly chi
        bool converged = false;
        const int max_iter = std::max(1, opt.max_iterations);
        int iterations = 0;
        for (int it = 0; it < max_iter; ++it)
        {
            ++iterations;
            const double z = alpha * chi * chi;
            const double C = stumpff_C(z);
            const double S = stumpff_S(z);

            const double F = (r0 * vr0 / sqrt_mu) * (chi * chi) * C + (1.0 - alpha * r0) * (chi * chi * chi) * S +
                             r0 * chi - sqrt_mu * dt_s;

            if (std::abs(F) <= opt.abs_tolerance)
            {
                converged = true;
                break;
            }

            const double dF = (r0 * vr0 / sqrt_mu) * chi * (1.0 - z * S) + (1.0 - alpha * r0) * (chi * chi) * C + r0;

            if (!(std::abs(dF) > 0.0) || !std::isfinite(dF))
            {
                break;
            }

            const double delta = F / dF;
            chi -= delta;
            if (!std::isfinite(chi))
            {
                break;
            }
            if (std::abs(delta) <= opt.abs_tolerance)
            {
                converged = true;
                break;
            }
        }

        out.converged = converged;
        out.iterations = iterations;

        const double z = alpha * chi * chi;
        const double C = stumpff_C(z);
        const double S = stumpff_S(z);

        const double f = 1.0 - (chi * chi / r0) * C;
        const double g = dt_s - (1.0 / sqrt_mu) * (chi * chi * chi) * S;

        const Vec3 r_m = (f * r0_m) + (g * v0_mps);
        const double r = glm::length(r_m);
        if (!(r > 0.0) || !std::isfinite(r))
        {
            return out;
        }

        const double fdot = (sqrt_mu / (r * r0)) * (z * S - 1.0) * chi;
        const double gdot = 1.0 - (chi * chi / r) * C;

        out.position_m = r_m;
        out.velocity_mps = (fdot * r0_m) + (gdot * v0_mps);
        return out;
    }

    /** @brief Overload accepting State instead of separate position/velocity. */
    inline KeplerStepResult propagate_kepler_universal(const double mu_m3_s2, const State &in, const double dt_s,
                                                       const KeplerOptions &opt = {})
    {
        return propagate_kepler_universal(mu_m3_s2, in.position_m, in.velocity_mps, dt_s, opt);
    }
} // namespace orbitsim
