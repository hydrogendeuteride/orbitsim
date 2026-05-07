#pragma once

#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace orbitsim
{
    enum class KeplerStatus
    {
        Ok,
        InvalidMu,
        InvalidTime,
        InvalidInitialState,
        StumpffOverflow,
        ZeroDerivative,
        DidNotConverge,
        InvalidFinalState,
    };

    enum class KeplerOrbitRegime
    {
        Unknown,
        Elliptic,
        NearParabolic,
        Hyperbolic,
    };

    struct KeplerPropagationOptions
    {
        int max_iterations{64};
        double abs_tolerance{1e-12};
        double rel_tolerance{1e-12};
        double step_tolerance{1e-12};
        double rel_step_tolerance{1e-11};
        double near_parabolic_alpha_tolerance{1e-12};
        double max_hyperbolic_stumpff_arg{700.0};
    };

    struct KeplerDiagnostics
    {
        KeplerStatus status{KeplerStatus::InvalidInitialState};
        KeplerOrbitRegime regime{KeplerOrbitRegime::Unknown};
        int iterations{0};
        double residual{std::numeric_limits<double>::quiet_NaN()};
        double chi{0.0};
        double alpha{std::numeric_limits<double>::quiet_NaN()};
        double z{std::numeric_limits<double>::quiet_NaN()};
        bool used_fallback{false};
    };

    struct KeplerPropagationResult
    {
        State state{};
        KeplerDiagnostics diagnostics{};

        [[nodiscard]] bool ok() const { return diagnostics.status == KeplerStatus::Ok; }
    };

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

    namespace detail
    {
        inline bool kepler_finite3_(const Vec3 &v)
        {
            return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
        }

        inline bool safe_stumpff_(const double z, const KeplerPropagationOptions &opt, double &out_C, double &out_S)
        {
            if (!std::isfinite(z))
            {
                return false;
            }

            if (z < 0.0)
            {
                const double arg = std::sqrt(-z);
                if (!std::isfinite(arg) || arg > opt.max_hyperbolic_stumpff_arg)
                {
                    return false;
                }
            }

            out_C = stumpff_C(z);
            out_S = stumpff_S(z);
            return std::isfinite(out_C) && std::isfinite(out_S);
        }
    } // namespace detail

    inline KeplerPropagationResult propagate_kepler_universal_safe(const double mu_m3_s2, const Vec3 &r0_m,
                                                                   const Vec3 &v0_mps, const double dt_s,
                                                                   const KeplerPropagationOptions &opt = {})
    {
        KeplerPropagationResult out;
        out.state = make_state(r0_m, v0_mps);

        if (!(mu_m3_s2 > 0.0) || !std::isfinite(mu_m3_s2))
        {
            out.diagnostics.status = KeplerStatus::InvalidMu;
            return out;
        }
        if (!std::isfinite(dt_s))
        {
            out.diagnostics.status = KeplerStatus::InvalidTime;
            return out;
        }
        if (!detail::kepler_finite3_(r0_m) || !detail::kepler_finite3_(v0_mps))
        {
            out.diagnostics.status = KeplerStatus::InvalidInitialState;
            return out;
        }

        const double r0 = glm::length(r0_m);
        const double v0 = glm::length(v0_mps);
        if (!(r0 > 0.0) || !std::isfinite(r0) || !std::isfinite(v0))
        {
            out.diagnostics.status = KeplerStatus::InvalidInitialState;
            return out;
        }

        if (dt_s == 0.0)
        {
            out.diagnostics.status = KeplerStatus::Ok;
            out.diagnostics.residual = 0.0;
            return out;
        }

        const double sqrt_mu = std::sqrt(mu_m3_s2);
        const double vr0 = glm::dot(r0_m, v0_mps) / r0;
        const double alpha = (2.0 / r0) - ((v0 * v0) / mu_m3_s2);
        out.diagnostics.alpha = alpha;

        const double alpha_tol = std::max(0.0, opt.near_parabolic_alpha_tolerance);
        if (alpha > alpha_tol)
        {
            out.diagnostics.regime = KeplerOrbitRegime::Elliptic;
        }
        else if (alpha < -alpha_tol)
        {
            out.diagnostics.regime = KeplerOrbitRegime::Hyperbolic;
        }
        else
        {
            out.diagnostics.regime = KeplerOrbitRegime::NearParabolic;
        }

        double chi = 0.0;
        if (std::abs(alpha) > alpha_tol)
        {
            chi = sqrt_mu * std::abs(alpha) * dt_s;
        }
        else
        {
            chi = sqrt_mu * dt_s / r0;
        }

        const int max_iter = std::max(1, opt.max_iterations);
        const double residual_scale = std::max(1.0, std::abs(sqrt_mu * dt_s));
        const double residual_tolerance =
                std::max(std::abs(opt.abs_tolerance), std::abs(opt.rel_tolerance) * residual_scale);
        const double abs_step_tolerance = std::abs(opt.step_tolerance);
        const double rel_step_tolerance = std::abs(opt.rel_step_tolerance);

        bool converged = false;
        for (int it = 0; it < max_iter; ++it)
        {
            ++out.diagnostics.iterations;
            const double z = alpha * chi * chi;
            out.diagnostics.z = z;

            double C = 0.0;
            double S = 0.0;
            if (!detail::safe_stumpff_(z, opt, C, S))
            {
                out.diagnostics.status = KeplerStatus::StumpffOverflow;
                out.diagnostics.chi = chi;
                return out;
            }

            const double F = (r0 * vr0 / sqrt_mu) * (chi * chi) * C +
                             (1.0 - alpha * r0) * (chi * chi * chi) * S + r0 * chi - sqrt_mu * dt_s;
            out.diagnostics.residual = std::abs(F);

            if (out.diagnostics.residual <= residual_tolerance)
            {
                converged = true;
                break;
            }

            const double dF = (r0 * vr0 / sqrt_mu) * chi * (1.0 - z * S) +
                              (1.0 - alpha * r0) * (chi * chi) * C + r0;
            if (!(std::abs(dF) > 0.0) || !std::isfinite(dF))
            {
                out.diagnostics.status = KeplerStatus::ZeroDerivative;
                out.diagnostics.chi = chi;
                return out;
            }

            const double delta = F / dF;
            chi -= delta;
            if (!std::isfinite(chi))
            {
                out.diagnostics.status = KeplerStatus::InvalidFinalState;
                out.diagnostics.chi = chi;
                return out;
            }
            const double relative_step_limit = rel_step_tolerance * std::max(1.0, std::abs(chi));
            if (std::abs(delta) <= std::max(abs_step_tolerance, relative_step_limit))
            {
                converged = true;
                break;
            }
        }

        out.diagnostics.chi = chi;
        if (!converged)
        {
            out.diagnostics.status = KeplerStatus::DidNotConverge;
        }

        const double z = alpha * chi * chi;
        out.diagnostics.z = z;
        double C = 0.0;
        double S = 0.0;
        if (!detail::safe_stumpff_(z, opt, C, S))
        {
            out.diagnostics.status = KeplerStatus::StumpffOverflow;
            return out;
        }

        const double f = 1.0 - (chi * chi / r0) * C;
        const double g = dt_s - (1.0 / sqrt_mu) * (chi * chi * chi) * S;

        const Vec3 r_m = (f * r0_m) + (g * v0_mps);
        const double r = glm::length(r_m);
        if (!(r > 0.0) || !std::isfinite(r))
        {
            out.diagnostics.status = KeplerStatus::InvalidFinalState;
            return out;
        }

        const double fdot = (sqrt_mu / (r * r0)) * (z * S - 1.0) * chi;
        const double gdot = 1.0 - (chi * chi / r) * C;
        const Vec3 v_mps = (fdot * r0_m) + (gdot * v0_mps);
        if (!detail::kepler_finite3_(r_m) || !detail::kepler_finite3_(v_mps))
        {
            out.diagnostics.status = KeplerStatus::InvalidFinalState;
            return out;
        }

        out.state.position_m = r_m;
        out.state.velocity_mps = v_mps;
        if (converged)
        {
            out.diagnostics.status = KeplerStatus::Ok;
        }
        return out;
    }

    inline KeplerPropagationResult propagate_kepler_universal_safe(const double mu_m3_s2, const State &in,
                                                                   const double dt_s,
                                                                   const KeplerPropagationOptions &opt = {})
    {
        KeplerPropagationResult out =
                propagate_kepler_universal_safe(mu_m3_s2, in.position_m, in.velocity_mps, dt_s, opt);
        out.state.spin = in.spin;
        return out;
    }

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
