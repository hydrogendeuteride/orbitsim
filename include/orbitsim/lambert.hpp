#pragma once

#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

namespace orbitsim
{

    /**
     * @brief Options for the Lambert solver.
     */
    struct LambertOptions
    {
        /// If true, select prograde transfer (orbit normal aligned with +Z).
        /// Uses sign of (r1 x r2).z to determine direction.
        bool prograde{true};

        /// If true, use short-path transfer (< 180 deg); otherwise long-path.
        bool short_path{true};

        /// Maximum revolutions to consider (0 = direct transfer only).
        /// Multi-rev solutions exist when TOF exceeds minimum-energy transfer time.
        int max_revolutions{0};

        /// Maximum iterations for bisection root-finder.
        int max_bisect_iters{96};

        /// Convergence tolerance for time-of-flight matching [s].
        double time_abs_tolerance_s{1e-9};

        /// Number of samples for root scanning (0 = auto, clamped to [2000, 40000]).
        /// Higher values find multi-rev roots more reliably at cost of performance.
        int scan_samples{0};
    };

    /**
     * @brief Result of a Lambert problem solution.
     */
    struct LambertSolution
    {
        Vec3 v1_mps{0.0, 0.0, 0.0};  ///< Velocity at departure position [m/s].
        Vec3 v2_mps{0.0, 0.0, 0.0};  ///< Velocity at arrival position [m/s].
        double z{0.0};               ///< Universal variable root parameter.
        bool converged{false};       ///< True if solution converged within tolerance.
        int iterations{0};           ///< Number of iterations used.
    };

    namespace detail
    {
        inline double clamp_cos_(const double c) { return std::clamp(c, -1.0, 1.0); }

        /**
         * @brief Computes the transfer angle between two position vectors.
         * @param r1 Departure position vector [m].
         * @param r2 Arrival position vector [m].
         * @param opt Lambert options (prograde/retrograde, short/long path).
         * @return Transfer angle in radians, or nullopt if degenerate (collinear vectors).
         */
        inline std::optional<double> transfer_angle_(const Vec3 &r1, const Vec3 &r2, const LambertOptions &opt)
        {
            const double r1n = glm::length(r1);
            const double r2n = glm::length(r2);
            if (!(r1n > 0.0) || !(r2n > 0.0) || !std::isfinite(r1n) || !std::isfinite(r2n))
            {
                return std::nullopt;
            }

            const Vec3 c = glm::cross(r1, r2);
            const double c2 = glm::dot(c, c);
            if (!(c2 > 0.0) || !std::isfinite(c2))
            {
                const double cosT = clamp_cos_(glm::dot(r1, r2) / (r1n * r2n));
                const double theta = std::acos(cosT);
                return theta;
            }

            const double cosT = clamp_cos_(glm::dot(r1, r2) / (r1n * r2n));
            const double theta = std::acos(cosT); // [0, pi]

            const bool alpha_pos = (c.z > 0.0);
            double dTheta = 0.0;
            if (opt.prograde)
            {
                dTheta = alpha_pos ? theta : (2.0 * std::acos(-1.0) - theta);
            }
            else
            {
                dTheta = (!alpha_pos) ? theta : (2.0 * std::acos(-1.0) - theta);
            }

            if (opt.short_path)
            {
                if (dTheta > std::acos(-1.0))
                {
                    dTheta = 2.0 * std::acos(-1.0) - dTheta;
                }
            }
            else
            {
                if (dTheta < std::acos(-1.0))
                {
                    dTheta = 2.0 * std::acos(-1.0) - dTheta;
                }
            }

            return dTheta;
        }

        /// Intermediate result from time-of-flight evaluation.
        struct TimeOfFlightEval
        {
            double dt_s{0.0};  ///< Computed time-of-flight [s].
            double y{0.0};     ///< Intermediate variable for velocity computation.
            double C{0.0};     ///< Stumpff C(z) value.
            double S{0.0};     ///< Stumpff S(z) value.
            bool valid{false}; ///< True if evaluation succeeded.
        };

        /**
         * @brief Evaluates time-of-flight for a given universal variable z.
         *
         * Uses Stumpff functions C(z) and S(z) to compute TOF for the universal
         * variable formulation of Lambert's problem.
         *
         * @param mu_m3_s2 Gravitational parameter [m^3/s^2].
         * @param r1_m Magnitude of departure position [m].
         * @param r2_m Magnitude of arrival position [m].
         * @param A Geometry parameter: sin(dTheta) * sqrt(r1*r2 / (1 - cos(dTheta))).
         * @param z Universal variable to evaluate.
         * @return TimeOfFlightEval with computed TOF and intermediate values.
         */
        inline TimeOfFlightEval time_of_flight_universal_(const double mu_m3_s2, const double r1_m, const double r2_m,
                                                          const double A, const double z)
        {
            TimeOfFlightEval out;
            if (!(mu_m3_s2 > 0.0) || !std::isfinite(mu_m3_s2))
            {
                return out;
            }

            const double C = stumpff_C(z);
            const double S = stumpff_S(z);
            if (!(C > 0.0) || !std::isfinite(C) || !std::isfinite(S))
            {
                return out;
            }

            const double sqrtC = std::sqrt(C);
            if (!(sqrtC > 0.0) || !std::isfinite(sqrtC))
            {
                return out;
            }

            const double y = r1_m + r2_m + (A * (z * S - 1.0)) / sqrtC;
            if (!(y > 0.0) || !std::isfinite(y))
            {
                return out;
            }

            const double x = std::sqrt(y / C);
            if (!(x > 0.0) || !std::isfinite(x))
            {
                return out;
            }

            const double sqrt_mu = std::sqrt(mu_m3_s2);
            const double dt_s = (x * x * x * S + A * std::sqrt(y)) / sqrt_mu;
            if (!std::isfinite(dt_s))
            {
                return out;
            }

            out.dt_s = dt_s;
            out.y = y;
            out.C = C;
            out.S = S;
            out.valid = true;
            return out;
        }

        /**
         * @brief Finds z root via bisection where TOF(z) = target_dt_s.
         * @param z_lo Lower bound of search interval.
         * @param z_hi Upper bound of search interval.
         * @param target_dt_s Desired time-of-flight [s].
         * @param out_iters If non-null, receives iteration count.
         * @return The z value satisfying TOF(z) ≈ target_dt_s, or nullopt if no root.
         */
        inline std::optional<double> bisect_root_(const double mu_m3_s2, const double r1_m, const double r2_m, const double A,
                                                  const double target_dt_s, const double z_lo, const double z_hi,
                                                  const LambertOptions &opt, int *out_iters)
        {
            if (out_iters != nullptr)
            {
                *out_iters = 0;
            }

            TimeOfFlightEval e_lo = time_of_flight_universal_(mu_m3_s2, r1_m, r2_m, A, z_lo);
            TimeOfFlightEval e_hi = time_of_flight_universal_(mu_m3_s2, r1_m, r2_m, A, z_hi);
            if (!e_lo.valid || !e_hi.valid)
            {
                return std::nullopt;
            }

            double f_lo = e_lo.dt_s - target_dt_s;
            double f_hi = e_hi.dt_s - target_dt_s;
            if (!std::isfinite(f_lo) || !std::isfinite(f_hi))
            {
                return std::nullopt;
            }
            if (f_lo == 0.0)
            {
                return z_lo;
            }
            if (f_hi == 0.0)
            {
                return z_hi;
            }
            if (f_lo * f_hi > 0.0)
            {
                return std::nullopt;
            }

            double a = z_lo;
            double b = z_hi;
            double fa = f_lo;
            double fb = f_hi;

            const int max_iters = std::max(1, opt.max_bisect_iters);
            for (int it = 0; it < max_iters; ++it)
            {
                if (out_iters != nullptr)
                {
                    *out_iters = it + 1;
                }

                const double m = 0.5 * (a + b);
                TimeOfFlightEval e_m = time_of_flight_universal_(mu_m3_s2, r1_m, r2_m, A, m);
                if (!e_m.valid)
                {
                    break;
                }

                const double fm = e_m.dt_s - target_dt_s;
                if (!std::isfinite(fm))
                {
                    break;
                }

                if (std::abs(fm) <= std::max(0.0, opt.time_abs_tolerance_s))
                {
                    return m;
                }

                const bool left = (fa <= 0.0 && fm >= 0.0) || (fa >= 0.0 && fm <= 0.0);
                if (left)
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
                    return 0.5 * (a + b);
                }
            }

            return 0.5 * (a + b);
        }

        /// Adds z to roots if no existing root is within eps of it.
        inline bool add_unique_root_(std::vector<double> &roots, const double z, const double eps)
        {
            for (double r: roots)
            {
                if (std::abs(r - z) <= eps)
                {
                    return false;
                }
            }
            roots.push_back(z);
            return true;
        }

    } // namespace detail

    /**
     * @brief Solves Lambert's problem using universal variable formulation.
     *
     * Given two position vectors and a time-of-flight, finds the velocities
     * that connect them via a conic trajectory under two-body dynamics.
     *
     * @param mu_m3_s2 Gravitational parameter of central body [m^3/s^2].
     * @param r1_m     Departure position vector [m].
     * @param r2_m     Arrival position vector [m].
     * @param dt_s     Time-of-flight [s], must be positive.
     * @param opt      Solver options (transfer type, multi-rev settings).
     * @return Vector of solutions (may be empty, or contain multiple for multi-rev).
     *
     * @note Returns empty vector for degenerate cases (collinear r1/r2, dt <= 0).
     * @note Multi-rev transfers can yield up to 2 solutions per revolution count.
     */
    inline std::vector<LambertSolution>
    solve_lambert_universal(const double mu_m3_s2, const Vec3 &r1_m, const Vec3 &r2_m, const double dt_s,
                            const LambertOptions &opt = {})
    {
        std::vector<LambertSolution> out;
        if (!(mu_m3_s2 > 0.0) || !std::isfinite(mu_m3_s2) || !(dt_s > 0.0) || !std::isfinite(dt_s))
        {
            return out;
        }

        const double r1n = glm::length(r1_m);
        const double r2n = glm::length(r2_m);
        if (!(r1n > 0.0) || !(r2n > 0.0) || !std::isfinite(r1n) || !std::isfinite(r2n))
        {
            return out;
        }

        const std::optional<double> dTheta_opt = detail::transfer_angle_(r1_m, r2_m, opt);
        if (!dTheta_opt.has_value())
        {
            return out;
        }
        const double dTheta = *dTheta_opt;

        const double cos_dT = std::cos(dTheta);
        const double sin_dT = std::sin(dTheta);
        const double denom = 1.0 - cos_dT;
        if (!(std::abs(sin_dT) > 0.0) || !(denom > 0.0) || !std::isfinite(denom))
        {
            return out;
        }

        const double A = sin_dT * std::sqrt((r1n * r2n) / denom);
        if (!(std::abs(A) > 0.0) || !std::isfinite(A))
        {
            return out;
        }

        const int max_revs = std::max(0, opt.max_revolutions);
        const double pi = std::acos(-1.0);
        const double z_min = -4.0 * pi * pi;
        const double z_max = std::pow(2.0 * pi * static_cast<double>(max_revs + 1), 2.0);

        int scan_samples = opt.scan_samples;
        if (scan_samples <= 0)
        {
            scan_samples = 4000 * (max_revs + 1);
        }
        scan_samples = std::clamp(scan_samples, 2000, 40000);

        const double dz = (z_max - z_min) / static_cast<double>(scan_samples);
        if (!(dz > 0.0) || !std::isfinite(dz))
        {
            return out;
        }

        std::vector<double> roots;
        roots.reserve(static_cast<std::size_t>(2 * (max_revs + 1) + 2));

        auto eval_f = [&](const double z) -> std::optional<double> {
            const detail::TimeOfFlightEval e = detail::time_of_flight_universal_(mu_m3_s2, r1n, r2n, A, z);
            if (!e.valid)
            {
                return std::nullopt;
            }
            return e.dt_s - dt_s;
        };

        double z_prev = z_min;
        std::optional<double> f_prev = eval_f(z_prev);

        for (int i = 1; i <= scan_samples; ++i)
        {
            const double z = z_min + dz * static_cast<double>(i);
            const std::optional<double> f = eval_f(z);

            if (!f_prev.has_value() || !f.has_value())
            {
                z_prev = z;
                f_prev = f;
                continue;
            }

            const double fa = *f_prev;
            const double fb = *f;
            if (!std::isfinite(fa) || !std::isfinite(fb))
            {
                z_prev = z;
                f_prev = f;
                continue;
            }

            if (fa == 0.0)
            {
                detail::add_unique_root_(roots, z_prev, 1e-10);
            }
            if (fa * fb < 0.0)
            {
                int iters = 0;
                const std::optional<double> z_root =
                        detail::bisect_root_(mu_m3_s2, r1n, r2n, A, dt_s, z_prev, z, opt, &iters);
                if (z_root.has_value())
                {
                    detail::add_unique_root_(roots, *z_root, 1e-10);
                }
            }

            z_prev = z;
            f_prev = f;
        }

        std::sort(roots.begin(), roots.end());

        out.reserve(roots.size());
        for (double z: roots)
        {
            const detail::TimeOfFlightEval e = detail::time_of_flight_universal_(mu_m3_s2, r1n, r2n, A, z);
            if (!e.valid)
            {
                continue;
            }

            const double y = e.y;
            const double g = A * std::sqrt(y / mu_m3_s2);
            if (!(std::abs(g) > 0.0) || !std::isfinite(g))
            {
                continue;
            }

            const double f = 1.0 - y / r1n;
            const double gdot = 1.0 - y / r2n;

            const Vec3 v1 = (r2_m - f * r1_m) / g;
            const Vec3 v2 = (gdot * r2_m - r1_m) / g;

            LambertSolution sol;
            sol.v1_mps = v1;
            sol.v2_mps = v2;
            sol.z = z;
            sol.converged = true;
            sol.iterations = 0;
            out.push_back(sol);
        }

        return out;
    }

} // namespace orbitsim

