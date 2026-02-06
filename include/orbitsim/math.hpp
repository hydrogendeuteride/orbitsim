#pragma once

#include "orbitsim/frame_utils.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace orbitsim
{

    inline constexpr double kGravitationalConstant_SI = 6.67430e-11; // m^3 / (kg s^2)
    inline constexpr double kStandardGravity_mps2 = 9.80665; // m / s^2

    inline Vec3 normalized_or(const Vec3 &v, const Vec3 &fallback_unit)
    {
        const double len2 = glm::dot(v, v);
        if (!(len2 > 0.0) || !std::isfinite(len2))
        {
            return fallback_unit;
        }
        return v / std::sqrt(len2);
    }

    inline double clamp01(const double x) { return std::clamp(x, 0.0, 1.0); }

    /**
     * @brief Cubic Hermite interpolation for position between two endpoints.
     * @param p0 Start position
     * @param v0_mps Velocity at start
     * @param p1 End position
     * @param v1_mps Velocity at end
     * @param dt_s Time span of the interval (p0 to p1)
     * @param tau01 Normalized parameter in [0,1]; clamped internally
     * @return Interpolated position
     */
    inline Vec3 hermite_position(const Vec3 &p0, const Vec3 &v0_mps, const Vec3 &p1, const Vec3 &v1_mps,
                                 const double dt_s, const double tau01)
    {
        const double t = clamp01(tau01);
        const double t2 = t * t;
        const double t3 = t2 * t;
        const double h = dt_s;

        const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        const double h10 = t3 - 2.0 * t2 + t;
        const double h01 = -2.0 * t3 + 3.0 * t2;
        const double h11 = t3 - t2;

        return (h00 * p0) + (h10 * (h * v0_mps)) + (h01 * p1) + (h11 * (h * v1_mps));
    }

    /**
     * @brief Derivative of Hermite interpolation, giving velocity at tau01.
     * @param p0 Start position
     * @param v0_mps Velocity at start
     * @param p1 End position
     * @param v1_mps Velocity at end
     * @param dt_s Time span of the interval
     * @param tau01 Normalized parameter in [0,1]
     * @return Interpolated velocity (m/s)
     */
    inline Vec3 hermite_velocity_mps(const Vec3 &p0, const Vec3 &v0_mps, const Vec3 &p1, const Vec3 &v1_mps,
                                     const double dt_s, const double tau01)
    {
        if (!(dt_s != 0.0) || !std::isfinite(dt_s))
        {
            return Vec3{0.0, 0.0, 0.0};
        }

        const double t = clamp01(tau01);
        const double t2 = t * t;
        const double h = dt_s;

        const double h00p = 6.0 * t2 - 6.0 * t;
        const double h10p = 3.0 * t2 - 4.0 * t + 1.0;
        const double h01p = -6.0 * t2 + 6.0 * t;
        const double h11p = 3.0 * t2 - 2.0 * t;

        const Vec3 dp_dtau = (h00p * p0) + (h10p * (h * v0_mps)) + (h01p * p1) + (h11p * (h * v1_mps));
        return dp_dtau / h;
    }

    struct OrbitalElements
    {
        double semi_major_axis_m{std::numeric_limits<double>::infinity()};
        double eccentricity{0.0};
        double inclination_rad{0.0};
        double raan_rad{0.0};
        double arg_periapsis_rad{0.0};
        double true_anomaly_rad{0.0};
        double mean_anomaly_rad{std::numeric_limits<double>::quiet_NaN()};
    };

    struct KeplerAnomalyResult
    {
        double mean_anomaly_rad{std::numeric_limits<double>::quiet_NaN()};
        double eccentric_anomaly_rad{std::numeric_limits<double>::quiet_NaN()}; // Elliptic only
        double hyperbolic_anomaly_rad{std::numeric_limits<double>::quiet_NaN()}; // Hyperbolic only
        bool elliptic{false};
        bool hyperbolic{false};
        bool valid{false};
    };

    /**
     * @brief Convert true anomaly to eccentric/hyperbolic and mean anomaly.
     * @param e Eccentricity (must be >= 0)
     * @param true_anomaly_rad True anomaly in radians
     * @return Anomaly values; check `elliptic` or `hyperbolic` flag for orbit type.
     *         Parabolic (e~1) returns invalid result.
     */
    inline KeplerAnomalyResult kepler_anomalies_from_true_anomaly(const double e, const double true_anomaly_rad)
    {
        KeplerAnomalyResult out;
        if (!(e >= 0.0) || !std::isfinite(e) || !std::isfinite(true_anomaly_rad))
        {
            return out;
        }

        const double nu = true_anomaly_rad;
        const double half_nu = 0.5 * nu;
        const double s = std::sin(half_nu);
        const double c = std::cos(half_nu);

        const double eps = 1e-12;
        if (e < 1.0 - eps)
        {
            // Elliptic: tan(nu/2) = sqrt((1+e)/(1-e)) * tan(E/2)
            const double k = std::sqrt((1.0 - e) / (1.0 + e));
            const double E = 2.0 * std::atan2(k * s, c);
            const double M = E - e * std::sin(E);
            out.eccentric_anomaly_rad = E;
            out.mean_anomaly_rad = M;
            out.elliptic = true;
            out.valid = std::isfinite(M);
            return out;
        }

        if (e > 1.0 + eps)
        {
            // Hyperbolic: tanh(H/2) = sqrt((e-1)/(e+1)) * tan(nu/2)
            const double t = std::tan(half_nu);
            const double k = std::sqrt((e - 1.0) / (e + 1.0));
            const double x = k * t;
            const double H = 2.0 * std::atanh(x);
            const double M = e * std::sinh(H) - H;
            out.hyperbolic_anomaly_rad = H;
            out.mean_anomaly_rad = M;
            out.hyperbolic = true;
            out.valid = std::isfinite(M);
            return out;
        }

        // Parabolic (e ~= 1): leave NaN for now.
        return out;
    }

    /**
     * @brief Compute classical orbital elements from relative state vector.
     *
     * Uses inertial Z-axis as reference for inclination and RAAN.
     *
     * @param mu_m3_s2 Gravitational parameter (G*M) of the primary
     * @param r_rel_m Position relative to primary (m)
     * @param v_rel_mps Velocity relative to primary (m/s)
     * @return Orbital elements; check individual fields for validity.
     */
    inline OrbitalElements orbital_elements_from_relative_state(const double mu_m3_s2, const Vec3 &r_rel_m,
                                                                const Vec3 &v_rel_mps)
    {
        OrbitalElements el;
        if (!(mu_m3_s2 > 0.0) || !std::isfinite(mu_m3_s2))
        {
            return el;
        }

        const double r = glm::length(r_rel_m);
        const double v2 = glm::dot(v_rel_mps, v_rel_mps);
        if (!(r > 0.0) || !std::isfinite(r) || !std::isfinite(v2))
        {
            return el;
        }

        const Vec3 h = glm::cross(r_rel_m, v_rel_mps);
        const double hmag = glm::length(h);
        if (!(hmag > 0.0) || !std::isfinite(hmag))
        {
            el.true_anomaly_rad = wrap_angle_0_2pi(std::atan2(r_rel_m.y, r_rel_m.x));
            return el;
        }

        const Vec3 k{0.0, 0.0, 1.0};
        const Vec3 n = glm::cross(k, h);
        const double nmag = glm::length(n);

        const Vec3 evec = (glm::cross(v_rel_mps, h) / mu_m3_s2) - (r_rel_m / r);
        const double e = glm::length(evec);
        el.eccentricity = (std::isfinite(e) ? e : 0.0);

        const double energy = 0.5 * v2 - mu_m3_s2 / r;
        if (std::abs(energy) > 1e-30 && std::isfinite(energy))
        {
            el.semi_major_axis_m = -mu_m3_s2 / (2.0 * energy);
        }

        el.inclination_rad = wrap_angle_0_2pi(std::acos(std::clamp(h.z / hmag, -1.0, 1.0)));

        if (nmag > 1e-24 && std::isfinite(nmag))
        {
            const double raan = std::acos(std::clamp(n.x / nmag, -1.0, 1.0));
            el.raan_rad = wrap_angle_0_2pi((n.y >= 0.0) ? raan : (2.0 * std::acos(-1.0) - raan));
        }
        else
        {
            el.raan_rad = 0.0;
        }

        const double eps_e = 1e-10;
        if (el.eccentricity > eps_e && nmag > 1e-24)
        {
            const double argp = std::acos(std::clamp(glm::dot(n, evec) / (nmag * el.eccentricity), -1.0, 1.0));
            el.arg_periapsis_rad = wrap_angle_0_2pi((evec.z >= 0.0) ? argp : (2.0 * std::acos(-1.0) - argp));
        }
        else
        {
            el.arg_periapsis_rad = 0.0;
        }

        if (el.eccentricity > eps_e)
        {
            const double ta = std::acos(std::clamp(glm::dot(evec, r_rel_m) / (el.eccentricity * r), -1.0, 1.0));
            el.true_anomaly_rad =
                    wrap_angle_0_2pi((glm::dot(r_rel_m, v_rel_mps) >= 0.0) ? ta : (2.0 * std::acos(-1.0) - ta));
        }
        else if (nmag > 1e-24)
        {
            const double u = std::acos(std::clamp(glm::dot(n, r_rel_m) / (nmag * r), -1.0, 1.0));
            el.true_anomaly_rad = wrap_angle_0_2pi((r_rel_m.z >= 0.0) ? u : (2.0 * std::acos(-1.0) - u));
        }
        else
        {
            el.true_anomaly_rad = wrap_angle_0_2pi(std::atan2(r_rel_m.y, r_rel_m.x));
        }

        const KeplerAnomalyResult anom = kepler_anomalies_from_true_anomaly(el.eccentricity, el.true_anomaly_rad);
        el.mean_anomaly_rad = anom.mean_anomaly_rad;

        return el;
    }

    /**
     * @brief Compute orbital elements using a custom reference axis for inclination.
     *
     * Useful when the primary's spin axis differs from the inertial Z-axis.
     *
     * @param mu_m3_s2 Gravitational parameter
     * @param r_rel_m Position relative to primary
     * @param v_rel_mps Velocity relative to primary
     * @param ref_axis_unit_i Reference axis (e.g., primary's spin axis) in inertial frame
     */
    inline OrbitalElements orbital_elements_from_relative_state_about_axis(const double mu_m3_s2, const Vec3 &r_rel_m,
                                                                           const Vec3 &v_rel_mps,
                                                                           const Vec3 &ref_axis_unit_i)
    {
        OrbitalElements el;
        if (!(mu_m3_s2 > 0.0) || !std::isfinite(mu_m3_s2))
        {
            return el;
        }

        const Vec3 k = normalized_or(ref_axis_unit_i, Vec3{0.0, 0.0, 0.0});
        const double k2 = glm::dot(k, k);
        if (!(k2 > 0.0) || !std::isfinite(k2))
        {
            return el;
        }

        const double r = glm::length(r_rel_m);
        const double v2 = glm::dot(v_rel_mps, v_rel_mps);
        if (!(r > 0.0) || !std::isfinite(r) || !std::isfinite(v2))
        {
            return el;
        }

        const Vec3 h = glm::cross(r_rel_m, v_rel_mps);
        const double hmag = glm::length(h);
        if (!(hmag > 0.0) || !std::isfinite(hmag))
        {
            return el;
        }
        const Vec3 hhat = h / hmag;

        const Vec3 evec = (glm::cross(v_rel_mps, h) / mu_m3_s2) - (r_rel_m / r);
        const double e = glm::length(evec);
        el.eccentricity = (std::isfinite(e) ? e : 0.0);

        const double energy = 0.5 * v2 - mu_m3_s2 / r;
        if (std::abs(energy) > 1e-30 && std::isfinite(energy))
        {
            el.semi_major_axis_m = -mu_m3_s2 / (2.0 * energy);
        }

        el.inclination_rad = wrap_angle_0_2pi(std::acos(std::clamp(glm::dot(hhat, k), -1.0, 1.0)));

        // Build a deterministic reference direction in the reference plane.
        Vec3 ref_x = Vec3{1.0, 0.0, 0.0} - glm::dot(Vec3{1.0, 0.0, 0.0}, k) * k;
        if (glm::dot(ref_x, ref_x) <= 1e-24 || !std::isfinite(glm::dot(ref_x, ref_x)))
        {
            ref_x = Vec3{0.0, 1.0, 0.0} - glm::dot(Vec3{0.0, 1.0, 0.0}, k) * k;
        }
        ref_x = normalized_or(ref_x, Vec3{1.0, 0.0, 0.0});
        const Vec3 ref_y = normalized_or(glm::cross(k, ref_x), Vec3{0.0, 1.0, 0.0});

        const Vec3 n = glm::cross(k, h);
        const double nmag = glm::length(n);

        const double eps = 1e-10;
        Vec3 nhat{0.0, 0.0, 0.0};
        if (nmag > 1e-24 && std::isfinite(nmag))
        {
            nhat = n / nmag;
            const double x = glm::dot(nhat, ref_x);
            const double y = glm::dot(nhat, ref_y);
            el.raan_rad = wrap_angle_0_2pi(std::atan2(y, x));
        }
        else
        {
            el.raan_rad = 0.0;
        }

        if (el.eccentricity > eps && nmag > 1e-24)
        {
            const Vec3 ehat = evec / el.eccentricity;
            const double x = glm::dot(ehat, nhat);
            const double y = glm::dot(glm::cross(nhat, ehat), hhat);
            el.arg_periapsis_rad = wrap_angle_0_2pi(std::atan2(y, x));

            const Vec3 rhat = r_rel_m / r;
            const double xt = glm::dot(rhat, ehat);
            const double yt = glm::dot(glm::cross(ehat, rhat), hhat);
            el.true_anomaly_rad = wrap_angle_0_2pi(std::atan2(yt, xt));
        }
        else if (nmag > 1e-24)
        {
            // Circular: use argument of latitude u measured from node.
            const Vec3 rhat = r_rel_m / r;
            const double x = glm::dot(rhat, nhat);
            const double y = glm::dot(glm::cross(nhat, rhat), hhat);
            el.arg_periapsis_rad = 0.0;
            el.true_anomaly_rad = wrap_angle_0_2pi(std::atan2(y, x));
        }
        else
        {
            // Equatorial circular: use inertial angle in plane.
            const Vec3 r_proj = r_rel_m - glm::dot(r_rel_m, k) * k;
            const double rproj = glm::length(r_proj);
            if (rproj > 0.0 && std::isfinite(rproj))
            {
                const Vec3 rhat = r_proj / rproj;
                const double x = glm::dot(rhat, ref_x);
                const double y = glm::dot(rhat, ref_y);
                el.raan_rad = 0.0;
                el.arg_periapsis_rad = 0.0;
                el.true_anomaly_rad = wrap_angle_0_2pi(std::atan2(y, x));
            }
        }

        const KeplerAnomalyResult anom = kepler_anomalies_from_true_anomaly(el.eccentricity, el.true_anomaly_rad);
        el.mean_anomaly_rad = anom.mean_anomaly_rad;
        return el;
    }

    inline OrbitalElements orbital_elements_from_state(const double mu_m3_s2, const State &primary_state_in,
                                                       const State &sc_state_in)
    {
        return orbital_elements_from_relative_state(mu_m3_s2, sc_state_in.position_m - primary_state_in.position_m,
                                                    sc_state_in.velocity_mps - primary_state_in.velocity_mps);
    }

    inline OrbitalElements orbital_elements_from_state_about_axis(const double mu_m3_s2, const State &primary_state_in,
                                                                  const State &sc_state_in, const Vec3 &ref_axis_unit_i)
    {
        return orbital_elements_from_relative_state_about_axis(
                mu_m3_s2, sc_state_in.position_m - primary_state_in.position_m,
                sc_state_in.velocity_mps - primary_state_in.velocity_mps, ref_axis_unit_i);
    }

    /**
     * @brief Convert orbital elements back to Cartesian state vector.
     *
     * Inverse of orbital_elements_from_relative_state. Returns position/velocity
     * relative to primary in inertial frame.
     */
    inline State relative_state_from_orbital_elements(const double mu_m3_s2, const OrbitalElements &el)
    {
        State out{};
        if (!(mu_m3_s2 > 0.0) || !std::isfinite(mu_m3_s2))
        {
            return out;
        }
        if (!std::isfinite(el.semi_major_axis_m) || !(std::abs(el.semi_major_axis_m) > 0.0))
        {
            return out;
        }
        const double e = std::max(0.0, el.eccentricity);
        const double a = el.semi_major_axis_m;
        const double p = a * (1.0 - e * e);
        if (!(p > 0.0) || !std::isfinite(p))
        {
            return out;
        }

        const double ta = el.true_anomaly_rad;
        if (!std::isfinite(ta))
        {
            return out;
        }

        const double cta = std::cos(ta);
        const double sta = std::sin(ta);
        const double denom = 1.0 + e * cta;
        if (!(denom > 0.0) || !std::isfinite(denom))
        {
            return out;
        }

        const double r = p / denom;
        if (!(r > 0.0) || !std::isfinite(r))
        {
            return out;
        }

        // Perifocal coordinates (PQW), z=0.
        const Vec3 r_pf{r * cta, r * sta, 0.0};
        const double sqrt_mu_over_p = std::sqrt(mu_m3_s2 / p);
        const Vec3 v_pf{-sqrt_mu_over_p * sta, sqrt_mu_over_p * (e + cta), 0.0};

        // Basis vectors of the perifocal frame expressed in inertial coordinates.
        const double cO = std::cos(el.raan_rad);
        const double sO = std::sin(el.raan_rad);
        const double ci = std::cos(el.inclination_rad);
        const double si = std::sin(el.inclination_rad);
        const double cw = std::cos(el.arg_periapsis_rad);
        const double sw = std::sin(el.arg_periapsis_rad);

        const Vec3 P{cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si};
        const Vec3 Q{-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si};

        out.position_m = r_pf.x * P + r_pf.y * Q;
        out.velocity_mps = v_pf.x * P + v_pf.y * Q;
        return out;
    }

    inline State state_from_orbital_elements(const double mu_m3_s2, const State &primary_state_in,
                                             const OrbitalElements &el)
    {
        State rel = relative_state_from_orbital_elements(mu_m3_s2, el);
        rel.position_m += primary_state_in.position_m;
        rel.velocity_mps += primary_state_in.velocity_mps;
        return rel;
    }

    struct OrbitScalars
    {
        double periapsis_radius_m{std::numeric_limits<double>::infinity()};
        double apoapsis_radius_m{std::numeric_limits<double>::infinity()};
        double period_s{std::numeric_limits<double>::infinity()};
        double mean_motion_radps{0.0};
        bool valid{false};
    };

    inline OrbitScalars orbit_scalars_from_elements(const double mu_m3_s2, const OrbitalElements &el)
    {
        OrbitScalars out;
        if (!(mu_m3_s2 > 0.0) || !std::isfinite(mu_m3_s2))
        {
            return out;
        }
        const double a = el.semi_major_axis_m;
        const double e = el.eccentricity;
        if (!std::isfinite(a) || !std::isfinite(e))
        {
            return out;
        }

        out.periapsis_radius_m = a * (1.0 - e);
        if (e < 1.0)
        {
            out.apoapsis_radius_m = a * (1.0 + e);
            if (a > 0.0)
            {
                const double n = std::sqrt(mu_m3_s2 / (a * a * a));
                out.mean_motion_radps = (std::isfinite(n) ? n : 0.0);
                const double two_pi = 2.0 * std::acos(-1.0);
                out.period_s = (out.mean_motion_radps > 0.0) ? (two_pi / out.mean_motion_radps)
                                                             : std::numeric_limits<double>::infinity();
            }
        }
        out.valid = std::isfinite(out.periapsis_radius_m);
        return out;
    }

    struct OrbitApsides
    {
        Vec3 periapsis_rel_m{0.0, 0.0, 0.0};
        Vec3 apoapsis_rel_m{0.0, 0.0, 0.0};
        double periapsis_radius_m{std::numeric_limits<double>::infinity()};
        double apoapsis_radius_m{std::numeric_limits<double>::infinity()};
        bool has_apoapsis{false};
        bool valid{false};
    };

    /**
     * @brief Compute periapsis and apoapsis positions from state.
     *
     * For hyperbolic orbits, only periapsis is valid (has_apoapsis = false).
     * For near-circular orbits, periapsis direction defaults to current position.
     */
    inline OrbitApsides apsides_from_relative_state(const double mu_m3_s2, const Vec3 &r_rel_m, const Vec3 &v_rel_mps)
    {
        OrbitApsides out;
        if (!(mu_m3_s2 > 0.0) || !std::isfinite(mu_m3_s2))
        {
            return out;
        }

        const double r = glm::length(r_rel_m);
        if (!(r > 0.0) || !std::isfinite(r))
        {
            return out;
        }

        const Vec3 h = glm::cross(r_rel_m, v_rel_mps);
        const double h2 = glm::dot(h, h);
        if (!(h2 > 0.0) || !std::isfinite(h2))
        {
            return out;
        }

        const double p = h2 / mu_m3_s2;
        if (!(p > 0.0) || !std::isfinite(p))
        {
            return out;
        }

        // Eccentricity vector points toward periapsis.
        const Vec3 evec = (glm::cross(v_rel_mps, h) / mu_m3_s2) - (r_rel_m / r);
        const double e = glm::length(evec);
        if (!std::isfinite(e))
        {
            return out;
        }

        const double rp = p / (1.0 + e);
        if (!(rp > 0.0) || !std::isfinite(rp))
        {
            return out;
        }

        // For near-circular orbits, periapsis direction is undefined; use current radius direction so markers behave.
        const double eps_e = 1e-10;
        const Vec3 ehat = (e > eps_e) ? normalized_or(evec, r_rel_m / r) : (r_rel_m / r);

        out.periapsis_radius_m = rp;
        out.periapsis_rel_m = rp * ehat;

        if (e < 1.0 - eps_e)
        {
            const double ra = p / (1.0 - e);
            if (ra > 0.0 && std::isfinite(ra))
            {
                out.has_apoapsis = true;
                out.apoapsis_radius_m = ra;
                out.apoapsis_rel_m = -ra * ehat;
            }
        }

        out.valid = true;
        return out;
    }

    inline OrbitApsides apsides_from_state(const double mu_m3_s2, const State &primary_state_in,
                                           const State &sc_state_in)
    {
        OrbitApsides out = apsides_from_relative_state(mu_m3_s2, sc_state_in.position_m - primary_state_in.position_m,
                                                       sc_state_in.velocity_mps - primary_state_in.velocity_mps);
        if (!out.valid)
        {
            return out;
        }
        out.periapsis_rel_m += primary_state_in.position_m;
        out.apoapsis_rel_m += primary_state_in.position_m;
        return out;
    }

    /**
     * @brief Stumpff function C(z) used in universal variable formulation.
     *
     * C(z) = (1 - cos(sqrt(z))) / z  for z > 0 (elliptic)
     * C(z) = (cosh(sqrt(-z)) - 1) / (-z)  for z < 0 (hyperbolic)
     * Series expansion near z = 0.
     */
    inline double stumpff_C(const double z)
    {
        if (std::abs(z) < 1e-8)
        {
            // Series: C(z) = 1/2 - z/24 + z^2/720 - ...
            const double z2 = z * z;
            return 0.5 - (z / 24.0) + (z2 / 720.0);
        }
        if (z > 0.0)
        {
            const double sz = std::sqrt(z);
            return (1.0 - std::cos(sz)) / z;
        }
        const double sz = std::sqrt(-z);
        return (std::cosh(sz) - 1.0) / (-z);
    }

    /**
     * @brief Stumpff function S(z) used in universal variable formulation.
     *
     * S(z) = (sqrt(z) - sin(sqrt(z))) / z^(3/2)  for z > 0 (elliptic)
     * S(z) = (sinh(sqrt(-z)) - sqrt(-z)) / (-z)^(3/2)  for z < 0 (hyperbolic)
     * Series expansion near z = 0.
     */
    inline double stumpff_S(const double z)
    {
        if (std::abs(z) < 1e-8)
        {
            // Series: S(z) = 1/6 - z/120 + z^2/5040 - ...
            const double z2 = z * z;
            return (1.0 / 6.0) - (z / 120.0) + (z2 / 5040.0);
        }
        if (z > 0.0)
        {
            const double sz = std::sqrt(z);
            return (sz - std::sin(sz)) / (sz * sz * sz);
        }
        const double sz = std::sqrt(-z);
        return (std::sinh(sz) - sz) / (sz * sz * sz);
    }

} // namespace orbitsim
