#pragma once

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

    struct RtnFrame
    {
        Vec3 R{1.0, 0.0, 0.0};
        Vec3 T{0.0, 1.0, 0.0};
        Vec3 N{0.0, 0.0, 1.0};
    };

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

    inline double wrap_angle_0_2pi(const double rad)
    {
        if (!std::isfinite(rad))
        {
            return 0.0;
        }
        const double two_pi = 2.0 * std::acos(-1.0);
        double x = std::fmod(rad, two_pi);
        if (x < 0.0)
        {
            x += two_pi;
        }
        return x;
    }

    struct OrbitalElements
    {
        double semi_major_axis_m{std::numeric_limits<double>::infinity()};
        double eccentricity{0.0};
        double inclination_rad{0.0};
        double raan_rad{0.0};
        double arg_periapsis_rad{0.0};
        double true_anomaly_rad{0.0};
    };

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

        return el;
    }

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
