#pragma once

#include <orbitsim/orbitsim.hpp>

#include <algorithm>
#include <cmath>

inline bool near_abs(const double a, const double b, const double abs_tol) { return std::abs(a - b) <= abs_tol; }

inline bool near_rel(const double a, const double b, const double rel_tol, const double abs_floor = 0.0)
{
    const double scale = std::max(abs_floor, std::max(std::abs(a), std::abs(b)));
    return std::abs(a - b) <= rel_tol * scale;
}

inline bool near_vec_abs(const orbitsim::Vec3 &a, const orbitsim::Vec3 &b, const double abs_tol)
{
    return near_abs(a.x, b.x, abs_tol) && near_abs(a.y, b.y, abs_tol) && near_abs(a.z, b.z, abs_tol);
}
