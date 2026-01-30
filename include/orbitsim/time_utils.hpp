#pragma once

namespace orbitsim
{

    /// @brief Convert seconds to seconds (identity, for consistency).
    inline constexpr double seconds(const double s) { return s; }

    /// @brief Convert minutes to seconds.
    inline constexpr double minutes(const double m) { return m * 60.0; }

    /// @brief Convert hours to seconds.
    inline constexpr double hours(const double h) { return h * 3600.0; }

    /// @brief Convert days to seconds.
    inline constexpr double days(const double d) { return d * 86400.0; }

} // namespace orbitsim
