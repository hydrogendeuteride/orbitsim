#pragma once

#include "orbitsim/types.hpp"

#include <functional>
#include <optional>

namespace orbitsim
{

    /// @brief Optional callback to look up a spacecraft inertial state at a given time.
    ///
    /// Used primarily for LVLH/relative-frame maneuvers and plotting. Implementations may provide a cached/interpolated
    /// state (recommended) to avoid repeatedly re-integrating the target trajectory.
    using SpacecraftStateLookup = std::function<std::optional<State>(SpacecraftId, double /*t_s*/)>;

} // namespace orbitsim

