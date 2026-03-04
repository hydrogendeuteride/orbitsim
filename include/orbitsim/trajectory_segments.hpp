#pragma once

#include "orbitsim/types.hpp"

#include <cstddef>
#include <cstdint>

namespace orbitsim
{

    struct TrajectorySegment
    {
        double t0_s{0.0};
        double dt_s{0.0};
        State start{};
        State end{};
        std::uint32_t flags{0};
    };

    struct TrajectorySegmentOptions
    {
        double duration_s{3600.0};
        std::size_t max_segments{2048};
        bool include_start{true};
        bool include_end{true};
        bool stop_on_impact{false};
        double lookup_dt_s{0.0};
    };

} // namespace orbitsim
