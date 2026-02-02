#pragma once

#include "orbitsim/types.hpp"

namespace orbitsim
{

    struct TrajectorySample
    {
        double t_s{0.0};
        Vec3 position_m{0.0, 0.0, 0.0};
        Vec3 velocity_mps{0.0, 0.0, 0.0};
    };

} // namespace orbitsim

