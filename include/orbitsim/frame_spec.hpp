#pragma once

#include "orbitsim/types.hpp"

namespace orbitsim
{

    /**
     * @brief Identifies a coordinate frame for representing states/trajectories.
     *
     * This is used for visualization, analysis, and maneuver planning; it does not change the simulation's
     * underlying dynamics integration, which always occurs in the library inertial frame.
     */
    enum class TrajectoryFrameType
    {
        Inertial,
        BodyCenteredInertial,
        BodyFixed,
        Synodic,
    };

    /**
     * @brief Parameterization for a coordinate frame.
     *
     * - Inertial: the simulation inertial frame.
     * - BodyCenteredInertial: translation to a body's instantaneous state (no rotation).
     * - BodyFixed: body-fixed rotating frame derived from the body's spin state.
     * - Synodic: two-body barycentric co-rotating frame for (A,B) with +X along A->B.
     */
    struct TrajectoryFrameSpec
    {
        TrajectoryFrameType type{TrajectoryFrameType::Inertial};

        /// Primary body id (used by BodyCenteredInertial, BodyFixed, and as body A in Synodic).
        BodyId primary_body_id{kInvalidBodyId};
        /// Secondary body id (used only by Synodic as body B).
        BodyId secondary_body_id{kInvalidBodyId};

        static TrajectoryFrameSpec inertial() { return {}; }

        static TrajectoryFrameSpec body_centered_inertial(const BodyId body_id)
        {
            return TrajectoryFrameSpec{.type = TrajectoryFrameType::BodyCenteredInertial, .primary_body_id = body_id};
        }

        static TrajectoryFrameSpec body_fixed(const BodyId body_id)
        {
            return TrajectoryFrameSpec{.type = TrajectoryFrameType::BodyFixed, .primary_body_id = body_id};
        }

        static TrajectoryFrameSpec synodic(const BodyId body_a_id, const BodyId body_b_id)
        {
            return TrajectoryFrameSpec{.type = TrajectoryFrameType::Synodic,
                                       .primary_body_id = body_a_id,
                                       .secondary_body_id = body_b_id};
        }
    };

} // namespace orbitsim

