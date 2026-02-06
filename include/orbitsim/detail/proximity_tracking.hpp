#pragma once

#include "orbitsim/events.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <vector>

namespace orbitsim::detail
{

    /** @brief Resolved proximity detection configuration. */
    struct ProximityConfig
    {
        bool enabled{false};
        SpacecraftId center_id{kInvalidSpacecraftId};
        double enter_radius_m{0.0};
        double exit_radius_m{0.0};
    };

    /** @brief Validate and normalize proximity options into a resolved config. */
    inline ProximityConfig resolve_proximity_config(const bool enable, const SpacecraftId center_spacecraft_id,
                                                    const double enter_radius_m, const double exit_radius_m)
    {
        ProximityConfig cfg;
        cfg.enabled = enable && (center_spacecraft_id != kInvalidSpacecraftId) && (enter_radius_m > 0.0) &&
                      std::isfinite(enter_radius_m) && (exit_radius_m > 0.0) && std::isfinite(exit_radius_m);

        if (!cfg.enabled)
        {
            return cfg;
        }

        cfg.center_id = center_spacecraft_id;
        cfg.enter_radius_m = enter_radius_m;
        cfg.exit_radius_m = std::max(enter_radius_m, exit_radius_m);
        return cfg;
    }

    /**
     * @brief Initialize the proximity active set from current spacecraft positions.
     *
     * Adds any spacecraft (other than the center) that is within exit_radius_m of the center.
     *
     * @param spacecraft All spacecraft in the simulation
     * @param center_pos_m Position of the center spacecraft
     * @param center_id ID of the center spacecraft
     * @param exit_radius_m Exit radius threshold
     * @param[out] active Set to populate with spacecraft within range
     */
    inline void initialize_proximity_active_set(const std::vector<Spacecraft> &spacecraft, const Vec3 &center_pos_m,
                                                const SpacecraftId center_id, const double exit_radius_m,
                                                std::unordered_set<SpacecraftId> &active)
    {
        for (const auto &sc : spacecraft)
        {
            if (sc.id == center_id)
            {
                continue;
            }
            const double d = glm::length(sc.state.position_m - center_pos_m);
            if (d <= exit_radius_m)
            {
                active.insert(sc.id);
            }
        }
    }

    /** @brief Update the proximity active set after a proximity event. */
    inline void update_proximity_active_set(const SpacecraftId spacecraft_id, const Crossing crossing,
                                            std::unordered_set<SpacecraftId> &active)
    {
        if (crossing == Crossing::Enter)
        {
            active.insert(spacecraft_id);
        }
        else
        {
            active.erase(spacecraft_id);
        }
    }

} // namespace orbitsim::detail
