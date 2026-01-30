#pragma once

#include <glm/glm.hpp>

#include <cstdint>
#include <vector>

namespace orbitsim
{

    using Vec3 = glm::dvec3;
    using BodyId = std::uint32_t;
    inline constexpr BodyId kInvalidBodyId = 0;
    using SpacecraftId = std::uint32_t;
    inline constexpr SpacecraftId kInvalidSpacecraftId = 0;

    struct SpinState
    {
        Vec3 axis{0.0, 1.0, 0.0}; // Unit vector
        double angle_rad{0.0}; // Radians
        double rate_rad_per_s{0.0}; // Radians / second
    };

    struct State
    {
        Vec3 position_m{0.0, 0.0, 0.0};
        Vec3 velocity_mps{0.0, 0.0, 0.0};
        SpinState spin{};
    };

    struct MassiveBody
    {
        double mass_kg{0.0};
        double radius_m{0.0};
        double atmosphere_top_height_m{0.0};
        double terrain_max_height_m{0.0};
        double soi_radius_m{0.0};
        State state{};
        BodyId id{kInvalidBodyId};
    };

    struct Engine
    {
        double max_thrust_N{0.0};
        double isp_s{0.0};
        double min_throttle_0_1{0.0};
    };

    struct Spacecraft
    {
        State state{};
        double dry_mass_kg{0.0};
        double prop_mass_kg{0.0};
        std::vector<Engine> engines{};
        SpacecraftId id{kInvalidSpacecraftId};

        inline double mass_kg() const { return dry_mass_kg + prop_mass_kg; }
    };

    // -------------------------------------------------------------------------
    // Factory helpers for concise object creation
    // -------------------------------------------------------------------------

    /// @brief Create a SpinState with the given parameters.
    /// @param axis Rotation axis (will be used as-is; caller should normalize if needed).
    /// @param rate_rad_per_s Rotation rate [rad/s].
    /// @param angle_rad Initial rotation angle [rad]. Default is 0.
    inline SpinState make_spin(const Vec3 &axis, const double rate_rad_per_s, const double angle_rad = 0.0)
    {
        return SpinState{.axis = axis, .angle_rad = angle_rad, .rate_rad_per_s = rate_rad_per_s};
    }

    /// @brief Create a State with position, velocity, and optional spin.
    /// @param position_m Position vector [m].
    /// @param velocity_mps Velocity vector [m/s].
    /// @param spin Optional spin state. Default is no spin.
    inline State make_state(const Vec3 &position_m, const Vec3 &velocity_mps, const SpinState &spin = {})
    {
        return State{.position_m = position_m, .velocity_mps = velocity_mps, .spin = spin};
    }

    /// @brief Create a State with position, velocity, and spin parameters.
    /// @param position_m Position vector [m].
    /// @param velocity_mps Velocity vector [m/s].
    /// @param spin_axis Spin axis vector.
    /// @param spin_rate_rad_per_s Spin rate [rad/s].
    /// @param spin_angle_rad Initial spin angle [rad]. Default is 0.
    inline State make_state(const Vec3 &position_m, const Vec3 &velocity_mps, const Vec3 &spin_axis,
                            const double spin_rate_rad_per_s, const double spin_angle_rad = 0.0)
    {
        return State{.position_m = position_m,
                     .velocity_mps = velocity_mps,
                     .spin = make_spin(spin_axis, spin_rate_rad_per_s, spin_angle_rad)};
    }

} // namespace orbitsim
