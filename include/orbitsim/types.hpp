#pragma once

#include <glm/glm.hpp>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace orbitsim
{

    using Vec3 = glm::dvec3;

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

    struct Pose
    {
        Vec3 position_m{0.0, 0.0, 0.0};
        Vec3 rotation_axis{0.0, 1.0, 0.0};
        double rotation_angle_rad{0.0};
    };

    struct MassiveBody
    {
        double mass_kg{0.0};
        double radius_m{0.0};
        double atmosphere_top_height_m{0.0};
        double terrain_max_height_m{0.0};
        double soi_radius_m{0.0};
        State state{};
    };

    // Convenience alias for "planet/moon" bodies (massive bodies).
    using CelestialBody = MassiveBody;

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

        inline double mass_kg() const { return dry_mass_kg + prop_mass_kg; }
    };

    inline Pose pose_from_state(const State &s)
    {
        return Pose{
                .position_m = s.position_m,
                .rotation_axis = s.spin.axis,
                .rotation_angle_rad = s.spin.angle_rad,
        };
    }

    inline Pose pose_from_body(const MassiveBody &b) { return pose_from_state(b.state); }
    inline Pose pose_from_spacecraft(const Spacecraft &sc) { return pose_from_state(sc.state); }

} // namespace orbitsim
