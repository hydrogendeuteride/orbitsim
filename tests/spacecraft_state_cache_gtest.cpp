#include <gtest/gtest.h>
#include "test_helpers.hpp"

TEST(SpacecraftStateCache, SplitsSegmentsAtImpulses)
{
    std::vector<orbitsim::MassiveBody> bodies;
    bodies.push_back(orbitsim::MassiveBody{
            .mass_kg = 1.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {0.0, 0.0, 0.0},
                            .velocity_mps = {0.0, 0.0, 0.0},
                    },
            .id = 1,
    });

    orbitsim::CelestialEphemerisSegment eph;
    eph.t0_s = 0.0;
    eph.dt_s = 10.0;
    eph.start = {bodies[0].state};
    eph.end = {bodies[0].state};

    constexpr orbitsim::SpacecraftId sc_id = 42;
    orbitsim::Spacecraft sc0{};
    sc0.id = sc_id;
    sc0.state.position_m = {1.0, 0.0, 0.0};
    sc0.state.velocity_mps = {0.0, 0.0, 0.0};
    sc0.dry_mass_kg = 1.0;
    sc0.prop_mass_kg = 0.0;

    orbitsim::ManeuverPlan plan;
    plan.impulses.push_back(orbitsim::impulse().time(5.0).prograde(10.0).primary(bodies[0].id).spacecraft(sc_id).build());

    orbitsim::DOPRI5Options integ{};
    integ.adaptive = true;
    integ.max_substeps = 256;

    const orbitsim::SpacecraftStateCache<orbitsim::CelestialEphemerisSegment> cache(
            bodies,
            eph,
            plan,
            0.0,
            0.0,
            integ,
            0.0,
            10.0,
            [&](const orbitsim::SpacecraftId id) -> const orbitsim::Spacecraft * { return (id == sc_id) ? &sc0 : nullptr; },
            orbitsim::SpacecraftStateCache<orbitsim::CelestialEphemerisSegment>::Options{.lookup_dt_s = 0.0});

    const orbitsim::SpacecraftStateLookup lookup = cache.lookup();
    const std::optional<orbitsim::State> s_before = lookup(sc_id, 4.0);
    const std::optional<orbitsim::State> s_after = lookup(sc_id, 6.0);
    ASSERT_TRUE(s_before.has_value());
    ASSERT_TRUE(s_after.has_value());

    EXPECT_TRUE(near_vec_abs(s_before->velocity_mps, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_after->velocity_mps, orbitsim::Vec3{0.0, 10.0, 0.0}, 1e-9));
}

