#include <gtest/gtest.h>
#include "test_helpers.hpp"

TEST(SoiRails, SelectLocalBodyAndUseHysteresis)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 1.0;
    cfg.enable_events = false;

    orbitsim::GameSimulation sim(cfg);

    const orbitsim::GameSimulation::BodyHandle planet_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 1.0e6,
            .radius_m = 1.0,
            .soi_radius_m = 1000.0,
            .state = {.position_m = {0.0, 0.0, 0.0}, .velocity_mps = {0.0, 0.0, 0.0}},
    });
    ASSERT_TRUE(planet_h.valid());

    const orbitsim::GameSimulation::BodyHandle moon_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 1.0e3,
            .radius_m = 1.0,
            .soi_radius_m = 100.0,
            .state = {.position_m = {500.0, 0.0, 0.0}, .velocity_mps = {0.0, 0.0, 0.0}},
    });
    ASSERT_TRUE(moon_h.valid());

    const orbitsim::CelestialEphemeris eph{}; // use live body states

    orbitsim::SoiSwitchOptions opt{};
    opt.enter_scale = 1.0;
    opt.exit_scale = 1.02;
    opt.prefer_smallest_soi = true;
    opt.fallback_to_max_accel = true;

    // Inside both SOIs -> choose the more local (moon).
    const orbitsim::Vec3 pos_inside_moon{550.0, 0.0, 0.0}; // 50m from moon center
    const orbitsim::BodyId primary0 =
            orbitsim::select_primary_body_id_rails(sim, eph, pos_inside_moon, 0.0, planet_h.id, opt);
    EXPECT_EQ(primary0, moon_h.id);

    // Slightly outside moon enter radius but inside exit hysteresis -> keep moon.
    const orbitsim::Vec3 pos_just_outside_enter{601.0, 0.0, 0.0}; // 101m from moon center
    const orbitsim::BodyId primary1 =
            orbitsim::select_primary_body_id_rails(sim, eph, pos_just_outside_enter, 0.0, moon_h.id, opt);
    EXPECT_EQ(primary1, moon_h.id);

    // Outside exit hysteresis -> switch back to planet.
    const orbitsim::Vec3 pos_outside_exit{603.0, 0.0, 0.0}; // 103m from moon center
    const orbitsim::BodyId primary2 =
            orbitsim::select_primary_body_id_rails(sim, eph, pos_outside_exit, 0.0, moon_h.id, opt);
    EXPECT_EQ(primary2, planet_h.id);
}

TEST(SoiRails, FallbackToMaxAccelWhenNoSoiConfigured)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 1.0;
    cfg.enable_events = false;

    orbitsim::GameSimulation sim(cfg);

    const orbitsim::GameSimulation::BodyHandle b0_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 100.0,
            .radius_m = 1.0,
            .soi_radius_m = 0.0,
            .state = {.position_m = {0.0, 0.0, 0.0}, .velocity_mps = {0.0, 0.0, 0.0}},
    });
    const orbitsim::GameSimulation::BodyHandle b1_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 1.0e6,
            .radius_m = 1.0,
            .soi_radius_m = 0.0,
            .state = {.position_m = {100.0, 0.0, 0.0}, .velocity_mps = {0.0, 0.0, 0.0}},
    });
    ASSERT_TRUE(b0_h.valid());
    ASSERT_TRUE(b1_h.valid());

    const orbitsim::CelestialEphemeris eph{};
    orbitsim::SoiSwitchOptions opt{};
    opt.fallback_to_max_accel = true;

    const orbitsim::Vec3 sc_pos{10.0, 0.0, 0.0};
    const orbitsim::BodyId primary =
            orbitsim::select_primary_body_id_rails(sim, eph, sc_pos, 0.0, orbitsim::kInvalidBodyId, opt);
    EXPECT_EQ(primary, b1_h.id);
}

