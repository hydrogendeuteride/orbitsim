#include <gtest/gtest.h>
#include "test_helpers.hpp"

TEST(Maneuvers, ImpulseChangesVelocity)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    const orbitsim::GameSimulation::BodyHandle body_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 0.0,
            .radius_m = 1.0,
            .state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}),
    });
    ASSERT_TRUE(body_h.valid());

    const orbitsim::GameSimulation::SpacecraftHandle sc_h = sim.create_spacecraft(orbitsim::Spacecraft{
            .state = orbitsim::make_state({10.0, 0.0, 0.0}, {0.0, 1.0, 0.0}),
            .dry_mass_kg = 1.0,
            .prop_mass_kg = 0.0,
    });
    ASSERT_TRUE(sc_h.valid());

    constexpr double t_impulse_s = 1.0;
    constexpr double dv_prograde_mps = 5.0;
    sim.maneuver_plan().impulses.push_back(
            orbitsim::impulse().time(t_impulse_s).prograde(dv_prograde_mps).primary(body_h).spacecraft(sc_h));

    sim.step(2.0);

    const orbitsim::Spacecraft *out = sim.spacecraft_by_id(sc_h.id);
    ASSERT_NE(out, nullptr);

    const orbitsim::Vec3 r_rel_m = orbitsim::Vec3{10.0, 0.0, 0.0} + orbitsim::Vec3{0.0, 1.0, 0.0} * t_impulse_s;
    const orbitsim::Vec3 v_rel_mps{0.0, 1.0, 0.0};
    const orbitsim::RtnFrame rtn = orbitsim::compute_rtn_frame(r_rel_m, v_rel_mps);
    const orbitsim::Vec3 expected_vel_mps = v_rel_mps + dv_prograde_mps * rtn.T;

    EXPECT_TRUE(near_vec_abs(out->state.velocity_mps, expected_vel_mps, 1e-9));
}

