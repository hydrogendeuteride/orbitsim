#include <gtest/gtest.h>
#include "test_helpers.hpp"

TEST(GameSimulation, TargetedBurnOnlyAffectsOneSpacecraft)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    const orbitsim::GameSimulation::BodyHandle body_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 0.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {0.0, 0.0, 0.0},
                            .velocity_mps = {0.0, 0.0, 0.0},
                    },
    });
    ASSERT_TRUE(body_h.valid());

    orbitsim::Spacecraft sc{};
    sc.state.position_m = {3.0, 0.0, 0.0};
    sc.state.velocity_mps = {0.0, 0.0, 0.0};
    sc.dry_mass_kg = 1.0;
    sc.prop_mass_kg = 0.1;
    sc.engines.push_back(orbitsim::Engine{.max_thrust_N = 10.0, .isp_s = 1.0, .min_throttle_0_1 = 0.0});
    const orbitsim::GameSimulation::SpacecraftHandle sc0_h = sim.create_spacecraft(sc);
    const orbitsim::GameSimulation::SpacecraftHandle sc1_h = sim.create_spacecraft(sc);
    ASSERT_TRUE(sc0_h.valid());
    ASSERT_TRUE(sc1_h.valid());

    sim.maneuver_plan().segments.push_back(orbitsim::BurnSegment{
            .t_start_s = 0.0,
            .t_end_s = 10.0,
            .primary_body_id = body_h,
            .dir_rtn_unit = {1.0, 0.0, 0.0},
            .throttle_0_1 = 1.0,
            .engine_index = 0,
            .spacecraft_id = sc1_h,
    });

    sim.step(1.0);

    const orbitsim::Spacecraft *out0 = sim.spacecraft_by_id(sc0_h.id);
    const orbitsim::Spacecraft *out1 = sim.spacecraft_by_id(sc1_h.id);
    ASSERT_NE(out0, nullptr);
    ASSERT_NE(out1, nullptr);
    EXPECT_TRUE(near_abs(out0->prop_mass_kg, 0.1, 1e-12));
    EXPECT_TRUE(near_abs(out0->state.velocity_mps.x, 0.0, 1e-12));
    EXPECT_GE(out1->prop_mass_kg, 0.0);
    EXPECT_TRUE(near_abs(out1->prop_mass_kg, 0.0, 1e-12));
    EXPECT_TRUE(std::isfinite(out1->state.velocity_mps.x));
    EXPECT_GT(out1->state.velocity_mps.x, 0.0);

}

TEST(GameSimulation, LVLHImpulseUsesTargetStateAtImpulseTime)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    const auto primary_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 1.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {0.0, 0.0, 0.0},
                            .velocity_mps = {0.0, 0.0, 0.0},
                    },
    });
    ASSERT_TRUE(primary_h.valid());

    orbitsim::Spacecraft target{};
    target.state.position_m = {1.0, 0.0, 0.0};
    target.state.velocity_mps = {0.0, 1.0, 0.0};
    target.dry_mass_kg = 1.0;
    target.prop_mass_kg = 0.0;
    const auto target_h = sim.create_spacecraft(target);
    ASSERT_TRUE(target_h.valid());

    orbitsim::Spacecraft chaser{};
    chaser.state.position_m = {10.0, 0.0, 0.0};
    chaser.state.velocity_mps = {0.0, 0.0, 0.0};
    chaser.dry_mass_kg = 1.0;
    chaser.prop_mass_kg = 0.0;
    const auto chaser_h = sim.create_spacecraft(chaser);
    ASSERT_TRUE(chaser_h.valid());

    sim.maneuver_plan().impulses.push_back(
            orbitsim::impulse()
                    .time(1.0)
                    .prograde(10.0)
                    .primary(primary_h)
                    .rtn_lvlh(target_h, primary_h)
                    .spacecraft(chaser_h)
                    .build());

    sim.step(2.0);

    const orbitsim::Spacecraft *chaser_out = sim.spacecraft_by_id(chaser_h);
    ASSERT_NE(chaser_out, nullptr);

    // Target at t=1: r=(1,1,0), v=(0,1,0). Prograde (+T) is not aligned with +Y at this time.
    const orbitsim::Vec3 r_rel{1.0, 1.0, 0.0};
    const orbitsim::Vec3 v_rel{0.0, 1.0, 0.0};
    const orbitsim::RtnFrame f = orbitsim::compute_rtn_frame(r_rel, v_rel);
    const orbitsim::Vec3 expected_dv = 10.0 * f.T;

    EXPECT_TRUE(near_vec_abs(chaser_out->state.velocity_mps, expected_dv, 1e-9));
    EXPECT_LT(chaser_out->state.velocity_mps.x, -1.0);
}

TEST(GameSimulation, LVLHImpulseHonorsPrimaryBodyIdInRtnLvlh)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    // Create bodies in an order that makes "auto-select" pick the wrong one when G=0.
    const auto moon_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 1.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {100.0, 0.0, 0.0},
                            .velocity_mps = {0.0, 0.0, 0.0},
                    },
    });
    ASSERT_TRUE(moon_h.valid());

    const auto earth_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 1.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {0.0, 0.0, 0.0},
                            .velocity_mps = {0.0, 0.0, 0.0},
                    },
    });
    ASSERT_TRUE(earth_h.valid());

    orbitsim::Spacecraft target{};
    target.state.position_m = {1.0, 0.0, 0.0};
    target.state.velocity_mps = {0.0, 1.0, 0.0};
    target.dry_mass_kg = 1.0;
    target.prop_mass_kg = 0.0;
    const auto target_h = sim.create_spacecraft(target);
    ASSERT_TRUE(target_h.valid());

    orbitsim::Spacecraft chaser{};
    chaser.state.position_m = {10.0, 0.0, 0.0};
    chaser.state.velocity_mps = {0.0, 0.0, 0.0};
    chaser.dry_mass_kg = 1.0;
    chaser.prop_mass_kg = 0.0;
    const auto chaser_h = sim.create_spacecraft(chaser);
    ASSERT_TRUE(chaser_h.valid());

    sim.maneuver_plan().impulses.push_back(
            orbitsim::impulse()
                    .time(1.0)
                    .prograde(10.0)
                    .rtn_lvlh(target_h, earth_h)
                    .spacecraft(chaser_h)
                    .build());

    sim.step(2.0);

    const orbitsim::Spacecraft *chaser_out = sim.spacecraft_by_id(chaser_h);
    ASSERT_NE(chaser_out, nullptr);

    // At t=1: target r=(1,1,0), v=(0,1,0) in inertial. Primary selection must use Earth, not Moon.
    const orbitsim::Vec3 target_pos{1.0, 1.0, 0.0};
    const orbitsim::Vec3 target_vel{0.0, 1.0, 0.0};

    const orbitsim::Vec3 r_rel_earth = target_pos - orbitsim::Vec3{0.0, 0.0, 0.0};
    const orbitsim::Vec3 r_rel_moon = target_pos - orbitsim::Vec3{100.0, 0.0, 0.0};
    const orbitsim::RtnFrame f_earth = orbitsim::compute_rtn_frame(r_rel_earth, target_vel);
    const orbitsim::RtnFrame f_moon = orbitsim::compute_rtn_frame(r_rel_moon, target_vel);

    const orbitsim::Vec3 expected_dv = 10.0 * f_earth.T;
    const orbitsim::Vec3 wrong_dv = 10.0 * f_moon.T;

    EXPECT_TRUE(near_vec_abs(chaser_out->state.velocity_mps, expected_dv, 1e-9));
    EXPECT_FALSE(near_vec_abs(chaser_out->state.velocity_mps, wrong_dv, 1e-6));
    EXPECT_LT(chaser_out->state.velocity_mps.x, -1.0);
}

TEST(GameSimulation, LVLHImpulseHonorsPrimaryBodyIdInTrajectoryFrameSpec)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    // Same setup as the previous test: Moon is created first so auto-select would pick it when G=0.
    const auto moon_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 1.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {100.0, 0.0, 0.0},
                            .velocity_mps = {0.0, 0.0, 0.0},
                    },
    });
    ASSERT_TRUE(moon_h.valid());

    const auto earth_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 1.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {0.0, 0.0, 0.0},
                            .velocity_mps = {0.0, 0.0, 0.0},
                    },
    });
    ASSERT_TRUE(earth_h.valid());

    orbitsim::Spacecraft target{};
    target.state.position_m = {1.0, 0.0, 0.0};
    target.state.velocity_mps = {0.0, 1.0, 0.0};
    target.dry_mass_kg = 1.0;
    target.prop_mass_kg = 0.0;
    const auto target_h = sim.create_spacecraft(target);
    ASSERT_TRUE(target_h.valid());

    orbitsim::Spacecraft chaser{};
    chaser.state.position_m = {10.0, 0.0, 0.0};
    chaser.state.velocity_mps = {0.0, 0.0, 0.0};
    chaser.dry_mass_kg = 1.0;
    chaser.prop_mass_kg = 0.0;
    const auto chaser_h = sim.create_spacecraft(chaser);
    ASSERT_TRUE(chaser_h.valid());

    sim.maneuver_plan().impulses.push_back(
            orbitsim::impulse()
                    .time(1.0)
                    .prograde(10.0)
                    .rtn_frame(orbitsim::TrajectoryFrameSpec::lvlh(target_h, earth_h))
                    .spacecraft(chaser_h)
                    .build());

    sim.step(2.0);

    const orbitsim::Spacecraft *chaser_out = sim.spacecraft_by_id(chaser_h);
    ASSERT_NE(chaser_out, nullptr);

    const orbitsim::Vec3 target_pos{1.0, 1.0, 0.0};
    const orbitsim::Vec3 target_vel{0.0, 1.0, 0.0};
    const orbitsim::Vec3 r_rel_earth = target_pos - orbitsim::Vec3{0.0, 0.0, 0.0};
    const orbitsim::RtnFrame f_earth = orbitsim::compute_rtn_frame(r_rel_earth, target_vel);
    const orbitsim::Vec3 expected_dv = 10.0 * f_earth.T;

    EXPECT_TRUE(near_vec_abs(chaser_out->state.velocity_mps, expected_dv, 1e-9));
    EXPECT_LT(chaser_out->state.velocity_mps.x, -1.0);
}

TEST(GameSimulation, StepReturnsImpactEvent)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = true;
    cfg.events.time_tol_s = 1e-6;
    cfg.events.dist_tol_m = 1e-6;
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
            .state = orbitsim::make_state({2.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}),
            .dry_mass_kg = 1.0,
            .prop_mass_kg = 0.0,
    });
    ASSERT_TRUE(sc_h.valid());

    std::vector<orbitsim::Event> events;
    // Crossing should be refined near t=1.
    sim.step(2.0, &events);

    ASSERT_FALSE(events.empty());
    EXPECT_EQ(events.front().type, orbitsim::EventType::Impact);
    EXPECT_EQ(events.front().body_id, body_h.id);
    EXPECT_EQ(events.front().spacecraft_id, sc_h.id);
    EXPECT_EQ(events.front().crossing, orbitsim::Crossing::Enter);
    EXPECT_TRUE(near_abs(events.front().t_event_s, 1.0, 1e-3));
}

TEST(GameSimulation, ProximityEnterExitWithHysteresis)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = true;
    cfg.events.time_tol_s = 1e-6;
    cfg.events.dist_tol_m = 1e-3;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    constexpr orbitsim::SpacecraftId center_id = 123u;
    constexpr orbitsim::SpacecraftId target_id = 456u;

    cfg.proximity.enable = true;
    cfg.proximity.center_spacecraft_id = center_id;
    cfg.proximity.enter_radius_m = 1.0e6;
    cfg.proximity.exit_radius_m = 1.2e6;

    orbitsim::GameSimulation sim(cfg);

    const auto center_h = sim.create_spacecraft_with_id(center_id, orbitsim::Spacecraft{
                                                                  .state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}),
                                                                  .dry_mass_kg = 1.0,
                                                                  .prop_mass_kg = 0.0,
                                                          });
    ASSERT_TRUE(center_h.valid());

    const auto target_h = sim.create_spacecraft_with_id(target_id, orbitsim::Spacecraft{
                                                                  .state = orbitsim::make_state({2.0e6, 0.0, 0.0}, {-500.0, 0.0, 0.0}),
                                                                  .dry_mass_kg = 1.0,
                                                                  .prop_mass_kg = 0.0,
                                                          });
    ASSERT_TRUE(target_h.valid());

    std::vector<orbitsim::Event> events;
    std::vector<orbitsim::Event> tmp;

    // Step in two chunks so each proximity crossing is bracketed by endpoints.
    sim.step(3000.0, &tmp); // Brackets enter at t=2000
    events.insert(events.end(), tmp.begin(), tmp.end());
    sim.step(4000.0, &tmp); // Brackets exit at t=6400
    events.insert(events.end(), tmp.begin(), tmp.end());

    std::vector<orbitsim::Event> prox;
    for (const auto &e: events)
    {
        if (e.type == orbitsim::EventType::Proximity)
        {
            prox.push_back(e);
        }
    }

    ASSERT_EQ(prox.size(), 2u);
    EXPECT_EQ(prox[0].crossing, orbitsim::Crossing::Enter);
    EXPECT_EQ(prox[0].spacecraft_id, target_id);
    EXPECT_EQ(prox[0].other_spacecraft_id, center_id);
    EXPECT_TRUE(near_abs(prox[0].t_event_s, 2000.0, 1e-2));

    EXPECT_EQ(prox[1].crossing, orbitsim::Crossing::Exit);
    EXPECT_EQ(prox[1].spacecraft_id, target_id);
    EXPECT_EQ(prox[1].other_spacecraft_id, center_id);
    EXPECT_TRUE(near_abs(prox[1].t_event_s, 6400.0, 1e-2));
}

TEST(GameSimulation, StepDetectsIntraStepImpactEnterAndExit)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = true;
    cfg.events.time_tol_s = 1e-6;
    cfg.events.dist_tol_m = 1e-6;
    cfg.events.crossing_scan_substeps = 32;
    cfg.events.max_event_splits_per_step = 16;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    const auto body_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 0.0,
            .radius_m = 1.0,
            .state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}),
    });
    ASSERT_TRUE(body_h.valid());

    const auto sc_h = sim.create_spacecraft(orbitsim::Spacecraft{
            .state = orbitsim::make_state({-2.0, 0.5, 0.0}, {4.0, 0.0, 0.0}),
            .dry_mass_kg = 1.0,
            .prop_mass_kg = 0.0,
    });
    ASSERT_TRUE(sc_h.valid());

    std::vector<orbitsim::Event> events;
    sim.step(1.0, &events);

    std::vector<orbitsim::Event> impacts;
    for (const auto &e : events)
    {
        if (e.type == orbitsim::EventType::Impact && e.body_id == body_h.id && e.spacecraft_id == sc_h.id)
        {
            impacts.push_back(e);
        }
    }

    ASSERT_GE(impacts.size(), 2u);
    EXPECT_EQ(impacts[0].crossing, orbitsim::Crossing::Enter);
    EXPECT_EQ(impacts[1].crossing, orbitsim::Crossing::Exit);
    EXPECT_TRUE(near_abs(impacts[0].t_event_s, 0.2834936, 5e-3));
    EXPECT_TRUE(near_abs(impacts[1].t_event_s, 0.7165064, 5e-3));
}

TEST(GameSimulation, StepDetectsIntraStepProximityEnterAndExit)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = true;
    cfg.events.time_tol_s = 1e-6;
    cfg.events.dist_tol_m = 1e-6;
    cfg.events.crossing_scan_substeps = 32;
    cfg.events.max_event_splits_per_step = 16;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    constexpr orbitsim::SpacecraftId center_id = 9001u;
    constexpr orbitsim::SpacecraftId target_id = 9002u;
    cfg.proximity.enable = true;
    cfg.proximity.center_spacecraft_id = center_id;
    cfg.proximity.enter_radius_m = 1.0;
    cfg.proximity.exit_radius_m = 1.0;

    orbitsim::GameSimulation sim(cfg);

    ASSERT_TRUE(sim.create_spacecraft_with_id(center_id, orbitsim::Spacecraft{
                                                              .state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}),
                                                              .dry_mass_kg = 1.0,
                                                              .prop_mass_kg = 0.0,
                                                      }).valid());

    ASSERT_TRUE(sim.create_spacecraft_with_id(target_id, orbitsim::Spacecraft{
                                                              .state = orbitsim::make_state({-2.0, 0.5, 0.0}, {4.0, 0.0, 0.0}),
                                                              .dry_mass_kg = 1.0,
                                                              .prop_mass_kg = 0.0,
                                                      }).valid());

    std::vector<orbitsim::Event> events;
    sim.step(1.0, &events);

    std::vector<orbitsim::Event> prox;
    for (const auto &e : events)
    {
        if (e.type == orbitsim::EventType::Proximity && e.spacecraft_id == target_id && e.other_spacecraft_id == center_id)
        {
            prox.push_back(e);
        }
    }

    ASSERT_GE(prox.size(), 2u);
    EXPECT_EQ(prox[0].crossing, orbitsim::Crossing::Enter);
    EXPECT_EQ(prox[1].crossing, orbitsim::Crossing::Exit);
    EXPECT_TRUE(near_abs(prox[0].t_event_s, 0.2834936, 5e-3));
    EXPECT_TRUE(near_abs(prox[1].t_event_s, 0.7165064, 5e-3));
}

TEST(GameSimulation, ZeroTimeToleranceStillAdvancesPastEvent)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = true;
    cfg.events.time_tol_s = 0.0;
    cfg.events.dist_tol_m = 1e-6;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);
    ASSERT_TRUE(sim.create_body(orbitsim::MassiveBody{
                                        .mass_kg = 0.0,
                                        .radius_m = 1.0,
                                        .state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}),
                                }).valid());
    ASSERT_TRUE(sim.create_spacecraft(orbitsim::Spacecraft{
                                              .state = orbitsim::make_state({2.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}),
                                              .dry_mass_kg = 1.0,
                                              .prop_mass_kg = 0.0,
                                      }).valid());

    std::vector<orbitsim::Event> events;
    sim.step(2.0, &events);

    EXPECT_TRUE(near_abs(sim.time_s(), 2.0, 1e-12));
    ASSERT_FALSE(events.empty());
}

TEST(GameSimulation, BackwardStepRewindsBurnAndPropellant)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    const auto body_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 0.0,
            .radius_m = 1.0,
            .state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}),
    });
    ASSERT_TRUE(body_h.valid());

    orbitsim::Spacecraft sc{};
    sc.state = orbitsim::make_state({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    sc.dry_mass_kg = 1.0;
    sc.prop_mass_kg = 0.2;
    sc.engines.push_back(orbitsim::Engine{.max_thrust_N = 10.0, .isp_s = 100.0, .min_throttle_0_1 = 0.0});
    const auto sc_h = sim.create_spacecraft(sc);
    ASSERT_TRUE(sc_h.valid());

    sim.maneuver_plan().segments.push_back(
        orbitsim::burn()
            .start(0.0)
            .duration(1.0)
            .radial_out()
            .primary(body_h)
            .engine(0)
            .spacecraft(sc_h)
            .build());

    sim.step(2.0);
    const orbitsim::Spacecraft *after_forward = sim.spacecraft_by_id(sc_h);
    ASSERT_NE(after_forward, nullptr);
    ASSERT_LT(after_forward->prop_mass_kg, 0.2);
    ASSERT_GT(after_forward->state.velocity_mps.x, 0.0);

    sim.step(-2.0);
    const orbitsim::Spacecraft *after_backward = sim.spacecraft_by_id(sc_h);
    ASSERT_NE(after_backward, nullptr);

    EXPECT_TRUE(near_abs(sim.time_s(), 0.0, 1e-9));
    EXPECT_TRUE(near_abs(after_backward->prop_mass_kg, 0.2, 1e-4));
    EXPECT_TRUE(near_vec_abs(after_backward->state.position_m, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-4));
    EXPECT_TRUE(near_vec_abs(after_backward->state.velocity_mps, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-4));
}

TEST(GameSimulation, BackwardStepReturnsReversedBoundaryCrossing)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = true;
    cfg.events.time_tol_s = 1e-6;
    cfg.events.dist_tol_m = 1e-6;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    const auto body_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 0.0,
            .radius_m = 1.0,
            .state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}),
    });
    ASSERT_TRUE(body_h.valid());

    const auto sc_h = sim.create_spacecraft(orbitsim::Spacecraft{
            .state = orbitsim::make_state({2.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}),
            .dry_mass_kg = 1.0,
            .prop_mass_kg = 0.0,
    });
    ASSERT_TRUE(sc_h.valid());

    std::vector<orbitsim::Event> forward_events;
    sim.step(2.0, &forward_events);
    ASSERT_FALSE(forward_events.empty());
    EXPECT_EQ(forward_events.front().crossing, orbitsim::Crossing::Enter);

    std::vector<orbitsim::Event> backward_events;
    sim.step(-2.0, &backward_events);
    ASSERT_FALSE(backward_events.empty());
    EXPECT_EQ(backward_events.front().type, orbitsim::EventType::Impact);
    EXPECT_EQ(backward_events.front().body_id, body_h.id);
    EXPECT_EQ(backward_events.front().crossing, orbitsim::Crossing::Exit);
    EXPECT_TRUE(near_abs(backward_events.front().t_event_s, 1.0, 1e-3));

    const orbitsim::Spacecraft *sc_back = sim.spacecraft_by_id(sc_h);
    ASSERT_NE(sc_back, nullptr);
    EXPECT_TRUE(near_abs(sim.time_s(), 0.0, 1e-9));
    EXPECT_TRUE(near_vec_abs(sc_back->state.position_m, orbitsim::Vec3{2.0, 0.0, 0.0}, 1e-6));
    EXPECT_TRUE(near_vec_abs(sc_back->state.velocity_mps, orbitsim::Vec3{-1.0, 0.0, 0.0}, 1e-6));
}
