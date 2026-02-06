#include <gtest/gtest.h>
#include "test_helpers.hpp"

TEST(Nodes, PredictEquatorialNodesLinearMotion)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;
    cfg.events.time_tol_s = 1e-6;
    cfg.events.dist_tol_m = 1e-6;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    orbitsim::MassiveBody body{};
    body.radius_m = 1.0;
    body.state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 0.0, 0.0);
    const auto body_h = sim.create_body(body);
    ASSERT_TRUE(body_h.valid());

    const auto sc_h = sim.create_spacecraft(orbitsim::Spacecraft{
            .state = orbitsim::make_state({0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}),
            .dry_mass_kg = 1.0,
            .prop_mass_kg = 0.0,
    });
    ASSERT_TRUE(sc_h.valid());

    orbitsim::TrajectoryOptions opt{};
    opt.duration_s = 3.0;
    opt.sample_dt_s = 1.0;
    opt.celestial_dt_s = 1.0;
    opt.max_samples = 64;

    const auto eph = orbitsim::build_celestial_ephemeris(sim, opt);
    const std::vector<orbitsim::NodeEvent> nodes = orbitsim::predict_equatorial_nodes(sim, eph, sc_h.id, body_h.id, opt, cfg.events);
    ASSERT_FALSE(nodes.empty());
    EXPECT_TRUE(near_abs(nodes.front().t_event_s, 1.0, 1e-3));
    EXPECT_EQ(nodes.front().crossing, orbitsim::NodeCrossing::Ascending);
}

TEST(Nodes, PredictTargetPlaneNodesLinearMotion)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;
    cfg.events.time_tol_s = 1e-6;
    cfg.events.dist_tol_m = 1e-6;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    orbitsim::MassiveBody body{};
    body.radius_m = 1.0;
    body.state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 0.0, 0.0);
    const auto body_h = sim.create_body(body);
    ASSERT_TRUE(body_h.valid());

    // Target defines the reference plane: choose motion in XY plane => normal roughly +Z.
    const auto target_h = sim.create_spacecraft(orbitsim::Spacecraft{
            .state = orbitsim::make_state({0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}),
            .dry_mass_kg = 1.0,
            .prop_mass_kg = 0.0,
    });
    ASSERT_TRUE(target_h.valid());

    // Spacecraft crosses the target plane at t=1 (z goes from +1 to 0).
    const auto sc_h = sim.create_spacecraft(orbitsim::Spacecraft{
            .state = orbitsim::make_state({0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}),
            .dry_mass_kg = 1.0,
            .prop_mass_kg = 0.0,
    });
    ASSERT_TRUE(sc_h.valid());

    orbitsim::TrajectoryOptions opt{};
    opt.duration_s = 3.0;
    opt.sample_dt_s = 1.0;
    opt.celestial_dt_s = 1.0;
    opt.max_samples = 64;

    const auto eph = orbitsim::build_celestial_ephemeris(sim, opt);
    const std::vector<orbitsim::NodeEvent> nodes =
            orbitsim::predict_target_plane_nodes(sim, eph, sc_h.id, target_h.id, body_h.id, opt, cfg.events);
    ASSERT_FALSE(nodes.empty());
    EXPECT_TRUE(near_abs(nodes.front().t_event_s, 1.0, 1e-3));
}

TEST(Nodes, PredictTargetPlaneNodesUsesTimeVaryingTargetPlane)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;
    cfg.events.time_tol_s = 1e-6;
    cfg.events.dist_tol_m = 1e-6;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    orbitsim::MassiveBody body{};
    body.radius_m = 1.0;
    body.state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 0.0, 0.0);
    const auto body_h = sim.create_body(body);
    ASSERT_TRUE(body_h.valid());

    const auto target_h = sim.create_spacecraft(orbitsim::Spacecraft{
            .state = orbitsim::make_state({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}),
            .dry_mass_kg = 1.0,
            .prop_mass_kg = 0.0,
    });
    ASSERT_TRUE(target_h.valid());

    // Rotate the target's orbital plane at t=1 by applying +N impulse.
    sim.maneuver_plan().impulses.push_back(
            orbitsim::impulse()
                    .time(1.0)
                    .normal(1.0)
                    .primary(body_h)
                    .spacecraft(target_h)
                    .build());

    // This spacecraft never crosses z=0, so fixed initial target plane (+Z) would miss it.
    // After the target impulse, the updated plane is tilted and crossed near t=2.
    const auto sc_h = sim.create_spacecraft(orbitsim::Spacecraft{
            .state = orbitsim::make_state({0.0, -3.0, 1.0}, {0.0, 2.0, 0.0}),
            .dry_mass_kg = 1.0,
            .prop_mass_kg = 0.0,
    });
    ASSERT_TRUE(sc_h.valid());

    orbitsim::TrajectoryOptions opt{};
    opt.duration_s = 4.0;
    opt.sample_dt_s = 0.25;
    opt.celestial_dt_s = 0.25;
    opt.spacecraft_lookup_dt_s = 0.25;
    opt.max_samples = 128;

    const auto eph = orbitsim::build_celestial_ephemeris(sim, opt);
    const std::vector<orbitsim::NodeEvent> nodes =
            orbitsim::predict_target_plane_nodes(sim, eph, sc_h.id, target_h.id, body_h.id, opt, cfg.events);
    ASSERT_FALSE(nodes.empty());
    EXPECT_TRUE(near_abs(nodes.front().t_event_s, 2.0, 5e-2));
}
