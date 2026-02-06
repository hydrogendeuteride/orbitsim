#include <gtest/gtest.h>
#include "test_helpers.hpp"

TEST(Trajectories, PredictSpacecraftTrajectoryLinearMotion)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;
    cfg.spacecraft_integrator.adaptive = true;
    cfg.spacecraft_integrator.max_substeps = 256;

    orbitsim::GameSimulation sim(cfg);

    const orbitsim::GameSimulation::BodyHandle origin_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 0.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {100.0, 0.0, 0.0},
                            .velocity_mps = {0.5, 0.0, 0.0},
                    },
    });
    ASSERT_TRUE(origin_h.valid());

    orbitsim::Spacecraft sc{};
    sc.state.position_m = {110.0, 0.0, 0.0};
    sc.state.velocity_mps = {1.0, 0.0, 0.0};
    sc.dry_mass_kg = 1.0;
    sc.prop_mass_kg = 0.0;
    const orbitsim::GameSimulation::SpacecraftHandle sc_h = sim.create_spacecraft(sc);
    ASSERT_TRUE(sc_h.valid());

    orbitsim::TrajectoryOptions opt{};
    opt.duration_s = 10.0;
    opt.sample_dt_s = 1.0;
    opt.max_samples = 64;
    opt.include_start = true;
    opt.include_end = true;

    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris(sim, opt);
    const std::vector<orbitsim::TrajectorySample> traj_inertial =
            orbitsim::predict_spacecraft_trajectory(sim, eph, sc_h.id, opt);
    ASSERT_EQ(traj_inertial.size(), 11u);

    const std::vector<orbitsim::TrajectorySample> traj = orbitsim::trajectory_to_frame_spec(
            traj_inertial, eph, sim.massive_bodies(), orbitsim::TrajectoryFrameSpec::body_centered_inertial(origin_h));
    ASSERT_EQ(traj.size(), traj_inertial.size());

    // Relative position should drift at (v_sc - v_body) = 0.5 m/s, starting from 10m.
    EXPECT_TRUE(near_abs(traj.front().t_s, 0.0, 1e-12));
    EXPECT_TRUE(near_vec_abs(traj.front().position_m, orbitsim::Vec3{10.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(traj.front().velocity_mps, orbitsim::Vec3{0.5, 0.0, 0.0}, 1e-12));

    EXPECT_TRUE(near_abs(traj.back().t_s, 10.0, 1e-12));
    EXPECT_TRUE(near_vec_abs(traj.back().position_m, orbitsim::Vec3{15.0, 0.0, 0.0}, 1e-9));
    EXPECT_TRUE(near_vec_abs(traj.back().velocity_mps, orbitsim::Vec3{0.5, 0.0, 0.0}, 1e-9));
}

TEST(Trajectories, CelestialEphemerisLinearInterpolationAndClamp)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;

    orbitsim::GameSimulation sim(cfg);

    const orbitsim::GameSimulation::BodyHandle moving_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 1.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {0.0, 0.0, 0.0},
                            .velocity_mps = {1.0, 0.0, 0.0},
                    },
    });
    ASSERT_TRUE(moving_h.valid());

    orbitsim::TrajectoryOptions opt{};
    opt.duration_s = 4.0;
    opt.sample_dt_s = 1.0;
    opt.celestial_dt_s = 2.0;
    opt.max_samples = 64;

    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris(sim, opt);
    ASSERT_FALSE(eph.empty());
    ASSERT_EQ(eph.segments.size(), 2u);

    const orbitsim::State s_t1 = eph.body_state_at_by_id(moving_h, 1.0);
    EXPECT_TRUE(near_vec_abs(s_t1.position_m, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_t1.velocity_mps, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));

    const orbitsim::State s_t3 = eph.body_state_at_by_id(moving_h, 3.0);
    EXPECT_TRUE(near_vec_abs(s_t3.position_m, orbitsim::Vec3{3.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_t3.velocity_mps, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));

    const orbitsim::State s_before = eph.body_state_at_by_id(moving_h, -1.0);
    EXPECT_TRUE(near_vec_abs(s_before.position_m, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_before.velocity_mps, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));

    const orbitsim::State s_after = eph.body_state_at_by_id(moving_h, 10.0);
    EXPECT_TRUE(near_vec_abs(s_after.position_m, orbitsim::Vec3{4.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_after.velocity_mps, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));
}

TEST(Trajectories, BodyTrajectoryRelativeToSelfIsZero)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;

    orbitsim::GameSimulation sim(cfg);

    const orbitsim::GameSimulation::BodyHandle self_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 1.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {5.0, 0.0, 0.0},
                            .velocity_mps = {2.0, 0.0, 0.0},
                    },
    });
    ASSERT_TRUE(self_h.valid());

    orbitsim::TrajectoryOptions opt{};
    opt.duration_s = 5.0;
    opt.sample_dt_s = 1.0;
    opt.celestial_dt_s = 1.0;
    opt.max_samples = 64;

    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris(sim, opt);
    const std::vector<orbitsim::TrajectorySample> traj_inertial = orbitsim::predict_body_trajectory(sim, eph, self_h, opt);
    ASSERT_EQ(traj_inertial.size(), 6u);

    const std::vector<orbitsim::TrajectorySample> rel = orbitsim::trajectory_to_frame_spec(
            traj_inertial, eph, sim.massive_bodies(), orbitsim::TrajectoryFrameSpec::body_centered_inertial(self_h));
    ASSERT_EQ(rel.size(), traj_inertial.size());

    for (const auto &s: rel)
    {
        EXPECT_TRUE(near_vec_abs(s.position_m, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-12));
        EXPECT_TRUE(near_vec_abs(s.velocity_mps, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-12));
    }
}

TEST(Trajectories, TrajectoryToBodyCenteredInertialIsRelative)
{
    orbitsim::MassiveBody body{};
    body.state = orbitsim::make_state({100.0, 0.0, 0.0}, {1.0, 0.0, 0.0});

    const orbitsim::CelestialEphemeris eph{};

    const std::vector<orbitsim::TrajectorySample> inertial = {
            orbitsim::TrajectorySample{.t_s = 0.0, .position_m = {110.0, 0.0, 0.0}, .velocity_mps = {2.0, 0.0, 0.0}},
    };

    const std::vector<orbitsim::TrajectorySample> bci = orbitsim::trajectory_to_body_centered_inertial(inertial, eph, body);
    ASSERT_EQ(bci.size(), 1u);

    EXPECT_TRUE(near_vec_abs(bci[0].position_m, orbitsim::Vec3{10.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(bci[0].velocity_mps, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));
}

TEST(Trajectories, TrajectoryToBodyFixedMakesFixedPointStationaryAndRoundTrips)
{
    const double pi = std::acos(-1.0);
    constexpr double a_m = 1000.0;
    constexpr double omega_radps = 2.0;
    const double theta = 0.5 * pi;

    orbitsim::MassiveBody body{};
    body.state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, omega_radps, theta);

    const orbitsim::CelestialEphemeris eph{};

    const orbitsim::Vec3 r_i{a_m * std::cos(theta), a_m * std::sin(theta), 0.0};
    const orbitsim::Vec3 omega_i{0.0, 0.0, omega_radps};
    const orbitsim::Vec3 v_i = glm::cross(omega_i, r_i);

    const std::vector<orbitsim::TrajectorySample> inertial = {
            orbitsim::TrajectorySample{.t_s = 0.0, .position_m = r_i, .velocity_mps = v_i},
    };

    const std::vector<orbitsim::TrajectorySample> ecef = orbitsim::trajectory_to_body_fixed(inertial, eph, body);
    ASSERT_EQ(ecef.size(), 1u);

    EXPECT_TRUE(near_vec_abs(ecef[0].position_m, orbitsim::Vec3{a_m, 0.0, 0.0}, 1e-9));
    EXPECT_TRUE(near_vec_abs(ecef[0].velocity_mps, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-9));

    const std::optional<orbitsim::RotatingFrame> frame = orbitsim::make_body_fixed_frame(body);
    ASSERT_TRUE(frame.has_value());
    const orbitsim::TrajectorySample back = orbitsim::frame_sample_to_inertial(ecef[0], *frame);

    EXPECT_TRUE(near_vec_abs(back.position_m, inertial[0].position_m, 1e-9));
    EXPECT_TRUE(near_vec_abs(back.velocity_mps, inertial[0].velocity_mps, 1e-9));
}

TEST(Trajectories, PredictSpacecraftEventsFindsImpact)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;
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

    orbitsim::TrajectoryOptions opt{};
    opt.duration_s = 5.0;
    opt.sample_dt_s = 1.0;
    opt.celestial_dt_s = 1.0;
    opt.max_samples = 64;

    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris(sim, opt);
    const std::vector<orbitsim::Event> events = orbitsim::predict_spacecraft_events(sim, eph, sc_h.id, opt, cfg.events);

    ASSERT_FALSE(events.empty());
    EXPECT_EQ(events.front().type, orbitsim::EventType::Impact);
    EXPECT_EQ(events.front().crossing, orbitsim::Crossing::Enter);
    EXPECT_TRUE(near_abs(events.front().t_event_s, 1.0, 1e-3));
}

