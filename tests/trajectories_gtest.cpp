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

TEST(Trajectories, PredictSpacecraftTrajectorySegmentsLinearMotion)
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

    orbitsim::TrajectorySegmentOptions seg_opt{};
    seg_opt.duration_s = 10.0;
    seg_opt.max_segments = 16;
    seg_opt.lookup_dt_s = 1.0;
    seg_opt.include_start = true;
    seg_opt.include_end = true;

    const std::vector<orbitsim::TrajectorySegment> segments =
            orbitsim::predict_spacecraft_trajectory_segments(sim, sc_h.id, seg_opt);
    ASSERT_EQ(segments.size(), 10u);

    EXPECT_TRUE(near_abs(segments.front().t0_s, 0.0, 1e-12));
    EXPECT_TRUE(near_abs(segments.front().dt_s, 1.0, 1e-12));
    EXPECT_TRUE(near_vec_abs(segments.front().start.position_m, orbitsim::Vec3{110.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(segments.front().end.position_m, orbitsim::Vec3{111.0, 0.0, 0.0}, 1e-9));

    EXPECT_TRUE(near_abs(segments.back().t0_s, 9.0, 1e-12));
    EXPECT_TRUE(near_abs(segments.back().dt_s, 1.0, 1e-12));
    EXPECT_TRUE(near_vec_abs(segments.back().start.position_m, orbitsim::Vec3{119.0, 0.0, 0.0}, 1e-9));
    EXPECT_TRUE(near_vec_abs(segments.back().end.position_m, orbitsim::Vec3{120.0, 0.0, 0.0}, 1e-9));

    const std::vector<orbitsim::TrajectorySample> sampled =
            orbitsim::sample_trajectory_segments_uniform_dt(segments, 1.0, 64, true, true);
    ASSERT_EQ(sampled.size(), 11u);
    EXPECT_TRUE(near_abs(sampled.front().t_s, 0.0, 1e-12));
    EXPECT_TRUE(near_abs(sampled.back().t_s, 10.0, 1e-12));
    EXPECT_TRUE(near_vec_abs(sampled.front().position_m, orbitsim::Vec3{110.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(sampled.back().position_m, orbitsim::Vec3{120.0, 0.0, 0.0}, 1e-9));
}

TEST(Trajectories, AdaptiveSegmentsLinearMotionUsesSingleSegment)
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
            .state = orbitsim::make_state({100.0, 0.0, 0.0}, {0.5, 0.0, 0.0}),
    });
    ASSERT_TRUE(origin_h.valid());

    orbitsim::Spacecraft sc{};
    sc.state.position_m = {110.0, 0.0, 0.0};
    sc.state.velocity_mps = {1.0, 0.0, 0.0};
    sc.dry_mass_kg = 1.0;
    const auto sc_h = sim.create_spacecraft(sc);
    ASSERT_TRUE(sc_h.valid());

    orbitsim::AdaptiveEphemerisOptions eph_opt{};
    eph_opt.duration_s = 10.0;
    eph_opt.min_dt_s = 0.25;
    eph_opt.max_dt_s = 10.0;
    eph_opt.soft_max_segments = 8;
    eph_opt.hard_max_segments = 32;

    orbitsim::AdaptiveEphemerisDiagnostics eph_diag{};
    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris_adaptive(sim, eph_opt, &eph_diag);
    ASSERT_FALSE(eph.empty());
    EXPECT_EQ(eph_diag.accepted_segments, 1u);

    orbitsim::AdaptiveSegmentOptions seg_opt{};
    seg_opt.duration_s = 10.0;
    seg_opt.min_dt_s = 0.25;
    seg_opt.max_dt_s = 10.0;
    seg_opt.lookup_max_dt_s = 1.0;
    seg_opt.soft_max_segments = 8;
    seg_opt.hard_max_segments = 32;

    orbitsim::AdaptiveSegmentDiagnostics seg_diag{};
    const std::vector<orbitsim::TrajectorySegment> segments =
            orbitsim::predict_spacecraft_trajectory_segments_adaptive(sim, eph, sc_h.id, seg_opt, &seg_diag);
    ASSERT_EQ(segments.size(), 1u);
    EXPECT_EQ(seg_diag.accepted_segments, 1u);
    EXPECT_TRUE(near_abs(segments.front().t0_s, 0.0, 1e-12));
    EXPECT_TRUE(near_abs(segments.front().dt_s, 10.0, 1e-12));
    EXPECT_TRUE(near_vec_abs(segments.front().start.position_m, orbitsim::Vec3{110.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(segments.front().end.position_m, orbitsim::Vec3{120.0, 0.0, 0.0}, 1e-9));
}

TEST(Trajectories, AdaptiveSegmentsSplitAtImpulseBoundary)
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

    const auto sc_h = sim.create_spacecraft(orbitsim::Spacecraft{
            .state = orbitsim::make_state({10.0, 0.0, 0.0}, {0.0, 0.0, 0.0}),
            .dry_mass_kg = 1.0,
            .prop_mass_kg = 0.0,
    });
    ASSERT_TRUE(sc_h.valid());

    sim.maneuver_plan().impulses.push_back(orbitsim::ImpulseSegment{
            .t_s = 5.0,
            .primary_body_id = body_h.id,
            .dv_rtn_mps = {0.0, 1.0, 0.0},
            .spacecraft_id = sc_h.id,
    });
    orbitsim::sort_impulses_by_time(sim.maneuver_plan());

    orbitsim::AdaptiveEphemerisOptions eph_opt{};
    eph_opt.duration_s = 10.0;
    eph_opt.min_dt_s = 0.25;
    eph_opt.max_dt_s = 10.0;
    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris_adaptive(sim, eph_opt);

    orbitsim::AdaptiveSegmentOptions seg_opt{};
    seg_opt.duration_s = 10.0;
    seg_opt.min_dt_s = 0.25;
    seg_opt.max_dt_s = 10.0;
    seg_opt.lookup_max_dt_s = 1.0;

    orbitsim::AdaptiveSegmentDiagnostics diag{};
    const std::vector<orbitsim::TrajectorySegment> segments =
            orbitsim::predict_spacecraft_trajectory_segments_adaptive(sim, eph, sc_h.id, seg_opt, &diag);
    ASSERT_EQ(segments.size(), 2u);
    EXPECT_EQ(diag.forced_boundary_splits, 1u);

    EXPECT_TRUE(near_abs(segments[0].t0_s, 0.0, 1e-12));
    EXPECT_TRUE(near_abs(segments[0].dt_s, 5.0, 1e-12));
    EXPECT_EQ(segments[0].flags, 0u);
    EXPECT_TRUE(near_vec_abs(segments[0].end.position_m, orbitsim::Vec3{10.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(segments[0].end.velocity_mps, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-12));

    EXPECT_TRUE(near_abs(segments[1].t0_s, 5.0, 1e-12));
    EXPECT_TRUE((segments[1].flags & orbitsim::kTrajectorySegmentFlagImpulseBoundary) != 0u);
    EXPECT_TRUE(near_vec_abs(segments[1].start.position_m, orbitsim::Vec3{10.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(segments[1].start.velocity_mps, orbitsim::Vec3{0.0, 1.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(segments[1].end.position_m, orbitsim::Vec3{10.0, 5.0, 0.0}, 1e-9));
}

TEST(Trajectories, PredictSpacecraftTrajectorySegmentsStopOnImpact)
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

    orbitsim::TrajectoryOptions eph_opt{};
    eph_opt.duration_s = 5.0;
    eph_opt.sample_dt_s = 1.0;
    eph_opt.celestial_dt_s = 1.0;
    eph_opt.max_samples = 64;
    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris(sim, eph_opt);

    orbitsim::TrajectorySegmentOptions seg_opt{};
    seg_opt.duration_s = 5.0;
    seg_opt.max_segments = 64;
    seg_opt.lookup_dt_s = 1.0;
    seg_opt.stop_on_impact = true;

    const std::vector<orbitsim::TrajectorySegment> segments =
            orbitsim::predict_spacecraft_trajectory_segments(sim, eph, sc_h.id, seg_opt);
    ASSERT_FALSE(segments.empty());
    EXPECT_TRUE(near_abs(segments.front().t0_s, 0.0, 1e-12));
    EXPECT_TRUE(near_abs(segments.back().t0_s + segments.back().dt_s, 1.0, 1e-3));
    EXPECT_TRUE(near_vec_abs(segments.back().end.position_m, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-3));
}

TEST(Trajectories, AdaptiveSegmentsStopOnImpactMarksTerminalSegment)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;
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

    orbitsim::AdaptiveEphemerisOptions eph_opt{};
    eph_opt.duration_s = 5.0;
    eph_opt.min_dt_s = 0.25;
    eph_opt.max_dt_s = 5.0;
    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris_adaptive(sim, eph_opt);

    orbitsim::AdaptiveSegmentOptions seg_opt{};
    seg_opt.duration_s = 5.0;
    seg_opt.min_dt_s = 0.25;
    seg_opt.max_dt_s = 0.5;
    seg_opt.lookup_max_dt_s = 1.0;
    seg_opt.stop_on_impact = true;

    const std::vector<orbitsim::TrajectorySegment> segments =
            orbitsim::predict_spacecraft_trajectory_segments_adaptive(sim, eph, sc_h.id, seg_opt);
    ASSERT_FALSE(segments.empty());
    EXPECT_TRUE((segments.back().flags & orbitsim::kTrajectorySegmentFlagImpactTerminal) != 0u);
    EXPECT_TRUE(near_abs(segments.back().t0_s + segments.back().dt_s, 1.0, 1e-3));
    EXPECT_TRUE(near_vec_abs(segments.back().end.position_m, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-3));
}

TEST(Trajectories, AdaptiveEphemerisLinearInterpolationAndClamp)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;

    orbitsim::GameSimulation sim(cfg);

    const auto moving_h = sim.create_body(orbitsim::MassiveBody{
            .mass_kg = 1.0,
            .radius_m = 1.0,
            .state = orbitsim::make_state({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}),
    });
    ASSERT_TRUE(moving_h.valid());

    orbitsim::AdaptiveEphemerisOptions opt{};
    opt.duration_s = 4.0;
    opt.min_dt_s = 0.25;
    opt.max_dt_s = 4.0;

    orbitsim::AdaptiveEphemerisDiagnostics diag{};
    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris_adaptive(sim, opt, &diag);
    ASSERT_FALSE(eph.empty());
    ASSERT_EQ(eph.segments.size(), 1u);
    EXPECT_EQ(diag.accepted_segments, 1u);

    const orbitsim::State s_t2 = eph.body_state_at_by_id(moving_h.id, 2.0);
    EXPECT_TRUE(near_vec_abs(s_t2.position_m, orbitsim::Vec3{2.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_t2.velocity_mps, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));

    const orbitsim::State s_after = eph.body_state_at_by_id(moving_h.id, 10.0);
    EXPECT_TRUE(near_vec_abs(s_after.position_m, orbitsim::Vec3{4.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_after.velocity_mps, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));
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
