#include <orbitsim/orbitsim.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <optional>
#include <vector>

namespace
{

    inline bool near_abs(const double a, const double b, const double abs_tol) { return std::abs(a - b) <= abs_tol; }

    inline bool near_rel(const double a, const double b, const double rel_tol, const double abs_floor = 0.0)
    {
        const double scale = std::max(abs_floor, std::max(std::abs(a), std::abs(b)));
        return std::abs(a - b) <= rel_tol * scale;
    }

    inline bool near_vec_abs(const orbitsim::Vec3 &a, const orbitsim::Vec3 &b, const double abs_tol)
    {
        return near_abs(a.x, b.x, abs_tol) && near_abs(a.y, b.y, abs_tol) && near_abs(a.z, b.z, abs_tol);
    }

} // namespace

TEST(Math, HermiteVelocityEndpoints)
{
    const orbitsim::Vec3 p0{0.0, 0.0, 0.0};
    const orbitsim::Vec3 v0{1.0, 2.0, 3.0};
    const orbitsim::Vec3 p1{10.0, -3.0, 5.0};
    const orbitsim::Vec3 v1{-4.0, 0.5, 1.25};
    const double dt_s = 2.0;

    const orbitsim::Vec3 dv0 = orbitsim::hermite_velocity_mps(p0, v0, p1, v1, dt_s, 0.0);
    const orbitsim::Vec3 dv1 = orbitsim::hermite_velocity_mps(p0, v0, p1, v1, dt_s, 1.0);
    EXPECT_TRUE(near_vec_abs(dv0, v0, 1e-12));
    EXPECT_TRUE(near_vec_abs(dv1, v1, 1e-12));
}

TEST(Math, RTNFrameOrthonormal)
{
    const orbitsim::Vec3 r{7'000'000.0, 0.0, 0.0};
    const orbitsim::Vec3 v{0.0, 7'500.0, 0.0};
    const orbitsim::RtnFrame f = orbitsim::compute_rtn_frame(r, v);
    EXPECT_TRUE(near_rel(glm::length(f.R), 1.0, 1e-12));
    EXPECT_TRUE(near_rel(glm::length(f.T), 1.0, 1e-12));
    EXPECT_TRUE(near_rel(glm::length(f.N), 1.0, 1e-12));
    EXPECT_LT(std::abs(glm::dot(f.R, f.T)), 1e-12);
    EXPECT_LT(std::abs(glm::dot(f.R, f.N)), 1e-12);
    EXPECT_LT(std::abs(glm::dot(f.T, f.N)), 1e-12);
}

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

    std::cout << "====== Test completed ======" << std::endl;
}

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
    opt.frame = orbitsim::TrajectoryFrameSpec::body_centered_inertial(origin_h);

    const std::vector<orbitsim::TrajectorySample> traj = orbitsim::predict_spacecraft_trajectory(sim, sc_h, opt);
    ASSERT_EQ(traj.size(), 11u);

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
    opt.frame = orbitsim::TrajectoryFrameSpec::body_centered_inertial(self_h);

    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris(sim, opt);
    const std::vector<orbitsim::TrajectorySample> traj = orbitsim::predict_body_trajectory(sim, eph, self_h, opt);
    ASSERT_EQ(traj.size(), 6u);
    for (const auto &s: traj)
    {
        EXPECT_TRUE(near_vec_abs(s.position_m, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-12));
        EXPECT_TRUE(near_vec_abs(s.velocity_mps, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-12));
    }
}

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

TEST(Frames, SynodicFrameOrthonormalAndBodyFixed)
{
    // Two-body circular orbit in barycentric inertial coordinates.
    constexpr double m_a_kg = 10.0;
    constexpr double m_b_kg = 1.0;
    constexpr double separation_m = 1000.0;

    const auto states = orbitsim::two_body_circular_barycentric(m_a_kg, m_b_kg, separation_m);
    const std::optional<orbitsim::SynodicFrame> frame =
            orbitsim::make_synodic_frame(states.state_a, m_a_kg, states.state_b, m_b_kg);
    ASSERT_TRUE(frame.has_value());
    ASSERT_TRUE(frame->valid());

    EXPECT_TRUE(near_rel(glm::length(frame->ex_i), 1.0, 1e-12));
    EXPECT_TRUE(near_rel(glm::length(frame->ey_i), 1.0, 1e-12));
    EXPECT_TRUE(near_rel(glm::length(frame->ez_i), 1.0, 1e-12));
    EXPECT_LT(std::abs(glm::dot(frame->ex_i, frame->ey_i)), 1e-12);
    EXPECT_LT(std::abs(glm::dot(frame->ex_i, frame->ez_i)), 1e-12);
    EXPECT_LT(std::abs(glm::dot(frame->ey_i, frame->ez_i)), 1e-12);
    EXPECT_GT(glm::dot(glm::cross(frame->ex_i, frame->ey_i), frame->ez_i), 0.999999999999);

    const double mu = m_b_kg / (m_a_kg + m_b_kg);

    const orbitsim::State a_rot = orbitsim::inertial_state_to_frame(states.state_a, *frame);
    const orbitsim::State b_rot = orbitsim::inertial_state_to_frame(states.state_b, *frame);

    EXPECT_TRUE(near_vec_abs(a_rot.position_m, orbitsim::Vec3{-mu * separation_m, 0.0, 0.0}, 1e-10));
    EXPECT_TRUE(near_vec_abs(b_rot.position_m, orbitsim::Vec3{(1.0 - mu) * separation_m, 0.0, 0.0}, 1e-10));

    // In the synodic frame, both primaries should be (nearly) stationary.
    EXPECT_TRUE(near_vec_abs(a_rot.velocity_mps, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-9));
    EXPECT_TRUE(near_vec_abs(b_rot.velocity_mps, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-9));

    // Round-trip transform sanity check.
    const orbitsim::State a_back = orbitsim::frame_state_to_inertial(a_rot, *frame);
    EXPECT_TRUE(near_vec_abs(a_back.position_m, states.state_a.position_m, 1e-9));
    EXPECT_TRUE(near_vec_abs(a_back.velocity_mps, states.state_a.velocity_mps, 1e-9));
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
    sim.step(5.0, &events);

    ASSERT_FALSE(events.empty());
    EXPECT_EQ(events.front().type, orbitsim::EventType::Impact);
    EXPECT_EQ(events.front().body_id, body_h.id);
    EXPECT_EQ(events.front().spacecraft_id, sc_h.id);
    EXPECT_EQ(events.front().crossing, orbitsim::Crossing::Enter);
    EXPECT_TRUE(near_abs(events.front().t_event_s, 1.0, 1e-3));
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

    sim.maneuver_plan().impulses.push_back(orbitsim::impulse().time(1.0).prograde(5.0).primary(body_h).spacecraft(sc_h));

    sim.step(2.0);

    const orbitsim::Spacecraft *out = sim.spacecraft_by_id(sc_h.id);
    ASSERT_NE(out, nullptr);
    EXPECT_TRUE(near_abs(out->state.velocity_mps.y, 6.0, 1e-9));
}

TEST(Geodesy, RoundTripGeodeticSpherical)
{
    orbitsim::MassiveBody body{};
    body.radius_m = 1'000.0;
    body.state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 1.0, 0.0);

    const std::optional<orbitsim::RotatingFrame> bf = orbitsim::make_body_fixed_frame(body);
    ASSERT_TRUE(bf.has_value());

    const orbitsim::GeodeticCoord c0{.latitude_rad = 0.3, .longitude_rad = -1.0, .altitude_m = 123.0};
    const orbitsim::Vec3 p_i = orbitsim::inertial_position_from_geodetic(*bf, c0, body.radius_m);
    const std::optional<orbitsim::GeodeticCoord> c1 = orbitsim::geodetic_from_inertial(*bf, p_i, body.radius_m);
    ASSERT_TRUE(c1.has_value());

    EXPECT_TRUE(near_abs(c1->latitude_rad, c0.latitude_rad, 1e-12));
    EXPECT_TRUE(near_abs(c1->longitude_rad, c0.longitude_rad, 1e-12));
    EXPECT_TRUE(near_abs(c1->altitude_m, c0.altitude_m, 1e-9));
}

TEST(Math, OrbitalElementsRoundTripToState)
{
    const double mu = 3.986004418e14; // Earth mu [m^3/s^2]

    orbitsim::OrbitalElements el{};
    el.semi_major_axis_m = 7'000'000.0;
    el.eccentricity = 0.1;
    el.inclination_rad = 0.5;
    el.raan_rad = 1.0;
    el.arg_periapsis_rad = 0.3;
    el.true_anomaly_rad = 0.9;

    const orbitsim::State s0 = orbitsim::relative_state_from_orbital_elements(mu, el);
    ASSERT_TRUE(std::isfinite(s0.position_m.x));
    ASSERT_TRUE(std::isfinite(s0.velocity_mps.x));

    const orbitsim::OrbitalElements el2 = orbitsim::orbital_elements_from_relative_state(mu, s0.position_m, s0.velocity_mps);
    const orbitsim::State s1 = orbitsim::relative_state_from_orbital_elements(mu, el2);

    EXPECT_TRUE(near_vec_abs(s1.position_m, s0.position_m, 1e-3));
    EXPECT_TRUE(near_vec_abs(s1.velocity_mps, s0.velocity_mps, 1e-6));
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
    sim.step(7000.0, &events);

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

TEST(Math, OrbitalElementsAboutAxisUsesSpinAxis)
{
    const double mu = 1.0e5;
    const orbitsim::Vec3 ref_axis = orbitsim::normalized_or(orbitsim::Vec3{0.0, 1.0, 0.0}, {0.0, 1.0, 0.0});

    // Simple orbit in XZ plane => angular momentum along +Y.
    const orbitsim::Vec3 r{1.0, 0.0, 0.0};
    const orbitsim::Vec3 v{0.0, 0.0, 1.0};

    const orbitsim::OrbitalElements el = orbitsim::orbital_elements_from_relative_state_about_axis(mu, r, v, ref_axis);
    EXPECT_TRUE(near_abs(el.inclination_rad, 0.0, 1e-9));
    EXPECT_TRUE(std::isfinite(el.mean_anomaly_rad) || std::isnan(el.mean_anomaly_rad));
}

TEST(Math, ApsidesFromElementsMatchRadii)
{
    const double mu = 3.986004418e14; // Earth mu [m^3/s^2]

    orbitsim::OrbitalElements el{};
    el.semi_major_axis_m = 7'000'000.0;
    el.eccentricity = 0.2;
    el.inclination_rad = 0.1;
    el.raan_rad = 0.4;
    el.arg_periapsis_rad = 0.7;
    el.true_anomaly_rad = 1.2;

    const orbitsim::State s = orbitsim::relative_state_from_orbital_elements(mu, el);
    const orbitsim::OrbitApsides aps = orbitsim::apsides_from_relative_state(mu, s.position_m, s.velocity_mps);
    ASSERT_TRUE(aps.valid);

    const double rp_expected = el.semi_major_axis_m * (1.0 - el.eccentricity);
    const double ra_expected = el.semi_major_axis_m * (1.0 + el.eccentricity);

    EXPECT_TRUE(near_rel(glm::length(aps.periapsis_rel_m), rp_expected, 1e-12, 1.0));
    ASSERT_TRUE(aps.has_apoapsis);
    EXPECT_TRUE(near_rel(glm::length(aps.apoapsis_rel_m), ra_expected, 1e-12, 1.0));

    // Periapsis and apoapsis should be opposite directions in the orbit plane.
    const double d = glm::dot(orbitsim::normalized_or(aps.periapsis_rel_m, {1.0, 0.0, 0.0}),
                              orbitsim::normalized_or(aps.apoapsis_rel_m, {-1.0, 0.0, 0.0}));
    EXPECT_LT(d, -0.999999);
}

TEST(Frames, Cr3bpLagrangePointMarkersConsistency)
{
    constexpr double m_a_kg = 4.0;
    constexpr double m_b_kg = 1.0;
    constexpr double separation_m = 2.0;
    const auto states = orbitsim::two_body_circular_barycentric(m_a_kg, m_b_kg, separation_m);
    const std::optional<orbitsim::SynodicFrame> frame =
            orbitsim::make_synodic_frame(states.state_a, m_a_kg, states.state_b, m_b_kg);
    ASSERT_TRUE(frame.has_value());

    const std::optional<orbitsim::Cr3bpLagrangePoints> pts = orbitsim::cr3bp_lagrange_points_m(*frame);
    ASSERT_TRUE(pts.has_value());

    const double mu = m_b_kg / (m_a_kg + m_b_kg);
    EXPECT_TRUE(near_vec_abs(pts->primary_m, orbitsim::Vec3{-mu * separation_m, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(pts->secondary_m, orbitsim::Vec3{(1.0 - mu) * separation_m, 0.0, 0.0}, 1e-12));

    // Triangular points are equidistant from both primaries by construction.
    const double dL4A = glm::length(pts->L4_m - pts->primary_m);
    const double dL4B = glm::length(pts->L4_m - pts->secondary_m);
    EXPECT_TRUE(near_rel(dL4A, separation_m, 1e-12));
    EXPECT_TRUE(near_rel(dL4B, separation_m, 1e-12));

    // Collinear points satisfy dΩ/dx = 0 in nondimensional CR3BP coordinates.
    const auto f = [&](const double x_nd) -> double {
        const double x1 = x_nd + mu;
        const double x2 = x_nd - (1.0 - mu);
        return x_nd - (1.0 - mu) * (x1 / std::pow(std::abs(x1), 3.0)) - mu * (x2 / std::pow(std::abs(x2), 3.0));
    };
    EXPECT_NEAR(f(pts->L1_m.x / separation_m), 0.0, 1e-12);
    EXPECT_NEAR(f(pts->L2_m.x / separation_m), 0.0, 1e-12);
    EXPECT_NEAR(f(pts->L3_m.x / separation_m), 0.0, 1e-12);
}

TEST(CoordinateFrames, BodyCenteredInertialRoundTrip)
{
    const orbitsim::State body_state = orbitsim::make_state({100.0, 200.0, 300.0}, {1.0, 2.0, 3.0});
    const orbitsim::RotatingFrame bci = orbitsim::make_body_centered_inertial_frame(body_state);

    const orbitsim::State s_in = orbitsim::make_state({110.0, 250.0, 310.0}, {5.0, -2.0, 1.0});
    const orbitsim::State s_frame = orbitsim::inertial_state_to_frame(s_in, bci);
    const orbitsim::State s_back = orbitsim::frame_state_to_inertial(s_frame, bci);

    EXPECT_TRUE(near_vec_abs(s_back.position_m, s_in.position_m, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_back.velocity_mps, s_in.velocity_mps, 1e-12));
}

TEST(CoordinateFrames, BodyFixedQuarterTurnAboutZ)
{
    orbitsim::State body_state{};
    body_state.position_m = {0.0, 0.0, 0.0};
    body_state.velocity_mps = {0.0, 0.0, 0.0};
    body_state.spin.axis = {0.0, 0.0, 1.0};
    body_state.spin.angle_rad = 0.5 * std::acos(-1.0);
    body_state.spin.rate_rad_per_s = 1.0;

    const std::optional<orbitsim::RotatingFrame> ecef = orbitsim::make_body_fixed_frame(body_state);
    ASSERT_TRUE(ecef.has_value());

    // With a +90deg rotation about +Z, inertial +Y aligns with body-fixed +X.
    const orbitsim::Vec3 v_in{0.0, 1.0, 0.0};
    const orbitsim::Vec3 v_frame = orbitsim::inertial_vector_to_frame(*ecef, v_in);
    EXPECT_TRUE(near_vec_abs(v_frame, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));

    // And inertial +X maps to body-fixed -Y.
    const orbitsim::Vec3 x_in{1.0, 0.0, 0.0};
    const orbitsim::Vec3 x_frame = orbitsim::inertial_vector_to_frame(*ecef, x_in);
    EXPECT_TRUE(near_vec_abs(x_frame, orbitsim::Vec3{0.0, -1.0, 0.0}, 1e-12));
}

TEST(CoordinateFrames, LvlhAxesAndOmega)
{
    const orbitsim::State primary = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    const orbitsim::State sc = orbitsim::make_state({7'000'000.0, 0.0, 0.0}, {0.0, 7'500.0, 0.0});

    const std::optional<orbitsim::RotatingFrame> lvlh = orbitsim::make_lvlh_frame(primary, sc);
    ASSERT_TRUE(lvlh.has_value());

    EXPECT_TRUE(near_vec_abs(lvlh->ex_i, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(lvlh->ey_i, orbitsim::Vec3{0.0, 1.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(lvlh->ez_i, orbitsim::Vec3{0.0, 0.0, 1.0}, 1e-12));

    const double expected_omega = 7'500.0 / 7'000'000.0;
    EXPECT_TRUE(near_vec_abs(lvlh->omega_inertial_radps, orbitsim::Vec3{0.0, 0.0, expected_omega}, 1e-12));
}

TEST(Lambert, ZeroRevSolutionPropagatesToTarget)
{
    const double mu = 1.0;

    const orbitsim::Vec3 r1{1.0, 0.0, 0.0};
    const double ang = std::acos(-1.0) / 3.0; // 60 deg
    const orbitsim::Vec3 r2{std::cos(ang), std::sin(ang), 0.0};
    const double dt_s = 1.0;

    orbitsim::LambertOptions opt{};
    opt.max_revolutions = 0;
    opt.short_path = true;
    opt.prograde = true;

    const std::vector<orbitsim::LambertSolution> sols = orbitsim::solve_lambert_universal(mu, r1, r2, dt_s, opt);
    ASSERT_FALSE(sols.empty());

    // For this configuration we expect a single 0-rev solution.
    ASSERT_EQ(sols.size(), 1u);

    const auto &sol = sols.front();
    const orbitsim::KeplerStepResult step = orbitsim::propagate_kepler_universal(mu, r1, sol.v1_mps, dt_s);
    ASSERT_TRUE(step.converged);

    EXPECT_TRUE(near_vec_abs(step.position_m, r2, 1e-6));
    EXPECT_TRUE(near_vec_abs(step.velocity_mps, sol.v2_mps, 1e-6));
}

TEST(Lambert, MultiRevSolutionsPropagateToTarget)
{
    const double mu = 1.0;

    const orbitsim::Vec3 r1{1.0, 0.0, 0.0};
    const double ang = std::acos(-1.0) / 3.0; // 60 deg
    const orbitsim::Vec3 r2{std::cos(ang), std::sin(ang), 0.0};

    // Longer time-of-flight allows multi-revolution solutions.
    const double dt_s = 10.0;

    orbitsim::LambertOptions opt{};
    opt.max_revolutions = 2;
    opt.short_path = true;
    opt.prograde = true;

    const std::vector<orbitsim::LambertSolution> sols = orbitsim::solve_lambert_universal(mu, r1, r2, dt_s, opt);
    ASSERT_GE(sols.size(), 2u);

    for (const auto &sol: sols)
    {
        const orbitsim::KeplerStepResult step = orbitsim::propagate_kepler_universal(mu, r1, sol.v1_mps, dt_s);
        ASSERT_TRUE(step.converged);
        EXPECT_TRUE(near_vec_abs(step.position_m, r2, 1e-5));
    }
}
