#include <gtest/gtest.h>
#include "test_helpers.hpp"

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

TEST(Frames, AdaptiveBodyFixedSegmentTransformResegmentsCurvedFrameMotion)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;

    orbitsim::GameSimulation sim(cfg);

    orbitsim::MassiveBody body{};
    body.mass_kg = 1.0;
    body.radius_m = 1.0;
    body.state = orbitsim::make_state({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 1.0, 0.0);
    const auto body_h = sim.create_body(body);
    ASSERT_TRUE(body_h.valid());

    orbitsim::AdaptiveEphemerisOptions eph_opt{};
    eph_opt.duration_s = 2.0;
    eph_opt.min_dt_s = 0.125;
    eph_opt.max_dt_s = 2.0;
    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris_adaptive(sim, eph_opt);
    ASSERT_FALSE(eph.empty());

    const orbitsim::TrajectorySegment inertial_segment{
            .t0_s = 0.0,
            .dt_s = 2.0,
            .start = orbitsim::make_state({1000.0, 0.0, 0.0}, {0.0, 0.0, 0.0}),
            .end = orbitsim::make_state({1000.0, 0.0, 0.0}, {0.0, 0.0, 0.0}),
            .flags = 0u,
    };

    orbitsim::FrameSegmentTransformOptions frame_opt{};
    frame_opt.min_dt_s = 0.125;
    frame_opt.max_dt_s = 2.0;
    frame_opt.soft_max_segments = 16;
    frame_opt.hard_max_segments = 64;
    frame_opt.tolerance.pos_near_m = 1.0e-3;
    frame_opt.tolerance.pos_far_m = 1.0e-3;
    frame_opt.tolerance.vel_near_mps = 1.0e-3;
    frame_opt.tolerance.vel_far_mps = 1.0e-3;

    orbitsim::FrameSegmentTransformDiagnostics diag{};
    const std::vector<orbitsim::TrajectorySegment> transformed =
            orbitsim::transform_trajectory_segments_to_frame_spec(
                    {inertial_segment},
                    eph,
                    sim.massive_bodies(),
                    orbitsim::TrajectoryFrameSpec::body_fixed(body_h.id),
                    frame_opt,
                    nullptr,
                    &diag);
    ASSERT_GT(transformed.size(), 1u);
    EXPECT_GT(diag.frame_resegmentation_count, 0u);

    bool saw_resegmented_flag = false;
    for (const auto &segment : transformed)
    {
        saw_resegmented_flag = saw_resegmented_flag ||
                               ((segment.flags & orbitsim::kTrajectorySegmentFlagFrameResegmented) != 0u);
    }
    EXPECT_TRUE(saw_resegmented_flag);

    const std::vector<orbitsim::TrajectorySample> exact_inertial = {
            {.t_s = 0.0, .position_m = {1000.0, 0.0, 0.0}, .velocity_mps = {0.0, 0.0, 0.0}},
            {.t_s = 1.0, .position_m = {1000.0, 0.0, 0.0}, .velocity_mps = {0.0, 0.0, 0.0}},
            {.t_s = 2.0, .position_m = {1000.0, 0.0, 0.0}, .velocity_mps = {0.0, 0.0, 0.0}},
    };
    const std::vector<orbitsim::TrajectorySample> exact_body_fixed = orbitsim::trajectory_to_frame_spec(
            exact_inertial,
            eph,
            sim.massive_bodies(),
            orbitsim::TrajectoryFrameSpec::body_fixed(body_h.id));
    ASSERT_EQ(exact_body_fixed.size(), 3u);

    const std::vector<orbitsim::TrajectorySample> sampled =
            orbitsim::sample_trajectory_segments_uniform_dt(transformed, 1.0, 8, true, true);
    ASSERT_EQ(sampled.size(), 3u);

    for (std::size_t i = 0; i < sampled.size(); ++i)
    {
        EXPECT_TRUE(near_abs(sampled[i].t_s, exact_body_fixed[i].t_s, 1e-12));
        EXPECT_TRUE(near_vec_abs(sampled[i].position_m, exact_body_fixed[i].position_m, 5.0e-3));
        EXPECT_TRUE(near_vec_abs(sampled[i].velocity_mps, exact_body_fixed[i].velocity_mps, 5.0e-3));
    }
}
