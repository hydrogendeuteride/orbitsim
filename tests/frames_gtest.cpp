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

