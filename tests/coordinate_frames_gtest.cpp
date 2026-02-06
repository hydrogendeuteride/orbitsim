#include <gtest/gtest.h>
#include "test_helpers.hpp"

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

