#include <gtest/gtest.h>
#include "test_helpers.hpp"

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

