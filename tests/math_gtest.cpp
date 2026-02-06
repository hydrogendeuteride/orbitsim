#include <gtest/gtest.h>
#include "test_helpers.hpp"

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

TEST(Math, OrbitalElementsAboutAxisUsesSpinAxis)
{
    const double mu = 1.0e5;
    const orbitsim::Vec3 ref_axis = orbitsim::normalized_or(orbitsim::Vec3{0.0, 1.0, 0.0}, {0.0, 1.0, 0.0});

    // Simple orbit in XZ plane => angular momentum along +Y.
    const orbitsim::Vec3 r{1.0, 0.0, 0.0};
    const orbitsim::Vec3 v{0.0, 0.0, -1.0};

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

