#include <gtest/gtest.h>
#include "test_helpers.hpp"

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

TEST(Lambert, ProgradeRetrogradeUsesReferenceNormalIn3D)
{
    const double mu = 1.0;

    const orbitsim::Vec3 r1{1.0, 0.0, 0.2};
    const orbitsim::Vec3 r2{0.2, 1.0, 0.4};
    const double dt_s = 20.0;
    const orbitsim::Vec3 ref_n = orbitsim::normalized_or(glm::cross(r1, r2), orbitsim::Vec3{0.0, 0.0, 1.0});

    orbitsim::LambertOptions prograde{};
    prograde.max_revolutions = 0;
    prograde.short_path = true;
    prograde.prograde = true;
    prograde.reference_normal_i = ref_n;

    orbitsim::LambertOptions retrograde = prograde;
    retrograde.prograde = false;

    const std::vector<orbitsim::LambertSolution> sols_pro =
            orbitsim::solve_lambert_universal(mu, r1, r2, dt_s, prograde);
    const std::vector<orbitsim::LambertSolution> sols_ret =
            orbitsim::solve_lambert_universal(mu, r1, r2, dt_s, retrograde);
    ASSERT_FALSE(sols_pro.empty());
    ASSERT_FALSE(sols_ret.empty());

    const orbitsim::Vec3 h_pro = glm::cross(r1, sols_pro.front().v1_mps);
    const orbitsim::Vec3 h_ret = glm::cross(r1, sols_ret.front().v1_mps);
    EXPECT_GT(glm::dot(h_pro, ref_n), 1e-8);
    EXPECT_LT(glm::dot(h_ret, ref_n), -1e-8);

    const orbitsim::KeplerStepResult step_pro =
            orbitsim::propagate_kepler_universal(mu, r1, sols_pro.front().v1_mps, dt_s);
    const orbitsim::KeplerStepResult step_ret =
            orbitsim::propagate_kepler_universal(mu, r1, sols_ret.front().v1_mps, dt_s);
    ASSERT_TRUE(step_pro.converged);
    ASSERT_TRUE(step_ret.converged);
    EXPECT_TRUE(near_vec_abs(step_pro.position_m, r2, 1e-5));
    EXPECT_TRUE(near_vec_abs(step_ret.position_m, r2, 1e-5));
}
