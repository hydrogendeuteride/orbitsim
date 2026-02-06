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
