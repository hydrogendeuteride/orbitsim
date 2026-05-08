#include <gtest/gtest.h>
#include "test_helpers.hpp"

#include <cmath>
#include <vector>

namespace
{
    constexpr double kEarthMu = 3.986004418e14;

    bool finite_state(const orbitsim::State &state)
    {
        return std::isfinite(state.position_m.x) && std::isfinite(state.position_m.y) &&
               std::isfinite(state.position_m.z) && std::isfinite(state.velocity_mps.x) &&
               std::isfinite(state.velocity_mps.y) && std::isfinite(state.velocity_mps.z);
    }
} // namespace

TEST(Kepler, CircularOrbitRoundTripsOnePeriod)
{
    constexpr double r_m = 7'000'000.0;
    const double v_mps = std::sqrt(kEarthMu / r_m);
    const double period_s = 2.0 * std::acos(-1.0) * std::sqrt((r_m * r_m * r_m) / kEarthMu);
    const orbitsim::State s0 = orbitsim::make_state({r_m, 0.0, 0.0}, {0.0, v_mps, 0.0});

    const orbitsim::KeplerPropagationResult step =
            orbitsim::propagate_kepler_universal_safe(kEarthMu, s0, period_s);

    ASSERT_TRUE(step.ok());
    EXPECT_EQ(step.diagnostics.regime, orbitsim::KeplerOrbitRegime::Elliptic);
    EXPECT_TRUE(near_vec_abs(step.state.position_m, s0.position_m, 1e-3));
    EXPECT_TRUE(near_vec_abs(step.state.velocity_mps, s0.velocity_mps, 1e-6));
}

TEST(Kepler, EccentricOrbitRoundTripsOnePeriod)
{
    orbitsim::OrbitalElements el{};
    el.semi_major_axis_m = 10'000'000.0;
    el.eccentricity = 0.4;
    el.inclination_rad = 0.5;
    el.raan_rad = 1.0;
    el.arg_periapsis_rad = 0.25;
    el.true_anomaly_rad = 1.3;

    const orbitsim::State s0 = orbitsim::relative_state_from_orbital_elements(kEarthMu, el);
    const double a = el.semi_major_axis_m;
    const double period_s = 2.0 * std::acos(-1.0) * std::sqrt((a * a * a) / kEarthMu);

    const orbitsim::KeplerPropagationResult step =
            orbitsim::propagate_kepler_universal_safe(kEarthMu, s0, period_s);

    ASSERT_TRUE(step.ok());
    EXPECT_TRUE(near_vec_abs(step.state.position_m, s0.position_m, 1e-2));
    EXPECT_TRUE(near_vec_abs(step.state.velocity_mps, s0.velocity_mps, 1e-5));
}

TEST(Kepler, HighEccentricOrbitBuildsDenseSamplesThroughOnePeriod)
{
    orbitsim::OrbitalElements el{};
    el.semi_major_axis_m = 100'000'000.0;
    el.eccentricity = 0.95;
    el.inclination_rad = 0.35;
    el.raan_rad = 0.4;
    el.arg_periapsis_rad = 0.9;
    el.true_anomaly_rad = 1.57;

    const orbitsim::State s0 = orbitsim::relative_state_from_orbital_elements(kEarthMu, el);
    const double period_s =
            2.0 * std::acos(-1.0) *
            std::sqrt((el.semi_major_axis_m * el.semi_major_axis_m * el.semi_major_axis_m) /
                      kEarthMu);

    orbitsim::KeplerArc arc{};
    arc.mu_m3_s2 = kEarthMu;
    arc.primary_body_id = 42;
    arc.t0_s = 0.0;
    arc.t1_s = period_s;
    arc.state0_relative = s0;

    orbitsim::KeplerTrajectoryOptions opt{};
    opt.sample_dt_s = period_s / 19'999.0;
    opt.max_samples = 20'000;

    orbitsim::KeplerTrajectoryDiagnostics diagnostics{};
    const std::vector<orbitsim::KeplerArcSample> samples =
            orbitsim::build_kepler_arc_samples(arc, opt, &diagnostics);

    ASSERT_EQ(samples.size(), 20'000u);
    EXPECT_EQ(diagnostics.accepted_samples, samples.size());
    EXPECT_EQ(diagnostics.first_failure, orbitsim::KeplerStatus::Ok);
    EXPECT_TRUE(near_abs(samples.front().t_s, 0.0, 1e-12));
    EXPECT_TRUE(near_abs(samples.back().t_s, period_s, 1e-6));
    for (const orbitsim::KeplerArcSample &sample : samples)
    {
        ASSERT_TRUE(sample.ok());
        EXPECT_TRUE(finite_state(sample.state_relative));
    }
}

TEST(Kepler, HyperbolicEscapePropagatesFinite)
{
    const orbitsim::State s0 = orbitsim::make_state({7'000'000.0, 0.0, 0.0}, {0.0, 12'000.0, 0.0});

    const orbitsim::KeplerPropagationResult step =
            orbitsim::propagate_kepler_universal_safe(kEarthMu, s0, 3600.0);

    ASSERT_TRUE(step.ok());
    EXPECT_EQ(step.diagnostics.regime, orbitsim::KeplerOrbitRegime::Hyperbolic);
    EXPECT_TRUE(finite_state(step.state));
    EXPECT_GT(glm::length(step.state.position_m), glm::length(s0.position_m));
}

TEST(Kepler, NearParabolicCasePropagatesFinite)
{
    constexpr double r_m = 7'000'000.0;
    const double v_mps = std::sqrt(2.0 * kEarthMu / r_m) * 0.999999;
    const orbitsim::State s0 = orbitsim::make_state({r_m, 0.0, 0.0}, {0.0, v_mps, 0.0});

    const orbitsim::KeplerPropagationResult step =
            orbitsim::propagate_kepler_universal_safe(kEarthMu, s0, 1200.0);

    ASSERT_TRUE(step.ok());
    EXPECT_EQ(step.diagnostics.regime, orbitsim::KeplerOrbitRegime::NearParabolic);
    EXPECT_TRUE(finite_state(step.state));
}

TEST(Kepler, NegativeDtRoundTrips)
{
    constexpr double r_m = 7'000'000.0;
    const double v_mps = std::sqrt(kEarthMu / r_m);
    const orbitsim::State s0 = orbitsim::make_state({r_m, 0.0, 0.0}, {0.0, v_mps, 0.0});

    const orbitsim::KeplerPropagationResult forward =
            orbitsim::propagate_kepler_universal_safe(kEarthMu, s0, 1800.0);
    ASSERT_TRUE(forward.ok());

    const orbitsim::KeplerPropagationResult back =
            orbitsim::propagate_kepler_universal_safe(kEarthMu, forward.state, -1800.0);
    ASSERT_TRUE(back.ok());

    EXPECT_TRUE(near_vec_abs(back.state.position_m, s0.position_m, 1e-3));
    EXPECT_TRUE(near_vec_abs(back.state.velocity_mps, s0.velocity_mps, 1e-6));
}

TEST(Kepler, ReportsStumpffOverflow)
{
    const orbitsim::State s0 = orbitsim::make_state({7'000'000.0, 0.0, 0.0}, {0.0, 12'000.0, 0.0});

    const orbitsim::KeplerPropagationResult step =
            orbitsim::propagate_kepler_universal_safe(kEarthMu, s0, 1.0e12);

    EXPECT_FALSE(step.ok());
    EXPECT_EQ(step.diagnostics.status, orbitsim::KeplerStatus::StumpffOverflow);
}

TEST(Kepler, BuildsArcSamplesIncludingEndpoints)
{
    constexpr double r_m = 7'000'000.0;
    const double v_mps = std::sqrt(kEarthMu / r_m);
    orbitsim::KeplerArc arc{};
    arc.mu_m3_s2 = kEarthMu;
    arc.primary_body_id = 42;
    arc.t0_s = 10.0;
    arc.t1_s = 40.0;
    arc.state0_relative = orbitsim::make_state({r_m, 0.0, 0.0}, {0.0, v_mps, 0.0});

    orbitsim::KeplerTrajectoryOptions opt{};
    opt.sample_dt_s = 10.0;
    opt.max_samples = 16;

    orbitsim::KeplerTrajectoryDiagnostics diagnostics{};
    const std::vector<orbitsim::KeplerArcSample> samples =
            orbitsim::build_kepler_arc_samples(arc, opt, &diagnostics);

    ASSERT_EQ(samples.size(), 4u);
    EXPECT_TRUE(near_abs(samples[0].t_s, 10.0, 1e-12));
    EXPECT_TRUE(near_abs(samples[1].t_s, 20.0, 1e-12));
    EXPECT_TRUE(near_abs(samples[2].t_s, 30.0, 1e-12));
    EXPECT_TRUE(near_abs(samples[3].t_s, 40.0, 1e-12));
    EXPECT_EQ(diagnostics.accepted_samples, samples.size());
    EXPECT_EQ(diagnostics.first_failure, orbitsim::KeplerStatus::Ok);

    opt.include_start = false;
    const std::vector<orbitsim::KeplerArcSample> samples_without_start =
            orbitsim::build_kepler_arc_samples(arc, opt);
    ASSERT_EQ(samples_without_start.size(), 3u);
    EXPECT_TRUE(near_abs(samples_without_start[0].t_s, 20.0, 1e-12));
    EXPECT_TRUE(near_abs(samples_without_start[1].t_s, 30.0, 1e-12));
    EXPECT_TRUE(near_abs(samples_without_start[2].t_s, 40.0, 1e-12));
}

TEST(Kepler, RtnImpulseAppliesTangentialDeltaV)
{
    constexpr double r_m = 7'000'000.0;
    const double v_mps = std::sqrt(kEarthMu / r_m);
    const orbitsim::State pre = orbitsim::make_state({r_m, 0.0, 0.0}, {0.0, v_mps, 0.0});

    const orbitsim::State post = orbitsim::apply_impulse_rtn(pre, {0.0, 10.0, 0.0});

    EXPECT_TRUE(near_vec_abs(post.position_m, pre.position_m, 1e-12));
    EXPECT_TRUE(near_vec_abs(post.velocity_mps, pre.velocity_mps + orbitsim::Vec3{0.0, 10.0, 0.0}, 1e-12));
}

TEST(Kepler, ManeuveredArcsSplitAtImpulse)
{
    constexpr double r_m = 7'000'000.0;
    const double v_mps = std::sqrt(kEarthMu / r_m);
    orbitsim::KeplerArc arc{};
    arc.mu_m3_s2 = kEarthMu;
    arc.primary_body_id = 42;
    arc.t0_s = 0.0;
    arc.t1_s = 1000.0;
    arc.state0_relative = orbitsim::make_state({r_m, 0.0, 0.0}, {0.0, v_mps, 0.0});

    const std::vector<orbitsim::KeplerImpulse> impulses{
            orbitsim::KeplerImpulse{.t_s = 100.0, .dv_rtn_mps = {0.0, 10.0, 0.0}},
    };

    orbitsim::KeplerManeuverDiagnostics diagnostics{};
    const std::vector<orbitsim::KeplerArc> arcs =
            orbitsim::build_maneuvered_kepler_arcs(arc, impulses, {}, &diagnostics);

    ASSERT_EQ(arcs.size(), 2u);
    EXPECT_TRUE(near_abs(arcs[0].t0_s, 0.0, 1e-12));
    EXPECT_TRUE(near_abs(arcs[0].t1_s, 100.0, 1e-12));
    EXPECT_TRUE(near_abs(arcs[1].t0_s, 100.0, 1e-12));
    EXPECT_TRUE(near_abs(arcs[1].t1_s, 1000.0, 1e-12));
    EXPECT_EQ(diagnostics.impulses_applied, 1u);
    EXPECT_EQ(diagnostics.arcs_built, 2u);

    const orbitsim::KeplerArcSample pre_impulse = orbitsim::sample_kepler_arc_state(arc, 100.0);
    ASSERT_TRUE(pre_impulse.ok());
    const orbitsim::State expected_post = orbitsim::apply_impulse_rtn(pre_impulse.state_relative, {0.0, 10.0, 0.0});
    EXPECT_TRUE(near_vec_abs(arcs[1].state0_relative.position_m, expected_post.position_m, 1e-6));
    EXPECT_TRUE(near_vec_abs(arcs[1].state0_relative.velocity_mps, expected_post.velocity_mps, 1e-9));
}
