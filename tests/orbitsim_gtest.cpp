#include <orbitsim/orbitsim.hpp>

#include <gtest/gtest.h>

#include <cmath>
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

    sim.massive_bodies().push_back(orbitsim::MassiveBody{
            .mass_kg = 0.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {0.0, 0.0, 0.0},
                            .velocity_mps = {0.0, 0.0, 0.0},
                    },
    });

    orbitsim::Spacecraft sc{};
    sc.state.position_m = {3.0, 0.0, 0.0};
    sc.state.velocity_mps = {0.0, 0.0, 0.0};
    sc.dry_mass_kg = 1.0;
    sc.prop_mass_kg = 0.1;
    sc.engines.push_back(orbitsim::Engine{.max_thrust_N = 10.0, .isp_s = 1.0, .min_throttle_0_1 = 0.0});
    sim.spacecraft().push_back(sc);
    sim.spacecraft().push_back(sc);

    sim.maneuver_plan().segments.push_back(orbitsim::BurnSegment{
            .t_start_s = 0.0,
            .t_end_s = 10.0,
            .primary_index = 0,
            .dir_rtn_unit = {1.0, 0.0, 0.0},
            .throttle_0_1 = 1.0,
            .engine_index = 0,
            .spacecraft_index = 1,
    });

    sim.step(1.0);

    const orbitsim::Spacecraft &out0 = sim.spacecraft()[0];
    const orbitsim::Spacecraft &out1 = sim.spacecraft()[1];
    EXPECT_TRUE(near_abs(out0.prop_mass_kg, 0.1, 1e-12));
    EXPECT_TRUE(near_abs(out0.state.velocity_mps.x, 0.0, 1e-12));
    EXPECT_GE(out1.prop_mass_kg, 0.0);
    EXPECT_TRUE(near_abs(out1.prop_mass_kg, 0.0, 1e-12));
    EXPECT_TRUE(std::isfinite(out1.state.velocity_mps.x));
    EXPECT_GT(out1.state.velocity_mps.x, 0.0);

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

    sim.massive_bodies().push_back(orbitsim::MassiveBody{
            .mass_kg = 0.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {100.0, 0.0, 0.0},
                            .velocity_mps = {0.5, 0.0, 0.0},
                    },
    });

    orbitsim::Spacecraft sc{};
    sc.state.position_m = {110.0, 0.0, 0.0};
    sc.state.velocity_mps = {1.0, 0.0, 0.0};
    sc.dry_mass_kg = 1.0;
    sc.prop_mass_kg = 0.0;
    sim.spacecraft().push_back(sc);

    orbitsim::TrajectoryOptions opt{};
    opt.duration_s = 10.0;
    opt.sample_dt_s = 1.0;
    opt.max_samples = 64;
    opt.include_start = true;
    opt.include_end = true;
    opt.origin_body_index = 0;

    const std::vector<orbitsim::TrajectorySample> traj = orbitsim::predict_spacecraft_trajectory(sim, 0, opt);
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

    sim.massive_bodies().push_back(orbitsim::MassiveBody{
            .mass_kg = 1.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {0.0, 0.0, 0.0},
                            .velocity_mps = {1.0, 0.0, 0.0},
                    },
    });

    orbitsim::TrajectoryOptions opt{};
    opt.duration_s = 4.0;
    opt.sample_dt_s = 1.0;
    opt.celestial_dt_s = 2.0;
    opt.max_samples = 64;

    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris(sim, opt);
    ASSERT_FALSE(eph.empty());
    ASSERT_EQ(eph.segments.size(), 2u);

    const orbitsim::State s_t1 = eph.body_state_at(0, 1.0);
    EXPECT_TRUE(near_vec_abs(s_t1.position_m, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_t1.velocity_mps, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));

    const orbitsim::State s_t3 = eph.body_state_at(0, 3.0);
    EXPECT_TRUE(near_vec_abs(s_t3.position_m, orbitsim::Vec3{3.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_t3.velocity_mps, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));

    const orbitsim::State s_before = eph.body_state_at(0, -1.0);
    EXPECT_TRUE(near_vec_abs(s_before.position_m, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_before.velocity_mps, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));

    const orbitsim::State s_after = eph.body_state_at(0, 10.0);
    EXPECT_TRUE(near_vec_abs(s_after.position_m, orbitsim::Vec3{4.0, 0.0, 0.0}, 1e-12));
    EXPECT_TRUE(near_vec_abs(s_after.velocity_mps, orbitsim::Vec3{1.0, 0.0, 0.0}, 1e-12));
}

TEST(Trajectories, BodyTrajectoryRelativeToSelfIsZero)
{
    orbitsim::GameSimulation::Config cfg{};
    cfg.gravitational_constant = 0.0;
    cfg.enable_events = false;

    orbitsim::GameSimulation sim(cfg);

    sim.massive_bodies().push_back(orbitsim::MassiveBody{
            .mass_kg = 1.0,
            .radius_m = 1.0,
            .state =
                    {
                            .position_m = {5.0, 0.0, 0.0},
                            .velocity_mps = {2.0, 0.0, 0.0},
                    },
    });

    orbitsim::TrajectoryOptions opt{};
    opt.duration_s = 5.0;
    opt.sample_dt_s = 1.0;
    opt.celestial_dt_s = 1.0;
    opt.max_samples = 64;
    opt.origin_body_index = 0;

    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris(sim, opt);
    const std::vector<orbitsim::TrajectorySample> traj = orbitsim::predict_body_trajectory(sim, eph, 0, opt);
    ASSERT_EQ(traj.size(), 6u);
    for (const auto &s: traj)
    {
        EXPECT_TRUE(near_vec_abs(s.position_m, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-12));
        EXPECT_TRUE(near_vec_abs(s.velocity_mps, orbitsim::Vec3{0.0, 0.0, 0.0}, 1e-12));
    }
}
