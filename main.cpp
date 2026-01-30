#include <orbitsim/orbitsim.hpp>

#include <cmath>
#include <cstdio>

int main()
{
    orbitsim::GameSimulation sim;

    constexpr double Ms = 1.98847e30; // kg
    constexpr double Me = 5.9722e24; // kg
    constexpr double AU = 1.495978707e11; // m

    const orbitsim::Vec3 r_rel{AU, 0.0, 0.0};
    const double mu = orbitsim::kGravitationalConstant_SI * (Ms + Me);
    const double v_rel_mag = std::sqrt(mu / glm::length(r_rel));
    const orbitsim::Vec3 v_rel{0.0, v_rel_mag, 0.0};

    // Barycentric initial conditions for a simple Sun/Earth circular orbit.
    const double m_tot = Ms + Me;
    orbitsim::MassiveBody sun{
            .mass_kg = Ms,
            .radius_m = 6.9634e8,
            .state =
                    {
                            .position_m = -(Me / m_tot) * r_rel,
                            .velocity_mps = -(Me / m_tot) * v_rel,
                            .spin = {.axis = {0.0, 1.0, 0.0}, .angle_rad = 0.0, .rate_rad_per_s = 2.9e-6},
                    },
    };
    orbitsim::MassiveBody earth{
            .mass_kg = Me,
            .radius_m = 6.371e6,
            .atmosphere_top_height_m = 1.0e5,
            .soi_radius_m = 9.25e8,
            .state =
                    {
                            .position_m = (Ms / m_tot) * r_rel,
                            .velocity_mps = (Ms / m_tot) * v_rel,
                            .spin = {.axis = {0.0, 1.0, 0.0}, .angle_rad = 0.0, .rate_rad_per_s = 7.2921159e-5},
                    },
    };

    sim.massive_bodies().push_back(sun);
    sim.massive_bodies().push_back(earth);

    // A "spacecraft" near Earth, initially co-moving.
    orbitsim::Spacecraft sc{
            .state =
                    {
                            .position_m = earth.state.position_m + orbitsim::Vec3{42'000'000.0, 0.0, 0.0},
                            .velocity_mps = earth.state.velocity_mps + orbitsim::Vec3{0.0, 0.0, 0.0},
                            .spin = {.axis = {0.0, 0.0, 1.0}, .angle_rad = 0.0, .rate_rad_per_s = 0.0},
                    },
            .dry_mass_kg = 800.0,
            .prop_mass_kg = 200.0,
            .engines =
                    {
                            {.max_thrust_N = 2'000.0, .isp_s = 300.0, .min_throttle_0_1 = 0.05},
                    },
    };
    sim.spacecraft().push_back(sc);

    // A short prograde burn in RTN (relative to Earth as primary).
    sim.maneuver_plan().segments.push_back(orbitsim::BurnSegment{
            .t_start_s = 3600.0,
            .t_end_s = 4200.0,
            .primary_index = 1,
            .dir_rtn_unit = {0.0, 1.0, 0.0},
            .throttle_0_1 = 1.0,
            .engine_index = 0,
            .spacecraft_index = 0,
    });

    constexpr double dt = 60.0; // s
    constexpr int steps = 24 * 60; // 1 day

    for (int i = 0; i < steps; ++i)
    {
        sim.step(dt);
    }

    // Trajectory sampling helpers (useful for rendering orbit lines):
    // - Build a celestial ephemeris once (massive bodies only)
    // - Sample any body/spacecraft trajectory against that ephemeris
    orbitsim::TrajectoryOptions traj_opt{};
    traj_opt.duration_s = 6.0 * 3600.0;
    traj_opt.sample_dt_s = 5.0 * 60.0;
    traj_opt.celestial_dt_s = 60.0;
    traj_opt.max_samples = 4096;

    const orbitsim::CelestialEphemeris eph = orbitsim::build_celestial_ephemeris(sim, traj_opt);

    // Earth orbit around Sun (relative to the Sun).
    traj_opt.origin_body_index = 0;
    const std::vector<orbitsim::TrajectorySample> earth_traj = orbitsim::predict_body_trajectory(sim, eph, 1, traj_opt);

    // Spacecraft orbit around Earth (relative to Earth).
    traj_opt.origin_body_index = 1;
    const std::vector<orbitsim::TrajectorySample> sc_traj =
            orbitsim::predict_spacecraft_trajectory(sim, eph, 0, traj_opt);

    std::printf("\n--- trajectory csv (t_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps) ---\n");
    std::printf("earth_rel_sun\n");
    for (const auto &s: earth_traj)
    {
        std::printf("%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s.t_s, s.position_m.x, s.position_m.y, s.position_m.z,
                    s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
    }
    std::printf("sc_rel_earth\n");
    for (const auto &s: sc_traj)
    {
        std::printf("%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s.t_s, s.position_m.x, s.position_m.y, s.position_m.z,
                    s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
    }

    const auto &earth_out = sim.massive_bodies()[1].state;
    const auto &sc_out = sim.spacecraft()[0].state;
    std::printf("t = %.0f s\n", sim.time_s());
    std::printf("earth pos (km): %.3f %.3f %.3f\n", earth_out.position_m.x / 1000.0, earth_out.position_m.y / 1000.0,
                earth_out.position_m.z / 1000.0);
    std::printf("sc pos (km):    %.3f %.3f %.3f\n", sc_out.position_m.x / 1000.0, sc_out.position_m.y / 1000.0,
                sc_out.position_m.z / 1000.0);
    std::printf("earth spin axis: %.3f %.3f %.3f, angle(rad)=%.6f\n", earth_out.spin.axis.x, earth_out.spin.axis.y,
                earth_out.spin.axis.z, earth_out.spin.angle_rad);

    const orbitsim::Vec3 r_sc_rel = sc_out.position_m - earth_out.position_m;
    const orbitsim::Vec3 v_sc_rel = sc_out.velocity_mps - earth_out.velocity_mps;
    const double mu_earth = orbitsim::kGravitationalConstant_SI * earth.mass_kg;
    const orbitsim::OrbitalElements el = orbitsim::orbital_elements_from_relative_state(mu_earth, r_sc_rel, v_sc_rel);
    std::printf("sc elements wrt Earth: a(km)=%.3f e=%.6f i(deg)=%.6f\n", el.semi_major_axis_m / 1000.0,
                el.eccentricity, el.inclination_rad * 180.0 / std::acos(-1.0));

    return 0;
}
