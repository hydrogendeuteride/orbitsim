#include <orbitsim/orbitsim.hpp>

#include <cmath>
#include <cstdio>

int main()
{
    using namespace orbitsim;

    GameSimulation::Config cfg{};
    cfg.spacecraft_integrator.max_step_s = 2.0;
    cfg.spacecraft_integrator.max_substeps = 512;
    GameSimulation sim(cfg);

    constexpr double Ms = 1.98847e30; // kg
    constexpr double Me = 5.9722e24; // kg
    constexpr double Mm = 7.342e22; // kg
    constexpr double AU = 1.495978707e11; // m

    const double m_em = Me + Mm;

    // Sun <-> Earth-Moon barycenter initial conditions (tilted for 3D visualization).
    const double inc_em_rad = glm::radians(10.0);
    const auto [sun_state, em_bary_state] = two_body_circular_barycentric(Ms, m_em, AU, inc_em_rad);

    MassiveBody sun{
            .mass_kg = Ms,
            .radius_m = 6.9634e8,
            .state = make_state(sun_state.position_m, sun_state.velocity_mps, {0.0, 1.0, 0.0}, 2.9e-6),
    };

    // Earth <-> Moon circular orbit around their barycenter (slightly inclined).
    constexpr double d_em_m = 384'400'000.0; // mean distance [m]
    const double inc_moon_rad = glm::radians(5.0);
    const auto [earth_rel, moon_rel] = two_body_circular_barycentric(Me, Mm, d_em_m, inc_moon_rad);

    MassiveBody earth{
            .mass_kg = Me,
            .radius_m = 6.371e6,
            .atmosphere_top_height_m = 1.0e5,
            .soi_radius_m = 9.25e8,
            .state = make_state(em_bary_state.position_m + earth_rel.position_m,
                                em_bary_state.velocity_mps + earth_rel.velocity_mps, {0.0, 1.0, 0.0}, 7.2921159e-5),
    };
    MassiveBody moon{
            .mass_kg = Mm,
            .radius_m = 1.7374e6,
            .soi_radius_m = 6.61e7,
            .state = make_state(em_bary_state.position_m + moon_rel.position_m,
                                em_bary_state.velocity_mps + moon_rel.velocity_mps, {0.0, 1.0, 0.0}, 2.66e-6),
    };

    const GameSimulation::BodyHandle sun_h = sim.create_body(std::move(sun));
    const GameSimulation::BodyHandle earth_h = sim.create_body(std::move(earth));
    const GameSimulation::BodyHandle moon_h = sim.create_body(std::move(moon));
    if (!sun_h.valid() || !earth_h.valid() || !moon_h.valid())
    {
        return 1;
    }
    const BodyId sun_id = sun_h.id;
    const BodyId earth_id = earth_h.id;
    const BodyId moon_id = moon_h.id;
    const MassiveBody *earth_ptr = sim.body_by_id(earth_id);
    const MassiveBody *moon_ptr = sim.body_by_id(moon_id);
    if (earth_ptr == nullptr || moon_ptr == nullptr)
    {
        return 1;
    }

    // A spacecraft in LEO with a few planned burns (Earth injection + mid-course correction + lunar capture).
    const double sc_altitude_m = 300'000.0;
    const double sc_inc_rad = glm::radians(28.5);
    const auto sc_orbit = circular_orbit_relative_state(earth_ptr->mass_kg, earth_ptr->radius_m + sc_altitude_m, sc_inc_rad);

    Spacecraft sc{
            .state = make_state(earth_ptr->state.position_m + sc_orbit.position_m,
                                earth_ptr->state.velocity_mps + sc_orbit.velocity_mps,
                                {0.0, 0.0, 1.0}, 0.0),
            .dry_mass_kg = 800.0,
            .prop_mass_kg = 200.0,
            .engines =
                    {
                            {.max_thrust_N = 2'000.0, .isp_s = 300.0, .min_throttle_0_1 = 0.05},
                    },
    };
    const GameSimulation::SpacecraftHandle sc_h = sim.create_spacecraft(std::move(sc));
    if (!sc_h.valid())
    {
        return 1;
    }
    const SpacecraftId sc_id = sc_h.id;

    // Injection burn: prograde in RTN (Earth primary).
    sim.maneuver_plan().segments.push_back(
            burn().start(hours(2.0))
                  .duration(minutes(20.0))
                  .prograde()
                  .primary(earth_id)
                  .spacecraft(sc_id));

    // Mid-course correction: small normal component in RTN (Earth primary).
    sim.maneuver_plan().segments.push_back(
            burn().start(days(2.0))
                  .duration(minutes(5.0))
                  .normal()
                  .throttle(0.4)
                  .primary(earth_id)
                  .spacecraft(sc_id));

    // Lunar capture attempt: retrograde in RTN (Moon primary).
    sim.maneuver_plan().segments.push_back(
            burn().start(days(4.5))
                  .duration(minutes(15.0))
                  .retrograde()
                  .throttle(0.7)
                  .primary(moon_id)
                  .spacecraft(sc_id));

    // Trajectory sampling helpers (useful for rendering orbit lines):
    // - Build a celestial ephemeris once (massive bodies only)
    // - Sample any body/spacecraft trajectory against that ephemeris
    const auto traj_opt = trajectory_options()
            .duration(days(20.0))
            .sample_dt(minutes(10.0))
            .celestial_dt(minutes(5.0))
            .max_samples(100'000);

    const CelestialEphemeris eph = build_celestial_ephemeris(sim, traj_opt);

    // Earth/Moon orbits around Sun (relative to the Sun).
    const auto traj_opt_sun = trajectory_options()
            .duration(days(20.0))
            .sample_dt(minutes(10.0))
            .celestial_dt(minutes(5.0))
            .max_samples(100'000)
            .origin(sun_id);
    const std::vector<TrajectorySample> earth_traj = predict_body_trajectory(sim, eph, earth_id, traj_opt_sun);
    const std::vector<TrajectorySample> moon_traj_sun = predict_body_trajectory(sim, eph, moon_id, traj_opt_sun);
    const std::vector<TrajectorySample> sc_traj_sun = predict_spacecraft_trajectory(sim, eph, sc_id, traj_opt_sun);

    // Moon + spacecraft around Earth (relative to Earth).
    const auto traj_opt_earth = trajectory_options()
            .duration(days(20.0))
            .sample_dt(minutes(10.0))
            .celestial_dt(minutes(5.0))
            .max_samples(100'000)
            .origin(earth_id);
    const std::vector<TrajectorySample> moon_traj_earth = predict_body_trajectory(sim, eph, moon_id, traj_opt_earth);
    const std::vector<TrajectorySample> sc_traj = predict_spacecraft_trajectory(sim, eph, sc_id, traj_opt_earth);

    // High-resolution spacecraft path near Earth for smooth LEO rendering.
    const auto earth_hi = trajectory_options()
            .duration(hours(8.0))
            .sample_dt(seconds(5.0))
            .spacecraft_sample_dt(seconds(5.0))
            .celestial_dt(seconds(30.0))
            .max_samples(50'000)
            .origin(earth_id);
    const std::vector<TrajectorySample> sc_traj_earth_hi = predict_spacecraft_trajectory(sim, eph, sc_id, earth_hi);

    // Spacecraft around Moon (relative to Moon).
    const auto traj_opt_moon = trajectory_options()
            .duration(days(20.0))
            .sample_dt(minutes(10.0))
            .celestial_dt(minutes(5.0))
            .max_samples(100'000)
            .origin(moon_id);
    const std::vector<TrajectorySample> sc_traj_moon = predict_spacecraft_trajectory(sim, eph, sc_id, traj_opt_moon);

    std::printf("\n--- trajectory csv (t_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps) ---\n");
    std::printf("earth_rel_sun\n");
    for (const auto &s: earth_traj)
    {
        std::printf("%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s.t_s, s.position_m.x, s.position_m.y, s.position_m.z,
                    s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
    }
    std::printf("moon_rel_sun\n");
    for (const auto &s: moon_traj_sun)
    {
        std::printf("%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s.t_s, s.position_m.x, s.position_m.y, s.position_m.z,
                    s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
    }
    std::printf("sc_rel_sun\n");
    for (const auto &s: sc_traj_sun)
    {
        std::printf("%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s.t_s, s.position_m.x, s.position_m.y, s.position_m.z,
                    s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
    }
    std::printf("moon_rel_earth\n");
    for (const auto &s: moon_traj_earth)
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
    std::printf("sc_rel_earth_hi\n");
    for (const auto &s: sc_traj_earth_hi)
    {
        std::printf("%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s.t_s, s.position_m.x, s.position_m.y, s.position_m.z,
                    s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
    }
    std::printf("sc_rel_moon\n");
    for (const auto &s: sc_traj_moon)
    {
        std::printf("%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s.t_s, s.position_m.x, s.position_m.y, s.position_m.z,
                    s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
    }

    return 0;
}
