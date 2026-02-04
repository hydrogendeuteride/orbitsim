#include <orbitsim/orbitsim.hpp>

#include <cmath>
#include <cstdio>
#include <limits>

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

    // -------------------------------------------------------------------------
    // Basic body metadata for downstream parsers/visualizers
    // -------------------------------------------------------------------------
    std::printf("\n--- bodies ---\n");
    std::printf("bodies\n");
    auto print_body = [&](BodyId id, const char *name) {
        const MassiveBody *b = sim.body_by_id(id);
        if (b == nullptr)
        {
            return;
        }
        std::printf("%u,%s,%.6e,%.6e,%.6e\n",
                    static_cast<unsigned>(b->id),
                    name,
                    b->radius_m,
                    b->atmosphere_top_height_m,
                    b->soi_radius_m);
    };
    print_body(sun_id, "sun");
    print_body(earth_id, "earth");
    print_body(moon_id, "moon");

    // A spacecraft in LEO with a few planned burns (Earth injection + mid-course correction + lunar capture).
    const double sc_altitude_m = 300'000.0;
    const double sc_inc_rad = glm::radians(28.5);
    const auto sc_orbit = circular_orbit_relative_state(earth_ptr->mass_kg, earth_ptr->radius_m + sc_altitude_m, sc_inc_rad);

    Spacecraft sc{
            .state = make_state(earth_ptr->state.position_m + sc_orbit.position_m,
                                earth_ptr->state.velocity_mps + sc_orbit.velocity_mps,
                                {0.0, 0.0, 1.0}, 0.0),
            .dry_mass_kg = 800.0,
            .prop_mass_kg = 2000.0,
            .engines =
                    {
                            {.max_thrust_N = 5'000.0, .isp_s = 300.0, .min_throttle_0_1 = 0.05},
                    },
    };
    const GameSimulation::SpacecraftHandle sc_h = sim.create_spacecraft(std::move(sc));
    if (!sc_h.valid())
    {
        return 1;
    }
    const SpacecraftId sc_id = sc_h.id;

    // -------------------------------------------------------------------------
    // Extra "event probe" spacecraft for testing event detection/visualization.
    // -------------------------------------------------------------------------

    // Reentry/impact probe near Earth (guaranteed to cross atmosphere + impact soon).
    GameSimulation::SpacecraftHandle sc_reentry_h{};
    {
        const double alt_m = 120'000.0;
        const Vec3 r0 = Vec3{earth_ptr->radius_m + alt_m, 0.0, 0.0};
        const Vec3 v0 = Vec3{-250.0, 0.0, 0.0}; // inward in Earth-centered inertial
        Spacecraft sc_reentry{
                .state = make_state(earth_ptr->state.position_m + r0, earth_ptr->state.velocity_mps + v0),
                .dry_mass_kg = 100.0,
                .prop_mass_kg = 0.0,
        };
        sc_reentry_h = sim.create_spacecraft(std::move(sc_reentry));
    }

    // Moon SOI probe: starts just outside Moon SOI and drifts inward, triggering SOI enter event.
    GameSimulation::SpacecraftHandle sc_moon_soi_h{};
    {
        const Vec3 r0 = Vec3{moon_ptr->soi_radius_m + 1.0e6, 0.0, 0.0};
        const Vec3 v0 = Vec3{-200.0, 0.0, 0.0};
        Spacecraft sc_moon_soi{
                .state = make_state(moon_ptr->state.position_m + r0, moon_ptr->state.velocity_mps + v0),
                .dry_mass_kg = 100.0,
                .prop_mass_kg = 0.0,
        };
        sc_moon_soi_h = sim.create_spacecraft(std::move(sc_moon_soi));
    }

    // Injection burn: prograde in RTN (Earth primary).
    sim.maneuver_plan().segments.push_back(
            burn().start(hours(2.0))
                  .duration(minutes(18.0))
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

    // Predict inertial trajectories once, then transform for plotting.
    const std::vector<TrajectorySample> earth_traj_inertial = predict_body_trajectory(sim, eph, earth_id, traj_opt);
    const std::vector<TrajectorySample> moon_traj_inertial = predict_body_trajectory(sim, eph, moon_id, traj_opt);
    const std::vector<TrajectorySample> sc_traj_inertial = predict_spacecraft_trajectory(sim, eph, sc_id, traj_opt);

    // Earth/Moon/spacecraft around Sun (sun-centered inertial coordinates).
    const TrajectoryFrameSpec sun_frame = TrajectoryFrameSpec::body_centered_inertial(sun_id);
    const std::vector<TrajectorySample> earth_traj =
            trajectory_to_frame_spec(earth_traj_inertial, eph, sim.massive_bodies(), sun_frame);
    const std::vector<TrajectorySample> moon_traj_sun =
            trajectory_to_frame_spec(moon_traj_inertial, eph, sim.massive_bodies(), sun_frame);
    const std::vector<TrajectorySample> sc_traj_sun =
            trajectory_to_frame_spec(sc_traj_inertial, eph, sim.massive_bodies(), sun_frame);

    // Moon + spacecraft around Earth (earth-centered inertial coordinates).
    const TrajectoryFrameSpec earth_frame = TrajectoryFrameSpec::body_centered_inertial(earth_id);
    const std::vector<TrajectorySample> moon_traj_earth =
            trajectory_to_frame_spec(moon_traj_inertial, eph, sim.massive_bodies(), earth_frame);
    const std::vector<TrajectorySample> sc_traj =
            trajectory_to_frame_spec(sc_traj_inertial, eph, sim.massive_bodies(), earth_frame);

    // High-resolution spacecraft path near Earth for smooth LEO rendering.
    const auto earth_hi = trajectory_options()
            .duration(hours(8.0))
            .sample_dt(seconds(5.0))
            .spacecraft_sample_dt(seconds(5.0))
            .celestial_dt(seconds(30.0))
            .max_samples(50'000);
    const std::vector<TrajectorySample> sc_traj_earth_hi_inertial = predict_spacecraft_trajectory(sim, eph, sc_id, earth_hi);
    const std::vector<TrajectorySample> sc_traj_earth_hi =
            trajectory_to_frame_spec(sc_traj_earth_hi_inertial, eph, sim.massive_bodies(), earth_frame);

    // Spacecraft around Moon (moon-centered inertial coordinates).
    const TrajectoryFrameSpec moon_frame = TrajectoryFrameSpec::body_centered_inertial(moon_id);
    const std::vector<TrajectorySample> sc_traj_moon =
            trajectory_to_frame_spec(sc_traj_inertial, eph, sim.massive_bodies(), moon_frame);

    // Short-horizon trajectories for event probes (keep output small and avoid singular behavior past impact).
    const auto reentry_short_inertial = trajectory_options()
            .duration(minutes(12.0))
            .sample_dt(seconds(1.0))
            .celestial_dt(seconds(1.0))
            .max_samples(20'000);

    const auto moon_soi_short_inertial = trajectory_options()
            .duration(hours(4.0))
            .sample_dt(seconds(10.0))
            .celestial_dt(seconds(10.0))
            .max_samples(20'000);

    std::vector<TrajectorySample> sc_reentry_rel_earth;
    std::vector<TrajectorySample> sc_moon_soi_rel_moon;
    if (sc_reentry_h.valid())
    {
        const std::vector<TrajectorySample> sc_reentry_inertial =
                predict_spacecraft_trajectory(sim, eph, sc_reentry_h.id, reentry_short_inertial);
        sc_reentry_rel_earth = trajectory_to_frame_spec(sc_reentry_inertial, eph, sim.massive_bodies(), earth_frame);
    }
    if (sc_moon_soi_h.valid())
    {
        const std::vector<TrajectorySample> sc_moon_soi_inertial =
                predict_spacecraft_trajectory(sim, eph, sc_moon_soi_h.id, moon_soi_short_inertial);
        sc_moon_soi_rel_moon = trajectory_to_frame_spec(sc_moon_soi_inertial, eph, sim.massive_bodies(), moon_frame);
    }

    // -------------------------------------------------------------------------
    // Lambert example: Earth-centered transfer targeting the Moon's predicted position
    // -------------------------------------------------------------------------
    //
    // This is a planning/initial-guess tool:
    // - Solve Lambert in an Earth-centered 2-body model (mu = G * M_earth)
    // - Convert the Lambert departure velocity into inertial coordinates
    // - Propagate the resulting spacecraft with the library's N-body gravity (via ephemeris interpolation)
    //
    // The result will generally NOT hit the Moon exactly (because N-body != 2-body),
    // but it provides a useful starting point.
    const double t_depart_s = sim.time_s();
    const double dt_transfer_s = days(9.0);
    const double t_arrive_s = t_depart_s + dt_transfer_s;

    const State earth_at_depart = eph.body_state_at_by_id(earth_id, t_depart_s);
    const State earth_at_arrive = eph.body_state_at_by_id(earth_id, t_arrive_s);
    const State moon_at_arrive = eph.body_state_at_by_id(moon_id, t_arrive_s);

    const Vec3 r1_rel_earth_m = (sim.spacecraft_by_id(sc_id)->state.position_m - earth_at_depart.position_m);
    const Vec3 v0_rel_earth_mps = (sim.spacecraft_by_id(sc_id)->state.velocity_mps - earth_at_depart.velocity_mps);
    const Vec3 r2_rel_earth_m = (moon_at_arrive.position_m - earth_at_arrive.position_m);

    const double mu_earth_m3_s2 = kGravitationalConstant_SI * Me;
    LambertOptions lam_opt{};
    lam_opt.prograde = true;
    lam_opt.short_path = true;
    lam_opt.max_revolutions = 2;

    const std::vector<LambertSolution> lam_solutions =
            solve_lambert_universal(mu_earth_m3_s2, r1_rel_earth_m, r2_rel_earth_m, dt_transfer_s, lam_opt);

    std::optional<LambertSolution> best;
    double best_dv = std::numeric_limits<double>::infinity();
    for (const auto &sol: lam_solutions)
    {
        const double dv = glm::length(sol.v1_mps - v0_rel_earth_mps);
        if (dv < best_dv)
        {
            best_dv = dv;
            best = sol;
        }
    }

    GameSimulation::SpacecraftHandle sc_lambert_h{};
    if (best.has_value())
    {
        const Spacecraft *base_sc = sim.spacecraft_by_id(sc_id);
        if (base_sc == nullptr)
        {
            return 1;
        }
        Spacecraft sc_lambert = *base_sc;
        sc_lambert.id = kInvalidSpacecraftId;
        sc_lambert.engines.clear();
        sc_lambert.prop_mass_kg = 0.0;

        sc_lambert.state.position_m = earth_at_depart.position_m + r1_rel_earth_m;
        sc_lambert.state.velocity_mps = earth_at_depart.velocity_mps + best->v1_mps;

        sc_lambert_h = sim.create_spacecraft(std::move(sc_lambert));
    }

    // Lambert output sections (custom parser in visualize_lambert_transfer.py)
    std::printf("\n--- lambert (earth-centered) ---\n");
    std::printf("lambert_info\n");
    std::printf("t_depart_s,%.6f\n", t_depart_s);
    std::printf("t_arrive_s,%.6f\n", t_arrive_s);
    std::printf("dt_s,%.6f\n", dt_transfer_s);
    std::printf("mu_earth_m3_s2,%.6e\n", mu_earth_m3_s2);
    std::printf("num_solutions,%zu\n", lam_solutions.size());
    if (best.has_value())
    {
        std::printf("best_z,%.12e\n", best->z);
        std::printf("best_dv_mps,%.6f\n", best_dv);
        std::printf("best_v1_mps,%.6e,%.6e,%.6e\n", best->v1_mps.x, best->v1_mps.y, best->v1_mps.z);
        std::printf("best_v2_mps,%.6e,%.6e,%.6e\n", best->v2_mps.x, best->v2_mps.y, best->v2_mps.z);
    }
    else
    {
        std::printf("best_z,nan\n");
        std::printf("best_dv_mps,nan\n");
    }

    if (sc_lambert_h.valid())
    {
        const auto lam_opt_traj = trajectory_options()
                .duration(dt_transfer_s)
                .sample_dt(minutes(30.0))
                .celestial_dt(minutes(10.0))
                .max_samples(50'000);

        const std::vector<TrajectorySample> lam_sc_inertial =
                predict_spacecraft_trajectory(sim, eph, sc_lambert_h.id, lam_opt_traj);
        const std::vector<TrajectorySample> lam_moon_inertial =
                predict_body_trajectory(sim, eph, moon_id, lam_opt_traj);

        const std::vector<TrajectorySample> lam_sc_rel_earth =
                trajectory_to_frame_spec(lam_sc_inertial, eph, sim.massive_bodies(), earth_frame);
        const std::vector<TrajectorySample> lam_moon_rel_earth =
                trajectory_to_frame_spec(lam_moon_inertial, eph, sim.massive_bodies(), earth_frame);

        std::printf("lambert_sc_rel_earth\n");
        for (const auto &s: lam_sc_rel_earth)
        {
            std::printf("%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s.t_s, s.position_m.x, s.position_m.y, s.position_m.z,
                        s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
        }

        std::printf("lambert_moon_rel_earth\n");
        for (const auto &s: lam_moon_rel_earth)
        {
            std::printf("%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s.t_s, s.position_m.x, s.position_m.y, s.position_m.z,
                        s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
        }
    }

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

    // -------------------------------------------------------------------------
    // Rotating frame (synodic frame) visualization
    // -------------------------------------------------------------------------
    std::printf("\n--- rotating frame (Earth-Moon synodic) ---\n");

    // Build synodic frame at t=0
    const std::optional<SynodicFrame> syn = make_synodic_frame(*earth_ptr, *moon_ptr);
    if (syn.has_value())
    {
        std::printf("synodic_frame_t0\n");
        std::printf("separation_m,%.6e\n", syn->separation_m);
        std::printf("mu,%.6e\n", syn->mu);
        std::printf("omega_radps,%.6e,%.6e,%.6e\n",
                    syn->omega_inertial_radps.x,
                    syn->omega_inertial_radps.y,
                    syn->omega_inertial_radps.z);

        // Compute Lagrange points
        const std::optional<Cr3bpLagrangePoints> lpts = cr3bp_lagrange_points_m(*syn);
        if (lpts.has_value())
        {
            std::printf("lagrange_points\n");
            std::printf("primary,%.6e,%.6e,%.6e\n", lpts->primary_m.x, lpts->primary_m.y, lpts->primary_m.z);
            std::printf("secondary,%.6e,%.6e,%.6e\n", lpts->secondary_m.x, lpts->secondary_m.y, lpts->secondary_m.z);
            std::printf("L1,%.6e,%.6e,%.6e\n", lpts->L1_m.x, lpts->L1_m.y, lpts->L1_m.z);
            std::printf("L2,%.6e,%.6e,%.6e\n", lpts->L2_m.x, lpts->L2_m.y, lpts->L2_m.z);
            std::printf("L3,%.6e,%.6e,%.6e\n", lpts->L3_m.x, lpts->L3_m.y, lpts->L3_m.z);
            std::printf("L4,%.6e,%.6e,%.6e\n", lpts->L4_m.x, lpts->L4_m.y, lpts->L4_m.z);
            std::printf("L5,%.6e,%.6e,%.6e\n", lpts->L5_m.x, lpts->L5_m.y, lpts->L5_m.z);
        }
    }

    // Convert spacecraft trajectory to synodic frame
    // Use inertial trajectory (not relative to Earth)
    const auto traj_opt_inertial = trajectory_options()
            .duration(days(20.0))
            .sample_dt(minutes(10.0))
            .celestial_dt(minutes(5.0))
            .max_samples(100'000);
    // Note: reuse the inertial trajectory sampled earlier; do not feed a pre-transformed (relative) trajectory
    // into trajectory_to_synodic().
    const std::vector<TrajectorySample> sc_traj_synodic = trajectory_to_synodic(sc_traj_inertial, eph, *earth_ptr, *moon_ptr);

    std::printf("sc_synodic\n");
    for (const auto &s: sc_traj_synodic)
    {
        std::printf("%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s.t_s, s.position_m.x, s.position_m.y, s.position_m.z,
                    s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
    }

    // Example: using the generic frame transform API with SynodicFrame (display/analysis only).
    // Show the Earth state expressed in the synodic frame at t=0.
    if (syn.has_value())
    {
        const State earth_in_syn = inertial_state_to_frame(earth_ptr->state, *syn);
        std::printf("earth_state_in_synodic_t0\n");
        std::printf("r_m,%.6e,%.6e,%.6e\n", earth_in_syn.position_m.x, earth_in_syn.position_m.y, earth_in_syn.position_m.z);
        std::printf("v_mps,%.6e,%.6e,%.6e\n", earth_in_syn.velocity_mps.x, earth_in_syn.velocity_mps.y, earth_in_syn.velocity_mps.z);
    }

    // -------------------------------------------------------------------------
    // Event and node prediction output (for visualization/debugging)
    // -------------------------------------------------------------------------
    std::printf("\n--- predicted events/nodes ---\n");

    auto event_type_str = [](const EventType t) -> const char * {
        switch (t)
        {
        case EventType::Impact:
            return "impact";
        case EventType::AtmosphereBoundary:
            return "atmosphere";
        case EventType::SoiBoundary:
            return "soi";
        default:
            return "unknown";
        }
    };
    auto crossing_str = [](const Crossing c) -> const char * {
        switch (c)
        {
        case Crossing::Enter:
            return "enter";
        case Crossing::Exit:
            return "exit";
        default:
            return "unknown";
        }
    };
    auto node_str = [](const NodeCrossing c) -> const char * {
        switch (c)
        {
        case NodeCrossing::Ascending:
            return "AN";
        case NodeCrossing::Descending:
            return "DN";
        default:
            return "?";
        }
    };

    // Predict boundary events (impact/atmosphere/SOI) on the N-body ephemeris.
    std::printf("events_sc\n");
    const std::vector<Event> ev = predict_spacecraft_events(sim, eph, sc_id, traj_opt_inertial);
    for (const auto &e: ev)
    {
        std::printf("%.6f,%s,%u,%s\n",
                    e.t_event_s,
                    event_type_str(e.type),
                    static_cast<unsigned>(e.body_id),
                    crossing_str(e.crossing));
    }

    // Event probe outputs (small CSV blocks)
    if (sc_reentry_h.valid())
    {
        std::printf("events_reentry_sc\n");
        const std::vector<Event> ev_reentry =
                predict_spacecraft_events(sim, eph, sc_reentry_h.id, reentry_short_inertial);
        for (const auto &e: ev_reentry)
        {
            std::printf("%.6f,%s,%u,%s\n",
                        e.t_event_s,
                        event_type_str(e.type),
                        static_cast<unsigned>(e.body_id),
                        crossing_str(e.crossing));
        }

        std::printf("reentry_sc_rel_earth\n");
        for (const auto &s: sc_reentry_rel_earth)
        {
            std::printf("%.6f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        s.t_s,
                        s.position_m.x,
                        s.position_m.y,
                        s.position_m.z,
                        s.velocity_mps.x,
                        s.velocity_mps.y,
                        s.velocity_mps.z);
        }
    }

    if (sc_moon_soi_h.valid())
    {
        std::printf("events_moon_soi_sc\n");
        const std::vector<Event> ev_soi =
                predict_spacecraft_events(sim, eph, sc_moon_soi_h.id, moon_soi_short_inertial);
        for (const auto &e: ev_soi)
        {
            std::printf("%.6f,%s,%u,%s\n",
                        e.t_event_s,
                        event_type_str(e.type),
                        static_cast<unsigned>(e.body_id),
                        crossing_str(e.crossing));
        }

        std::printf("moon_soi_sc_rel_moon\n");
        for (const auto &s: sc_moon_soi_rel_moon)
        {
            std::printf("%.6f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        s.t_s,
                        s.position_m.x,
                        s.position_m.y,
                        s.position_m.z,
                        s.velocity_mps.x,
                        s.velocity_mps.y,
                        s.velocity_mps.z);
        }
    }

    // Equatorial nodes: crossings of the primary body's equator plane (spin axis normal).
    std::printf("equatorial_nodes_sc_earth\n");
    const std::vector<NodeEvent> eq_nodes = predict_equatorial_nodes(sim, eph, sc_id, earth_id, traj_opt_inertial);
    for (const auto &n: eq_nodes)
    {
        std::printf("%.6f,%s\n", n.t_event_s, node_str(n.crossing));
    }

    // Target plane nodes: crossings of the target spacecraft's current orbit plane about the same primary.
    // (In a KSP-like UI, this is the "AN/DN relative to target orbit" you show in the map view.)
    if (sc_lambert_h.valid())
    {
        std::printf("target_nodes_sc_to_lambert_earth\n");
        const std::vector<NodeEvent> tgt_nodes =
                predict_target_plane_nodes(sim, eph, sc_id, sc_lambert_h.id, earth_id, traj_opt_inertial);
        for (const auto &n: tgt_nodes)
        {
            std::printf("%.6f,%s\n", n.t_event_s, node_str(n.crossing));
        }
    }

    // Apsides markers from osculating Earth-centered elements at t=0 (relative to Earth).
    {
        const State earth0 = eph.body_state_at_by_id(earth_id, sim.time_s());
        const State sc0 = sim.spacecraft_by_id(sc_id)->state;
        const Vec3 r_rel = sc0.position_m - earth0.position_m;
        const Vec3 v_rel = sc0.velocity_mps - earth0.velocity_mps;
        const OrbitApsides aps = apsides_from_relative_state(mu_earth_m3_s2, r_rel, v_rel);
        std::printf("apsides_sc_rel_earth_t0\n");
        if (aps.valid)
        {
            std::printf("peri,%.6e,%.6e,%.6e\n", aps.periapsis_rel_m.x, aps.periapsis_rel_m.y, aps.periapsis_rel_m.z);
            if (aps.has_apoapsis)
            {
                std::printf("apo,%.6e,%.6e,%.6e\n", aps.apoapsis_rel_m.x, aps.apoapsis_rel_m.y, aps.apoapsis_rel_m.z);
            }
        }
    }

    return 0;
}
