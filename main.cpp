#include <orbitsim/orbitsim.hpp>

#include <cmath>
#include <cstdio>
#include <limits>

using namespace orbitsim;

// =============================================================================
// CSV Output Helpers (for Python visualization scripts)
// =============================================================================

inline void print_trajectory(const char *label, const std::vector<TrajectorySample> &traj)
{
    std::printf("%s\n", label);
    for (const auto &s : traj)
    {
        std::printf("%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                    s.t_s,
                    s.position_m.x, s.position_m.y, s.position_m.z,
                    s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
    }
}

inline void print_trajectory_hires(const char *label, const std::vector<TrajectorySample> &traj)
{
    std::printf("%s\n", label);
    for (const auto &s : traj)
    {
        std::printf("%.6f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                    s.t_s,
                    s.position_m.x, s.position_m.y, s.position_m.z,
                    s.velocity_mps.x, s.velocity_mps.y, s.velocity_mps.z);
    }
}

inline const char *event_type_str(EventType t)
{
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
}

inline const char *crossing_str(Crossing c)
{
    switch (c)
    {
    case Crossing::Enter:
        return "enter";
    case Crossing::Exit:
        return "exit";
    default:
        return "unknown";
    }
}

inline const char *node_str(NodeCrossing c)
{
    switch (c)
    {
    case NodeCrossing::Ascending:
        return "AN";
    case NodeCrossing::Descending:
        return "DN";
    default:
        return "?";
    }
}

// =============================================================================
// Physical Constants
// =============================================================================

constexpr double Ms = 1.98847e30;      // Sun mass [kg]
constexpr double Me = 5.9722e24;       // Earth mass [kg]
constexpr double Mm = 7.342e22;        // Moon mass [kg]
constexpr double AU = 1.495978707e11;  // Astronomical unit [m]
constexpr double d_em_m = 384'400'000.0;  // Earth-Moon mean distance [m]

// =============================================================================
// Section 1: Body Metadata Output
// Purpose: visualize_orbit.py - body info for drawing spheres
// =============================================================================

void output_body_metadata(const GameSimulation &sim, BodyId sun_id, BodyId earth_id, BodyId moon_id)
{
    std::printf("\n--- bodies ---\n");
    std::printf("bodies\n");

    auto print_body = [&](BodyId id, const char *name) {
        const MassiveBody *b = sim.body_by_id(id);
        if (b == nullptr) return;
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
}

// =============================================================================
// Section 2: Orbit Trajectories Output
// Purpose: visualize_orbit.py - orbital paths in various reference frames
// =============================================================================

void output_orbit_trajectories(
    const std::vector<TrajectorySample> &earth_traj,
    const std::vector<TrajectorySample> &moon_traj_sun,
    const std::vector<TrajectorySample> &sc_traj_sun,
    const std::vector<TrajectorySample> &moon_traj_earth,
    const std::vector<TrajectorySample> &sc_traj,
    const std::vector<TrajectorySample> &sc_traj_earth_hi,
    const std::vector<TrajectorySample> &sc_traj_moon)
{
    std::printf("\n--- trajectory csv (t_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps) ---\n");

    print_trajectory("earth_rel_sun", earth_traj);
    print_trajectory("moon_rel_sun", moon_traj_sun);
    print_trajectory("sc_rel_sun", sc_traj_sun);
    print_trajectory("moon_rel_earth", moon_traj_earth);
    print_trajectory("sc_rel_earth", sc_traj);
    print_trajectory("sc_rel_earth_hi", sc_traj_earth_hi);
    print_trajectory("sc_rel_moon", sc_traj_moon);
}

// =============================================================================
// Section 3: Lambert Transfer Output
// Purpose: visualize_lambert_transfer.py - Earth-Moon transfer solution
// =============================================================================

void output_lambert_section(
    const GameSimulation &sim,
    const CelestialEphemeris &eph,
    SpacecraftId sc_id,
    BodyId earth_id,
    BodyId moon_id,
    const TrajectoryFrameSpec &earth_frame)
{
    const double t_depart_s = sim.time_s();
    const double dt_transfer_s = days(9.0);
    const double t_arrive_s = t_depart_s + dt_transfer_s;

    const State earth_at_depart = eph.body_state_at_by_id(earth_id, t_depart_s);
    const State earth_at_arrive = eph.body_state_at_by_id(earth_id, t_arrive_s);
    const State moon_at_arrive = eph.body_state_at_by_id(moon_id, t_arrive_s);

    const Vec3 r1_rel_earth_m = sim.spacecraft_by_id(sc_id)->state.position_m - earth_at_depart.position_m;
    const Vec3 v0_rel_earth_mps = sim.spacecraft_by_id(sc_id)->state.velocity_mps - earth_at_depart.velocity_mps;
    const Vec3 r2_rel_earth_m = moon_at_arrive.position_m - earth_at_arrive.position_m;

    const double mu_earth_m3_s2 = kGravitationalConstant_SI * Me;

    LambertOptions lam_opt{};
    lam_opt.prograde = true;
    lam_opt.short_path = true;
    lam_opt.max_revolutions = 2;

    const std::vector<LambertSolution> lam_solutions =
        solve_lambert_universal(mu_earth_m3_s2, r1_rel_earth_m, r2_rel_earth_m, dt_transfer_s, lam_opt);

    std::optional<LambertSolution> best;
    double best_dv = std::numeric_limits<double>::infinity();
    for (const auto &sol : lam_solutions)
    {
        const double dv = glm::length(sol.v1_mps - v0_rel_earth_mps);
        if (dv < best_dv)
        {
            best_dv = dv;
            best = sol;
        }
    }

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

        // Create Lambert spacecraft for trajectory visualization
        const Spacecraft *base_sc = sim.spacecraft_by_id(sc_id);
        if (base_sc != nullptr)
        {
            // Make a copy of the simulation for Lambert trajectory
            GameSimulation sim_lam = sim;
            Spacecraft sc_lambert = *base_sc;
            sc_lambert.id = kInvalidSpacecraftId;
            sc_lambert.engines.clear();
            sc_lambert.prop_mass_kg = 0.0;
            sc_lambert.state.position_m = earth_at_depart.position_m + r1_rel_earth_m;
            sc_lambert.state.velocity_mps = earth_at_depart.velocity_mps + best->v1_mps;

            const auto sc_lambert_h = sim_lam.create_spacecraft(std::move(sc_lambert));
            if (sc_lambert_h.valid())
            {
                const auto lam_opt_traj = trajectory_options()
                    .duration(dt_transfer_s)
                    .sample_dt(minutes(30.0))
                    .celestial_dt(minutes(10.0))
                    .max_samples(50'000);

                const std::vector<TrajectorySample> lam_sc_inertial =
                    predict_spacecraft_trajectory(sim_lam, eph, sc_lambert_h.id, lam_opt_traj);
                const std::vector<TrajectorySample> lam_moon_inertial =
                    predict_body_trajectory(sim_lam, eph, moon_id, lam_opt_traj);

                const std::vector<TrajectorySample> lam_sc_rel_earth =
                    trajectory_to_frame_spec(lam_sc_inertial, eph, sim_lam.massive_bodies(), earth_frame);
                const std::vector<TrajectorySample> lam_moon_rel_earth =
                    trajectory_to_frame_spec(lam_moon_inertial, eph, sim_lam.massive_bodies(), earth_frame);

                print_trajectory("lambert_sc_rel_earth", lam_sc_rel_earth);
                print_trajectory("lambert_moon_rel_earth", lam_moon_rel_earth);
            }
        }
    }
    else
    {
        std::printf("best_z,nan\n");
        std::printf("best_dv_mps,nan\n");
    }
}

// =============================================================================
// Section 4: Synodic Frame (Rotating Frame) Output
// Purpose: visualize_synodic_frame.py - CR3BP visualization
// =============================================================================

void output_synodic_frame(
    const MassiveBody &earth,
    const MassiveBody &moon,
    const CelestialEphemeris &eph,
    const std::vector<TrajectorySample> &sc_traj_inertial,
    const std::vector<TrajectorySample> &earth_traj_inertial,
    const std::vector<TrajectorySample> &moon_traj_inertial)
{
    std::printf("\n--- rotating frame (Earth-Moon synodic) ---\n");

    const std::optional<SynodicFrame> syn = make_synodic_frame(earth, moon);
    if (!syn.has_value()) return;

    std::printf("synodic_frame_t0\n");
    std::printf("separation_m,%.6e\n", syn->separation_m);
    std::printf("mu,%.6e\n", syn->mu);
    std::printf("omega_radps,%.6e,%.6e,%.6e\n",
                syn->omega_inertial_radps.x,
                syn->omega_inertial_radps.y,
                syn->omega_inertial_radps.z);

    // Lagrange points (at t=0)
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

    // Spacecraft trajectory in synodic frame
    const std::vector<TrajectorySample> sc_traj_synodic =
        trajectory_to_synodic(sc_traj_inertial, eph, earth, moon);
    print_trajectory("sc_synodic", sc_traj_synodic);

    // Earth and Moon trajectories in synodic frame (shows variable separation)
    const std::vector<TrajectorySample> earth_traj_synodic =
        trajectory_to_synodic(earth_traj_inertial, eph, earth, moon);
    const std::vector<TrajectorySample> moon_traj_synodic =
        trajectory_to_synodic(moon_traj_inertial, eph, earth, moon);

    print_trajectory("earth_synodic", earth_traj_synodic);
    print_trajectory("moon_synodic", moon_traj_synodic);

    // Earth state in synodic frame at t=0
    const State earth_in_syn = inertial_state_to_frame(earth.state, *syn);
    std::printf("earth_state_in_synodic_t0\n");
    std::printf("r_m,%.6e,%.6e,%.6e\n",
                earth_in_syn.position_m.x, earth_in_syn.position_m.y, earth_in_syn.position_m.z);
    std::printf("v_mps,%.6e,%.6e,%.6e\n",
                earth_in_syn.velocity_mps.x, earth_in_syn.velocity_mps.y, earth_in_syn.velocity_mps.z);
}

// =============================================================================
// Section 5: Event Detection Output
// Purpose: visualize_events.py - SOI/atmosphere/impact events
// =============================================================================

void output_events_section(
    const GameSimulation &sim,
    const CelestialEphemeris &eph,
    SpacecraftId sc_id,
    BodyId earth_id,
    const TrajectoryOptions &traj_opt)
{
    std::printf("\n--- predicted events/nodes ---\n");

    const double mu_earth_m3_s2 = kGravitationalConstant_SI * Me;

    // Boundary events (impact/atmosphere/SOI)
    std::printf("events_sc\n");
    const std::vector<Event> ev = predict_spacecraft_events(sim, eph, sc_id, traj_opt);
    for (const auto &e : ev)
    {
        std::printf("%.6f,%s,%u,%s\n",
                    e.t_event_s,
                    event_type_str(e.type),
                    static_cast<unsigned>(e.body_id),
                    crossing_str(e.crossing));
    }

    // Equatorial nodes
    std::printf("equatorial_nodes_sc_earth\n");
    const std::vector<NodeEvent> eq_nodes = predict_equatorial_nodes(sim, eph, sc_id, earth_id, traj_opt);
    for (const auto &n : eq_nodes)
    {
        std::printf("%.6f,%s\n", n.t_event_s, node_str(n.crossing));
    }

    // Apsides markers
    const State earth0 = eph.body_state_at_by_id(earth_id, sim.time_s());
    const State sc0 = sim.spacecraft_by_id(sc_id)->state;
    const Vec3 r_rel = sc0.position_m - earth0.position_m;
    const Vec3 v_rel = sc0.velocity_mps - earth0.velocity_mps;
    const OrbitApsides aps = apsides_from_relative_state(mu_earth_m3_s2, r_rel, v_rel);

    std::printf("apsides_sc_rel_earth_t0\n");
    if (aps.valid)
    {
        std::printf("peri,%.6e,%.6e,%.6e\n",
                    aps.periapsis_rel_m.x, aps.periapsis_rel_m.y, aps.periapsis_rel_m.z);
        if (aps.has_apoapsis)
        {
            std::printf("apo,%.6e,%.6e,%.6e\n",
                        aps.apoapsis_rel_m.x, aps.apoapsis_rel_m.y, aps.apoapsis_rel_m.z);
        }
    }
}

// =============================================================================
// Section 6: Event Probe Spacecraft Output
// Purpose: visualize_events.py - reentry/SOI test cases
// =============================================================================

void output_event_probes(
    const GameSimulation &sim,
    const CelestialEphemeris &eph,
    const std::optional<GameSimulation::SpacecraftHandle> &sc_reentry_h,
    const std::optional<GameSimulation::SpacecraftHandle> &sc_moon_soi_h,
    const TrajectoryFrameSpec &earth_frame,
    const TrajectoryFrameSpec &moon_frame)
{
    // Reentry probe
    if (sc_reentry_h.has_value() && sc_reentry_h->valid())
    {
        const auto reentry_opt = trajectory_options()
            .duration(minutes(12.0))
            .sample_dt(seconds(1.0))
            .celestial_dt(seconds(1.0))
            .max_samples(20'000);

        std::printf("events_reentry_sc\n");
        const std::vector<Event> ev_reentry =
            predict_spacecraft_events(sim, eph, sc_reentry_h->id, reentry_opt);
        for (const auto &e : ev_reentry)
        {
            std::printf("%.6f,%s,%u,%s\n",
                        e.t_event_s,
                        event_type_str(e.type),
                        static_cast<unsigned>(e.body_id),
                        crossing_str(e.crossing));
        }

        const std::vector<TrajectorySample> sc_reentry_inertial =
            predict_spacecraft_trajectory(sim, eph, sc_reentry_h->id, reentry_opt);
        const std::vector<TrajectorySample> sc_reentry_rel_earth =
            trajectory_to_frame_spec(sc_reentry_inertial, eph, sim.massive_bodies(), earth_frame);
        print_trajectory_hires("reentry_sc_rel_earth", sc_reentry_rel_earth);
    }

    // Moon SOI probe
    if (sc_moon_soi_h.has_value() && sc_moon_soi_h->valid())
    {
        const auto moon_soi_opt = trajectory_options()
            .duration(hours(4.0))
            .sample_dt(seconds(10.0))
            .celestial_dt(seconds(10.0))
            .max_samples(20'000);

        std::printf("events_moon_soi_sc\n");
        const std::vector<Event> ev_soi =
            predict_spacecraft_events(sim, eph, sc_moon_soi_h->id, moon_soi_opt);
        for (const auto &e : ev_soi)
        {
            std::printf("%.6f,%s,%u,%s\n",
                        e.t_event_s,
                        event_type_str(e.type),
                        static_cast<unsigned>(e.body_id),
                        crossing_str(e.crossing));
        }

        const std::vector<TrajectorySample> sc_moon_soi_inertial =
            predict_spacecraft_trajectory(sim, eph, sc_moon_soi_h->id, moon_soi_opt);
        const std::vector<TrajectorySample> sc_moon_soi_rel_moon =
            trajectory_to_frame_spec(sc_moon_soi_inertial, eph, sim.massive_bodies(), moon_frame);
        print_trajectory_hires("moon_soi_sc_rel_moon", sc_moon_soi_rel_moon);
    }
}

// =============================================================================
// Section 7: Spacecraft Trajectory Comparison (ECI vs Synodic)
// Purpose: visualize_frame_comparison.py - coordinate frame demonstration
// =============================================================================

void output_frame_comparison_example(
    GameSimulation &sim,
    const CelestialEphemeris &eph,
    const MassiveBody &earth,
    const MassiveBody &moon,
    BodyId earth_id,
    BodyId moon_id)
{
    std::printf("\n--- frame comparison (ECI vs Synodic) ---\n");

    // Create two spacecraft in the Earth-Moon system:
    // 1. A spacecraft in a high elliptical orbit (will show different behavior in each frame)
    // 2. A spacecraft near L1 point (quasi-stable in synodic frame)

    // Spacecraft 1: High elliptical orbit around Earth
    // This will appear as a Keplerian ellipse in ECI, but a complex path in synodic
    {
        const double peri_alt_m = 400'000.0;   // 400 km periapsis
        const double apo_alt_m = 300'000'000.0; // ~300,000 km apoapsis (near Moon distance)
        const double inc_rad = glm::radians(15.0);

        const double rp = earth.radius_m + peri_alt_m;
        const double ra = earth.radius_m + apo_alt_m;
        const double a = (rp + ra) / 2.0;
        const double e = (ra - rp) / (ra + rp);
        const double mu_earth = kGravitationalConstant_SI * earth.mass_kg;
        const double v_peri = std::sqrt(mu_earth * (2.0 / rp - 1.0 / a));

        // Start at periapsis, velocity tangent
        const Vec3 r0_eci = Vec3{rp * std::cos(inc_rad), 0.0, rp * std::sin(inc_rad)};
        const Vec3 v0_eci = Vec3{0.0, v_peri, 0.0};

        Spacecraft sc_elliptical{
            .state = make_state(earth.state.position_m + r0_eci,
                                earth.state.velocity_mps + v0_eci),
            .dry_mass_kg = 500.0,
            .prop_mass_kg = 0.0,
        };

        const auto sc_elliptical_h = sim.create_spacecraft(std::move(sc_elliptical));

        if (sc_elliptical_h.valid())
        {
            // Output info
            std::printf("elliptical_orbit_info\n");
            std::printf("periapsis_alt_m,%.6e\n", peri_alt_m);
            std::printf("apoapsis_alt_m,%.6e\n", apo_alt_m);
            std::printf("semi_major_m,%.6e\n", a);
            std::printf("eccentricity,%.6f\n", e);
            std::printf("inclination_deg,%.2f\n", glm::degrees(inc_rad));

            // Predict trajectory (one orbital period ~ 27 days for this orbit)
            const double period_s = 2.0 * M_PI * std::sqrt(a * a * a / mu_earth);
            const auto traj_opt = trajectory_options()
                .duration(std::min(period_s * 1.2, days(30.0)))
                .sample_dt(minutes(30.0))
                .celestial_dt(minutes(10.0))
                .max_samples(100'000);

            const std::vector<TrajectorySample> sc_inertial =
                predict_spacecraft_trajectory(sim, eph, sc_elliptical_h.id, traj_opt);

            // ECI frame (Earth-centered inertial)
            const TrajectoryFrameSpec eci_frame = TrajectoryFrameSpec::body_centered_inertial(earth_id);
            const std::vector<TrajectorySample> sc_eci =
                trajectory_to_frame_spec(sc_inertial, eph, sim.massive_bodies(), eci_frame);

            // Synodic frame (Earth-Moon rotating)
            const std::vector<TrajectorySample> sc_synodic =
                trajectory_to_synodic(sc_inertial, eph, earth, moon);

            print_trajectory("elliptical_eci", sc_eci);
            print_trajectory("elliptical_synodic", sc_synodic);
        }
    }

    // Spacecraft 2: Near L1 Lagrange point
    // In synodic frame, this should stay near a fixed point
    // In ECI, it will trace a complex path
    {
        const std::optional<SynodicFrame> syn = make_synodic_frame(earth, moon);
        if (syn.has_value())
        {
            const std::optional<Cr3bpLagrangePoints> lpts = cr3bp_lagrange_points_m(*syn);
            if (lpts.has_value())
            {
                // Position slightly displaced from L1 (to show oscillation)
                const Vec3 l1_synodic = lpts->L1_m;
                const Vec3 offset_synodic = Vec3{5'000'000.0, 1'000'000.0, 500'000.0}; // 5000 km offset
                const Vec3 pos_synodic = l1_synodic + offset_synodic;

                // Convert synodic position to inertial
                const State synodic_state{.position_m = pos_synodic, .velocity_mps = Vec3{0.0}};
                const State inertial_state = frame_state_to_inertial(synodic_state, *syn);

                Spacecraft sc_l1{
                    .state = make_state(inertial_state.position_m, inertial_state.velocity_mps),
                    .dry_mass_kg = 200.0,
                    .prop_mass_kg = 0.0,
                };

                const auto sc_l1_h = sim.create_spacecraft(std::move(sc_l1));

                if (sc_l1_h.valid())
                {
                    std::printf("l1_spacecraft_info\n");
                    std::printf("l1_synodic_m,%.6e,%.6e,%.6e\n",
                                l1_synodic.x, l1_synodic.y, l1_synodic.z);
                    std::printf("offset_synodic_m,%.6e,%.6e,%.6e\n",
                                offset_synodic.x, offset_synodic.y, offset_synodic.z);

                    // Predict trajectory (30 days)
                    const auto traj_opt = trajectory_options()
                        .duration(days(30.0))
                        .sample_dt(minutes(30.0))
                        .celestial_dt(minutes(10.0))
                        .max_samples(100'000);

                    const std::vector<TrajectorySample> sc_inertial =
                        predict_spacecraft_trajectory(sim, eph, sc_l1_h.id, traj_opt);

                    // ECI frame
                    const TrajectoryFrameSpec eci_frame = TrajectoryFrameSpec::body_centered_inertial(earth_id);
                    const std::vector<TrajectorySample> sc_eci =
                        trajectory_to_frame_spec(sc_inertial, eph, sim.massive_bodies(), eci_frame);

                    // Synodic frame
                    const std::vector<TrajectorySample> sc_synodic =
                        trajectory_to_synodic(sc_inertial, eph, earth, moon);

                    print_trajectory("l1_eci", sc_eci);
                    print_trajectory("l1_synodic", sc_synodic);
                }
            }
        }
    }

    // Spacecraft 3: Circular orbit for reference
    // Shows how a simple circular orbit transforms between frames
    {
        const double alt_m = 100'000'000.0;  // 100,000 km (about 1/4 Moon distance)
        const double inc_rad = glm::radians(5.0);
        const auto circ_state = circular_orbit_relative_state(earth.mass_kg, earth.radius_m + alt_m, inc_rad);

        Spacecraft sc_circ{
            .state = make_state(earth.state.position_m + circ_state.position_m,
                                earth.state.velocity_mps + circ_state.velocity_mps),
            .dry_mass_kg = 300.0,
            .prop_mass_kg = 0.0,
        };

        const auto sc_circ_h = sim.create_spacecraft(std::move(sc_circ));

        if (sc_circ_h.valid())
        {
            const double mu_earth = kGravitationalConstant_SI * earth.mass_kg;
            const double r = earth.radius_m + alt_m;
            const double period_s = 2.0 * M_PI * std::sqrt(r * r * r / mu_earth);

            std::printf("circular_orbit_info\n");
            std::printf("altitude_m,%.6e\n", alt_m);
            std::printf("period_s,%.6f\n", period_s);
            std::printf("period_days,%.3f\n", period_s / 86400.0);

            // Predict trajectory (2 periods or 30 days max)
            const auto traj_opt = trajectory_options()
                .duration(std::min(period_s * 2.0, days(30.0)))
                .sample_dt(minutes(20.0))
                .celestial_dt(minutes(10.0))
                .max_samples(100'000);

            const std::vector<TrajectorySample> sc_inertial =
                predict_spacecraft_trajectory(sim, eph, sc_circ_h.id, traj_opt);

            // ECI frame
            const TrajectoryFrameSpec eci_frame = TrajectoryFrameSpec::body_centered_inertial(earth_id);
            const std::vector<TrajectorySample> sc_eci =
                trajectory_to_frame_spec(sc_inertial, eph, sim.massive_bodies(), eci_frame);

            // Synodic frame
            const std::vector<TrajectorySample> sc_synodic =
                trajectory_to_synodic(sc_inertial, eph, earth, moon);

            print_trajectory("circular_eci", sc_eci);
            print_trajectory("circular_synodic", sc_synodic);
        }
    }
}

// =============================================================================
// Main Entry Point
// =============================================================================

int main()
{
    // -------------------------------------------------------------------------
    // Simulation Setup
    // -------------------------------------------------------------------------

    GameSimulation::Config cfg{};
    cfg.spacecraft_integrator.max_step_s = 2.0;
    cfg.spacecraft_integrator.max_substeps = 512;
    GameSimulation sim(cfg);

    const double m_em = Me + Mm;

    // Sun <-> Earth-Moon barycenter initial conditions (tilted for 3D visualization)
    const double inc_em_rad = glm::radians(10.0);
    const auto [sun_state, em_bary_state] = two_body_circular_barycentric(Ms, m_em, AU, inc_em_rad);

    MassiveBody sun{
        .mass_kg = Ms,
        .radius_m = 6.9634e8,
        .state = make_state(sun_state.position_m, sun_state.velocity_mps, {0.0, 1.0, 0.0}, 2.9e-6),
    };

    // Earth <-> Moon circular orbit around their barycenter
    const double inc_moon_rad = glm::radians(5.0);
    const auto [earth_rel, moon_rel] = two_body_circular_barycentric(Me, Mm, d_em_m, inc_moon_rad);

    MassiveBody earth{
        .mass_kg = Me,
        .radius_m = 6.371e6,
        .atmosphere_top_height_m = 1.0e5,
        .soi_radius_m = 9.25e8,
        .state = make_state(em_bary_state.position_m + earth_rel.position_m,
                            em_bary_state.velocity_mps + earth_rel.velocity_mps,
                            {0.0, 1.0, 0.0}, 7.2921159e-5),
    };

    MassiveBody moon{
        .mass_kg = Mm,
        .radius_m = 1.7374e6,
        .soi_radius_m = 6.61e7,
        .state = make_state(em_bary_state.position_m + moon_rel.position_m,
                            em_bary_state.velocity_mps + moon_rel.velocity_mps,
                            {0.0, 1.0, 0.0}, 2.66e-6),
    };

    // -------------------------------------------------------------------------
    // Create Massive Bodies
    // -------------------------------------------------------------------------

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
    // Create Main Spacecraft (LEO with maneuvers)
    // -------------------------------------------------------------------------

    const double sc_altitude_m = 300'000.0;
    const double sc_inc_rad = glm::radians(28.5);
    const auto sc_orbit = circular_orbit_relative_state(earth_ptr->mass_kg,
                                                        earth_ptr->radius_m + sc_altitude_m,
                                                        sc_inc_rad);

    Spacecraft sc{
        .state = make_state(earth_ptr->state.position_m + sc_orbit.position_m,
                            earth_ptr->state.velocity_mps + sc_orbit.velocity_mps,
                            {0.0, 0.0, 1.0}, 0.0),
        .dry_mass_kg = 800.0,
        .prop_mass_kg = 2000.0,
        .engines = {{.max_thrust_N = 5'000.0, .isp_s = 300.0, .min_throttle_0_1 = 0.05}},
    };

    const GameSimulation::SpacecraftHandle sc_h = sim.create_spacecraft(std::move(sc));
    if (!sc_h.valid())
    {
        return 1;
    }
    const SpacecraftId sc_id = sc_h.id;

    // -------------------------------------------------------------------------
    // Create Event Probe Spacecraft
    // -------------------------------------------------------------------------

    // Reentry probe (guaranteed to cross atmosphere + impact)
    std::optional<GameSimulation::SpacecraftHandle> sc_reentry_h;
    {
        const double alt_m = 120'000.0;
        const Vec3 r0 = Vec3{earth_ptr->radius_m + alt_m, 0.0, 0.0};
        const Vec3 v0 = Vec3{-250.0, 0.0, 0.0};

        Spacecraft sc_reentry{
            .state = make_state(earth_ptr->state.position_m + r0,
                                earth_ptr->state.velocity_mps + v0),
            .dry_mass_kg = 100.0,
            .prop_mass_kg = 0.0,
        };
        auto h = sim.create_spacecraft(std::move(sc_reentry));
        if (h.valid()) sc_reentry_h = h;
    }

    // Moon SOI probe (starts just outside Moon SOI, drifts inward)
    std::optional<GameSimulation::SpacecraftHandle> sc_moon_soi_h;
    {
        const Vec3 r0 = Vec3{moon_ptr->soi_radius_m + 1.0e6, 0.0, 0.0};
        const Vec3 v0 = Vec3{-200.0, 0.0, 0.0};

        Spacecraft sc_moon_soi{
            .state = make_state(moon_ptr->state.position_m + r0,
                                moon_ptr->state.velocity_mps + v0),
            .dry_mass_kg = 100.0,
            .prop_mass_kg = 0.0,
        };
        auto h = sim.create_spacecraft(std::move(sc_moon_soi));
        if (h.valid()) sc_moon_soi_h = h;
    }

    // -------------------------------------------------------------------------
    // Maneuver Plan
    // -------------------------------------------------------------------------

    // Injection burn: prograde in RTN (Earth primary)
    sim.maneuver_plan().segments.push_back(
        burn().start(hours(2.0))
              .duration(minutes(18.0))
              .prograde()
              .primary(earth_id)
              .spacecraft(sc_id));

    // Mid-course correction: small normal component
    sim.maneuver_plan().segments.push_back(
        burn().start(days(2.0))
              .duration(minutes(5.0))
              .normal()
              .throttle(0.4)
              .primary(earth_id)
              .spacecraft(sc_id));

    // Lunar capture: retrograde (Moon primary)
    sim.maneuver_plan().segments.push_back(
        burn().start(days(4.5))
              .duration(minutes(15.0))
              .retrograde()
              .throttle(0.7)
              .primary(moon_id)
              .spacecraft(sc_id));

    // -------------------------------------------------------------------------
    // Build Ephemeris and Trajectories
    // -------------------------------------------------------------------------

    const auto traj_opt = trajectory_options()
        .duration(days(20.0))
        .sample_dt(minutes(10.0))
        .celestial_dt(minutes(5.0))
        .max_samples(100'000);

    const CelestialEphemeris eph = build_celestial_ephemeris(sim, traj_opt);

    // Predict inertial trajectories
    const std::vector<TrajectorySample> earth_traj_inertial = predict_body_trajectory(sim, eph, earth_id, traj_opt);
    const std::vector<TrajectorySample> moon_traj_inertial = predict_body_trajectory(sim, eph, moon_id, traj_opt);
    const std::vector<TrajectorySample> sc_traj_inertial = predict_spacecraft_trajectory(sim, eph, sc_id, traj_opt);

    // Frame specifications
    const TrajectoryFrameSpec sun_frame = TrajectoryFrameSpec::body_centered_inertial(sun_id);
    const TrajectoryFrameSpec earth_frame = TrajectoryFrameSpec::body_centered_inertial(earth_id);
    const TrajectoryFrameSpec moon_frame = TrajectoryFrameSpec::body_centered_inertial(moon_id);

    // Transform to Sun-centered frame
    const std::vector<TrajectorySample> earth_traj =
        trajectory_to_frame_spec(earth_traj_inertial, eph, sim.massive_bodies(), sun_frame);
    const std::vector<TrajectorySample> moon_traj_sun =
        trajectory_to_frame_spec(moon_traj_inertial, eph, sim.massive_bodies(), sun_frame);
    const std::vector<TrajectorySample> sc_traj_sun =
        trajectory_to_frame_spec(sc_traj_inertial, eph, sim.massive_bodies(), sun_frame);

    // Transform to Earth-centered frame
    const std::vector<TrajectorySample> moon_traj_earth =
        trajectory_to_frame_spec(moon_traj_inertial, eph, sim.massive_bodies(), earth_frame);
    const std::vector<TrajectorySample> sc_traj =
        trajectory_to_frame_spec(sc_traj_inertial, eph, sim.massive_bodies(), earth_frame);

    // High-resolution LEO trajectory
    const auto earth_hi = trajectory_options()
        .duration(hours(8.0))
        .sample_dt(seconds(5.0))
        .spacecraft_sample_dt(seconds(5.0))
        .celestial_dt(seconds(30.0))
        .max_samples(50'000);

    const std::vector<TrajectorySample> sc_traj_earth_hi_inertial =
        predict_spacecraft_trajectory(sim, eph, sc_id, earth_hi);
    const std::vector<TrajectorySample> sc_traj_earth_hi =
        trajectory_to_frame_spec(sc_traj_earth_hi_inertial, eph, sim.massive_bodies(), earth_frame);

    // Moon-centered spacecraft trajectory
    const std::vector<TrajectorySample> sc_traj_moon =
        trajectory_to_frame_spec(sc_traj_inertial, eph, sim.massive_bodies(), moon_frame);

    // =========================================================================
    // Output Sections (for Python visualization scripts)
    // =========================================================================

    // Section 1: Body metadata
    output_body_metadata(sim, sun_id, earth_id, moon_id);

    // Section 2: Orbit trajectories
    output_orbit_trajectories(earth_traj, moon_traj_sun, sc_traj_sun,
                              moon_traj_earth, sc_traj, sc_traj_earth_hi, sc_traj_moon);

    // Section 3: Lambert transfer
    output_lambert_section(sim, eph, sc_id, earth_id, moon_id, earth_frame);

    // Section 4: Synodic frame
    output_synodic_frame(*earth_ptr, *moon_ptr, eph, sc_traj_inertial,
                         earth_traj_inertial, moon_traj_inertial);

    // Section 5: Events
    output_events_section(sim, eph, sc_id, earth_id, traj_opt);

    // Section 6: Event probes
    output_event_probes(sim, eph, sc_reentry_h, sc_moon_soi_h, earth_frame, moon_frame);

    // Section 7: Frame comparison (ECI vs Synodic) - NEW
    output_frame_comparison_example(sim, eph, *earth_ptr, *moon_ptr, earth_id, moon_id);

    return 0;
}
