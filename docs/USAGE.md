# orbitsim Usage Guide

A header-only C++20 orbit simulation library for KSP-style games.

## Table of Contents

1. [Getting Started](#getting-started)
2. [Core Concepts](#core-concepts)
3. [Basic Types](#basic-types)
4. [Simulation Setup](#simulation-setup)
5. [Maneuver System](#maneuver-system)
6. [Event Detection](#event-detection)
7. [Trajectory Prediction](#trajectory-prediction)
8. [Coordinate Frame Transformations](#coordinate-frame-transformations)
9. [Orbital Mechanics Utilities](#orbital-mechanics-utilities)
10. [Lambert Solver](#lambert-solver)

---

## Getting Started

### Installation

This is a header-only library. Simply add the `include/` directory to your project.

```cpp
#include <orbitsim/orbitsim.hpp>
```

### Building

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

---

## Core Concepts

### Dual-Simulation Strategy

orbitsim integrates massive bodies and spacecraft differently:

- Massive Bodies (MassiveBody): 4th-order Yoshida symplectic integrator (long-term energy conservation)
- Spacecraft: Adaptive DOPRI5(4) integrator (high-precision local dynamics)

This separation allows spacecraft to have accurate trajectories without affecting massive body dynamics.

### Coordinate System

- Inertial barycentric frame
- SI units: meters, seconds, kilograms
- Angles: radians
- Vectors: `Vec3 = glm::dvec3` (double precision)

---

## Basic Types

### SpinState - Rotation State

```cpp
using namespace orbitsim;

// Create spin state
SpinState spin = make_spin(
    Vec3{0.0, 0.0, 1.0},  // Rotation axis (unit vector)
    7.29e-5,               // Rotation rate [rad/s]
    0.0                    // Initial angle [rad]
);
```

### State - Position/Velocity/Spin

```cpp
State state = make_state(
    Vec3{1.5e11, 0.0, 0.0},     // Position [m]
    Vec3{0.0, 29780.0, 0.0},    // Velocity [m/s]
    spin                         // Spin state (optional)
);
```

### MassiveBody - Celestial Body

```cpp
MassiveBody earth;
earth.mass_kg = 5.972e24;
earth.radius_m = 6.371e6;
earth.atmosphere_top_height_m = 100000.0;  // Atmosphere height
earth.terrain_max_height_m = 8848.0;       // Maximum terrain height
earth.soi_radius_m = 9.24e8;               // Sphere of influence radius
earth.state = make_state(
    Vec3{1.5e11, 0.0, 0.0},
    Vec3{0.0, 29780.0, 0.0}
);
```

### Engine

```cpp
Engine engine;
engine.max_thrust_N = 10000.0;  // Maximum thrust [N]
engine.isp_s = 320.0;            // Specific impulse [s]
engine.min_throttle_0_1 = 0.1;   // Minimum throttle [0-1]
```

### Spacecraft

```cpp
Spacecraft ship;
ship.state = make_state(
    Vec3{1.5e11 + 6.771e6, 0.0, 0.0},  // LEO position
    Vec3{0.0, 29780.0 + 7700.0, 0.0}    // LEO velocity
);
ship.dry_mass_kg = 1000.0;
ship.prop_mass_kg = 500.0;
ship.engines.push_back(engine);
```

---

## Simulation Setup

### Creating GameSimulation

```cpp
GameSimulation::Config cfg;
cfg.gravitational_constant = kGravitationalConstant_SI;  // 6.67430e-11
cfg.softening_length_m = 0.0;                             // N-body softening
cfg.enable_events = true;                                 // Enable event detection

GameSimulation sim(cfg);
```

### Adding Bodies

```cpp
// Automatic ID assignment (recommended)
auto earth_h = sim.create_body(earth);
if (!earth_h.valid()) {
    // handle error
}
BodyId earth_id = earth_h.id;

// Manual ID assignment
auto moon_h = sim.create_body_with_id(100, moon);
if (!moon_h.valid()) {
    // handle error
}
BodyId moon_id = moon_h.id;
```

### Adding Spacecraft

```cpp
auto ship_h = sim.create_spacecraft(ship);
if (!ship_h.valid()) {
    // handle error
}
SpacecraftId ship_id = ship_h.id;
```

### State Setters (Rails/Time-Jump)

For rails or time-jump workflows where states are advanced externally:

- `set_time_s(t)`: Sets the simulation clock without modifying any body/spacecraft states.
- `set_body_state(id, state)`: Replaces the inertial state of a massive body.
- `set_spacecraft_state(id, state)`: Replaces the inertial state of a spacecraft; also resets proximity tracking.

```cpp
// Set simulation time (does not modify body/spacecraft states)
sim.set_time_s(hours(10.0));

// Set body state
State new_earth_state = make_state(
    Vec3{1.5e11, 0.0, 0.0},
    Vec3{0.0, 29780.0, 0.0}
);
bool ok = sim.set_body_state(earth_id, new_earth_state);

// Set spacecraft state
State new_ship_state = make_state(
    Vec3{1.5e11 + 7e6, 0.0, 0.0},
    Vec3{0.0, 29780.0 + 7500.0, 0.0}
);
ok = sim.set_spacecraft_state(ship_id, new_ship_state);
```

### Running the Simulation

```cpp
// Simple step (no events)
sim.step(10.0);  // Advance 10 seconds

// Collect events
std::vector<Event> events;
sim.step(10.0, &events);

// Handle events
for (const Event& e : events) {
    if (e.type == EventType::SoiBoundary && e.crossing == Crossing::Enter) {
        std::cout << "SOI entry: body " << e.body_id << std::endl;
    }
}
```

### Time Utilities

```cpp
#include <orbitsim/time_utils.hpp>

double t1 = seconds(60.0);    // 60 seconds
double t2 = minutes(5.0);     // 300 seconds
double t3 = hours(2.0);       // 7200 seconds
double t4 = days(1.0);        // 86400 seconds

sim.step(hours(1.0));  // Simulate 1 hour
```

---

## Maneuver System

### RTN Coordinate Frame

Maneuvers are defined in the RTN (Radial-Tangential-Normal) frame:
- R (Radial): From central body toward spacecraft
- T (Tangential): Velocity direction (prograde)
- N (Normal): Orbital plane normal

### Direction Constants

```cpp
kPrograde     // +T: velocity direction
kRetrograde   // -T: opposite to velocity
kRadialOut    // +R: away from center
kRadialIn     // -R: toward center
kNormal       // +N: orbital plane normal
kAntiNormal   // -N: opposite to normal
```

### BurnSegment - Finite Duration Burn

```cpp
// Method 1: Direct creation
BurnSegment burn;
burn.t_start_s = hours(1.0);
burn.t_end_s = hours(1.0) + minutes(5.0);
burn.dir_rtn_unit = kPrograde;
burn.throttle_0_1 = 1.0;
burn.engine_index = 0;
burn.spacecraft_id = ship_id;  // or kAllSpacecraft

// Method 2: Factory functions
auto burn2 = prograde_burn(hours(2.0), minutes(3.0), ship_id);
auto burn3 = retrograde_burn(hours(3.0), minutes(2.0));
auto burn4 = normal_burn(hours(4.0), minutes(1.0));

// Method 3: Builder pattern
auto burn5 = burn()
    .start(hours(5.0))
    .duration(minutes(10.0))
    .prograde()
    .full_throttle()
    .spacecraft(ship_id)
    .build();

// Optional: compute the RTN basis in a chosen reference frame (default: inertial).
// Example: RTN about Earth, but with the basis computed in the Earth-Moon synodic frame.
auto burn6 = burn()
    .start(hours(5.5))
    .duration(minutes(10.0))
    .prograde()
    .primary(earth_id)
    .rtn_synodic(earth_id, moon_id)
    .spacecraft(ship_id)
    .build();
```

### ImpulseSegment - Instantaneous Delta-v

```cpp
// Direct creation
ImpulseSegment imp;
imp.t_s = hours(1.0);
imp.dv_rtn_mps = Vec3{0.0, 100.0, 0.0};  // 100 m/s prograde
imp.spacecraft_id = ship_id;

// Builder pattern
auto imp2 = impulse()
    .time(hours(2.0))
    .prograde(150.0)     // 150 m/s
    .spacecraft(ship_id)
    .build();

// Optional: compute the RTN basis in a chosen reference frame (default: inertial).
auto imp4 = impulse()
    .time(hours(2.5))
    .prograde(150.0)
    .primary(earth_id)
    .rtn_synodic(earth_id, moon_id)
    .spacecraft(ship_id)
    .build();

auto imp3 = impulse()
    .time(hours(3.0))
    .dv_rtn(Vec3{10.0, 50.0, -5.0})  // Combined direction
    .build();
```

### ManeuverPlan

```cpp
ManeuverPlan& plan = sim.maneuver_plan();

// Add burn segments
plan.segments.push_back(prograde_burn(hours(1.0), minutes(5.0)));
plan.segments.push_back(retrograde_burn(hours(3.0), minutes(2.0)));

// Add impulses
plan.impulses.push_back(
    impulse().time(hours(2.0)).prograde(100.0).build()
);

// Sort by time
sort_plan(plan);
```

---

## Event Detection

### Event Types

```cpp
enum class EventType {
    Impact,              // Surface collision
    AtmosphereBoundary,  // Atmosphere entry/exit
    SoiBoundary,         // Sphere of influence crossing
    Proximity,           // Spacecraft proximity
};

enum class Crossing {
    Enter,  // Entry
    Exit,   // Exit
};
```

### Event Options

```cpp
EventOptions opt;
opt.time_tol_s = 1e-3;             // Time tolerance
opt.dist_tol_m = 1e-2;             // Distance tolerance
opt.max_bisect_iters = 64;         // Maximum bisection iterations
opt.max_event_splits_per_step = 8; // Maximum splits per step
```

### Proximity Event Configuration

```cpp
GameSimulation::Config cfg;
cfg.proximity.enable = true;
cfg.proximity.center_spacecraft_id = main_ship_id;  // Center spacecraft
cfg.proximity.enter_radius_m = 1000.0;              // Entry radius
cfg.proximity.exit_radius_m = 1200.0;               // Exit radius (hysteresis)
```

### Event Handling Example

```cpp
std::vector<Event> events;
sim.step(hours(1.0), &events);

for (const Event& e : events) {
    switch (e.type) {
        case EventType::Impact:
            if (e.crossing == Crossing::Enter) {
                std::cout << "Impact! Spacecraft " << e.spacecraft_id
                          << " -> Body " << e.body_id << std::endl;
            }
            break;
        case EventType::AtmosphereBoundary:
            std::cout << (e.crossing == Crossing::Enter ? "Atmosphere entry" : "Atmosphere exit")
                      << " at t=" << e.t_event_s << std::endl;
            break;
        case EventType::SoiBoundary:
            std::cout << "SOI " << (e.crossing == Crossing::Enter ? "entry" : "exit")
                      << std::endl;
            break;
        case EventType::Proximity:
            std::cout << "Proximity event: " << e.spacecraft_id
                      << " <-> " << e.other_spacecraft_id << std::endl;
            break;
    }
}
```

---

## Trajectory Prediction

### Creating Ephemeris

```cpp
TrajectoryOptions opts;
opts.duration_s = days(10.0);        // Prediction duration
opts.sample_dt_s = minutes(10.0);    // Sampling interval
opts.celestial_dt_s = minutes(5.0);  // Celestial integration interval
opts.max_samples = 10000;            // Maximum samples

// Or use builder
auto opts2 = trajectory_options()
    .duration(days(20.0))
    .sample_dt(minutes(10.0))
    .celestial_dt(minutes(5.0))
    .max_samples(100000)
    .build();

// Build ephemeris
CelestialEphemeris eph = build_celestial_ephemeris(sim, opts);
```

### Predicting Body Trajectories

```cpp
// Sample future trajectory of a celestial body
std::vector<TrajectorySample> earth_traj =
    predict_body_trajectory(sim, eph, earth_id, opts);

for (const auto& sample : earth_traj) {
    // Use sample.t_s, sample.position_m, sample.velocity_mps
}

// Trajectory prediction always returns inertial samples.
auto earth_traj_inertial = predict_body_trajectory(sim, eph, earth_id, opts);

// Convert for display/analysis (e.g., Sun-centered inertial coordinates).
TrajectoryFrameSpec sun_frame = TrajectoryFrameSpec::body_centered_inertial(sun_id);
auto earth_rel_traj = trajectory_to_frame_spec(earth_traj_inertial, eph, sim.massive_bodies(), sun_frame);

// Convert to a synodic (co-rotating) frame (Earth-Moon).
TrajectoryFrameSpec em_syn = TrajectoryFrameSpec::synodic(earth_id, moon_id);
auto ship_traj_inertial = predict_spacecraft_trajectory(sim, eph, ship_id, opts);
auto ship_syn_traj = trajectory_to_frame_spec(ship_traj_inertial, eph, sim.massive_bodies(), em_syn);
```

### Predicting Spacecraft Trajectories

```cpp
std::vector<TrajectorySample> ship_traj =
    predict_spacecraft_trajectory(sim, eph, ship_id, opts);

// Convert for plotting in Earth-centered inertial coordinates
TrajectoryFrameSpec earth_frame = TrajectoryFrameSpec::body_centered_inertial(earth_id);
auto ship_rel = trajectory_to_frame_spec(ship_traj, eph, sim.massive_bodies(), earth_frame);
```

### Predicting Spacecraft Events

```cpp
std::vector<Event> future_events =
    predict_spacecraft_events(sim, eph, ship_id, opts);

for (const Event& e : future_events) {
    std::cout << "Expected event at t=" << e.t_event_s << std::endl;
}
```

### Predicting Orbital Nodes

```cpp
// Equatorial plane crossing nodes (ascending/descending)
std::vector<NodeEvent> eq_nodes =
    predict_equatorial_nodes(sim, eph, ship_id, earth_id, opts);

for (const NodeEvent& n : eq_nodes) {
    if (n.crossing == NodeCrossing::Ascending) {
        std::cout << "Ascending node at t=" << n.t_event_s << std::endl;
    } else {
        std::cout << "Descending node at t=" << n.t_event_s << std::endl;
    }
}

// Target spacecraft orbital plane crossing nodes
std::vector<NodeEvent> target_nodes =
    predict_target_plane_nodes(sim, eph, ship_id, target_id, earth_id, opts);
```

## Coordinate Frame Transformations

### RotatingFrame

```cpp
// Body-centered inertial frame (ECI)
RotatingFrame eci = make_body_centered_inertial_frame(earth);

// Body-fixed frame (ECEF)
auto ecef_opt = make_body_fixed_frame(earth);
if (ecef_opt.has_value()) {
    RotatingFrame ecef = *ecef_opt;
}

// LVLH frame (spacecraft-centered)
auto lvlh_opt = make_lvlh_frame(earth.state, ship.state);
if (lvlh_opt.has_value()) {
    RotatingFrame lvlh = *lvlh_opt;
}
```

### Coordinate Transformations

```cpp
// Inertial -> Rotating frame
Vec3 pos_frame = inertial_position_to_frame(ecef, pos_inertial);
Vec3 vec_frame = inertial_vector_to_frame(ecef, vec_inertial);
State state_frame = inertial_state_to_frame(state_in, ecef);

// Rotating -> Inertial frame
Vec3 pos_inertial = frame_position_to_inertial(ecef, pos_frame);
Vec3 vec_inertial = frame_vector_to_inertial(ecef, vec_frame);
State state_inertial = frame_state_to_inertial(state_frame, ecef);
```

### SynodicFrame - Synodic Rotating Frame

A rotating frame for two-body systems.

```cpp
// Create synodic frame (Earth-Moon)
auto synodic_opt = make_synodic_frame(earth, moon);
if (synodic_opt.has_value()) {
    SynodicFrame frame = *synodic_opt;
    double separation = frame.separation_m;  // Earth-Moon distance
    double mu = frame.mu;                    // Mass ratio
}

// Time-varying synodic frame
auto frame_t = make_synodic_frame_at(eph, earth, moon, t_s);

// Convert trajectory to synodic coordinates
auto synodic_traj = trajectory_to_synodic(inertial_traj, eph, earth, moon);
```

### Trajectory Frame Conversion (ECI/ECEF, etc.)

```cpp
// Note: input samples must be inertial (not already origin-offset).

// Body-centered inertial frame (often called "ECI" for Earth)
auto eci_traj = trajectory_to_body_centered_inertial(inertial_traj, eph, earth);
// Alias
auto eci_traj2 = trajectory_to_eci(inertial_traj, eph, earth);

// Body-fixed rotating frame (often called "ECEF" for Earth)
auto ecef_traj = trajectory_to_body_fixed(inertial_traj, eph, earth);
// Alias
auto ecef_traj2 = trajectory_to_ecef(inertial_traj, eph, earth);
```

### Lagrange Point Calculation

```cpp
auto lp_opt = cr3bp_lagrange_points_m(synodic_frame);
if (lp_opt.has_value()) {
    Cr3bpLagrangePoints lp = *lp_opt;
    Vec3 L1 = lp.L1_m;  // L1 position in synodic frame
    Vec3 L2 = lp.L2_m;
    Vec3 L3 = lp.L3_m;
    Vec3 L4 = lp.L4_m;  // Triangular Lagrange points
    Vec3 L5 = lp.L5_m;
}
```

### Geodetic Coordinates

```cpp
auto ecef_opt = make_body_fixed_frame(earth);
if (ecef_opt.has_value()) {
    // Inertial position -> Geodetic coordinates
    auto geo = geodetic_from_inertial(*ecef_opt, pos_inertial, earth.radius_m);
    if (geo.has_value()) {
        double lat = geo->latitude_rad;   // Latitude
        double lon = geo->longitude_rad;  // Longitude
        double alt = geo->altitude_m;     // Altitude
    }

    // Geodetic coordinates -> Inertial position
    GeodeticCoord coord{0.5, 2.0, 400000.0};  // lat, lon, alt
    Vec3 pos = inertial_position_from_geodetic(*ecef_opt, coord, earth.radius_m);

    // Inertial velocity of a surface-fixed point
    Vec3 surface_vel = inertial_velocity_of_fixed_point(*ecef_opt, pos);
}
```

### NED Frame

```cpp
// North-East-Down local frame
auto ned_opt = make_ned_frame(*ecef_opt, coord, earth.radius_m);
if (ned_opt.has_value()) {
    RotatingFrame ned = *ned_opt;
    // ned.ex_i = North, ned.ey_i = East, ned.ez_i = Down
}
```

---

## Orbital Mechanics Utilities

### Circular Orbit State Computation

```cpp
// Circular orbit relative state
RelativeOrbitState rel = circular_orbit_relative_state(
    earth.mass_kg,     // Central body mass
    6.771e6,           // Orbital radius
    0.5,               // Inclination [rad]
    1.0                // Argument of latitude [rad]
);

// Two-body barycentric circular orbit
TwoBodyBarycentricStates bary = two_body_circular_barycentric(
    earth.mass_kg,
    moon.mass_kg,
    3.844e8,           // Separation distance
    0.089,             // Inclination
    0.0                // Initial argument of latitude
);
```

### Keplerian Orbital Elements

```cpp
// State -> Orbital elements
double mu = kGravitationalConstant_SI * earth.mass_kg;
OrbitalElements el = orbital_elements_from_relative_state(mu, r_rel, v_rel);

// Or from absolute states
OrbitalElements el2 = orbital_elements_from_state(mu, earth.state, ship.state);

// Orbital elements structure
double a = el.semi_major_axis_m;      // Semi-major axis
double e = el.eccentricity;           // Eccentricity
double i = el.inclination_rad;        // Inclination
double raan = el.raan_rad;            // Right ascension of ascending node
double argp = el.arg_periapsis_rad;   // Argument of periapsis
double ta = el.true_anomaly_rad;      // True anomaly
double ma = el.mean_anomaly_rad;      // Mean anomaly
```

### State from Orbital Elements

```cpp
OrbitalElements el;
el.semi_major_axis_m = 7000000.0;
el.eccentricity = 0.01;
el.inclination_rad = 0.5;
el.raan_rad = 1.0;
el.arg_periapsis_rad = 0.5;
el.true_anomaly_rad = 0.0;

State rel_state = relative_state_from_orbital_elements(mu, el);
State abs_state = state_from_orbital_elements(mu, earth.state, el);
```

### Orbital Scalars

```cpp
OrbitScalars scalars = orbit_scalars_from_elements(mu, el);
double periapsis = scalars.periapsis_radius_m;  // Periapsis radius
double apoapsis = scalars.apoapsis_radius_m;    // Apoapsis radius
double period = scalars.period_s;                // Orbital period
double n = scalars.mean_motion_radps;           // Mean motion
```

### Apsides Calculation

```cpp
OrbitApsides apsides = apsides_from_relative_state(mu, r_rel, v_rel);
if (apsides.valid) {
    Vec3 pe_pos = apsides.periapsis_rel_m;       // Periapsis position
    double pe_r = apsides.periapsis_radius_m;    // Periapsis radius

    if (apsides.has_apoapsis) {
        Vec3 ap_pos = apsides.apoapsis_rel_m;    // Apoapsis position
        double ap_r = apsides.apoapsis_radius_m; // Apoapsis radius
    }
}
```

### Kepler Propagator

Used for fast two-body orbit propagation.

```cpp
KeplerOptions opt;
opt.max_iterations = 64;
opt.abs_tolerance = 1e-12;

KeplerStepResult result = propagate_kepler_universal(
    mu,          // Gravitational parameter
    r0,          // Initial position
    v0,          // Initial velocity
    3600.0,      // Time interval [s]
    opt
);

if (result.converged) {
    Vec3 r_new = result.position_m;
    Vec3 v_new = result.velocity_mps;
}
```

### SOI-Based Primary Selection (Rails Mode)

For rails/patched-conics timewarp, select the appropriate primary body based on SOI:

```cpp
#include <orbitsim/soi.hpp>

// Configure SOI switching options
SoiSwitchOptions soi_opt;
soi_opt.enter_scale = 1.0;           // Enter SOI at 1.0x radius
soi_opt.exit_scale = 1.02;           // Exit SOI at 1.02x radius (hysteresis)
soi_opt.prefer_smallest_soi = true;  // Prefer most local body (moon over planet)
soi_opt.fallback_to_max_accel = true; // Fall back to max gravitational acceleration

// Select primary body for a spacecraft
BodyId current_primary = earth_id;
BodyId new_primary = select_primary_body_id_rails(
    sim, eph, ship_id, t_s, current_primary, soi_opt
);

// Or with position directly
BodyId primary = select_primary_body_id_rails(
    sim, eph, sc_pos_m, t_s, current_primary, soi_opt
);
```

Behavior:
- Immediately switches into smaller/local SOIs (e.g., Earth → Moon)
- Uses hysteresis when exiting to prevent boundary thrashing
- Falls back to max-acceleration selection when no SOI contains the spacecraft

---

## Lambert Solver

Computes the required velocities given two positions and time-of-flight.

### Basic Usage

```cpp
LambertOptions opt;
opt.prograde = true;       // Prograde orbit
opt.short_path = true;     // Short path
opt.max_revolutions = 0;   // Zero-revolution solutions only
opt.max_bisect_iters = 96;
opt.time_abs_tolerance_s = 1e-9;

std::vector<LambertSolution> solutions = solve_lambert_universal(
    mu,          // Gravitational parameter
    r1,          // Start position
    r2,          // Target position
    tof,         // Time of flight [s]
    opt
);

for (const auto& sol : solutions) {
    if (sol.converged) {
        Vec3 v1 = sol.v1_mps;  // Departure velocity
        Vec3 v2 = sol.v2_mps;  // Arrival velocity

        // Compute required delta-v
        Vec3 dv_departure = v1 - current_velocity;
        Vec3 dv_arrival = target_velocity - v2;
    }
}
```

### Multi-Revolution Solutions

```cpp
LambertOptions opt;
opt.max_revolutions = 3;  // Up to 3 revolutions

auto solutions = solve_lambert_universal(mu, r1, r2, tof, opt);
// Each revolution count may have 0-2 solutions
```

---

## Example: Complete Simulation

```cpp
#include <orbitsim/orbitsim.hpp>
#include <iostream>

using namespace orbitsim;

int main() {
    // Simulation configuration
    GameSimulation::Config cfg;
    cfg.enable_events = true;
    GameSimulation sim(cfg);

    // Add Earth
    MassiveBody earth;
    earth.mass_kg = 5.972e24;
    earth.radius_m = 6.371e6;
    earth.atmosphere_top_height_m = 100000.0;
    earth.soi_radius_m = 9.24e8;
    earth.state = make_state(Vec3{0, 0, 0}, Vec3{0, 0, 0});
    auto earth_h = sim.create_body(earth);
    if (!earth_h.valid()) {
        return 1;
    }
    BodyId earth_id = earth_h.id;

    // Add spacecraft (LEO)
    Spacecraft ship;
    double orbit_r = earth.radius_m + 400000.0;
    double v_circ = std::sqrt(kGravitationalConstant_SI * earth.mass_kg / orbit_r);
    ship.state = make_state(Vec3{orbit_r, 0, 0}, Vec3{0, v_circ, 0});
    ship.dry_mass_kg = 1000.0;
    ship.prop_mass_kg = 500.0;
    ship.engines.push_back(Engine{10000.0, 320.0, 0.1});
    auto ship_h = sim.create_spacecraft(ship);
    if (!ship_h.valid()) {
        return 1;
    }
    SpacecraftId ship_id = ship_h.id;

    // Add prograde burn
    sim.maneuver_plan().segments.push_back(
        prograde_burn(minutes(30.0), minutes(2.0), ship_id)
    );

    // Predict trajectory
    auto opts = trajectory_options()
        .duration(hours(2.0))
        .sample_dt(seconds(30.0))
        .build();

    CelestialEphemeris eph = build_celestial_ephemeris(sim, opts);
    auto traj_inertial = predict_spacecraft_trajectory(sim, eph, ship_id, opts);

    TrajectoryFrameSpec earth_frame = TrajectoryFrameSpec::body_centered_inertial(earth_id);
    auto traj = trajectory_to_frame_spec(traj_inertial, eph, sim.massive_bodies(), earth_frame);

    std::cout << "Trajectory samples: " << traj.size() << std::endl;

    // Run simulation
    std::vector<Event> events;
    while (sim.time_s() < hours(2.0)) {
        sim.step(seconds(10.0), &events);

        for (const Event& e : events) {
            std::cout << "Event at t=" << e.t_event_s << "s" << std::endl;
        }
        events.clear();
    }

    // Output final state
    const Spacecraft* final_ship = sim.spacecraft_by_id(ship_id);
    if (final_ship) {
        std::cout << "Final position: "
                  << final_ship->state.position_m.x << ", "
                  << final_ship->state.position_m.y << ", "
                  << final_ship->state.position_m.z << std::endl;
    }

    return 0;
}
```

---

## API Reference Summary

| Header | Key Features |
|--------|--------------|
| `types.hpp` | Basic types (Vec3, State, MassiveBody, Spacecraft) |
| `game_sim.hpp` | GameSimulation class |
| `maneuvers.hpp` | BurnSegment, ImpulseSegment, ManeuverPlan |
| `events.hpp` | Event, EventType, EventOptions |
| `trajectory_types.hpp` | TrajectorySample types |
| `trajectories.hpp` | Trajectory prediction functions |
| `trajectory_transforms.hpp` | Trajectory/frame conversion (ECI/ECEF/Synodic, etc.) |
| `ephemeris.hpp` | CelestialEphemeris, CelestialEphemerisSegment |
| `coordinate_frames.hpp` | RotatingFrame, coordinate transformations |
| `synodic.hpp` | SynodicFrame, Lagrange points |
| `geodesy.hpp` | GeodeticCoord, NED frame |
| `math.hpp` | Math utilities, OrbitalElements |
| `kepler.hpp` | Kepler propagator |
| `lambert.hpp` | Lambert solver |
| `orbit_utils.hpp` | Circular orbit utilities |
| `nodes.hpp` | NodeEvent, NodeCrossing |
| `soi.hpp` | Rails primary selection (SOI + hysteresis) |
| `integrators.hpp` | Integrators (symplectic4, DOPRI5) |
| `time_utils.hpp` | Time conversions (seconds, minutes, hours, days) |
