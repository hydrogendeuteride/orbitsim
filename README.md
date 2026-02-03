# orbitsim (header-only)

Small header-only orbit simulation library intended for games:

- 2-body Kepler propagator (universal variables) for quick/debug orbits
- Restricted full N-body simulation:
  - Massive bodies (planets/moons): symplectic 4th order (Yoshida composition)
  - Massless spacecraft: DOPRI5(4) with interpolation of massive body states
- Trajectory sampling helpers (future prediction) for orbit-line rendering

Assumptions / conventions:

- Units: SI (`m`, `s`, `kg`)
- Frame: inertial barycentric
- Angles: radians
- Rotation: fixed spin axis + constant spin rate (no torques)

## Usage
See [docs/usage.md](docs/USAGE.md).
