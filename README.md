# orbitsim (header-only)

Small C++20 header-only orbit simulation library intended for games:

- **2-body Kepler propagator** (universal variables) for quick/debug orbits
- **Restricted full N-body** simulation:
  - **Massive bodies** (planets/moons): symplectic 4th order (Yoshida composition)
  - **Massless spacecraft**: DOPRI5(4) with interpolation of massive body states

Assumptions / conventions:

- **Units**: SI (`m`, `s`, `kg`)
- **Frame**: inertial barycentric
- **Angles**: radians
- **Rotation**: fixed spin axis + constant spin rate (no torques)

## Build the example

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/orbitsim_example
```

## Use in your project

Include headers from `include/` and GLM from `third_party/glm`.

```cpp
#include <orbitsim/orbitsim.hpp>

orbitsim::NBodySimulation sim;
sim.massive_bodies().push_back(orbitsim::MassiveBody{ .mass_kg = 1.0e30 });
sim.spacecraft().push_back(orbitsim::Spacecraft{});
sim.step(60.0);
```

