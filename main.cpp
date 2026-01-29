#include <orbitsim/orbitsim.hpp>

#include <glm/glm.hpp>

#include <cstdio>

int main()
{
    orbitsim::NBodySimulation sim;

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
            .state =
                    {
                            .position_m = -(Me / m_tot) * r_rel,
                            .velocity_mps = -(Me / m_tot) * v_rel,
                            .spin = {.axis = {0.0, 1.0, 0.0}, .angle_rad = 0.0, .rate_rad_per_s = 2.9e-6},
                    },
    };
    orbitsim::MassiveBody earth{
            .mass_kg = Me,
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
    };
    sim.spacecraft().push_back(sc);

    constexpr double dt = 60.0; // s
    constexpr int steps = 24 * 60; // 1 day

    for (int i = 0; i < steps; ++i)
    {
        sim.step(dt);
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

    return 0;
}
