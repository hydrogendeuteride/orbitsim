#pragma once

#include "orbitsim/game_sim.hpp"

#include <utility>
#include <vector>

namespace orbitsim
{

    class NBodySimulation
    {
    public:
        struct Config
        {
            double gravitational_constant{orbitsim::kGravitationalConstant_SI};
            double softening_length_m{0.0};
            DOPRI5Options spacecraft_integrator{};
        };

        NBodySimulation() : sim_(make_game_cfg_()) {}
        explicit NBodySimulation(Config cfg) : sim_(make_game_cfg_(std::move(cfg))) {}

        double time_s() const { return sim_.time_s(); }

        std::vector<MassiveBody> &massive_bodies() { return sim_.massive_bodies(); }
        const std::vector<MassiveBody> &massive_bodies() const { return sim_.massive_bodies(); }

        std::vector<Spacecraft> &spacecraft() { return sim_.spacecraft(); }
        const std::vector<Spacecraft> &spacecraft() const { return sim_.spacecraft(); }

        void step(const double dt_s) { sim_.step(dt_s); }

    private:
        static GameSimulation::Config make_game_cfg_()
        {
            GameSimulation::Config cfg;
            cfg.enable_events = false;
            return cfg;
        }

        static GameSimulation::Config make_game_cfg_(Config cfg_in)
        {
            GameSimulation::Config cfg;
            cfg.gravitational_constant = cfg_in.gravitational_constant;
            cfg.softening_length_m = cfg_in.softening_length_m;
            cfg.spacecraft_integrator = cfg_in.spacecraft_integrator;
            cfg.enable_events = false;
            return cfg;
        }

        GameSimulation sim_;
    };

} // namespace orbitsim
