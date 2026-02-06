#pragma once

#include "orbitsim/detail/spacecraft_propagation.hpp"
#include "orbitsim/ephemeris.hpp"
#include "orbitsim/integrators.hpp"
#include "orbitsim/spacecraft_state_cache.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace orbitsim::detail
{

    /** @brief Create ephemeris segment from start/end body states. */
    inline CelestialEphemerisSegment make_ephemeris_segment(const std::vector<State> &start,
                                                            const std::vector<State> &end, const double t0_s,
                                                            const double dt_s)
    {
        CelestialEphemerisSegment eph;
        eph.t0_s = t0_s;
        eph.dt_s = dt_s;
        eph.start = start;
        eph.end = end;
        return eph;
    }

    /** @brief Preview step for massive bodies only (for event search). */
    inline void preview_step_massive_bodies(std::vector<MassiveBody> &massive, double &t_s, const double dt_s,
                                            const double gravitational_constant, const double softening_length_m,
                                            CelestialEphemerisSegment *out_eph)
    {
        if (out_eph != nullptr)
        {
            *out_eph = {};
        }
        if (!(dt_s != 0.0) || !std::isfinite(dt_s))
        {
            return;
        }

        std::vector<State> start_states;
        snapshot_states(massive, &start_states);

        symplectic4_step(massive, dt_s, gravitational_constant, softening_length_m);

        std::vector<State> end_states;
        snapshot_states(massive, &end_states);

        CelestialEphemerisSegment eph = make_ephemeris_segment(start_states, end_states, t_s, dt_s);
        if (out_eph != nullptr)
        {
            *out_eph = eph;
        }

        t_s += dt_s;
    }

    /** @brief Propagate spacecraft through ephemeris segment with SpacecraftStateCache wrapper. */
    inline Spacecraft propagate_spacecraft_with_cache(const Spacecraft &sc0,
                                                      const std::vector<MassiveBody> &massive,
                                                      const CelestialEphemerisSegment &eph, const ManeuverPlan &plan,
                                                      const double gravitational_constant,
                                                      const double softening_length_m,
                                                      const DOPRI5Options &spacecraft_integrator, const double t0_s,
                                                      const double dt_s,
                                                      const std::vector<Spacecraft> &spacecraft_at_t0)
    {
        if (dt_s < 0.0)
        {
            return propagate_spacecraft_in_ephemeris(sc0, massive, eph, plan, gravitational_constant,
                                                     softening_length_m, spacecraft_integrator, t0_s, dt_s, nullptr);
        }

        const double t_end_s = t0_s + dt_s;
        const double lookup_dt = [&]() -> double {
            const double max_step = spacecraft_integrator.max_step_s;
            if (max_step > 0.0 && std::isfinite(max_step))
            {
                return max_step;
            }
            const double abs_dt = std::abs(dt_s);
            if (abs_dt > 0.0 && std::isfinite(abs_dt))
            {
                return std::min(abs_dt, 1.0);
            }
            return 0.0;
        }();

        SpacecraftStateCache<CelestialEphemerisSegment> sc_cache(
                massive, eph, plan, gravitational_constant, softening_length_m, spacecraft_integrator, t0_s, t_end_s,
                [&](const SpacecraftId id) -> const Spacecraft * {
                    for (const auto &sc : spacecraft_at_t0)
                    {
                        if (sc.id == id)
                        {
                            return &sc;
                        }
                    }
                    return nullptr;
                },
                SpacecraftStateCache<CelestialEphemerisSegment>::Options{.lookup_dt_s = lookup_dt});
        const SpacecraftStateLookup sc_lookup = sc_cache.lookup();

        return propagate_spacecraft_in_ephemeris(sc0, massive, eph, plan, gravitational_constant, softening_length_m,
                                                  spacecraft_integrator, t0_s, dt_s, sc_lookup);
    }

    /** @brief Execute a single timestep without event subdivision. */
    inline void step_no_events(std::vector<MassiveBody> &massive, std::vector<Spacecraft> &spacecraft, double &time_s,
                               const double dt_s, const ManeuverPlan &plan, const double gravitational_constant,
                               const double softening_length_m, const DOPRI5Options &spacecraft_integrator)
    {
        if (!(dt_s != 0.0) || !std::isfinite(dt_s))
        {
            return;
        }

        std::vector<State> start_states;
        snapshot_states(massive, &start_states);

        symplectic4_step(massive, dt_s, gravitational_constant, softening_length_m);

        std::vector<State> end_states;
        snapshot_states(massive, &end_states);

        const CelestialEphemerisSegment eph = make_ephemeris_segment(start_states, end_states, time_s, dt_s);

        const std::vector<Spacecraft> spacecraft_start = spacecraft;
        for (std::size_t sc_index = 0; sc_index < spacecraft.size(); ++sc_index)
        {
            spacecraft[sc_index] = propagate_spacecraft_with_cache(spacecraft_start[sc_index], massive, eph, plan,
                                                                    gravitational_constant, softening_length_m,
                                                                    spacecraft_integrator, time_s, dt_s,
                                                                    spacecraft_start);
        }

        time_s += dt_s;
    }

} // namespace orbitsim::detail
