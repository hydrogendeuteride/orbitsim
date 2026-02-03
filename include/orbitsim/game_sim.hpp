#pragma once

#include "orbitsim/detail/spacecraft_propagation.hpp"
#include "orbitsim/ephemeris.hpp"
#include "orbitsim/events.hpp"
#include "orbitsim/integrators.hpp"
#include "orbitsim/maneuvers.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace orbitsim
{

    /**
     * @brief Main simulation class managing massive bodies and spacecraft.
     *
     * Uses a dual-simulation strategy:
     * - Massive bodies: 4th-order symplectic integrator (energy-conserving)
     * - Spacecraft: adaptive DOPRI5 integrator (high precision)
     *
     * Supports event detection (SOI crossings, impacts, etc.) with automatic
     * timestep subdivision, maneuver planning, and proximity tracking.
     */
    class GameSimulation
    {
    public:
        /** @brief Opaque handle returned when creating a massive body. */
        struct BodyHandle
        {
            BodyId id{kInvalidBodyId};

            inline bool valid() const { return id != kInvalidBodyId; }
            inline operator BodyId() const { return id; }
        };

        /** @brief Opaque handle returned when creating a spacecraft. */
        struct SpacecraftHandle
        {
            SpacecraftId id{kInvalidSpacecraftId};

            inline bool valid() const { return id != kInvalidSpacecraftId; }
            inline operator SpacecraftId() const { return id; }
        };

        /** @brief Simulation configuration parameters. */
        struct Config
        {
            double gravitational_constant{orbitsim::kGravitationalConstant_SI}; ///< G constant (m^3/kg/s^2)
            double softening_length_m{0.0}; ///< Softening to avoid singularities
            DOPRI5Options spacecraft_integrator{}; ///< Spacecraft integrator settings
            EventOptions events{}; ///< Event detection settings
            bool enable_events{true}; ///< Enable event-aware timestepping

            /** @brief Options for spacecraft-to-spacecraft proximity detection. */
            struct ProximityOptions
            {
                bool enable{false};
                SpacecraftId center_spacecraft_id{kInvalidSpacecraftId}; ///< Reference spacecraft
                double enter_radius_m{0.0}; ///< Radius to trigger "enter" event
                double exit_radius_m{0.0}; ///< Radius to trigger "exit" (use >= enter for hysteresis)
            };

            ProximityOptions proximity{};
        };

        GameSimulation() = default;
        explicit GameSimulation(Config cfg) : cfg_(std::move(cfg)) {}

        double time_s() const { return time_s_; }
        const Config &config() const { return cfg_; }

        const std::vector<MassiveBody> &massive_bodies() const { return massive_; }

        BodyId allocate_body_id()
        {
            BodyId id = next_body_id_++;
            if (id == kInvalidBodyId)
            {
                id = next_body_id_++;
            }
            while (body_id_to_index_.contains(id) || id == kInvalidBodyId)
            {
                id = next_body_id_++;
            }
            return id;
        }

        /** @brief Add a massive body to the simulation. Auto-assigns ID if invalid. */
        BodyHandle create_body(MassiveBody body)
        {
            if (body.id == kInvalidBodyId)
            {
                body.id = allocate_body_id();
            }
            if (body_id_to_index_.contains(body.id))
            {
                return BodyHandle{};
            }

            body_id_to_index_[body.id] = massive_.size();
            massive_.push_back(std::move(body));
            return BodyHandle{.id = massive_.back().id};
        }

        BodyHandle create_body_with_id(const BodyId id, MassiveBody body)
        {
            if (id == kInvalidBodyId)
            {
                return BodyHandle{};
            }
            body.id = id;
            return create_body(std::move(body));
        }

        MassiveBody *body_by_id(const BodyId id)
        {
            auto it = body_id_to_index_.find(id);
            if (it == body_id_to_index_.end())
            {
                return nullptr;
            }
            return &massive_[it->second];
        }

        const MassiveBody *body_by_id(const BodyId id) const
        {
            auto it = body_id_to_index_.find(id);
            if (it == body_id_to_index_.end())
            {
                return nullptr;
            }
            return &massive_[it->second];
        }

        bool has_body(const BodyId id) const { return body_id_to_index_.contains(id); }

        /** @brief Remove a massive body by ID. Returns false if not found. */
        bool remove_body(const BodyId id)
        {
            auto it = body_id_to_index_.find(id);
            if (it == body_id_to_index_.end())
            {
                return false;
            }

            const std::size_t index = it->second;
            const std::size_t last = massive_.empty() ? 0 : (massive_.size() - 1);
            if (index != last)
            {
                massive_[index] = std::move(massive_[last]);
                body_id_to_index_[massive_[index].id] = index;
            }

            massive_.pop_back();
            body_id_to_index_.erase(it);
            return true;
        }

        const std::vector<Spacecraft> &spacecraft() const { return spacecraft_; }

        SpacecraftId allocate_spacecraft_id()
        {
            SpacecraftId id = next_spacecraft_id_++;
            if (id == kInvalidSpacecraftId || id == kAllSpacecraft)
            {
                id = next_spacecraft_id_++;
            }
            while (spacecraft_id_to_index_.contains(id) || id == kInvalidSpacecraftId || id == kAllSpacecraft)
            {
                id = next_spacecraft_id_++;
            }
            return id;
        }

        /** @brief Add a spacecraft to the simulation. Auto-assigns ID if invalid. */
        SpacecraftHandle create_spacecraft(Spacecraft sc)
        {
            if (sc.id == kInvalidSpacecraftId || sc.id == kAllSpacecraft)
            {
                sc.id = allocate_spacecraft_id();
            }
            if (spacecraft_id_to_index_.contains(sc.id))
            {
                return SpacecraftHandle{};
            }

            spacecraft_id_to_index_[sc.id] = spacecraft_.size();
            spacecraft_.push_back(std::move(sc));
            return SpacecraftHandle{.id = spacecraft_.back().id};
        }

        SpacecraftHandle create_spacecraft_with_id(const SpacecraftId id, Spacecraft sc)
        {
            if (id == kInvalidSpacecraftId || id == kAllSpacecraft)
            {
                return SpacecraftHandle{};
            }
            sc.id = id;
            return create_spacecraft(std::move(sc));
        }

        bool has_spacecraft(const SpacecraftId id) const { return spacecraft_id_to_index_.contains(id); }

        Spacecraft *spacecraft_by_id(const SpacecraftId id)
        {
            auto it = spacecraft_id_to_index_.find(id);
            if (it == spacecraft_id_to_index_.end())
            {
                return nullptr;
            }
            return &spacecraft_[it->second];
        }

        const Spacecraft *spacecraft_by_id(const SpacecraftId id) const
        {
            auto it = spacecraft_id_to_index_.find(id);
            if (it == spacecraft_id_to_index_.end())
            {
                return nullptr;
            }
            return &spacecraft_[it->second];
        }

        /** @brief Remove a spacecraft by ID. Returns false if not found. */
        bool remove_spacecraft(const SpacecraftId id)
        {
            auto it = spacecraft_id_to_index_.find(id);
            if (it == spacecraft_id_to_index_.end())
            {
                return false;
            }

            const std::size_t index = it->second;
            const std::size_t last = spacecraft_.empty() ? 0 : (spacecraft_.size() - 1);
            if (index != last)
            {
                spacecraft_[index] = std::move(spacecraft_[last]);
                spacecraft_id_to_index_[spacecraft_[index].id] = index;
            }

            spacecraft_.pop_back();
            spacecraft_id_to_index_.erase(it);
            return true;
        }

        ManeuverPlan &maneuver_plan() { return plan_; }
        const ManeuverPlan &maneuver_plan() const { return plan_; }

        /**
         * @brief Find the massive body exerting the strongest gravity on a spacecraft.
         * @return Index into massive_bodies() of the dominant body.
         */
        std::size_t select_primary_by_max_accel(const Spacecraft &sc) const
        {
            if (massive_.empty())
            {
                return 0;
            }
            const Vec3 p = sc.state.position_m;
            const double eps2 = cfg_.softening_length_m * cfg_.softening_length_m;

            std::size_t best = 0;
            double best_a = -1.0;
            for (std::size_t i = 0; i < massive_.size(); ++i)
            {
                const Vec3 dr = massive_[i].state.position_m - p;
                const double r2 = glm::dot(dr, dr) + eps2;
                if (!(r2 > 0.0) || !std::isfinite(r2))
                {
                    continue;
                }
                const double amag = (cfg_.gravitational_constant * massive_[i].mass_kg) / r2;
                if (amag > best_a)
                {
                    best_a = amag;
                    best = i;
                }
            }
            return best;
        }

        /** @brief Advance simulation by dt_s seconds. */
        void step(const double dt_s) { step(dt_s, nullptr); }

        /**
         * @brief Advance simulation with event detection.
         *
         * When events are enabled, the timestep is automatically subdivided
         * at event boundaries (SOI crossings, impacts, periapsis, etc.).
         * Detected events are appended to out_events if provided.
         *
         * @param dt_s Time to advance (positive = forward, negative = backward without events)
         * @param out_events Optional output vector for detected events
         */
        void step(const double dt_s, std::vector<Event> *out_events)
        {
            if (out_events != nullptr)
            {
                out_events->clear();
            }
            if (!(dt_s != 0.0) || !std::isfinite(dt_s))
            {
                return;
            }
            if (!cfg_.enable_events || cfg_.events.max_event_splits_per_step <= 0)
            {
                do_step_no_events_(dt_s);
                return;
            }
            if (dt_s < 0.0)
            {
                // Backwards integration with event splitting is not supported.
                do_step_no_events_(dt_s);
                return;
            }

            sort_plan(plan_);

            const bool proximity_enabled =
                    cfg_.proximity.enable && (cfg_.proximity.center_spacecraft_id != kInvalidSpacecraftId) &&
                    (cfg_.proximity.enter_radius_m > 0.0) && std::isfinite(cfg_.proximity.enter_radius_m) &&
                    (cfg_.proximity.exit_radius_m > 0.0) && std::isfinite(cfg_.proximity.exit_radius_m);

            const double proximity_enter_m = proximity_enabled ? cfg_.proximity.enter_radius_m : 0.0;
            const double proximity_exit_m =
                    proximity_enabled ? std::max(cfg_.proximity.enter_radius_m, cfg_.proximity.exit_radius_m) : 0.0;

            // Initialize/refresh proximity active set when enabled or center changes.
            if (proximity_enabled)
            {
                if (proximity_center_id_ != cfg_.proximity.center_spacecraft_id)
                {
                    proximity_center_id_ = cfg_.proximity.center_spacecraft_id;
                    proximity_active_.clear();
                    proximity_initialized_ = false;
                }

                const Spacecraft *center_ptr = spacecraft_by_id(proximity_center_id_);
                if (center_ptr == nullptr)
                {
                    proximity_active_.clear();
                    proximity_initialized_ = false;
                }
                else if (!proximity_initialized_)
                {
                    // Treat anything within the exit radius as active at initialization.
                    for (const auto &sc: spacecraft_)
                    {
                        if (sc.id == proximity_center_id_)
                        {
                            continue;
                        }
                        const double d = glm::length(sc.state.position_m - center_ptr->state.position_m);
                        if (d <= proximity_exit_m)
                        {
                            proximity_active_.insert(sc.id);
                        }
                    }
                    proximity_initialized_ = true;
                }
            }
            else
            {
                proximity_center_id_ = kInvalidSpacecraftId;
                proximity_active_.clear();
                proximity_initialized_ = false;
            }

            double remaining = dt_s;
            int splits = 0;
            while (remaining > 0.0 && splits++ < std::max(1, cfg_.events.max_event_splits_per_step))
            {
                // Preview a full remaining step to build a single ephemeris segment for event search.
                std::vector<MassiveBody> massive_preview = massive_;
                std::vector<Spacecraft> spacecraft_preview = spacecraft_;
                double t_preview = time_s_;

                CelestialEphemerisSegment eph_preview{};
                preview_step_no_events_(massive_preview, spacecraft_preview, t_preview, remaining, &eph_preview);

                std::optional<Event> best;

                // Boundary events (impact/atmosphere/SOI) for all spacecraft.
                for (const auto &sc: spacecraft_)
                {
                    auto propagate_sc = [&](const Spacecraft &sc_start, const double t0_s,
                                            const double dt_sc_s) -> Spacecraft {
                        return propagate_spacecraft_(sc_start, eph_preview, t0_s, dt_sc_s);
                    };
                    const std::optional<Event> e = find_earliest_event_in_interval(
                            massive_, eph_preview, sc, time_s_, remaining, plan_, cfg_.events, propagate_sc);
                    if (e.has_value() && (!best.has_value() || e->t_event_s < best->t_event_s))
                    {
                        best = e;
                    }
                }

                // Proximity events: one center spacecraft vs all targets.
                if (proximity_enabled && proximity_center_id_ != kInvalidSpacecraftId)
                {
                    const Spacecraft *center_ptr = spacecraft_by_id(proximity_center_id_);
                    if (center_ptr != nullptr)
                    {
                        const Spacecraft center0 = *center_ptr;
                        auto propagate_sc = [&](const Spacecraft &sc_start, const double t0_s,
                                                const double dt_sc_s) -> Spacecraft {
                            return propagate_spacecraft_(sc_start, eph_preview, t0_s, dt_sc_s);
                        };

                        for (const auto &target: spacecraft_)
                        {
                            if (target.id == proximity_center_id_)
                            {
                                continue;
                            }

                            const bool active = proximity_active_.contains(target.id);
                            const double threshold = active ? proximity_exit_m : proximity_enter_m;
                            const std::optional<Event> e = find_earliest_proximity_event_in_interval(
                                    eph_preview, center0, target, time_s_, remaining, threshold, cfg_.events,
                                    propagate_sc);
                            if (e.has_value() && (!best.has_value() || e->t_event_s < best->t_event_s))
                            {
                                best = e;
                            }
                        }
                    }
                }

                if (!best.has_value())
                {
                    // No events found: apply preview result for massive bodies to avoid recomputation,
                    // but still need to propagate spacecraft since preview_step_no_events_ skips them.
                    massive_ = std::move(massive_preview);

                    for (std::size_t sc_index = 0; sc_index < spacecraft_.size(); ++sc_index)
                    {
                        spacecraft_[sc_index] =
                                propagate_spacecraft_(spacecraft_[sc_index], eph_preview, time_s_, remaining);
                    }

                    time_s_ = t_preview;
                    return;
                }

                double dt_event = best->t_event_s - time_s_;
                const double min_step = std::max(0.0, cfg_.events.time_tol_s);
                if (!(dt_event > min_step) || !std::isfinite(dt_event))
                {
                    dt_event = min_step;
                }
                if (dt_event > remaining)
                {
                    dt_event = remaining;
                }

                if (out_events != nullptr)
                {
                    Event e = *best;
                    e.t_event_s = time_s_ + dt_event;
                    out_events->push_back(e);
                }

                do_step_no_events_(dt_event);

                if (best->type == EventType::Proximity)
                {
                    if (best->crossing == Crossing::Enter)
                    {
                        proximity_active_.insert(best->spacecraft_id);
                    }
                    else
                    {
                        proximity_active_.erase(best->spacecraft_id);
                    }
                }
                remaining -= dt_event;
            }

            if (remaining > 0.0)
            {
                do_step_no_events_(remaining);
            }
        }

    private:
        /** @brief Create ephemeris segment from start/end body states. */
        inline CelestialEphemerisSegment make_segment_(const std::vector<State> &start, const std::vector<State> &end,
                                                       const double t0_s, const double dt_s) const
        {
            CelestialEphemerisSegment eph;
            eph.t0_s = t0_s;
            eph.dt_s = dt_s;
            eph.start = start;
            eph.end = end;
            return eph;
        }

        /** @brief Propagate spacecraft through ephemeris segment with maneuvers. */
        inline Spacecraft propagate_spacecraft_(const Spacecraft &sc0, const CelestialEphemerisSegment &eph,
                                                const double t0_s, const double dt_s) const
        {
            return detail::propagate_spacecraft_in_ephemeris(sc0, massive_, eph, plan_, cfg_.gravitational_constant,
                                                             cfg_.softening_length_m, cfg_.spacecraft_integrator, t0_s,
                                                             dt_s);
        }

        /** @brief Preview step for massive bodies only (for event search). */
        inline void preview_step_no_events_(std::vector<MassiveBody> &massive, std::vector<Spacecraft> &spacecraft,
                                            double &t_s, const double dt_s, CelestialEphemerisSegment *out_eph) const
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
            detail::snapshot_states(massive, &start_states);

            symplectic4_step(massive, dt_s, cfg_.gravitational_constant, cfg_.softening_length_m);

            std::vector<State> end_states;
            detail::snapshot_states(massive, &end_states);

            CelestialEphemerisSegment eph = make_segment_(start_states, end_states, t_s, dt_s);
            if (out_eph != nullptr)
            {
                *out_eph = eph;
            }

            // Note: spacecraft preview ignores the maneuver plan mass/prop changes in the output container; event
            // search re-propagates from the true start state via propagate_spacecraft_ anyway.
            (void) spacecraft;
            t_s += dt_s;
        }

        /** @brief Execute a single timestep without event subdivision. */
        inline void do_step_no_events_(const double dt_s)
        {
            if (!(dt_s != 0.0) || !std::isfinite(dt_s))
            {
                return;
            }

            std::vector<State> start_states;
            detail::snapshot_states(massive_, &start_states);

            symplectic4_step(massive_, dt_s, cfg_.gravitational_constant, cfg_.softening_length_m);

            std::vector<State> end_states;
            detail::snapshot_states(massive_, &end_states);

            const CelestialEphemerisSegment eph = make_segment_(start_states, end_states, time_s_, dt_s);

            for (std::size_t sc_index = 0; sc_index < spacecraft_.size(); ++sc_index)
            {
                spacecraft_[sc_index] = propagate_spacecraft_(spacecraft_[sc_index], eph, time_s_, dt_s);
            }

            time_s_ += dt_s;
        }

        Config cfg_{};
        double time_s_{0.0}; ///< Current simulation time
        std::vector<MassiveBody> massive_{}; ///< All massive bodies
        std::vector<Spacecraft> spacecraft_{}; ///< All spacecraft
        BodyId next_body_id_{1};
        std::unordered_map<BodyId, std::size_t> body_id_to_index_{};
        SpacecraftId next_spacecraft_id_{1};
        std::unordered_map<SpacecraftId, std::size_t> spacecraft_id_to_index_{};
        ManeuverPlan plan_{}; ///< Scheduled maneuvers

        // Proximity tracking state
        SpacecraftId proximity_center_id_{kInvalidSpacecraftId};
        std::unordered_set<SpacecraftId> proximity_active_{}; ///< Spacecraft currently within exit radius
        bool proximity_initialized_{false};
    };

} // namespace orbitsim
