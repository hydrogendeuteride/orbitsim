#pragma once

#include "orbitsim/detail/proximity_tracking.hpp"
#include "orbitsim/detail/simulation_stepping.hpp"
#include "orbitsim/detail/spacecraft_propagation.hpp"
#include "orbitsim/ephemeris.hpp"
#include "orbitsim/events.hpp"
#include "orbitsim/integrators.hpp"
#include "orbitsim/maneuvers.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/spacecraft_state_cache.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
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
        /**
         * @brief Set the simulation time (seconds).
         *
         * This only updates the internal clock; it does not modify body/spacecraft states.
         * Intended for "rails"/time-jump workflows where the caller has already advanced
         * states externally (e.g., via ephemeris/Kepler propagation) and needs to keep
         * the simulation clock consistent.
         *
         * Resets proximity-tracking state to avoid stale enter/exit events after a jump.
         *
         * @return true if the time was finite and applied.
         */
        bool set_time_s(const double t_s)
        {
            if (!std::isfinite(t_s))
            {
                return false;
            }
            time_s_ = t_s;
            proximity_initialized_ = false;
            proximity_active_.clear();
            return true;
        }

        const Config &config() const { return cfg_; }

        const std::vector<MassiveBody> &massive_bodies() const { return massive_; }

        /**
         * @brief Replace the inertial state of a massive body by ID.
         *
         * Intended for rails/time-jump workflows where the caller advances states externally and then injects the
         * results back into the simulation.
         *
         * @return true if the body exists and the state was finite/applied.
         */
        bool set_body_state(const BodyId id, const State &state)
        {
            MassiveBody *b = body_by_id(id);
            if (b == nullptr)
            {
                return false;
            }
            if (!is_finite_state_(state))
            {
                return false;
            }
            b->state = state;
            return true;
        }

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
         * @brief Replace the inertial state of a spacecraft by ID.
         *
         * Intended for rails/time-jump workflows where the caller advances states externally and then injects the
         * results back into the simulation.
         *
         * Resets proximity-tracking state to avoid stale enter/exit events after a jump.
         *
         * @return true if the spacecraft exists and the state was finite/applied.
         */
        bool set_spacecraft_state(const SpacecraftId id, const State &state)
        {
            Spacecraft *sc = spacecraft_by_id(id);
            if (sc == nullptr)
            {
                return false;
            }
            if (!is_finite_state_(state))
            {
                return false;
            }

            sc->state = state;
            proximity_initialized_ = false;
            proximity_active_.clear();
            return true;
        }

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
         * @param dt_s Time to advance (positive = forward, negative = backward)
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
                step_no_events_(dt_s);
                if (dt_s < 0.0)
                {
                    reset_proximity_tracking_();
                }
                return;
            }

            if (dt_s < 0.0)
            {
                // Rewind state first, then replay forward in a temporary copy to recover
                // backwards-visible events with inverted crossing directions.
                GameSimulation rewound = *this;
                rewound.step_no_events_(dt_s);

                if (out_events != nullptr)
                {
                    std::vector<Event> forward_events;
                    GameSimulation replay = rewound;
                    replay.step(-dt_s, &forward_events);
                    out_events->reserve(forward_events.size());
                    for (auto it = forward_events.rbegin(); it != forward_events.rend(); ++it)
                    {
                        Event e = *it;
                        e.crossing = invert_crossing_(e.crossing);
                        out_events->push_back(e);
                    }
                }

                *this = std::move(rewound);
                reset_proximity_tracking_();
                return;
            }

            sort_plan(plan_);

            const detail::ProximityConfig prox = resolve_proximity_config_();
            prepare_proximity_tracking_(prox);

            double remaining = dt_s;
            int splits = 0;
            const int max_splits = std::max(1, cfg_.events.max_event_splits_per_step);
            while (remaining > 0.0 && splits++ < max_splits)
            {
                StepPreview preview = preview_and_find_events_(remaining, prox);
                if (!preview.best_event.has_value())
                {
                    // No events found: apply preview result for massive bodies to avoid recomputation,
                    // but still need to propagate spacecraft since preview_step_massive_bodies skips them.
                    massive_ = std::move(preview.massive_preview);
                    spacecraft_ = std::move(preview.spacecraft_end);
                    time_s_ = preview.t_preview_s;
                    return;
                }

                const double dt_event = event_step_dt_(*preview.best_event, remaining);

                if (out_events != nullptr)
                {
                    Event e = *preview.best_event;
                    e.t_event_s = time_s_ + dt_event;
                    out_events->push_back(e);
                }

                step_no_events_(dt_event);

                if (preview.best_event->type == EventType::Proximity)
                {
                    detail::update_proximity_active_set(preview.best_event->spacecraft_id,
                                                        preview.best_event->crossing, proximity_active_);
                }
                remaining -= dt_event;
            }

            if (remaining > 0.0)
            {
                step_no_events_(remaining);
            }
        }

    private:
        struct StepPreview
        {
            std::vector<MassiveBody> massive_preview{};
            std::vector<Spacecraft> spacecraft_end{};
            std::optional<Event> best_event{};
            double t_preview_s{0.0};
        };

        static bool is_finite_vec3_(const Vec3 &v)
        {
            return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
        }

        static bool is_finite_state_(const State &state)
        {
            return is_finite_vec3_(state.position_m) && is_finite_vec3_(state.velocity_mps) &&
                   is_finite_vec3_(state.spin.axis) && std::isfinite(state.spin.angle_rad) &&
                   std::isfinite(state.spin.rate_rad_per_s);
        }

        void step_no_events_(const double dt_s)
        {
            detail::step_no_events(massive_, spacecraft_, time_s_, dt_s, plan_, cfg_.gravitational_constant,
                                   cfg_.softening_length_m, cfg_.spacecraft_integrator);
        }

        static Crossing invert_crossing_(const Crossing c)
        {
            return (c == Crossing::Enter) ? Crossing::Exit : Crossing::Enter;
        }

        void reset_proximity_tracking_()
        {
            proximity_initialized_ = false;
            proximity_active_.clear();
        }

        detail::ProximityConfig resolve_proximity_config_() const
        {
            return detail::resolve_proximity_config(cfg_.proximity.enable, cfg_.proximity.center_spacecraft_id,
                                                    cfg_.proximity.enter_radius_m, cfg_.proximity.exit_radius_m);
        }

        void prepare_proximity_tracking_(const detail::ProximityConfig &prox)
        {
            if (!prox.enabled)
            {
                proximity_center_id_ = kInvalidSpacecraftId;
                proximity_active_.clear();
                proximity_initialized_ = false;
                return;
            }

            if (proximity_center_id_ != prox.center_id)
            {
                proximity_center_id_ = prox.center_id;
                proximity_active_.clear();
                proximity_initialized_ = false;
            }

            const Spacecraft *center_ptr = spacecraft_by_id(proximity_center_id_);
            if (center_ptr == nullptr)
            {
                proximity_active_.clear();
                proximity_initialized_ = false;
                return;
            }

            if (!proximity_initialized_)
            {
                detail::initialize_proximity_active_set(spacecraft_, center_ptr->state.position_m,
                                                        proximity_center_id_, prox.exit_radius_m, proximity_active_);
                proximity_initialized_ = true;
            }
        }

        StepPreview preview_and_find_events_(const double remaining, const detail::ProximityConfig &prox) const
        {
            StepPreview out;
            out.massive_preview = massive_;
            out.t_preview_s = time_s_;
            out.spacecraft_end.resize(spacecraft_.size());

            CelestialEphemerisSegment eph_preview{};
            detail::preview_step_massive_bodies(out.massive_preview, out.t_preview_s, remaining,
                                                cfg_.gravitational_constant, cfg_.softening_length_m, &eph_preview);

            auto propagate_sc = [&](const Spacecraft &sc_start, const double t0_s, const double dt_sc_s) -> Spacecraft {
                return detail::propagate_spacecraft_with_cache(sc_start, massive_, eph_preview, plan_,
                                                               cfg_.gravitational_constant, cfg_.softening_length_m,
                                                               cfg_.spacecraft_integrator, t0_s, dt_sc_s, spacecraft_);
            };

            // Boundary events (impact/atmosphere/SOI) for all spacecraft.
            for (std::size_t sc_index = 0; sc_index < spacecraft_.size(); ++sc_index)
            {
                const Spacecraft &sc = spacecraft_[sc_index];
                Spacecraft sc1{};
                const std::optional<Event> e = find_earliest_event_in_interval(
                        massive_, eph_preview, sc, time_s_, remaining, plan_, cfg_.events, propagate_sc, &sc1);
                out.spacecraft_end[sc_index] = sc1;
                if (e.has_value() && (!out.best_event.has_value() || e->t_event_s < out.best_event->t_event_s))
                {
                    out.best_event = e;
                }
            }

            if (!prox.enabled || proximity_center_id_ == kInvalidSpacecraftId)
            {
                return out;
            }

            const Spacecraft *center_ptr = spacecraft_by_id(proximity_center_id_);
            if (center_ptr == nullptr)
            {
                return out;
            }

            const Spacecraft center0 = *center_ptr;

            // Proximity events: one center spacecraft vs all targets.
            for (const auto &target : spacecraft_)
            {
                if (target.id == proximity_center_id_)
                {
                    continue;
                }

                const bool active = proximity_active_.contains(target.id);
                const double threshold = active ? prox.exit_radius_m : prox.enter_radius_m;

                const std::optional<Event> pe = find_earliest_proximity_event_in_interval(
                        eph_preview, center0, target, time_s_, remaining, threshold, cfg_.events, propagate_sc);
                if (pe.has_value() && (!out.best_event.has_value() || pe->t_event_s < out.best_event->t_event_s))
                {
                    out.best_event = pe;
                }
            }

            return out;
        }

        double event_step_dt_(const Event &best_event, const double remaining) const
        {
            double dt_event = best_event.t_event_s - time_s_;
            const double min_step = min_event_progress_dt_(remaining);
            if (!(dt_event > min_step) || !std::isfinite(dt_event))
            {
                dt_event = min_step;
            }
            if (dt_event > remaining)
            {
                dt_event = remaining;
            }
            return dt_event;
        }

        double min_event_progress_dt_(const double remaining) const
        {
            const double abs_remaining = std::abs(remaining);
            if (!(abs_remaining > 0.0) || !std::isfinite(abs_remaining))
            {
                return 0.0;
            }
            const double tol = std::max(0.0, cfg_.events.time_tol_s);
            const double hard_floor = std::max(1e-9, abs_remaining * 1e-12);
            return std::min(abs_remaining, std::max(tol, hard_floor));
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
