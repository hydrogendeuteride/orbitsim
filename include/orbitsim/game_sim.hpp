#pragma once

#include "orbitsim/ephemeris.hpp"
#include "orbitsim/detail/spacecraft_propagation.hpp"
#include "orbitsim/events.hpp"
#include "orbitsim/integrators.hpp"
#include "orbitsim/maneuvers.hpp"
#include "orbitsim/math.hpp"
#include "orbitsim/types.hpp"

#include <algorithm>
#include <cstddef>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace orbitsim
{

    class GameSimulation
    {
    public:
        struct BodyHandle
        {
            BodyId id{kInvalidBodyId};

            inline bool valid() const { return id != kInvalidBodyId; }
            inline operator BodyId() const { return id; }
        };

        struct SpacecraftHandle
        {
            SpacecraftId id{kInvalidSpacecraftId};

            inline bool valid() const { return id != kInvalidSpacecraftId; }
            inline operator SpacecraftId() const { return id; }
        };

        struct Config
        {
            double gravitational_constant{orbitsim::kGravitationalConstant_SI};
            double softening_length_m{0.0};
            DOPRI5Options spacecraft_integrator{};
            EventOptions events{};
            bool enable_events{true};

            struct ProximityOptions
            {
                bool enable{false};
                SpacecraftId center_spacecraft_id{kInvalidSpacecraftId};
                double enter_radius_m{0.0};
                double exit_radius_m{0.0}; // Use >= enter_radius_m for hysteresis.
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

        BodyId add_body(MassiveBody body)
        {
            if (body.id == kInvalidBodyId)
            {
                body.id = allocate_body_id();
            }
            if (body_id_to_index_.contains(body.id))
            {
                return kInvalidBodyId;
            }

            body_id_to_index_[body.id] = massive_.size();
            massive_.push_back(std::move(body));
            return massive_.back().id;
        }

        BodyId add_body_with_id(const BodyId id, MassiveBody body)
        {
            if (id == kInvalidBodyId)
            {
                return kInvalidBodyId;
            }
            body.id = id;
            return add_body(std::move(body));
        }

        BodyHandle create_body(MassiveBody body)
        {
            const BodyId id = add_body(std::move(body));
            return BodyHandle{.id = id};
        }

        BodyHandle create_body_with_id(const BodyId id, MassiveBody body)
        {
            const BodyId out = add_body_with_id(id, std::move(body));
            return BodyHandle{.id = out};
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

        SpacecraftId add_spacecraft(Spacecraft sc)
        {
            if (sc.id == kInvalidSpacecraftId || sc.id == kAllSpacecraft)
            {
                sc.id = allocate_spacecraft_id();
            }
            if (spacecraft_id_to_index_.contains(sc.id))
            {
                return kInvalidSpacecraftId;
            }

            spacecraft_id_to_index_[sc.id] = spacecraft_.size();
            spacecraft_.push_back(std::move(sc));
            return spacecraft_.back().id;
        }

        SpacecraftId add_spacecraft_with_id(const SpacecraftId id, Spacecraft sc)
        {
            if (id == kInvalidSpacecraftId || id == kAllSpacecraft)
            {
                return kInvalidSpacecraftId;
            }
            sc.id = id;
            return add_spacecraft(std::move(sc));
        }

        SpacecraftHandle create_spacecraft(Spacecraft sc)
        {
            const SpacecraftId id = add_spacecraft(std::move(sc));
            return SpacecraftHandle{.id = id};
        }

        SpacecraftHandle create_spacecraft_with_id(const SpacecraftId id, Spacecraft sc)
        {
            const SpacecraftId out = add_spacecraft_with_id(id, std::move(sc));
            return SpacecraftHandle{.id = out};
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

        void step(const double dt_s)
        {
            step(dt_s, nullptr);
        }

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
                                    eph_preview, center0, target, time_s_, remaining, threshold, cfg_.events, propagate_sc);
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
                        spacecraft_[sc_index] = propagate_spacecraft_(spacecraft_[sc_index], eph_preview, time_s_, remaining);
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

	        inline Spacecraft propagate_spacecraft_(const Spacecraft &sc0, const CelestialEphemerisSegment &eph,
	                                                const double t0_s, const double dt_s) const
	        {
	            return detail::propagate_spacecraft_in_ephemeris(
	                    sc0, massive_, eph, plan_, cfg_.gravitational_constant, cfg_.softening_length_m,
	                    cfg_.spacecraft_integrator, t0_s, dt_s);
	        }

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
        double time_s_{0.0};
        std::vector<MassiveBody> massive_{};
        std::vector<Spacecraft> spacecraft_{};
        BodyId next_body_id_{1};
        std::unordered_map<BodyId, std::size_t> body_id_to_index_{};
        SpacecraftId next_spacecraft_id_{1};
        std::unordered_map<SpacecraftId, std::size_t> spacecraft_id_to_index_{};
        ManeuverPlan plan_{};

        SpacecraftId proximity_center_id_{kInvalidSpacecraftId};
        std::unordered_set<SpacecraftId> proximity_active_{};
        bool proximity_initialized_{false};
    };

} // namespace orbitsim
