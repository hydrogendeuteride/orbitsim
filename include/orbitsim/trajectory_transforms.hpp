#pragma once

#include "orbitsim/coordinate_frames.hpp"
#include "orbitsim/frame_spec.hpp"
#include "orbitsim/spacecraft_lookup.hpp"
#include "orbitsim/synodic.hpp"
#include "orbitsim/trajectory_types.hpp"

#include <cmath>
#include <optional>
#include <vector>

namespace orbitsim
{

    // -------------------------------------------------------------------------
    // Trajectory <-> frame transforms (display/analysis only; does not affect dynamics).
    // -------------------------------------------------------------------------

    /// @brief Convert a single inertial trajectory sample into a (possibly rotating) frame.
    inline TrajectorySample inertial_sample_to_frame(const TrajectorySample &sample_in, const RotatingFrame &frame)
    {
        const State s_in = make_state(sample_in.position_m, sample_in.velocity_mps);
        const State s_out = inertial_state_to_frame(s_in, frame);
        return TrajectorySample{.t_s = sample_in.t_s, .position_m = s_out.position_m, .velocity_mps = s_out.velocity_mps};
    }

    /// @brief Convert a single trajectory sample expressed in a frame back to inertial coordinates.
    inline TrajectorySample frame_sample_to_inertial(const TrajectorySample &sample_frame, const RotatingFrame &frame)
    {
        const State s_frame = make_state(sample_frame.position_m, sample_frame.velocity_mps);
        const State s_out = frame_state_to_inertial(s_frame, frame);
        return TrajectorySample{.t_s = sample_frame.t_s, .position_m = s_out.position_m, .velocity_mps = s_out.velocity_mps};
    }

    /// @brief Convert a single trajectory sample from one frame to another (via inertial).
    inline TrajectorySample frame_sample_to_frame(const TrajectorySample &sample_in, const RotatingFrame &frame_from,
                                                  const RotatingFrame &frame_to)
    {
        const State s_from = make_state(sample_in.position_m, sample_in.velocity_mps);
        const State s_to = frame_state_to_frame(s_from, frame_from, frame_to);
        return TrajectorySample{.t_s = sample_in.t_s, .position_m = s_to.position_m, .velocity_mps = s_to.velocity_mps};
    }

    /// @brief Convert inertial trajectory samples into a constant frame.
    inline std::vector<TrajectorySample> trajectory_to_frame(const std::vector<TrajectorySample> &samples_in,
                                                             const RotatingFrame &frame)
    {
        std::vector<TrajectorySample> out;
        out.reserve(samples_in.size());
        for (const auto &s: samples_in)
        {
            out.push_back(inertial_sample_to_frame(s, frame));
        }
        return out;
    }

    /// @brief Convert trajectory samples expressed in a constant frame back to inertial coordinates.
    inline std::vector<TrajectorySample> trajectory_from_frame(const std::vector<TrajectorySample> &samples_in_frame,
                                                               const RotatingFrame &frame)
    {
        std::vector<TrajectorySample> out;
        out.reserve(samples_in_frame.size());
        for (const auto &s: samples_in_frame)
        {
            out.push_back(frame_sample_to_inertial(s, frame));
        }
        return out;
    }

    /// @brief Convert inertial trajectory samples into the time-varying body-centered inertial frame of a body (ECI/BCI).
    /// @note The input samples must be inertial (i.e., not already expressed in another frame).
    inline std::vector<TrajectorySample> trajectory_to_body_centered_inertial(const std::vector<TrajectorySample> &samples_in,
                                                                              const CelestialEphemeris &eph,
                                                                              const MassiveBody &body)
    {
        std::vector<TrajectorySample> out;
        out.reserve(samples_in.size());

        for (const auto &s: samples_in)
        {
            const RotatingFrame frame = make_body_centered_inertial_frame_at(eph, body, s.t_s);
            out.push_back(inertial_sample_to_frame(s, frame));
        }

        return out;
    }

    /// @brief Convert inertial trajectory samples into the time-varying body-fixed rotating frame of a body (ECEF).
    /// @note The input samples must be inertial (i.e., not already expressed in another frame).
    inline std::vector<TrajectorySample> trajectory_to_body_fixed(const std::vector<TrajectorySample> &samples_in,
                                                                  const CelestialEphemeris &eph,
                                                                  const MassiveBody &body)
    {
        std::vector<TrajectorySample> out;
        out.reserve(samples_in.size());

        for (const auto &s: samples_in)
        {
            const std::optional<RotatingFrame> frame = make_body_fixed_frame_at(eph, body, s.t_s);
            if (!frame.has_value())
            {
                return {};
            }
            out.push_back(inertial_sample_to_frame(s, *frame));
        }

        return out;
    }

    /// @brief Alias for trajectory_to_body_centered_inertial (commonly called ECI when the body is Earth).
    inline std::vector<TrajectorySample> trajectory_to_eci(const std::vector<TrajectorySample> &samples_in,
                                                           const CelestialEphemeris &eph,
                                                           const MassiveBody &body)
    {
        return trajectory_to_body_centered_inertial(samples_in, eph, body);
    }

    /// @brief Alias for trajectory_to_body_fixed (commonly called ECEF when the body is Earth).
    inline std::vector<TrajectorySample> trajectory_to_ecef(const std::vector<TrajectorySample> &samples_in,
                                                            const CelestialEphemeris &eph,
                                                            const MassiveBody &body)
    {
        return trajectory_to_body_fixed(samples_in, eph, body);
    }

    /// @brief Convert a single inertial trajectory sample into a synodic rotating frame.
    /// @note The input sample must be in the same inertial frame as the body states used to construct the SynodicFrame.
    inline TrajectorySample inertial_sample_to_synodic(const TrajectorySample &sample_in, const SynodicFrame &frame)
    {
        return inertial_sample_to_frame(sample_in, frame);
    }

    /// @brief Convert inertial trajectory samples into the time-varying synodic frame of bodies (A,B).
    /// @note The input samples must be inertial (i.e., not already expressed in another frame).
    inline std::vector<TrajectorySample> trajectory_to_synodic(const std::vector<TrajectorySample> &samples_in,
                                                               const CelestialEphemeris &eph,
                                                               const MassiveBody &body_a,
                                                               const MassiveBody &body_b)
    {
        std::vector<TrajectorySample> out;
        out.reserve(samples_in.size());

        for (const auto &s: samples_in)
        {
            const std::optional<SynodicFrame> frame = make_synodic_frame_at(eph, body_a, body_b, s.t_s);
            if (!frame.has_value())
            {
                return {};
            }
            out.push_back(inertial_sample_to_synodic(s, *frame));
        }

        return out;
    }

    /// @brief Convert inertial trajectory samples into a time-varying frame specified by TrajectoryFrameSpec.
    ///
    /// @note The input samples must be inertial (i.e., not already expressed in another frame).
    inline std::vector<TrajectorySample> trajectory_to_frame_spec(const std::vector<TrajectorySample> &samples_in,
                                                                  const CelestialEphemeris &eph,
                                                                  const std::vector<MassiveBody> &bodies,
                                                                  const TrajectoryFrameSpec &frame,
                                                                  const SpacecraftStateLookup &sc_lookup = nullptr)
    {
        if (samples_in.empty())
        {
            return samples_in;
        }

        const auto body_by_id = [&](const BodyId id) -> const MassiveBody * {
            if (id == kInvalidBodyId)
            {
                return nullptr;
            }
            for (const auto &b: bodies)
            {
                if (b.id == id)
                {
                    return &b;
                }
            }
            return nullptr;
        };

        switch (frame.type)
        {
        case TrajectoryFrameType::Inertial:
            return samples_in;
        case TrajectoryFrameType::BodyCenteredInertial:
        {
            const MassiveBody *body = body_by_id(frame.primary_body_id);
            if (body == nullptr)
            {
                return {};
            }
            return trajectory_to_body_centered_inertial(samples_in, eph, *body);
        }
        case TrajectoryFrameType::BodyFixed:
        {
            const MassiveBody *body = body_by_id(frame.primary_body_id);
            if (body == nullptr)
            {
                return {};
            }
            return trajectory_to_body_fixed(samples_in, eph, *body);
        }
        case TrajectoryFrameType::Synodic:
        {
            const MassiveBody *body_a = body_by_id(frame.primary_body_id);
            const MassiveBody *body_b = body_by_id(frame.secondary_body_id);
            if (body_a == nullptr || body_b == nullptr)
            {
                return {};
            }
            return trajectory_to_synodic(samples_in, eph, *body_a, *body_b);
        }
        case TrajectoryFrameType::LVLH:
        {
            if (!sc_lookup)
            {
                return {};
            }
            if (frame.target_spacecraft_id == kInvalidSpacecraftId)
            {
                return {};
            }
            if (bodies.empty())
            {
                return {};
            }

            const MassiveBody *primary_fixed = body_by_id(frame.primary_body_id);

            std::vector<std::optional<std::size_t>> eph_indices;
            if (!eph.empty())
            {
                eph_indices.reserve(bodies.size());
                for (const auto &b: bodies)
                {
                    std::size_t idx = 0;
                    if (eph.body_index_for_id(b.id, &idx))
                    {
                        eph_indices.push_back(idx);
                    }
                    else
                    {
                        eph_indices.push_back(std::nullopt);
                    }
                }
            }

            std::vector<TrajectorySample> out;
            out.reserve(samples_in.size());

            for (const auto &s: samples_in)
            {
                const std::optional<State> target_state = sc_lookup(frame.target_spacecraft_id, s.t_s);
                if (!target_state.has_value())
                {
                    return {};
                }

                const MassiveBody *primary = primary_fixed;
                std::optional<std::size_t> primary_eph_index;

                if (primary == nullptr)
                {
                    const std::size_t best = auto_select_primary_index(
                            bodies, target_state->position_m,
                            [&](std::size_t i) -> Vec3 {
                                if (!eph.empty() && i < eph_indices.size() && eph_indices[i].has_value())
                                    return eph.body_position_at(*eph_indices[i], s.t_s);
                                return bodies[i].state.position_m;
                            });

                    primary = &bodies[best];
                    if (!eph.empty() && best < eph_indices.size())
                    {
                        primary_eph_index = eph_indices[best];
                    }
                }
                else if (!eph.empty())
                {
                    for (std::size_t i = 0; i < bodies.size(); ++i)
                    {
                        if (bodies[i].id == primary->id)
                        {
                            if (i < eph_indices.size())
                            {
                                primary_eph_index = eph_indices[i];
                            }
                            break;
                        }
                    }
                }

                const State primary_state = [&]() -> State {
                    if (eph.empty() || !primary_eph_index.has_value())
                    {
                        return primary->state;
                    }
                    return eph.body_state_at(*primary_eph_index, s.t_s);
                }();
                const std::optional<RotatingFrame> lvlh = make_lvlh_frame(primary_state, *target_state);
                if (!lvlh.has_value())
                {
                    return {};
                }
                out.push_back(inertial_sample_to_frame(s, *lvlh));
            }

            return out;
        }
        }

        return samples_in;
    }

} // namespace orbitsim
