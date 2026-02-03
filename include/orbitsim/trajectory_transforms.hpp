#pragma once

#include "orbitsim/coordinate_frames.hpp"
#include "orbitsim/synodic.hpp"
#include "orbitsim/trajectory_types.hpp"

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

} // namespace orbitsim
