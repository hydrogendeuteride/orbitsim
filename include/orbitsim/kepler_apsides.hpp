#pragma once

#include "orbitsim/kepler_trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace orbitsim
{
    enum class KeplerApsisKind : uint8_t
    {
        Periapsis,
        Apoapsis,
    };

    struct KeplerApsisEvent
    {
        KeplerApsisKind kind{KeplerApsisKind::Periapsis};
        double t_s{0.0};
        State state_relative{};
        double radius_m{0.0};
    };

    struct KeplerApsisEventOptions
    {
        double min_eccentricity{1.0e-6};
        double near_parabolic_epsilon{1.0e-8};
        double time_epsilon_s{1.0e-6};
        KeplerPropagationOptions propagation{};
    };

    enum class KeplerNodeKind : uint8_t
    {
        Ascending,
        Descending,
    };

    struct KeplerNodeEvent
    {
        KeplerNodeKind kind{KeplerNodeKind::Ascending};
        double t_s{0.0};
        State state_relative{};
    };

    struct KeplerNodeEventOptions
    {
        double min_inclination_rad{1.0e-6};
        double near_parabolic_epsilon{1.0e-8};
        double time_epsilon_s{1.0e-6};
        KeplerPropagationOptions propagation{};
    };

    namespace detail
    {
        inline bool time_in_arc_(const KeplerArc &arc,
                                        const double t_s,
                                        const double eps_s)
        {
            if (!std::isfinite(t_s) || !std::isfinite(eps_s))
            {
                return false;
            }
            if (arc.t1_s >= arc.t0_s)
            {
                return t_s + eps_s >= arc.t0_s && t_s <= arc.t1_s + eps_s;
            }
            return t_s - eps_s <= arc.t0_s && t_s >= arc.t1_s - eps_s;
        }

        inline bool push_apsis_(std::vector<KeplerApsisEvent> &out,
                                      const KeplerArc &arc,
                                      const KeplerApsisKind kind,
                                      const double t_s,
                                      const KeplerPropagationOptions &propagation)
        {
            const KeplerArcSample sample = sample_kepler_arc_state(arc, t_s, propagation);
            if (!sample.ok() || !kepler_finite3_(sample.state_relative.position_m))
            {
                return false;
            }

            const double radius_m = glm::length(sample.state_relative.position_m);
            if (!(radius_m > 0.0) || !std::isfinite(radius_m))
            {
                return false;
            }

            out.push_back(KeplerApsisEvent{
                    .kind = kind,
                    .t_s = sample.t_s,
                    .state_relative = sample.state_relative,
                    .radius_m = radius_m,
            });
            return true;
        }

        inline bool push_node_(std::vector<KeplerNodeEvent> &out,
                                     const KeplerArc &arc,
                                     const KeplerNodeKind kind,
                                     const double t_s,
                                     const KeplerPropagationOptions &propagation)
        {
            const KeplerArcSample sample = sample_kepler_arc_state(arc, t_s, propagation);
            if (!sample.ok() || !kepler_finite3_(sample.state_relative.position_m))
            {
                return false;
            }

            out.push_back(KeplerNodeEvent{
                    .kind = kind,
                    .t_s = sample.t_s,
                    .state_relative = sample.state_relative,
            });
            return true;
        }
    } // namespace detail

    inline std::vector<KeplerApsisEvent> find_kepler_apsides(
            const KeplerArc &arc,
            const KeplerApsisEventOptions &options = {})
    {
        std::vector<KeplerApsisEvent> out{};
        if (!kepler_arc_valid(arc))
        {
            return out;
        }

        const OrbitalElements elements =
                orbital_elements_from_relative_state(arc.mu_m3_s2,
                                                     arc.state0_relative.position_m,
                                                     arc.state0_relative.velocity_mps);
        const double e = elements.eccentricity;
        if (!(e > options.min_eccentricity) ||
            !std::isfinite(e) ||
            !std::isfinite(elements.semi_major_axis_m) ||
            !std::isfinite(elements.mean_anomaly_rad) ||
            std::abs(e - 1.0) <= options.near_parabolic_epsilon)
        {
            return out;
        }

        const double safe_eps_s = std::max(0.0, options.time_epsilon_s);
        if (e < 1.0 - options.near_parabolic_epsilon)
        {
            const OrbitScalars scalars = orbit_scalars_from_elements(arc.mu_m3_s2, elements);
            if (!(scalars.mean_motion_radps > 0.0) ||
                !std::isfinite(scalars.mean_motion_radps))
            {
                return out;
            }

            const auto event_time = [&](const double event_true_anomaly_rad) {
                const KeplerAnomalyResult anomaly =
                        kepler_anomalies_from_true_anomaly(e, event_true_anomaly_rad);
                if (!anomaly.valid || !std::isfinite(anomaly.mean_anomaly_rad))
                {
                    return std::numeric_limits<double>::quiet_NaN();
                }

                double dM = wrap_angle_0_2pi(anomaly.mean_anomaly_rad -
                                             elements.mean_anomaly_rad);
                if (std::abs(dM) <= 1.0e-12)
                {
                    dM = 0.0;
                }
                return arc.t0_s + dM / scalars.mean_motion_radps;
            };

            const double periapsis_t_s = event_time(0.0);
            if (detail::time_in_arc_(arc, periapsis_t_s, safe_eps_s))
            {
                detail::push_apsis_(out,
                                          arc,
                                          KeplerApsisKind::Periapsis,
                                          std::clamp(periapsis_t_s,
                                                     std::min(arc.t0_s, arc.t1_s),
                                                     std::max(arc.t0_s, arc.t1_s)),
                                          options.propagation);
            }

            const double apoapsis_t_s = event_time(std::acos(-1.0));
            if (detail::time_in_arc_(arc, apoapsis_t_s, safe_eps_s))
            {
                detail::push_apsis_(out,
                                          arc,
                                          KeplerApsisKind::Apoapsis,
                                          std::clamp(apoapsis_t_s,
                                                     std::min(arc.t0_s, arc.t1_s),
                                                     std::max(arc.t0_s, arc.t1_s)),
                                          options.propagation);
            }
        }
        else if (e > 1.0 + options.near_parabolic_epsilon)
        {
            const double abs_a = -elements.semi_major_axis_m;
            if (!(abs_a > 0.0) || !std::isfinite(abs_a))
            {
                return out;
            }

            const double mean_motion_radps = std::sqrt(arc.mu_m3_s2 /
                                                       (abs_a * abs_a * abs_a));
            if (!(mean_motion_radps > 0.0) || !std::isfinite(mean_motion_radps))
            {
                return out;
            }

            const double periapsis_t_s =
                    arc.t0_s - elements.mean_anomaly_rad / mean_motion_radps;
            if (detail::time_in_arc_(arc, periapsis_t_s, safe_eps_s))
            {
                detail::push_apsis_(out,
                                          arc,
                                          KeplerApsisKind::Periapsis,
                                          std::clamp(periapsis_t_s,
                                                     std::min(arc.t0_s, arc.t1_s),
                                                     std::max(arc.t0_s, arc.t1_s)),
                                          options.propagation);
            }
        }

        std::sort(out.begin(),
                  out.end(),
                  [](const KeplerApsisEvent &a, const KeplerApsisEvent &b) {
                      if (a.t_s == b.t_s)
                      {
                          return static_cast<uint8_t>(a.kind) <
                                 static_cast<uint8_t>(b.kind);
                      }
                      return a.t_s < b.t_s;
                  });
        return out;
    }

    inline std::vector<KeplerNodeEvent> find_kepler_nodes(
            const KeplerArc &arc,
            const Vec3 &reference_axis_unit_i,
            const KeplerNodeEventOptions &options = {})
    {
        std::vector<KeplerNodeEvent> out{};
        if (!kepler_arc_valid(arc))
        {
            return out;
        }

        const OrbitalElements elements =
                orbital_elements_from_relative_state_about_axis(arc.mu_m3_s2,
                                                                arc.state0_relative.position_m,
                                                                arc.state0_relative.velocity_mps,
                                                                reference_axis_unit_i);
        const double e = elements.eccentricity;
        if (!std::isfinite(e) ||
            !std::isfinite(elements.semi_major_axis_m) ||
            !std::isfinite(elements.mean_anomaly_rad) ||
            !std::isfinite(elements.inclination_rad) ||
            !std::isfinite(elements.arg_periapsis_rad) ||
            elements.inclination_rad <= options.min_inclination_rad ||
            std::abs(e - 1.0) <= options.near_parabolic_epsilon ||
            !(e < 1.0 - options.near_parabolic_epsilon))
        {
            return out;
        }

        const OrbitScalars scalars = orbit_scalars_from_elements(arc.mu_m3_s2, elements);
        if (!(scalars.mean_motion_radps > 0.0) ||
            !std::isfinite(scalars.mean_motion_radps))
        {
            return out;
        }

        const double safe_eps_s = std::max(0.0, options.time_epsilon_s);
        const auto event_time = [&](const double event_true_anomaly_rad) {
            const KeplerAnomalyResult anomaly =
                    kepler_anomalies_from_true_anomaly(e, event_true_anomaly_rad);
            if (!anomaly.valid || !std::isfinite(anomaly.mean_anomaly_rad))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }

            double dM = wrap_angle_0_2pi(anomaly.mean_anomaly_rad -
                                         elements.mean_anomaly_rad);
            if (std::abs(dM) <= 1.0e-12)
            {
                dM = 0.0;
            }
            return arc.t0_s + dM / scalars.mean_motion_radps;
        };

        const double pi = std::acos(-1.0);
        const double node_true_anomaly[2]{
                wrap_angle_0_2pi(-elements.arg_periapsis_rad),
                wrap_angle_0_2pi(pi - elements.arg_periapsis_rad),
        };
        const KeplerNodeKind node_kind[2]{
                KeplerNodeKind::Ascending,
                KeplerNodeKind::Descending,
        };

        for (int i = 0; i < 2; ++i)
        {
            const double t_s = event_time(node_true_anomaly[i]);
            if (!detail::time_in_arc_(arc, t_s, safe_eps_s))
            {
                continue;
            }

            detail::push_node_(out,
                                     arc,
                                     node_kind[i],
                                     std::clamp(t_s,
                                                std::min(arc.t0_s, arc.t1_s),
                                                std::max(arc.t0_s, arc.t1_s)),
                                     options.propagation);
        }

        std::sort(out.begin(),
                  out.end(),
                  [](const KeplerNodeEvent &a, const KeplerNodeEvent &b) {
                      if (a.t_s == b.t_s)
                      {
                          return static_cast<uint8_t>(a.kind) <
                                 static_cast<uint8_t>(b.kind);
                      }
                      return a.t_s < b.t_s;
                  });
        return out;
    }
} // namespace orbitsim
