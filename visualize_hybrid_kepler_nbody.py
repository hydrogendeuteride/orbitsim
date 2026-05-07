#!/usr/bin/env python3
import argparse
import math
from pathlib import Path

import numpy as np

G = 6.67430e-11
M_EARTH = 5.97219e24
M_MOON = 7.34767309e22
MU_EARTH = G * M_EARTH
EARTH_MOON_R = 384_400_000.0
EARTH_RADIUS = 6_371_000.0
MOON_RADIUS = 1_737_400.0


def earth_moon_initial_state() -> np.ndarray:
    total_mass = M_EARTH + M_MOON
    earth_x = -EARTH_MOON_R * M_MOON / total_mass
    moon_x = EARTH_MOON_R * M_EARTH / total_mass
    omega = math.sqrt(G * total_mass / EARTH_MOON_R**3)

    return np.array(
        [
            earth_x, 0.0, 0.0,
            0.0, omega * earth_x, 0.0,
            moon_x, 0.0, 0.0,
            0.0, omega * moon_x, 0.0,
        ],
        dtype=float,
    )


def nbody_derivative(y: np.ndarray) -> np.ndarray:
    earth_r, earth_v = y[0:3], y[3:6]
    moon_r, moon_v = y[6:9], y[9:12]
    dr = moon_r - earth_r
    inv_r3 = 1.0 / np.linalg.norm(dr) ** 3

    dy = np.zeros_like(y)
    dy[0:3] = earth_v
    dy[3:6] = G * M_MOON * dr * inv_r3
    dy[6:9] = moon_v
    dy[9:12] = -G * M_EARTH * dr * inv_r3
    return dy


def rk4_step(y: np.ndarray, dt: float) -> np.ndarray:
    k1 = nbody_derivative(y)
    k2 = nbody_derivative(y + 0.5 * dt * k1)
    k3 = nbody_derivative(y + 0.5 * dt * k2)
    k4 = nbody_derivative(y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def solve_kepler_elliptic(mean_anomaly: float, eccentricity: float) -> float:
    mean_anomaly = math.fmod(mean_anomaly, 2.0 * math.pi)
    eccentric_anomaly = mean_anomaly
    for _ in range(12):
        f = eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly) - mean_anomaly
        df = 1.0 - eccentricity * math.cos(eccentric_anomaly)
        eccentric_anomaly -= f / df
    return eccentric_anomaly


def spacecraft_kepler_relative(t: float) -> np.ndarray:
    a = 70_000_000.0
    e = 0.45
    n = math.sqrt(MU_EARTH / a**3)
    E = solve_kepler_elliptic(n * t, e)
    x = a * (math.cos(E) - e)
    y = a * math.sqrt(1.0 - e * e) * math.sin(E)
    return np.array([x, y, 0.0], dtype=float)


def simulate(days: float, dt: float):
    steps = int(math.ceil(days * 86_400.0 / dt))
    times = np.arange(steps + 1, dtype=float) * dt

    earth = np.zeros((steps + 1, 3))
    moon = np.zeros((steps + 1, 3))
    ship_rel = np.zeros((steps + 1, 3))
    ship_world = np.zeros((steps + 1, 3))

    y = earth_moon_initial_state()
    for i, t in enumerate(times):
        earth[i] = y[0:3]
        moon[i] = y[6:9]
        ship_rel[i] = spacecraft_kepler_relative(t)
        ship_world[i] = earth[i] + ship_rel[i]
        if i < steps:
            y = rk4_step(y, dt)

    return times, earth, moon, ship_rel, ship_world


def plot(times, earth, moon, ship_rel, ship_world, out: Path, show: bool) -> None:
    if not show:
        import matplotlib

        matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    scale = 1.0e6
    moon_rel = moon - earth

    fig, (ax_world, ax_earth) = plt.subplots(1, 2, figsize=(14, 6))

    ax_world.plot(earth[:, 0] / scale, earth[:, 1] / scale, label="Earth N-body", color="tab:blue")
    ax_world.plot(moon[:, 0] / scale, moon[:, 1] / scale, label="Moon N-body", color="0.4")
    ax_world.plot(ship_world[:, 0] / scale, ship_world[:, 1] / scale, label="Ship Kepler + Earth", color="tab:orange")
    ax_world.set_title("Barycentric inertial")
    ax_world.set_xlabel("x [Mm]")
    ax_world.set_ylabel("y [Mm]")
    ax_world.grid(True, alpha=0.3)
    ax_world.axis("equal")
    ax_world.legend(fontsize=8)

    ax_earth.plot(moon_rel[:, 0] / scale, moon_rel[:, 1] / scale, label="Moon relative to Earth", color="0.4")
    ax_earth.plot(ship_rel[:, 0] / scale, ship_rel[:, 1] / scale, label="Ship Kepler relative to Earth", color="tab:orange")
    ax_earth.add_patch(Circle((0.0, 0.0), 3.0 * EARTH_RADIUS / scale, color="tab:blue", alpha=0.5, label="Earth x3"))
    ax_earth.add_patch(
        Circle(
            (moon_rel[0, 0] / scale, moon_rel[0, 1] / scale),
            6.0 * MOON_RADIUS / scale,
            color="0.5",
            alpha=0.5,
            label="Moon x6 at t0",
        )
    )
    ax_earth.scatter([ship_rel[0, 0] / scale], [ship_rel[0, 1] / scale], color="green", s=30, label="Ship start")
    ax_earth.set_title("Earth-centered hybrid view")
    ax_earth.set_xlabel("x [Mm]")
    ax_earth.set_ylabel("y [Mm]")
    ax_earth.grid(True, alpha=0.3)
    ax_earth.axis("equal")
    ax_earth.legend(fontsize=8)

    fig.suptitle(
        f"Earth-Moon N-body + Earth-primary Kepler spacecraft ({times[-1] / 86_400.0:.1f} days)",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}")
    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--days", type=float, default=7.0)
    parser.add_argument("--dt", type=float, default=600.0)
    parser.add_argument("--out", type=Path, default=Path("hybrid_kepler_nbody.png"))
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    times, earth, moon, ship_rel, ship_world = simulate(args.days, args.dt)
    plot(times, earth, moon, ship_rel, ship_world, args.out, args.show)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
