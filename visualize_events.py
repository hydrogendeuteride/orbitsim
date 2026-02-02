#!/usr/bin/env python3
"""
Lightweight visualization for event-related features emitted by `orbitsim_example`.

It expects the example to print these sections:
  - bodies                    (id,name,radius_m,atmosphere_top_height_m,soi_radius_m)
  - events_reentry_sc          (t_s,type,body_id,crossing)
  - reentry_sc_rel_earth       (t_s,x,y,z,vx,vy,vz)  [Earth-relative]
  - events_moon_soi_sc         (t_s,type,body_id,crossing)
  - moon_soi_sc_rel_moon       (t_s,x,y,z,vx,vy,vz)  [Moon-relative]

Usage:
  python3 visualize_events.py
"""

import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle


@dataclass(frozen=True)
class BodyInfo:
    body_id: int
    name: str
    radius_m: float
    atmosphere_top_height_m: float
    soi_radius_m: float


@dataclass(frozen=True)
class EventRow:
    t_s: float
    type: str
    body_id: int
    crossing: str


def _is_section_header(line: str) -> bool:
    if not line:
        return False
    if line.startswith("---"):
        return False
    if "," in line:
        return False
    # section names in this repo are snake_case-ish
    return all(c.isalnum() or c in "_-" for c in line)


def parse_output(output: str) -> Tuple[Dict[str, np.ndarray], Dict[int, BodyInfo], Dict[str, List[EventRow]], Dict[str, np.ndarray]]:
    lines = output.splitlines()
    traj: Dict[str, np.ndarray] = {}
    bodies: Dict[int, BodyInfo] = {}
    events: Dict[str, List[EventRow]] = {}
    apsides: Dict[str, np.ndarray] = {}

    current: Optional[str] = None
    rows_floats: List[List[float]] = []

    def flush_trajectory():
        nonlocal current, rows_floats
        if current and rows_floats:
            traj[current] = np.array(rows_floats, dtype=float)
        rows_floats = []

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("---"):
            continue

        if _is_section_header(line):
            flush_trajectory()
            current = line
            continue

        if current == "bodies":
            parts = line.split(",")
            if len(parts) != 5:
                continue
            try:
                body_id = int(parts[0])
                name = parts[1]
                radius_m = float(parts[2])
                atmos_m = float(parts[3])
                soi_m = float(parts[4])
            except ValueError:
                continue
            bodies[body_id] = BodyInfo(body_id=body_id, name=name, radius_m=radius_m, atmosphere_top_height_m=atmos_m, soi_radius_m=soi_m)
            continue

        if current and current.startswith("events_"):
            parts = line.split(",")
            if len(parts) != 4:
                continue
            try:
                t_s = float(parts[0])
                typ = parts[1].strip()
                body_id = int(parts[2])
                crossing = parts[3].strip()
            except ValueError:
                continue
            events.setdefault(current, []).append(EventRow(t_s=t_s, type=typ, body_id=body_id, crossing=crossing))
            continue

        if current == "apsides_sc_rel_earth_t0":
            parts = line.split(",")
            if len(parts) != 4:
                continue
            try:
                label = parts[0].strip()
                xyz = np.array([float(parts[1]), float(parts[2]), float(parts[3])], dtype=float)
            except ValueError:
                continue
            apsides[label] = xyz
            continue

        # Default: try to parse 7-float trajectory rows.
        parts = line.split(",")
        if len(parts) == 7:
            try:
                rows_floats.append([float(x) for x in parts])
            except ValueError:
                flush_trajectory()
                current = None
            continue

    flush_trajectory()
    return traj, bodies, events, apsides


def vec_interp_at_t(traj: np.ndarray, t_s: float) -> np.ndarray:
    t = traj[:, 0]
    x = np.interp(t_s, t, traj[:, 1])
    y = np.interp(t_s, t, traj[:, 2])
    z = np.interp(t_s, t, traj[:, 3])
    return np.array([x, y, z], dtype=float)


def plot_event_panel_xy(ax, traj_xyz_m: np.ndarray, body_radius_m: float, title: str, events: List[EventRow]):
    x_km = traj_xyz_m[:, 0] / 1000.0
    y_km = traj_xyz_m[:, 1] / 1000.0
    ax.plot(x_km, y_km, lw=1.0, color="tab:blue", alpha=0.8)
    ax.scatter([x_km[0]], [y_km[0]], c="g", s=40, marker="o", label="start")
    ax.scatter([x_km[-1]], [y_km[-1]], c="r", s=40, marker="x", label="end")

    r_km = body_radius_m / 1000.0
    ax.add_patch(Circle((0.0, 0.0), r_km, color="gray", alpha=0.15, label="body"))

    for e in events:
        pos = vec_interp_at_t(traj_xyz_m, e.t_s) / 1000.0
        if e.type == "impact":
            ax.scatter([pos[0]], [pos[1]], c="red", s=60, marker="x", label="impact" if "impact" not in ax.get_legend_handles_labels()[1] else None)
        elif e.type == "atmosphere":
            ax.scatter([pos[0]], [pos[1]], c="orange", s=60, marker="^", label="atmosphere" if "atmosphere" not in ax.get_legend_handles_labels()[1] else None)
        elif e.type == "soi":
            ax.scatter([pos[0]], [pos[1]], c="cyan", s=60, marker="s", label="soi" if "soi" not in ax.get_legend_handles_labels()[1] else None)

    ax.set_title(title)
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.grid(True, alpha=0.3)
    ax.axis("equal")


def plot_event_panel_r(ax, traj_xyz_m: np.ndarray, body_info: BodyInfo, title: str, events: List[EventRow], extra_thresholds: List[Tuple[str, float, str]]):
    t_s = traj_xyz_m[:, 0]
    r_m = np.linalg.norm(traj_xyz_m[:, 1:4], axis=1)

    ax.plot(t_s, r_m / 1000.0, lw=1.0, color="tab:blue")
    ax.axhline(body_info.radius_m / 1000.0, color="gray", ls="--", lw=1.0, alpha=0.8, label="radius")

    for name, val_m, color in extra_thresholds:
        if val_m > 0.0 and np.isfinite(val_m):
            ax.axhline(val_m / 1000.0, color=color, ls=":", lw=1.0, alpha=0.9, label=name)

    for e in events:
        color = {"impact": "red", "atmosphere": "orange", "soi": "cyan"}.get(e.type, "k")
        ax.axvline(e.t_s, color=color, lw=1.0, alpha=0.7)

    ax.set_title(title)
    ax.set_xlabel("t [s]")
    ax.set_ylabel("|r| [km]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc="best")


def main() -> int:
    print("Building and running orbitsim_example...")

    try:
        subprocess.run(["cmake", "--build", "build", "-j"], check=True, cwd=".")
    except subprocess.CalledProcessError:
        print("Error: build failed. Run: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release", file=sys.stderr)
        return 1

    try:
        result = subprocess.run(["./build/orbitsim_example"], capture_output=True, text=True, check=True, cwd=".")
    except subprocess.CalledProcessError as e:
        print(f"Error running orbitsim_example: {e}", file=sys.stderr)
        return 1

    traj, bodies_by_id, events_by_section, aps = parse_output(result.stdout)

    # Lookup by name where possible (the example uses stable names).
    bodies_by_name = {b.name: b for b in bodies_by_id.values()}
    earth = bodies_by_name.get("earth", BodyInfo(body_id=2, name="earth", radius_m=6.371e6, atmosphere_top_height_m=1.0e5, soi_radius_m=0.0))
    moon = bodies_by_name.get("moon", BodyInfo(body_id=3, name="moon", radius_m=1.7374e6, atmosphere_top_height_m=0.0, soi_radius_m=6.61e7))

    if "reentry_sc_rel_earth" not in traj:
        print("Error: missing section 'reentry_sc_rel_earth' in output.", file=sys.stderr)
        print(f"Found trajectory sections: {list(traj.keys())}", file=sys.stderr)
        return 1
    if "moon_soi_sc_rel_moon" not in traj:
        print("Error: missing section 'moon_soi_sc_rel_moon' in output.", file=sys.stderr)
        print(f"Found trajectory sections: {list(traj.keys())}", file=sys.stderr)
        return 1

    reentry_traj = traj["reentry_sc_rel_earth"]
    moon_soi_traj = traj["moon_soi_sc_rel_moon"]

    reentry_events = events_by_section.get("events_reentry_sc", [])
    moon_soi_events = events_by_section.get("events_moon_soi_sc", [])

    # Clip reentry samples at impact to avoid point-mass singular behavior if someone increases duration later.
    impact_times = [e.t_s for e in reentry_events if e.type == "impact" and e.crossing == "enter"]
    if impact_times:
        t_imp = min(impact_times)
        reentry_traj = reentry_traj[reentry_traj[:, 0] <= (t_imp + 5.0)]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    plot_event_panel_xy(
        axes[0, 0],
        reentry_traj,
        earth.radius_m,
        "Reentry probe (Earth-relative XY)",
        reentry_events,
    )
    if earth.atmosphere_top_height_m > 0.0:
        axes[0, 0].add_patch(Circle((0.0, 0.0), (earth.radius_m + earth.atmosphere_top_height_m) / 1000.0, color="orange", fill=False, ls=":", lw=1.0, alpha=0.6))

    plot_event_panel_r(
        axes[0, 1],
        reentry_traj,
        earth,
        "Reentry probe |r| vs time",
        reentry_events,
        extra_thresholds=[
            ("atmosphere top", earth.radius_m + earth.atmosphere_top_height_m, "orange"),
        ],
    )

    plot_event_panel_xy(
        axes[1, 0],
        moon_soi_traj,
        moon.radius_m,
        "Moon SOI probe (Moon-relative XY)",
        moon_soi_events,
    )
    if moon.soi_radius_m > 0.0 and np.isfinite(moon.soi_radius_m):
        axes[1, 0].add_patch(Circle((0.0, 0.0), moon.soi_radius_m / 1000.0, color="cyan", fill=False, ls=":", lw=1.0, alpha=0.6))

    plot_event_panel_r(
        axes[1, 1],
        moon_soi_traj,
        moon,
        "Moon SOI probe |r| vs time",
        moon_soi_events,
        extra_thresholds=[
            ("SOI", moon.soi_radius_m, "cyan"),
        ],
    )

    if aps:
        # Optional markers for osculating apsides at t=0 in Earth-relative frame.
        if "peri" in aps:
            p = aps["peri"] / 1000.0
            axes[0, 0].scatter([p[0]], [p[1]], c="tab:purple", s=60, marker="*", label="peri (t0)")
        if "apo" in aps:
            a = aps["apo"] / 1000.0
            axes[0, 0].scatter([a[0]], [a[1]], c="tab:pink", s=60, marker="*", label="apo (t0)")
        axes[0, 0].legend(fontsize=9, loc="best")

    plt.tight_layout()
    out = "events.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}", file=sys.stderr)
    plt.show()  # comment out for headless environments
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

