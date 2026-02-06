#!/usr/bin/env python3
"""
Lightweight visualization for event-related features emitted by `orbitsim_example`.

It expects the example to print these sections:
  - bodies                    (id,name,radius_m,atmosphere_top_height_m,soi_radius_m)
  - event_probe_catalog       (key,title,body_name,events_section,traj_section)
  - events_*                  (t_s,type,body_id,crossing)
  - *_rel_*                   (t_s,x,y,z,vx,vy,vz)

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


@dataclass(frozen=True)
class EventProbeSpec:
    key: str
    title: str
    body_name: str
    events_section: str
    traj_section: str


def _is_section_header(line: str) -> bool:
    if not line:
        return False
    if line.startswith("---"):
        return False
    if "," in line:
        return False
    return all(c.isalnum() or c in "_-" for c in line)


def parse_output(output: str) -> Tuple[
    Dict[str, np.ndarray],
    Dict[int, BodyInfo],
    Dict[str, List[EventRow]],
    Dict[str, np.ndarray],
    List[EventProbeSpec],
]:
    lines = output.splitlines()
    traj: Dict[str, np.ndarray] = {}
    bodies: Dict[int, BodyInfo] = {}
    events: Dict[str, List[EventRow]] = {}
    apsides: Dict[str, np.ndarray] = {}
    probe_specs: List[EventProbeSpec] = []

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
            bodies[body_id] = BodyInfo(
                body_id=body_id,
                name=name,
                radius_m=radius_m,
                atmosphere_top_height_m=atmos_m,
                soi_radius_m=soi_m,
            )
            continue

        if current == "event_probe_catalog":
            parts = [p.strip() for p in line.split(",")]
            if len(parts) != 5:
                continue
            probe_specs.append(
                EventProbeSpec(
                    key=parts[0],
                    title=parts[1],
                    body_name=parts[2].lower(),
                    events_section=parts[3],
                    traj_section=parts[4],
                )
            )
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
            events.setdefault(current, []).append(
                EventRow(t_s=t_s, type=typ, body_id=body_id, crossing=crossing)
            )
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

        parts = line.split(",")
        if len(parts) == 7:
            try:
                rows_floats.append([float(x) for x in parts])
            except ValueError:
                flush_trajectory()
                current = None
            continue

    flush_trajectory()
    return traj, bodies, events, apsides, probe_specs


def vec_interp_at_t(traj: np.ndarray, t_s: float) -> np.ndarray:
    t = traj[:, 0]
    x = np.interp(t_s, t, traj[:, 1])
    y = np.interp(t_s, t, traj[:, 2])
    z = np.interp(t_s, t, traj[:, 3])
    return np.array([x, y, z], dtype=float)


def clip_trajectory_after_impact(traj_xyz_m: np.ndarray, events: List[EventRow]) -> np.ndarray:
    impact_times = [e.t_s for e in events if e.type == "impact" and e.crossing == "enter"]
    if not impact_times:
        return traj_xyz_m
    t_impact = min(impact_times)
    return traj_xyz_m[traj_xyz_m[:, 0] <= (t_impact + 5.0)]


def threshold_is_relevant(threshold_m: float, r_m: np.ndarray) -> bool:
    if not (threshold_m > 0.0 and np.isfinite(threshold_m)):
        return False
    r_min = float(np.min(r_m))
    r_max = float(np.max(r_m))
    span = max(r_max - r_min, 1.0)
    lo = r_min - 0.25 * span
    hi = r_max + 0.50 * span
    return lo <= threshold_m <= hi


def plot_event_panel_xy(
    ax,
    traj_xyz_m: np.ndarray,
    body_radius_m: float,
    title: str,
    events: List[EventRow],
    extra_rings: List[Tuple[str, float, str]],
):
    x_km = traj_xyz_m[:, 1] / 1000.0
    y_km = traj_xyz_m[:, 2] / 1000.0
    ax.plot(x_km, y_km, lw=1.0, color="tab:blue", alpha=0.8)
    ax.scatter([x_km[0]], [y_km[0]], c="g", s=40, marker="o", label="start")
    ax.scatter([x_km[-1]], [y_km[-1]], c="r", s=40, marker="x", label="end")

    r_km = body_radius_m / 1000.0
    ax.add_patch(Circle((0.0, 0.0), r_km, color="gray", alpha=0.15, label="body"))

    for label, radius_m, color in extra_rings:
        ring_km = radius_m / 1000.0
        labels = ax.get_legend_handles_labels()[1]
        ax.add_patch(
            Circle(
                (0.0, 0.0),
                ring_km,
                color=color,
                fill=False,
                ls=":",
                lw=1.0,
                alpha=0.7,
                label=label if label not in labels else None,
            )
        )

    for e in events:
        pos = vec_interp_at_t(traj_xyz_m, e.t_s) / 1000.0
        labels = ax.get_legend_handles_labels()[1]
        if e.type == "impact":
            ax.scatter([pos[0]], [pos[1]], c="red", s=60, marker="x", label="impact" if "impact" not in labels else None)
        elif e.type == "atmosphere":
            ax.scatter(
                [pos[0]],
                [pos[1]],
                c="orange",
                s=60,
                marker="^",
                label="atmosphere" if "atmosphere" not in labels else None,
            )
        elif e.type == "soi":
            ax.scatter([pos[0]], [pos[1]], c="cyan", s=60, marker="s", label="soi" if "soi" not in labels else None)

    ax.set_title(title)
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.grid(True, alpha=0.3)
    ax.axis("equal")
    ax.legend(fontsize=8, loc="best")


def plot_event_panel_r(
    ax,
    traj_xyz_m: np.ndarray,
    body_info: BodyInfo,
    title: str,
    events: List[EventRow],
    extra_thresholds: List[Tuple[str, float, str]],
):
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
    ax.legend(fontsize=8, loc="best")


def default_probe_specs() -> List[EventProbeSpec]:
    return [
        EventProbeSpec(
            key="reentry_sc",
            title="Earth reentry impact",
            body_name="earth",
            events_section="events_reentry_sc",
            traj_section="reentry_sc_rel_earth",
        ),
        EventProbeSpec(
            key="moon_soi_sc",
            title="Moon SOI ingress",
            body_name="moon",
            events_section="events_moon_soi_sc",
            traj_section="moon_soi_sc_rel_moon",
        ),
    ]


def resolve_body_info(
    spec: EventProbeSpec,
    bodies_by_name: Dict[str, BodyInfo],
    bodies_by_id: Dict[int, BodyInfo],
    probe_events: List[EventRow],
    earth_fallback: BodyInfo,
    moon_fallback: BodyInfo,
) -> BodyInfo:
    by_name = bodies_by_name.get(spec.body_name.lower())
    if by_name is not None:
        return by_name
    if probe_events:
        by_id = bodies_by_id.get(probe_events[0].body_id)
        if by_id is not None:
            return by_id
    if spec.body_name.lower() == "moon":
        return moon_fallback
    return earth_fallback


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

    traj, bodies_by_id, events_by_section, aps, probe_specs = parse_output(result.stdout)

    bodies_by_name = {b.name.lower(): b for b in bodies_by_id.values()}
    earth = bodies_by_name.get(
        "earth",
        BodyInfo(
            body_id=2,
            name="earth",
            radius_m=6.371e6,
            atmosphere_top_height_m=1.0e5,
            soi_radius_m=0.0,
        ),
    )
    moon = bodies_by_name.get(
        "moon",
        BodyInfo(
            body_id=3,
            name="moon",
            radius_m=1.7374e6,
            atmosphere_top_height_m=0.0,
            soi_radius_m=6.61e7,
        ),
    )

    if not probe_specs:
        probe_specs = default_probe_specs()

    available_specs: List[EventProbeSpec] = []
    for spec in probe_specs:
        if spec.traj_section in traj:
            available_specs.append(spec)

    if not available_specs:
        print("Error: no event probe trajectory sections were found.", file=sys.stderr)
        print(f"Probe catalog entries: {[s.traj_section for s in probe_specs]}", file=sys.stderr)
        print(f"Found trajectory sections: {list(traj.keys())}", file=sys.stderr)
        return 1

    fig, axes = plt.subplots(len(available_specs), 2, figsize=(14, 5.0 * len(available_specs)), squeeze=False)

    earth_xy_ax = None
    for row, spec in enumerate(available_specs):
        probe_events = events_by_section.get(spec.events_section, [])
        probe_traj = clip_trajectory_after_impact(traj[spec.traj_section], probe_events)
        probe_body = resolve_body_info(spec, bodies_by_name, bodies_by_id, probe_events, earth, moon)

        r_m = np.linalg.norm(probe_traj[:, 1:4], axis=1)
        extra_thresholds: List[Tuple[str, float, str]] = []
        extra_rings: List[Tuple[str, float, str]] = []

        atmosphere_top_m = probe_body.radius_m + probe_body.atmosphere_top_height_m
        if probe_body.atmosphere_top_height_m > 0.0 and threshold_is_relevant(atmosphere_top_m, r_m):
            extra_thresholds.append(("atmosphere top", atmosphere_top_m, "orange"))
            extra_rings.append(("atmosphere", atmosphere_top_m, "orange"))

        if probe_body.soi_radius_m > 0.0 and threshold_is_relevant(probe_body.soi_radius_m, r_m):
            extra_thresholds.append(("SOI", probe_body.soi_radius_m, "cyan"))
            extra_rings.append(("SOI", probe_body.soi_radius_m, "cyan"))

        plot_event_panel_xy(
            axes[row, 0],
            probe_traj,
            probe_body.radius_m,
            f"{spec.title} ({probe_body.name}-relative XY)",
            probe_events,
            extra_rings,
        )
        plot_event_panel_r(
            axes[row, 1],
            probe_traj,
            probe_body,
            f"{spec.title} |r| vs time",
            probe_events,
            extra_thresholds,
        )

        if earth_xy_ax is None and probe_body.name.lower() == "earth":
            earth_xy_ax = axes[row, 0]

    if aps and earth_xy_ax is not None:
        if "peri" in aps:
            p = aps["peri"] / 1000.0
            earth_xy_ax.scatter([p[0]], [p[1]], c="tab:purple", s=60, marker="*", label="peri (t0)")
        if "apo" in aps:
            a = aps["apo"] / 1000.0
            earth_xy_ax.scatter([a[0]], [a[1]], c="tab:pink", s=60, marker="*", label="apo (t0)")
        earth_xy_ax.legend(fontsize=8, loc="best")

    plt.tight_layout()
    out = "events.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}", file=sys.stderr)

    plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
