#!/usr/bin/env python3
"""
Visualize the Earth-centered Lambert transfer example emitted by `orbitsim_example`.

It expects the example to print these sections:
  - lambert_info           (key,value lines)
  - lambert_sc_rel_earth   (t_s,x,y,z,vx,vy,vz)
  - lambert_moon_rel_earth (t_s,x,y,z,vx,vy,vz)

Usage:
  python3 visualize_lambert_transfer.py
"""

import sys
import subprocess
from typing import Dict, Tuple

import numpy as np
import matplotlib.pyplot as plt


def parse_output(output: str) -> Tuple[Dict[str, np.ndarray], Dict[str, str]]:
    lines = output.splitlines()
    data: Dict[str, np.ndarray] = {}
    info: Dict[str, str] = {}

    current = None
    rows = []

    def flush():
        nonlocal current, rows
        if current and rows:
            data[current] = np.array(rows, dtype=float)
        rows = []

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("---"):
            continue

        # Section headers (exact names we emit)
        if line in ("lambert_info", "lambert_sc_rel_earth", "lambert_moon_rel_earth"):
            flush()
            current = line
            continue

        if current == "lambert_info":
            parts = line.split(",")
            if len(parts) >= 2:
                info[parts[0]] = ",".join(parts[1:])
            continue

        if current in ("lambert_sc_rel_earth", "lambert_moon_rel_earth"):
            parts = line.split(",")
            if len(parts) != 7:
                # Any non-data line terminates the current trajectory block, since the example prints many
                # other trajectory sections afterwards (also 7-column CSV).
                flush()
                current = None
                continue
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                flush()
                current = None
                continue

    flush()
    return data, info


def main() -> int:
    print("Building and running orbitsim_example...")

    try:
        subprocess.run(["cmake", "--build", "build", "-j"], check=True, cwd=".")
    except subprocess.CalledProcessError:
        print("Error: Failed to build. Run: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release", file=sys.stderr)
        return 1

    try:
        result = subprocess.run(["./build/orbitsim_example"], capture_output=True, text=True, check=True, cwd=".")
    except subprocess.CalledProcessError as e:
        print(f"Error running orbitsim_example: {e}", file=sys.stderr)
        return 1

    data, info = parse_output(result.stdout)
    if "lambert_sc_rel_earth" not in data or "lambert_moon_rel_earth" not in data:
        print("Error: Lambert sections not found in output.", file=sys.stderr)
        print(f"Found sections: {list(data.keys())}", file=sys.stderr)
        return 1

    sc = data["lambert_sc_rel_earth"]
    moon = data["lambert_moon_rel_earth"]

    # Convert to km for display
    sc_xy_km = sc[:, 1:3] / 1000.0
    moon_xy_km = moon[:, 1:3] / 1000.0

    t_s = sc[:, 0]
    dt_s = float(info.get("dt_s", "nan"))
    best_dv = float(info.get("best_dv_mps", "nan"))
    num_solutions = info.get("num_solutions", "?")

    fig = plt.figure(figsize=(14, 6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.plot(moon_xy_km[:, 0], moon_xy_km[:, 1], color="gray", lw=1.0, alpha=0.7, label="Moon (rel Earth)")
    ax1.plot(sc_xy_km[:, 0], sc_xy_km[:, 1], color="tab:blue", lw=1.2, label="Lambert-initialized SC (N-body)")
    ax1.scatter([sc_xy_km[0, 0]], [sc_xy_km[0, 1]], c="g", s=40, marker="o", label="Start")
    ax1.scatter([sc_xy_km[-1, 0]], [sc_xy_km[-1, 1]], c="r", s=40, marker="x", label="End")
    ax1.scatter([moon_xy_km[-1, 0]], [moon_xy_km[-1, 1]], c="k", s=30, marker=".", label="Moon @ arrival")

    ax1.set_title("Earth-centered view (XY plane)")
    ax1.set_xlabel("x [km]")
    ax1.set_ylabel("y [km]")
    ax1.grid(True, alpha=0.3)
    ax1.axis("equal")
    ax1.legend(fontsize=9, loc="best")

    r_sc_km = np.linalg.norm(sc[:, 1:4], axis=1) / 1000.0
    r_moon_km = np.linalg.norm(moon[:, 1:4], axis=1) / 1000.0
    ax2.plot(t_s / 86400.0, r_sc_km, color="tab:blue", lw=1.2, label="SC |r|")
    ax2.plot(t_s / 86400.0, r_moon_km, color="gray", lw=1.0, alpha=0.7, label="Moon |r|")
    ax2.set_title("Radius vs time")
    ax2.set_xlabel("t [days]")
    ax2.set_ylabel("|r| [km]")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9, loc="best")

    title = f"Lambert transfer example: dt={dt_s/86400.0:.2f} days, solutions={num_solutions}"
    if np.isfinite(best_dv):
        title += f", best Δv≈{best_dv:.0f} m/s"
    fig.suptitle(title, fontsize=12)

    out = "lambert_transfer.png"
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}", file=sys.stderr)
    plt.show()  # comment out for headless environments
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
