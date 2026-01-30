#!/usr/bin/env python3
from __future__ import annotations

import argparse
import io
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class Trajectory:
    name: str
    t_s: List[float]
    x_m: List[float]
    y_m: List[float]
    z_m: List[float]
    vx_mps: List[float]
    vy_mps: List[float]
    vz_mps: List[float]


def _is_csv_sample_line(line: str) -> bool:
    # Expected: t_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps (floats, scientific notation ok).
    parts = line.split(",")
    if len(parts) != 7:
        return False
    try:
        float(parts[0])
        float(parts[1])
        float(parts[2])
        float(parts[3])
        float(parts[4])
        float(parts[5])
        float(parts[6])
    except ValueError:
        return False
    return True


def parse_trajectory_dump(text: str) -> Dict[str, Trajectory]:
    trajectories: Dict[str, Trajectory] = {}

    current_name: Optional[str] = None
    buf: List[Tuple[float, float, float, float, float, float, float]] = []

    def flush() -> None:
        nonlocal current_name, buf
        if current_name is None:
            return
        if not buf:
            current_name = None
            return
        t_s, x_m, y_m, z_m, vx_mps, vy_mps, vz_mps = (list(v) for v in zip(*buf))
        trajectories[current_name] = Trajectory(
            name=current_name,
            t_s=t_s,
            x_m=x_m,
            y_m=y_m,
            z_m=z_m,
            vx_mps=vx_mps,
            vy_mps=vy_mps,
            vz_mps=vz_mps,
        )
        current_name = None
        buf = []

    started = False
    for raw in io.StringIO(text):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("--- trajectory csv"):
            started = True
            continue
        if not started:
            continue

        if _is_csv_sample_line(line):
            parts = [float(p) for p in line.split(",")]
            buf.append((parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6]))
            continue

        # Section name (e.g. "earth_rel_sun").
        flush()
        current_name = line

    flush()
    return trajectories


def run_and_capture(cmd: List[str]) -> str:
    proc = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return proc.stdout


def _require_numpy_matplotlib():
    try:
        import numpy as np  # noqa: F401
        import matplotlib.pyplot as plt  # noqa: F401
    except Exception as e:
        raise RuntimeError(
            "This script requires numpy + matplotlib.\n"
            "If you're on Ubuntu/Debian: `sudo apt-get install python3-numpy python3-matplotlib`.\n"
            "Or use your venv/conda environment."
        ) from e


def _set_axes_equal_3d(ax, xs, ys, zs) -> None:
    import numpy as np

    x_min, x_max = float(np.min(xs)), float(np.max(xs))
    y_min, y_max = float(np.min(ys)), float(np.max(ys))
    z_min, z_max = float(np.min(zs)), float(np.max(zs))
    max_range = max(x_max - x_min, y_max - y_min, z_max - z_min)
    if not (max_range > 0.0):
        return
    x_mid = 0.5 * (x_min + x_max)
    y_mid = 0.5 * (y_min + y_max)
    z_mid = 0.5 * (z_min + z_max)

    ax.set_xlim(x_mid - 0.5 * max_range, x_mid + 0.5 * max_range)
    ax.set_ylim(y_mid - 0.5 * max_range, y_mid + 0.5 * max_range)
    ax.set_zlim(z_mid - 0.5 * max_range, z_mid + 0.5 * max_range)

    # Newer matplotlib supports set_box_aspect; if unavailable, axes limits above still help.
    if hasattr(ax, "set_box_aspect"):
        ax.set_box_aspect((1.0, 1.0, 1.0))


def _extent_mid_and_range_2d(xs, ys) -> Tuple[Tuple[float, float], float]:
    import numpy as np

    x_min, x_max = float(np.min(xs)), float(np.max(xs))
    y_min, y_max = float(np.min(ys)), float(np.max(ys))
    max_range = max(x_max - x_min, y_max - y_min)
    x_mid = 0.5 * (x_min + x_max)
    y_mid = 0.5 * (y_min + y_max)
    return (x_mid, y_mid), max_range


def _unit_scale_and_label(units: str) -> Tuple[float, str]:
    if units == "m":
        return 1.0, "m"
    if units == "km":
        return 1e-3, "km"
    if units == "Mm":
        return 1e-6, "Mm"
    if units == "Gm":
        return 1e-9, "Gm"
    if units == "AU":
        return 1.0 / 1.495978707e11, "AU"
    raise RuntimeError(f"Unsupported units: {units}")


def _match_names(items: List[Tuple[str, Trajectory]], pattern: str) -> List[Tuple[str, Trajectory]]:
    rx = re.compile(pattern)
    return [(k, v) for (k, v) in items if rx.search(k)]


def plot_trajectories(
    trajectories: Dict[str, Trajectory],
    select_pattern: Optional[str],
    focus_pattern: Optional[str],
    mode: str,
    save_path: Optional[str],
    title: Optional[str],
    equal_axes: bool,
    free_aspect_2d: bool,
    units: str,
) -> None:
    _require_numpy_matplotlib()
    import numpy as np
    import matplotlib.pyplot as plt

    items = list(trajectories.items())
    if select_pattern:
        items = _match_names(items, select_pattern)

    if not items:
        available = "\n".join(f"  - {k}" for k in sorted(trajectories.keys()))
        raise RuntimeError(f"No trajectories selected.\nAvailable sections:\n{available}")

    focus_items = items
    if focus_pattern:
        focus_items = _match_names(items, focus_pattern)
        if not focus_items:
            available = "\n".join(f"  - {k}" for k, _ in items)
            raise RuntimeError(f"No focus trajectories matched.\nPlotted sections:\n{available}")

    scale, unit_label = _unit_scale_and_label(units)

    if mode == "3d":
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        for name, traj in items:
            xs = scale * np.asarray(traj.x_m, dtype=float)
            ys = scale * np.asarray(traj.y_m, dtype=float)
            zs = scale * np.asarray(traj.z_m, dtype=float)
            ax.plot(xs, ys, zs, label=name, linewidth=1.2)
        ax.set_xlabel(f"x [{unit_label}]")
        ax.set_ylabel(f"y [{unit_label}]")
        ax.set_zlabel(f"z [{unit_label}]")
        if title:
            ax.set_title(title)
        ax.legend(loc="best")
        if equal_axes:
            all_x = np.concatenate([scale * np.asarray(traj.x_m, dtype=float) for _, traj in focus_items])
            all_y = np.concatenate([scale * np.asarray(traj.y_m, dtype=float) for _, traj in focus_items])
            all_z = np.concatenate([scale * np.asarray(traj.z_m, dtype=float) for _, traj in focus_items])
            _set_axes_equal_3d(ax, all_x, all_y, all_z)
    else:
        fig, ax = plt.subplots()
        if mode == "xy":
            for name, traj in items:
                ax.plot(scale * np.asarray(traj.x_m, dtype=float), scale * np.asarray(traj.y_m, dtype=float),
                        label=name, linewidth=1.2)
            ax.set_xlabel(f"x [{unit_label}]")
            ax.set_ylabel(f"y [{unit_label}]")
        elif mode == "xz":
            for name, traj in items:
                ax.plot(scale * np.asarray(traj.x_m, dtype=float), scale * np.asarray(traj.z_m, dtype=float),
                        label=name, linewidth=1.2)
            ax.set_xlabel(f"x [{unit_label}]")
            ax.set_ylabel(f"z [{unit_label}]")
        elif mode == "yz":
            for name, traj in items:
                ax.plot(scale * np.asarray(traj.y_m, dtype=float), scale * np.asarray(traj.z_m, dtype=float),
                        label=name, linewidth=1.2)
            ax.set_xlabel(f"y [{unit_label}]")
            ax.set_ylabel(f"z [{unit_label}]")
        else:
            raise RuntimeError(f"Unsupported mode: {mode}")
        if not free_aspect_2d:
            ax.set_aspect("equal", adjustable="box")
        if focus_pattern:
            if mode == "xy":
                focus_x = np.concatenate([scale * np.asarray(traj.x_m, dtype=float) for _, traj in focus_items])
                focus_y = np.concatenate([scale * np.asarray(traj.y_m, dtype=float) for _, traj in focus_items])
            elif mode == "xz":
                focus_x = np.concatenate([scale * np.asarray(traj.x_m, dtype=float) for _, traj in focus_items])
                focus_y = np.concatenate([scale * np.asarray(traj.z_m, dtype=float) for _, traj in focus_items])
            else:  # yz
                focus_x = np.concatenate([scale * np.asarray(traj.y_m, dtype=float) for _, traj in focus_items])
                focus_y = np.concatenate([scale * np.asarray(traj.z_m, dtype=float) for _, traj in focus_items])
            (x_mid, y_mid), max_range = _extent_mid_and_range_2d(focus_x, focus_y)
            if max_range > 0.0:
                ax.set_xlim(x_mid - 0.5 * max_range, x_mid + 0.5 * max_range)
                ax.set_ylim(y_mid - 0.5 * max_range, y_mid + 0.5 * max_range)
        if title:
            ax.set_title(title)
        ax.grid(True, linewidth=0.3)
        ax.legend(loc="best")

    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        plt.savefig(save_path, dpi=160, bbox_inches="tight")
    else:
        plt.show()


def main(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Visualize orbits from the orbitsim_example trajectory CSV dump.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = parser.add_mutually_exclusive_group(required=False)
    src.add_argument("--input", help="Read trajectory dump from a text file ('orbitsim_example' stdout).")
    src.add_argument(
        "--cmd",
        "--run",
        dest="cmd",
        help="Run a command and parse its stdout (e.g. --cmd ./build/orbitsim_example).",
    )
    parser.add_argument(
        "--select",
        default=None,
        help="Regex to select which sections to plot (matches section name).",
    )
    parser.add_argument(
        "--focus",
        default=None,
        help="Regex to compute axis limits from (still plots all selected). Useful when mixing very different scales.",
    )
    parser.add_argument("--mode", choices=["3d", "xy", "xz", "yz"], default="3d", help="Plot type.")
    parser.add_argument("--title", default=None, help="Optional plot title.")
    parser.add_argument("--equal", action="store_true", help="Force equal axis scaling (3D only).")
    parser.add_argument("--free-aspect", action="store_true", help="Allow non-equal aspect ratio in 2D.")
    parser.add_argument("--units", choices=["m", "km", "Mm", "Gm", "AU"], default="m", help="Axis units.")
    parser.add_argument("--save", default=None, help="Write image instead of showing a window (e.g. out.png).")
    parser.add_argument("--list", action="store_true", help="List available sections and exit.")
    args = parser.parse_args(argv)

    if args.cmd is not None:
        cmd = shlex.split(args.cmd)
        if not cmd:
            raise RuntimeError("--cmd requires a command, e.g. --cmd ./build/orbitsim_example")
        text = run_and_capture(cmd)
    elif args.input is not None:
        with open(args.input, "r", encoding="utf-8") as f:
            text = f.read()
    else:
        default_exe = os.path.join("build", "orbitsim_example")
        if os.name == "nt":
            default_exe += ".exe"
        if os.path.exists(default_exe):
            text = run_and_capture([default_exe])
        else:
            raise RuntimeError(
                "No input provided and default executable not found.\n"
                "Build the example first:\n"
                "  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release\n"
                "  cmake --build build -j\n"
                "Then run:\n"
                "  python3 tools/plot_orbits.py --cmd ./build/orbitsim_example\n"
                "Or save stdout and use --input."
            )

    trajs = parse_trajectory_dump(text)
    if not trajs:
        raise RuntimeError(
            "No trajectories found in input.\n"
            "Make sure the input includes the line starting with '--- trajectory csv'."
        )

    if args.list:
        for k in sorted(trajs.keys()):
            print(k)
        return 0

    plot_trajectories(
        trajectories=trajs,
        select_pattern=args.select,
        focus_pattern=args.focus,
        mode=args.mode,
        save_path=args.save,
        title=args.title,
        equal_axes=bool(args.equal),
        free_aspect_2d=bool(args.free_aspect),
        units=args.units,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except subprocess.CalledProcessError as e:
        print(e.stdout or str(e), file=sys.stderr)
        return_code = e.returncode if e.returncode is not None else 1
        raise SystemExit(return_code)
    except Exception as e:
        print(f"error: {e}", file=sys.stderr)
        raise SystemExit(1)
