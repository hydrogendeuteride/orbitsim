#!/usr/bin/env python3
"""
Visualize Earth and Moon motion in the Synodic (rotating) frame.

Demonstrates that the Earth-Moon distance is VARIABLE in the N-body simulation:
  - Earth oscillates along the x-axis (x = -mu * separation)
  - Moon oscillates along the x-axis (x = (1-mu) * separation)
  - By definition of synodic frame, y = z = 0 for both bodies

Reads CSV output from orbitsim_example:
  - earth_synodic (Earth trajectory in synodic frame)
  - moon_synodic (Moon trajectory in synodic frame)
  - sc_synodic (Spacecraft trajectory for reference)
  - lagrange_points (L1-L5 positions at t=0)

Usage:
  python3 visualize_synodic_bodies.py
"""

import sys
import subprocess
from typing import Dict, List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def parse_output(output: str) -> Tuple[Dict[str, np.ndarray], Dict[str, List[float]]]:
    """Parse orbitsim_example output."""
    lines = output.splitlines()
    trajectories: Dict[str, np.ndarray] = {}
    lagrange_points: Dict[str, List[float]] = {}

    traj_sections = {'earth_synodic', 'moon_synodic', 'sc_synodic'}
    current = None
    rows = []

    def flush():
        nonlocal current, rows
        if current and rows:
            trajectories[current] = np.array(rows, dtype=float)
        rows = []

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith('---'):
            continue

        if line in traj_sections:
            flush()
            current = line
            continue

        if line == 'lagrange_points':
            flush()
            current = 'lagrange_points'
            continue

        if current == 'lagrange_points':
            parts = line.split(',')
            if len(parts) == 4 and parts[0] in ('L1', 'L2', 'L3', 'L4', 'L5', 'primary', 'secondary'):
                try:
                    lagrange_points[parts[0]] = [float(x) for x in parts[1:4]]
                except ValueError:
                    pass
            continue

        if current in traj_sections:
            parts = line.split(',')
            if len(parts) != 7:
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
    return trajectories, lagrange_points


def main() -> int:
    print("Building and running orbitsim_example...")

    try:
        subprocess.run(['cmake', '--build', 'build', '-j'], check=True, cwd='.')
    except subprocess.CalledProcessError:
        print("Error: Failed to build.", file=sys.stderr)
        return 1

    try:
        result = subprocess.run(['./build/orbitsim_example'],
                               capture_output=True, text=True, check=True, cwd='.')
    except subprocess.CalledProcessError as e:
        print(f"Error running orbitsim_example: {e}", file=sys.stderr)
        return 1

    trajectories, lagrange_points = parse_output(result.stdout)

    if 'earth_synodic' not in trajectories or 'moon_synodic' not in trajectories:
        print("Error: earth_synodic or moon_synodic not found in output.", file=sys.stderr)
        return 1

    earth = trajectories['earth_synodic']
    moon = trajectories['moon_synodic']
    sc = trajectories.get('sc_synodic')

    t_days = earth[:, 0] / 86400.0
    earth_x_mm = earth[:, 1] / 1e6  # Mm
    moon_x_mm = moon[:, 1] / 1e6

    # Calculate Earth-Moon separation
    separation_mm = moon_x_mm - earth_x_mm

    print(f"\n=== Analysis from C++ Output ===")
    print(f"Time span: {t_days[-1]:.1f} days, {len(earth)} points")
    print(f"\nEarth X in Synodic: {earth_x_mm.min():.4f} to {earth_x_mm.max():.4f} Mm")
    print(f"  Variation: {(earth_x_mm.max() - earth_x_mm.min()) * 1000:.1f} km")
    print(f"\nMoon X in Synodic: {moon_x_mm.min():.2f} to {moon_x_mm.max():.2f} Mm")
    print(f"  Variation: {(moon_x_mm.max() - moon_x_mm.min()) * 1000:.0f} km")
    print(f"\nEarth-Moon Separation: {separation_mm.min():.2f} to {separation_mm.max():.2f} Mm")
    print(f"  Variation: {(separation_mm.max() - separation_mm.min()) * 1000:.0f} km ({(separation_mm.max() - separation_mm.min()) / separation_mm.mean() * 100:.2f}%)")

    # Create visualization
    fig = plt.figure(figsize=(16, 12))

    # Plot 1: Earth-Moon separation over time
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(t_days, separation_mm, 'b-', lw=1.5)
    ax1.axhline(separation_mm.mean(), color='r', linestyle='--', alpha=0.5,
                label=f'Mean: {separation_mm.mean():.2f} Mm')
    ax1.fill_between(t_days, separation_mm.min(), separation_mm, alpha=0.2)
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Earth-Moon Distance (Mm)')
    ax1.set_title('Earth-Moon Separation Over Time\n(Variable due to N-body dynamics)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Plot 2: Earth and Moon X-position in synodic frame
    ax2 = fig.add_subplot(2, 2, 2)
    ax2_earth = ax2
    ax2_moon = ax2.twinx()

    line1, = ax2_earth.plot(t_days, earth_x_mm, 'b-', lw=1.5, label='Earth')
    line2, = ax2_moon.plot(t_days, moon_x_mm, 'gray', lw=1.5, label='Moon')

    ax2_earth.set_xlabel('Time (days)')
    ax2_earth.set_ylabel('Earth X (Mm)', color='blue')
    ax2_moon.set_ylabel('Moon X (Mm)', color='gray')
    ax2_earth.tick_params(axis='y', labelcolor='blue')
    ax2_moon.tick_params(axis='y', labelcolor='gray')
    ax2.set_title('Earth & Moon X-Position in Synodic Frame\n(Both oscillate along x-axis)')
    ax2.grid(True, alpha=0.3)
    ax2.legend([line1, line2], ['Earth', 'Moon'], loc='best')

    # Plot 3: XY plane view with Earth/Moon oscillation
    ax3 = fig.add_subplot(2, 2, 3)

    # Earth trajectory (oscillates along x-axis only)
    scatter_earth = ax3.scatter(earth_x_mm, np.zeros_like(earth_x_mm) - 5,
                                c=t_days, cmap='Blues', s=5, alpha=0.7, label='Earth')

    # Moon trajectory (oscillates along x-axis only)
    scatter_moon = ax3.scatter(moon_x_mm, np.zeros_like(moon_x_mm) + 5,
                               c=t_days, cmap='Oranges', s=5, alpha=0.7, label='Moon')

    # Barycenter
    ax3.axvline(0, color='k', linestyle='--', alpha=0.3)
    ax3.scatter([0], [0], c='black', s=100, marker='+', zorder=10, label='Barycenter')

    # Lagrange points at t=0
    if 'L1' in lagrange_points:
        ax3.axvline(lagrange_points['L1'][0]/1e6, color='red', linestyle=':', alpha=0.5)
        ax3.text(lagrange_points['L1'][0]/1e6, 15, 'L1', ha='center', fontsize=9, color='red')

    ax3.set_xlabel('X (Mm)')
    ax3.set_ylabel('Offset for visualization')
    ax3.set_title('Earth & Moon Oscillation in Synodic Frame\n(y=0 by definition, offset here for visibility)')
    ax3.set_ylim(-20, 25)
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc='upper right')
    plt.colorbar(scatter_moon, ax=ax3, label='Time (days)', shrink=0.8)

    # Plot 4: Full synodic view with spacecraft (if available)
    ax4 = fig.add_subplot(2, 2, 4)

    if sc is not None:
        sc_x = sc[:, 1] / 1e6
        sc_y = sc[:, 2] / 1e6
        scatter_sc = ax4.scatter(sc_x, sc_y, c=sc[:, 0]/86400, cmap='viridis', s=2, alpha=0.7)
        plt.colorbar(scatter_sc, ax=ax4, label='Time (days)', shrink=0.8)

    # Earth and Moon at start/end positions
    ax4.scatter([earth_x_mm[0]], [0], c='blue', s=200, marker='o', edgecolors='black',
                linewidths=2, zorder=10, label=f'Earth (start)')
    ax4.scatter([earth_x_mm[-1]], [0], c='lightblue', s=200, marker='o', edgecolors='blue',
                linewidths=2, zorder=10, label=f'Earth (end)')

    ax4.scatter([moon_x_mm[0]], [0], c='gray', s=100, marker='o', edgecolors='black',
                linewidths=2, zorder=10, label=f'Moon (start)')
    ax4.scatter([moon_x_mm[-1]], [0], c='lightgray', s=100, marker='o', edgecolors='gray',
                linewidths=2, zorder=10, label=f'Moon (end)')

    # Lagrange points
    lp_colors = {'L1': 'red', 'L2': 'orange', 'L3': 'purple', 'L4': 'cyan', 'L5': 'magenta'}
    for name, coords in lagrange_points.items():
        if name in lp_colors:
            ax4.scatter([coords[0]/1e6], [coords[1]/1e6], c=lp_colors[name], s=80, marker='*',
                       zorder=5, label=name)

    ax4.axvline(0, color='k', linestyle='--', alpha=0.3)
    ax4.scatter([0], [0], c='black', s=50, marker='+', zorder=10)

    ax4.set_xlabel('X (Mm)')
    ax4.set_ylabel('Y (Mm)')
    ax4.set_title('Full Synodic Frame View\n(Spacecraft trajectory + Earth/Moon positions)')
    ax4.grid(True, alpha=0.3)
    ax4.axis('equal')
    ax4.legend(loc='best', fontsize=7)

    # Add overall title
    fig.suptitle('Earth-Moon System in Synodic (Rotating) Frame\n'
                 'Demonstrates VARIABLE Earth-Moon distance from N-body simulation',
                 fontsize=14, fontweight='bold', y=1.0)

    plt.tight_layout()

    out_file = 'synodic_bodies.png'
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {out_file}")
    plt.show()

    return 0


if __name__ == '__main__':
    sys.exit(main())
