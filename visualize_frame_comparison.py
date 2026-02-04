#!/usr/bin/env python3
"""
Visualize spacecraft trajectories in ECI (Earth-Centered Inertial) vs Synodic frame.

Demonstrates how the same trajectory appears differently in:
  - ECI: Earth-centered, non-rotating (Keplerian orbits are closed)
  - Synodic: Earth-Moon rotating frame (complex patterns, Lagrange points fixed)

Reads CSV output from orbitsim_example:
  - elliptical_eci / elliptical_synodic (high elliptical orbit)
  - l1_eci / l1_synodic (near L1 Lagrange point)
  - circular_eci / circular_synodic (circular orbit at ~100,000 km)
  - lagrange_points (for reference in synodic plots)
  - elliptical_orbit_info, l1_spacecraft_info, circular_orbit_info

Usage:
  python3 visualize_frame_comparison.py
"""

import sys
import subprocess
from typing import Dict, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def parse_output(output: str) -> Tuple[Dict[str, np.ndarray], Dict[str, dict], Dict[str, list]]:
    """Parse the orbitsim_example output into trajectories and info."""
    lines = output.splitlines()
    trajectories: Dict[str, np.ndarray] = {}
    info: Dict[str, dict] = {}
    lagrange_points: Dict[str, list] = {}

    current = None
    rows = []

    # Recognized trajectory sections
    traj_sections = {
        'elliptical_eci', 'elliptical_synodic',
        'l1_eci', 'l1_synodic',
        'circular_eci', 'circular_synodic',
    }

    # Recognized info sections
    info_sections = {
        'elliptical_orbit_info',
        'l1_spacecraft_info',
        'circular_orbit_info',
    }

    def flush():
        nonlocal current, rows
        if current and rows:
            trajectories[current] = np.array(rows, dtype=float)
        rows = []

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith('---'):
            continue

        # Check for trajectory section headers
        if line in traj_sections:
            flush()
            current = line
            continue

        # Check for info section headers
        if line in info_sections:
            flush()
            current = line
            info[current] = {}
            continue

        # Check for lagrange_points section
        if line == 'lagrange_points':
            flush()
            current = 'lagrange_points'
            continue

        # Parse lagrange points (only L1-L5, primary, secondary)
        if current == 'lagrange_points':
            parts = line.split(',')
            if len(parts) == 4:
                name = parts[0]
                # Only accept known Lagrange point names
                if name in ('L1', 'L2', 'L3', 'L4', 'L5', 'primary', 'secondary'):
                    try:
                        coords = [float(x) for x in parts[1:4]]
                        lagrange_points[name] = coords
                    except ValueError:
                        pass
                else:
                    # Unknown name means we've moved to a different section
                    current = None
            else:
                # Non-4-column line means section ended
                current = None
            continue

        # Parse info sections (key,value format)
        if current in info_sections:
            parts = line.split(',')
            if len(parts) >= 2:
                key = parts[0]
                try:
                    # Try to parse as number
                    info[current][key] = float(parts[1])
                except ValueError:
                    info[current][key] = ','.join(parts[1:])
            continue

        # Parse trajectory data (7 columns: t,x,y,z,vx,vy,vz)
        if current in traj_sections:
            parts = line.split(',')
            if len(parts) != 7:
                # Non-data line terminates the block
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
    return trajectories, info, lagrange_points


def plot_comparison(ax_eci, ax_syn, eci_data: np.ndarray, syn_data: np.ndarray,
                    lagrange_points: Dict[str, list], title: str, color: str = 'tab:blue'):
    """Plot ECI and Synodic trajectories side by side with time coloring."""

    # Convert to thousands of km for readability
    scale = 1e6  # to Mm (megameters)
    scale_label = 'Mm'

    eci_x = eci_data[:, 1] / scale
    eci_y = eci_data[:, 2] / scale
    syn_x = syn_data[:, 1] / scale
    syn_y = syn_data[:, 2] / scale

    # Time in days for coloring
    t_days = eci_data[:, 0] / 86400.0

    # ECI plot - use scatter with time coloring to show evolution
    scatter_eci = ax_eci.scatter(eci_x, eci_y, c=t_days, s=2, cmap='viridis', alpha=0.8)
    ax_eci.scatter([eci_x[0]], [eci_y[0]], c='lime', s=80, marker='o', zorder=5,
                   edgecolors='black', linewidths=1, label='Start')
    ax_eci.scatter([eci_x[-1]], [eci_y[-1]], c='red', s=80, marker='X', zorder=5,
                   edgecolors='black', linewidths=1, label='End')

    # Earth at origin (ECI)
    earth_radius_mm = 6.371e6 / scale
    earth_eci = Circle((0, 0), earth_radius_mm * 3, color='blue', alpha=0.6)
    ax_eci.add_patch(earth_eci)
    ax_eci.scatter([0], [0], c='blue', s=50, marker='o', label='Earth')

    ax_eci.set_xlabel(f'X ({scale_label})')
    ax_eci.set_ylabel(f'Y ({scale_label})')
    ax_eci.set_title(f'{title}\nECI (Earth-Centered Inertial)')
    ax_eci.grid(True, alpha=0.3)
    ax_eci.axis('equal')
    ax_eci.legend(fontsize=8, loc='best')

    # Add colorbar for ECI
    cbar_eci = plt.colorbar(scatter_eci, ax=ax_eci, shrink=0.8)
    cbar_eci.set_label('Time (days)', fontsize=8)

    # Synodic plot - use scatter with time coloring
    t_days_syn = syn_data[:, 0] / 86400.0
    scatter_syn = ax_syn.scatter(syn_x, syn_y, c=t_days_syn, s=2, cmap='viridis', alpha=0.8)
    ax_syn.scatter([syn_x[0]], [syn_y[0]], c='lime', s=80, marker='o', zorder=5,
                   edgecolors='black', linewidths=1, label='Start')
    ax_syn.scatter([syn_x[-1]], [syn_y[-1]], c='red', s=80, marker='X', zorder=5,
                   edgecolors='black', linewidths=1, label='End')

    # Plot Lagrange points (only L1 for cleaner view)
    lp_colors = {'L1': 'red', 'L2': 'orange'}
    for name, coords in lagrange_points.items():
        if name in lp_colors:
            ax_syn.scatter([coords[0]/scale], [coords[1]/scale],
                          c=lp_colors[name], s=100, marker='*', zorder=4, label=name)

    # Plot Earth and Moon positions in synodic frame
    if 'primary' in lagrange_points:
        earth_pos = lagrange_points['primary']
        earth_syn = Circle((earth_pos[0]/scale, earth_pos[1]/scale),
                          earth_radius_mm * 3, color='blue', alpha=0.6)
        ax_syn.add_patch(earth_syn)
        ax_syn.scatter([earth_pos[0]/scale], [earth_pos[1]/scale],
                      c='blue', s=50, marker='o', zorder=4, label='Earth')

    if 'secondary' in lagrange_points:
        moon_pos = lagrange_points['secondary']
        moon_radius_mm = 1.737e6 / scale
        moon_syn = Circle((moon_pos[0]/scale, moon_pos[1]/scale),
                         moon_radius_mm * 10, color='gray', alpha=0.6)
        ax_syn.add_patch(moon_syn)
        ax_syn.scatter([moon_pos[0]/scale], [moon_pos[1]/scale],
                      c='gray', s=40, marker='o', zorder=4, label='Moon')

    # Mark barycenter
    ax_syn.scatter([0], [0], c='black', s=30, marker='+', zorder=4, label='Barycenter')

    ax_syn.set_xlabel(f'X ({scale_label})')
    ax_syn.set_ylabel(f'Y ({scale_label})')
    ax_syn.set_title(f'{title}\nSynodic (Earth-Moon Rotating)')
    ax_syn.grid(True, alpha=0.3)
    ax_syn.axis('equal')
    ax_syn.legend(fontsize=7, loc='best')

    # Add colorbar for synodic
    cbar_syn = plt.colorbar(scatter_syn, ax=ax_syn, shrink=0.8)
    cbar_syn.set_label('Time (days)', fontsize=8)


def main() -> int:
    print("Building and running orbitsim_example...")

    # Build the project
    try:
        subprocess.run(['cmake', '--build', 'build', '-j'], check=True, cwd='.')
    except subprocess.CalledProcessError:
        print("Error: Failed to build. Run: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release", file=sys.stderr)
        return 1

    # Run the example
    try:
        result = subprocess.run(['./build/orbitsim_example'],
                               capture_output=True, text=True, check=True, cwd='.')
    except subprocess.CalledProcessError as e:
        print(f"Error running orbitsim_example: {e}", file=sys.stderr)
        return 1

    # Parse output
    trajectories, info, lagrange_points = parse_output(result.stdout)

    print(f"Found trajectories: {list(trajectories.keys())}")
    print(f"Found info sections: {list(info.keys())}")
    print(f"Found Lagrange points: {list(lagrange_points.keys())}")

    # Create figure with 3 rows (elliptical, L1, circular) x 2 columns (ECI, synodic)
    fig, axes = plt.subplots(3, 2, figsize=(14, 18))

    # Row 0: Elliptical orbit
    if 'elliptical_eci' in trajectories and 'elliptical_synodic' in trajectories:
        orbit_info = info.get('elliptical_orbit_info', {})
        peri = orbit_info.get('periapsis_alt_m', 0) / 1000
        apo = orbit_info.get('apoapsis_alt_m', 0) / 1000
        ecc = orbit_info.get('eccentricity', 0)
        title = f'High Elliptical Orbit\n(Pe={peri:.0f} km, Ap={apo:.0f} km, e={ecc:.3f})'

        plot_comparison(axes[0, 0], axes[0, 1],
                       trajectories['elliptical_eci'],
                       trajectories['elliptical_synodic'],
                       lagrange_points, title, color='tab:blue')
    else:
        axes[0, 0].text(0.5, 0.5, 'elliptical data not found', ha='center', va='center')
        axes[0, 1].text(0.5, 0.5, 'elliptical data not found', ha='center', va='center')

    # Row 1: L1 spacecraft
    if 'l1_eci' in trajectories and 'l1_synodic' in trajectories:
        l1_info = info.get('l1_spacecraft_info', {})
        title = 'Spacecraft near L1 Lagrange Point\n(displaced ~5000 km from L1)'

        plot_comparison(axes[1, 0], axes[1, 1],
                       trajectories['l1_eci'],
                       trajectories['l1_synodic'],
                       lagrange_points, title, color='tab:orange')
    else:
        axes[1, 0].text(0.5, 0.5, 'L1 data not found', ha='center', va='center')
        axes[1, 1].text(0.5, 0.5, 'L1 data not found', ha='center', va='center')

    # Row 2: Circular orbit
    if 'circular_eci' in trajectories and 'circular_synodic' in trajectories:
        circ_info = info.get('circular_orbit_info', {})
        alt = circ_info.get('altitude_m', 0) / 1000
        period_days = circ_info.get('period_days', 0)
        title = f'Circular Orbit\n(alt={alt:.0f} km, period={period_days:.2f} days)'

        plot_comparison(axes[2, 0], axes[2, 1],
                       trajectories['circular_eci'],
                       trajectories['circular_synodic'],
                       lagrange_points, title, color='tab:green')
    else:
        axes[2, 0].text(0.5, 0.5, 'circular data not found', ha='center', va='center')
        axes[2, 1].text(0.5, 0.5, 'circular data not found', ha='center', va='center')

    # Overall title
    fig.suptitle('Coordinate Frame Comparison: ECI vs Synodic (Earth-Moon)\n'
                 'Same trajectory, different reference frames',
                 fontsize=14, fontweight='bold', y=0.995)

    # Add explanation text
    explanation = (
        "ECI (left): Non-rotating frame centered on Earth. Keplerian orbits appear as ellipses.\n"
        "Synodic (right): Rotating frame where Earth-Moon line is fixed. "
        "Lagrange points are stationary; simple orbits become complex curves."
    )
    fig.text(0.5, 0.01, explanation, ha='center', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout(rect=[0, 0.03, 1, 0.98])

    out_file = 'frame_comparison.png'
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {out_file}")
    plt.show()

    return 0


if __name__ == '__main__':
    sys.exit(main())
