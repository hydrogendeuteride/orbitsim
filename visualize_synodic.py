#!/usr/bin/env python3
"""
Visualize orbitsim trajectories in the Earth-Moon synodic (rotating) frame.
Reads CSV output from orbitsim_example and plots:
  - Spacecraft trajectory in rotating coordinates
  - Lagrange points (L1-L5)
  - Primary (Earth) and Secondary (Moon) positions
"""

import sys
import subprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np


def parse_output(output):
    """Parse the orbitsim_example output into dictionaries."""
    lines = output.strip().split('\n')
    data = {}
    current_key = None
    current_data = []

    lagrange_points = {}
    synodic_params = {}

    for line in lines:
        line = line.strip()
        if not line:
            continue

        # Check for section headers
        if line in ['sc_synodic', 'lagrange_points', 'synodic_frame_t0']:
            if current_key and current_data:
                data[current_key] = np.array(current_data)
            current_key = line
            current_data = []
            continue

        # Parse Lagrange points
        if current_key == 'lagrange_points':
            parts = line.split(',')
            if len(parts) == 4:
                name = parts[0]
                coords = [float(x) for x in parts[1:4]]
                lagrange_points[name] = coords
            continue

        # Parse synodic frame parameters
        if current_key == 'synodic_frame_t0':
            parts = line.split(',')
            if parts[0] == 'separation_m':
                synodic_params['separation_m'] = float(parts[1])
            elif parts[0] == 'mu':
                synodic_params['mu'] = float(parts[1])
            elif parts[0] == 'omega_radps':
                synodic_params['omega'] = [float(x) for x in parts[1:4]]
            continue

        # Parse trajectory data
        if current_key == 'sc_synodic':
            try:
                values = [float(x) for x in line.split(',')]
                if len(values) == 7:
                    current_data.append(values)
            except ValueError:
                continue

    # Save the last section
    if current_key and current_data:
        data[current_key] = np.array(current_data)

    return data, lagrange_points, synodic_params


def plot_synodic_frame(data, lagrange_points, synodic_params):
    """Create a 2D visualization of the synodic frame."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Get spacecraft trajectory in synodic frame
    if 'sc_synodic' not in data:
        print("Error: sc_synodic trajectory not found in output")
        return

    sc_synodic = data['sc_synodic']

    # Convert to km for better readability
    x_km = sc_synodic[:, 1] / 1000.0
    y_km = sc_synodic[:, 2] / 1000.0
    z_km = sc_synodic[:, 3] / 1000.0

    # Plot 1: XY plane (orbital plane view)
    ax1.plot(x_km, y_km, 'b-', linewidth=0.8, label='Spacecraft', alpha=0.7)
    ax1.plot(x_km[0], y_km[0], 'go', markersize=8, label='Start')
    ax1.plot(x_km[-1], y_km[-1], 'ro', markersize=8, label='End')

    # Plot Lagrange points
    colors = {'L1': 'red', 'L2': 'orange', 'L3': 'purple', 'L4': 'cyan', 'L5': 'magenta'}
    for name, coords in lagrange_points.items():
        if name in colors:
            ax1.plot(coords[0]/1000.0, coords[1]/1000.0, 'o',
                    color=colors[name], markersize=10, label=name)
            ax1.text(coords[0]/1000.0, coords[1]/1000.0 + 10000, name,
                    ha='center', fontsize=10, fontweight='bold')

    # Plot primary (Earth) and secondary (Moon)
    if 'primary' in lagrange_points and 'secondary' in lagrange_points:
        earth_pos = lagrange_points['primary']
        moon_pos = lagrange_points['secondary']

        # Draw Earth (larger blue circle)
        earth_radius_km = 6371.0
        earth_circle = Circle((earth_pos[0]/1000.0, earth_pos[1]/1000.0),
                             earth_radius_km, color='blue', alpha=0.6, label='Earth')
        ax1.add_patch(earth_circle)

        # Draw Moon (smaller gray circle)
        moon_radius_km = 1737.4
        moon_circle = Circle((moon_pos[0]/1000.0, moon_pos[1]/1000.0),
                            moon_radius_km, color='gray', alpha=0.6, label='Moon')
        ax1.add_patch(moon_circle)

    ax1.set_xlabel('X (km)', fontsize=12)
    ax1.set_ylabel('Y (km)', fontsize=12)
    ax1.set_title('Earth-Moon Synodic Frame (XY Plane)', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axis('equal')
    ax1.legend(loc='best', fontsize=9)

    # Plot 2: XZ plane (out-of-plane view)
    ax2.plot(x_km, z_km, 'b-', linewidth=0.8, label='Spacecraft', alpha=0.7)
    ax2.plot(x_km[0], z_km[0], 'go', markersize=8, label='Start')
    ax2.plot(x_km[-1], z_km[-1], 'ro', markersize=8, label='End')

    # Plot Lagrange points in XZ
    for name, coords in lagrange_points.items():
        if name in colors:
            ax2.plot(coords[0]/1000.0, coords[2]/1000.0, 'o',
                    color=colors[name], markersize=10, label=name)

    # Plot primary and secondary in XZ
    if 'primary' in lagrange_points and 'secondary' in lagrange_points:
        earth_pos = lagrange_points['primary']
        moon_pos = lagrange_points['secondary']

        earth_circle = Circle((earth_pos[0]/1000.0, earth_pos[2]/1000.0),
                             earth_radius_km, color='blue', alpha=0.6, label='Earth')
        ax2.add_patch(earth_circle)

        moon_circle = Circle((moon_pos[0]/1000.0, moon_pos[2]/1000.0),
                            moon_radius_km, color='gray', alpha=0.6, label='Moon')
        ax2.add_patch(moon_circle)

    ax2.set_xlabel('X (km)', fontsize=12)
    ax2.set_ylabel('Z (km)', fontsize=12)
    ax2.set_title('Earth-Moon Synodic Frame (XZ Plane)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.axis('equal')
    ax2.legend(loc='best', fontsize=9)

    # Add info text
    if synodic_params:
        info_text = f"Separation: {synodic_params.get('separation_m', 0)/1e6:.2f} Mm\n"
        info_text += f"μ: {synodic_params.get('mu', 0):.6f}\n"
        if 'omega' in synodic_params:
            omega_mag = np.linalg.norm(synodic_params['omega'])
            info_text += f"ω: {omega_mag:.6e} rad/s"
        fig.text(0.5, 0.02, info_text, ha='center', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig('synodic_frame.png', dpi=150, bbox_inches='tight')
    print(f"Saved visualization to synodic_frame.png")
    plt.show()  # Comment out for headless environments


def main():
    print("Building and running orbitsim_example...")

    # Build the project
    try:
        subprocess.run(['cmake', '--build', 'build', '-j'], check=True, cwd='.')
    except subprocess.CalledProcessError:
        print("Error: Failed to build orbitsim_example")
        print("Make sure to run: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release")
        return 1

    # Run the example and capture output
    try:
        result = subprocess.run(['./build/orbitsim_example'],
                               capture_output=True, text=True, check=True, cwd='.')
        output = result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running orbitsim_example: {e}")
        return 1

    # Parse and visualize
    data, lagrange_points, synodic_params = parse_output(output)

    if not data:
        print("Error: No data parsed from output")
        return 1

    print(f"Parsed {len(data)} trajectory sections")
    print(f"Found {len(lagrange_points)} Lagrange points")

    plot_synodic_frame(data, lagrange_points, synodic_params)

    return 0


if __name__ == '__main__':
    sys.exit(main())
