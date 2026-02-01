#!/usr/bin/env python3
"""Plot orbits from orbitsim CSV output. Usage: python3 plot_orbits.py"""

import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt

R_EARTH = 6.371e6

def parse(lines):
    trajs = {}
    name, data = None, []
    for line in lines:
        line = line.strip()
        if not line or line.startswith('---'):
            continue
        if line.count(',') < 6:
            if name and data:
                trajs[name] = np.array(data)
            name, data = line, []
        else:
            vals = [float(x) for x in line.split(',')]
            if len(vals) == 7:
                data.append(vals)
    if name and data:
        trajs[name] = np.array(data)
    return trajs

def plot(trajs):
    fig = plt.figure(figsize=(12, 5))

    # 3D view
    ax1 = fig.add_subplot(121, projection='3d')
    u, v = np.linspace(0, 2*np.pi, 30), np.linspace(0, np.pi, 30)
    ax1.plot_surface(R_EARTH * np.outer(np.cos(u), np.sin(v)),
                     R_EARTH * np.outer(np.sin(u), np.sin(v)),
                     R_EARTH * np.outer(np.ones(30), np.cos(v)), color='blue', alpha=0.5)

    if 'moon_rel_earth' in trajs:
        m = trajs['moon_rel_earth']
        ax1.plot(m[:,1], m[:,2], m[:,3], 'gray', lw=0.5, alpha=0.5)

    if 'sc_rel_earth' in trajs:
        s = trajs['sc_rel_earth']
        ax1.plot(s[:,1], s[:,2], s[:,3], 'r', lw=1)
        ax1.scatter(*s[0,1:4], c='g', s=50, marker='o')
        ax1.scatter(*s[-1,1:4], c='r', s=50, marker='x')

    ax1.set_xlabel('X [m]')
    ax1.set_ylabel('Y [m]')
    ax1.set_zlabel('Z [m]')
    ax1.set_box_aspect([1,1,1])

    # 2D distance plot
    ax2 = fig.add_subplot(122)
    if 'sc_rel_earth' in trajs:
        s = trajs['sc_rel_earth']
        t = s[:,0] / 86400
        r = np.linalg.norm(s[:,1:4], axis=1)
        ax2.plot(t, r/1e6, 'r', lw=1)

    ax2.axhline(R_EARTH/1e6, color='b', ls=':', lw=1)
    ax2.set_xlabel('Time [days]')
    ax2.set_ylabel('Distance [Mm]')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('orbit.png', dpi=100)
    print("Saved: orbit.png", file=sys.stderr)
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
        output_lines = result.stdout.splitlines()
    except subprocess.CalledProcessError as e:
        print(f"Error running orbitsim_example: {e}")
        return 1

    # Parse and visualize
    trajs = parse(output_lines)
    print(f"Loaded: {list(trajs.keys())}")

    if 'sc_rel_earth' in trajs:
        plot(trajs)
    else:
        print("Error: No 'sc_rel_earth' trajectory found in output")
        return 1

    return 0

if __name__ == '__main__':
    sys.exit(main())
