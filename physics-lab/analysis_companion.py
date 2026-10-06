"""
Graviton Engine: Numerical Integrator Benchmark & Hamiltonian Drift Analysis
Author: Muto Syuleyman (BSc Physics, Nottingham Trent University)
GitHub: @muto12292 | Email: mutoerduhan05@gmail.com

This module benchmarks three classical numerical integrators on the Gravitational Three-Body Problem:
1. Explicit Euler (1st order, non-symplectic)
2. Velocity Verlet (2nd order, symplectic, phase-space preserving)
3. Runge-Kutta 4th Order (4th order, high-precision classical integrator)

Demonstrates Hamiltonian Energy conservation: H = T + V
and long-term secular energy drift Delta E / E_0.
"""

import numpy as np
import time

G = 1000.0
EPSILON = 5.0  # Plummer softening parameter


def compute_accelerations(positions, masses):
    """
    Computes Newtonian gravitational accelerations with Plummer softening:
    a_i = sum_{j != i} G * m_j * (r_j - r_i) / (|r_j - r_i|^2 + eps^2)^(3/2)
    """
    n = len(masses)
    acc = np.zeros_like(positions)
    for i in range(n):
        for j in range(i + 1, n):
            diff = positions[j] - positions[i]
            dist_sq = np.sum(diff ** 2) + EPSILON ** 2
            inv_dist_cube = 1.0 / (dist_sq * np.sqrt(dist_sq))
            force_ij = G * diff * inv_dist_cube
            acc[i] += force_ij * masses[j]
            acc[j] -= force_ij * masses[i]
    return acc


def compute_energy(positions, velocities, masses):
    """
    Calculates Hamiltonian total mechanical energy:
    H = Kinetic (T) + Potential (V)
    """
    n = len(masses)
    # Kinetic energy: sum 0.5 * m_i * v_i^2
    T = 0.5 * np.sum(masses[:, None] * (velocities ** 2))

    # Potential energy: - sum_{i < j} G * m_i * m_j / sqrt(r_ij^2 + eps^2)
    V = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            dist = np.sqrt(np.sum((positions[j] - positions[i]) ** 2) + EPSILON ** 2)
            V -= (G * masses[i] * masses[j]) / dist
    return T + V


def simulate(integrator_name, init_pos, init_vel, masses, dt, total_steps):
    """Runs numerical integration and records relative energy drift."""
    pos = np.copy(init_pos)
    vel = np.copy(init_vel)
    n = len(masses)

    e0 = compute_energy(pos, vel, masses)
    drift_history = np.zeros(total_steps)

    for step in range(total_steps):
        if integrator_name == 'euler':
            # Explicit Euler
            acc = compute_accelerations(pos, masses)
            pos += vel * dt
            vel += acc * dt

        elif integrator_name == 'verlet':
            # Velocity Verlet (Symplectic)
            acc_old = compute_accelerations(pos, masses)
            pos += vel * dt + 0.5 * acc_old * (dt ** 2)
            acc_new = compute_accelerations(pos, masses)
            vel += 0.5 * (acc_old + acc_new) * dt

        elif integrator_name == 'rk4':
            # Runge-Kutta 4th Order
            # State vector: [pos, vel]
            k1_v = compute_accelerations(pos, masses)
            k1_x = vel

            pos_k2 = pos + 0.5 * dt * k1_x
            vel_k2 = vel + 0.5 * dt * k1_v
            k2_v = compute_accelerations(pos_k2, masses)
            k2_x = vel_k2

            pos_k3 = pos + 0.5 * dt * k2_x
            vel_k3 = vel + 0.5 * dt * k2_v
            k3_v = compute_accelerations(pos_k3, masses)
            k3_x = vel_k3

            pos_k4 = pos + dt * k3_x
            vel_k4 = vel + dt * k3_v
            k4_v = compute_accelerations(pos_k4, masses)
            k4_x = vel_k4

            pos += (dt / 6.0) * (k1_x + 2 * k2_x + 2 * k3_x + k4_x)
            vel += (dt / 6.0) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)

        # Track relative energy error: (E(t) - E_0) / |E_0|
        current_e = compute_energy(pos, vel, masses)
        drift_history[step] = (current_e - e0) / abs(e0)

    return drift_history


def main():
    print("=" * 65)
    print(" AetherPhys: Numerical Integrator Benchmark & Symplectic Analysis")
    print("=" * 65)

    # Initial condition: Chenciner-Montgomery Figure-8 Choreography
    masses = np.array([150.0, 150.0, 150.0])
    x1 = 150.0 * 0.97000436
    y1 = -150.0 * 0.24308753
    vx3 = 2.0 * -0.93240737 * 1.05
    vy3 = 2.0 * -0.86473146 * 1.05

    init_pos = np.array([
        [x1, y1],
        [-x1, -y1],
        [0.0, 0.0]
    ])
    init_vel = np.array([
        [-vx3 / 2, -vy3 / 2],
        [-vx3 / 2, -vy3 / 2],
        [vx3, vy3]
    ])

    dt = 0.015
    steps = 4000
    time_array = np.arange(steps) * dt

    integrators = ['euler', 'verlet', 'rk4']
    results = {}

    for name in integrators:
        start_t = time.perf_counter()
        drift = simulate(name, init_pos, init_vel, masses, dt, steps)
        elapsed = time.perf_counter() - start_t
        results[name] = drift
        print(f"[{name.upper():<6}] Max Relative Drift: {np.max(np.abs(drift)) * 100:.4f}% | Execution: {elapsed:.3f}s")

    print("\nBenchmark Summary:")
    print("-" * 55)
    for name in integrators:
        print(f"{name.upper():<8}: Max Drift = {np.max(np.abs(results[name])) * 100:.5f}% | Final Drift = {results[name][-1] * 100:.5f}%")

    try:
        import matplotlib.pyplot as plt
        print("\nPlotting energy drift comparative curves with Matplotlib...")
        plt.figure(figsize=(10, 5))
        plt.plot(time_array, results['euler'] * 100, label='Explicit Euler (O(h)) - Secular Divergence', color='crimson')
        plt.plot(time_array, results['verlet'] * 100, label='Velocity Verlet (Symplectic O(h²)) - Bounded Drift', color='teal', linewidth=1.8)
        plt.plot(time_array, results['rk4'] * 100, label='Runge-Kutta 4th Order (O(h⁴)) - Minimal Truncation Error', color='royalblue')

        plt.title('Hamiltonian Energy Conservation Error: $\\Delta E / |E_0|$ (%)', fontsize=13)
        plt.xlabel('Simulation Time (seconds)', fontsize=11)
        plt.ylabel('Relative Energy Drift (%)', fontsize=11)
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig('integrator_benchmark_drift.png', dpi=200)
        print("Benchmark figure saved as 'integrator_benchmark_drift.png'.")
    except ImportError:
        print("\n[Notice] matplotlib is not installed. Exporting benchmark results to 'integrator_benchmark.csv'...")
        with open('integrator_benchmark.csv', 'w') as f:
            f.write("time_sec,euler_drift_pct,verlet_drift_pct,rk4_drift_pct\n")
            for i in range(steps):
                f.write(f"{time_array[i]:.4f},{results['euler'][i]*100:.6f},{results['verlet'][i]*100:.6f},{results['rk4'][i]*100:.6f}\n")
        print("Benchmark data saved successfully to 'integrator_benchmark.csv'.")


if __name__ == '__main__':
    main()
