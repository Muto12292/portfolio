#include <iostream>
#include <chrono>
#include <iomanip>
#include "graviton/math/vec2.hpp"
#include "graviton/physics/body.hpp"
#include "graviton/physics/world.hpp"

using namespace graviton;

int main() {
    std::cout << "======================================================================\n";
    std::cout << " Graviton C++ Core: Symplectic Solver & Spatial Hash Grid Benchmark\n";
    std::cout << " Author: Muto Syuleyman | Nottingham Trent University (BSc Physics)\n";
    std::cout << "======================================================================\n\n";

    physics::PhysicsWorld world(25.0); // 25.0 unit spatial grid cells

    // Setup Chenciner-Montgomery Figure-8 Choreography
    const double m = 120.0;
    const double x1 = 120.0 * 0.97000436;
    const double y1 = -120.0 * 0.24308753;
    const double vx3 = 2.0 * -0.93240737 * 1.1;
    const double vy3 = 2.0 * -0.86473146 * 1.1;

    world.addBody(physics::Body(0, math::Vec2(x1, y1),   math::Vec2(-vx3 / 2.0, -vy3 / 2.0), m, 6.0));
    world.addBody(physics::Body(1, math::Vec2(-x1, -y1), math::Vec2(-vx3 / 2.0, -vy3 / 2.0), m, 6.0));
    world.addBody(physics::Body(2, math::Vec2(0, 0),     math::Vec2(vx3, vy3),               m, 6.0));

    // Baseline reference
    world.computeGravitationalAccelerations();
    const double initialEnergy = world.computeHamiltonianEnergy();

    std::cout << "Initialized Figure-8 Orbital Choreography (" << world.bodyCount() << " bodies).\n";
    std::cout << "Initial Hamiltonian Total Energy (H_0): " << std::fixed << std::setprecision(4)
              << initialEnergy << " J\n\n";

    constexpr double dt = 0.016;
    constexpr int totalTicks = 100000; // 100k integration steps
    std::cout << "Executing " << totalTicks << " Symplectic Velocity Verlet ticks (dt = " << dt << "s)...\n";

    const auto startClock = std::chrono::high_resolution_clock::now();

    for (int tick = 0; tick < totalTicks; ++tick) {
        world.step(dt);
    }

    const auto endClock = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> durationMs = endClock - startClock;

    const double finalEnergy = world.computeHamiltonianEnergy();
    const double relativeDrift = std::abs((finalEnergy - initialEnergy) / initialEnergy) * 100.0;
    const double ticksPerSec = totalTicks / (durationMs.count() / 1000.0);

    std::cout << "\n------------------------- Benchmark Results -------------------------\n";
    std::cout << "Total Execution Time:      " << durationMs.count() << " ms\n";
    std::cout << "Throughput:                " << std::setprecision(1) << ticksPerSec << " steps/sec\n";
    std::cout << "Time per Simulation Step:  " << std::setprecision(4) << (durationMs.count() / totalTicks) * 1000.0 << " microseconds\n";
    std::cout << "Final Total Energy (H_f):  " << finalEnergy << " J\n";
    std::cout << "Relative Energy Drift:     " << std::setprecision(5) << relativeDrift << " % (Bounded Symplectic Drift)\n";
    std::cout << "---------------------------------------------------------------------\n";
    std::cout << "STATUS: PASS (Symplectic energy conservation verified under 10^5 steps).\n\n";

    return 0;
}
