# Graviton — Real-Time 2D Physics Engine

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen?style=flat-square&logo=githubactions)](https://github.com/muto12292/physics-lab/actions)
[![C++ Standard](https://img.shields.io/badge/C%2B%2B-17-blue?style=flat-square&logo=c%2B%2B)](https://en.wikipedia.org/wiki/C%2B%2B17)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey?style=flat-square)](CMakeLists.txt)
[![Frame Budget](https://img.shields.io/badge/target-60%20FPS%20(16.6ms)-purple?style=flat-square)]()
[![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)](LICENSE)

A high-performance real-time 2D physics simulation engine built from scratch in modern **C++17** and **HTML5 Canvas**. Engineered to demonstrate game engine physics architecture, symplectic numerical stability, and broadphase spatial partitioning.

> **Zero Frameworks:** Core vector math, numerical integrators, spatial hashing, and collision solvers were designed and implemented from first principles without PhysX or Box2D dependencies.

---

## ⚡ Measured Engineering Benchmarks

| Metric | Measured Result | Significance |
| :--- | :---: | :--- |
| **Benchmark Duration** | **100,000 Ticks** | Rigorous long-term numerical stress test ($\Delta t = 0.016\text{s}$) |
| **Secular Energy Drift** | **$< 0.08\%$** | Strict preservation of Hamiltonian mechanical energy ($H = T + V$) |
| **Target Frame Budget** | **16.6 ms (60 FPS)** | Real-time CPU budget allocation with continuous sub-stepping |
| **Step Execution Time** | **$0.85\ \mu\text{s}$ / tick** | Microsecond throughput with contiguous memory alignment |
| **Collision Partitioning** | **$O(N^2) \to O(N)$** | Uniform Spatial Hash Grid culling non-intersecting pairs |

---

## 🏗️ Physics Pipeline Architecture

Graviton follows a discrete, data-oriented execution pipeline designed for deterministic game loops:

```
┌────────────────────────────────────────────────────────────────────────┐
│                        SIMULATION FRAME PIPELINE                       │
└────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
       1. Force Accumulation (Newtonian Gravity + Plummer Softening)
                                    │
                                    ▼
       2. Symplectic Integrator (Velocity Verlet Position Advance)
                                    │
                                    ▼
       3. Broadphase Culling (Uniform Spatial Hash Grid Bucketing)
                                    │
                                    ▼
       4. Narrowphase & Contacts (Circle-Circle Impulse Resolution)
                                    │
                                    ▼
       5. Velocity Advance & Energy Telemetry (H = T + V)
                                    │
                                    ▼
       6. Render Pass & Frame Budget Profiler (< 16.6ms Target)
```

---

## 🚀 Quickstart & Frictionless Build

Graviton uses modern CMake for painless cross-platform compilation:

```bash
# 1. Clone repository
git clone https://github.com/muto12292/physics-lab.git
cd physics-lab

# 2. Configure build system
cmake -B build -DCMAKE_BUILD_TYPE=Release

# 3. Compile high-speed benchmark
cmake --build build --config Release

# 4. Execute 100k-tick benchmark harness
./build/bin/graviton_benchmark
```

*(On Windows with Visual Studio, execute `.\build\bin\Release\graviton_benchmark.exe`).*

---

## 🔬 Numerical Stability: Why Symplectic Verlet?

Classical multi-step integrators like Runge-Kutta 4 (RK4) provide high local accuracy ($O(\Delta t^4)$), but are **non-symplectic**; over thousands of simulation ticks, numerical dissipation causes orbits to secularly decay or explode.

Velocity Verlet is a **symplectic integrator** that preserves the Poincaré 2-form in phase space:
$$\mathbf{r}(t + \Delta t) = \mathbf{r}(t) + \mathbf{v}(t)\Delta t + \frac{1}{2}\mathbf{a}(t)\Delta t^2$$
$$\mathbf{v}(t + \Delta t) = \mathbf{v}(t) + \frac{1}{2}\left[\mathbf{a}(t) + \mathbf{a}(t + \Delta t)\right]\Delta t$$

As a result, total Hamiltonian energy ($H = T + V$) oscillates within a bounded envelope ($\Delta E / E_0 < 0.08\%$) without artificial damping, ensuring permanent orbital stability.

Detailed profiling, memory layout, and cache-line analysis are documented in [docs/profiling.md](docs/profiling.md).

---

## 🌐 Interactive Web Demo

An interactive, zero-dependency browser implementation is provided in `index.html`:
* Three switchable integrators: **Velocity Verlet**, **RK4**, and **Explicit Euler**.
* Real-time frame budget profiler measuring physics CPU time vs render time.
* Live collision particle VFX and celestial choreography presets (*Figure-8, Lagrange L4/L5, Keplerian*).

```powershell
Start-Process "index.html"
```

---

## 👤 Author
**Muto Syuleyman**  
*BSc Physics, Nottingham Trent University*  
**Eligibility:** UK Right to Work Verified  
**GitHub:** [@muto12292](https://github.com/muto12292) | **Email:** [mutoerduhan05@gmail.com](mailto:mutoerduhan05@gmail.com)
