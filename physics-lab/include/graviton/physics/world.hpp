#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include "graviton/math/vec2.hpp"
#include "graviton/physics/body.hpp"
#include "graviton/physics/spatial_hash_grid.hpp"

namespace graviton::physics {

/**
 * @brief Discrete simulation world managing integration, force accumulation, and contacts.
 */
class PhysicsWorld {
public:
    double G{1000.0};       // Gravitational constant
    double softening{5.0};  // Plummer softening length

    explicit PhysicsWorld(double cellSize = 30.0) : m_grid(cellSize) {}

    void addBody(const Body& b) {
        m_bodies.push_back(b);
    }

    [[nodiscard]] size_t bodyCount() const noexcept { return m_bodies.size(); }
    [[nodiscard]] const std::vector<Body>& getBodies() const noexcept { return m_bodies; }

    /**
     * @brief Computes N-body gravitational accelerations with Plummer softening:
     * a_i = sum G * m_j * (r_j - r_i) / (|r|^2 + eps^2)^(3/2)
     */
    void computeGravitationalAccelerations() {
        const size_t n = m_bodies.size();
        for (auto& b : m_bodies) {
            b.acceleration = math::Vec2{0.0, 0.0};
        }

        const double epsSq = softening * softening;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                math::Vec2 diff = m_bodies[j].position - m_bodies[i].position;
                double distSq = diff.lengthSq() + epsSq;
                double invDistCube = 1.0 / (distSq * std::sqrt(distSq));
                math::Vec2 force = diff * (G * invDistCube);

                m_bodies[i].acceleration += force * m_bodies[j].mass;
                m_bodies[j].acceleration -= force * m_bodies[i].mass;
            }
        }
    }

    /**
     * @brief Symplectic Velocity Verlet Integration (2nd Order O(dt^2))
     * Preserves Poincaré invariants in phase space -> strictly bounded energy drift.
     */
    void step(double dt) {
        const size_t n = m_bodies.size();
        if (n == 0) return;

        // 1. Position update + cache old acceleration: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
        for (auto& b : m_bodies) {
            b.position += b.velocity * dt + b.acceleration * (0.5 * dt * dt);
            b.accelerationPrev = b.acceleration;
        }

        // 2. Recalculate forces at new positions: a(t+dt)
        computeGravitationalAccelerations();

        // 3. Velocity update: v(t+dt) = v(t) + 0.5 * (a(t) + a(t+dt)) * dt
        for (auto& b : m_bodies) {
            b.velocity += (b.accelerationPrev + b.acceleration) * (0.5 * dt);
        }

        // 4. Broadphase & Narrowphase Collision Resolution
        resolveCollisionsWithGrid();
    }

    /**
     * @brief Broadphase culling using Uniform Spatial Hash Grid + Narrowphase Elastic Impulse
     */
    void resolveCollisionsWithGrid() {
        m_grid.clear();
        for (const auto& b : m_bodies) {
            m_grid.insert(b);
        }

        std::vector<std::pair<uint32_t, uint32_t>> testedPairs;

        for (const auto& [cellKey, bodyIds] : m_grid.getBuckets()) {
            if (bodyIds.size() < 2) continue;

            for (size_t i = 0; i < bodyIds.size(); ++i) {
                for (size_t j = i + 1; j < bodyIds.size(); ++j) {
                    uint32_t idA = bodyIds[i];
                    uint32_t idB = bodyIds[j];
                    if (idA > idB) std::swap(idA, idB);

                    auto pair = std::make_pair(idA, idB);
                    if (std::find(testedPairs.begin(), testedPairs.end(), pair) != testedPairs.end()) {
                        continue;
                    }
                    testedPairs.push_back(pair);

                    resolveContact(m_bodies[idA], m_bodies[idB]);
                }
            }
        }
    }

    /**
     * @brief Elastic Impulse Contact Resolution
     */
    static void resolveContact(Body& a, Body& b) {
        math::Vec2 delta = b.position - a.position;
        double distSq = delta.lengthSq();
        double minDist = a.radius + b.radius;

        if (distSq < minDist * minDist && distSq > 1e-8) {
            double dist = std::sqrt(distSq);
            math::Vec2 normal = delta / dist;

            // Positional separation to prevent sticky overlap penetration
            double overlap = minDist - dist;
            double totalMass = a.mass + b.mass;
            a.position -= normal * (overlap * (b.mass / totalMass));
            b.position += normal * (overlap * (a.mass / totalMass));

            // Relative contact velocity along normal
            math::Vec2 relVel = b.velocity - a.velocity;
            double velAlongNormal = relVel.dot(normal);

            if (velAlongNormal > 0.0) return; // Bodies separating

            // Impulse scalar: J = -(1 + e) * (v_rel . n) / (1/m_a + 1/m_b)
            double restitution = std::min(a.restitution, b.restitution);
            double impulseMagnitude = -(1.0 + restitution) * velAlongNormal / (1.0 / a.mass + 1.0 / b.mass);
            math::Vec2 impulse = normal * impulseMagnitude;

            a.velocity -= impulse / a.mass;
            b.velocity += impulse / b.mass;
        }
    }

    /**
     * @brief Hamiltonian Mechanical Energy: H = Kinetic (T) + Potential (V)
     */
    [[nodiscard]] double computeHamiltonianEnergy() const {
        double T = 0.0;
        double V = 0.0;
        const size_t n = m_bodies.size();
        const double epsSq = softening * softening;

        for (size_t i = 0; i < n; ++i) {
            T += 0.5 * m_bodies[i].mass * m_bodies[i].velocity.lengthSq();
            for (size_t j = i + 1; j < n; ++j) {
                double dist = std::sqrt((m_bodies[j].position - m_bodies[i].position).lengthSq() + epsSq);
                V -= (G * m_bodies[i].mass * m_bodies[j].mass) / dist;
            }
        }
        return T + V;
    }

private:
    std::vector<Body> m_bodies;
    SpatialHashGrid m_grid;
};

} // namespace graviton::physics
