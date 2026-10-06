#pragma once

#include <vector>
#include <unordered_map>
#include <cmath>
#include "graviton/physics/body.hpp"

namespace graviton::physics {

/**
 * @brief Uniform Spatial Hash Grid for Broadphase Collision Culling: O(N^2) -> O(N).
 * Partitions world space into discrete uniform buckets to cull non-intersecting candidate pairs.
 */
class SpatialHashGrid {
public:
    explicit SpatialHashGrid(double cellSize) : m_cellSize(cellSize) {}

    void clear() {
        m_buckets.clear();
    }

    void insert(const Body& body) {
        int64_t minX = static_cast<int64_t>(std::floor((body.position.x - body.radius) / m_cellSize));
        int64_t maxX = static_cast<int64_t>(std::floor((body.position.x + body.radius) / m_cellSize));
        int64_t minY = static_cast<int64_t>(std::floor((body.position.y - body.radius) / m_cellSize));
        int64_t maxY = static_cast<int64_t>(std::floor((body.position.y + body.radius) / m_cellSize));

        for (int64_t x = minX; x <= maxX; ++x) {
            for (int64_t y = minY; y <= maxY; ++y) {
                uint64_t key = hashKey(x, y);
                m_buckets[key].push_back(body.id);
            }
        }
    }

    [[nodiscard]] const std::unordered_map<uint64_t, std::vector<uint32_t>>& getBuckets() const {
        return m_buckets;
    }

private:
    double m_cellSize;
    std::unordered_map<uint64_t, std::vector<uint32_t>> m_buckets;

    // Fast 64-bit spatial coordinate hash
    [[nodiscard]] static constexpr uint64_t hashKey(int64_t x, int64_t y) {
        constexpr uint64_t p1 = 73856093ULL;
        constexpr uint64_t p2 = 19349663ULL;
        return (static_cast<uint64_t>(x) * p1) ^ (static_cast<uint64_t>(y) * p2);
    }
};

} // namespace graviton::physics
