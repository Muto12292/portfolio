#pragma once

#include <cstdint>
#include "graviton/math/vec2.hpp"

namespace graviton::physics {

/**
 * @brief Discrete rigid particle with state and material attributes.
 */
struct Body {
    uint32_t id{0};
    math::Vec2 position{0.0, 0.0};
    math::Vec2 velocity{0.0, 0.0};
    math::Vec2 acceleration{0.0, 0.0};
    math::Vec2 accelerationPrev{0.0, 0.0};
    double mass{1.0};
    double radius{5.0};
    double restitution{0.8}; // Coefficient of restitution [0.0 = inelastic, 1.0 = elastic]

    constexpr Body() = default;

    Body(uint32_t id_, math::Vec2 pos, math::Vec2 vel, double mass_, double radius_, double rest_ = 0.8)
        : id(id_), position(pos), velocity(vel), mass(mass_), radius(radius_), restitution(rest_) {}
};

} // namespace graviton::physics
