#pragma once

#include <cmath>

namespace graviton::math {

/**
 * @brief 2D Vector structure with contiguous memory layout and SIMD-friendly alignment.
 */
struct Vec2 {
    double x{0.0};
    double y{0.0};

    constexpr Vec2() = default;
    constexpr Vec2(double x_, double y_) : x(x_), y(y_) {}

    constexpr Vec2 operator+(const Vec2& rhs) const { return {x + rhs.x, y + rhs.y}; }
    constexpr Vec2 operator-(const Vec2& rhs) const { return {x - rhs.x, y - rhs.y}; }
    constexpr Vec2 operator*(double s) const { return {x * s, y * s}; }
    constexpr Vec2 operator/(double s) const { return {x / s, y / s}; }

    Vec2& operator+=(const Vec2& rhs) { x += rhs.x; y += rhs.y; return *this; }
    Vec2& operator-=(const Vec2& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
    Vec2& operator*=(double s) { x *= s; y *= s; return *this; }

    [[nodiscard]] constexpr double dot(const Vec2& rhs) const { return x * rhs.x + y * rhs.y; }
    [[nodiscard]] constexpr double lengthSq() const { return x * x + y * y; }
    [[nodiscard]] double length() const { return std::sqrt(lengthSq()); }

    [[nodiscard]] Vec2 normalized() const {
        double len = length();
        return (len > 1e-9) ? (*this / len) : Vec2{0.0, 0.0};
    }
};

} // namespace graviton::math
