#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>

#if __cplusplus < 201703L
namespace std {
template <class T>
constexpr const T& clamp(const T& v, const T& lo, const T& hi) {
    return v < lo ? lo : (hi < v ? hi : v);
}
}
#endif

namespace rt {

constexpr double EPSILON = 1e-8;
constexpr double INF = std::numeric_limits<double>::infinity();
constexpr double PI = 3.14159265358979323846;

struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double v) : x(v), y(v), z(v) {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
    double& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }

    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(const Vec3& v) const { return Vec3(x * v.x, y * v.y, z * v.z); }
    Vec3 operator/(const Vec3& v) const { return Vec3(x / v.x, y / v.y, z / v.z); }

    Vec3 operator*(double t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator/(double t) const { return Vec3(x / t, y / t, z / t); }

    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    Vec3& operator*=(const Vec3& v) { x *= v.x; y *= v.y; z *= v.z; return *this; }
    Vec3& operator*=(double t) { x *= t; y *= t; z *= t; return *this; }
    Vec3& operator/=(double t) { x /= t; y /= t; z /= t; return *this; }

    bool operator==(const Vec3& v) const { return x == v.x && y == v.y && z == v.z; }
    bool operator!=(const Vec3& v) const { return !(*this == v); }

    double length_squared() const { return x*x + y*y + z*z; }
    double length() const { return std::sqrt(length_squared()); }

    Vec3 normalized() const {
        double len = length();
        return len > EPSILON ? *this / len : Vec3(0);
    }

    bool near_zero() const {
        return (std::fabs(x) < EPSILON) && (std::fabs(y) < EPSILON) && (std::fabs(z) < EPSILON);
    }

    static Vec3 min(const Vec3& a, const Vec3& b) {
        return Vec3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
    }

    static Vec3 max(const Vec3& a, const Vec3& b) {
        return Vec3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
    }
};

inline Vec3 operator*(double t, const Vec3& v) { return v * t; }

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

inline Vec3 normalize(const Vec3& v) {
    return v.normalized();
}

inline Vec3 reflect(const Vec3& v, const Vec3& n) {
    return v - 2.0 * dot(v, n) * n;
}

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

// Type aliases
using Point3 = Vec3;
using Color = Vec3;

// Color utilities
inline Color clamp_color(const Color& c) {
    return Color(
        std::clamp(c.x, 0.0, 1.0),
        std::clamp(c.y, 0.0, 1.0),
        std::clamp(c.z, 0.0, 1.0)
    );
}

inline double to_srgb_component(double c) {
    c = std::clamp(c, 0.0, 1.0);
    return c > 0.0031308 ? 1.055 * std::pow(c, 1.0/2.4) - 0.055 : 12.92 * c;
}

inline Color to_srgb(const Color& c) {
    return Color(
        to_srgb_component(c.x),
        to_srgb_component(c.y),
        to_srgb_component(c.z)
    );
}

/** Struct to represent Mat3s for transformations */
struct Mat3 {
    double m[3][3];

    Mat3() {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m[i][j] = (i == j) ? 1.0 : 0.0;
    }

    Mat3(double m00, double m01, double m02,
         double m10, double m11, double m12,
         double m20, double m21, double m22) {
        m[0][0] = m00; m[0][1] = m01; m[0][2] = m02;
        m[1][0] = m10; m[1][1] = m11; m[1][2] = m12;
        m[2][0] = m20; m[2][1] = m21; m[2][2] = m22;
    }

    Vec3 operator*(const Vec3& v) const {
        return Vec3(
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z
        );
    }

    Mat3 operator*(const Mat3& other) const {
        Mat3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m[i][j] = 0;
                for (int k = 0; k < 3; ++k) {
                    result.m[i][j] += m[i][k] * other.m[k][j];
                }
            }
        }
        return result;
    }

    static Mat3 identity() { return Mat3(); }

    static Mat3 scale(double s) {
        return Mat3(s, 0, 0,  0, s, 0,  0, 0, s);
    }

    static Mat3 scale(double sx, double sy, double sz) {
        return Mat3(sx, 0, 0,  0, sy, 0,  0, 0, sz);
    }

    static Mat3 rotate_x(double angle_deg) {
        double r = angle_deg * PI / 180.0;
        double c = std::cos(r), s = std::sin(r);
        return Mat3(1, 0, 0,  0, c, -s,  0, s, c);
    }

    static Mat3 rotate_y(double angle_deg) {
        double r = angle_deg * PI / 180.0;
        double c = std::cos(r), s = std::sin(r);
        return Mat3(c, 0, s,  0, 1, 0,  -s, 0, c);
    }

    static Mat3 rotate_z(double angle_deg) {
        double r = angle_deg * PI / 180.0;
        double c = std::cos(r), s = std::sin(r);
        return Mat3(c, -s, 0,  s, c, 0,  0, 0, 1);
    }

    static Mat3 rotate(double rx, double ry, double rz) {
        return rotate_z(rz) * rotate_y(ry) * rotate_x(rx);
    }

    // Transpose (useful for transforming normals)
    Mat3 transposed() const {
        return Mat3(
            m[0][0], m[1][0], m[2][0],
            m[0][1], m[1][1], m[2][1],
            m[0][2], m[1][2], m[2][2]
        );
    }
};

} // namespace rt
