#pragma once

#include "vec3.hpp"
#include "material.hpp"

namespace rt {

struct Hit {
    double t;
    Point3 point;
    Vec3 normal;
    MaterialPtr material;
    bool valid;

    Hit() : t(INF), point(0), normal(0), material(nullptr), valid(false) {}

    Hit(double t, const Point3& point, const Vec3& normal, MaterialPtr material)
        : t(t), point(point), normal(normal), material(material), valid(true) {}

    static Hit no_hit() {
        return Hit();
    }

    bool is_hit() const { return valid && t < INF; }

    operator bool() const { return is_hit(); }
};

// Global no_hit constant
inline const Hit NO_HIT = Hit::no_hit();

} // namespace rt
