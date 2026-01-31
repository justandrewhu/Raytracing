#pragma once

#include "vec3.hpp"

namespace rt {

struct Ray {
    Point3 origin;
    Vec3 direction;
    double t_min;
    double t_max;

    Ray() : origin(0), direction(0, 0, -1), t_min(0), t_max(INF) {}

    Ray(const Point3& origin, const Vec3& direction, double t_min = 0.0, double t_max = INF)
        : origin(origin), direction(direction), t_min(t_min), t_max(t_max) {}

    Point3 at(double t) const {
        return origin + t * direction;
    }
};

} // namespace rt
