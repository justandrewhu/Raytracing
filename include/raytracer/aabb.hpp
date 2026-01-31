#pragma once

#include "vec3.hpp"
#include "ray.hpp"
#include <algorithm>

namespace rt {

struct AABB {
    Point3 min_corner;
    Point3 max_corner;

    AABB() : min_corner(INF), max_corner(-INF) {}

    AABB(const Point3& min_corner, const Point3& max_corner)
        : min_corner(min_corner), max_corner(max_corner) {}

    static AABB from_points(const Point3* points, size_t count) {
        Point3 min_c(INF);
        Point3 max_c(-INF);
        for (size_t i = 0; i < count; ++i) {
            min_c = Vec3::min(min_c, points[i]);
            max_c = Vec3::max(max_c, points[i]);
        }
        return AABB(min_c, max_c);
    }

    static AABB merge(const AABB& a, const AABB& b) {
        return AABB(
            Vec3::min(a.min_corner, b.min_corner),
            Vec3::max(a.max_corner, b.max_corner)
        );
    }

    Point3 centroid() const {
        return (min_corner + max_corner) * 0.5;
    }

    Vec3 extent() const {
        return max_corner - min_corner;
    }

    int longest_axis() const {
        Vec3 e = extent();
        if (e.x > e.y && e.x > e.z) return 0;
        if (e.y > e.z) return 1;
        return 2;
    }

    double surface_area() const {
        Vec3 e = extent();
        return 2.0 * (e.x * e.y + e.y * e.z + e.z * e.x);
    }

    bool intersect(const Ray& ray, double& t_near, double& t_far) const {
        t_near = ray.t_min;
        t_far = ray.t_max;

        for (int axis = 0; axis < 3; ++axis) {
            double dir_comp = ray.direction[axis];
            double origin_comp = ray.origin[axis];
            double min_bound = min_corner[axis];
            double max_bound = max_corner[axis];

            if (std::fabs(dir_comp) < 1e-12) {
                // Ray parallel to this slab
                if (origin_comp < min_bound || origin_comp > max_bound) {
                    return false;
                }
            } else {
                double inv_dir = 1.0 / dir_comp;
                double t_to_min = (min_bound - origin_comp) * inv_dir;
                double t_to_max = (max_bound - origin_comp) * inv_dir;

                if (t_to_min > t_to_max) {
                    std::swap(t_to_min, t_to_max);
                }

                t_near = std::max(t_near, t_to_min);
                t_far = std::min(t_far, t_to_max);

                if (t_far < t_near) {
                    return false;
                }
            }
        }
        return true;
    }

    bool intersect(const Ray& ray) const {
        double t_near, t_far;
        return intersect(ray, t_near, t_far);
    }
};

} // namespace rt
