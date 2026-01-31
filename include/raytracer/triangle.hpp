#pragma once

#include "shape.hpp"
#include <array>

namespace rt {

class Triangle : public Shape {
public:
    std::array<Point3, 3> vertices;
    MaterialPtr material;

    Triangle(const Point3& v0, const Point3& v1, const Point3& v2, MaterialPtr material)
        : vertices({v0, v1, v2}), material(material) {}

    Triangle(const std::array<Point3, 3>& vs, MaterialPtr material)
        : vertices(vs), material(material) {}

    Vec3 normal() const {
        Vec3 edge1 = vertices[1] - vertices[0];
        Vec3 edge2 = vertices[2] - vertices[0];
        return normalize(cross(edge1, edge2));
    }

    Hit intersect(const Ray& ray) const override {
        // Moller-Trumbore intersection algorithm
        const Point3& a = vertices[0];
        const Point3& b = vertices[1];
        const Point3& c = vertices[2];

        Vec3 edge1 = b - a;
        Vec3 edge2 = c - a;
        Vec3 h = cross(ray.direction, edge2);
        double det = dot(edge1, h);

        // Check if ray is parallel to triangle
        if (std::fabs(det) < EPSILON) {
            return NO_HIT;
        }

        double inv_det = 1.0 / det;
        Vec3 s = ray.origin - a;
        double u = inv_det * dot(s, h);

        if (u < 0.0 || u > 1.0) {
            return NO_HIT;
        }

        Vec3 q = cross(s, edge1);
        double v = inv_det * dot(ray.direction, q);

        if (v < 0.0 || u + v > 1.0) {
            return NO_HIT;
        }

        double t = inv_det * dot(edge2, q);

        if (t < ray.t_min || t > ray.t_max) {
            return NO_HIT;
        }

        Point3 hit_point = ray.at(t);
        Vec3 n = normal();

        return Hit(t, hit_point, n, material);
    }

    AABB bounds() const override {
        return AABB::from_points(vertices.data(), 3);
    }
};

} // namespace rt
