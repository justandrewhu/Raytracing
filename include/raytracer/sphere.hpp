#pragma once

#include "shape.hpp"
#include <cmath>

namespace rt {

class Sphere : public Shape {
public:
    Point3 center;
    double radius;
    MaterialPtr material;

    Sphere(const Point3& center, double radius, MaterialPtr material)
        : center(center), radius(radius), material(material) {}

    Hit intersect(const Ray& ray) const override {
        Vec3 oc = ray.origin - center;

        double a = dot(ray.direction, ray.direction);
        double b = 2.0 * dot(oc, ray.direction);
        double c = dot(oc, oc) - radius * radius;

        double discriminant = b * b - 4.0 * a * c;

        if (discriminant < 0) {
            return NO_HIT;
        }

        double sqrt_disc = std::sqrt(discriminant);
        double t1 = (-b - sqrt_disc) / (2.0 * a);
        double t2 = (-b + sqrt_disc) / (2.0 * a);

        double t = t1;
        if (t < ray.t_min || t > ray.t_max) {
            t = t2;
            if (t < ray.t_min || t > ray.t_max) {
                return NO_HIT;
            }
        }

        Point3 hit_point = ray.at(t);
        Vec3 normal = normalize(hit_point - center);

        return Hit(t, hit_point, normal, material);
    }

    std::vector<Hit> all_hits(const Ray& ray) const override {
        Vec3 oc = ray.origin - center;

        double a = dot(ray.direction, ray.direction);
        double b = 2.0 * dot(oc, ray.direction);
        double c = dot(oc, oc) - radius * radius;

        double discriminant = b * b - 4.0 * a * c;

        if (discriminant <= 0) {
            return {};
        }

        double sqrt_disc = std::sqrt(discriminant);
        double t_near = (-b - sqrt_disc) / (2.0 * a);
        double t_far = (-b + sqrt_disc) / (2.0 * a);

        std::vector<Hit> hits;

        if (t_near >= ray.t_min && t_near <= ray.t_max) {
            Point3 p = ray.at(t_near);
            hits.emplace_back(t_near, p, normalize(p - center), material);
        }

        if (t_far >= ray.t_min && t_far <= ray.t_max) {
            Point3 p = ray.at(t_far);
            hits.emplace_back(t_far, p, normalize(p - center), material);
        }

        return hits;
    }

    AABB bounds() const override {
        return AABB(
            center - Vec3(radius),
            center + Vec3(radius)
        );
    }
};

} // namespace rt
