#pragma once

#include "vec3.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include <random>
#include <cmath>

namespace rt {

// Base class for area lights that can be importance sampled
class AreaLight {
public:
    Color emission;
    int samples;

    AreaLight(const Color& emission, int samples = 16)
        : emission(emission), samples(samples) {}

    virtual ~AreaLight() = default;

    // Sample a point on the light surface, returns point and normal at that point
    virtual void sample_point(Point3& point, Vec3& normal, std::mt19937& rng) const = 0;

    // Sample a direction from a given point toward the light
    virtual Vec3 sample_direction_from_point(const Point3& point, std::mt19937& rng) const = 0;

    // PDF of sampling a direction from point toward the light (in solid angle measure)
    virtual double pdf_direction(const Point3& point, const Vec3& direction) const = 0;

    // Check if ray hits this light, returns t and hit point
    virtual bool intersect(const Point3& origin, const Vec3& direction, double& t, Point3& hit_point) const = 0;

    // Get the surface area
    virtual double area() const = 0;
};

using AreaLightPtr = std::shared_ptr<AreaLight>;

class SphereAreaLight : public AreaLight {
public:
    Point3 center;
    double radius;

    SphereAreaLight(const Point3& center, double radius, const Color& emission, int samples = 16)
        : AreaLight(emission, samples), center(center), radius(radius) {}

    void sample_point(Point3& point, Vec3& normal, std::mt19937& rng) const override {
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        double z = 2.0 * dist(rng) - 1.0;
        double r = std::sqrt(std::max(0.0, 1.0 - z * z));
        double phi = 2.0 * PI * dist(rng);

        normal = Vec3(r * std::cos(phi), r * std::sin(phi), z);
        point = center + radius * normal;
    }

    Vec3 sample_direction_from_point(const Point3& point, std::mt19937& rng) const override {
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        Vec3 to_center = center - point;
        double dist_to_center = to_center.length();

        if (dist_to_center < radius) {
            double z = 2.0 * dist(rng) - 1.0;
            double r = std::sqrt(std::max(0.0, 1.0 - z * z));
            double phi = 2.0 * PI * dist(rng);
            return normalize(Vec3(r * std::cos(phi), r * std::sin(phi), z));
        }

        Vec3 direction_to_center = to_center / dist_to_center;

        double sin_theta_max = radius / dist_to_center;
        double cos_theta_max = std::sqrt(std::max(0.0, 1.0 - sin_theta_max * sin_theta_max));

        double cos_theta = 1.0 - dist(rng) * (1.0 - cos_theta_max);
        double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
        double phi = 2.0 * PI * dist(rng);

        Vec3 arbitrary = (std::fabs(direction_to_center.z) > 0.9) ? Vec3(1, 0, 0) : Vec3(0, 0, 1);
        Vec3 tangent = normalize(cross(arbitrary, direction_to_center));
        Vec3 bitangent = cross(direction_to_center, tangent);

        return normalize(
            cos_theta * direction_to_center +
            sin_theta * std::cos(phi) * tangent +
            sin_theta * std::sin(phi) * bitangent
        );
    }

    double pdf_direction(const Point3& point, const Vec3& direction) const override {
        Vec3 to_center = center - point;
        double dist_to_center = to_center.length();

        if (dist_to_center < radius) {
            return 1.0 / (4.0 * PI);
        }

        double sin_theta_max = radius / dist_to_center;
        double cos_theta_max = std::sqrt(std::max(0.0, 1.0 - sin_theta_max * sin_theta_max));
        double solid_angle = 2.0 * PI * (1.0 - cos_theta_max);

        return solid_angle > 0 ? 1.0 / solid_angle : 0.0;
    }

    bool intersect(const Point3& origin, const Vec3& direction, double& t, Point3& hit_point) const override {
        Vec3 oc = origin - center;
        double a = dot(direction, direction);
        double b = 2.0 * dot(oc, direction);
        double c = dot(oc, oc) - radius * radius;
        double disc = b * b - 4.0 * a * c;

        if (disc < 0) return false;

        double sqrt_disc = std::sqrt(disc);
        t = (-b - sqrt_disc) / (2.0 * a);
        if (t < EPSILON) {
            t = (-b + sqrt_disc) / (2.0 * a);
            if (t < EPSILON) return false;
        }

        hit_point = origin + t * direction;
        return true;
    }

    double area() const override {
        return 4.0 * PI * radius * radius;
    }
};

class RectAreaLight : public AreaLight {
public:
    Point3 corner;
    Vec3 edge1, edge2;
    Vec3 normal_;
    double area_;

    RectAreaLight(const Point3& corner, const Vec3& edge1, const Vec3& edge2,
                  const Color& emission, int samples = 16)
        : AreaLight(emission, samples), corner(corner), edge1(edge1), edge2(edge2) {
        Vec3 c = cross(edge1, edge2);
        area_ = c.length();
        normal_ = normalize(c);
    }

    void sample_point(Point3& point, Vec3& normal, std::mt19937& rng) const override {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(rng);
        double v = dist(rng);
        point = corner + u * edge1 + v * edge2;
        normal = normal_;
    }

    Vec3 sample_direction_from_point(const Point3& point, std::mt19937& rng) const override {
        Point3 light_point;
        Vec3 light_normal;
        sample_point(light_point, light_normal, rng);
        return normalize(light_point - point);
    }

    double pdf_direction(const Point3& point, const Vec3& direction) const override {
        double t;
        Point3 hit_point;
        if (!intersect(point, direction, t, hit_point)) {
            return 0.0;
        }

        double dist_sq = (hit_point - point).length_squared();
        double cos_light = std::fabs(dot(normal_, -direction));

        if (cos_light < EPSILON) return 0.0;

        return dist_sq / (area_ * cos_light);
    }

    bool intersect(const Point3& origin, const Vec3& direction, double& t, Point3& hit_point) const override {
        double denom = dot(direction, normal_);
        if (std::fabs(denom) < EPSILON) return false;

        t = dot(corner - origin, normal_) / denom;
        if (t < EPSILON) return false;

        hit_point = origin + t * direction;
        Vec3 v = hit_point - corner;

        // Project onto edges
        double u1 = dot(v, edge1) / dot(edge1, edge1);
        double u2 = dot(v, edge2) / dot(edge2, edge2);

        if (u1 >= 0 && u1 <= 1 && u2 >= 0 && u2 <= 1) {
            return true;
        }
        return false;
    }

    double area() const override {
        return area_;
    }
};

class BoxAreaLight : public AreaLight {
public:
    Point3 min_corner, max_corner;
    std::vector<RectAreaLight> faces;
    double total_area_;

    BoxAreaLight(const Point3& min_corner, const Point3& max_corner,
                 const Color& emission, int samples = 16)
        : AreaLight(emission, samples), min_corner(min_corner), max_corner(max_corner) {

        Vec3 size = max_corner - min_corner;


        faces.emplace_back(min_corner, Vec3(size.x, 0, 0), Vec3(0, 0, size.z), emission, samples);

        faces.emplace_back(Point3(min_corner.x, max_corner.y, min_corner.z),
                          Vec3(size.x, 0, 0), Vec3(0, 0, size.z), emission, samples);

        faces.emplace_back(min_corner, Vec3(size.x, 0, 0), Vec3(0, size.y, 0), emission, samples);

        faces.emplace_back(Point3(min_corner.x, min_corner.y, max_corner.z),
                          Vec3(size.x, 0, 0), Vec3(0, size.y, 0), emission, samples);

        faces.emplace_back(min_corner, Vec3(0, size.y, 0), Vec3(0, 0, size.z), emission, samples);

        faces.emplace_back(Point3(max_corner.x, min_corner.y, min_corner.z),
                          Vec3(0, size.y, 0), Vec3(0, 0, size.z), emission, samples);

        total_area_ = 0;
        for (const auto& face : faces) {
            total_area_ += face.area();
        }
    }

    void sample_point(Point3& point, Vec3& normal, std::mt19937& rng) const override {
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        double r = dist(rng) * total_area_;
        double cumulative = 0;

        for (const auto& face : faces) {
            cumulative += face.area();
            if (r <= cumulative) {
                face.sample_point(point, normal, rng);
                return;
            }
        }

        // Fallback to last face
        faces.back().sample_point(point, normal, rng);
    }

    Vec3 sample_direction_from_point(const Point3& point, std::mt19937& rng) const override {
        Point3 light_point;
        Vec3 light_normal;
        sample_point(light_point, light_normal, rng);
        return normalize(light_point - point);
    }

    double pdf_direction(const Point3& point, const Vec3& direction) const override {
        // Find which face we hit
        double closest_t = INF;
        Point3 closest_hit;
        const RectAreaLight* hit_face = nullptr;

        for (const auto& face : faces) {
            double t;
            Point3 hit_point;
            if (face.intersect(point, direction, t, hit_point) && t < closest_t) {
                closest_t = t;
                closest_hit = hit_point;
                hit_face = &face;
            }
        }

        if (!hit_face) return 0.0;

        double dist_sq = (closest_hit - point).length_squared();
        double cos_light = std::fabs(dot(hit_face->normal_, -direction));

        if (cos_light < EPSILON) return 0.0;

        // PDF = 1/total_area in area measure, convert to solid angle
        return dist_sq / (total_area_ * cos_light);
    }

    bool intersect(const Point3& origin, const Vec3& direction, double& t, Point3& hit_point) const override {
        double closest_t = INF;
        Point3 closest_hit;
        bool found = false;

        for (const auto& face : faces) {
            double face_t;
            Point3 face_hit;
            if (face.intersect(origin, direction, face_t, face_hit) && face_t < closest_t) {
                closest_t = face_t;
                closest_hit = face_hit;
                found = true;
            }
        }

        if (found) {
            t = closest_t;
            hit_point = closest_hit;
        }
        return found;
    }

    double area() const override {
        return total_area_;
    }
};

inline Vec3 sample_cosine_hemisphere(const Vec3& normal, std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double u1 = dist(rng);
    double u2 = dist(rng);

    double r = std::sqrt(u1);
    double theta = 2.0 * PI * u2;

    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    double z = std::sqrt(std::max(0.0, 1.0 - u1));

    Vec3 arbitrary = (std::fabs(normal.z) > 0.9) ? Vec3(1, 0, 0) : Vec3(0, 0, 1);
    Vec3 tangent = normalize(cross(arbitrary, normal));
    Vec3 bitangent = cross(normal, tangent);

    return normalize(x * tangent + y * bitangent + z * normal);
}

inline double cosine_hemisphere_pdf(double cos_theta) {
    return cos_theta > 0 ? cos_theta / PI : 0.0;
}

} // namespace rt
