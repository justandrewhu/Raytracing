#pragma once

#include "shape.hpp"
#include "light.hpp"
#include "bvh.hpp"
#include "camera.hpp"
#include <vector>

namespace rt {

class Scene {
public:
    std::vector<ShapePtr> shapes;
    std::vector<LightPtr> lights;
    Color background_color;
    std::shared_ptr<BVHNode> bvh_root;

    Scene(const Color& bg_color = Color(0.2, 0.3, 0.5))
        : background_color(bg_color) {}

    void add_shape(ShapePtr shape) {
        shapes.push_back(std::move(shape));
    }

    void add_light(LightPtr light) {
        lights.push_back(std::move(light));
    }

    void build_bvh() {
        if (!shapes.empty()) {
            std::vector<ShapePtr> shapes_copy = shapes;
            bvh_root = rt::build_bvh(shapes_copy);
        }
    }

    Hit intersect(const Ray& ray) const {
        if (bvh_root) {
            return bvh_root->intersect(ray);
        }

        // Linear fallback
        Hit closest = NO_HIT;
        double closest_t = ray.t_max;

        for (const auto& shape : shapes) {
            Hit h = shape->intersect(ray);
            if (h.is_hit() && h.t >= ray.t_min && h.t < closest_t) {
                closest = h;
                closest_t = h.t;
            }
        }

        return closest;
    }
};

// Now implement the light illuminate methods (need Scene defined)

inline Color PointLight::illuminate(const Ray& ray, const Hit& hit, const Scene& scene) const {
    Point3 point = hit.point;
    Vec3 normal = hit.normal;
    Vec3 to_light = position - point;
    double dist_squared = dot(to_light, to_light);
    double dist = std::sqrt(dist_squared);
    Vec3 L = to_light / dist;

    // Shadow ray
    Ray shadow_ray(point + normal * EPSILON * 10, to_light, 0.0, 1.0);
    Hit shadow_hit = scene.intersect(shadow_ray);

    Color attenuation(1.0);

    if (shadow_hit.is_hit() && shadow_hit.t < 1.0) {
        if (shadow_hit.material && shadow_hit.material->has_transmission()) {
            attenuation = shadow_hit.material->k_t;
            if (attenuation.x < EPSILON && attenuation.y < EPSILON && attenuation.z < EPSILON) {
                return Color(0);
            }
        } else {
            return Color(0);  // Fully shadowed
        }
    }

    // Irradiance
    Color irradiance = intensity / dist_squared;

    // View direction
    Vec3 view = normalize(-ray.direction);

    // Material properties
    Color k_d = hit.material->k_d;
    Color k_s = hit.material->k_s;
    double p = hit.material->p;

    // Diffuse (Lambertian)
    double cos_theta = std::max(0.0, dot(normal, L));
    Color diffuse = k_d * cos_theta * irradiance;

    // Specular (Blinn-Phong)
    Vec3 H = normalize(L + view);
    double cos_alpha = std::max(0.0, dot(normal, H));
    Color specular = k_s * std::pow(cos_alpha, p) * cos_theta * irradiance;

    return attenuation * (diffuse + specular);
}

inline Color AmbientLight::illuminate(const Ray& ray, const Hit& hit, const Scene& scene) const {
    (void)ray;
    (void)scene;
    Color k_a = hit.material->k_a;
    return k_a * intensity;
}

} // namespace rt
