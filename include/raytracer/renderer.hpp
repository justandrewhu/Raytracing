#pragma once

#include "scene.hpp"
#include "camera.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace rt {

constexpr int MAX_DEPTH = 4;

// Fresnel-Schlick approximation for fresnel-reflection
inline double fresnel_schlick(const Vec3& normal, const Vec3& view, double ior1, double ior2) {
    double r0 = (ior1 - ior2) / (ior1 + ior2);
    r0 = r0 * r0;
    double cos_theta = std::fabs(dot(normal, -view));
    return r0 + (1.0 - r0) * std::pow(1.0 - cos_theta, 5.0);
}

/** computes the refraction based on Snell's law
 * 
 * @return true if the ray will refract
*/
inline bool refract(const Vec3& incident, const Vec3& normal, double ior_ratio, Vec3& refracted) {
    double cos_i = -dot(incident, normal);
    double sin_t2 = ior_ratio * ior_ratio * (1.0 - cos_i * cos_i);

    if (sin_t2 > 1.0) {
        return false;  // Total internal reflection
    }

    double cos_t = std::sqrt(1.0 - sin_t2);
    refracted = normalize(ior_ratio * incident + (ior_ratio * cos_i - cos_t) * normal);
    return true;
}

class Renderer {
public:
    int samples_per_pixel = 1;

    /** Recursive ray-tracing shading function */
    Color shade(const Ray& ray, const Hit& hit, const Scene& scene, int depth, bool inside) const {
        if (!hit.is_hit()) {
            return scene.background_color;
        }

        if (hit.material && hit.material->has_emission() && depth > 0) {
            return hit.material->k_e;
        }

        Color color(0);

        for (const auto& light : scene.lights) {
            color += light->illuminate(ray, hit, scene);
        }

        if (depth < MAX_DEPTH && hit.material) {
            Vec3 normal = normalize(hit.normal);
            Vec3 view_dir = normalize(ray.direction);

            double ior1 = inside ? hit.material->ior : 1.0;
            double ior2 = inside ? 1.0 : hit.material->ior;
            double fresnel = fresnel_schlick(normal, view_dir, ior1, ior2);

            // Reflection
            Color reflection_coef = hit.material->k_m;
            Color total_reflection = reflection_coef + fresnel * (Color(1.0) - reflection_coef);

            if (total_reflection.x > EPSILON || total_reflection.y > EPSILON || total_reflection.z > EPSILON) {
                Vec3 reflect_dir = reflect(view_dir, normal);
                Point3 offset_origin = hit.point + normal * EPSILON * 10;
                Ray reflect_ray(offset_origin, reflect_dir);
                Hit reflect_hit = scene.intersect(reflect_ray);

                Color reflected_color = shade(reflect_ray, reflect_hit, scene, depth + 1, inside);
                color += total_reflection * reflected_color;
            }

            Color transmission_coef = hit.material->k_t;
            if (transmission_coef.x > EPSILON || transmission_coef.y > EPSILON || transmission_coef.z > EPSILON) {
                Vec3 refract_normal = inside ? -normal : normal;
                Vec3 refract_dir;

                if (refract(view_dir, refract_normal, ior1 / ior2, refract_dir)) {
                    Vec3 offset_normal = -refract_normal;
                    Point3 offset_origin = hit.point + offset_normal * EPSILON * 10;
                    Ray refract_ray(offset_origin, refract_dir);
                    Hit refract_hit = scene.intersect(refract_ray);

                    Color refracted_color = shade(refract_ray, refract_hit, scene, depth + 1, !inside);
                    color += transmission_coef * (1.0 - fresnel) * refracted_color;
                } else {
                    Vec3 reflect_dir = reflect(view_dir, normal);
                    Point3 offset_origin = hit.point + reflect_dir * EPSILON * 10;
                    Ray reflect_ray(offset_origin, reflect_dir);
                    Hit reflect_hit = scene.intersect(reflect_ray);

                    Color tir_color = shade(reflect_ray, reflect_hit, scene, depth + 1, inside);
                    color += transmission_coef * tir_color;
                }
            }
        }

        return clamp_color(color);
    }

    std::vector<Color> render(const Camera& camera, const Scene& scene, int width, int height) const {
        std::vector<Color> framebuffer(width * height);

        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (int j = 0; j < height; ++j) {
            for (int i = 0; i < width; ++i) {
                Color pixel_color(0);

                if (samples_per_pixel > 1) {
                    for (int s = 0; s < samples_per_pixel; ++s) {
                        double x = (i + dist(rng)) / width;
                        double y = (j + dist(rng)) / height;

                        Ray ray = camera.generate_ray(x, y);
                        Hit hit = scene.intersect(ray);
                        pixel_color += shade(ray, hit, scene, 0, false);
                    }
                    pixel_color = pixel_color / static_cast<double>(samples_per_pixel);
                } else {
                    double x = (i + 0.5) / width;
                    double y = (j + 0.5) / height;

                    Ray ray = camera.generate_ray(x, y);
                    Hit hit = scene.intersect(ray);
                    pixel_color = shade(ray, hit, scene, 0, false);
                }

                framebuffer[j * width + i] = pixel_color;
            }

            if (j % 10 == 0) {
                std::cout << "Rendering row " << j << "/" << height << std::endl;
            }
        }

        return framebuffer;
    }
};

} // namespace rt
