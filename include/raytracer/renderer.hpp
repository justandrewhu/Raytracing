#pragma once

#include "scene.hpp"
#include "camera.hpp"
#include "area_light.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace rt {

constexpr int MAX_DEPTH = 4;
constexpr int PATH_TRACE_MAX_DEPTH = 8;


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

// blinn phong brdf
inline Color eval_brdf(const Material& mat, const Vec3& normal, const Vec3& view, const Vec3& light_dir) {
    Color f_diffuse = mat.k_d / PI;

    Vec3 H = normalize(light_dir + view);
    double cos_h = std::max(0.0, dot(normal, H));
    Color f_specular = mat.k_s * ((mat.p + 2.0) / (2.0 * PI)) * std::pow(cos_h, mat.p);

    return f_diffuse + f_specular;
}

// path tracing w/ multiple importance sampling
class PathTracer {
public:
    int samples_per_pixel = 64;
    int max_depth = PATH_TRACE_MAX_DEPTH;
    double light_sample_prob = 0.5;  // Probability of sampling light vs BRDF
    std::vector<AreaLightPtr> area_lights;

    void add_area_light(AreaLightPtr light) {
        area_lights.push_back(std::move(light));
    }

    // Check if a ray hits any area light, returns the light and emission
    bool hit_area_light(const Point3& origin, const Vec3& direction,
                        double max_t, Color& emission, double& t) const {
        double closest_t = max_t;
        bool found = false;

        for (const auto& light : area_lights) {
            double light_t;
            Point3 hit_point;
            if (light->intersect(origin, direction, light_t, hit_point) && light_t < closest_t) {
                closest_t = light_t;
                emission = light->emission;
                found = true;
            }
        }

        if (found) {
            t = closest_t;
        }
        return found;
    }

    // MIS weight using balance heuristic
    double mis_weight(double pdf1, double pdf2) const {
        return pdf1 / (pdf1 + pdf2 + 1e-10);
    }

    // Path trace a single ray with MIS
    Color path_trace(Ray ray, const Scene& scene, std::mt19937& rng) const {
        Color color(0);
        Color throughput(1);
        bool inside = false;

        for (int depth = 0; depth < max_depth; ++depth) {
            Hit hit = scene.intersect(ray);

            // Check if we hit an area light directly
            Color light_emission;
            double light_t;
            bool hit_light = hit_area_light(ray.origin, ray.direction,
                                            hit.is_hit() ? hit.t : INF,
                                            light_emission, light_t);

            if (hit_light && (!hit.is_hit() || light_t < hit.t)) {
                if (depth == 0) {
                    color += throughput * light_emission;
                }
                break;
            }

            if (!hit.is_hit()) {
                color += throughput * scene.background_color;
                break;
            }

            if (!hit.material) break;

            // Add emission from emissive materials
            if (hit.material->has_emission()) {
                color += throughput * hit.material->k_e;
                if (depth > 0) break;  // Don't continue after hitting emissive surface
            }

            Point3 point = hit.point;
            Vec3 normal = normalize(hit.normal);
            Vec3 view = normalize(-ray.direction);

            if (dot(normal, view) < 0) {
                normal = -normal;
            }

            const Material& mat = *hit.material;

            // Handle transmission
            if (mat.has_transmission()) {
                double ior1 = inside ? mat.ior : 1.0;
                double ior2 = inside ? 1.0 : mat.ior;
                double fresnel = fresnel_schlick(normal, -view, ior1, ior2);

                std::uniform_real_distribution<double> dist(0.0, 1.0);

                if (dist(rng) < fresnel) {
                    // Reflect
                    Vec3 reflect_dir = reflect(-view, normal);
                    ray = Ray(point + normal * EPSILON * 10, reflect_dir);
                    throughput *= mat.k_t;
                } else {
                    // Refract
                    Vec3 refract_normal = inside ? -normal : normal;
                    Vec3 refract_dir;
                    if (refract(-view, refract_normal, ior1 / ior2, refract_dir)) {
                        ray = Ray(point - refract_normal * EPSILON * 10, refract_dir);
                        throughput *= mat.k_t;
                        inside = !inside;
                    } else {
                        // Total internal reflection
                        Vec3 reflect_dir = reflect(-view, normal);
                        ray = Ray(point + normal * EPSILON * 10, reflect_dir);
                        throughput *= mat.k_t;
                    }
                }
                continue;
            }

            if (depth > 2) {
                double p_continue = std::min(0.95, std::max({throughput.x, throughput.y, throughput.z}));
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                if (dist(rng) > p_continue) {
                    break;
                }
                throughput = throughput / p_continue;
            }

            std::uniform_real_distribution<double> dist(0.0, 1.0);

            for (const auto& area_light : area_lights) {
                int n_light_samples = std::max(1, area_light->samples / 2);
                Color light_contrib(0);

                for (int s = 0; s < n_light_samples; ++s) {
                    Vec3 L = area_light->sample_direction_from_point(point, rng);
                    double cos_o = dot(normal, L);

                    if (cos_o <= 0) continue;

                    Point3 shadow_origin = point + normal * EPSILON * 10;

                    double light_t;
                    Point3 light_hit;
                    if (!area_light->intersect(shadow_origin, L, light_t, light_hit)) {
                        continue;
                    }

                    Ray shadow_ray(shadow_origin, L, 0, light_t - EPSILON);
                    Hit shadow_hit = scene.intersect(shadow_ray);

                    if (shadow_hit.is_hit() && shadow_hit.t < light_t - EPSILON * 10) {
                        if (!shadow_hit.material || !shadow_hit.material->has_transmission()) {
                            continue;  
                        }
                    }

                    Color f_r = eval_brdf(mat, normal, view, L);

                    double pdf_light = area_light->pdf_direction(point, L);
                    double pdf_brdf = cosine_hemisphere_pdf(cos_o);

                    double weight = mis_weight(pdf_light * light_sample_prob,
                                               pdf_brdf * (1.0 - light_sample_prob));

                    Color contrib = weight * f_r * cos_o * area_light->emission / (pdf_light + 1e-10);
                    light_contrib += contrib;
                }

                color += throughput * light_contrib / static_cast<double>(n_light_samples);

                // BRDF SAMPLING for this light
                int n_brdf_samples = std::max(1, area_light->samples - n_light_samples);
                Color brdf_contrib(0);

                for (int s = 0; s < n_brdf_samples; ++s) {
                    Vec3 L = sample_cosine_hemisphere(normal, rng);
                    double cos_o = dot(normal, L);

                    if (cos_o <= 0) continue;

                    // Check if we hit this area light
                    Point3 ray_origin = point + normal * EPSILON * 10;
                    double light_t;
                    Point3 light_hit;

                    if (!area_light->intersect(ray_origin, L, light_t, light_hit)) {
                        continue;
                    }

                    // Check for occlusion
                    Ray shadow_ray(ray_origin, L, 0, light_t - EPSILON);
                    Hit shadow_hit = scene.intersect(shadow_ray);

                    if (shadow_hit.is_hit() && shadow_hit.t < light_t - EPSILON * 10) {
                        if (!shadow_hit.material || !shadow_hit.material->has_transmission()) {
                            continue;
                        }
                    }

                    Color f_r = eval_brdf(mat, normal, view, L);

                    double pdf_brdf = cosine_hemisphere_pdf(cos_o);
                    double pdf_light = area_light->pdf_direction(point, L);

                    double weight = mis_weight(pdf_brdf * (1.0 - light_sample_prob),
                                               pdf_light * light_sample_prob);

                    Color contrib = weight * f_r * cos_o * area_light->emission / (pdf_brdf + 1e-10);
                    brdf_contrib += contrib;
                }

                color += throughput * brdf_contrib / static_cast<double>(n_brdf_samples);
            }

            // Continue path with BRDF sampling
            Vec3 new_dir = sample_cosine_hemisphere(normal, rng);
            double cos_o = dot(normal, new_dir);

            if (cos_o <= 0) break;

            Color f_r = eval_brdf(mat, normal, view, new_dir);
            double pdf = cosine_hemisphere_pdf(cos_o);

            throughput *= f_r * cos_o / (pdf + 1e-10);

            ray = Ray(point + normal * EPSILON * 10, new_dir);
        }

        return clamp_color(color);
    }

    std::vector<Color> render(const Camera& camera, const Scene& scene, int width, int height) const {
        std::vector<Color> framebuffer(width * height);

        #pragma omp parallel for schedule(dynamic, 1)
        for (int j = 0; j < height; ++j) {
            std::mt19937 rng(j * width + 42);  // Per-row RNG for reproducibility
            std::uniform_real_distribution<double> dist(0.0, 1.0);

            for (int i = 0; i < width; ++i) {
                Color pixel_color(0);

                for (int s = 0; s < samples_per_pixel; ++s) {
                    double x = (i + dist(rng)) / width;
                    double y = (j + dist(rng)) / height;

                    Ray ray = camera.generate_ray(x, y);
                    pixel_color += path_trace(ray, scene, rng);
                }

                framebuffer[j * width + i] = pixel_color / static_cast<double>(samples_per_pixel);
            }

            if (j % 10 == 0) {
                #pragma omp critical
                std::cout << "Path tracing row " << j << "/" << height << std::endl;
            }
        }

        return framebuffer;
    }
};

} // namespace rt

