#pragma once

#include "vec3.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include <memory>

namespace rt {

// Forward declaration
class Scene;

class Light {
public:
    virtual ~Light() = default;

    // Calculate the illumination at a hit point
    // Returns the color contribution from this light
    virtual Color illuminate(const Ray& ray, const Hit& hit, const Scene& scene) const = 0;
};

using LightPtr = std::shared_ptr<Light>;

class PointLight : public Light {
public:
    Point3 position;
    Color intensity;

    PointLight(const Point3& position, const Color& intensity)
        : position(position), intensity(intensity) {}

    PointLight(const Point3& position, double intensity)
        : position(position), intensity(Color(intensity)) {}

    Color illuminate(const Ray& ray, const Hit& hit, const Scene& scene) const override;
};

class AmbientLight : public Light {
public:
    Color intensity;

    AmbientLight(const Color& intensity) : intensity(intensity) {}
    AmbientLight(double intensity) : intensity(Color(intensity)) {}

    Color illuminate(const Ray& ray, const Hit& hit, const Scene& scene) const override;
};

} // namespace rt
