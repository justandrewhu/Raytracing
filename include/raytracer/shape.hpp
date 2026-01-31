#pragma once

#include "ray.hpp"
#include "hit.hpp"
#include "aabb.hpp"
#include <vector>
#include <memory>

namespace rt {

// Forward declaration
class Shape;
using ShapePtr = std::shared_ptr<Shape>;

class Shape {
public:
    virtual ~Shape() = default;

    // Returns the first intersection with the ray
    virtual Hit intersect(const Ray& ray) const = 0;

    // Returns all intersections, used for CSG
    virtual std::vector<Hit> all_hits(const Ray& ray) const {
        Hit h = intersect(ray);
        if (h.is_hit()) {
            return {h};
        }
        return {};
    }

    // Returns the bounding box
    virtual AABB bounds() const = 0;
};

} // namespace rt
