#pragma once

#include "shape.hpp"
#include <algorithm>

namespace rt {

enum class CSGOperation {
    Union,
    Intersection,
    Difference
};

class CSGNode : public Shape {
public:
    ShapePtr left;
    ShapePtr right;
    CSGOperation op;
    MaterialPtr material;  // Optional override material
    AABB bbox;

    CSGNode(ShapePtr left, ShapePtr right, CSGOperation op, MaterialPtr material = nullptr)
        : left(std::move(left)), right(std::move(right)), op(op), material(material)
    {
        bbox = AABB::merge(this->left->bounds(), this->right->bounds());
    }

    bool inside_helper(bool inside_left, bool inside_right) const {
        switch (op) {
            case CSGOperation::Union:
                return inside_left || inside_right;
            case CSGOperation::Intersection:
                return inside_left && inside_right;
            case CSGOperation::Difference:
                return inside_left && !inside_right;
        }
        return false;
    }

    std::vector<Hit> all_hits(const Ray& ray) const override {
        // Collect hits from both children
        auto hits_left = left->all_hits(ray);
        auto hits_right = right->all_hits(ray);

        // Create boundary events
        struct Event {
            double t;
            bool is_left;
            Hit hit;
        };

        std::vector<Event> events;
        events.reserve(hits_left.size() + hits_right.size());

        for (const auto& h : hits_left) {
            events.push_back({h.t, true, h});
        }
        for (const auto& h : hits_right) {
            events.push_back({h.t, false, h});
        }

        if (events.empty()) {
            return {};
        }

        std::sort(events.begin(), events.end(),
                  [](const Event& a, const Event& b) { return a.t < b.t; });

        bool inside_left = false;
        bool inside_right = false;
        bool was_inside = inside_helper(inside_left, inside_right);

        std::vector<Hit> result;
        Vec3 ray_dir = ray.direction;

        for (const auto& event : events) {
            if (event.t < ray.t_min || event.t > ray.t_max) {
                continue;
            }

            // Toggle inside state
            if (event.is_left) {
                inside_left = !inside_left;
            } else {
                inside_right = !inside_right;
            }

            bool is_inside = inside_helper(inside_left, inside_right);

            // Boundary crossing
            if (was_inside != is_inside) {
                Vec3 normal = event.hit.normal;

                // Ensure normal points correctly
                if (is_inside) {
                    // Entering: normal should point against ray
                    if (dot(normal, ray_dir) > 0) {
                        normal = -normal;
                    }
                } else {
                    // Leaving: normal should point with ray (outward)
                    if (dot(normal, ray_dir) < 0) {
                        normal = -normal;
                    }
                }

                MaterialPtr hit_mat = material ? material : event.hit.material;
                result.emplace_back(event.t, event.hit.point, normal, hit_mat);
            }

            was_inside = is_inside;
        }

        std::sort(result.begin(), result.end(),
                  [](const Hit& a, const Hit& b) { return a.t < b.t; });

        return result;
    }

    Hit intersect(const Ray& ray) const override {
        auto hits = all_hits(ray);
        for (const auto& h : hits) {
            if (h.t >= ray.t_min && h.t <= ray.t_max) {
                return h;
            }
        }
        return NO_HIT;
    }

    AABB bounds() const override {
        return bbox;
    }
};

// Convenience classes
class CSGUnion : public CSGNode {
public:
    CSGUnion(ShapePtr left, ShapePtr right, MaterialPtr material = nullptr)
        : CSGNode(std::move(left), std::move(right), CSGOperation::Union, material) {}
};

class CSGIntersection : public CSGNode {
public:
    CSGIntersection(ShapePtr left, ShapePtr right, MaterialPtr material = nullptr)
        : CSGNode(std::move(left), std::move(right), CSGOperation::Intersection, material) {}
};

class CSGDifference : public CSGNode {
public:
    CSGDifference(ShapePtr left, ShapePtr right, MaterialPtr material = nullptr)
        : CSGNode(std::move(left), std::move(right), CSGOperation::Difference, material) {}
};

// Factory functions
inline ShapePtr csg_union(ShapePtr left, ShapePtr right, MaterialPtr material = nullptr) {
    return std::make_shared<CSGUnion>(std::move(left), std::move(right), material);
}

inline ShapePtr csg_intersection(ShapePtr left, ShapePtr right, MaterialPtr material = nullptr) {
    return std::make_shared<CSGIntersection>(std::move(left), std::move(right), material);
}

inline ShapePtr csg_difference(ShapePtr left, ShapePtr right, MaterialPtr material = nullptr) {
    return std::make_shared<CSGDifference>(std::move(left), std::move(right), material);
}

} // namespace rt
