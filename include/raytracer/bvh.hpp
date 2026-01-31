#pragma once

#include "shape.hpp"
#include <algorithm>
#include <memory>

namespace rt {

class BVHNode : public Shape {
public:
    AABB bbox;
    std::shared_ptr<Shape> left;
    std::shared_ptr<Shape> right;
    std::vector<ShapePtr> primitives;  // Only for leaf nodes

    BVHNode() = default;

    bool is_leaf() const {
        return !primitives.empty();
    }

    Hit intersect(const Ray& ray) const override {
        double t_near, t_far;
        if (!bbox.intersect(ray, t_near, t_far)) {
            return NO_HIT;
        }

        if (is_leaf()) {
            Hit closest = NO_HIT;
            double closest_t = ray.t_max;

            for (const auto& prim : primitives) {
                Hit h = prim->intersect(ray);
                if (h.is_hit() && h.t >= ray.t_min && h.t < closest_t) {
                    closest = h;
                    closest_t = h.t;
                }
            }
            return closest;
        }

        // Check children in order of potential intersection
        struct ChildHit {
            double t_min;
            Shape* child;
        };

        std::vector<ChildHit> children;

        if (left) {
            double t_min_left, t_max_left;
            if (left->bounds().intersect(ray, t_min_left, t_max_left)) {
                children.push_back({t_min_left, left.get()});
            }
        }

        if (right) {
            double t_min_right, t_max_right;
            if (right->bounds().intersect(ray, t_min_right, t_max_right)) {
                children.push_back({t_min_right, right.get()});
            }
        }

        if (children.empty()) {
            return NO_HIT;
        }

        // Sort by t_min to check closer child first
        std::sort(children.begin(), children.end(),
                  [](const ChildHit& a, const ChildHit& b) { return a.t_min < b.t_min; });

        Hit closest = NO_HIT;
        double closest_t = ray.t_max;

        for (const auto& ch : children) {
            if (ch.t_min > closest_t) {
                break;  // Early termination
            }

            Ray limited_ray(ray.origin, ray.direction, ray.t_min, closest_t);
            Hit h = ch.child->intersect(limited_ray);

            if (h.is_hit() && h.t >= ray.t_min && h.t < closest_t) {
                closest = h;
                closest_t = h.t;
            }
        }

        return closest;
    }

    std::vector<Hit> all_hits(const Ray& ray) const override {
        double t_near, t_far;
        if (!bbox.intersect(ray, t_near, t_far)) {
            return {};
        }

        if (is_leaf()) {
            std::vector<Hit> hits;
            for (const auto& prim : primitives) {
                auto prim_hits = prim->all_hits(ray);
                hits.insert(hits.end(), prim_hits.begin(), prim_hits.end());
            }
            std::sort(hits.begin(), hits.end(),
                      [](const Hit& a, const Hit& b) { return a.t < b.t; });
            return hits;
        }

        std::vector<Hit> hits;
        if (left) {
            auto left_hits = left->all_hits(ray);
            hits.insert(hits.end(), left_hits.begin(), left_hits.end());
        }
        if (right) {
            auto right_hits = right->all_hits(ray);
            hits.insert(hits.end(), right_hits.begin(), right_hits.end());
        }

        std::sort(hits.begin(), hits.end(),
                  [](const Hit& a, const Hit& b) { return a.t < b.t; });
        return hits;
    }

    AABB bounds() const override {
        return bbox;
    }
};

// Build BVH from a list of shapes
inline std::shared_ptr<BVHNode> build_bvh(std::vector<ShapePtr>& shapes,
                                           int max_leaf_size = 4,
                                           int depth = 0) {
    if (shapes.empty()) {
        return nullptr;
    }

    auto node = std::make_shared<BVHNode>();

    // Compute bounding box for all shapes
    node->bbox = shapes[0]->bounds();
    for (size_t i = 1; i < shapes.size(); ++i) {
        node->bbox = AABB::merge(node->bbox, shapes[i]->bounds());
    }

    // Create leaf if few shapes or max depth reached
    if (static_cast<int>(shapes.size()) <= max_leaf_size || depth > 32) {
        node->primitives = shapes;
        return node;
    }

    // Find split axis (longest extent)
    int split_axis = node->bbox.longest_axis();

    // Sort shapes by centroid along split axis
    std::sort(shapes.begin(), shapes.end(),
              [split_axis](const ShapePtr& a, const ShapePtr& b) {
                  return a->bounds().centroid()[split_axis] < b->bounds().centroid()[split_axis];
              });

    // Split at midpoint
    size_t mid = shapes.size() / 2;

    std::vector<ShapePtr> left_shapes(shapes.begin(), shapes.begin() + mid);
    std::vector<ShapePtr> right_shapes(shapes.begin() + mid, shapes.end());

    // Handle degenerate case
    if (left_shapes.empty() || right_shapes.empty()) {
        node->primitives = shapes;
        return node;
    }

    node->left = build_bvh(left_shapes, max_leaf_size, depth + 1);
    node->right = build_bvh(right_shapes, max_leaf_size, depth + 1);

    return node;
}

} // namespace rt
