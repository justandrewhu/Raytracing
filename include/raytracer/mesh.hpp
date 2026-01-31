#pragma once

#include "shape.hpp"
#include "triangle.hpp"
#include "bvh.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

namespace rt {

class Mesh : public Shape {
public:
    std::vector<ShapePtr> triangles;
    MaterialPtr material;
    std::shared_ptr<BVHNode> bvh_root;
    AABB bbox;

    Mesh(const std::vector<std::array<Point3, 3>>& faces, MaterialPtr material)
        : material(material)
    {
        for (const auto& face : faces) {
            triangles.push_back(std::make_shared<Triangle>(face, material));
        }
        build();
    }

    Mesh(std::vector<ShapePtr>&& tris, MaterialPtr material)
        : triangles(std::move(tris)), material(material)
    {
        build();
    }

    void build() {
        if (triangles.empty()) {
            bbox = AABB(Vec3(0), Vec3(0));
            return;
        }

        // Build BVH
        std::vector<ShapePtr> shapes_copy = triangles;
        bvh_root = build_bvh(shapes_copy);

        if (bvh_root) {
            bbox = bvh_root->bounds();
        } else {
            bbox = triangles[0]->bounds();
            for (size_t i = 1; i < triangles.size(); ++i) {
                bbox = AABB::merge(bbox, triangles[i]->bounds());
            }
        }
    }

    Hit intersect(const Ray& ray) const override {
        if (bvh_root) {
            Hit h = bvh_root->intersect(ray);
            if (h.is_hit()) {
                // Override material with mesh material
                return Hit(h.t, h.point, h.normal, material);
            }
            return NO_HIT;
        }

        // Fallback: linear search
        Hit closest = NO_HIT;
        double closest_t = ray.t_max;

        for (const auto& tri : triangles) {
            Hit h = tri->intersect(ray);
            if (h.is_hit() && h.t >= ray.t_min && h.t < closest_t) {
                closest = Hit(h.t, h.point, h.normal, material);
                closest_t = h.t;
            }
        }

        return closest;
    }

    std::vector<Hit> all_hits(const Ray& ray) const override {
        std::vector<Hit> hits;

        if (bvh_root) {
            hits = bvh_root->all_hits(ray);
        } else {
            for (const auto& tri : triangles) {
                Hit h = tri->intersect(ray);
                if (h.is_hit() && h.t >= ray.t_min && h.t <= ray.t_max) {
                    hits.push_back(h);
                }
            }
        }

        // Override material
        for (auto& h : hits) {
            h.material = material;
        }

        std::sort(hits.begin(), hits.end(),
                  [](const Hit& a, const Hit& b) { return a.t < b.t; });

        return hits;
    }

    AABB bounds() const override {
        return bbox;
    }

    // Transform all vertices (scale, rotate, translate)
    void transform(const Mat3& rotation, const Vec3& scale, const Vec3& translation) {
        std::vector<std::array<Point3, 3>> new_faces;

        for (const auto& tri_ptr : triangles) {
            auto* tri = dynamic_cast<Triangle*>(tri_ptr.get());
            if (tri) {
                std::array<Point3, 3> new_verts;
                for (int i = 0; i < 3; ++i) {
                    // Apply scale, then rotation, then translation
                    Vec3 v = tri->vertices[i];
                    v = Vec3(v.x * scale.x, v.y * scale.y, v.z * scale.z);
                    v = rotation * v;
                    v = v + translation;
                    new_verts[i] = v;
                }
                new_faces.push_back(new_verts);
            }
        }

        // Rebuild mesh with transformed vertices
        triangles.clear();
        for (const auto& face : new_faces) {
            triangles.push_back(std::make_shared<Triangle>(face, material));
        }
        build();
    }

    // Convenience methods
    void scale(double s) {
        transform(Mat3::identity(), Vec3(s, s, s), Vec3(0));
    }

    void scale(double sx, double sy, double sz) {
        transform(Mat3::identity(), Vec3(sx, sy, sz), Vec3(0));
    }

    void translate(const Vec3& t) {
        transform(Mat3::identity(), Vec3(1, 1, 1), t);
    }

    void rotate(double rx, double ry, double rz) {
        transform(Mat3::rotate(rx, ry, rz), Vec3(1, 1, 1), Vec3(0));
    }

    // Load from OBJ file
    static std::shared_ptr<Mesh> load_obj(const std::string& filename, MaterialPtr material) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open OBJ file: " + filename);
        }

        std::vector<Point3> positions;
        std::vector<std::array<int, 3>> face_indices;

        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string prefix;
            iss >> prefix;

            if (prefix == "v") {
                double x, y, z;
                iss >> x >> y >> z;
                positions.emplace_back(x, y, z);
            } else if (prefix == "f") {
                std::vector<int> indices;
                std::string vertex;
                while (iss >> vertex) {
                    // Handle "v", "v/vt", "v/vt/vn", "v//vn" formats
                    size_t slash_pos = vertex.find('/');
                    int idx;
                    if (slash_pos != std::string::npos) {
                        idx = std::stoi(vertex.substr(0, slash_pos));
                    } else {
                        idx = std::stoi(vertex);
                    }
                    indices.push_back(idx - 1);
                }

                for (size_t i = 1; i + 1 < indices.size(); ++i) {
                    face_indices.push_back({indices[0], indices[i], indices[i + 1]});
                }
            }
        }

        std::vector<std::array<Point3, 3>> faces;
        for (const auto& fi : face_indices) {
            faces.push_back({positions[fi[0]], positions[fi[1]], positions[fi[2]]});
        }

        return std::make_shared<Mesh>(faces, material);
    }
};

// Helper to create a box mesh
inline std::shared_ptr<Mesh> make_box(const Point3& min_corner, const Point3& max_corner, MaterialPtr material) {
    Point3 p0 = min_corner;
    Point3 p7 = max_corner;
    Point3 p1(p7.x, p0.y, p0.z);
    Point3 p2(p7.x, p0.y, p7.z);
    Point3 p3(p0.x, p0.y, p7.z);
    Point3 p4(p0.x, p7.y, p0.z);
    Point3 p5(p7.x, p7.y, p0.z);
    Point3 p6(p7.x, p7.y, p7.z);
    Point3 p_4(p0.x, p7.y, p7.z);  // renamed from p7 collision

    std::vector<std::array<Point3, 3>> faces = {
        // Bottom
        {p0, p1, p2}, {p0, p2, p3},
        // Top
        {p4, p6, p5}, {p4, p_4, p6},
        // Front
        {p0, p1, p5}, {p0, p5, p4},
        // Back
        {p3, p6, p2}, {p3, p_4, p6},
        // Left
        {p0, p3, p_4}, {p0, p_4, p4},
        // Right
        {p1, p2, p6}, {p1, p6, p5}
    };

    return std::make_shared<Mesh>(faces, material);
}

} // namespace rt
