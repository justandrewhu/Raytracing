#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include "raytracer/vec3.hpp"
#include "raytracer/ray.hpp"
#include "raytracer/material.hpp"
#include "raytracer/hit.hpp"
#include "raytracer/camera.hpp"
#include "raytracer/aabb.hpp"
#include "raytracer/shape.hpp"
#include "raytracer/sphere.hpp"
#include "raytracer/triangle.hpp"
#include "raytracer/mesh.hpp"
#include "raytracer/bvh.hpp"
#include "raytracer/csg.hpp"
#include "raytracer/light.hpp"
#include "raytracer/scene.hpp"
#include "raytracer/renderer.hpp"
#include "raytracer/scene_loader.hpp"
#include "raytracer/image_output.hpp"

using namespace rt;

SceneLoader::SceneData create_default_scene() {
    SceneLoader::SceneData data;

    // Materials
    auto tan = std::make_shared<Material>(Color(0.7, 0.7, 0.4), 0.6);
    auto gray = std::make_shared<Material>(Color(0.2, 0.2, 0.2));

    // Shapes
    data.scene.add_shape(std::make_shared<Sphere>(Point3(0, 0, 0), 0.5, tan));
    data.scene.add_shape(std::make_shared<Sphere>(Point3(0, -40, 0), 39.5, gray));

    // Lights
    data.scene.add_light(std::make_shared<PointLight>(Point3(12, 10, 5), Color(300, 300, 300)));
    data.scene.add_light(std::make_shared<AmbientLight>(0.1));

    // Camera
    data.camera = Camera(Point3(3, 1.7, 5), Point3(0, 0, 0), Vec3(0, 1, 0), 25, 16.0 / 9.0);

    // Build BVH
    data.scene.build_bvh();

    return data;
}

SceneLoader::SceneData create_three_spheres_scene() {
    SceneLoader::SceneData data;

    auto tan = std::make_shared<Material>(Color(0.4, 0.4, 0.2), 0.3, 90.0, 0.3);
    auto blue = std::make_shared<Material>(Color(0.2, 0.2, 0.5));
    blue->set_mirror(0.5);
    auto gray = std::make_shared<Material>(Color(0.2, 0.2, 0.2));
    gray->set_mirror(0.4);

    data.scene.add_shape(std::make_shared<Sphere>(Point3(-0.7, 0, 0), 0.5, tan));
    data.scene.add_shape(std::make_shared<Sphere>(Point3(0.7, 0, 0), 0.5, blue));
    data.scene.add_shape(std::make_shared<Sphere>(Point3(0, -40, 0), 39.5, gray));

    data.scene.add_light(std::make_shared<PointLight>(Point3(12, 10, 5), Color(300, 300, 300)));
    data.scene.add_light(std::make_shared<AmbientLight>(0.1));

    data.camera = Camera(Point3(3, 1.2, 5), Point3(0, -0.4, 0), Vec3(0, 1, 0), 24, 16.0 / 9.0);

    data.scene.build_bvh();

    return data;
}

SceneLoader::SceneData create_glass_scene() {
    SceneLoader::SceneData data;

    auto glass = std::make_shared<Material>(Color(0.0, 0.0, 0.0));
    glass->set_transmission(Color(0.9, 0.9, 0.9));
    glass->set_ior(1.5);
    glass->set_mirror(0.1);

    auto gray = std::make_shared<Material>(Color(0.3, 0.3, 0.3));
    gray->set_mirror(0.2);

    auto red = std::make_shared<Material>(Color(0.7, 0.2, 0.2), 0.3);

    data.scene.add_shape(std::make_shared<Sphere>(Point3(0, 0.5, 0), 0.5, glass));
    data.scene.add_shape(std::make_shared<Sphere>(Point3(-1.2, 0.3, 0.5), 0.3, red));
    data.scene.add_shape(std::make_shared<Sphere>(Point3(0, -40, 0), 39.5, gray));

    data.scene.add_light(std::make_shared<PointLight>(Point3(5, 8, 5), Color(200, 200, 200)));
    data.scene.add_light(std::make_shared<PointLight>(Point3(-5, 5, 3), Color(100, 100, 150)));
    data.scene.add_light(std::make_shared<AmbientLight>(0.1));

    data.camera = Camera(Point3(3, 2, 4), Point3(0, 0, 0), Vec3(0, 1, 0), 30, 16.0 / 9.0);

    data.scene.build_bvh();
    return data;
}

SceneLoader::SceneData create_csg_scene() {
    SceneLoader::SceneData data;

    auto white = std::make_shared<Material>(Color(0.9, 0.9, 0.9), 0.3);
    auto blue = std::make_shared<Material>(Color(0.3, 0.3, 0.8), 0.5);
    auto gray = std::make_shared<Material>(Color(0.2, 0.2, 0.2));

    auto main_sphere = std::make_shared<Sphere>(Point3(0, 0.5, 0), 0.5, white);
    auto hole1 = std::make_shared<Sphere>(Point3(0.3, 0.6, 0.3), 0.25, blue);
    auto hole2 = std::make_shared<Sphere>(Point3(-0.3, 0.4, 0.3), 0.2, blue);

    auto csg_shape = csg_difference(
        csg_difference(main_sphere, hole1, white),
        hole2,
        white
    );

    data.scene.add_shape(csg_shape);

    data.scene.add_shape(std::make_shared<Sphere>(Point3(0, -40, 0), 39.5, gray));

    data.scene.add_light(std::make_shared<PointLight>(Point3(5, 8, 5), Color(300, 300, 300)));
    data.scene.add_light(std::make_shared<AmbientLight>(0.15));

    data.camera = Camera(Point3(2, 1.5, 3), Point3(0, 0.3, 0), Vec3(0, 1, 0), 35, 16.0 / 9.0);

    data.scene.build_bvh();
    return data;
}

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n"
              << "Options:\n"
              << "  -s, --scene <file>    Load scene from JSON file\n"
              << "  -o, --output <file>   Output image file (default: output.png)\n"
              << "  -w, --width <n>       Image width (default: 800)\n"
              << "  -h, --height <n>      Image height (default: 450)\n"
              << "  --spp <n>             Samples per pixel (default: 1)\n"
              << "  --builtin <name>      Use built-in scene:\n"
              << "                        - two_spheres (default)\n"
              << "                        - three_spheres\n"
              << "                        - glass\n"
              << "                        - csg\n"
              << "  --help                Show this help message\n";
}

int main(int argc, char* argv[]) {
    std::string scene_file;
    std::string output_file = "output.png";
    std::string builtin_scene = "two_spheres";
    int width = 800;
    int height = 450;
    int spp = 1;
    bool use_path_tracing = false;

    std::vector<std::string> _args;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-s" || arg == "--scene") && i + 1 < argc) {
            scene_file = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_file = argv[++i];
        } else if ((arg == "-w" || arg == "--width") && i + 1 < argc) {
            width = std::stoi(argv[++i]);
        } else if ((arg == "-h" || arg == "--height") && i + 1 < argc) {
            height = std::stoi(argv[++i]);
        } else if (arg == "--spp" && i + 1 < argc) {
            spp = std::stoi(argv[++i]);
        } else if (arg == "--path-trace") {
            use_path_tracing = true;
        } else if (arg == "--builtin" && i + 1 < argc) {
            builtin_scene = argv[++i];
        } else if (arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg[0] != '-') {
            _args.push_back(arg);
        }
    }

    if (!_args.empty() && scene_file.empty()) {
        scene_file = _args[0];
    }
    if (_args.size() >= 2 && output_file == "output.png") {
        output_file = _args[1];
    }

    std::cout << "Resolution: " << width << "x" << height << "\n";
    std::cout << "Samples per pixel: " << spp << "\n";

    SceneLoader::SceneData scene_data;

    try {
        if (!scene_file.empty()) {
            std::ifstream test(scene_file);
            if (!test.good()) {
                std::cerr << "Error: Scene file not found: " << scene_file << "\n";
                return 1;
            }
            test.close();


            std::cout << "Loading scene from: " << scene_file << "\n";
            scene_data = SceneLoader::load(scene_file);
        } else {
            std::cout << "Rendering example scene: " << builtin_scene << "\n";
            if (builtin_scene == "two_spheres") {
                scene_data = create_default_scene();
            } else if (builtin_scene == "three_spheres") {
                scene_data = create_three_spheres_scene();
            } else if (builtin_scene == "glass") {
                scene_data = create_glass_scene();
            } else if (builtin_scene == "csg") {
                scene_data = create_csg_scene();
            } else {
                std::cerr << "unknown scene" << builtin_scene << "\n";
                return 1;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error loading scene: " << e.what() << "\n";
        return 1;
    }

    scene_data.camera.aspect = static_cast<double>(width) / height;
    scene_data.camera = Camera(
        scene_data.camera.eye,
        scene_data.camera.eye + scene_data.camera.w_forward,
        scene_data.camera.v_up,
        std::atan(scene_data.camera.vertical_half_size) * 2.0 * 180.0 / PI,
        scene_data.camera.aspect
    );


    std::vector<Color> framebuffer;
    auto start = std::chrono::high_resolution_clock::now();

    if (use_path_tracing || !scene_data.area_lights.empty()) {
        std::cout << "# of area lights: " << scene_data.area_lights.size() << "\n";

        PathTracer path_tracer;
        path_tracer.samples_per_pixel = spp > 1 ? spp : 64;  // Default to 64 spp for path tracing

        for (const auto& al : scene_data.area_lights) {
            path_tracer.add_area_light(al);
        }

        framebuffer = path_tracer.render(scene_data.camera, scene_data.scene, width, height);
    } else {
        Renderer renderer;
        renderer.samples_per_pixel = spp;
        framebuffer = renderer.render(scene_data.camera, scene_data.scene, width, height);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Render time: " << duration.count() / 1000.0 << " seconds\n";

    if (write_image(output_file, framebuffer, width, height)) {
        std::cout << "Saved to: " << output_file << "\n";
    } else {
        std::cerr << "Failed to save image\n";
        return 1;
    }

    return 0;
}
