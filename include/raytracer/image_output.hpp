#pragma once

#include "vec3.hpp"
#include "png_writer.hpp"
#include <vector>
#include <fstream>
#include <string>
#include <cstdint>

namespace rt {

// PNG output (no external dependencies)
inline bool write_png(const std::string& filename, const std::vector<Color>& framebuffer,
                      int width, int height, bool srgb = true) {
    std::vector<uint8_t> pixels(width * height * 3);

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            int idx = j * width + i;
            Color c = framebuffer[idx];
            if (srgb) {
                c = to_srgb(c);
            }

            pixels[idx * 3 + 0] = static_cast<uint8_t>(std::clamp(c.x * 255.0, 0.0, 255.0));
            pixels[idx * 3 + 1] = static_cast<uint8_t>(std::clamp(c.y * 255.0, 0.0, 255.0));
            pixels[idx * 3 + 2] = static_cast<uint8_t>(std::clamp(c.z * 255.0, 0.0, 255.0));
        }
    }

    return rt::write_png(filename, pixels.data(), width, height);
}

// PPM output (no external dependencies)
inline bool write_ppm(const std::string& filename, const std::vector<Color>& framebuffer,
                      int width, int height, bool srgb = true) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }

    file << "P6\n" << width << " " << height << "\n255\n";

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            Color c = framebuffer[j * width + i];
            if (srgb) {
                c = to_srgb(c);
            }

            uint8_t r = static_cast<uint8_t>(std::clamp(c.x * 255.0, 0.0, 255.0));
            uint8_t g = static_cast<uint8_t>(std::clamp(c.y * 255.0, 0.0, 255.0));
            uint8_t b = static_cast<uint8_t>(std::clamp(c.z * 255.0, 0.0, 255.0));

            file.write(reinterpret_cast<char*>(&r), 1);
            file.write(reinterpret_cast<char*>(&g), 1);
            file.write(reinterpret_cast<char*>(&b), 1);
        }
    }

    return true;
}


// Unified write function
inline bool write_image(const std::string& filename, const std::vector<Color>& framebuffer,
                        int width, int height, bool srgb = true) {
    // Check extension
    size_t dot_pos = filename.rfind('.');
    if (dot_pos != std::string::npos) {
        std::string ext = filename.substr(dot_pos);
        for (auto& c : ext) c = std::tolower(c);

        if (ext == ".png") {
            return write_png(filename, framebuffer, width, height, srgb);
        }
        if (ext == ".ppm") {
            return write_ppm(filename, framebuffer, width, height, srgb);
        }
    }

    // Default to PNG
    return write_png(filename + ".png", framebuffer, width, height, srgb);
}

} // namespace rt
