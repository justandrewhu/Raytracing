#pragma once

#include <vector>
#include <fstream>
#include <cstdint>
#include <string>

namespace rt {

namespace png {

inline uint32_t crc32(const uint8_t* data, size_t len) {
    static uint32_t table[256];
    static bool initialized = false;

    if (!initialized) {
        for (uint32_t i = 0; i < 256; i++) {
            uint32_t c = i;
            for (int k = 0; k < 8; k++) {
                c = (c & 1) ? 0xedb88320 ^ (c >> 1) : c >> 1;
            }
            table[i] = c;
        }
        initialized = true;
    }

    uint32_t crc = 0xffffffff;
    for (size_t i = 0; i < len; i++) {
        crc = table[(crc ^ data[i]) & 0xff] ^ (crc >> 8);
    }
    return crc ^ 0xffffffff;
}

inline uint32_t adler32(const uint8_t* data, size_t len) {
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < len; i++) {
        a = (a + data[i]) % 65521;
        b = (b + a) % 65521;
    }
    return (b << 16) | a;
}

inline void write_be32(std::vector<uint8_t>& out, uint32_t val) {
    out.push_back((val >> 24) & 0xff);
    out.push_back((val >> 16) & 0xff);
    out.push_back((val >> 8) & 0xff);
    out.push_back(val & 0xff);
}

inline void write_be16(std::vector<uint8_t>& out, uint16_t val) {
    out.push_back((val >> 8) & 0xff);
    out.push_back(val & 0xff);
}

inline void write_chunk(std::vector<uint8_t>& out, const char* type, const std::vector<uint8_t>& data) {
    write_be32(out, static_cast<uint32_t>(data.size()));

    size_t type_start = out.size();
    for (int i = 0; i < 4; i++) out.push_back(type[i]);
    for (auto b : data) out.push_back(b);

    uint32_t crc = crc32(&out[type_start], 4 + data.size());
    write_be32(out, crc);
}

}

inline bool write_png(const std::string& filename, const uint8_t* pixels, int width, int height) {
    std::vector<uint8_t> out;

    // PNG signature
    const uint8_t sig[] = {137, 80, 78, 71, 13, 10, 26, 10};
    for (auto b : sig) out.push_back(b);

    // IHDR chunk
    {
        std::vector<uint8_t> ihdr;
        png::write_be32(ihdr, width);
        png::write_be32(ihdr, height);
        ihdr.push_back(8);  
        ihdr.push_back(2);  
        ihdr.push_back(0);  
        ihdr.push_back(0); 
        ihdr.push_back(0);
        png::write_chunk(out, "IHDR", ihdr);
    }

    // IDAT chunk 
    {
        // Build raw image data with filter bytes
        std::vector<uint8_t> raw;
        for (int y = 0; y < height; y++) {
            raw.push_back(0);
            for (int x = 0; x < width; x++) {
                int idx = (y * width + x) * 3;
                raw.push_back(pixels[idx + 0]);
                raw.push_back(pixels[idx + 1]);
                raw.push_back(pixels[idx + 2]);
            }
        }

        // Zlib wrapper with uncompressed deflate blocks
        std::vector<uint8_t> zlib;
        zlib.push_back(0x78);  
        zlib.push_back(0x01); 

        // Split into 65535-byte blocks (max uncompressed block size)
        size_t pos = 0;
        while (pos < raw.size()) {
            size_t block_size = std::min(raw.size() - pos, size_t(65535));
            bool last = (pos + block_size >= raw.size());

            zlib.push_back(last ? 0x01 : 0x00);  // BFINAL + BTYPE=00
            uint16_t len = static_cast<uint16_t>(block_size);
            uint16_t nlen = ~len;
            zlib.push_back(len & 0xff);
            zlib.push_back((len >> 8) & 0xff);
            zlib.push_back(nlen & 0xff);
            zlib.push_back((nlen >> 8) & 0xff);

            for (size_t i = 0; i < block_size; i++) {
                zlib.push_back(raw[pos + i]);
            }
            pos += block_size;
        }

        // Adler-32 checksum
        uint32_t adler = png::adler32(raw.data(), raw.size());
        png::write_be32(zlib, adler);

        png::write_chunk(out, "IDAT", zlib);
    }

    // IEND chunk
    png::write_chunk(out, "IEND", {});

    // Write to file
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) return false;
    file.write(reinterpret_cast<char*>(out.data()), out.size());
    return true;
}

} // namespace rt
