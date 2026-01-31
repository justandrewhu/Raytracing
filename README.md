# C++ Raytracer

This is a feature-complete raytracer, originally written in Python but ported to C++17. It supports meshes, CSG operations, BVH acceleration, reflection, refraction,   e and JSON scene loading.

## Features

- **Primitives**: Spheres, Triangles, Meshes (OBJ loading)
- **Materials**: Diffuse, specular, mirror reflection, refraction with Fresnel
- **Acceleration**: BVH (Bounding Volume Hierarchy)
- **CSG**: Union, Intersection, Difference operations
- **Lighting**: Point lights, ambient lights
- **Output**: PPM format (PNG with optional stb_image_write)
- **Scene Loading**: JSON scene files

## Building

### Using CMake (recommended)
```bash
mkdir build && cd build
cmake ..
make
```

### Using Make directly
```bash
make
```

### Manual compilation
```bash
g++ -std=c++17 -O3 -I include -o raytracer src/main.cpp
```

## Usage

```bash
# Render built-in scenes
./raytracer --builtin two_spheres -w 800 -h 450 -o output.ppm
./raytracer --builtin three_spheres -w 800 -h 450 -o output.ppm
./raytracer --builtin glass -w 800 -h 450 -o output.ppm
./raytracer --builtin csg -w 800 -h 450 -o output.ppm

# Render from JSON scene file
./raytracer -s scenes/example.json -w [width_resolution] -h [height resolution] -o output.ppm

# Multi-sample anti-aliasing
./raytracer --builtin glass -w 800 -h 450 --spp 4 -o output.ppm
```

### Command Line Options

| Option | Description |
|--------|-------------|
| `-s, --scene <file>` | Load scene from JSON file |
| `-o, --output <file>` | Output image file (default: output.ppm) |
| `-w, --width <n>` | Image width (default: 800) |
| `-h, --height <n>` | Image height (default: 450) |
| `--spp <n>` | Samples per pixel (default: 1) |
| `--builtin <name>` | Use built-in scene |
| `--help` | Show help message |

### Built-in Scenes

- `two_spheres` - Simple scene with one sphere on a ground plane
- `three_spheres` - Three spheres with mirror reflections
- `glass` - Transparent glass sphere with refraction
- `csg` - CSG difference operation demo

## JSON Scene Format

```json
{
    "background": [0.2, 0.3, 0.5],

    "camera": {
        "eye": [3, 1.7, 5],
        "target": [0, 0, 0],
        "up": [0, 1, 0],
        "vfov": 25,
        "aspect": 1.777
    },

    "materials": {
        "myMaterial": {
            "diffuse": [0.7, 0.7, 0.4],
            "specular": 0.6,
            "exponent": 20,
            "mirror": 0.0,
            "transmission": [0, 0, 0],
            "ior": 1.0
        }
    },

    "shapes": [
        {
            "type": "sphere",
            "center": [0, 0, 0],
            "radius": 0.5,
            "material": "myMaterial"
        },
        {
            "type": "triangle",
            "v0": [0, 0, 0],
            "v1": [1, 0, 0],
            "v2": [0, 1, 0],
            "material": "myMaterial"
        },
        {
            "type": "mesh",
            "file": "model.obj",
            "material": "myMaterial"
        },
        {
            "type": "csg_difference",
            "left": { "type": "sphere", "center": [0, 0, 0], "radius": 1 },
            "right": { "type": "sphere", "center": [0.5, 0, 0], "radius": 0.5 },
            "material": "myMaterial"
        }
    ],

    "lights": [
        {
            "type": "point",
            "position": [12, 10, 5],
            "intensity": [300, 300, 300]
        },
        {
            "type": "ambient",
            "intensity": 0.1
        }
    ]
}
```

## Adding in Custom Meshes
You can add custom obj files into the models folder. Then, upon defining an object in a json scene, use `type: "mesh"` and set the path to be `models.name_of_mesh.obj`. 

## Extending

### Adding New Shapes

1. Create a new header in `include/raytracer/`
2. Inherit from `Shape` and implement:
   - `Hit intersect(const Ray& ray) const`
   - `AABB bounds() const`
   - Optionally: `std::vector<Hit> all_hits(const Ray& ray) const` for CSG

### Adding PNG Output

Download `stb_image_write.h` from [stb](https://github.com/nothings/stb) and place in `external/`. Then compile with `-DUSE_STB_IMAGE_WRITE`.

## License

MIT License
