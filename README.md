# Raytracing

This is a custom raytracer, originally written in Python but ported to C++17. It supports meshes, CSG operations, BVH acceleration, reflection, refraction, multiple importance sampling, area lights, (and eventually more).

![gif](examples/glassbunnyanim1/glassbunny.gif)
A render of the stanford bunny with an ice-like material. Also in this image is a good example of refraction (the crystal-like object) and fresnel reflection (the glass ball). If you can't tell I like crystal/ice-like effects.



## Beginnings
There are 2 main primers, 3D transformation math, and the physics behind natural lighting phenomenom. Much of the transfomration math is used a lot in graphics to define the geometry needed to actually conceptualize 3D space. 

Raytracing is an illumination system modeled on the simple fact that to see things, light must bounce from that thing into our eye. For each pixel, we can look up ray intersections with geometry. For each intersection we also want to generate rays to light sources, and figure out how the light illuminates that surface. 


## Refraction and Fresnel Reflection
Snell's law, IORs, and Fresnel-Schlick approximation



## Acceleration
While it can be computationally expensive to find the nearest triangle for every single ray, this can be sped up in many ways. A popular method is the Bounding Volume Hierchy, where triangles in the scene are giving some bounding volume (cube, sphere, axis aligned or non-axis aligned prism), and before we test intersections with the triangle itself we test for the bounding volume which is usually a simpler computation (ideally small cost for bv test, larger cost for obj test). In addition, the boxes are split up into recursively structured groups, where bonding volumes can contain other bounding volumes, and at the lowest level a list of objects/triangles, hence the "Hierchy" part. It is split in a way such that if you miss some bounding volume, you also miss all of its bounding volume children, to make a tree like structure that is more efficient (log(n) vs. n) to traverse rather than iterively checking.

### BVH tradeoffs and BVH construction








## Annoying build stuff

### Using CMake
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

You can create your own scenes just through defining objects and materials in a json file. Materials support different types of reflection models, including diffuse, specular, blinn-phong. There is also transparency controlled by transmission, an index of refraction you can tune to create different glass/prism-like effects, and even a mirror coefficient that doesn't look trippy at all when thrown onto a sphere. Certain **primitive** shapes are supported, but if you want to get really fancy, I'd suggest downloading 3D model(s) in the form of obj files (or making your own through software like Blender) and try to achieve an aesthetic looking scene using them.

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

### Adding New Shapes (obj files are easier though...)

1. Create a new header in `include/raytracer/`
2. Inherit from `Shape` and implement:
   - `Hit intersect(const Ray& ray) const`
   - `AABB bounds() const`
   - Optionally: `std::vector<Hit> all_hits(const Ray& ray) const` for CSG

### Adding PNG Output

Download `stb_image_write.h` from [stb](https://github.com/nothings/stb) and place in `external/`. Then compile with `-DUSE_STB_IMAGE_WRITE`.

## License

MIT License
