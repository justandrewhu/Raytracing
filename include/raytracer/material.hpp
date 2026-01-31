#pragma once

#include "vec3.hpp"
#include <memory>

namespace rt {

struct Material {
    Color k_d;      // Diffuse coefficient
    Color k_s;      // Specular coefficient
    double p;       // specular exponent 
    Color k_m;      // mirror reflection coefficient
    Color k_a;      // ambient coefficient
    double ior;     // ior
    Color k_t;      // Transmission coefficient
    Color k_e;      // Emission coefficient

    Material()
        : k_d(0.5), k_s(0.0), p(20.0), k_m(0.0), k_a(0.5),
          ior(1.0), k_t(0.0), k_e(0.0) {}

    Material(const Color& diffuse, double specular = 0.0, double exponent = 20.0,
             double mirror = 0.0, const Color& ambient = Color(-1))
        : k_d(diffuse), k_s(specular), p(exponent), k_m(mirror),
          k_a(ambient.x < 0 ? diffuse : ambient), ior(1.0), k_t(0.0), k_e(0.0) {}

    // Builder pattern methods for cleaner material construction
    Material& set_specular(const Color& spec) { k_s = spec; return *this; }
    Material& set_specular(double spec) { k_s = Color(spec); return *this; }
    Material& set_exponent(double exp) { p = exp; return *this; }
    Material& set_mirror(const Color& m) { k_m = m; return *this; }
    Material& set_mirror(double m) { k_m = Color(m); return *this; }
    Material& set_ior(double i) { ior = i; return *this; }
    Material& set_transmission(const Color& t) { k_t = t; return *this; }
    Material& set_transmission(double t) { k_t = Color(t); return *this; }
    Material& set_emission(const Color& e) { k_e = e; return *this; }
    Material& set_emission(double e) { k_e = Color(e); return *this; }
    Material& set_ambient(const Color& a) { k_a = a; return *this; }

    bool has_reflection() const {
        return k_m.x > EPSILON || k_m.y > EPSILON || k_m.z > EPSILON;
    }

    bool has_transmission() const {
        return k_t.x > EPSILON || k_t.y > EPSILON || k_t.z > EPSILON;
    }

    bool has_emission() const {
        return k_e.x > EPSILON || k_e.y > EPSILON || k_e.z > EPSILON;
    }
};

using MaterialPtr = std::shared_ptr<Material>;

inline MaterialPtr make_material(const Color& diffuse, double specular = 0.0, double exponent = 20.0,
                                  double mirror = 0.0) {
    return std::make_shared<Material>(diffuse, specular, exponent, mirror);
}

} // namespace rt
