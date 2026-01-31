#pragma once

#include "vec3.hpp"
#include "ray.hpp"
#include <cmath>

namespace rt {

class Camera {
public:
    Point3 eye;
    Vec3 w_forward;  // forward direction
    Vec3 u_right;    // right direction
    Vec3 v_up;       // up direction
    double f;        // Focal length
    double vertical_half_size;
    double horizontal_half_size;
    double aspect;

    Camera(const Point3& eye = Point3(0, 0, 0),
           const Point3& target = Point3(0, 0, -1),
           const Vec3& up = Vec3(0, 1, 0),
           double vfov = 90.0,
           double aspect = 1.0)
        : eye(eye), aspect(aspect), f(1.0)
    {
        double theta = vfov * PI / 180.0;
        vertical_half_size = std::tan(theta * 0.5) * f;
        horizontal_half_size = aspect * vertical_half_size;

        w_forward = normalize(target - eye);
        u_right = normalize(cross(w_forward, up));
        v_up = normalize(cross(u_right, w_forward));
    }

    Ray generate_ray(double img_x, double img_y) const {
        // img_x, img_y are in [0,1] where (0,0) is upper left
        double screen_x = 2.0 * img_x - 1.0;
        double screen_y = 1.0 - 2.0 * img_y;

        Vec3 direction = w_forward * f
                       + u_right * (screen_x * horizontal_half_size)
                       + v_up * (screen_y * vertical_half_size);

        return Ray(eye, direction, 0.0, INF);
    }
};

} // namespace rt
