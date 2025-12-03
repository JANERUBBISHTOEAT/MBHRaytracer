#include "solve.h"
#include "ray.h"
#include "vec3.h"
#include <utility>

string ray_iterator::fmt() {
    stringstream s;
    s << "m_r: " << m_r.fmt() << "\n";
    s << "bh count: " << holes.size();
    return s.str();
}

solve_ret ray_iterator::iter(ray *r) {
    point3 pos = m_r.at(0);
    vec3 dir = unit_vector(m_r.direction());

    if (!disable_bh) {
        vec3 total_accel(0, 0, 0);
        for (const auto &bh : holes) {
            vec3 rel = pos - bh.origin;
            double dist = rel.length();

            if (dist <= bh.rs) {
                return S_SUCC;
            }

            if (dist <= 0) {
                return S_ERROR;
            }

            double inv_dist_cubed = 1.0 / (dist * dist * dist);
            total_accel += (-1 * bh.mass * rel) * inv_dist_cubed;
        }
        dir = unit_vector(dir + total_accel * epsilon);
    }

    point3 next_pos = pos + dir * epsilon;

    if (!freedom && next_pos.length() > 200) {
        fprintf(stderr, "freedom denied\n");
        return S_ERROR;
    }

    m_r = ray(next_pos, dir);
    *r = m_r;
    return S_GOOD;
}

ray_iterator::ray_iterator(std::vector<black_hole> holes, ray initial_ray,
                           double epsilon, bool prevent_freedom,
                           bool disable_bh)
    : holes(std::move(holes)), epsilon(epsilon), freedom(prevent_freedom),
      disable_bh(disable_bh) {
    m_r = ray(initial_ray.origin(), unit_vector(initial_ray.direction()));

    for (auto &hole : this->holes) {
        double G = 6.6743;
        hole.rs = 7 * G * hole.mass / 10;
    }
}
