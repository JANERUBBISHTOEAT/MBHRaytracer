#ifndef SOLVE_H
#define SOLVE_H

#include "common.h"
#include "ray.h"
#include "vec3.h"
#include <vector>

/** The possible return values from a solution attempt */
enum solve_ret {
    S_GOOD,
    S_ERROR,
    S_SUCC,
};

/** A struct which owns a ray as it travels through space, applying
 * gravitational distortion */
struct ray_iterator {
    struct black_hole {
        double mass;
        double rs;
        point3 origin;
    };

    std::vector<black_hole> holes;
    /** Step size */
    double epsilon;
    /** Debugging parameter used for `tgsl` */
    bool freedom;
    /** Debugging parameter used for `tgsl` */
    bool disable_bh;

    /** Convert the state of the iterator to a string for debugging */
    string fmt();

    /** March the ray forwards, returns the updated ray as an out parameter */
    solve_ret iter(ray *r);

    /**
     * @param holes The black holes in the scene
     * @param initial_ray The initial light ray which will become distorted
     * @param epsilon The epsilon of the simulation
     * @param prevent_freedom Used for `tgsl` to stop the ray once it is out of
     * a certain range
     * @param disable_bh Another parameter for `tgsl` which disables distortion
     */
    ray_iterator(std::vector<black_hole> holes, ray initial_ray,
                 double epsilon, bool prevent_freedom, bool disable_bh);

    /** Dtor */
    ~ray_iterator() = default;

  private:
    ray m_r;                 /* Current ray */
};

#endif
