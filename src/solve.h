#ifndef SOLVE_H
#define SOLVE_H

#include "common.h"
#include "ray.h"
#include "vec3.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
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
    
    struct multi_bh_params {
        std::vector<black_hole> *holes;
    };

    /** BH mass (for single BH mode with GSL) */
    double mass;
    /** schwartzschield radius (for single BH mode) */
    double rs;
    /** BH origin (for single BH mode) */
    point3 origin;
    /** Black holes (for multi-BH mode) */
    std::vector<black_hole> holes;
    /** Step size */
    double epsilon;
    /** Debugging parameter used for `tgsl` */
    bool freedom;
    /** Debugging parameter used for `tgsl` */
    bool disable_bh;
    /** Use simple integration for multi-BH */
    bool use_simple_integration;
    /** Use GSL high-precision integration for multi-BH */
    bool use_gsl_integration;

    /** Convert the state of the iterator to a string for debugging */
    string fmt();

    /** March the ray forwards, returns the updated ray as an out parameter */
    solve_ret iter(ray *r);

    /**
     * @param mass The mass of the black hole
     * @param initial_ray The initial light ray which will become distorted
     * @param origin The origin of the black hole
     * @param epsilon The epsilon of the simulation
     * @param prevent_freedom Used for `tgsl` to stop the ray once it is out of
     * a certain range
     * @param disable_bh Another parameter for `tgsl` which disables distortion
     */
    ray_iterator(double mass, ray initial_ray, point3 origin, double epsilon,
                 bool prevent_freedom, bool disable_bh);

    /**
     * Multi-black hole constructor using simple integration
     * @param holes The black holes in the scene
     * @param initial_ray The initial light ray which will become distorted
     * @param epsilon The epsilon of the simulation
     * @param prevent_freedom Used for `tgsl` to stop the ray once it is out of
     * a certain range
     * @param disable_bh Another parameter for `tgsl` which disables distortion
     * @param use_gsl If true, use GSL Runge-Kutta solver for high precision
     */
    ray_iterator(std::vector<black_hole> holes, ray initial_ray,
                 double epsilon, bool prevent_freedom, bool disable_bh,
                 bool use_gsl = false);

    /** Dtor */
    ~ray_iterator() { 
        // Free GSL driver if it was allocated
        // For single BH: use_simple_integration=false, use_gsl_integration=false, but m_d is allocated
        // For multi-BH GSL: use_simple_integration=false, use_gsl_integration=true, m_d is allocated
        // For multi-BH simple: use_simple_integration=true, use_gsl_integration=false, m_d is nullptr
        if (m_d != nullptr) {
            gsl_odeiv2_driver_free(m_d);
            m_d = nullptr;
        }
    }

  private:
    point3 transfer_out(point3 pt);
    double m_z_rot;
    ray m_r;                 /* Current ray */
    double m_t;              /* alias for current phi */
    double m_y[2];           /* (u,phi) for single BH */
    double m_state[6];       /* (x,y,z,vx,vy,vz) for multi-BH GSL */
    gsl_odeiv2_system m_sys; /*= {func, jac, 2, &mass}; */
    gsl_odeiv2_driver *m_d;  /* gsl_odeiv2_driver_apply (d, &t, new_phi, y); */
    gsl_odeiv2_system m_multi_sys; /* For multi-BH GSL */
    multi_bh_params m_multi_params;
};

#endif
