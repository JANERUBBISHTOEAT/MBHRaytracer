#include "solve.h"
#include "ray.h"
#include "vec3.h"
#include <utility>

// ODE function for single black hole (spherical coordinates)
static int func(double t, const double y[], double f[], void *params) {
    (void)(t); /* avoid unused parameter warning */
    double mass = *(double *)params;
    // we cannot traverse through the wormhole.
    if (y[0] < 0.0) {
        return GSL_FAILURE;
    }
    f[0] = y[1];
    f[1] = 3 * mass * y[0] * y[0] - y[0];
    return GSL_SUCCESS;
}

// ODE function for multiple black holes (Cartesian coordinates)
// State: [x, y, z, vx, vy, vz]
static int func_multi(double t, const double y[], double f[], void *params) {
    (void)(t);
    ray_iterator::multi_bh_params *p = (ray_iterator::multi_bh_params *)params;
    
    // Position
    point3 pos(y[0], y[1], y[2]);
    
    // Velocity derivatives (acceleration)
    // Calculate gravitational acceleration from all black holes
    vec3 total_accel(0, 0, 0);
    for (const auto &bh : *(p->holes)) {
        vec3 rel = pos - bh.origin;
        double dist = rel.length();
        
        // Single BH func() only checks y[0] < 0.0, NOT event horizon
        // We do the same - only check for invalid distance, let iter() check event horizon
        if (dist <= 0) {
            return GSL_FAILURE;
        }
        
        // Newtonian gravitational acceleration: a = -GM * r / r^3
        // For light rays, this is an approximation (should use GR, but simpler for multi-BH)
        // Near event horizon, use higher precision calculation
        double dist_sq = dist * dist;
        double inv_dist_cubed = 1.0 / (dist_sq * dist);
        total_accel += (-1 * bh.mass * rel) * inv_dist_cubed;
    }
    
    // d/dt [x, y, z] = [vx, vy, vz]
    f[0] = y[3];  // dx/dt = vx
    f[1] = y[4];  // dy/dt = vy
    f[2] = y[5];  // dz/dt = vz
    
    // d/dt [vx, vy, vz] = acceleration
    f[3] = total_accel.x();  // dvx/dt = ax
    f[4] = total_accel.y();  // dvy/dt = ay
    f[5] = total_accel.z();  // dvz/dt = az
    
    return GSL_SUCCESS;
}

static int jac(double t, const double y[], double *dfdy, double dfdt[],
               void *params) {
    (void)(t); /* avoid unused parameter warning */
    double mass = *(double *)params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;
    gsl_matrix_set(m, 0, 0, 0.0);
    gsl_matrix_set(m, 0, 1, 1.0);
    gsl_matrix_set(m, 1, 0, 6 * mass * y[0] - 1);
    gsl_matrix_set(m, 1, 1, 0);
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}

static double find_z_rot(point3 start) {
    double y = start.y();
    double x = start.x();
    y *= -1;
    // fix atan stupidity returning -pi when not wanted or needed
    if (y == -0.0) {
        return 0;
    }
    double theta = atan2(y, x);
    return theta;
}

point3 ray_iterator::transfer_out(point3 pt) {
    assert(fabs(pt.e[2]) < 0.000001);
    pt.e[2] = -1 * pt.e[1];
    pt.e[1] = 0.0;
    return pt.rotz(-1 * m_z_rot) + origin;
}

string ray_iterator::fmt() {
    stringstream s;
    if (use_simple_integration) {
        s << "m_r: " << m_r.fmt() << "\n";
        s << "bh count: " << holes.size();
    } else {
        s << "m_r: " << m_r.fmt() << "\nu: " << m_y[0] << "\nu': " << m_y[1]
          << "\nm_t:" << m_t;
    }
    return s.str();
}

solve_ret ray_iterator::iter(ray *r) {
    if (use_simple_integration && !use_gsl_integration) {
        // Simple Euler integration for multiple black holes (first-order, fast but less accurate)
        point3 pos = m_r.at(0);
        vec3 dir = unit_vector(m_r.direction());

        if (!disable_bh) {
            vec3 total_accel(0, 0, 0);
            for (const auto &bh : holes) {
                vec3 rel = pos - bh.origin;
                double dist = rel.length();

                // Check for event horizon - use exact radius like single BH
                if (dist <= bh.rs) {
                    return S_SUCC;
                }

                if (dist <= 0) {
                    return S_ERROR;
                }

                // Calculate acceleration with higher precision
                double dist_sq = dist * dist;
                double inv_dist_cubed = 1.0 / (dist_sq * dist);
                total_accel += (-1 * bh.mass * rel) * inv_dist_cubed;
            }
            // Euler method: v_new = v_old + a * dt, then normalize (light speed is constant)
            dir = unit_vector(dir + total_accel * epsilon);
        }

        // Update position: x_new = x_old + v * dt
        point3 next_pos = pos + dir * epsilon;

        if (!freedom && next_pos.length() > 200) {
            fprintf(stderr, "freedom denied\n");
            return S_ERROR;
        }

        m_r = ray(next_pos, dir);
        *r = m_r;
        return S_GOOD;
    }
    
    if (use_gsl_integration) {
        // GSL Runge-Kutta high-precision integration for multiple black holes
        point3 pos = m_r.at(0);
        vec3 dir = unit_vector(m_r.direction());
        
        // Single BH method: check event horizon AFTER GSL call, not before
        // We do the same - let GSL integrate first, then check capture
        if (!disable_bh) {
            double t = 0.0;
            double t1 = epsilon;
            
            // Current state: [x, y, z, vx, vy, vz]
            // Store position and velocity (direction * speed, where speed = 1 for light)
            m_state[0] = pos.x();
            m_state[1] = pos.y();
            m_state[2] = pos.z();
            m_state[3] = dir.x();  // vx (light speed normalized to 1)
            m_state[4] = dir.y();  // vy
            m_state[5] = dir.z();  // vz
            
            // Single BH method: no nmax set, allowing unlimited internal steps
            // We do the same - don't set nmax to match single BH behavior exactly
            int status = gsl_odeiv2_driver_apply(m_d, &t, t1, m_state);
            
            // Single BH method: check event horizon AFTER GSL call: if (1 / m_y[0] <= rs)
            // We do the same - check AFTER GSL integration, using exact radius
            point3 next_pos(m_state[0], m_state[1], m_state[2]);
            vec3 next_dir(m_state[3], m_state[4], m_state[5]);
            
            // Check for black hole capture - use exact radius like single BH
            // Single BH uses: if (1 / m_y[0] <= rs), we use: if (dist <= bh.rs)
            for (const auto &bh : holes) {
                vec3 rel = next_pos - bh.origin;
                double dist = rel.length();
                if (dist <= bh.rs) {
                    return S_SUCC;
                }
            }
            
            // Single BH method: if GSL error, check capture then return error
            // We do the same - simple error handling
            if (status != GSL_SUCCESS) {
                // Check if captured despite error (shouldn't happen, but check anyway)
                for (const auto &bh : holes) {
                    vec3 rel = next_pos - bh.origin;
                    double dist = rel.length();
                    if (dist <= bh.rs) {
                        return S_SUCC;
                    }
                }
                // Like single BH: return error
                fprintf(stderr, "error, return value=%d\n", status);
                return S_ERROR;
            }
            
            // Single BH method: calculate new_direction from iter_pos - pos
            // iter_pos comes from GSL result in spherical coords: spher3(1/m_y[0], ...)
            // This gives more accurate direction near event horizon
            // We do the same - use GSL result position to calculate direction
            // This ensures we capture the full distortion from GSL integration
            vec3 new_direction = unit_vector(next_pos - pos);
            
            // Single BH also checks: if (1 / m_y[0] > 200 && !freedom) return error
            // We do similar check for distance from origin
            if (!freedom && next_pos.length() > 200) {
                fprintf(stderr, "freedom denied\n");
                return S_ERROR;
            }
            
            // Single BH uses: m_r = ray(iter_pos, new_direction)
            // We do the same - use GSL result position and calculated direction
            m_r = ray(next_pos, new_direction);
            *r = m_r;
            return S_GOOD;
        } else {
            // No black hole effects
            point3 next_pos = pos + dir * epsilon;
            m_r = ray(next_pos, dir);
            *r = m_r;
            return S_GOOD;
        }
    }
    
    // Original GSL-based method for single black hole
    double new_phi;
    point3 pos = m_r.at(0);
    assert(m_r.orig.z() == 0);
    assert(m_r.dir.z() == 0);
    assert(pos.length() == pos.length());

    point3 next_pos = m_r.at(epsilon);
    spher3 s_next_pos = next_pos.to_spher3();

    new_phi = s_next_pos.z();

    if (!disable_bh) {
        int status;
        if (new_phi < m_t) {
            fprintf(stderr, "weird integration bounds\n");
            return S_ERROR;
        }
        status = gsl_odeiv2_driver_apply(m_d, &m_t, new_phi, m_y);

        if (status != GSL_SUCCESS) {
            fprintf(stderr, "error, return value=%d\n", status);
            return S_ERROR;
        }
        // update current sim location
        m_t = new_phi;

        if (1 / m_y[0] <= rs) {
            return S_SUCC;
        } else if (1 / m_y[0] > 200 && !freedom) {
            fprintf(stderr, "freedom denied\n");
            return S_ERROR;
        }
    }

    point3 iter_pos;
    point3 new_direction;
    if (!disable_bh) {
        assert(fabs(s_next_pos.y() - pi / 2) < 0.00001);
        iter_pos = spher3(1 / m_y[0], s_next_pos.y(), m_t).to_cartesian();
        assert(fabs(iter_pos.z()) < 0.00001);
        new_direction = unit_vector(iter_pos - pos);
        // new_direction = iter_pos - pos;
        // nan
        assert(new_direction.y() == new_direction.y());
    } else {
        // needed to avoid epsilon influencing direction magnitude
        new_direction = unit_vector(next_pos - pos);
        iter_pos = next_pos;
    }
    point3 ret_dir;
    ret_dir = new_direction;
    ret_dir.e[2] = -1 * ret_dir.e[1];
    ret_dir.e[1] = 0.0;
    *r = ray(transfer_out(iter_pos), ret_dir.rotz(-1 * m_z_rot));
    assert(iter_pos.length() == iter_pos.length() &&
           new_direction.length() == new_direction.length());
    assert(fabs(iter_pos.z()) < 0.00001);
    assert(fabs(new_direction.z()) < 0.00001);
    iter_pos.e[2] = 0;
    new_direction.e[2] = 0;
    m_r = ray(iter_pos, new_direction);
    return S_GOOD;
}

ray_iterator::ray_iterator(double mass, ray initial_ray, point3 origin,
                           double epsilon, bool prevent_freedom,
                           bool disable_bh)
    : mass(mass), origin(origin), epsilon(epsilon), freedom(prevent_freedom),
      disable_bh(disable_bh), use_simple_integration(false), 
      use_gsl_integration(false), m_d(nullptr) {
    point3 start = initial_ray.origin() - origin;
    // cannot start a light ray at the center of a black hole.
    assert(start.length() > 0);
    vec3 direction = initial_ray.direction();
    m_z_rot = find_z_rot(direction);
    direction = direction.rotz(m_z_rot);
    start = start.rotz(m_z_rot);
    assert(fabs(start.y()) < 0.0000001);
    assert(fabs(direction.y()) < 0.0000001);
    start.e[1] = -1 * start.e[2];
    start.e[2] = 0;
    direction.e[1] = -1 * direction.e[2];
    direction.e[2] = 0;

    // checking that our system actually works for extrapolated points
    assert((transfer_out(start + direction) - initial_ray.at(1)).length() <
           0.001);

    m_r = ray(start, direction);

    assert(fabs(start.z()) < 0.000001);
    assert(fabs(direction.z()) < 0.000001);

    point3 next = m_r.at(epsilon);
    m_r.orig = next;

    spher3 s_start = start.to_spher3();
    spher3 s_next = next.to_spher3();

    double initial_r = s_start.x();
    double next_r = s_next.x();
    double initial_phi = s_start.z();
    double next_phi = s_next.z();

    double initial_velo =
        (1 / next_r - 1 / initial_r) / (next_phi - initial_phi);
    assert(2 * initial_velo != initial_velo); /* Check for nan and inf */

    m_sys = {func, jac, 2, &this->mass};
    // hstart is initial step size. the other two arguments are error.
    //                                                                hstart
    //                                                                epsabs
    //                                                                epsrel
    m_d = gsl_odeiv2_driver_alloc_y_new(&m_sys, gsl_odeiv2_step_rk4, 1e-4, 1e-4,
                                        0.0);

    m_t = initial_phi;
    // if (m_t < 0)
    //   m_t += 2*pi;
    m_y[0] = 1 / initial_r;
    m_y[1] = initial_velo;

    double G = 6.6743;
    rs = 7 * G * mass / 10;
}

// Multi-black hole constructor
ray_iterator::ray_iterator(std::vector<black_hole> holes, ray initial_ray,
                           double epsilon, bool prevent_freedom,
                           bool disable_bh, bool use_gsl)
    : holes(std::move(holes)), epsilon(epsilon), freedom(prevent_freedom),
      disable_bh(disable_bh), use_simple_integration(!use_gsl), 
      use_gsl_integration(use_gsl), m_d(nullptr) {
    m_r = ray(initial_ray.origin(), unit_vector(initial_ray.direction()));

    double G = 6.6743;
    for (auto &hole : this->holes) {
        hole.rs = 7 * G * hole.mass / 10;
    }
    
    if (use_gsl_integration) {
        // Initialize GSL ODE solver for multi-BH
        m_multi_params.holes = &this->holes;
        m_multi_sys = {func_multi, nullptr, 6, &m_multi_params};
        // Use EXACT same method as single BH: rk4 with same tolerances
        // Single BH uses: gsl_odeiv2_step_rk4, 1e-4, 1e-4, 0.0
        // We use the same to match single BH behavior exactly
        // The key is that single BH uses spherical coords with u=1/r variable
        // which has better numerical properties near event horizon
        // For Cartesian coords, we need to rely on GSL's adaptive step size
        // with the same tolerances to achieve similar precision
        m_d = gsl_odeiv2_driver_alloc_y_new(&m_multi_sys, gsl_odeiv2_step_rk4, 
                                           1e-4, 1e-4, 0.0);
        m_t = 0.0;
    } else {
        m_d = nullptr;
    }
}
