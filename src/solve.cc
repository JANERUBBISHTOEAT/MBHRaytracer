#include "solve.h"
#include "ray.h"
#include "vec3.h"
#include <utility>

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

static int multi_func(double t, const double y[], double f[], void *params) {
    (void)(t);
    auto *p = static_cast<ray_iterator::multi_params *>(params);

    point3 pos(y[0], y[1], y[2]);
    vec3 vel(y[3], y[4], y[5]);
    vec3 accel(0, 0, 0);

    for (const auto &bh : *p->holes) {
        vec3 rel = pos - bh.origin;
        double dist = rel.length();

        if (dist <= bh.rs) {
            dist = bh.rs;
        }

        double inv_dist_cubed = 1.0 / (dist * dist * dist);
        accel += (-1 * bh.mass * rel) * inv_dist_cubed;
    }

    f[0] = vel.x();
    f[1] = vel.y();
    f[2] = vel.z();
    f[3] = accel.x();
    f[4] = accel.y();
    f[5] = accel.z();

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
    if (mode == Mode::SingleGsl) {
        s << "m_r: " << m_r.fmt() << "\nu: " << m_y[0] << "\nu': " << m_y[1]
          << "\nm_t:" << m_t;
    } else {
        s << "m_r: " << m_r.fmt() << "\nm_t:" << m_t;
        s << "\nbh count: " << holes.size();
    }
    return s.str();
}

solve_ret ray_iterator::iter(ray *r) {
    if (mode == Mode::MultiGsl) {
        point3 pos(m_state[0], m_state[1], m_state[2]);
        vec3 vel(m_state[3], m_state[4], m_state[5]);

        // Check if we're close to any black hole (within threshold)
        // If so, use Schwarzschild-corrected acceleration
        const double schwarz_threshold_factor = 10.0;
        bool use_schwarz_correction = false;
        double schwarz_correction_factor = 1.0;
        
        for (const auto &bh : holes) {
            double dist = (pos - bh.origin).length();
            if (dist <= bh.rs) {
                return S_SUCC;
            }
            // Check if within Schwarzschild threshold
            if (dist < schwarz_threshold_factor * bh.rs) {
                // Apply Schwarzschild correction: a_schwarz ≈ a_newton * (1 + 3*rs/r)
                // This is a first-order correction for weak fields
                double correction = 1.0 + 3.0 * bh.rs / dist;
                schwarz_correction_factor *= correction;
                use_schwarz_correction = true;
            }
        }

        // Use modified multi_func that applies Schwarzschild correction
        if (use_schwarz_correction) {
            // Manually compute acceleration with Schwarzschild correction
            vec3 accel(0, 0, 0);
            for (const auto &bh : holes) {
                vec3 rel = pos - bh.origin;
                double dist = rel.length();
                
                if (dist <= bh.rs) {
                    dist = bh.rs;
                }
                
                double inv_dist_cubed = 1.0 / (dist * dist * dist);
                double correction = 1.0;
                
                // Apply Schwarzschild correction if close
                if (dist < schwarz_threshold_factor * bh.rs) {
                    correction = 1.0 + 3.0 * bh.rs / dist;
                }
                
                accel += (-1 * bh.mass * rel) * inv_dist_cubed * correction;
            }
            
            // Update state using Euler method with corrected acceleration
            point3 new_pos = pos + vel * epsilon;
            vec3 new_vel = vel + accel * epsilon;
            
            // Normalize velocity direction (light speed is constant)
            new_vel = unit_vector(new_vel);
            
            // Update state
            m_state[0] = new_pos.x();
            m_state[1] = new_pos.y();
            m_state[2] = new_pos.z();
            m_state[3] = new_vel.x();
            m_state[4] = new_vel.y();
            m_state[5] = new_vel.z();
            m_t += epsilon;
            
            pos = new_pos;
            
            if (!freedom && pos.length() > 200) {
                fprintf(stderr, "freedom denied\n");
                return S_ERROR;
            }
            
            for (const auto &bh : holes) {
                double dist = (pos - bh.origin).length();
                if (dist <= bh.rs) {
                    return S_SUCC;
                }
            }
            
            m_r = ray(pos, new_vel);
            *r = m_r;
            return S_GOOD;
        }

        // Otherwise, use Newtonian method (original code)
        double new_t = m_t + epsilon;

        int status = gsl_odeiv2_driver_apply(m_d, &m_t, new_t, m_state);

        if (status != GSL_SUCCESS) {
            fprintf(stderr, "GSL error, return value=%d\n", status);
            return S_ERROR;
        }

        pos = point3(m_state[0], m_state[1], m_state[2]);

        if (!freedom && pos.length() > 200) {
            fprintf(stderr, "freedom denied\n");
            return S_ERROR;
        }

        for (const auto &bh : holes) {
            double dist = (pos - bh.origin).length();
            if (dist <= bh.rs) {
                return S_SUCC;
            }
        }

        vec3 dir(m_state[3], m_state[4], m_state[5]);
        dir = unit_vector(dir);
        m_state[3] = dir.x();
        m_state[4] = dir.y();
        m_state[5] = dir.z();

        m_r = ray(pos, dir);
        *r = m_r;
        return S_GOOD;
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
      disable_bh(disable_bh), mode(Mode::SingleGsl) {
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

// Multi-black hole constructor with simple integration
ray_iterator::ray_iterator(std::vector<black_hole> holes, ray initial_ray,
                           double epsilon, bool prevent_freedom,
                           bool disable_bh)
    : holes(std::move(holes)), epsilon(epsilon), freedom(prevent_freedom),
      disable_bh(disable_bh), mode(Mode::MultiGsl) {
    vec3 initial_dir = unit_vector(initial_ray.direction());
    m_r = ray(initial_ray.origin(), initial_dir);

    double G = 6.6743;
    for (auto &hole : this->holes) {
        hole.rs = 7 * G * hole.mass / 10;
    }

    m_state[0] = m_r.origin().x();
    m_state[1] = m_r.origin().y();
    m_state[2] = m_r.origin().z();
    m_state[3] = initial_dir.x();
    m_state[4] = initial_dir.y();
    m_state[5] = initial_dir.z();

    m_multi_params.holes = &this->holes;
    m_sys = {multi_func, nullptr, 6, &m_multi_params};
    m_d = gsl_odeiv2_driver_alloc_y_new(&m_sys, gsl_odeiv2_step_rkf45, 1e-6,
                                        1e-6, 0.0);
    m_t = 0.0;
}
