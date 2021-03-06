#include <omp.h>

#include "kinematics.h"

// Updates particles involved in all valid collisions in collision array
void resolveValidCollisions(collision_t** collisionArray, int* numCollisions,
        float L, float r) {
    collision_t* curCollision;
    // ===== PARALLELISE: COMPUTE NEW POSITIONS AFTER COLLISION =====
    // Each collision involves particles that are updated independently
    #pragma omp parallel private(curCollision)
    {
        int chunk_size = (*numCollisions / omp_get_num_threads()) + 1;
        #pragma omp for schedule(static, chunk_size) 
        for (int i = 0; i < *numCollisions; i++) {
            curCollision = collisionArray[i];
            settleCollision(curCollision, L, r);
            free_collision(curCollision);
        }
    }
    // End parallel for-loop region
}

// Updates particles not involved in any collision
void updateParticles(particle_t** Array, int n, bool* hasCollided) {
    particle_t* curParticle;
    // ===== DO NOT PARALLELISE =====
    // As each particle is updated independently, it it safe to parallelise, but the
    // overhead is much greater than the time for 1 thread to complete the loop
    for (int i = 0; i < n; i++) {
        curParticle = Array[i];
        if (!hasCollided[i]) {
            // Advance particle by its velocity
            curParticle->x += curParticle->v_x;
            curParticle->y += curParticle->v_y;
        } else {
            // Particle had collided -> reset its collision status for next step
            hasCollided[i] = false;
        }
    }
}

// Moves particles involved in a collision to their rightful place after the timestep
void settleCollision(collision_t* curCollision, float L, float r) {
    particle_t* A;
    particle_t* B;
    // Particles A and B (null if wall collision) in this collision
    A = curCollision->p;
    B = curCollision->q;
    float time = curCollision->time;

    // Advance A by the fractional time step dT until collision occurs
    A->x += time * A->v_x;
    A->y += time * A->v_y;

    // If the collision is against the wall, toggle directions
    if (B == NULL) {
        // Add to wall collision counter of A
        A->w_collisions += 1;
        // printf("Step %.14lf: particle %d collided with wall\n", time, A->id);
        if (A->x <= r + EDGE_TOLERANCE || A->x >= L - r - EDGE_TOLERANCE)
            A->v_x *= -1;
        if (A->y <= r + EDGE_TOLERANCE || A->y >= L - r - EDGE_TOLERANCE)
            A->v_y *= -1;
    }
    // If collision is against another particle
    else {
        // Add to particle collision counters of both A and B
        A->p_collisions += 1;
        B->p_collisions += 1;
        // printf("Step %.14lf: particle %d collided with particle %d\n",
        //        time, A->id, B->id);
        // Advance B by dT until collision occurs
        B->x += time * B->v_x;
        B->y += time * B->v_y;

        // Compute distance between A, B
        float distance = sqrt(pow(B->x - A->x, 2) + pow(B->y - A->y, 2));

        // Compute normal and tangent unit vectors along x-, y-axes
        float n_x = (B->x - A->x) / distance;
        float n_y = (B->y - A->y) / distance;
        float t_x = -n_y;
        float t_y = n_x;

        // Compute new normal and tangent unit vectors for particles A, B
        float v_an = n_x * A->v_x + n_y * A->v_y;
        float v_at = t_x * A->v_x + t_y * A->v_y;
        float v_bn = n_x * B->v_x + n_y * B->v_y;
        float v_bt = t_x * B->v_x + t_y * B->v_y;

        // printf("n_x = %.14f, n_y = %.14f\n", n_x, n_y);
        // printf("t_x = %.14f, t_y = %.14f\n", t_x, t_y);
        // printf("v_an = %.14f, v_at = %.14f\n", v_an, v_at);
        // printf("v_bn = %.14f, v_bt = %.14f\n", v_bn, v_bt);

        // printf("Pre-collision velocities: %.14f, %.14f, %.14f, %.14f\n",
        //    A->v_x, A->v_y, B->v_x, B->v_y);

        // Update resultant velocities along x- and y-axes for particles A, B
        A->v_x = v_bn * n_x + v_at * t_x;
        A->v_y = v_bn * n_y + v_at * t_y;
        B->v_x = v_an * n_x + v_bt * t_x;
        B->v_y = v_an * n_y + v_bt * t_y;

        // printf("Post-collision velocities: %.14f, %.14f, %.14f, %.14f\n",
        //    A->v_x, A->v_y, B->v_x, B->v_y);

        // If particle B will collide against the wall, check when it will collide 
        // with the nearest wall and take that time
        float time_bx = 1 - time, time_by = 1 - time;
        if (B->v_x != 0) {
            if (B->x + time_bx * B->v_x < r) time_bx = -(B->x - r) / B->v_x;
            else if (B->x + time_bx * B->v_x > L - r)
                time_bx = (L - r - B->x) / B->v_x;
        }

        if (B->v_y != 0) {
            if (B->y + time_by * B->v_y < r) time_by = -(B->y - r) / B->v_y;
            else if (B->y + time_by * B->v_y > L - r)
                time_by = (L - r - B->y) / B->v_y;
        }

        // If B collides with two walls after colliding with A, take lesser of
        // two times
        float time_b = (time_bx < time_by) ? time_bx : time_by;

        B->x += time_b * B->v_x;
        B->y += time_b * B->v_y;
    }

    // If particle A will collide against the wall, check when it will collide
    // with the nearest wall and take that time
    float time_ax = 1 - time, time_ay = 1 - time;
    
    if (A->v_x != 0) {
        if (A->x + time_ax * A->v_x < r) time_ax = -(A->x - r) / A->v_x;
        else if (A->x + time_ax * A->v_x > L - r) time_ax = (L - r - A->x) / A->v_x;
    }

    if (A->v_y != 0) {
        if (A->y + time_ay * A->v_y < r) time_ay = -(A->y - r)/ A->v_y;
        else if (A->y + time_ay * A->v_y > L - r) time_ay = (L - r - A->y) / A->v_y;
    }

    // If A collides with another wall after colliding, take lesser of two times
    float time_a = (time_ax < time_ay) ? time_ax : time_ay;

    A->x += time_a * A->v_x;
    A->y += time_a * A->v_y;
}

