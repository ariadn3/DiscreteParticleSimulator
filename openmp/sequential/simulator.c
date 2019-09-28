#include <math.h>

#include "filter.h"
#include "io.h"
#include "kinematics.h"
#include "particle.h"

#define DEBUG_LEVEL 0
#define SLOW_FACTOR 1
#define NO_COLLISION 2

void simulate(void);
double checkWallCollision(double, double, particle_t*);
double checkCollision(double, particle_t*, particle_t*);

int main() {
    simulate();
    return 0;
}

void simulate() {
    params_t* params = read_file(SLOW_FACTOR);

    if (DEBUG_LEVEL > 3) {
        printf("%d %lf %lf %d\n", params->n, params->l, params->r, params->s);
        printf("Simulation printing: %d\n", params->willPrint);
    }

    int n = params->n;
    double l = params->l;
    double r = params->r;
    int s = params->s * SLOW_FACTOR;
    bool willPrint = params->willPrint;
    particle_t** ps = params->particles;

    printAll(n, 0, ps);

    int* numCollisions = (int*) malloc(sizeof(int));
    bool* states = (bool*) malloc(sizeof(bool) * n);
    for (int i = 0; i < n; i++) {
        states[i] = false;
    }
    collision_t** cs = (collision_t**) malloc(sizeof(collision_t*) * n * n / 2); 

    for (int step = 1; step <= s; step++) {
        if (DEBUG_LEVEL > 3) printf("Step %d\n", step);
        *numCollisions = 0;

        // ===== CHECKING AND ADDING COLLISION CANDIDATES =====
        for (int p = 0; p < n; p++) {
            if (DEBUG_LEVEL > 2) printf("Particle %d is p\n", p);
            double wallTime = checkWallCollision(r, l, ps[p]);
            if (DEBUG_LEVEL > 2)
                printf("Particle %d collides with wall at %lf\n", p, wallTime);

            if (wallTime != NO_COLLISION) {
                collision_t* candidate = build_collision(ps[p], NULL, wallTime);
                // #pragma CS
                cs[*numCollisions] = candidate;
                (*numCollisions)++;
                // #end CS
            }

            for (int q = p + 1; q < n; q++) {
                if (DEBUG_LEVEL > 2) printf("Particle %d is q\n", q);
                double time = checkCollision(r, ps[p], ps[q]);

                if (time != NO_COLLISION) {
                    collision_t* candidate = build_collision(ps[p], ps[q], time);
                    // #pragma CS
                    cs[*numCollisions] = candidate;
                    (*numCollisions)++;
                    // #end CS
                }
            }
        }

        if (DEBUG_LEVEL > 1) printf("%d collisions for step %d", *numCollisions, step);

        // ===== FILTER COLLISION CANDIDATES TO VALID COLLISION =====
        filterCollisions(cs, states, numCollisions);
        if (DEBUG_LEVEL > 3) printf("FILTER\n");

        // ===== RESOLVE VALID COLLISIONS =====
        resolveValidCollisions(cs, numCollisions, l, r);

        if (DEBUG_LEVEL > 2) {
            for (int i = 0; i < *numCollisions; i++) {
                printf("%s\n", collision_string(cs[i]));
            }
        }

        updateParticles(ps, n, states);
        if (DEBUG_LEVEL > 3) printf("UPDATE PARTICLES\n");

        // ===== PRINT =====
        if (willPrint || step == s) printAll(n, step, ps);
    }

    if (DEBUG_LEVEL > 3) printf("SIMULATION COMPLETE\n");
}

double checkWallCollision(double r, double l, particle_t* p) {    
    // Collision times with vertical and horizontal walls
    double x_time = NO_COLLISION;
    double y_time = NO_COLLISION;

    double margin = r + EDGE_TOLERANCE;
    // Particle's position after 1 time step
    double x1 = p->x + p->v_x;
    double y1 = p->y + p->v_y;

    // Check if particle would intersect a vertical wall after 1 time step
    // If yes -> compute the time this would happen
    // Also check: if x-velocity is 0 but particle collides with wall
    // -> moving along horizontal wall -> don't try to divide by 0
    if (p->v_x != 0) {
        if (x1 < margin) {
            x_time = (p->x - r) / -(p->v_x); 
        } else if (x1 > l - margin) {
            x_time = (l - r - p->x) / (p->v_x);
        }
    }

    // Check if particle would intersect a horizontal wall after 1 time step
    // If yes -> compute the time this would happen
    // Also check: if y-velocity is 0 but particle collides with wall
    // -> moving along vertical wall -> don't try to divide by 0
    if (p->v_y != 0) {
        if (y1 < margin) {
            y_time = (p->y - r) / -(p->v_y);
        } else if (y1 > l - margin) {
            y_time = (l - r - p->y) / (p->v_y);
        }
    }

    // printf("%lf %lf %lf %lf\n", x_time, y_time, x1, y1);

    // Pick earlier of two times the particle would collide with a wall
    return x_time < y_time ? x_time : y_time;
}

double checkCollision(double r, particle_t* p, particle_t* q) {
    // Difference in X and Y positions and velocities of particles P, Q
    double dX = q->x - p->x;
    double dY = q->y - p->y;
    double dVx = q->v_x - p->v_x;
    double dVy = q->v_y - p->v_y;

    // 0 <= dT <= 1 is the fraction of a time step
    // A, B, C are the coefficients of the (dT)^2, dT and 0-th order terms in
    // the quadratic equation describing distance between particles P, Q at time dT
    double A = dVx * dVx + dVy * dVy;
    double B = 2 * (dX * dVx + dY * dVy);
    double C = dX * dX + dY * dY - 4 * r * r;

    double discriminant = B * B - 4 * A * C;

    if (discriminant <= 0) {
        return NO_COLLISION;
    }
    
    // Distance curve y = d(t) is concave up and intersects y = 2r at two points
    // First intersect (root) is at smaller dT and we only compute this

    // Possible that two particles are currently phasing through (i.e. d(0) < 2r)
    // since only 1 collision was computed per particle -> we ignore any first roots
    // that are dT < 0
    double dT = (-B - sqrt(discriminant)) / 2 / A;

    if (dT >= 0 && dT <= 1) {
        return dT;
    } else {
        return NO_COLLISION;
    }
}

