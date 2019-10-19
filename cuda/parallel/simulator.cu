#include <math.h>

#include "init.h"
#include "kernels.h"
#include "structs.h"

#define DEBUG_LEVEL 0
#define SLOW_FACTOR 1
#define NO_COLLISION 2

__host__ void simulate(void);
__host__ void resolveValidCollisions(collision_t**, int*, double, double);
__host__ void filterCollisions(collision_t**, bool*, int*);
__host__ int cmpCollision(const void*, const void*);

__host__ int main() {
    simulate();
    return 0;
}

__host__ void simulate() {
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

    printAll(false, n, 0, ps);

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

        // ===== PRINT SIMULATION DETAILS =====
        if (step == s) printAll(true, n, step, ps);
        else if (willPrint) printAll(false, n, step, ps);
    }

    if (DEBUG_LEVEL > 3) printf("SIMULATION COMPLETE\n");
}

// Filters the collisions according to the time that it took place
__host__ void filterCollisions(collision_t** collisionArray, bool* hasCollided,
        int* numCollisions) {
    // Quicksort all collision candidates with the comparator function
    qsort(collisionArray, *numCollisions, sizeof(collision_t*), &cmpCollision);

    int saveIndex = 0;
    collision_t* curCollision;
    for (int curIndex = 0; curIndex < *numCollisions; curIndex++) {
        curCollision = collisionArray[curIndex];
        
        // printf("=== Particle %d and %d collided ===\n", curCollision->p->id,
        //         curCollision->q == NULL ? -1 : curCollision->q->id);
        if (hasCollided[curCollision->p->id]
                || (curCollision->q != NULL && hasCollided[curCollision->q->id])) {
            // Particle p has already collided OR particle q has already collided
            // -> discard this colision candidate
            free_collision(curCollision);
        } else {
            // Collision candidate is valid - marked p, q as collided
            hasCollided[curCollision->p->id] = true;

            if (curCollision-> q != NULL) hasCollided[curCollision->q->id] = true;
            // Re-use collision candidates array to store valid collisions
            collisionArray[saveIndex] = collisionArray[curIndex];
            saveIndex++;
        }
    }

    *numCollisions = saveIndex;
}

// Comparator for sorting collisions, earlier time then smaller particle 'p' id
__host__ int cmpCollision(const void* collisionA, const void* collisionB) {
    collision_t* firstCollision = *(collision_t**) collisionA;
    collision_t* secondCollision = *(collision_t**) collisionB;
    
    if (firstCollision->time == secondCollision->time) {
        // If both collisions involve the same first particle
        // Then prioritize wall collision, otherwise prioritize lower 2nd particle ID
        if (firstCollision->p->id == secondCollision->p->id) {
            if (firstCollision->q == NULL) return -1;
            else if (secondCollision->q == NULL) return 1;
            else return (firstCollision->q->id < secondCollision->q->id) ? -1 : 1;
        }
        // If two collisions occur at exactly the same time
        // Then prioritise the one which involves the particle P with lower ID
        return (firstCollision->p->id < secondCollision->p->id) ? -1 : 1;
    } else {
        // Otherwise prioritise the collision occurring at an earlier time
        return (firstCollision->time < secondCollision->time) ? -1 : 1;
    }
}

// Updates particles involved in all valid collisions in collision array
__host__ void resolveValidCollisions(collision_t** collisionArray, int* numCollisions,
        double L, double r) {
    collision_t* curCollision;
    for (int i = 0; i < *numCollisions; i++) {
        curCollision = collisionArray[i];
        settleCollision(curCollision, L, r);
        free_collision(curCollision);
    }
}

