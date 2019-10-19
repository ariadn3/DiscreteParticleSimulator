#include <math.h>

#include "init.h"
#include "kernels.h"
#include "structs.h"

#define DEBUG_LEVEL 0
#define SLOW_FACTOR 1
#define NO_COLLISION 2

__host__ void simulate();
__host__ void printAll(bool, int, int, particle_t**);
__host__ void resolveValidCollisions(collision_t**, int*, double, double);
__host__ void filterCollisions(collision_t**, bool*, int*);
__host__ int cmpCollision(const void*, const void*);

cudaError_t allocStatus;

// Shared simulation parameters
__constant__ int n, s;
__constant__ double l, r;

// Shared data
__managed__ int* numCollisions;
__managed__ particle_t** ps;
__managed__ bool* states;
__managed__ collision_t** cs;

__host__ void assertMallocSuccess(char* buff) {
    if (allocStatus != cudaSuccess) {
        printf("Failed to dynamically allocate memory for %s\n", buff);
        printf("%s\n", cudaGetErrorString(allocStatus));
        exit(1);
    }
}

__host__ int main(int argc, char** argv) {
    int hostN, hostL, hostR, hostS;
    bool willPrint;

    // Read in N, L, r, S and finally simulation mode
    scanf("%d\n%lf\n%lf\n%d\n", &hostN, &hostL, &hostR, &hostS);
    char* buffer = (char*) malloc(sizeof(char) * 140);
    scanf("%s\n", buffer);

    // Determine if this simulation will run in 'print' or 'perf' mode
    if(strcmp(buffer, "print") == 0) {
        willPrint = true;
    } else if (strcmp(buffer, "perf") == 0) {
        willPrint = false;
    } else {
        printf("Neither 'print' or 'perf' words are present. Exiting...\n");
        exit(1);
    }
    
    // Determine if there is a need to randomise particles
    int i;
    double x, y, v_x, v_y;
    bool isInitialised = false;
    allocStatus = cudaMallocManaged((void**) &ps, hostN * sizeof(particle_t*));
    assertMallocSuccess("particle_t** ps");

    // If initial positions and velocities of particles are provided, read them
    while (fgets(buffer, 140, stdin) != EOF) {
        isInitialised = true;
        sscanf(buffer, "%d %lf %lf %lf %lf", &i, &x, &y, &v_x, &v_y);
        particles[i] = build_particle(i, x, y, v_x / slowFactor, v_y / slowFactor);
    }

    // Otherwise randomise the initial positions and velocities
    if (!isInitialised) randomiseParticles(particles, slowFactor, p->n, p->l, p->r);
    free(buffer);

    // Copy to GPU constant memory
    cudaMemcpyToSymbol(n, &hostN, sizeof(n));
    cudaMemcpyToSymbol(l, &hostL, sizeof(l));
    cudaMemcpyToSymbol(r, &hostR, sizeof(r));
    cudaMemcpyToSymbol(s, &hostS, sizeof(s));

    // Initialise global collision counter
    allocStatus = cudaMallocManaged((void**) &numCollisions, sizeof(int));
    assertMallocSuccess("int* numCollisions");

    // Initialise global particle collision state array
    allocStatus = cudaMallocManaged((void**) &states, hostN * sizeof(bool));
    assertMallocSuccess("bool* states");

    for (int i = 0; i < hostN; i++) {
        states[i] = false;
    }
    
    // Initialise global collisions array - keep up to 8N collision candidates
    allocStatus = cudaMallocManaged((void**) &cs, 8 * hostN * sizeof(collision_t*));
    assertMallocSuccess("collision_t** cs");

    simulate();
    
    return 0;
}

__host__ void simulate() {
    // Unconditionally print the starting state of the simulation
    printAll(false, hostN, 0, ps);
    
    for (int step = 1; step <= s; step++) {
        *numCollisions = 0;

        // ===== CHECKING AND ADDING COLLISION CANDIDATES =====
        for (int p = 0; p < n; p++) {
            double wallTime = checkWallCollision(r, l, ps[p]);
            if (wallTime != NO_COLLISION) {
                collision_t* candidate = build_collision(ps[p], NULL, wallTime);
                // #pragma CS
                cs[*numCollisions] = candidate;
                (*numCollisions)++;
                // #end CS
            }

            for (int q = p + 1; q < n; q++) {
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

        // ===== FILTER COLLISION CANDIDATES TO VALID COLLISION =====
        filterCollisions(cs, states, numCollisions);
        
        // ===== RESOLVE VALID COLLISIONS =====
        resolveValidCollisions(cs, numCollisions, l, r);

        updateParticles(ps, n, states);

        // ===== PRINT SIMULATION DETAILS =====
        if (step == s) printAll(true, n, step, ps);
        else if (willPrint) printAll(false, n, step, ps);
    }
    
    return 0;
}

__host__ void printAll(bool includeCollisions, int n, int step,
        particle_t** particles) {
    // Parallelise this
    for (int i = 0; i < n; i++) {
        char* details;
        if (includeCollisions) {
            details = particle_string_full(particles[i]);
        } else {
            details = particle_string(particles[i]);
        }
        printf("%d %s", step, details);
        free(details);
    }
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

