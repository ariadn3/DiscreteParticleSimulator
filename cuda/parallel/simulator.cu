#include <math.h>

#include "init.h"
#include "structs.h"

#define DEBUG_LEVEL 0
#define SLOW_FACTOR 1
#define NO_COLLISION 2
#define EDGE_TOLERANCE 1e-14

__host__ void simulate();
__host__ void printAll(bool, int, int, particle_t*);
__host__ void resolveValidCollisions(collision_t*, int*, double, double);
__host__ void filterCollisions(collision_t*, bool*, int*);
__host__ int cmpCollision(const void*, const void*);

__global__ void checkWallCollision(double, double, particle_t*);
__global__ void checkCollision(double, particle_t*, particle_t*);
__global__ void updateParticles(particle_t**, int, bool*);
__global__ void settleCollision(collision_t*, double, double);

int hostN, hostS;
double hostL, hostR;
bool willPrint;

cudaError_t allocStatus;

// Shared simulation parameters
__constant__ int n, s;
__constant__ double l, r;

// Shared data
__managed__ int numCollisions;
__managed__ particle_t* ps;
__managed__ bool* states;
__managed__ collision_t* cs;

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
    allocStatus = cudaMallocManaged((void**) &ps, hostN * sizeof(particle_t));
    assertMallocSuccess("particle_t** ps");

    // If initial positions and velocities of particles are provided, read them
    while (fgets(buffer, 140, stdin) != EOF) {
        isInitialised = true;
        sscanf(buffer, "%d %lf %lf %lf %lf", &i, &x, &y, &v_x, &v_y);
        particles[i] = build_particle(i, x, y, v_x / slowFactor, v_y / slowFactor);
    }

    // Otherwise randomise the initial positions and velocities
    if (!isInitialised) randomiseParticles(particles, slowFactor, p.n, p.l, p.r);
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
    allocStatus = cudaMallocManaged((void**) &cs, 8 * hostN * sizeof(collision_t));
    assertMallocSuccess("collision_t** cs");

    simulate();
    
    return 0;
}

__host__ void simulate() {
    // Unconditionally print the starting state of the simulation
    printAll(false, hostN, 0, ps);
    
    int pwChunkSize = 32;
    dim3 pwGrid((hostN+pwChunkSize-1)/pwChunkSize);
    dim3 pwThread(pwChunkSize);

    int ppChunkSize = 32
    dim3 ppGrid((hostN+1)/2, (hostN+ppChunkSize-1)/ppChunkSize);
    dim3 ppThread(pwChunkSize);

    int resolveChunkSize = 32;

    int updateChunkSize = 32;
    dim3 updateGrid((hostN+updateChunkSize-1)/updateChunkSize);
    dim3 updateThread(updateChunkSize);

    for (int step = 1; step <= s; step++) {
        numCollisions = 0;

        // ===== CHECKING AND ADDING COLLISION CANDIDATES =====
        checkWallCollision<<<pwGrid, pwThread>>();
        settleCollision<<<ppGrid, ppChunkSize>>>();
        cudaDeviceSynchronize();

        // ===== FILTER COLLISION CANDIDATES TO VALID COLLISION =====
        filterCollisions(cs, states, numCollisions);
        
        // ===== RESOLVE VALID COLLISIONS =====
        dim3 resolveGrid((numCollisions+resolveChunkSize-1)/resolveChunkSize);
        dim3 resolveThread(resolveChunkSize);

        settleCollision<<<resolveGrid, resolveThread>>>();
        updateParticles<<<updateGrid, updateThread>>>();
        cudaDeviceSynchronize();

        // ===== PRINT SIMULATION DETAILS =====
        if (step == s) printAll(true, n, step, ps);
        else if (willPrint) printAll(false, n, step, ps);
    }
    
    return 0;
}

__host__ void printAll(bool includeCollisions, int step, particle_t* particles)
{
    for (int i = 0; i < hostN; i++) {
        char* details;
        if (includeCollisions) {
            details = particle_string_full(&particles[i]);
        } else {
            details = particle_string(&particles[i]);
        }
        printf("%d %s", step, details);
        free(details);
    }
}

// Filters the collisions according to the time that it took place
__host__ void filterCollisions(collision_t* collisionArray, bool* hasCollided) {
    // Quicksort all collision candidates with the comparator function
    qsort(collisionArray, numCollisions, sizeof(collision_t), cmpCollision);

    int saveIndex = 0;
    collision_t curCollision;

    // Walk down collision array and retain valid collisions
    for (int curIndex = 0; curIndex < numCollisions; curIndex++) {
        curCollision = collisionArray[curIndex];
        
        if (hasCollided[curCollision.p.id]
                || (curCollision.q != NULL && hasCollided[curCollision.q.id])) {
            // Particle p has already collided OR particle q has already collided
            // -> discard this colision candidate
            // DO NOTHING (allow this struct to be overwritten later)
        } else {
            // Collision candidate is valid - marked p, q as collided
            hasCollided[curCollision.p.id] = true;

            if (curCollision. q != NULL) hasCollided[curCollision.q.id] = true;
            // Re-use collision candidates array to store valid collisions
            collisionArray[saveIndex] = collisionArray[curIndex];
            saveIndex++;
        }
    }

    numCollisions = saveIndex;
}

// Comparator for sorting collisions, earlier time then smaller particle 'p' id
__host__ int cmpCollision(const void* collisionA, const void* collisionB) {
    collision_t firstCollision = *(collision_t*) collisionA;
    collision_t secondCollision = *(collision_t*) collisionB;
    
    if (firstCollision.time == secondCollision.time) {
        // If both collisions involve the same first particle
        // Then prioritize wall collision, otherwise prioritize lower 2nd particle ID
        if (firstCollision.p.id == secondCollision.p.id) {
            if (firstCollision.q == NULL) return -1;
            else if (secondCollision.q == NULL) return 1;
            else return (firstCollision.q.id < secondCollision.q.id) ? -1 : 1;
        }
        // If two collisions occur at exactly the same time
        // Then prioritise the one which involves the particle P with lower ID
        return (firstCollision.p.id < secondCollision.p.id) ? -1 : 1;
    } else {
        // Otherwise prioritise the collision occurring at an earlier time
        return (firstCollision.time < secondCollision.time) ? -1 : 1;
    }
}


        // for (int p = 0; p < n; p++) {
        //     double wallTime = checkWallCollision(r, l, ps[p]);
        //     if (wallTime != NO_COLLISION) {
        //         collision_t* candidate = build_collision(ps[p], NULL, wallTime);
        //         // #pragma CS
        //         cs[*numCollisions] = candidate;
        //         (*numCollisions)++;
        //         // #end CS
        //     }

        //     for (int q = p + 1; q < n; q++) {
        //         double time = checkCollision(r, ps[p], ps[q]);

        //         if (time != NO_COLLISION) {
        //             collision_t* candidate = build_collision(ps[p], ps[q], time);
        //             // #pragma CS
        //             cs[*numCollisions] = candidate;
        //             (*numCollisions)++;
        //             // #end CS
        //         }
        //     }
        // }

__global__ void checkWallCollision() {
    int index = blockIdx.x*gridDim.x + threadIdx.x;
    if (index >= n)
        return;
    particle_t p = ps[index];

    // Collision times with vertical and horizontal walls
    double x_time = NO_COLLISION;
    double y_time = NO_COLLISION;

    double margin = r + EDGE_TOLERANCE;
    // Particle's position after 1 time step
    double x1 = p.x + p.v_x;
    double y1 = p.y + p.v_y;

    // Check if particle would intersect a vertical wall after 1 time step
    // If yes -> compute the time this would happen
    // Also check: if x-velocity is 0 but particle collides with wall
    // -> moving along horizontal wall -> don't try to divide by 0
    if (p.v_x != 0) {
        if (x1 < margin) {
            x_time = (p.x - r) / -(p.v_x); 
        } else if (x1 > l - margin) {
            x_time = (l - r - p.x) / (p.v_x);
        }
    }

    // Check if particle would intersect a horizontal wall after 1 time step
    // If yes -> compute the time this would happen
    // Also check: if y-velocity is 0 but particle collides with wall
    // -> moving along vertical wall -> don't try to divide by 0
    if (p.v_y != 0) {
        if (y1 < margin) {
            y_time = (p.y - r) / -(p.v_y);
        } else if (y1 > l - margin) {
            y_time = (l - r - p.y) / (p.v_y);
        }
    }

    // printf("%lf %lf %lf %lf\n", x_time, y_time, x1, y1);

    // Pick earlier of two times the particle would collide with a wall
    return x_time < y_time ? x_time : y_time;
}

__global__ void checkCollision() {
    int pIndex = blockIdx.x;
    int qIndex = blockDim.x*blockIdx.y + threadIdx.x;
    particle_t p, q;
    if (qIndex > pIndex) {
        p = ps[pIndex];
        q = ps[qIndex];
    }
    else if (gridDim.x%2 == 0 || pIndex != gridDim.x-1) {
        p = ps[n - 1 - qIndex];
        q = ps[n - 1 - pIndex];
    }
    else
        return;

    // Difference in X and Y positions and velocities of particles P, Q
    double dX = q.x - p.x;
    double dY = q.y - p.y;
    double dVx = q.v_x - p.v_x;
    double dVy = q.v_y - p.v_y;

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

// Moves particles involved in a collision to their rightful place after the timestep
__global__ void settleCollision() {
    int collIndex = blockIdx.x*blockDim.x + threadIdx.x;
    if (collIndex >= numCollisions)
        return;
    collision_t curCollision = cs[collIndex];

    // Particles A and B (null if wall collision) in this collision
    particle_t A = curCollision.p;
    particle_t B = curCollision.q;
    double time = curCollision.time;

    // Advance A by the fractional time step dT until collision occurs
    A.x += time * A.v_x;
    A.y += time * A.v_y;

    // If the collision is against the wall, toggle directions
    if (B == NULL) {
        // Add to wall collision counter of A
        A.w_collisions += 1;
        // printf("Step %.14lf: particle %d collided with wall\n", time, A.id);
        if (A.x <= r + EDGE_TOLERANCE || A.x >= L - r - EDGE_TOLERANCE)
            A.v_x *= -1;
        if (A.y <= r + EDGE_TOLERANCE || A.y >= L - r - EDGE_TOLERANCE)
            A.v_y *= -1;
    }
    // If collision is against another particle
    else {
        // Add to particle collision counters of both A and B
        A.p_collisions += 1;
        B.p_collisions += 1;
        // printf("Step %.14lf: particle %d collided with particle %d\n",
        //        time, A.id, B.id);
        // Advance B by dT until collision occurs
        B.x += time * B.v_x;
        B.y += time * B.v_y;

        // Compute distance between A, B
        double distance = sqrt(pow(B.x - A.x, 2) + pow(B.y - A.y, 2));

        // Compute normal and tangent unit vectors along x-, y-axes
        double n_x = (B.x - A.x) / distance;
        double n_y = (B.y - A.y) / distance;
        double t_x = -n_y;
        double t_y = n_x;

        // Compute new normal and tangent unit vectors for particles A, B
        double v_an = n_x * A.v_x + n_y * A.v_y;
        double v_at = t_x * A.v_x + t_y * A.v_y;
        double v_bn = n_x * B.v_x + n_y * B.v_y;
        double v_bt = t_x * B.v_x + t_y * B.v_y;

        // printf("n_x = %.14f, n_y = %.14f\n", n_x, n_y);
        // printf("t_x = %.14f, t_y = %.14f\n", t_x, t_y);
        // printf("v_an = %.14f, v_at = %.14f\n", v_an, v_at);
        // printf("v_bn = %.14f, v_bt = %.14f\n", v_bn, v_bt);

        // printf("Pre-collision velocities: %.14f, %.14f, %.14f, %.14f\n",
        //    A.v_x, A.v_y, B.v_x, B.v_y);

        // Update resultant velocities along x- and y-axes for particles A, B
        A.v_x = v_bn * n_x + v_at * t_x;
        A.v_y = v_bn * n_y + v_at * t_y;
        B.v_x = v_an * n_x + v_bt * t_x;
        B.v_y = v_an * n_y + v_bt * t_y;

        // printf("Post-collision velocities: %.14f, %.14f, %.14f, %.14f\n",
        //    A.v_x, A.v_y, B.v_x, B.v_y);

        // If particle B will collide against the wall, check when it will collide 
        // with the nearest wall and take that time
        double time_bx = 1 - time, time_by = 1 - time;
        if (B.v_x != 0) {
            if (B.x + time_bx * B.v_x < r) time_bx = -(B.x - r) / B.v_x;
            else if (B.x + time_bx * B.v_x > L - r)
                time_bx = (L - r - B.x) / B.v_x;
        }

        if (B.v_y != 0) {
            if (B.y + time_by * B.v_y < r) time_by = -(B.y - r) / B.v_y;
            else if (B.y + time_by * B.v_y > L - r)
                time_by = (L - r - B.y) / B.v_y;
        }

        // If B collides with two walls after colliding with A, take lesser of
        // two times
        double time_b = (time_bx < time_by) ? time_bx : time_by;

        B.x += time_b * B.v_x;
        B.y += time_b * B.v_y;
    }

    // If particle A will collide against the wall, check when it will collide
    // with the nearest wall and take that time
    double time_ax = 1 - time, time_ay = 1 - time;
    if (A.v_x != 0) {
        if (A.x + time_ax * A.v_x < r) time_ax = -(A.x - r) / A.v_x;
        else if (A.x + time_ax * A.v_x > L - r) time_ax = (L - r - A.x) / A.v_x;
    }

    if (A.v_y != 0) {
        if (A.y + time_ay * A.v_y < r) time_ay = -(A.y - r)/ A.v_y;
        else if (A.y + time_ay * A.v_y > L - r) time_ay = (L - r - A.y) / A.v_y;
    }

    // If A collides with another wall after colliding, take lesser of two times
    double time_a = (time_ax < time_ay) ? time_ax : time_ay;

    A.x += time_a * A.v_x;
    A.y += time_a * A.v_y;
}

// Updates particles not involved in any collision
__global__ void updateParticles() {
    int index = blockIdx.x*gridDim.x + threadIdx.x;
    if (index >= n)
        return;

    particle_t curParticle = ps[index];
    if (!states[i]) {
        // Advance particle by its velocity
        curParticle.x += curParticle.v_x;
        curParticle.y += curParticle.v_y;
    } else {
        // Particle had collided -> reset its collision status for next step
        states[i] = false;
    }
}

