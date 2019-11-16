#include <math.h>

#include "init.h"
#include "structs.h"

#define SLOW_FACTOR 1
#define NO_COLLISION 2
#define EDGE_TOLERANCE 1e-8

#define PERF 0
#define PRINT 1
#define NOT_COLLIDED 0
#define COLLIDED 1

#define MASTER_ID 0

// MASTER routines
void MASTER_main();
void MASTER_init();
void MASTER_printAll(bool, int);
void MASTER_filterCollisions();
int cmpCollision(const void*, const void*);

// SLAVE routines
void SLAVE_main();
void SLAVE_init();

// COMMON routines
void ALL_assertMalloc();

// MASTER computation functions
void MASTER_checkWallCollision();
void MASTER_updateParticles();

// SLAVE Computation functions
void SLAVE_checkCollision();
void SLAVE_settleCollision();

// Simulation parameters common to all MPI processes, master or slave
int n, s;
double l, r;

// Data common to all MPI processes, master or slave
int numCollisions = 0;
particle_t* ps;
int* states;
collision_t* cs;
double minMargin, maxMargin, max;

// MPI information
int numSlaves;
int rank;           // Unique to each MPI process

int main(int argc, char** argv) {

    // ======== BEGIN MPI ========

    int numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WOLRD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int numSlaves = numProcs - 1;

    // Master (rank 0) process executes master routine and initialises simulation
    if (rank == MASTER_ID) {
        MASTER_main();
    } else {
        // Slave processes execute slave routine
        // Waits for data from master to perform computations
        SLAVE_main();
    }

    MPI_Finalize();
    
    // ======== END MPI ========

    return 0;
}

/**
 * Asserts that a malloc is successful, otherwise terminates the entire computation
 * (all MPI processes).
 */
void ALL_assertMalloc(void* ptr, char* buff) {
    if (ptr == NULL) {
        printf("Failed to malloc %s in MPI process with rank %d\n", buff, rank);
        MPI_Abort();
    }
}

/**
 * Master main routine begins by reading in simulation parameters and initialising the
 * simulation state, before beginning simulation.
 */
void MASTER_main() {

    // Determines if the master will print every step of the simulation
    int printMode;

    // Array of integers, each representing the number of data items to retrieve from
    // each MPI process (master itself and all slaves) with MPI_Gatherv
    int* numItems;

    MASTER_init(&printMode, numItems);

    // Unconditionally print the starting state of the simulation
    MASTER_printAll(false, 0);
    
    for (int step = 1; step <= s; step++) {
        numCollisions = 0;

        // ===== CHECKING AND ADDING COLLISION CANDIDATES =====
        // Both kernels have threads that only read particle data, but add collisions
        // to a common array (using atomicAdd to protect each index of this array)
        checkWallCollision<<<pwGrid, pwBlock>>>();
        // You know, we accidentally launched the settleCollision kernel here instead
        // of checkCollision and wondered why particles were colliding 3000 times in
        // 1 step - we wasted 3 hours on this :')
        checkCollision<<<ppGrid, ppBlock>>>();

        cudaDeviceSynchronize();

        // ===== FILTER COLLISION CANDIDATES TO VALID COLLISION =====
        filterCollisions();
        
        // ===== RESOLVE VALID COLLISIONS =====
        dim3 resolveGrid((numCollisions + resolveChunkSize - 1) / resolveChunkSize);
        dim3 resolveBlock(resolveChunkSize);

        // Both kernels have threads that update different particles
        // Safe to run concurrently
        settleCollision<<<resolveGrid, resolveBlock>>>();
        updateParticles<<<updateGrid, updateBlock>>>();
        
        cudaDeviceSynchronize();
        
        // ===== PRINT SIMULATION DETAILS =====
        if (step == s) printAll(true, step);
        else if (willPrint) printAll(false, step);
    }


    // Main body of the simulation
    
    return 0;
}

/**
 * ======== EXECUTED ONLY BY MASTER MPI PROCESS ========
 * Initialises the simulation parameters and the memory required by the master process
 * to store the simulation state.
 */
void MASTER_init(int* mode, int* numItems) {
    
    // Read in N, L, r, S and finally simulation mode
    scanf("%d\n%lf\n%lf\n%d\n", &n, &l, &r, &s);
    char* buffer = (char*) malloc(sizeof(char) * 140);
    scanf("%s\n", buffer);

    // Determine if this simulation will run in 'print' or 'perf' mode
    if(strcmp(buffer, "print") == 0) {
        mode = PRINT;
    } else if (strcmp(buffer, "perf") == 0) {
        mode = PERF;
    } else {
        printf("Neither 'print' or 'perf' words are present. Exiting...\n");
        exit(1);
    }
    
    // Determine if there is a need to randomise particles
    int i;
    double x, y, v_x, v_y;
    bool isInitialised = false;

    // MASTER: allocates its own array of particle structs
    ps = (particle_t*) malloc(n * sizeof(particle_t));
    sprintf(buffer, "particle_t* ps");
    ALL_assertMalloc((void*) ps, buffer);

    // If initial positions and velocities of particles are provided, read them
    while (fgets(buffer, 140, stdin) != NULL) {
        isInitialised = true;
        sscanf(buffer, "%d %lf %lf %lf %lf", &i, &x, &y, &v_x, &v_y);
        ps[i].id = i;
        ps[i].x = x;
        ps[i].y = y;
        ps[i].v_x = v_x / SLOW_FACTOR;
        ps[i].v_y = v_y / SLOW_FACTOR;
        ps[i].w_collisions = 0;
        ps[i].p_collisions = 0;
    }

    // Otherwise randomise the initial positions and velocities
    if (!isInitialised) randomiseParticles(ps, SLOW_FACTOR, n, l, r);

    minMargin = r + EDGE_TOLERANCE;
    maxMargin = l - r - EDGE_TOLERANCE;
    max = l - r;

    // MASTER: allocates its own array of particle collision states
    states = (int*) calloc(n * sizeof(int));
    sprintf(buffer, "int* states");
    ALL_assertMalloc((void*) ps, buffer);

    // MASTER: allocates its own array of collision candidates
    cs = (collision_t*) malloc(n * n * sizeof(collision_t) / 2);
    sprintf(buffer, "collision_t* cs");
    ALL_assertMalloc((void*) cs, buffer);

    // MASTER: allocates array of integers to use with MPI_Gatherv
    numItems = (int*) malloc((numSlaves + 1) * sizeof(int));
    sprintf(buffer, "int* numItems");
    ALL_assertMalloc((void*) numItems, buffer);

}

/**
 * ======== EXECUTED ONLY BY MASTER MPI PROCESS ========
 * Prints the state of the simulation at the current time step.
 */
void MASTER_printAll(bool includeCollisions, int step) {
    char* details = malloc(200 * sizeof(char));
    
    for (int i = 0; i < n; i++) {
        if (includeCollisions) {
            details = particle_string_full(details, &ps[i]);
        } else {
            details = particle_string(details, &ps[i]);
        }
        printf("%d %s", step, details);
    }
    
    free(details);
}

/**
 * ======== EXECUTED ONLY BY MASTER MPI PROCESS ========
 * Filters the array of collision candidates by the time they occurred and the
 * particles involved. Uses the comparator cmpCollision.
 */
void MASTER_filterCollisions() {
    // Quicksort all collision candidates with the comparator function
    qsort(cs, numCollisions, sizeof(collision_t), cmpCollision);

    int saveIndex = 0;
    collision_t curCollision;

    // Walk down collision array and retain valid collisions
    for (int curIndex = 0; curIndex < numCollisions; curIndex++) {
        curCollision = cs[curIndex];
        
        if (states[curCollision.pId]
                || (curCollision.qId != NULL && states[curCollision.qId])) {
            // Particle p has already collided OR particle q has already collided
            // -> discard this colision candidate
            // DO NOTHING (allow this struct to be overwritten later)
        } else {
            // Collision candidate is valid - marked p, q as collided
            states[curCollision.pId] = COLLIDED;

            if (curCollision.qId != NULL) states[curCollision.q->id] = COLLIDED;
            // Re-use collision candidates array to store valid collisions
            cs[saveIndex] = cs[curIndex];
            saveIndex++;
        }
    }

    numCollisions = saveIndex;
}

/**
 * Comparator for sorting collisions, prioritising by earlier time first, then lower
 * particle ID of p.
 */
int cmpCollision(const void* collisionA, const void* collisionB) {
    collision_t firstCollision = *(collision_t*) collisionA;
    collision_t secondCollision = *(collision_t*) collisionB;
   
    if (firstCollision.time == secondCollision.time) {
        // If both collisions involve the same first particle
        // Then prioritize wall collision, otherwise prioritize lower 2nd particle ID
        if (firstCollision.pId == secondCollision.pId) {
            if (firstCollision.qId == NULL) return -1;
            else if (secondCollision.qId == NULL) return 1;
            else return (firstCollision.qId < secondCollision.qId) ? -1 : 1;
        }
        // If two collisions occur at exactly the same time
        // Then prioritise the one which involves the particle P with lower ID
        return (firstCollision.pId < secondCollision.pId) ? -1 : 1;
    } else {
        // Otherwise prioritise the collision occurring at an earlier time
        return (firstCollision.time < secondCollision.time) ? -1 : 1;
    }
}

/**
 * ======== EXECUTED ONLY BY MASTER MPI PROCESS ========
 * Checks if the specified particle collides with a wall during this time step.
 */
void MASTER_checkWallCollision() {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (index >= n)
        return;
   
    particle_t p = ps[index];
    
    // printf("Thread %d checking particle %d\n", index, p.id);

    // Collision times with vertical and horizontal walls
    double x_time = NO_COLLISION;
    double y_time = NO_COLLISION;

    // Particle's position after 1 time step
    double x1 = p.x + p.v_x;
    double y1 = p.y + p.v_y;

    // Check if particle would intersect a vertical wall after 1 time step
    // If yes -> compute the time this would happen
    // Also check: if x-velocity is 0 but particle collides with wall
    // -> moving along horizontal wall -> don't try to divide by 0
    if (p.v_x != 0) {
        if (x1 < minPosMargin) {
            x_time = (p.x - r) / -(p.v_x); 
        } else if (x1 > maxPosMargin) {
            x_time = (maxPos - p.x) / (p.v_x);
        }
    }

    // Check if particle would intersect a horizontal wall after 1 time step
    // If yes -> compute the time this would happen
    // Also check: if y-velocity is 0 but particle collides with wall
    // -> moving along vertical wall -> don't try to divide by 0
    if (p.v_y != 0) {
        if (y1 < minPosMargin) {
            y_time = (p.y - r) / -(p.v_y);
        } else if (y1 > maxPosMargin) {
            y_time = (maxPos - p.y) / (p.v_y);
        }
    }

    // printf("%lf %lf %lf %lf\n", x_time, y_time, x1, y1);

    // Pick earlier of two times the particle would collide with a wall
    double wall_time = x_time < y_time ? x_time : y_time;
    
    if (wall_time != NO_COLLISION) {
        // atomicAdd returns the previous value of that address - we use this as a
        // ticket for this thread to write a collision to that specific index
        // Implicitly serves as a critical section
        int i = atomicAdd(&numCollisions, 1);
        // printf("CS%d: added by thread %d\n", i, index); 

        cs[i].p = &ps[p.id];
        cs[i].q = NULL;
        cs[i].time = wall_time;
    }
}

/**
 * ======== EXECUTED ONLY BY MASTER MPI PROCESS ========
 * Updates the positions of all particles that did not collide during this time step.
 */
void MASTER_updateParticles() {
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index >= n)
        return;

    particle_t* curParticle = &ps[index];
    
    if (!states[index]) {
        // Advance particle by its velocity
        curParticle->x += curParticle->v_x;
        curParticle->y += curParticle->v_y;
    } else {
        // Particle had collided -> reset its collision status for next step
        states[index] = false;
    }
}

/**
 * ======== EXECUTED ONLY BY SLAVE MPI PROCESSES  ========
 * Initialises memory required by slave processes to store data required for their
 * computations.
 */
void SLAVE_init() {
}


void SLAVE_checkCollision() {

    // printf("%d %d % %d\n", gridDim.x, gridDim.y, blockDim.x, threadIdx.x);
    
    int pIndex = blockIdx.x;
    int qIndex = blockDim.x * blockIdx.y + threadIdx.x;

    particle_t p, q;
    
    // Ignore excess threads beyond the row of computation
    if (qIndex >= n) return;

    // Compute upper half of triangle
    if (qIndex > pIndex) {
        p = ps[pIndex];
        q = ps[qIndex];
    } else if (n % 2 == 1 || pIndex != gridDim.x - 1) {
        // N is odd -> reflected lower half folds correctly to form rows of length N
        // with no excess
        // Compute reflected lower half of triangle folded to form a row
        p = ps[n - 2 - pIndex];
        q = ps[n - 1 - qIndex];
    } else {
        // Catch case where n is even -> odd number of rows of computation when folded
        // Ignore excess threads beyond the middle column
        return;
    }

    // printf("Checking array computation (%d, %d)\n", p.id, q.id);

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
        return;
    }
    
    // Distance curve y = d(t) is concave up and intersects y = 2r at two points
    // First intersect (root) is at smaller dT and we only compute this

    // Possible that two particles are currently phasing through (i.e. d(0) < 2r)
    // since only 1 collision was computed per particle -> we ignore any first roots
    // that are dT < 0
    double dT = (-B - sqrt(discriminant)) / 2 / A;

    // Add a collision candidate if P, Q would collide during this time step
    if (dT >= 0 && dT <= 1) {
        // atomicAdd returns the previous value of that address - we use this as a
        // ticket for this thread to write a collision to that specific index
        // Implicitly serves as a critical section
        int i = atomicAdd(&numCollisions, 1);
        // printf("CS%d: p-p added by block %d thread %d\n", i, pIndex, qIndex); 
        
        cs[i].p = &ps[p.id];
        cs[i].q = &ps[q.id];
        cs[i].time = dT;
    }
}

// Moves particles involved in a collision to their rightful place after the timestep
void SLAVE_settleCollision() {
    int collIndex = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (collIndex >= numCollisions)
        return;
    
    collision_t curCollision = cs[collIndex];

    // Particles A and B (null if wall collision) in this collision
    particle_t* A = curCollision.p;
    particle_t* B = curCollision.q;
    double time = curCollision.time;

    // Advance A by the fractional time step dT until collision occurs
    A->x += time * A->v_x;
    A->y += time * A->v_y;

    // If the collision is against the wall, toggle directions
    if (B == NULL) {
        // Add to wall collision counter of A
        A->w_collisions += 1;
        // printf("Step %.14lf: particle %d collided with wall\n", time, A->id);
        if (A->x <= minPosMargin || A->x >= maxPosMargin)
            A->v_x = A->v_x * -1;
        if (A->y <= minPosMargin || A->y >= maxPosMargin)
            A->v_y = A->v_y * -1;
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
        double distance = sqrt(pow(B->x - A->x, 2) + pow(B->y - A->y, 2));

        // Compute normal and tangent unit vectors along x-, y-axes
        double n_x = (B->x - A->x) / distance;
        double n_y = (B->y - A->y) / distance;
        double t_x = -n_y;
        double t_y = n_x;

        // Compute new normal and tangent unit vectors for particles A, B
        double v_an = n_x * A->v_x + n_y * A->v_y;
        double v_at = t_x * A->v_x + t_y * A->v_y;
        double v_bn = n_x * B->v_x + n_y * B->v_y;
        double v_bt = t_x * B->v_x + t_y * B->v_y;

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
        double time_bx = 1 - time, time_by = 1 - time;
        if (B->v_x != 0) {
            if (B->x + time_bx * B->v_x < r) time_bx = -(B->x - r) / B->v_x;
            else if (B->x + time_bx * B->v_x > maxPos)
                time_bx = (maxPos - B->x) / B->v_x;
        }

        if (B->v_y != 0) {
            if (B->y + time_by * B->v_y < r) time_by = -(B->y - r) / B->v_y;
            else if (B->y + time_by * B->v_y > maxPos)
                time_by = (maxPos - B->y) / B->v_y;
        }

        // If B collides with two walls after colliding with A, take lesser of
        // two times
        double time_b = (time_bx < time_by) ? time_bx : time_by;

        B->x += time_b * B->v_x;
        B->y += time_b * B->v_y;
    }

    // If particle A will collide against the wall, check when it will collide
    // with the nearest wall and take that time
    double time_ax = 1 - time;
    double time_ay = 1 - time;
    
    if (A->v_x != 0) {
        if (A->x + time_ax * A->v_x < r) time_ax = -(A->x - r) / A->v_x;
        else if (A->x + time_ax * A->v_x > maxPos) time_ax = (maxPos - A->x) / A->v_x;
    }

    if (A->v_y != 0) {
        if (A->y + time_ay * A->v_y < r) time_ay = -(A->y - r)/ A->v_y;
        else if (A->y + time_ay * A->v_y > maxPos) time_ay = (maxPos - A->y) / A->v_y;
    }

    // If A collides with another wall after colliding, take lesser of two times
    double time_a = (time_ax < time_ay) ? time_ax : time_ay;

    A->x += time_a * A->v_x;
    A->y += time_a * A->v_y;
}

