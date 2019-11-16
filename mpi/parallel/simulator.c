#include <math.h>
#include <mpi.h>

#include "init.h"
#include "structs.h"

#define SLOW_FACTOR 1
#define NO_COLLISION 2
#define EDGE_TOLERANCE 1e-8

#define FALSE 0
#define TRUE 1
#define PERF 0
#define PRINT 1
#define NOT_COLLIDED 0
#define COLLIDED 1

#define MASTER_ID 0
#define NUM_SLAVES numProcs - 1

// MASTER routines
void MASTER_main();
void MASTER_init(int*, int*, int*, int*);
void MASTER_printAll(int, int);
void MASTER_buildDisplacement(int*, int*);
void MASTER_divideCollisions(int*, int*);
void MASTER_filterCollisions(int*);
void MASTER_mergeResolvedParticles(particle_t*);
int cmpCollision(const void*, const void*);

// SLAVE routines
void SLAVE_main();
void SLAVE_init();

// COMMON routines
void ALL_assertMalloc();

// MASTER computation functions
void MASTER_checkWallCollisions();
void MASTER_updateParticles(int*);

// SLAVE Computation functions
void SLAVE_checkCollisions();
void SLAVE_settleCollisions(int*);

// Simulation parameters common to all MPI processes, master or slave
int n, s;
double l, r;

// Data common to all MPI processes, master or slave
int numCollisions = 0;
particle_t* ps;
particle_t* pBuffer;
collision_t* cs;
double minMargin, maxMargin, max;

// MPI information
int zero = 0;
int numProcs;
int rank;           // Unique to each MPI process
MPI_Datatype MPI_PARTICLE;
MPI_Datatype MPI_COLLISION;

int main(int argc, char** argv) {

    // ======== BEGIN MPI ========

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create the MPI struct representations for a particle and collision
    int particleBlockLengths[3] = {1, 4, 2};
    MPI_Aint particleBlockOffsets[3] = {0, 1, 5};
    MPI_Datatype particleBlockTypes[3] = {MPI_INT, MPI_DOUBLE, MPI_INT};
    MPI_Type_create_struct(3, particleBlockLengths,
            particleBlockOffsets, particleBlockTypes, &MPI_PARTICLE);

    int collisionBlockLengths[2] = {2, 1};
    MPI_Aint collisionBlockOffsets[2] = {0, 2};
    MPI_Datatype collisionBlockTypes[2] = {MPI_INT, MPI_DOUBLE};
    MPI_Type_create_struct(1, collisionBlockLengths,
            collisionBlockOffsets, collisionBlockTypes, &MPI_COLLISION);

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
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

/**
 * Master main routine begins by reading in simulation parameters and initialising the
 * simulation state, before beginning simulation.
 */
void MASTER_main() {
    // Determines if the master will print every step of the simulation
    int printMode;

    // Array storing collision states of particles for each time step
    int* states;

    // Array of integers, each representing the number of data items to retrieve from
    // each MPI process (master itself and all slaves) with MPI_Gatherv
    int* numItems;

    // Array of integers, each representing the displacement from the start of a buffer
    // at which to place the received data from the process with ith rank
    int* displc;

    MASTER_init(&printMode, states, numItems, displc);

    MASTER_printAll(FALSE, 0);

    // Synchronise to ensure all slaves have completed initialisation
    MPI_Barrier(MPI_COMM_WORLD);

    // ======== MAIN BODY OF THE SIMULATION ========
    for (int step = 1; step <= s; step++) {

        // ===== BROADCAST UPDATED PARTICLE DATA TO ALL SLAVES =====
        // Slaves will divide all the work of computing particle-particle collisions
        // among themselves; master will concurrently compute particle-wall collisions
        MPI_Bcast(ps, n, MPI_FLOAT, MASTER_ID, MPI_COMM_WORLD);

        // Master concurrently computes particle-wall collisions
        numCollisions = 0;
        MASTER_checkWallCollisions();

        // Ensure all slaves have finished computing before attempting to gather
        MPI_Barrier(MPI_COMM_WORLD);

        // Collect the collision candidates from all MPI processes (master included)
        // back into the master's array of collision candidates
        MPI_Gather(&numCollisions, 1, MPI_INT,
                numItems, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
        MASTER_buildDisplacement(numItems, displc);
        MPI_Gatherv(cs, numCollisions, MPI_COLLISION,
                cs, numItems, displc, MPI_COLLISION, MASTER_ID, MPI_COMM_WORLD);

        // ===== FILTER COLLISION CANDIDATES TO VALID COLLISION =====
        MASTER_filterCollisions(states);

        // ===== SCATTER VALID COLLISIONS AMONGST SLAVES ====
        // Scatter valid collisions amongst slave processes; master will concurrently
        // update positions of all particles that did not collide in this step
        // NOTE: this requires use of MPI_IN_PLACE to prevent master from being a
        // recipient of the scatter
        MASTER_divideCollisions(numItems, displc);
        MPI_Scatter(&numItems, 1, MPI_INT,
                &numCollisions, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
        MASTER_divideCollisions(numItems, displc);
        MPI_Scatterv(cs, numItems, displc, MPI_COLLISION,
                cs, 0, MPI_COLLISION, MASTER_ID, MPI_COMM_WORLD);

        MASTER_updateParticles(states);

        // Ensure all slaves have finished resolving collisions before gathering
        MPI_Barrier(MPI_COMM_WORLD);

        // Collect all updated particles back into the master's buffer of particles
        // Note the &zero is used here since the master has 0 updated particles to
        // send to itself
        MPI_Gather(&zero, 1, MPI_INT,
                numItems, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
        MASTER_buildDisplacement(numItems, displc);
        MPI_Gatherv(pBuffer, 0, MPI_PARTICLE,
                pBuffer, numItems, displc, MPI_PARTICLE, MASTER_ID, MPI_COMM_WORLD);

        MASTER_mergeResolvedParticles(pBuffer);

        // ===== PRINT SIMULATION DETAILS =====
        if (step == s) MASTER_printAll(TRUE, step);
        else if (printMode == PRINT) MASTER_printAll(FALSE, step);
    }
}

/**
 * ======== EXECUTED ONLY BY MASTER MPI PROCESS ========
 * Given an array describing number of items to receive from the ith process,
 * builds the array of displacements from the start of the buffer to place the items
 * received from the ith process.
 */
void MASTER_buildDisplacement(int* items, int* displc) {
    displc[0] = 0;
    for (int i = 1; i < numProcs; i++) {
        displc[i] = displc[i - 1] + items[i - 1];
    }
}

/**
 * ======== EXECUTED ONLY BY MASTER MPI PROCESS ========
 * Given a number of valid collisions, constructs the array describing the number of
 * collisions to assign to each slave process, and the displacement of the contiguous
 * block of collisions in the collisions array for each process.
 */
void MASTER_divideCollisions(int* items, int* displc) {
    // Master will not receive any collision structs to resolve
    items[0] = 0;
    displc[0] = 0;

    // Divide evenly and take the quotient
    int chunkSize = numCollisions / NUM_SLAVES;
    // Distribute the remainder R collisions amongst the first R slaves
    int remainder = numCollisions - (chunkSize * NUM_SLAVES);

    for (int i = 1; i < numProcs; i++) {
        items[i] = chunkSize;
        if (i <= remainder) items[i]++;
        displc[i] = displc[i - 1] + items[i - 1];
    }
}

/**
 * ======== EXECUTED ONLY BY MASTER MPI PROCESS ========
 * Initialises the simulation parameters and the memory required by the master process
 * to store the simulation state.
 */
void MASTER_init(int* mode, int* states, int* numItems, int* displc) {
    // Read in N, L, r, S and finally simulation mode
    scanf("%d\n%lf\n%lf\n%d\n", &n, &l, &r, &s);
    char* buffer = (char*) malloc(sizeof(char) * 140);
    scanf("%s\n", buffer);

    // Determine if this simulation will run in 'print' or 'perf' mode
    if(strcmp(buffer, "print") == 0) {
        *mode = PRINT;
    } else if (strcmp(buffer, "perf") == 0) {
        *mode = PERF;
    } else {
        printf("Neither 'print' or 'perf' words are present. Exiting...\n");
        exit(1);
    }

    // Determine if there is a need to randomise particles
    int i;
    double x, y, v_x, v_y;
    int isInitialised = 0;

    // MASTER: allocates its own array of particle structs
    ps = (particle_t*) malloc(n * sizeof(particle_t));
    sprintf(buffer, "particle_t* ps");
    ALL_assertMalloc((void*) ps, buffer);

    // If initial positions and velocities of particles are provided, read them
    while (fgets(buffer, 140, stdin) != NULL) {
        isInitialised = 1;
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
    if (isInitialised > 0) randomiseParticles(ps, SLOW_FACTOR, n, l, r);

    minMargin = r + EDGE_TOLERANCE;
    maxMargin = l - r - EDGE_TOLERANCE;
    max = l - r;

    // MASTER: allocates its own BUFFER array of particle structs
    pBuffer = (particle_t*) malloc(n * sizeof(particle_t));
    sprintf(buffer, "particle_t* pBuffer");
    ALL_assertMalloc((void*) pBuffer, buffer);

    // MASTER: allocates its own array of particle collision states
    states = (int*) calloc(n, sizeof(int));
    sprintf(buffer, "int* states");
    ALL_assertMalloc((void*) ps, buffer);

    // MASTER: allocates its own array of collision candidates
    cs = (collision_t*) malloc(n * n * sizeof(collision_t) / 2);
    sprintf(buffer, "collision_t* cs");
    ALL_assertMalloc((void*) cs, buffer);

    // MASTER: allocates array of integers to use with MPI_Gatherv
    numItems = (int*) malloc((int) numProcs * sizeof(int));
    sprintf(buffer, "int* numItems");
    ALL_assertMalloc((void*) numItems, buffer);

    // MASTER: allocates displacement array of itnegers to use with MPI_Gatherv
    displc = (int*) malloc((int) numProcs * sizeof(int));
    sprintf(buffer, "int* displc");
    ALL_assertMalloc((void*) displc, buffer);

    free(buffer);
}

/**
 * ======== EXECUTED ONLY BY MASTER MPI PROCESS ========
 * Prints the state of the simulation at the current time step.
 */
void MASTER_printAll(int includeCollisions, int step) {
    char* details = malloc(200 * sizeof(char));

    for (int i = 0; i < n; i++) {
        if (includeCollisions == TRUE) {
            particle_string_full(details, &ps[i]);
        } else {
            particle_string(details, &ps[i]);
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
void MASTER_filterCollisions(int* states) {
    // Quicksort all collision candidates with the comparator function
    qsort(cs, numCollisions, sizeof(collision_t), cmpCollision);

    int saveIndex = 0;
    collision_t curCollision;

    // Walk down collision array and retain valid collisions
    for (int curIndex = 0; curIndex < numCollisions; curIndex++) {
        curCollision = cs[curIndex];

        if (states[curCollision.pId]
                || (curCollision.qId != WALL && states[curCollision.qId])) {
            // Particle p has already collided OR particle q has already collided
            // -> discard this colision candidate
            // DO NOTHING (allow this struct to be overwritten later)
        } else {
            // Collision candidate is valid - marked p, q as collided
            states[curCollision.pId] = COLLIDED;

            if (curCollision.qId != WALL) states[curCollision.qId] = COLLIDED;
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
            if (firstCollision.qId == WALL) return -1;
            else if (secondCollision.qId == WALL) return 1;
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
void MASTER_checkWallCollisions() {
    for (int i = 0; i < n; i++) {
        particle_t p = ps[i];

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
            if (x1 < minMargin) {
                x_time = (p.x - r) / -(p.v_x); 
            } else if (x1 > maxMargin) {
                x_time = (max - p.x) / (p.v_x);
            }
        }

        // Check if particle would intersect a horizontal wall after 1 time step
        // If yes -> compute the time this would happen
        // Also check: if y-velocity is 0 but particle collides with wall
        // -> moving along vertical wall -> don't try to divide by 0
        if (p.v_y != 0) {
            if (y1 < minMargin) {
                y_time = (p.y - r) / -(p.v_y);
            } else if (y1 > maxMargin) {
                y_time = (max - p.y) / (p.v_y);
            }
        }

        // Pick earlier of two times the particle would collide with a wall
        double wall_time = x_time < y_time ? x_time : y_time;

        if (wall_time != NO_COLLISION) {
            // Add a new collision candidate to the master's collisions array
            cs[numCollisions].pId = i;
            cs[numCollisions].qId = WALL;
            cs[numCollisions].time = wall_time;
            numCollisions++;
        }
    }
}

/**
 * ======== EXECUTED ONLY BY MASTER MPI PROCESS ========
 * Updates the positions of all particles that did not collide during this time step.
 */
void MASTER_updateParticles(int* states) {
    for (int i = 0; i < 0; i++) {
        if (states[i] != COLLIDED) {
            // Advance particle by its velocity
            ps[i].x += ps[i].v_x;
            ps[i].y += ps[i].v_y;
        } else {
            // Particle had collided -> reset its collision status for next time step
            states[i] = NOT_COLLIDED;
        }
    }
}

/**
 * ======== EXECUTED ONLY BY MASTER MPI PROCESS ========
 * Merges the updated state of particles involved in collisions computed by slave
 * processes into the master's particle array.
 */
void MASTER_mergeResolvedParticles(particle_t* newParticles) {
    // The number of updated particles should be 2x the number of valid collisions
    for (int i = 0; i < 2 * numCollisions; i++) {
        ps[newParticles[i].id] = newParticles[i];
    }
}

/**
 * Slave main routine begins by having tea
 */
void SLAVE_main() {

    int numParticlesUpdated;

    SLAVE_init();

    for (int step = 1; step <= s; step++) {
        // Update its copy of particles from master's broadcast
        MPI_Bcast(ps, n, MPI_PARTICLE, MASTER_ID, MPI_COMM_WORLD);

        // Check and add collision candidates in its own area based on rank
        SLAVE_checkCollisions();

        // Synchronise all MPI processes prior to gather
        MPI_Barrier(MPI_COMM_WORLD);

        // Send number of collision candidates to master
        MPI_Gather(&numCollisions, 1, MPI_INT,
                NULL, 0, MPI_INT, MASTER_ID, MPI_COMM_WORLD);

        // Gatherv the collision candidates back to master
        MPI_Gatherv(cs, numCollisions, MPI_COLLISION,
                NULL, NULL, NULL, MPI_COLLISION, MASTER_ID, MPI_COMM_WORLD);

        // Obtain number of collisions to settle
        MPI_Scatter(NULL, 0, MPI_INT,
                &numCollisions, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);

        // Obtain list of collisions to settle
        MPI_Scatterv(NULL, NULL, NULL, MPI_COLLISION,
                cs, numCollisions, MPI_COLLISION, MASTER_ID, MPI_COMM_WORLD);

        // Settle them one at a time
        SLAVE_settleCollisions(&numParticlesUpdated);

        // Synchronise all MPI processes prior to gather
        MPI_Barrier(MPI_COMM_WORLD);

        // Send number of updated particles back to master
        MPI_Gather(&numParticlesUpdated, 1, MPI_INT,
                NULL, 0, MPI_INT, MASTER_ID, MPI_COMM_WORLD);

        // Gatherv the particles back to master
        MPI_Gatherv(pBuffer, numParticlesUpdated, MPI_PARTICLE,
                NULL, NULL, NULL, MPI_PARTICLE, MASTER_ID, MPI_COMM_WORLD);
    }
}

/**
 * ======== EXECUTED ONLY BY SLAVE MPI PROCESSES  ========
 * Initialises memory required by slave processes to store data required for their
 * computations.
 */
void SLAVE_init() {
    char* buffer = (char*) malloc(sizeof(char) * 140);
    
    // SLAVE: allocates its own array of particle structs
    ps = (particle_t*) malloc(n * sizeof(particle_t));
    sprintf(buffer, "particle_t* ps");
    ALL_assertMalloc((void*) ps, buffer);

    // SLAVE: allocates its own BUFFER array of particle structs
    pBuffer = (particle_t*) malloc(n * sizeof(particle_t));
    sprintf(buffer, "particle_t* pBuffer");
    ALL_assertMalloc((void*) pBuffer, buffer);

    // SLAVE: allocates its own array of collision candidates (smaller than master)
    cs = (collision_t*) malloc(n * n * sizeof(collision_t) / 2 / NUM_SLAVES);
    sprintf(buffer, "collision_t* cs");
    ALL_assertMalloc((void*) cs, buffer);

    free(buffer);
}

void SLAVE_checkCollisions() {
    int chunkSize = (n + NUM_SLAVES - 1)/NUM_SLAVES;
    numCollisions = 0;

    // printf("%d %d % %d\n", gridDim.x, gridDim.y, blockDim.x, threadIdx.x);

    for (int i = (rank - 1) * chunkSize; i < rank * chunkSize; i++) {
        for (int j = 0; j < n - 1; j++) {

            if (i > (n + 1) / 2 || (n % 2 == 0 && j < n/2) ) continue;

            int pIndex = i;
            int qIndex = j;

            particle_t p, q;

            // Compute upper half of triangle
            if (qIndex > pIndex) {
                p = ps[pIndex];
                q = ps[qIndex];
            } else {
                // N is odd -> reflected lower half folds correctly to form rows 
                // of length N with no excess
                // Compute reflected lower half of triangle folded to form a row
                p = ps[n - 2 - pIndex];
                q = ps[n - 1 - qIndex];
            }

            // Difference in X and Y positions and velocities of particles P, Q
            double dX = q.x - p.x;
            double dY = q.y - p.y;
            double dVx = q.v_x - p.v_x;
            double dVy = q.v_y - p.v_y;

            // 0 <= dT <= 1 is the fraction of a time step
            // A, B, C are the coefficients of the (dT)^2, dT and 0-th order terms in
            // the quadratic equation describing distance between particles P, Q at
            // time dT
            double A = dVx * dVx + dVy * dVy;
            double B = 2 * (dX * dVx + dY * dVy);
            double C = dX * dX + dY * dY - 4 * r * r;

            double discriminant = B * B - 4 * A * C;

            if (discriminant <= 0) {
                return;
            }

            // Distance curve y = d(t) is concave up and intersects y = 2r at 2 points
            // First intersect (root) is at smaller dT and we only compute this
            // Possible 2 particles are currently phasing through (i.e. d(0) < 2r)
            // since only 1 collision was computed per particle -> we ignore any first
            // roots that are dT < 0
            double dT = (-B - sqrt(discriminant)) / 2 / A;

            // Add a collision candidate if P, Q would collide during this time step
            if (dT >= 0 && dT <= 1) {
                cs[numCollisions].pId = p.id;
                cs[numCollisions].qId = q.id;
                cs[numCollisions].time = dT;
                numCollisions++;
            }
        }
    }
}

// Moves particles involved in a collision to their rightful place after the timestep
void SLAVE_settleCollisions(int* numParticlesUpdated) {
    *numParticlesUpdated = 0;
    for (int collIndex = 0; collIndex < numCollisions; collIndex++) {

        collision_t curCollision = cs[collIndex];

        // Particles A and B (null if wall collision) in this collision
        particle_t* A = &ps[curCollision.pId];
        particle_t* B = NULL;
        if (curCollision.qId != WALL) B = &ps[curCollision.qId];
        double time = curCollision.time;

        // Advance A by the fractional time step dT until collision occurs
        A->x += time * A->v_x;
        A->y += time * A->v_y;

        // If the collision is against the wall, toggle directions
        if (B == NULL) {
            // Add to wall collision counter of A
            A->w_collisions += 1;
            // printf("Step %.14lf: particle %d collided with wall\n", time, A->id);
            if (A->x <= minMargin || A->x >= maxMargin)
                A->v_x = A->v_x * -1;
            if (A->y <= minMargin || A->y >= maxMargin)
                A->v_y = A->v_y * -1;
        }
        // If collision is against another particle
        else {
            // Add to particle collision counters of both A and B
            A->p_collisions += 1;
            B->p_collisions += 1;
            
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

            // Update resultant velocities along x- and y-axes for particles A, B
            A->v_x = v_bn * n_x + v_at * t_x;
            A->v_y = v_bn * n_y + v_at * t_y;
            B->v_x = v_an * n_x + v_bt * t_x;
            B->v_y = v_an * n_y + v_bt * t_y;

            // If particle B will collide against the wall, check when it will collide 
            // with the nearest wall and take that time
            double time_bx = 1 - time, time_by = 1 - time;
            if (B->v_x != 0) {
                if (B->x + time_bx * B->v_x < r) time_bx = -(B->x - r) / B->v_x;
                else if (B->x + time_bx * B->v_x > max)
                    time_bx = (max - B->x) / B->v_x;
            }

            if (B->v_y != 0) {
                if (B->y + time_by * B->v_y < r) time_by = -(B->y - r) / B->v_y;
                else if (B->y + time_by * B->v_y > max)
                    time_by = (max - B->y) / B->v_y;
            }

            // If B collides with two walls after colliding with A, take lesser of
            // two times
            double time_b = (time_bx < time_by) ? time_bx : time_by;

            B->x += time_b * B->v_x;
            B->y += time_b * B->v_y;

            pBuffer[numCollisions] = *B;
            numCollisions++;
        }

        // If particle A will collide against the wall, check when it will collide
        // with the nearest wall and take that time
        double time_ax = 1 - time;
        double time_ay = 1 - time;

        if (A->v_x != 0) {
            if (A->x + time_ax * A->v_x < r) time_ax = -(A->x - r) / A->v_x;
            else if (A->x + time_ax * A->v_x > max) time_ax = (max - A->x) / A->v_x;
        }

        if (A->v_y != 0) {
            if (A->y + time_ay * A->v_y < r) time_ay = -(A->y - r)/ A->v_y;
            else if (A->y + time_ay * A->v_y > max) time_ay = (max - A->y) / A->v_y;
        }

        // If A collides with another wall after colliding, take lesser of two times
        double time_a = (time_ax < time_ay) ? time_ax : time_ay;

        A->x += time_a * A->v_x;
        A->y += time_a * A->v_y;

        pBuffer[numCollisions] = *A;
        numCollisions++;
    }
}

