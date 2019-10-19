#include "init.h"

// Read in inputs from file and return initial simulation parameters as params_t
__host__ params_t* read_file(int slowFactor) {
    params_t* p = params();

    // Read in N, L, r, S and finally simulation mode
    scanf("%d\n%lf\n%lf\n%d\n", &(p->n), &(p->l), &(p->r), &(p->s));
    char* buffer = (char*) malloc(sizeof(char) * 140);
    scanf("%s\n", buffer);

    // Determine if this simulation will run in 'print' or 'perf' mode
    if(strcmp(buffer, "print") == 0) {
        p->willPrint = true;
    } else if (strcmp(buffer, "perf") == 0) {
        p->willPrint = false;
    } else {
        printf("Neither 'print' or 'perf' words are present. Exiting...\n");
        exit(1);
    }

    int i;
    double x, y, v_x, v_y;
    bool isInitialised = false;
    particle_t** particles = (particle_t**) malloc(p->n * sizeof(particle_t));

    // If initial positions and velocities of particles are provided, read them
    while (fgets(buffer, 140, stdin) != NULL) {
        isInitialised = true;
        sscanf(buffer, "%d %lf %lf %lf %lf", &i, &x, &y, &v_x, &v_y);
        particles[i] = build_particle(i, x, y, v_x / slowFactor, v_y / slowFactor);
    }

    // Otherwise randomise the initial positions and velocities
    if (!isInitialised) {
        randomiseParticles(particles, slowFactor, p->n, p->l, p->r);
    }

    p->particles = particles;
    return p;
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

// Randomly generates an array of particles
__host__ void randomiseParticles(particle_t** particleArray, int slowFactor, int n,
        double L, double r) {
    double* posArray = generatePosition(n, L, r);
    double* veloArray = generateVelocity(slowFactor, n, L, r);
    for (int i = 0; i < n; i++) {
        particleArray[i] = build_particle(i, posArray[2 * i], posArray[2 * i + 1],
                veloArray[2 * i], veloArray[2 * i + 1]);
    }
    free(posArray);
    free(veloArray);
}

// Generates an array of position values
__host__ double* generatePosition(int n, double L, double r) {
    static double* posArray;
    posArray = (double*) malloc(n * 2 * sizeof(double));
    srand(SEED);
    
    // Checks pre-conditions (read our report)
    if (L < 2 * r) {
        printf("Assumption violated: L < (2 * r)\n");
        exit(1);
    } else if (n * r * r > L * L) {
        printf("Assumption violated: (n * r * r) > (L * L)\n");
        exit(1);
    }
    
    // Bounds for generating positions
    double leftLimit = r, rightLimit = L - r, lenDiff = rightLimit - leftLimit;
    for (int i = 0; i < n; i++) {
        while (true) {
            posArray[2 * i] = leftLimit + lenDiff * (rand() / (double)RAND_MAX);
            posArray[2 * i +1] = leftLimit + lenDiff * (rand() / (double)RAND_MAX);
            bool flag = true;
            for (int j = 0; j < i; j++) {
                if (2 * r > sqrt(pow(posArray[2 * i] - posArray[2 * j], 2)
                            + pow(posArray[2 * i + 1] - posArray[2 * j + 1], 2))) {
                    flag = false;
                    break;
                }
            }
            if (flag) break;
        }
    }

    return posArray;
}

// Generates an array of velocity values
__host__ double* generateVelocity(int slowFactor, int n, double L, double r) {
    static double* veloArray;
    veloArray = (double*) malloc(n * 2 * sizeof(double));
    
    // Bounds for generating velocities
    double veloLeftLimit = L / (8 * r * slowFactor);
    double veloRightLimit = L / (4 * slowFactor);
    double veloDiff = veloRightLimit - veloLeftLimit;
    double angleLeftLimit = 0, angleRightLimit = 2 * M_PI;
    double angleDiff = angleRightLimit - angleLeftLimit;
    double velo, angle;

    for (int i = 0; i < n; i++) {
        velo = veloLeftLimit + veloDiff * (rand() / (double)RAND_MAX);
        angle = angleLeftLimit + angleDiff * (rand() / (double)RAND_MAX);
        veloArray[2 * i] = velo * cos(angle);
        veloArray[2 * i + 1] = velo * sin(angle);
    }
    return veloArray;
}

