#include "random.h"

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

