#include "random.h"

// Randomly generates an array of particles
void randomiseParticles(particle_t** particleArray, int slowFactor, int n, float L,
        float r) {
    float* posArray = generatePosition(n, L, r);
    float* veloArray = generateVelocity(slowFactor, n, L, r);
    // ===== DO NOT PARALLELISE =====
    // Work is too trivial that runtime gets destroyed by overhead of //isation 
    for (int i = 0; i < n; i++) {
        particleArray[i] = build_particle(i, posArray[2 * i], posArray[2 * i + 1],
                veloArray[2 * i], veloArray[2 * i + 1]);
    }
    free(posArray);
    free(veloArray);
}

// Generates an array of position values
float* generatePosition(int n, float L, float r) {
    static float* posArray;
    posArray = (float*) malloc(n * 2 * sizeof(float));
    srand(SEED);
    
    // Check pre-conditions (read our report)
    if (L < 2 * r) {
        // No way to place any particle with diameter larger than the box
        printf("Assumption violated: L < (2 * r)\n");
        exit(1);
    } else if (n * r * r > L * L) {
        // No way to place n particles safely in a box of length L
        printf("Assumption violated: (n * r * r) > (L * L)\n");
        exit(1);
    }
    
    // Bounds for generating positions
    float leftLimit = r, rightLimit = L - r, lenDiff = rightLimit - leftLimit;
    // ===== CANNOT PARALLELISE =====
    // Significant work is done here attempting to place all n particles and checking
    // that each particle does not overlap with any previously placed particle
    // Not possible to //ise due to sequential nature of placement
    for (int i = 0; i < n; i++) {
        while (true) {
            posArray[2 * i] = leftLimit + lenDiff * (rand() / (float)RAND_MAX);
            posArray[2 * i +1] = leftLimit + lenDiff * (rand() / (float)RAND_MAX);
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
float* generateVelocity(int slowFactor, int n, float L, float r) {
    static float* veloArray;
    veloArray = (float*) malloc(n * 2 * sizeof(float));
    
    // Bounds for generating velocities
    float veloLeftLimit = L / (8 * r * slowFactor);
    float veloRightLimit = L / (4 * slowFactor);
    float veloDiff = veloRightLimit - veloLeftLimit;
    float angleLeftLimit = 0, angleRightLimit = 2 * M_PI;
    float angleDiff = angleRightLimit - angleLeftLimit;
    float velo, angle;

    // ===== DO NOT PARALLELISE =====
    // Work is too trivial that runtime gets destroyed by overhead of //isation 
    for (int i = 0; i < n; i++) {
        velo = veloLeftLimit + veloDiff * (rand() / (float)RAND_MAX);
        angle = angleLeftLimit + angleDiff * (rand() / (float)RAND_MAX);
        veloArray[2 * i] = velo * cos(angle);
        veloArray[2 * i + 1] = velo * sin(angle);
    }
    return veloArray;
}

