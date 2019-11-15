#include "init.h"

// Randomly generates an array of particles
__host__ void randomiseParticles(particle_t* particleArray, int slowFactor, int n,
        double L, double r) {
    double* posArray = generatePosition(n, L, r);
    double* veloArray = generateVelocity(slowFactor, n, L, r);
    
    for (int i = 0; i < n; i++) {
        particleArray[i].id = i;
        particleArray[i].x = posArray[2 * i];
        particleArray[i].y = posArray[2 * i + 1];
        particleArray[i].v_x = veloArray[2 * i];
        particleArray[i].v_y = veloArray[2 * i + 1];
        particleArray[i].w_collisions = 0;
        particleArray[i].p_collisions = 0;
    }

    free(posArray);
    free(veloArray);
}

// Generates an array of non-overlapping position values
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
    double minPos = r, maxPos = L - r, posRange = maxPos - minPos;
    for (int i = 0; i < n; i++) {
        while (true) {
            posArray[2 * i] = minPos + posRange * (rand() / (double)RAND_MAX);
            posArray[2 * i +1] = minPos + posRange * (rand() / (double)RAND_MAX);
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
    double minVelocity = L / (8 * r * slowFactor);
    double maxVelocity = L / (4 * slowFactor);
    double velocityRange = maxVelocity - minVelocity;
    double minPolarAngle = 0, maxPolarAngle = 2 * M_PI;
    double angleRange = maxPolarAngle - minPolarAngle;
    double v, theta;

    for (int i = 0; i < n; i++) {
        v = minVelocity + velocityRange * (rand() / (double)RAND_MAX);
        theta = minPolarAngle + angleRange * (rand() / (double)RAND_MAX);
        veloArray[2 * i] = v * cos(theta);
        veloArray[2 * i + 1] = v * sin(theta);
    }
    return veloArray;
}

