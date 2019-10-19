#include <math.h>

#include "structs.h"

#define SEED 3210

__host__ void randomiseParticles(particle_t**, int, int, double, double);
__host__ double* generatePosition(int, double, double);
__host__ double* generateVelocity(int, int, double, double);

