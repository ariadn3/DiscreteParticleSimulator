#include <math.h>

#include "structs.h"

#define SEED 3210

__host__ void randomiseParticles(particle_t*, int, int, float, float);
__host__ float* generatePosition(int, float, float);
__host__ float* generateVelocity(int, int, float, float);

