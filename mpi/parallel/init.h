#include <math.h>

#include "structs.h"

#define SEED 3210

void randomiseParticles(particle_t*, int, int, double, double);
double* generatePosition(int, double, double);
double* generateVelocity(int, int, double, double);

