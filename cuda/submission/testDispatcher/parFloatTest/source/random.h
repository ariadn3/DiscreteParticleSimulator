#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "particle.h"

#define SEED 3210

void randomiseParticles(particle_t**, int, int, float, float);
float* generatePosition(int, float, float);
float* generateVelocity(int, int, float, float);

