#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "particle.h"

#define SEED 3210

void randomiseParticles(particle_t**, int, int, double, double);
double* generatePosition(int, double, double);
double* generateVelocity(int, int, double, double);

