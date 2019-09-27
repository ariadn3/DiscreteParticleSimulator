#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "particle.h"

#define _SEED 3210

void randomiseParticles(particle_t** particleArray, int n, double L, double r);
double* generatePosition(int n, double L, double r);
double* generateVelocity(int n, double L, double r);
