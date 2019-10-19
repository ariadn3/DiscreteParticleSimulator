#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "particle.h"
#include "collision.h"

#define EDGE_TOLERANCE 1e-14

__host__void resolveValidCollisions(collision_t**, int*, double, double);
__global__ void updateParticles(particle_t**, int, bool*);
__global__ void settleCollision(collision_t*, double, double);

