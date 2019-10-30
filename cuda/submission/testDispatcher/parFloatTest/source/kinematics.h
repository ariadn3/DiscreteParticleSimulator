#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "particle.h"
#include "collision.h"

#define EDGE_TOLERANCE 1e-8

void resolveValidCollisions(collision_t**, int*, float, float);
void updateParticles(particle_t**, int, bool*);
void settleCollision(collision_t*, float, float);

