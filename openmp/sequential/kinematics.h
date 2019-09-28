#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "particle.h"
#include "collision.h"

#define EDGE_TOLERANCE 1e-14

void resolveValidCollisions(collision_t**, int*, double, double);
void updateParticles(particle_t**, int, bool*);
void settleCollision(collision_t*, double, double);

