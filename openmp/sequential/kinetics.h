#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "particle.h"
#include "collision.h"

#define _EDGE_TOLERANCE 1e-14

void resolveValidCollisions(collision_t** collisionArray, int* numCollisions, double L, double r);
void updateParticles(particle_t** particleArray, int N, bool* hasCollided);
void settleCollision(collision_t* curCollision, double L, double r);