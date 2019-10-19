#include <math.h>

#include "structs.h"

#define EDGE_TOLERANCE 1e-14

__global__ void checkWallCollision(double, double, particle_t*);
__global__ void checkCollision(double, particle_t*, particle_t*);
__global__ void updateParticles(particle_t**, int, bool*);
__global__ void settleCollision(collision_t*, double, double);

