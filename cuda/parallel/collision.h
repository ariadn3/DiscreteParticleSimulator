#ifndef COLLISION_H
#define COLLISION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "particle.h"

typedef struct collision_t {
    particle_t* p;
    particle_t* q;
    double time;
} collision_t;

__device__ collision_t* build_collision(particle_t*, particle_t*, double);
__device__ void free_collision(collision_t*);
__host__ char* collision_string(collision_t*);

#endif

