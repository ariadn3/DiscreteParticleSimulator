#ifndef COLLISION_H
#define COLLISION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "particle.h"

typedef struct collision_t {
    particle_t* p;
    particle_t* q;
    float time;
} collision_t;

collision_t* build_collision(particle_t*, particle_t*, float);
void free_collision(collision_t*);
char* collision_string(collision_t*);

#endif

