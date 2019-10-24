#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct particle_t {
    int id;
    float x;
    float y;
    float v_x;
    float v_y;
    int w_collisions;
    int p_collisions;
} particle_t;

particle_t* build_particle(int, float, float, float, float);
void free_particle(particle_t*);
char* particle_string(particle_t*);
char* particle_string_full(particle_t*);

#endif

