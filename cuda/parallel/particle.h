#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct particle_t {
    int id;
    double x;
    double y;
    double v_x;
    double v_y;
    int w_collisions;
    int p_collisions;
} particle_t;

__host__ particle_t* build_particle(int, double, double, double, double);
__host__ void free_particle(particle_t*);
__host__ char* particle_string(particle_t*);
__host__ char* particle_string_full(particle_t*);

#endif

