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
    int cell_x;
    int cell_y;
} particle_t;

particle_t* build_particle(int, int, double, double, double, double);
void free_particle(particle_t*);
char* particle_string(particle_t*);
char* particle_string_full(particle_t*);

#endif

