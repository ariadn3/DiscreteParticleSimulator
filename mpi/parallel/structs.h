#ifndef STRUCTS_H
#define STRUCTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Collisions with a wall do not have a number
#define WALL -1

typedef struct {
    int id;
    double x;
    double y;
    double v_x;
    double v_y;
    int w_collisions;
    int p_collisions;
} particle_t;

typedef struct {
    int pId;
    int qId;
    double time;
} collision_t;

void collision_string(char*, collision_t*);
void particle_string(char*, particle_t*);
void particle_string_full(char*, particle_t*);

#endif

