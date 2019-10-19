#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

typedef struct particle_t {
    int id;
    double x;
    double y;
    double v_x;
    double v_y;
    int w_collisions;
    int p_collisions;
} particle_t;

typedef struct collision_t {
    particle_t* p;
    particle_t* q;
    double time;
} collision_t;

__device__ collision_t* build_collision(particle_t*, particle_t*, double);
__device__ void free_collision(collision_t*);
__host__ char* collision_string(collision_t*);

__host__ particle_t* build_particle(int, double, double, double, double);
__host__ void free_particle(particle_t*);
__host__ char* particle_string(particle_t*);
__host__ char* particle_string_full(particle_t*);

