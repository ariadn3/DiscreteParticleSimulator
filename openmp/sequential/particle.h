#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct particle_t {
    int id;
    double x;
    double y;
    double v_x;
    double v_y;
} particle_t;

particle_t* build_particle(int id, double x, double y, double v_x, double v_y);
void free_particle(particle_t* particle);
char* toString(particle_t* particle);
