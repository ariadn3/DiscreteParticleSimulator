#include "particle.h"

typedef struct params_t {
    int n;
    double l;
    double r;
    int s;
    char* type;
    particle_t** particles;
} params_t;

params_t* params(int, double, double, int, char*, particle_t**);
void free_params(params_t*);
