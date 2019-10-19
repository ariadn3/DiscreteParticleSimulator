#include <stdbool.h>

#include "particle.h"

typedef struct params_t {
    int n;
    double l;
    double r;
    int s;
    bool willPrint;
    particle_t** particles;
} params_t;

__host__ params_t* params();
__host__ void free_params(params_t*);

