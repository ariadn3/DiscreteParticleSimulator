#include <stdbool.h>

#include "particle.h"

typedef struct params_t {
    int n;
    float l;
    float r;
    int s;
    bool willPrint;
    particle_t** particles;
} params_t;

params_t* params();
void free_params(params_t*);

