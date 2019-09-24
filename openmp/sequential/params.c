#include "params.h"

// Initialises a new parameter struct to store simulation parameters
params_t* params(int n, double l, double r, int s, char* type, particle_t** particles) {
    return (params_t*) malloc(sizeof(params_t));
}

// Destroys a built parameter struct
void free_params(params_t* params) {
    params->type = NULL;
    free(params);
}
