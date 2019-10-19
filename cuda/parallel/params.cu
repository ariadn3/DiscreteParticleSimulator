#include "params.h"

// Initialises a new parameter struct to store simulation parameters
__host__ params_t* params() {
    return (params_t*) malloc(sizeof(params_t));
}

// Destroys a built parameter struct
__host__ void free_params(params_t* params) {
    free(params);
}

