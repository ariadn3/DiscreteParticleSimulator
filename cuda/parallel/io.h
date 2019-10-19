#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "params.h"
#include "particle.h"
#include "random.h"

__host__ params_t* read_file(int slowFactor);
__host__ void printAll(bool, int, int, particle_t**);

