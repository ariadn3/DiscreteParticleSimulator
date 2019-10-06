#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "params.h"
#include "particle.h"
#include "random.h"

params_t* read_file(int slowFactor, int g);
void printAll(bool, int, int, particle_t**);

