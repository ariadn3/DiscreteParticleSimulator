#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "particle.h"
#include "collision.h"

__host__ void filterCollisions(collision_t**, bool*, int*);
__host__ int cmpCollision(const void*, const void*);

