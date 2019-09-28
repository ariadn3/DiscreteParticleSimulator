#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "particle.h"
#include "collision.h"

void filterCollisions(collision_t**, bool*, int*);
int cmpCollision(const void*, const void*);
