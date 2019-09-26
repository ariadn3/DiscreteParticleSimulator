#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "particle.h"
#include "collision.h"

void filterCollisions(collision_t** collisionArray, bool* hasCollided, int* numCollisions);
int cmpCollision(const void* collisionA, const void* collisionB);