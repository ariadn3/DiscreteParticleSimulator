#include "filter.h"

// Filters the collisions according to the time that it took place
void filterCollisions(collision_t** collisionArray, bool* hasCollided,
        int* numCollisions) {
    qsort(collisionArray, *numCollisions, sizeof(collision_t*), &cmpCollision);

    int saveIndex = 0;
    collision_t* curCollision;
    for (int curIndex = 0; curIndex < *numCollisions; curIndex++) {
        curCollision = collisionArray[curIndex];

        if (hasCollided[curCollision->p->id]
                || (curCollision->q != NULL && hasCollided[curCollision->q->id])) {
            // Particle p has already collided OR particle q has already collided
            // -> discard this colision candidate
            free_collision(curCollision);
        } else {
            // Collision candidate is valid - marked p, q as collided
            hasCollided[curCollision->p->id] = true;

            if (curCollision-> q != NULL) hasCollided[curCollision->q->id] = true;
            // Re-use collision candidates array to store valid collisions
            collisionArray[saveIndex] = collisionArray[curIndex];
            saveIndex++;
        }
    }

    *numCollisions = saveIndex;
}

// Comparator for sorting collisions, earlier time then smaller particle 'p' id
int cmpCollision(const void* collisionA, const void* collisionB) {
    collision_t* firstCollision = *(collision_t**) collisionA;
    collision_t* secondCollision = *(collision_t**) collisionB;

    if (firstCollision->time == secondCollision->time)
        return (firstCollision->p->id < secondCollision->p->id) ? -1 : 1;
    else
        return (firstCollision->time < secondCollision->time) ? -1 : 1;
}
