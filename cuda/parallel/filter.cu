#include "filter.h"

// Filters the collisions according to the time that it took place
__host__ void filterCollisions(collision_t** collisionArray, bool* hasCollided,
        int* numCollisions) {
    // Quicksort all collision candidates with the comparator function
    qsort(collisionArray, *numCollisions, sizeof(collision_t*), &cmpCollision);

    int saveIndex = 0;
    collision_t* curCollision;
    for (int curIndex = 0; curIndex < *numCollisions; curIndex++) {
        curCollision = collisionArray[curIndex];
        
        // printf("=== Particle %d and %d collided ===\n", curCollision->p->id,
        //         curCollision->q == NULL ? -1 : curCollision->q->id);
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
__host__ int cmpCollision(const void* collisionA, const void* collisionB) {
    collision_t* firstCollision = *(collision_t**) collisionA;
    collision_t* secondCollision = *(collision_t**) collisionB;
    
    if (firstCollision->time == secondCollision->time) {
        // If both collisions involve the same first particle
        // Then prioritize wall collision, otherwise prioritize lower 2nd particle ID
        if (firstCollision->p->id == secondCollision->p->id) {
            if (firstCollision->q == NULL) return -1;
            else if (secondCollision->q == NULL) return 1;
            else return (firstCollision->q->id < secondCollision->q->id) ? -1 : 1;
        }
        // If two collisions occur at exactly the same time
        // Then prioritise the one which involves the particle P with lower ID
        return (firstCollision->p->id < secondCollision->p->id) ? -1 : 1;
    } else {
        // Otherwise prioritise the collision occurring at an earlier time
        return (firstCollision->time < secondCollision->time) ? -1 : 1;
    }
}

