#include "collision.h"

// Initialises a new collision
__device__ collision_t* build_collision(particle_t* p, particle_t* q, double time) {
    collision_t* collision = (collision_t*) malloc(sizeof(collision_t));

    collision->p = p;
    collision->q = q;
    collision->time = time;

    return collision;
}

// Destroys a built collision
__device__ void free_collision(collision_t* collision) {
    collision->p = NULL;
    collision->q = NULL;
    free(collision);
}

// Returns a C string with the details of this collision
__host__ char* collision_string(collision_t* c) {
    char* details = (char*) (malloc(sizeof(char) * 120));

    sprintf(details, "Collision between particles with ids: %d and %d @ %.14f",
            c->p->id, c->q->id, c->time);
    return details;
}

