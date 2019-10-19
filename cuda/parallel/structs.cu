#include "structs.h"

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

// Initialises a new particle
__host__ particle_t* build_particle(int id, double x, double y, double v_x,
        double v_y) {
    particle_t* particle = (particle_t*) malloc(sizeof(particle_t));

    particle->id = id;
    particle->x = x;
    particle->y = y;
    particle->v_x = v_x;
    particle->v_y = v_y;
    particle->w_collisions = 0;
    particle->p_collisions = 0;

    return particle;
}

// Destroys a built particle
__host__ void free_particle(particle_t* particle) {
    free(particle);
}

// Returns a C string with the details of this particle
__host__ char* particle_string(particle_t* p) {
    char* details = (char*) (malloc(sizeof(char) * 140));

    sprintf(details, "%d %10.8lf %10.8lf %10.8lf %10.8lf\n",
            p->id, p->x, p->y, p->v_x, p->v_y);
    return details;
}

// Returns a C string with full details of this particle
__host__ char* particle_string_full(particle_t* p) {
    char* details = (char*) (malloc(sizeof(char) * 160));

    sprintf(details, "%d %10.8f %10.8lf %10.8lf %10.8lf %d %d\n",
            p->id, p->x, p->y, p->v_x, p->v_y, p->p_collisions, p->w_collisions);
    return details;
}

