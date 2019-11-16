#include "structs.h"

// Returns a C string with the details of this collision
__host__ char* collision_string(collision_t* c) {
    char* details = (char*) (malloc(sizeof(char) * 120));
    
    // Be careful of null-pointer exceptions since Q may be NULL for wall collisions
    sprintf(details, "Collision between particles with ids: %d and %d @ %.14f",
            c->p->id, c->q == NULL ? -1 : c->q->id, c->time);
    return details;
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

