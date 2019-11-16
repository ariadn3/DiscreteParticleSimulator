#include "structs.h"

// Returns a C string with the details of this collision
void collision_string(char* buff, collision_t* c) {
    // Be careful of null-pointer exceptions since Q may be NULL for wall collisions
    sprintf(details, "Collision between particles with ids: %d and %d @ %.14f",
            c->pId, c->qId == NULL ? -1 : c->qId, c->time);
}

// Returns a C string with the details of this particle
void particle_string(char* buff, particle_t* p) {
    sprintf(details, "%d %10.8lf %10.8lf %10.8lf %10.8lf\n",
            p->id, p->x, p->y, p->v_x, p->v_y);
}

// Returns a C string with full details of this particle
void particle_string_full(char* buff, particle_t* p) {
    sprintf(details, "%d %10.8f %10.8lf %10.8lf %10.8lf %d %d\n",
            p->id, p->x, p->y, p->v_x, p->v_y, p->p_collisions, p->w_collisions);
}

