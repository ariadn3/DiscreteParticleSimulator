#include "particle.h"

// Initialises a new particle
particle_t* build_particle(int id, double x, double y, double v_x, double v_y) {
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
void free_particle(particle_t* particle) {
    free(particle);
}

// Returns a C string with the details of this particle
char* particle_string(particle_t* p) {
    char* details = (char*) (malloc(sizeof(char) * 140));

    sprintf(details, "%d %10.8lf %10.8lf %10.8lf %10.8lf\n",
            p->id, p->x, p->y, p->v_x, p->v_y);
    return details;
}

// Returns a C string with full details of this particle
char* particle_string_full(particle_t* p) {
    char* details = (char*) (malloc(sizeof(char) * 160));

    sprintf(details, "%d %10.8f %10.8lf %10.8lf %10.8lf %d %d\n",
            p->id, p->x, p->y, p->v_x, p->v_y, p->p_collisions, p->w_collisions);
    return details;
}

