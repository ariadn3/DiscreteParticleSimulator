#include "particle.h"

// Initialises a new particle
particle_t* build_particle(int id, double x, double y, double v_x, double v_y) {
    particle_t* particle = (particle_t*) malloc(sizeof(particle_t));
    
    particle->id = id;
    particle->x = x;
    particle->y = y;
    particle->v_x = v_x;
    particle->v_y = v_y;

    return particle;
}

// Destroys a built particle
void free_particle(particle_t* particle) {
    free(particle);
}

// Returns a C string with the details of this particle
char* toString(particle_t* p) {
    char* details = (char*) (malloc(sizeof(char) * 80));
    
    sprintf(details, "%d %.14f %.14lf %.14f %.14f\n", p->id, p->x, p->y, p->v_x, p->v_y);
    return details;
}
