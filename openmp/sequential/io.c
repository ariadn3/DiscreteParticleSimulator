#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"
#include "random.h"

// Read in inputs from file and return initial simulation parameters as params_t
params_t* read_file() {
    params_t* p = params();

    // Read in N, L, r, S respectively
    scanf("%d\n%f\n%f\n%d\n", &(p->n), &(p->l), &(p->r), &(p->l));
    
    char* buffer = (char*) malloc(sizeof(char) * 140);
    scanf("%s\n", &buffer);

    // Determine if this simulation will run in 'print' or 'perf' mode
    if(strcmp(buffer, "print")) {
        p->willPrint = true;
    } else if (strcmp(buffer, "perf")) {
        p->willPrint = false;
    } else {
        printf("Neither 'print' or 'perf' words are present in input. Exiting...\n");
        sleep(1);
        exit(1);
    }

    int i;
    double x, y, v_x, v_y;
    bool isInitialised = false;
    particle_t** particles = (particle_t**) malloc(p->n * sizeof(particle_t));
    
    // If initial positions and velocities of particles are provided, read them
    while (fgets(buffer, 140, stdin) != NULL) {
        isInitialised = true;
        sscanf(buffer, "%d %f %f %f %f", i, x, y, v_x, v_y);
        particles[i] = build_particle(i, x, y, v_x, v_y);
    }

    // Otherwise randomise the initial positions and velocities
    if (!isInitialised) {
        randomiseParticles(particles, p->n, p->l, p->r);
    }

    return p;
}

void write_file() {
    return;
}

int main() {
    return 1;
}
