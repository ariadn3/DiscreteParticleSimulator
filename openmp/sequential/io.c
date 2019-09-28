#include "io.h"

// Read in inputs from file and return initial simulation parameters as params_t
params_t* read_file(int slowFactor) {
    params_t* p = params();

    // Read in N, L, r, S and finally simulation mode
    scanf("%d\n%lf\n%lf\n%d\n", &(p->n), &(p->l), &(p->r), &(p->s));
    char* buffer = (char*) malloc(sizeof(char) * 140);
    scanf("%s\n", buffer);

    // Determine if this simulation will run in 'print' or 'perf' mode
    if(strcmp(buffer, "print") == 0) {
        p->willPrint = true;
    } else if (strcmp(buffer, "perf") == 0) {
        p->willPrint = false;
    } else {
        printf("Neither 'print' or 'perf' words are present. Exiting...\n");
        exit(1);
    }

    int i;
    double x, y, v_x, v_y;
    bool isInitialised = false;
    particle_t** particles = (particle_t**) malloc(p->n * sizeof(particle_t));

    // If initial positions and velocities of particles are provided, read them
    while (fgets(buffer, 140, stdin) != NULL) {
        isInitialised = true;
        sscanf(buffer, "%d %lf %lf %lf %lf", &i, &x, &y, &v_x, &v_y);
        particles[i] = build_particle(i, x, y, v_x / slowFactor, v_y / slowFactor);
    }

    // Otherwise randomise the initial positions and velocities
    if (!isInitialised) {
        randomiseParticles(particles, slowFactor, p->n, p->l, p->r);
    }

    p->particles = particles;
    return p;
}

void printAll(int n, int step, particle_t** particles) {
    // Parallelise this
    for (int i = 0; i < n; i++) {
        char* details = particle_string(particles[i]);
        printf("%d %s", step, details);
        free(details);
    }
}
