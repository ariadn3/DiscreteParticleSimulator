#include "params.h"

params_t* read_file(void);
void simulate(void);

int main() {
    simulate();
}

void simulate() {
    params_t* params = read_file();
    
    printf("%d %lf %lf %d\n", params->n, params->l, params->r, params->s);

    for(int time = 0; time < params->s; time++) {
    }
}
