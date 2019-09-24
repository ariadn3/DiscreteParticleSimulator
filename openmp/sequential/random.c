#include "random.h"

// Randomly generates an array of particles
particle_t** randomiseParticles(particle_t** particleArray, int n, double L, double r) {
	double* posArray = generatePosition(n, L, r);
	double* veloArray = generateVelocity(n, L, r);
	for (int i = 0; i < n; i++)
		particleArray[i] = build_particle(i, posArray[2*i], posArray[2*i+1], veloArray[2*i], veloArray[2*i+1]);
	free(posArray);
	free(veloArray);

	return particleArray;
}

// Generates an array of position values
double* generatePosition(int n, double L, double r) {
	static double* posArray;
	posArray = (double*) malloc(n*2*sizeof(double));
	srand(_SEED);

	if (L < 2*r) {
		printf("Invalid parameters (L<2r)\n");
		exit(1);
	}

	double leftLimit = r, rightLimit = L-r, lenDiff = rightLimit - leftLimit;
	for (int i = 0; i < 2*n; i++)
		posArray[i] = leftLimit + lenDiff*(rand()/(double)RAND_MAX);

	return posArray;
}

// Generates an array of velocity values
double* generateVelocity(int n, double L, double r) {
	static double* veloArray;
	veloArray = (double*) malloc(n*2*sizeof(double));

	double veloLeftLimit = L/(8*r), veloRightLimit = L/4, veloDiff = veloRightLimit - veloLeftLimit;
	double angleLeftLimit = 0, angleRightLimit = 2*M_PI, angleDiff = angleRightLimit - angleLeftLimit;
	double velo, angle;
	for (int i = 0; i < n; i++) {
		velo = veloLeftLimit + veloDiff*(rand()/(double)RAND_MAX);
		angle = angleLeftLimit + angleDiff*(rand()/(double)RAND_MAX);
		veloArray[2*i] = velo*cos(angle);
		veloArray[2*i+1] = velo*sin(angle);
	}
	return veloArray;
}