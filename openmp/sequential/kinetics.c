#include "kinetics.h"

// Moves all valid collisions in the collision array
void resolveValidCollisions(collision_t** collisionArray, int* numCollisions, double L, double r) {
	collision_t* curCollision;
	for (int i = 0; i < *numCollisions; i++) {
        curCollision = collisionArray[i];
		settleCollision(curCollision, L, r);
		free_collision(curCollision);
	}
	return;
}

// Updates particles not involved in any collision
void updateParticles(particle_t** particleArray, int N, bool* hasCollided) {
	particle_t* curParticle;
	for (int i = 0; i < N; i++) {
		curParticle = particleArray[i];
        if (!hasCollided[i]) {
			curParticle->x += curParticle->v_x;
			curParticle->y += curParticle->v_y;
		} else {
            hasCollided[i] = false;
        }
	}
	return;
}

// Moves two particles involved in a collision to their rightful place after the timestep
void settleCollision(collision_t* curCollision, double L, double r) {
	particle_t* particleA;
	particle_t* particleB;
	particleA = curCollision->p;
	particleB = curCollision->q;
	double time = curCollision->time;
    
	particleA->x += time * particleA->v_x;
	particleA->y += time * particleA->v_y;
	// If the collision is against the wall, toggle directions
	if (particleB == NULL) {
		if (particleA->x <= r + _EDGE_TOLERANCE || particleA->x >= L - r - _EDGE_TOLERANCE) {
            particleA->v_x *= -1;
        }
		if (particleA->y <= r + _EDGE_TOLERANCE || particleA->y >= L - r - _EDGE_TOLERANCE)
			particleA->v_y *= -1;
	}
	// If collision is against another particle
	else {
		particleB->x += time * particleB->v_x;
		particleB->y += time * particleB->v_y;
	
		// Normal and tangent unit vectors
		double distance = sqrt(pow(particleB->x - particleA->x, 2) + pow(particleB->y - particleA->y, 2));
		double n_x = (particleB->x - particleA->x)/distance, n_y = (particleB->y - particleA->y)/distance;
		double t_x = -n_y, t_y = n_x;

		double v_an = n_x * particleA->v_x + n_y * particleA->v_y, v_at = t_x * particleA->v_x + t_y * particleA->v_y;
		double v_bn = n_x * particleB->v_x + n_y * particleB->v_y, v_bt = t_x * particleB->v_x + t_y * particleB->v_y;

		// printf("n_x = %.14f, n_y = %.14f\n", n_x, n_y);
		// printf("t_x = %.14f, t_y = %.14f\n", t_x, t_y);
		// printf("v_an = %.14f, v_at = %.14f\n", v_an, v_at);
		// printf("v_bn = %.14f, v_bt = %.14f\n", v_bn, v_bt);

		// printf("Pre-collision velocities: %.14f, %.14f, %.14f, %.14f\n", particleA->v_x, particleA->v_y, particleB->v_x, particleB->v_y);

		// Update resultant velocities
		particleA->v_x = v_bn*n_x + v_at*t_x;
		particleA->v_y = v_bn*n_y + v_at*t_y;
		particleB->v_x = v_an*n_x + v_bt*t_x;
		particleB->v_y = v_an*n_y + v_bt*t_y;

		// printf("Post-collision velocities: %.14f, %.14f, %.14f, %.14f\n", particleA->v_x, particleA->v_y, particleB->v_x, particleB->v_y);

		// If particle B will collide against the wall, check when it will collide with the nearest wall and take that time
		double time_bx = 1-time, time_by = 1-time;
		if (particleB->x + time_bx*particleB->v_x < r)
			time_bx = -(particleB->x-r)/particleB->v_x;
		else if (particleB->x + time_bx*particleB->v_x > L-r)
			time_bx = (L-r-particleB->x)/particleB->v_x;
		if (particleB->y + time_by*particleB->v_y < r)
			time_by = -(particleB->y-r)/particleB->v_y;
		else if (particleB->y + time_by*particleB->v_y > L-r)
			time_by = (L-r-particleB->y)/particleB->v_y;
		double time_b = (time_bx < time_by)? time_bx : time_by;

		particleB->x += time_b * particleB->v_x;
		particleB->y += time_b * particleB->v_y;
	}

	// If particle A will collide against the wall, check when it will collide with the nearest wall and take that time
	double time_ax = 1-time, time_ay = 1-time;
	if (particleA->x + time_ax*particleA->v_x < r)
		time_ax = -(particleA->x-r)/particleA->v_x;
	else if (particleA->x + time_ax*particleA->v_x > L-r)
		time_ax = (L-r-particleA->x)/particleA->v_x;
	if (particleA->y + time_ay*particleA->v_y < r)
		time_ay = -(particleA->y-r)/particleA->v_y;
	else if (particleA->y + time_ay*particleA->v_y > L-r)
		time_ay = (L-r-particleA->y)/particleA->v_y;
	double time_a = (time_ax < time_ay)? time_ax : time_ay;
	particleA->x += time_a * particleA->v_x;
	particleA->y += time_a * particleA->v_y;
    return;
}
