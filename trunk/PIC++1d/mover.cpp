#include "stdlib.h"
#include "stdio.h"
#include "cmath"
#include <omp.h>

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"

void Simulation::moveParticles(){
	printf("moving particles\n");
	int i = 0;
#pragma omp parallel for private(i) 
	for(i = 0; i < particles.size(); ++i){
		moveParticle(particles[i]);
	}
}

void Simulation::moveParticle(Particle* particle){
	//Vector3d oldE = correlationTempEfield(particle)*fieldScale
	Vector3d E = correlationTempEfield(particle)*fieldScale;
	Vector3d B = correlationBfield(particle)*fieldScale;

	Vector3d velocity = particle->velocity(speed_of_light_normalized);
	Vector3d newVelocity = velocity;
	Vector3d middleVelocity = velocity;

	//see Noguchi
	double beta = 0.5*particle->charge*deltaT/particle->mass;
	double Gamma = (beta*(E.scalarMult(velocity))/speed_of_light_normalized_sqr) + particle->gammaFactor(speed_of_light_normalized);
	double betaShift = beta/Gamma;

	int particleIterations = 20;


	Particle tempParticle = *particle;

	tempParticle.momentum += (E + (velocity.vectorMult(B)/speed_of_light_normalized))*particle->charge*deltaT;

	newVelocity = tempParticle.velocity(speed_of_light_normalized);

	tempParticle.x += ((1 - eta)*velocity.x + eta*newVelocity.x)*eta*deltaT;
	tempParticle.y += ((1 - eta)*velocity.y + eta*newVelocity.y)*eta*deltaT;
	tempParticle.z += ((1 - eta)*velocity.z + eta*newVelocity.z)*eta*deltaT;


	Vector3d prevVelocity = velocity; 
	int i = 0;

	while((prevVelocity - newVelocity).norm() > 1E-16*velocity.norm() && i < particleIterations){
		++i;
		prevVelocity = newVelocity;

		tempParticle = *particle;
		Vector3d rotatedE = tempParticle.rotationTensor*E;

		//tempParticle.momentum += (E + ((velocity + newVelocity).vectorMult(B)/(2.0*speed_of_light_normalized)))*particle->charge*deltaT;

		//mistake in noguchi - he writes betashift!
		middleVelocity = (tempParticle.rotationTensor*tempParticle.gammaFactor(speed_of_light_normalized)*velocity) + rotatedE*beta;

		tempParticle.x += (middleVelocity.x*theta*deltaT);
		correctParticlePosition(tempParticle);

		E = correlationTempEfield(tempParticle)*fieldScale;
		B = correlationBfield(tempParticle)*fieldScale;

		tempParticle.momentum += (E + (middleVelocity.vectorMult(B)/speed_of_light_normalized))*particle->charge*deltaT;

		newVelocity = tempParticle.velocity(speed_of_light_normalized);

	}

	particle->momentum = tempParticle.momentum;
	particle->momentum.x = 0;

	//particle->x += middleVelocity.x*deltaT;

	particle->y += middleVelocity.y*deltaT;
	particle->z += middleVelocity.z*deltaT;

	correctParticlePosition(particle);
}

void Simulation::correctParticlePosition(Particle* particle) {
	correctParticlePosition(*particle);
}

void Simulation::correctParticlePosition(Particle& particle) {
		if(particle.x < 0) {
			particle.x = particle.x + xsize;
		}
		if(particle.x > xsize) {
			particle.x = particle.x - xsize;
		}		
}

void Simulation::evaluateParticlesRotationTensor()
{
	for(int i = 0; i < particles.size(); ++i)
	{
		Particle* particle = particles[i];
		double beta = 0.5*particle->charge*deltaT/particle->mass;
		Vector3d velocity = particle->velocity(speed_of_light_normalized);

		Vector3d oldE = correlationEfield(particle)*fieldScale;
		Vector3d oldB = correlationBfield(particle)*fieldScale;

		particle->rotationTensor = evaluateAlphaRotationTensor(beta, velocity, oldE, oldB);

	}
}

Matrix3d Simulation::evaluateAlphaRotationTensor(double beta, Vector3d velocity, Vector3d EField, Vector3d BField) {
	Matrix3d result;

	beta = beta / speed_of_light_normalized;

	double gamma_factor = 1 / sqrt(1 - (velocity.x * velocity.x + velocity.y * velocity.y + velocity.z * velocity.z) / speed_of_light_normalized_sqr);
	double G = ((beta * (EField.scalarMult(velocity)) / speed_of_light_normalized) + gamma_factor);
	beta = beta / G;
	double denominator = G * (1 + beta * beta * BField.scalarMult(BField));

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			result.matrix[i][j] = Kronecker.matrix[i][j] + beta * beta * BField[i] * BField[j];
			for (int k = 0; k < 3; ++k) {
				for (int l = 0; l < 3; ++l) {
					result.matrix[i][j] -= beta * LeviCivita[j][k][l] * Kronecker.matrix[i][k] * BField[l];
				}
			}

			result.matrix[i][j] /= denominator;
		}
	}

	return result;
}
