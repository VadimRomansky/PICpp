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

	if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
		removeEscapedParticles();
	}
}

void Simulation::removeEscapedParticles(){
	std::vector<Particle*>::iterator it = particles.end();
	it = it - 1;
	while(it != particles.begin()){
		Particle* particle = *it;
		std::vector<Particle*>::iterator prev = it - 1;
		if(particle->escaped){
			particles.erase(it);
		}
		it = prev;
	}
}

/*void Simulation::moveParticle(Particle* particle){
	Vector3d E = correlationEfield(particle)*fieldScale;
	Vector3d B = correlationBfield(particle)*fieldScale;

	Vector3d velocity = particle->velocity(speed_of_light_normalized);

	if(B.norm() > 0){
		double omega = -1.0/(particle->gammaFactor(speed_of_light_normalized)*particle->mass*speed_of_light_normalized/(particle->charge*B.norm()));
		double deltaPhi = omega*deltaT;
		Matrix3d* rotation = Matrix3d::createBasisByOneVector(B);
		Matrix3d* inverce = rotation->Inverse();

		Matrix3d lorentzRotation = Matrix3d(cos(deltaPhi), -sin(deltaPhi), 0, sin(deltaPhi), cos(deltaPhi), 0, 0, 0, 1);

		Vector3d newMomentum = (*rotation)*(lorentzRotation*((*inverce)*(particle->momentum)));
		//Vector3d newMomentum = (*inverce)*(lorentzRotation*((*rotation)*(particle->momentum)));

		delete rotation;
		delete inverce;

		particle->momentum = newMomentum;
	}

	//Vector3d velocity = particle->velocity(speed_of_light_normalized);
	Vector3d lorentzForce = velocity.vectorMult(B)/speed_of_light_normalized;
	//particle->momenltatum += (E + (velocity.vectorMult(B)/speed_of_light_normalized))*particle->charge*deltaT;

	double kw = 2*pi/xsize;

	double realEy = Eyamplitude*cos(kw*particle->x - omega*time);
	double realVzprotton = VzamplitudeProton*cos(kw*particle->x - omega*time);
	double realVzelectron = VzamplitudeElectron*cos(kw*particle->x - omega*time);

	double force1 = E.y + B0.x*velocity.z/speed_of_light_normalized;
	double force5 = Eyamplitude*cos(kw*particle->x) + B0.x*velocity.z/speed_of_light_normalized;
	double force2 = (E + lorentzForce).y;
	double force3 = (Eyamplitude +B0.x*VzamplitudeElectron/speed_of_light_normalized)*cos(kw*particle->x - omega*time);
	double force4 = (Eyamplitude +B0.x*VzamplitudeProton/speed_of_light_normalized)*cos(kw*particle->x - omega*time);

	double acceleration = particle->charge*(E + lorentzForce).y/particle->mass;
	double derVe = omega*sqrt((VyamplitudeElectron - velocity.y)*(VyamplitudeElectron + velocity.y));
	double derVp = omega*sqrt((VyamplitudeProton - velocity.y)*(VyamplitudeProton + velocity.y));

	//particle->momentum += (E + lorentzForce)*particle->charge*deltaT;
	particle->momentum += E*particle->charge*deltaT;
	Vector3d newVelocity = particle->velocity(speed_of_light_normalized);

	//particle->x += 0.5*(velocity.x + newVelocity.x)*deltaT;
	particle->momentum.x = 0;

	particle->y += 0.5*(velocity.y + newVelocity.y)*deltaT;
	particle->z += 0.5*(velocity.z + newVelocity.z)*deltaT;
}*/

void Simulation::moveParticle(Particle* particle){
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

	tempParticle.coordinates.x += ((1 - eta)*velocity.x + eta*newVelocity.x)*eta*deltaT;
	tempParticle.coordinates.y += ((1 - eta)*velocity.y + eta*newVelocity.y)*eta*deltaT;
	tempParticle.coordinates.z += ((1 - eta)*velocity.z + eta*newVelocity.z)*eta*deltaT;

if(boundaryConditionType != PERIODIC){
		if(tempParticle.coordinates.x < xgrid[0]){
			if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
				particle->coordinates.x = 2*xgrid[0] - tempParticle.coordinates.x;
				particle->coordinates.y = tempParticle.coordinates.y;
				particle->coordinates.z = tempParticle.coordinates.z;
				newVelocity.x = -newVelocity.x;
				particle->setMomentumByV(newVelocity, speed_of_light_normalized);
				return;
			} else {
				escapedParticles.push_back(particle);
				particle->coordinates.x = xgrid[0];
				particle->coordinates.y = tempParticle.coordinates.y;
				particle->coordinates.z = tempParticle.coordinates.z;
				particle->setMomentumByV(newVelocity, speed_of_light_normalized);
				particle->escaped = true;
				return;
			}
		}
		if(tempParticle.coordinates.x > xgrid[xnumber]){
			escapedParticles.push_back(particle);
			particle->coordinates.x = xgrid[xnumber];
			particle->coordinates.y = tempParticle.coordinates.y;
			particle->coordinates.z = tempParticle.coordinates.z;
			particle->setMomentumByV(newVelocity, speed_of_light_normalized);
			particle->escaped = true;
			return;
		}
	}


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

		tempParticle.coordinates.x += (middleVelocity.x*eta*deltaT);
		if(boundaryConditionType != PERIODIC){
			//todo more accurate speed!!
			if(tempParticle.coordinates.x < xgrid[0]){
				if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
					particle->coordinates.x = 2*xgrid[0] - tempParticle.coordinates.x + fabs(middleVelocity.x*(1 - eta)*deltaT);
					particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y*(1 - eta)*deltaT;
					particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z*(1 - eta)*deltaT;
					newVelocity = middleVelocity;
					newVelocity.x = -newVelocity.x;
					particle->setMomentumByV(newVelocity, speed_of_light_normalized);
					return;
				} else {
					escapedParticles.push_back(particle);
					particle->coordinates.x = xgrid[0];
					particle->coordinates.y = tempParticle.coordinates.y;
					particle->coordinates.z = tempParticle.coordinates.z;
					particle->setMomentumByV(middleVelocity, speed_of_light_normalized);
					particle->escaped = true;
					return;
				}
			}
			if(tempParticle.coordinates.x > xgrid[xnumber]){
				escapedParticles.push_back(particle);
				particle->coordinates.x = xgrid[xnumber];
				particle->coordinates.y = tempParticle.coordinates.y;
				particle->coordinates.z = tempParticle.coordinates.z;
				particle->setMomentumByV(middleVelocity, speed_of_light_normalized);
				particle->escaped = true;
				return;
			}
		}
		correctParticlePosition(tempParticle);

		E = correlationTempEfield(tempParticle)*fieldScale;
		B = correlationBfield(tempParticle)*fieldScale;

		tempParticle.momentum += (E + (middleVelocity.vectorMult(B)/speed_of_light_normalized))*particle->charge*deltaT;

		newVelocity = tempParticle.velocity(speed_of_light_normalized);

	}

	particle->momentum = tempParticle.momentum;
	//particle->momentum.x = 0;

	particle->coordinates.x += middleVelocity.x*deltaT;

	particle->coordinates.y += middleVelocity.y*deltaT;
	particle->coordinates.z += middleVelocity.z*deltaT;

	correctParticlePosition(particle);
}

void Simulation::correctParticlePosition(Particle* particle) {
	if(particle->coordinates.y < ygrid[0]) {
		particle->coordinates.y = particle->coordinates.y + ysize;
	}
	if(particle->coordinates.y > ygrid[ynumber]) {
		particle->coordinates.y = particle->coordinates.y - ysize;
	}		

	if(particle->coordinates.z < zgrid[0]) {
		particle->coordinates.z = particle->coordinates.z + zsize;
	}
	if(particle->coordinates.z > zgrid[znumber]) {
		particle->coordinates.z = particle->coordinates.z - zsize;
	}		

	if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
		if(particle->coordinates.x < xgrid[0]){
			particle->coordinates.x = 2*xgrid[0] - particle->coordinates.x;
			particle->momentum.x = -particle->momentum.x;
			return;
		}
		if(particle->coordinates.x > xgrid[xnumber]){
			escapedParticles.push_back(particle);
			particle->coordinates.x = xsize;
			particle->escaped = true;
			return;
		}
	}
	if(boundaryConditionType == PERIODIC){
		if(particle->coordinates.x < xgrid[0]) {
			particle->coordinates.x = particle->coordinates.x + xsize;
		}
		if(particle->coordinates.x > xgrid[xnumber]) {
			particle->coordinates.x = particle->coordinates.x - xsize;
		}		
	}
}

void Simulation::correctParticlePosition(Particle& particle) {
	if(particle.coordinates.y < ygrid[0]) {
		particle.coordinates.y = particle.coordinates.y + ysize;
	}
	if(particle.coordinates.y > ygrid[ynumber]) {
		particle.coordinates.y = particle.coordinates.y - ysize;
	}

	if(particle.coordinates.z < zgrid[0]) {
		particle.coordinates.z = particle.coordinates.z + zsize;
	}
	if(particle.coordinates.z > zgrid[znumber]) {
		particle.coordinates.z = particle.coordinates.z - zsize;
	}	

	if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
		if(particle.coordinates.x < xgrid[0]){
			particle.coordinates.x = 2*xgrid[0] - particle.coordinates.x;
			particle.momentum.x = -particle.momentum.x;
			return;
		}
		if(particle.coordinates.x > xgrid[xnumber]){
			escapedParticles.push_back(&particle);
			particle.coordinates.x = xsize;
			particle.escaped = true;
			return;
		}
	}
	if(boundaryConditionType == PERIODIC){
		if(particle.coordinates.x < xgrid[0]) {
			particle.coordinates.x = particle.coordinates.x + xsize;
		}
		if(particle.coordinates.x > xgrid[xnumber]) {
			particle.coordinates.x = particle.coordinates.x - xsize;
		}		
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

	beta = beta;

	double gamma_factor = 1 / sqrt(1 - (velocity.x * velocity.x + velocity.y * velocity.y + velocity.z * velocity.z) / speed_of_light_normalized_sqr);
	double G = ((beta * (EField.scalarMult(velocity)) / speed_of_light_normalized_sqr) + gamma_factor);
	beta = beta / G;
	double denominator = G * (1 + beta * beta * BField.scalarMult(BField)/speed_of_light_normalized_sqr);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			result.matrix[i][j] = Kronecker.matrix[i][j] + (beta * beta * BField[i] * BField[j]/speed_of_light_normalized_sqr);
			for (int k = 0; k < 3; ++k) {
				for (int l = 0; l < 3; ++l) {
					result.matrix[i][j] -= (beta * LeviCivita[j][k][l] * Kronecker.matrix[i][k] * BField[l]/speed_of_light_normalized);
				}
			}

			result.matrix[i][j] /= denominator;
		}
	}

	return result;
}

void Simulation::injectNewParticles(int count, ParticleTypeContainer typeContainer, double length){
	printf("inject new particles\n");
	double concentration = density/(massProton + massElectron);

	int n = particles.size();
	//Particle* tempParticle = particles[0];

	double x = xgrid[xnumber] - length;

	if(typeContainer.type == ELECTRON && preserveCharge){
		return;
	}

	for(int j = 0; j < ynumber; ++j){
		for(int k = 0; k < znumber; ++k){
			double weight = (typeContainer.concentration / typeContainer.particlesPerBin) * volumeB(xnumber - 1, j, k);
			for (int l = 0; l < 2 * count; ++l) {
				ParticleTypes type = typeContainer.type;
				Particle* particle = createParticle(n, xnumber - 1, j, k, weight, type, typeContainer, temperature);
				n++;
				particle->coordinates.x = x;
				particle->addVelocity(V0, speed_of_light_normalized);
				particle->initialMomentum = particle->momentum;
				particles.push_back(particle);
				double en = particle->energy(speed_of_light_normalized)*particle->weight*sqr(gyroradius/plasma_period);
				theoreticalEnergy += particle->energy(speed_of_light_normalized)*particle->weight*sqr(gyroradius/plasma_period);
				theoreticalMomentum += particle->momentum*particle->weight*gyroradius/plasma_period;
				if(preserveCharge){
					int necessaryElectrons = 1;
					if(type == ALPHA){
						necessaryElectrons = 2;
					}
					while(necessaryElectrons > 0){
						particle = createParticle(n, xnumber - 1, j, k, weight, ELECTRON, types[0], temperature);
						n++;
						particle->coordinates.x = x;
						particle->addVelocity(V0, speed_of_light_normalized);
						particle->initialMomentum = particle->momentum;
						particles.push_back(particle);
						en = particle->energy(speed_of_light_normalized)*particle->weight*sqr(gyroradius/plasma_period);
						theoreticalEnergy += particle->energy(speed_of_light_normalized)*particle->weight*sqr(gyroradius/plasma_period);
						theoreticalMomentum += particle->momentum*particle->weight*gyroradius/plasma_period;
						necessaryElectrons--;
					}
				}
			}
		}
	}
}
