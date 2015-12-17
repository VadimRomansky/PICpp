#include "stdlib.h"
#include "stdio.h"
#include "cmath"
#include <omp.h>

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"
#include "random.h"

void Simulation::moveParticles(){
	printf("moving particles\n");
	int i = 0;
#pragma omp parallel for private(i) 
	for(int i = 0; i < particles.size(); ++i){
		moveParticle(particles[i]);
	}

	/*for(int i = 0; i < particles.size(); ++i){
		scatterParticle(particles[i]);
	}*/
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

	Particle* tempParticle = *particles.begin();
	if(tempParticle->escaped){
		particles.erase(particles.begin());
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
	//Vector3d E = Vector3d(0, 0, 0);
	//Vector3d B = Vector3d(0, 0, 0);

	Vector3d velocity = particle->velocity(speed_of_light_normalized);
	Vector3d newVelocity = velocity;
	Vector3d middleVelocity = velocity;

	Vector3d electricForce = E*particle->charge*deltaT;
	Vector3d lorentzForce = (velocity.vectorMult(B)/speed_of_light_normalized)*particle->charge*deltaT;

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

	if(boundaryConditionType != PERIODIC){
		if(tempParticle.x < xgrid[0]){
			if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
				particle->x = 2*xgrid[0] - tempParticle.x;
				particle->y = tempParticle.y;
				particle->z = tempParticle.z;
				newVelocity.x = -newVelocity.x;
				particle->setMomentumByV(newVelocity, speed_of_light_normalized);
				return;
			} else {
				escapedParticles.push_back(particle);
				particle->x = xgrid[0];
				particle->y = tempParticle.y;
				particle->z = tempParticle.z;
				particle->setMomentumByV(newVelocity, speed_of_light_normalized);
				particle->escaped = true;
				return;
			}
		}
		if(tempParticle.x > xgrid[xnumber]){
			escapedParticles.push_back(particle);
			particle->x = xgrid[xnumber];
			particle->y = tempParticle.y;
			particle->z = tempParticle.z;
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

		tempParticle.x += (middleVelocity.x*eta*deltaT);
		if(boundaryConditionType != PERIODIC){
			//todo more accurate speed!!
			if(tempParticle.x < xgrid[0]){
				if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
					particle->x = 2*xgrid[0] - tempParticle.x + fabs(middleVelocity.x*(1 - eta)*deltaT);
					particle->y = tempParticle.y + middleVelocity.y*(1 - eta)*deltaT;
					particle->z = tempParticle.z + middleVelocity.z*(1 - eta)*deltaT;
					newVelocity = middleVelocity;
					newVelocity.x = -newVelocity.x;
					particle->setMomentumByV(newVelocity, speed_of_light_normalized);
					return;
				} else {
					escapedParticles.push_back(particle);
					particle->x = xgrid[0];
					particle->y = tempParticle.y;
					particle->z = tempParticle.z;
					particle->setMomentumByV(middleVelocity, speed_of_light_normalized);
					particle->escaped = true;
					return;
				}
			}
			if(tempParticle.x > xgrid[xnumber]){
				escapedParticles.push_back(particle);
				particle->x = xgrid[xnumber];
				particle->y = tempParticle.y;
				particle->z = tempParticle.z;
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

	particle->x += middleVelocity.x*deltaT;

	particle->y += middleVelocity.y*deltaT;
	particle->z += middleVelocity.z*deltaT;

	correctParticlePosition(particle);
}

void Simulation::correctParticlePosition(Particle* particle) {
	if(boundaryConditionType != PERIODIC){
		if(particle->x < xgrid[0]){
			if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
				particle->x = 2*xgrid[0] - particle->x;
				particle->momentum.x = -particle->momentum.x;
				return;
			} else {
				escapedParticles.push_back(particle);
				particle->x = xgrid[0];
				particle->escaped = true;
				return;
			}
		}
		if(particle->x > xgrid[xnumber]){
			escapedParticles.push_back(particle);
			particle->x = xgrid[xnumber];
			particle->escaped = true;
			return;
		}
	}
	if(boundaryConditionType == PERIODIC){
		if(particle->x < xgrid[0]) {
			particle->x = particle->x + xsize;
		}
		if(particle->x > xgrid[xnumber]) {
			particle->x = particle->x - xsize;
		}		
	}
}

void Simulation::correctParticlePosition(Particle& particle) {
	if(boundaryConditionType != PERIODIC){
		if(particle.x < xgrid[0]){
			if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
				particle.x = 2*xgrid[0] - particle.x;
				particle.momentum.x = -particle.momentum.x;
				return;
			} else {
				escapedParticles.push_back(&particle);
				particle.x = xgrid[0];
				particle.escaped = true;
				return;
			}
		}
		if(particle.x > xgrid[xnumber]){
			escapedParticles.push_back(&particle);
			particle.x = xgrid[xnumber];
			particle.escaped = true;
			return;
		}
	}
	if(boundaryConditionType == PERIODIC){
		if(particle.x < xgrid[0]) {
			particle.x = particle.x + xsize;
		}
		if(particle.x > xgrid[xnumber]) {
			particle.x = particle.x - xsize;
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
			//alertNaNOrInfinity(result.matrix[i][j], "rotation tensot = NaN\n");
		}
	}

	return result;
}

void Simulation::injectNewParticles(int count, double length){
	printf("inject new particles\n");
	double concentration = density/(massProton + massElectron);

	int n = particles.size();

	double weight = (concentration / particlesPerBin) * volumeB(xnumber - 1);
	double x = xgrid[xnumber] - length;
	Particle* lastProton = getLastProton();
	Particle* lastElectron = getLastElectron();
	for (int l = 0; l < 2 * count; ++l) {
		ParticleTypes type;
		if (l % 2 == 0) {
			type = PROTON;
		} else {
			type = ELECTRON;
		}
		Particle* particle = createParticle(n, xnumber - 1, weight, type, temperature);
		n++;
		particle->x = x;
		particle->addVelocity(V0, speed_of_light_normalized);
		/*if(type == PROTON){
			particle->momentum = lastProton->momentum;
		} else {
			particle->momentum = lastElectron->momentum;
		}*/
		particles.push_back(particle);
	}
}

void Simulation::scatterParticle(Particle* particle){
	double B = correlationBfield(particle).norm()*fieldScale;
	double lambda = fabs(particle->momentum.norm()*speed_of_light_normalized/(particle->charge*B));

	int i = (particle->x - xgrid[0])/deltaX;
	Vector3d oldMomentum = particle->momentum;
	Vector3d reverseBulkVelocity;
	Vector3d bulkVelocity;
	if(particle->type == PROTON || particle->type == ALPHA){
		bulkVelocity = velocityBulkProton[i];
		reverseBulkVelocity = Vector3d(-velocityBulkProton[i].x, -velocityBulkProton[i].y, -velocityBulkProton[i].z);
	} else if(particle->type == ELECTRON || particle->type == POSITRON){
		bulkVelocity = velocityBulkElectron[i];
		reverseBulkVelocity = Vector3d(-velocityBulkElectron[i].x, -velocityBulkElectron[i].y, -velocityBulkElectron[i].z);
	}
	particle->addVelocity(reverseBulkVelocity, speed_of_light_normalized);
	Vector3d particleVelocity = particle->velocity(speed_of_light_normalized);
	double probability = particleVelocity.norm()*deltaT/lambda;
	double maxTheta = 2*pi*lambda/particle->velocity(speed_of_light_normalized).norm();
	if(maxTheta > pi){
		maxTheta = pi;
	}

	//if(uniformDistribution() > (1 - probability)){
		Vector3d newVelocity;
		Matrix3d* rotation = Matrix3d::createBasisByOneVector(particleVelocity);
		
		//newVelocity.x = particleVelocity.norm()*2*(uniformDistribution()-0.5);
		newVelocity.z = particleVelocity.norm()*(1 - uniformDistribution()*(1 - cos(maxTheta)));
		double perpendicularVelocity = sqrt(particleVelocity.x* particleVelocity.x + 
											particleVelocity.y* particleVelocity.y +
											particleVelocity.z* particleVelocity.z -
											newVelocity.z*newVelocity.z);
		double phi = 2*pi*uniformDistribution();
		newVelocity.y = perpendicularVelocity*cos(phi);
		newVelocity.x = perpendicularVelocity*sin(phi);
		newVelocity = (*rotation)*newVelocity;
		particle->setMomentumByV(newVelocity, speed_of_light_normalized);

		delete rotation;
	//}

	particle->addVelocity(bulkVelocity, speed_of_light_normalized);
	Vector3d newMomentum = particle->momentum;
}
