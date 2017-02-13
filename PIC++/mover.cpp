#include "stdlib.h"
#include "stdio.h"
#include "cmath"
#include <omp.h>
#include <time.h>
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "util.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "random.h"
#include "simulation.h"
#include "mpi_util.h"

void Simulation::moveParticles() {
	double procTime = 0;
	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	if ((rank == 0) && (verbosity > 0)) printf("moving particles\n");
	if ((rank == 0) && (verbosity > 0)) printLog("moving particles\n");
	int i = 0;

	for (i = 0; i < particles.size(); ++i) {
		if (i % 1000 == 0) {
			if (verbosity > 2) {
				printf("move particle number %d rank %d\n", i, rank);
			}
		}
		moveParticle(particles[i]);
		//particles[i]->momentum.x = 0;
	}

	//if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
	//removeEscapedParticles();
	//}
	if ((rank == 0) && (verbosity > 0)) printf("end moving particles\n");
	if ((rank == 0) && (verbosity > 0)) printLog("end moving particles\n");
	//MPI_Barrier(MPI_COMM_WORLD);
	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("moving particles time = %g sec\n", procTime/CLOCKS_PER_SEC);
	}
}

void Simulation::removeEscapedParticles() {
	double procTime= 0;
	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock();
	}
	if (nprocs > 1) {
		for (int i = 0; i < escapedParticlesLeft.size(); ++i) {
			Particle* particle = escapedParticlesLeft[i];
			delete particle;
		}
	}
	escapedParticlesLeft.clear();

	if (nprocs > 1) {
		for (int i = 0; i < escapedParticlesRight.size(); ++i) {
			Particle* particle = escapedParticlesRight[i];
			delete particle;
		}
	}
	escapedParticlesRight.clear();

	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("removing escaped particles time = %g sec\n", procTime/CLOCKS_PER_SEC);
	}
}

void Simulation::eraseEscapedPaticles() {
	double procTime = 0;
	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if (particles.size() > 0) {
		std::vector<Particle*>::iterator it = particles.end();
		it = it - 1;
		while (it != particles.begin()) {
			Particle* particle = *it;
			std::vector<Particle*>::iterator prev = it - 1;
			if (particle->escaped) {
				chargeBalance -= particle->chargeCount;
				particles.erase(it);
			}
			it = prev;
		}
		Particle* particle = particles[0];
		if (particle->escaped) {
			chargeBalance -= particle->chargeCount;
			particles.erase(particles.begin());
		}
	}
	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("erasing escaped particles time = %g sec\n", procTime/CLOCKS_PER_SEC);
	}
}

void Simulation::moveParticle(Particle* particle) {
 	updateCorrelationMaps(particle);
	Vector3d E = correlationTempEfield(particle);
	//Vector3d E = correlationNewEfield(particle) * fieldScale;
	Vector3d B = correlationBfield(particle);
	//printf("E = %g %g %g\n", E.x, E.y, E.z);
	//printf("B = %g %g %g\n", B.x, B.y, B.z);

	Vector3d velocity = particle->getVelocity(speed_of_light_normalized);
	Vector3d newVelocity = velocity;
	Vector3d middleVelocity = velocity;


	//see Noguchi
	double beta = 0.5 * particle->charge * deltaT / particle->mass;

	int particleIterations = 20;


	Particle tempParticle = *particle;

	tempParticle.addMomentum( (E + (velocity.vectorMult(B) / speed_of_light_normalized)) * particle->charge * deltaT);
	//if (debugMode) alertNaNOrInfinity(tempParticle.momentum.x, "p.x = naN\n");
	//if (debugMode) alertNaNOrInfinity(tempParticle.momentum.y, "p.y = naN\n");
	//if (debugMode) alertNaNOrInfinity(tempParticle.momentum.z, "p.z = naN\n");


	newVelocity = tempParticle.getVelocity(speed_of_light_normalized);

	middleVelocity = velocity * (1 - eta) + newVelocity * eta;

	tempParticle.coordinates.x += middleVelocity.x * eta * deltaT;
	tempParticle.coordinates.y += middleVelocity.y * eta * deltaT;
	tempParticle.coordinates.z += middleVelocity.z * eta * deltaT;


	//if (boundaryConditionType != PERIODIC) {
	if ((tempParticle.coordinates.x < xgrid[1]) && (boundaryConditionType == SUPER_CONDUCTOR_LEFT) && (rank == 0)) {
		particle->coordinates.x = 2 * xgrid[1] - tempParticle.coordinates.x + fabs(
			middleVelocity.x * (1 - eta) * deltaT);
		particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * (1 - eta) * deltaT;
		particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * (1 - eta) * deltaT;

		if (particle->coordinates.y > ygrid[ynumber]) {
			particle->coordinates.y -= ysize;
		}
		if (particle->coordinates.y < ygrid[0]) {
			particle->coordinates.y += ysize;
		}

		if (particle->coordinates.z > zgrid[znumber]) {
			particle->coordinates.z -= zsize;
		}
		if (particle->coordinates.z < zgrid[0]) {
			particle->coordinates.z += zsize;
		}
		newVelocity.x = -newVelocity.x;
		particle->setMomentumByV(newVelocity, speed_of_light_normalized);

		//if (debugMode) alertNaNOrInfinity(particle->momentum.x, "p.x = NaN in simulation::moveParticle\n");
		//if (debugMode) alertNaNOrInfinity(particle->momentum.y, "p.y = NaN in simulation::moveParticle\n");
		//if (debugMode) alertNaNOrInfinity(particle->momentum.z, "p.z = NaN in simulation::moveParticle\n");
		return;
	}

	//if (tempParticle.coordinates.x < xgrid[0]) {
	if (tempParticle.coordinates.x < xgrid[1]) {

		escapedParticlesLeft.push_back(particle);
		particle->coordinates.x = tempParticle.coordinates.x + middleVelocity.x * (1 - eta) * deltaT;

		particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * (1 - eta) * deltaT;
		if (particle->coordinates.y > ygrid[ynumber]) {
			particle->coordinates.y -= ysize;
		}
		if (particle->coordinates.y < ygrid[0]) {
			particle->coordinates.y += ysize;
		}

		particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * (1 - eta) * deltaT;
		if (particle->coordinates.z > zgrid[znumber]) {
			particle->coordinates.z -= zsize;
		}
		if (particle->coordinates.z < zgrid[0]) {
			particle->coordinates.z += zsize;
		}
		particle->setMomentumByV(newVelocity, speed_of_light_normalized);
		particle->escaped = true;
		particle->crossBoundaryCount++;
		//if (debugMode) alertNaNOrInfinity(particle->momentum.x, "p.x = NaN in simulation::moveParticle\n");
		//if (debugMode) alertNaNOrInfinity(particle->momentum.y, "p.y = NaN in simulation::moveParticle\n");
		//if (debugMode) alertNaNOrInfinity(particle->momentum.z, "p.z = NaN in simulation::moveParticle\n");
		return;

	}
	//if (tempParticle.coordinates.x > xgrid[xnumber+1]) {
	if (tempParticle.coordinates.x > xgrid[xnumber]) {
		escapedParticlesRight.push_back(particle);
		particle->coordinates.x = tempParticle.coordinates.x + middleVelocity.x * (1 - eta) * deltaT;

		particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * (1 - eta) * deltaT;
		if (particle->coordinates.y > ygrid[ynumber]) {
			particle->coordinates.y -= ysize;
		}
		if (particle->coordinates.y < ygrid[0]) {
			particle->coordinates.y += ysize;
		}

		particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * (1 - eta) * deltaT;
		if (particle->coordinates.z > zgrid[znumber]) {
			particle->coordinates.z -= zsize;
		}
		if (particle->coordinates.z < zgrid[0]) {
			particle->coordinates.z += zsize;
		}
		particle->setMomentumByV(newVelocity, speed_of_light_normalized);
		particle->escaped = true;
		particle->crossBoundaryCount++;
		//if (debugMode) alertNaNOrInfinity(particle->momentum.x, "p.x = NaN in simulation::moveParticle\n");
		//if (debugMode) alertNaNOrInfinity(particle->momentum.y, "p.y = NaN in simulation::moveParticle\n");
		//if (debugMode) alertNaNOrInfinity(particle->momentum.z, "p.z = NaN in simulation::moveParticle\n");
		return;
	}
	//}


	Vector3d prevVelocity = velocity;
	int i = 0;
	Vector3d velocityHat = (tempParticle.rotationTensor * tempParticle.gammaFactor(
			speed_of_light_normalized) * velocity);

	double etaDeltaT = eta*deltaT;
	double restEtaDeltaT = (1.0 - eta)*deltaT;
	double velocityNorm = velocity.norm();

	double error = (prevVelocity - newVelocity).norm();
	while (error > particleVelocityErrorLevel * velocityNorm && i < particleIterations) {
		++i;
		prevVelocity = newVelocity;

		tempParticle = *particle;
		Vector3d rotatedE = tempParticle.rotationTensor * E;

		//tempParticle.momentum += (E + ((getVelocity + newVelocity).vectorMult(B)/(2.0*speed_of_light_normalized)))*particle->charge*deltaT;

		//mistake in noguchi - he writes betashift!
		
		middleVelocity = velocityHat + rotatedE * beta;
		/*if(middleVelocity.norm() > speed_of_light_normalized) {
		    printf("middleVelocity = %g\n", middleVelocity.norm());
		    printf("speed of light = %g\n", speed_of_light_normalized);
		    printf("particle number %d\n", particle->number);
		}*/

		tempParticle.coordinates.x += (middleVelocity.x * etaDeltaT);
		tempParticle.coordinates.y += (middleVelocity.y * etaDeltaT);
		tempParticle.coordinates.z += (middleVelocity.z * etaDeltaT);
		//if (boundaryConditionType != PERIODIC) {
		//todo more accurate speed!!
		if ((tempParticle.coordinates.x < xgrid[1]) && (boundaryConditionType == SUPER_CONDUCTOR_LEFT) && (rank == 0)) {
			particle->coordinates.x = 2 * xgrid[1] - tempParticle.coordinates.x + fabs(
				middleVelocity.x * restEtaDeltaT);
			particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * restEtaDeltaT;
			particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * restEtaDeltaT;
			if (particle->coordinates.y > ygrid[ynumber]) {
				particle->coordinates.y -= ysize;
			}
			if (particle->coordinates.y < ygrid[0]) {
				particle->coordinates.y += ysize;
			}
			if (particle->coordinates.z > zgrid[znumber]) {
				particle->coordinates.z -= zsize;
			}
			if (particle->coordinates.z < zgrid[0]) {
				particle->coordinates.z += zsize;
			}
			newVelocity.x = -newVelocity.x;
			particle->setMomentumByV(newVelocity, speed_of_light_normalized);
			//if (debugMode) alertNaNOrInfinity(particle->momentum.x, "p.x = NaN in simulation::moveParticle\n");
			//if (debugMode) alertNaNOrInfinity(particle->momentum.y, "p.y = NaN in simulation::moveParticle\n");
			//if (debugMode) alertNaNOrInfinity(particle->momentum.z, "p.z = NaN in simulation::moveParticle\n");
			//particle->momentum.x = -particle->momentum.x;
			return;

		}
		//if (tempParticle.coordinates.x < xgrid[0]) {
		if (tempParticle.coordinates.x < xgrid[1]) {
			//printf("particle number %d escaped to left\n", particle->number);
			escapedParticlesLeft.push_back(particle);
			particle->coordinates.x = tempParticle.coordinates.x + middleVelocity.x * restEtaDeltaT;

			particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * restEtaDeltaT;
			if (particle->coordinates.y > ygrid[ynumber]) {
				particle->coordinates.y -= ysize;
			}
			if (particle->coordinates.y < ygrid[0]) {
				particle->coordinates.y += ysize;
			}

			particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * restEtaDeltaT;
			if (particle->coordinates.z > zgrid[znumber]) {
				particle->coordinates.z -= zsize;
			}
			if (particle->coordinates.z < zgrid[0]) {
				particle->coordinates.z += zsize;
			}
			particle->setMomentumByV(newVelocity, speed_of_light_normalized);
			particle->escaped = true;
			particle->crossBoundaryCount++;
			//if (debugMode) alertNaNOrInfinity(particle->momentum.x, "p.x = NaN in simulation::moveParticle\n");
			//if (debugMode) alertNaNOrInfinity(particle->momentum.y, "p.y = NaN in simulation::moveParticle\n");
			//if (debugMode) alertNaNOrInfinity(particle->momentum.z, "p.z = NaN in simulation::moveParticle\n");
			return;
		}

		//if (tempParticle.coordinates.x > xgrid[xnumber+1]) {
		if (tempParticle.coordinates.x > xgrid[xnumber]) {
			escapedParticlesRight.push_back(particle);
			particle->coordinates.x = tempParticle.coordinates.x + middleVelocity.x * restEtaDeltaT;

			particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * restEtaDeltaT;
			if (particle->coordinates.y > ygrid[ynumber]) {
				particle->coordinates.y -= ysize;
			}
			if (particle->coordinates.y < ygrid[0]) {
				particle->coordinates.y += ysize;
			}
			particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * restEtaDeltaT;
			if (particle->coordinates.z > zgrid[znumber]) {
				particle->coordinates.z -= zsize;
			}
			if (particle->coordinates.z < zgrid[0]) {
				particle->coordinates.z += zsize;
			}
			//if (debugMode) alertNaNOrInfinity(particle->momentum.x, "p.x = NaN in simulation::moveParticle\n");
			//if (debugMode) alertNaNOrInfinity(particle->momentum.y, "p.y = NaN in simulation::moveParticle\n");
			//if (debugMode) alertNaNOrInfinity(particle->momentum.z, "p.z = NaN in simulation::moveParticle\n");
			//particle->setMomentumByV(newVelocity, speed_of_light_normalized);
			particle->escaped = true;
			particle->crossBoundaryCount++;
			return;
		}
		//}
		correctParticlePosition(tempParticle);
		updateCorrelationMaps(tempParticle);
		//checkParticleInBox(tempParticle);

		E = correlationTempEfield(tempParticle);
		//E = correlationNewEfield(tempParticle) * fieldScale;
		B = correlationBfield(tempParticle);
		//printf("E = %g %g %g\n", E.x, E.y, E.z);
		//printf("B = %g %g %g\n", B.x, B.y, B.z);

		tempParticle.addMomentum((E + (middleVelocity.vectorMult(
			B) / speed_of_light_normalized)) * (particle->charge * deltaT));
		//if (debugMode) alertNaNOrInfinity(tempParticle.momentum.x, "p.x = NaN in simulation::moveParticle\n");
		//if (debugMode) alertNaNOrInfinity(tempParticle.momentum.y, "p.y = NaN in simulation::moveParticle\n");
		//if (debugMode) alertNaNOrInfinity(tempParticle.momentum.z, "p.z = NaN in simulation::moveParticle\n");
		newVelocity = tempParticle.getVelocity(speed_of_light_normalized);
		error = (prevVelocity - newVelocity).norm();
	}

	//Vector3d prevCoordinates = particle->coordinates;

	particle->copyMomentum(tempParticle);
	//if (debugMode) alertNaNOrInfinity(particle->momentum.x, "p.x = NaN in simulation::moveParticle\n");
	//if (debugMode) alertNaNOrInfinity(particle->momentum.y, "p.y = NaN in simulation::moveParticle\n");
	//if (debugMode) alertNaNOrInfinity(particle->momentum.z, "p.z = NaN in simulation::moveParticle\n");

	/*if (particle->momentum.x > 1E100) {
		printf("particle->momentum.x > 1E100\n");
	}*/
	//particle->momentum.x = 0;

	particle->coordinates.x += middleVelocity.x * deltaT;

	particle->coordinates.y += middleVelocity.y * deltaT;
	particle->coordinates.z += middleVelocity.z * deltaT;

	correctParticlePosition(particle);
	//updateCorrelationMaps(particle);

	if (particle->coordinates.x < xgrid[1]) {
		if (boundaryConditionType == SUPER_CONDUCTOR_LEFT && rank == 0) {
			particle->coordinates.x = 2 * xgrid[1] - particle->coordinates.x;
			particle->reflectMomentumX();
			return;
		} else {
			//printf("particle number %d escaped to left\n", particle->number);
			particle->escaped = true;
			particle->crossBoundaryCount++;
			escapedParticlesLeft.push_back(particle);
			//particle->coordinates.x = xsize;
			return;
		}
	}
	if (particle->coordinates.x > xgrid[xnumber]) {
		//printf("particle number %d escaped to right\n", particle->number);
		//printLog("particle escaped to right\n");
		escapedParticlesRight.push_back(particle);
		particle->crossBoundaryCount++;
		//particle->coordinates.x = xsize;
		particle->escaped = true;
		return;
	}
}

void Simulation::correctParticlePosition(Particle* particle) {
	if (particle->coordinates.y < ygrid[0]) {
		particle->coordinates.y = particle->coordinates.y + ysize;
	}
	if (particle->coordinates.y > ygrid[ynumber]) {
		particle->coordinates.y = particle->coordinates.y - ysize;
	}

	if (particle->coordinates.z < zgrid[0]) {
		particle->coordinates.z = particle->coordinates.z + zsize;
	}
	if (particle->coordinates.z > zgrid[znumber]) {
		particle->coordinates.z = particle->coordinates.z - zsize;
	}


	/*if (particle->coordinates.x < xgrid[1]) {
	    if (boundaryConditionType == SUPER_CONDUCTOR_LEFT && rank == 0) {
	        particle->coordinates.x = 2 * xgrid[1] - particle->coordinates.x;
	        particle->momentum.x = -particle->momentum.x;
	        return;
	    } else {
	        escapedParticlesLeft.push_back(particle);
	        particle->coordinates.x = xsize;
	        particle->escaped = true;
	        return;
	    }
	}
	if (particle->coordinates.x > xgrid[xnumber]) {
	    escapedParticlesRight.push_back(particle);
	    particle->coordinates.x = xsize;
	    particle->escaped = true;
	    return;
	}*/
	/*if (boundaryConditionType == PERIODIC) {
	    if (particle.coordinates.x < xgrid[0]) {
	        particle.coordinates.x = particle.coordinates.x + xsize;
	    }
	    if (particle.coordinates.x > xgrid[xnumber]) {
	        particle.coordinates.x = particle.coordinates.x - xsize;
	    }
	}*/
}

void Simulation::correctParticlePosition(Particle& particle) {
	if (particle.coordinates.y < ygrid[0]) {
		particle.coordinates.y = particle.coordinates.y + ysize;
	}
	if (particle.coordinates.y > ygrid[ynumber]) {
		particle.coordinates.y = particle.coordinates.y - ysize;
	}

	if (particle.coordinates.z < zgrid[0]) {
		particle.coordinates.z = particle.coordinates.z + zsize;
	}
	if (particle.coordinates.z > zgrid[znumber]) {
		particle.coordinates.z = particle.coordinates.z - zsize;
	}
	/*if (boundaryConditionType == PERIODIC) {
	    if (particle.coordinates.x < xgrid[0]) {
	        particle.coordinates.x = particle.coordinates.x + xsize;
	    }
	    if (particle.coordinates.x > xgrid[xnumber]) {
	        particle.coordinates.x = particle.coordinates.x - xsize;
	    }
	}*/
}

void Simulation::evaluateParticlesRotationTensor() {
	double procTime = 0;
	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	for (int i = 0; i < particles.size(); ++i) {
		Particle* particle = particles[i];
		double beta = 0.5 * particle->charge * deltaT / particle->mass;
		double gamma = particle->gammaFactor(speed_of_light_normalized);
		Vector3d velocity = particle->getVelocity(speed_of_light_normalized);

		Vector3d oldE = correlationEfield(particle);
		Vector3d oldB = correlationBfield(particle);

		particle->rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);

	}
	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("evaluating ParticlesRotationTensor time = %g sec\n", procTime/CLOCKS_PER_SEC);
	}
}

Matrix3d Simulation::evaluateAlphaRotationTensor(double beta, Vector3d& velocity, double& gamma, Vector3d& EField, Vector3d& BField) {
	Matrix3d result;

	double G = ((beta * (EField.scalarMult(velocity)) / speed_of_light_normalized_sqr) + gamma);
	beta = beta / G;
	double denominator = G * (1 + beta * beta * BField.scalarMult(BField) / speed_of_light_normalized_sqr);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			result.matrix[i][j] = Kronecker.matrix[i][j] + (beta * beta * BField[i] * BField[j] / speed_of_light_normalized_sqr);
			//for (int k = 0; k < 3; ++k) {
				for (int l = 0; l < 3; ++l) {
					if(LeviCivita[j][i][l] != 0){
						result.matrix[i][j] -= (beta * LeviCivita[j][i][l] * BField[l] / speed_of_light_normalized);
					//result.matrix[i][j] += (beta * LeviCivita[j][k][l] * Kronecker.matrix[i][k] * BField[l] / speed_of_light_normalized);
					}
				}
			//}

			result.matrix[i][j] /= denominator;
			if (debugMode) alertNaNOrInfinity(result.matrix[i][j], "rotation tensor = NaN");
		}
	}

	return result;
}

void Simulation::injectNewParticles(int count, ParticleTypeContainer typeContainer, double length) {
	if ((rank == 0) && (verbosity > 0)) printf("inject new particles\n");

	//Particle* tempParticle = particles[0];

	double x = xgrid[xnumber] - length;
	double tempDeltaY = deltaY * uniformDistribution();
	double tempDeltaZ = deltaZ * uniformDistribution();

	/*if (typeContainer.type == ELECTRON && preserveChargeLocal) {
	    return;
	}*/

	for (int j = 0; j < ynumber; ++j) {
		for (int k = 0; k < znumber; ++k) {
			double weight = (typeContainer.concentration / typeContainer.particlesPerBin) * volumeB(xnumber - 1, j,
			                                                                                        k);
			for (int l = 0; l < count; ++l) {
				ParticleTypes type = typeContainer.type;
                if(verbosity > 1){
                    printf("inject particle number = %d\n", particlesNumber);
                }
				Particle* particle = createParticle(particlesNumber, xnumber - 1, j, k, weight, type, typeContainer,
				                                    typeContainer.temperatureX, typeContainer.temperatureY,
				                                    typeContainer.temperatureZ);
				particlesNumber++;
				particle->coordinates.x = x;
				//particle->coordinates.y = ygrid[j] + tempDeltaY;
				//particle->coordinates.z = zgrid[k] + tempDeltaZ;
				particle->coordinates.y = middleYgrid[j];
				particle->coordinates.z = middleZgrid[k];
				double y = particle->coordinates.y;
				double z = particle->coordinates.z;

				particle->addVelocity(V0, speed_of_light_normalized);
				Vector3d momentum = particle->getMomentum();
				particle->initialMomentum = momentum;
				//particle->prevMomentum = momentum;
				particles.push_back(particle);
				double en = particle->energy(speed_of_light_normalized) * particle->weight * sqr(
					scaleFactor / plasma_period);
				theoreticalEnergy += particle->energy(speed_of_light_normalized) * particle->weight * sqr(
					scaleFactor / plasma_period);
				theoreticalMomentum += particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
				/*if (preserveChargeLocal) {
				    int necessaryElectrons = 1;
				    if (type == ALPHA) {
				        necessaryElectrons = 2;
				    }
				    while (necessaryElectrons > 0) {
				        particle = createParticle(n, xnumber - 1, j, k, weight, ELECTRON, types[0],
				                                  types[0].temperatureX, types[0].temperatureY,
				                                  types[0].temperatureZ);
				        n++;
				        particle->coordinates.x = x;
				        particle->coordinates.y = y;
				        particle->coordinates.z = z;
				        particle->addVelocity(V0, speed_of_light_normalized);
				        particle->initialMomentum = particle->momentum;
				        particles.push_back(particle);
				        en = particle->energy(speed_of_light_normalized) * particle->weight * sqr(
				            scaleFactor / plasma_period);
				        theoreticalEnergy += particle->energy(speed_of_light_normalized) * particle->weight * sqr(
				            scaleFactor / plasma_period);
				        theoreticalMomentum += particle->momentum * particle->weight * scaleFactor / plasma_period;
				        necessaryElectrons--;
				    }
				}*/
			}
		}
	}
}

void Simulation::exchangeParticles() {
	double procTime = 0;
	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if (nprocs == 1) {
		if (boundaryConditionType == PERIODIC) {
			for (int i = 0; i < escapedParticlesLeft.size(); ++i) {
				Particle* particle = escapedParticlesLeft[i];
				particle->coordinates.x += xsizeGeneral;
				particle->escaped = false;
				particles.push_back(particle);
			}
			for (int i = 0; i < escapedParticlesRight.size(); ++i) {
				Particle* particle = escapedParticlesRight[i];
				particle->coordinates.x -= xsizeGeneral;
				particle->escaped = false;
				particles.push_back(particle);
			}
		}
	} else {
		if (rank == 0 && boundaryConditionType == PERIODIC) {
			for (int i = 0; i < escapedParticlesLeft.size(); ++i) {
				Particle* particle = escapedParticlesLeft[i];
				particle->coordinates.x += xsizeGeneral;
			}
		}
		if (rank == nprocs - 1 && boundaryConditionType == PERIODIC) {
			for (int i = 0; i < escapedParticlesRight.size(); ++i) {
				Particle* particle = escapedParticlesRight[i];
				particle->coordinates.x -= xsizeGeneral;
			}
		}
		if(verbosity > 2) printf("send particles left rank = %d\n", rank);
        sendLeftReceiveRightParticles(escapedParticlesLeft, particles, types, typesNumber,
                                      boundaryConditionType == PERIODIC, verbosity);
        MPI_Barrier(MPI_COMM_WORLD);
		if(verbosity > 2) printf("send particles right rank = %d\n", rank);
        sendRightReceiveLeftParticles(escapedParticlesRight, particles, types, typesNumber,
                                      boundaryConditionType == PERIODIC, verbosity);
	}
	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("exchange particles time = %g sec\n", procTime/CLOCKS_PER_SEC);
	}
}
