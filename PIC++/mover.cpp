#include "stdlib.h"
#include "stdio.h"
#include "cmath"
#include <omp.h>
#include <mpi.h>
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
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	//MPI_Barrier(cartComm);
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
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("moving particles time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::sortParticleToEscaped(Particle* particle) {
	if (particle->coordinates.x < xgrid[1 + additionalBinNumber]) {
		escapedParticlesLeft.push_back(particle);
		particle->escaped = true;
		particle->crossBoundaryCount++;
		return;
	}
	if (particle->coordinates.x > xgrid[xnumberAdded - 1 - additionalBinNumber]) {
		escapedParticlesRight.push_back(particle);
		particle->escaped = true;
		particle->crossBoundaryCount++;
		return;
	}

	if (particle->coordinates.y < ygrid[1 + additionalBinNumber]) {
		escapedParticlesFront.push_back(particle);
		particle->escaped = true;
		particle->crossBoundaryCount++;
		return;
	}
	if (particle->coordinates.y > ygrid[ynumberAdded - 1 - additionalBinNumber]) {
		escapedParticlesBack.push_back(particle);
		particle->escaped = true;
		particle->crossBoundaryCount++;
		return;
	}

	if (particle->coordinates.z < zgrid[1 + additionalBinNumber]) {
		escapedParticlesBottom.push_back(particle);
		particle->escaped = true;
		particle->crossBoundaryCount++;
		return;
	}
	if (particle->coordinates.z > zgrid[znumberAdded - 1 - additionalBinNumber]) {
		escapedParticlesTop.push_back(particle);
		particle->escaped = true;
		particle->crossBoundaryCount++;
	}
}

void Simulation::removeEscapedParticles() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((cartDim[0] > 1) || (boundaryConditionType == FREE_BOTH)) {
		for (int i = 0; i < escapedParticlesLeft.size(); ++i) {
			Particle* particle = escapedParticlesLeft[i];
			delete particle;
		}
	}
	escapedParticlesLeft.clear();

	if ((cartDim[0] > 1) || (boundaryConditionType != PERIODIC)) {
		for (int i = 0; i < escapedParticlesRight.size(); ++i) {
			Particle* particle = escapedParticlesRight[i];
			delete particle;
		}
	}
	escapedParticlesRight.clear();

	if (cartDim[1] > 1) {
		for (int i = 0; i < escapedParticlesFront.size(); ++i) {
			Particle* particle = escapedParticlesFront[i];
			delete particle;
		}
	}
	escapedParticlesFront.clear();

	if (cartDim[1] > 1) {
		for (int i = 0; i < escapedParticlesBack.size(); ++i) {
			Particle* particle = escapedParticlesBack[i];
			delete particle;
		}
	}
	escapedParticlesBack.clear();

	if (cartDim[2] > 1) {
		for (int i = 0; i < escapedParticlesBottom.size(); ++i) {
			Particle* particle = escapedParticlesBottom[i];
			delete particle;
		}
	}
	escapedParticlesBottom.clear();

	if (cartDim[2] > 1) {
		for (int i = 0; i < escapedParticlesTop.size(); ++i) {
			Particle* particle = escapedParticlesTop[i];
			delete particle;
		}
	}
	escapedParticlesTop.clear();

	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("removing escaped particles time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::eraseEscapedPaticles() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	std::vector<Particle*> tempParticles;
	if (particles.size() > 0) {
		tempParticles.reserve(particles.size());
		for (int i = 0; i < particles.size(); ++i) {
			Particle* particle = particles[i];
			if (particle->escaped) {
				chargeBalance -= particle->chargeCount;
			} else {
				tempParticles.push_back(particle);
			}
		}
	}
	particles.clear();
	particles = tempParticles;
	tempParticles.clear();
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("erasing escaped particles time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::moveParticle(Particle* particle) {
	updateCorrelationMaps(particle);
	particle->prevMomentum = particle->getMomentum();
	Vector3d E;
	Vector3d B;

	if(solverType == BUNEMAN){
		E = (correlationBunemanEfield(particle) + correlationBunemanNewEfield(particle))/2.0;
		B = (correlationBunemanBfield(particle) + correlationBunemanNewBfield(particle))/2.0;
	} else {
		E = correlationTempEfield(particle);
		//E = correlationNewEfield(particle);
		//B = correlationBfield(particle)*(1-theta) + correlationNewBfield(particle)*theta;
		B = correlationBfield(particle);
	}
	//B = B0;
	//printf("E = %g %g %g\n", E.x, E.y, E.z);
	//printf("B = %g %g %g\n", B.x, B.y, B.z);

	Vector3d velocity = particle->getVelocity(speed_of_light_normalized);
	Vector3d newVelocity = velocity;
	Vector3d middleVelocity = velocity;


	//see Noguchi
	double gamma = particle->gammaFactor(speed_of_light_normalized);
	double beta = 0.5 * particle->charge * deltaT / particle->mass;

	int particleIterations = 50;


	Particle tempParticle = *particle;

	//tempParticle.addMomentum((E + (velocity.vectorMult(B) / speed_of_light_normalized)) * particle->charge * deltaT);

	//boris
	bool boris = false;
	if (boris) {
		moveParticleBoris(particle);
		return;
	}
	//
	tempParticle.addMomentum((E + (velocity.vectorMult(B) / speed_of_light_normalized)) * particle->charge * deltaT);
	alertNaNOrInfinity(E.x, "E.x = Nan in move particle\n");

	newVelocity = tempParticle.getVelocity(speed_of_light_normalized);

	middleVelocity = velocity * (1 - eta) + newVelocity * eta;

	tempParticle.coordinates.x += middleVelocity.x * eta * deltaT;
	tempParticle.coordinates.y += middleVelocity.y * eta * deltaT;
	tempParticle.coordinates.z += middleVelocity.z * eta * deltaT;

	if ((tempParticle.coordinates.x < xgrid[1 + additionalBinNumber]) && (boundaryConditionType == SUPER_CONDUCTOR_LEFT) && (cartCoord[0] == 0)) {
		particle->coordinates.x = 2 * xgrid[1 + additionalBinNumber] - tempParticle.coordinates.x + fabs(
			middleVelocity.x * (1 - eta) * deltaT);
		particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * (1 - eta) * deltaT;
		particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * (1 - eta) * deltaT;
		newVelocity.x = fabs(newVelocity.x);
		particle->setMomentumByV(newVelocity, speed_of_light_normalized);
		theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
		sortParticleToEscaped(particle);
		return;
	}


	Vector3d prevVelocity = velocity;
	int i = 0;
	Vector3d velocityHat = (particle->rotationTensor * particle->gammaFactor(
		speed_of_light_normalized) * velocity);
	Vector3d Eperp = E - velocity*(velocity.scalarMult(E)/(velocity.scalarMult(velocity)));
	Vector3d electricVelocityShift = (Eperp * (2 * eta * beta / gamma));
	//velocityHat += electricVelocityShift;

	if (velocityHat.norm() > speed_of_light_normalized) {
		//printf("velocity Hat norm > c\n");
		//MPI_Finalize();
		//exit(0);
	}
	double a = tempParticle.gammaFactor(speed_of_light_normalized);
	double etaDeltaT = eta * deltaT;
	double restEtaDeltaT = (1.0 - eta) * deltaT;
	double velocityNorm = velocity.norm();

	Vector3d prevMomentum = particle->getMomentum();
	Vector3d momentum = tempParticle.getMomentum();
	Vector3d newMomentum = tempParticle.getMomentum();
	double momentumNorm = momentum.norm();

	//double error = (prevVelocity - newVelocity).norm();
	double error = (prevMomentum - newMomentum).norm();
	//error = 0;


	while (error > particleVelocityErrorLevel * momentumNorm && i < particleIterations) {
		++i;
		prevVelocity = newVelocity;
		prevMomentum = newMomentum;

		tempParticle = *particle;
		Vector3d rotatedE = particle->rotationTensor * E;

		middleVelocity = velocityHat + rotatedE * beta;

		tempParticle.coordinates.x += (middleVelocity.x * etaDeltaT);
		tempParticle.coordinates.y += (middleVelocity.y * etaDeltaT);
		tempParticle.coordinates.z += (middleVelocity.z * etaDeltaT);

		if ((tempParticle.coordinates.x < xgrid[1 + additionalBinNumber]) && (boundaryConditionType == SUPER_CONDUCTOR_LEFT) && (cartCoord[0] == 0)) {
			particle->coordinates.x = 2 * xgrid[1 + additionalBinNumber] - tempParticle.coordinates.x + fabs(
				middleVelocity.x * restEtaDeltaT);
			particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * restEtaDeltaT;
			particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * restEtaDeltaT;

			newVelocity.x = fabs(newVelocity.x);
			particle->setMomentumByV(newVelocity, speed_of_light_normalized);
			theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
			sortParticleToEscaped(particle);
			return;

		}


		//correctParticlePosition(tempParticle);
		updateCorrelationMaps(tempParticle);

		if(solverType == BUNEMAN){
			E = (correlationBunemanEfield(tempParticle) + correlationBunemanNewEfield(tempParticle))/2.0;
			B = (correlationBunemanBfield(tempParticle) + correlationBunemanNewBfield(tempParticle))/2.0;
		} else {
			E = correlationTempEfield(tempParticle);
			//E = correlationNewEfield(tempParticle);
			//B = correlationBfield(tempParticle)*(1-theta) + correlationNewBfield(tempParticle)*theta;
			B = correlationBfield(tempParticle);
		}
		//E = E0;
		//B = B0;

		tempParticle.addMomentum((E + (middleVelocity.vectorMult(B) / speed_of_light_normalized)) * (particle->charge * deltaT));
		newVelocity = tempParticle.getVelocity(speed_of_light_normalized);
		newMomentum = tempParticle.getMomentum();
		//error = (prevVelocity - newVelocity).norm();
		error = (prevMomentum - newMomentum).norm();
	}

	particle->copyMomentum(tempParticle);

	particle->coordinates.x += middleVelocity.x * deltaT;

	particle->coordinates.y += middleVelocity.y * deltaT;
	particle->coordinates.z += middleVelocity.z * deltaT;

	//correctParticlePosition(particle);
	/*if(particle->coordinates.x > xgrid[xnumberAdded - additionalBinNumber]){
		printf("aaa\n");
	}*/

	if (particle->coordinates.x < xgrid[1 + additionalBinNumber]) {
		if (boundaryConditionType == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
			particle->coordinates.x = 2 * xgrid[1 + additionalBinNumber] - particle->coordinates.x;
			particle->reflectMomentumX();
			theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
			//return;
		} else {
			particle->escaped = true;
			particle->crossBoundaryCount++;
			escapedParticlesLeft.push_back(particle);
			return;
		}
	}
	if (particle->coordinates.x > xgrid[xnumberAdded - 1 - additionalBinNumber]) {
		escapedParticlesRight.push_back(particle);
		particle->crossBoundaryCount++;
		particle->escaped = true;
		return;
	}

	if (particle->coordinates.y < ygrid[1 + additionalBinNumber]) {
		escapedParticlesFront.push_back(particle);
		particle->crossBoundaryCount++;
		particle->escaped = true;
		return;
	}

	if (particle->coordinates.y > ygrid[ynumberAdded - 1 - additionalBinNumber]) {
		escapedParticlesBack.push_back(particle);
		particle->crossBoundaryCount++;
		particle->escaped = true;
		return;
	}

	if (particle->coordinates.z < zgrid[1 + additionalBinNumber]) {
		escapedParticlesBottom.push_back(particle);
		particle->crossBoundaryCount++;
		particle->escaped = true;
		return;
	}

	if (particle->coordinates.z > zgrid[znumberAdded - 1 - additionalBinNumber]) {
		escapedParticlesTop.push_back(particle);
		particle->crossBoundaryCount++;
		particle->escaped = true;
	}
}

void Simulation::moveParticleBoris(Particle* particle) {
	updateCorrelationMaps(particle);

	Vector3d E1;
	Vector3d E2;
	Vector3d B;

	if(solverType == BUNEMAN){
		E1 = correlationBunemanEfield(particle);// + correlationBunemanNewEfield(particle);
		E2 = correlationBunemanNewEfield(particle);
		B = correlationBunemanBfield(particle);// + correlationBunemanNewBfield(particle);
	} else {
		E1 = correlationEfield(particle);
		E2 = correlationNewEfield(particle);
		//E = correlationNewEfield(particle) * fieldScale;
		B = correlationBfield(particle)*(1-theta) + correlationNewBfield(particle)*theta;
		//B = correlationBfield(particle);
	}
	//B = B0;
	//printf("E = %g %g %g\n", E.x, E.y, E.z);
	//printf("B = %g %g %g\n", B.x, B.y, B.z);


	double Bnorm = B.norm();
	double Bnorm2 = B.scalarMult(B);
	double beta = particle->charge *deltaT/2.0;

	Vector3d p1 = particle->getMomentum() + E1*beta;
	double tempGamma = sqrt(1.0 + (p1.scalarMult(p1)/(particle->mass*particle->mass*speed_of_light_normalized_sqr)));
	double omega = particle->charge*Bnorm/(particle->mass*speed_of_light*tempGamma);
	double a1 = tan(omega*deltaT/2.0)/Bnorm;
	if(Bnorm <= 0){
		a1 = particle->charge/(2.0*particle->mass*speed_of_light*tempGamma);
	}
	//double a1 = particle->charge/(2.0*particle->mass*speed_of_light*tempGamma);
	double a2 = 2*a1/(1 + a1*a1*Bnorm2);
	Vector3d p3 = p1 + (p1.vectorMult(B))*a1;
	Vector3d p2 = p1 + (p3.vectorMult(B))*a2;

	Vector3d newMomentum = p2 + E2*beta;

	particle->setMomentum(newMomentum);
	Vector3d newVelocity = particle->getVelocity(speed_of_light_normalized);
	particle->coordinates += newVelocity*deltaT;

	if (particle->coordinates.x < xgrid[1 + additionalBinNumber]) {
		if (boundaryConditionType == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
			particle->coordinates.x = 2 * xgrid[1 + additionalBinNumber] - particle->coordinates.x;
			particle->reflectMomentumX();
			theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
			//return;
		} else {
			particle->escaped = true;
			particle->crossBoundaryCount++;
			escapedParticlesLeft.push_back(particle);
			return;
		}
	}
	if (particle->coordinates.x > xgrid[xnumberAdded - 1 - additionalBinNumber]) {
		escapedParticlesRight.push_back(particle);
		particle->crossBoundaryCount++;
		particle->escaped = true;
		return;
	}

	if (particle->coordinates.y < ygrid[1 + additionalBinNumber]) {
		escapedParticlesFront.push_back(particle);
		particle->crossBoundaryCount++;
		particle->escaped = true;
		return;
	}

	if (particle->coordinates.y > ygrid[ynumberAdded - 1 - additionalBinNumber]) {
		escapedParticlesBack.push_back(particle);
		particle->crossBoundaryCount++;
		particle->escaped = true;
		return;
	}

	if (particle->coordinates.z < zgrid[1 + additionalBinNumber]) {
		escapedParticlesBottom.push_back(particle);
		particle->crossBoundaryCount++;
		particle->escaped = true;
		return;
	}

	if (particle->coordinates.z > zgrid[znumberAdded - 1 - additionalBinNumber]) {
		escapedParticlesTop.push_back(particle);
		particle->crossBoundaryCount++;
		particle->escaped = true;
	}
}

void Simulation::evaluateParticlesRotationTensor() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	for (int i = 0; i < particles.size(); ++i) {
		Particle* particle = particles[i];
		double beta = 0.5 * particle->charge * deltaT / particle->mass;
		double gamma = particle->gammaFactor(speed_of_light_normalized);
		Vector3d velocity = particle->getVelocity(speed_of_light_normalized);

		Vector3d oldE;
		Vector3d oldB;

		if(solverType == BUNEMAN){
			oldE = correlationBunemanEfield(particle);
			oldB = correlationBunemanBfield(particle);
		} else {
			oldE = correlationEfield(particle);
			//E = correlationNewEfield(particle) * fieldScale;
			//B = correlationBfield(particle)*(1-theta) + correlationNewBfield(particle)*theta;
			oldB = correlationBfield(particle);
		}
		particle->rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);

	}
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("evaluating ParticlesRotationTensor time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

Matrix3d Simulation::evaluateAlphaRotationTensor(double beta, Vector3d& velocity, double& gamma, Vector3d& EField, Vector3d& BField) {
	Matrix3d result = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);

	double G = ((beta * (EField.scalarMult(velocity)) / speed_of_light_normalized_sqr) + gamma);
	beta = beta / G;
	double beta2c = beta * beta / speed_of_light_normalized_sqr;
	double denominator = G * (1 + beta2c * BField.scalarMult(BField));

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			result.matrix[i][j] = Kronecker.matrix[i][j] + (beta * beta * BField[i] * BField[j] / speed_of_light_normalized_sqr);
			//for (int k = 0; k < 3; ++k) {
			for (int l = 0; l < 3; ++l) {
				if (LeviCivita[j][i][l] != 0) {
					result.matrix[i][j] -= (beta * LeviCivita[j][i][l] * BField[l] / speed_of_light_normalized);
					//result.matrix[i][j] += (beta * LeviCivita[j][k][l] * Kronecker.matrix[i][k] * BField[l] / speed_of_light_normalized);
				}
			}
			//}

			result.matrix[i][j] /= denominator;
			if (debugMode) alertNaNOrInfinity(result.matrix[i][j], "rotation tensor = NaN");
		}
	}

	/*result.matrix[0][0] = 1.0 + beta2c*BField[0]*BField[0];
	result.matrix[0][1] = beta2c*BField[0]*BField[1];
	result.matrix[0][2] = beta2c*BField[0]*BField[2];
	result.matrix[1][1] = 1.0 + beta2c*BField[1]*BField[1];
	result.matrix[1][2] = beta2c*BField[1]*BField[2];
	result.matrix[2][2] = 1.0 + beta2c*BField[2]*BField[2];
	result.matrix[1][0] = result.matrix[0][1];
	result.matrix[2][0] = result.matrix[0][2];
	result.matrix[2][1] = result.matrix[1][2];

	result.matrix[0][1] += beta * BField[2] / speed_of_light_normalized;
	result.matrix[0][2] -= beta * BField[1] / speed_of_light_normalized;
	result.matrix[1][0] -= beta * BField[2] / speed_of_light_normalized;
	result.matrix[1][2] += beta * BField[0] / speed_of_light_normalized;
	result.matrix[2][0] += beta * BField[1] / speed_of_light_normalized;
	result.matrix[2][1] -= beta * BField[0] / speed_of_light_normalized;

	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 3; ++j) {
			result.matrix[i][j] /= denominator;
			if (debugMode) alertNaNOrInfinity(result.matrix[i][j], "rotation tensor = NaN");
		}
	}*/

	return result;
}

void Simulation::injectNewParticles(int count, ParticleTypeContainer typeContainer, double length) {
	if ((rank == 0) && (verbosity > 0)) printf("inject new particles\n");

	//Particle* tempParticle = particles[0];

	double x = xgrid[xnumberAdded - 1 - additionalBinNumber] - length;
	double tempDeltaY = deltaY * uniformDistribution();
	double tempDeltaZ = deltaZ * uniformDistribution();

	/*if (typeContainer.type == ELECTRON && preserveChargeLocal) {
	    return;
	}*/

	for (int j = 1 + additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
		for (int k = 1 + additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
			double weight = (typeContainer.concentration / typeContainer.particlesPerBin) * volumeB();
			for (int l = 0; l < count; ++l) {
				ParticleTypes type = typeContainer.type;
				if (verbosity > 1) {
					printf("inject particle number = %d\n", particlesNumber);
				}
				Particle* particle = createParticle(particlesNumber, xnumberAdded - 1 - additionalBinNumber, j, k, weight, type, typeContainer,
				                                    typeContainer.temperatureX, typeContainer.temperatureY,
				                                    typeContainer.temperatureZ);
				particlesNumber++;
				particle->coordinates.x = x;
				particle->coordinates.y = ygrid[j] + tempDeltaY;
				particle->coordinates.z = zgrid[k] + tempDeltaZ;
				//particle->coordinates.y = middleYgrid[j];
				//particle->coordinates.z = middleZgrid[k];
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
			}
		}
	}
}

void Simulation::exchangeParticles() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if (cartDim[0] == 1) {
		if (boundaryConditionType == PERIODIC) {
			for (int i = 0; i < escapedParticlesLeft.size(); ++i) {
				Particle* particle = escapedParticlesLeft[i];
				particle->coordinates.x += xsizeGeneral;
				particle->escaped = false;
				tempParticles.push_back(particle);
			}
			for (int i = 0; i < escapedParticlesRight.size(); ++i) {
				Particle* particle = escapedParticlesRight[i];
				particle->coordinates.x -= xsizeGeneral;
				particle->escaped = false;
				tempParticles.push_back(particle);
			}
		}
	} else {
		if (cartCoord[0] == 0 && boundaryConditionType == PERIODIC) {
			for (int i = 0; i < escapedParticlesLeft.size(); ++i) {
				Particle* particle = escapedParticlesLeft[i];
				particle->coordinates.x += xsizeGeneral;
			}
		}
		if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionType == PERIODIC) {
			for (int i = 0; i < escapedParticlesRight.size(); ++i) {
				Particle* particle = escapedParticlesRight[i];
				particle->coordinates.x -= xsizeGeneral;
			}
		}
		if (verbosity > 2) printf("send particles left rank = %d\n", rank);
		sendLeftReceiveRightParticles(escapedParticlesLeft, tempParticles, types, typesNumber,
		                              boundaryConditionType == PERIODIC, verbosity, cartComm, rank, leftRank, rightRank);
		MPI_Barrier(cartComm);
		if (verbosity > 2) printf("send particles right rank = %d\n", rank);
		sendRightReceiveLeftParticles(escapedParticlesRight, tempParticles, types, typesNumber,
		                              boundaryConditionType == PERIODIC, verbosity, cartComm, rank, leftRank, rightRank);
	}

	for (int pcount = 0; pcount < tempParticles.size(); ++pcount) {
		Particle* particle = tempParticles[pcount];
		if (particle->coordinates.y < ygrid[1 + additionalBinNumber]) {
			escapedParticlesFront.push_back(particle);
			particle->escaped = true;
			particle->crossBoundaryCount++;
		} else if (particle->coordinates.y > ygrid[ynumberAdded - 1 - additionalBinNumber]) {
			escapedParticlesBack.push_back(particle);
			particle->escaped = true;
			particle->crossBoundaryCount++;
		} else if (particle->coordinates.z < zgrid[1 + additionalBinNumber]) {
			escapedParticlesBottom.push_back(particle);
			particle->escaped = true;
			particle->crossBoundaryCount++;
		} else if (particle->coordinates.z > zgrid[znumberAdded - 1 - additionalBinNumber]) {
			escapedParticlesTop.push_back(particle);
			particle->escaped = true;
			particle->crossBoundaryCount++;
		} else {
			particles.push_back(particle);
		}
	}
	tempParticles.clear();

	MPI_Barrier(cartComm);

	if (cartDim[1] == 1) {
		for (int i = 0; i < escapedParticlesFront.size(); ++i) {
			Particle* particle = escapedParticlesFront[i];
			particle->coordinates.y += ysizeGeneral;
			particle->escaped = false;
			tempParticles.push_back(particle);
		}
		for (int i = 0; i < escapedParticlesBack.size(); ++i) {
			Particle* particle = escapedParticlesBack[i];
			particle->coordinates.y -= ysizeGeneral;
			particle->escaped = false;
			tempParticles.push_back(particle);
		}

	} else {
		if (cartCoord[1] == 0) {
			for (int i = 0; i < escapedParticlesFront.size(); ++i) {
				Particle* particle = escapedParticlesFront[i];
				particle->coordinates.y += ysizeGeneral;
			}
		}
		if (cartCoord[1] == cartDim[1] - 1) {
			for (int i = 0; i < escapedParticlesBack.size(); ++i) {
				Particle* particle = escapedParticlesBack[i];
				particle->coordinates.y -= ysizeGeneral;
			}
		}
		if (verbosity > 2) printf("send particles front rank = %d\n", rank);
		sendFrontReceiveBackParticles(escapedParticlesFront, tempParticles, types, typesNumber,
		                              boundaryConditionType == PERIODIC, verbosity, cartComm, rank, frontRank, backRank);
		MPI_Barrier(cartComm);
		if (verbosity > 2) printf("send particles back rank = %d\n", rank);
		sendBackReceiveFrontParticles(escapedParticlesBack, tempParticles, types, typesNumber,
		                              boundaryConditionType == PERIODIC, verbosity, cartComm, rank, frontRank, backRank);
	}

	for (int pcount = 0; pcount < tempParticles.size(); ++pcount) {
		Particle* particle = tempParticles[pcount];
		if (particle->coordinates.z < zgrid[1 + additionalBinNumber]) {
			escapedParticlesBottom.push_back(particle);
			particle->escaped = true;
			particle->crossBoundaryCount++;
		} else if (particle->coordinates.z > zgrid[znumberAdded - 1 - additionalBinNumber]) {
			escapedParticlesTop.push_back(particle);
			particle->escaped = true;
			particle->crossBoundaryCount++;
		} else {
			particles.push_back(particle);
		}
	}
	tempParticles.clear();


	MPI_Barrier(cartComm);

	if (cartDim[2] == 1) {
		for (int i = 0; i < escapedParticlesBottom.size(); ++i) {
			Particle* particle = escapedParticlesBottom[i];
			particle->coordinates.z += zsizeGeneral;
			particle->escaped = false;
			particles.push_back(particle);
		}
		for (int i = 0; i < escapedParticlesTop.size(); ++i) {
			Particle* particle = escapedParticlesTop[i];
			particle->coordinates.z -= zsizeGeneral;
			particle->escaped = false;
			particles.push_back(particle);
		}

	} else {
		if (cartCoord[2] == 0) {
			for (int i = 0; i < escapedParticlesBottom.size(); ++i) {
				Particle* particle = escapedParticlesBottom[i];
				particle->coordinates.z += zsizeGeneral;
			}
		}
		if (cartCoord[2] == cartDim[2] - 1) {
			for (int i = 0; i < escapedParticlesTop.size(); ++i) {
				Particle* particle = escapedParticlesTop[i];
				particle->coordinates.z -= zsizeGeneral;
			}
		}
		if (verbosity > 2) printf("send particles bottom rank = %d\n", rank);
		sendBottomReceiveTopParticles(escapedParticlesBottom, particles, types, typesNumber,
		                              boundaryConditionType == PERIODIC, verbosity, cartComm, rank, bottomRank, topRank);
		MPI_Barrier(cartComm);
		if (verbosity > 2) printf("send particles top rank = %d\n", rank);
		sendTopReceiveBottomParticles(escapedParticlesTop, particles, types, typesNumber,
		                              boundaryConditionType == PERIODIC, verbosity, cartComm, rank, bottomRank, topRank);
	}

	MPI_Barrier(cartComm);
	tempParticles.clear();
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("exchange particles time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}
