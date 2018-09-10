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
#include "paths.h"

void Simulation::moveParticles() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	//MPI_Barrier(cartComm);
	if ((rank == 0) && (verbosity > 0)) printf("moving particles\n");
	if ((rank == 0) && (verbosity > 0)) printLog("moving particles\n");
	int a[10];
	if(solverType == BUNEMAN){
		for (int i = 0; i < particles.size(); ++i) {
			Particle* particle = particles[i];
		//for (auto particle : particles) {
			moveParticleTristan(particle);
		}
	} else {
		for (int i = 0; i < particles.size(); ++i) {
		//Particle* particle = particles[i];
		//for (auto particle : particles) {
			moveParticle(particles[i]);
		}
	}


	//if (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
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

	if (ynumberGeneral > 1) {
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
	}

	if (znumberGeneral > 1) {
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
}

void Simulation::removeEscapedParticles() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((cartDim[0] > 1) || (boundaryConditionTypeX == FREE_BOTH)) {
		for (int i = 0; i < escapedParticlesLeft.size(); ++i) {
			Particle* particle = escapedParticlesLeft[i];
			//reservedParticles.push_back(particle);
			delete particle;
		}
	}
	escapedParticlesLeft.clear();

	if ((cartDim[0] > 1) || (boundaryConditionTypeX != PERIODIC)) {
		for (int i = 0; i < escapedParticlesRight.size(); ++i) {
			Particle* particle = escapedParticlesRight[i];
			//reservedParticles.push_back(particle);
			delete particle;
		}
	}
	escapedParticlesRight.clear();

	if (cartDim[1] > 1) {
		for (int i = 0; i < escapedParticlesFront.size(); ++i) {
			Particle* particle = escapedParticlesFront[i];
			//reservedParticles.push_back(particle);
			delete particle;
		}
	}
	escapedParticlesFront.clear();

	if (cartDim[1] > 1) {
		for (int i = 0; i < escapedParticlesBack.size(); ++i) {
			Particle* particle = escapedParticlesBack[i];
			//reservedParticles.push_back(particle);
			delete particle;
		}
	}
	escapedParticlesBack.clear();

	if (cartDim[2] > 1) {
		for (int i = 0; i < escapedParticlesBottom.size(); ++i) {
			Particle* particle = escapedParticlesBottom[i];
			//reservedParticles.push_back(particle);
			delete particle;
		}
	}
	escapedParticlesBottom.clear();

	if (cartDim[2] > 1) {
		for (int i = 0; i < escapedParticlesTop.size(); ++i) {
			Particle* particle = escapedParticlesTop[i];
			//reservedParticles.push_back(particle);
			delete particle;
		}
	}
	escapedParticlesTop.clear();

	//MPI_Barrier(cartComm);
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
	std::vector < Particle* > tempParticles1;
	if (particles.size() > 0) {
		tempParticles1.reserve(particles.size());
		for (int i = 0; i < particles.size(); ++i) {
			Particle* particle = particles[i];
			if (particle->escaped) {
				chargeBalance -= particle->chargeCount;
			} else {
				tempParticles1.push_back(particle);
			}
		}
	}
	particles.clear();
	particles = tempParticles1;
	tempParticles1.clear();
	//MPI_Barrier(cartComm);
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

	double oldGamma = particle->gammaFactor();

	if (solverType == BUNEMAN) {
		E = (correlationBunemanEfield(particle) + correlationBunemanNewEfield(particle)) / 2.0;
		B = (correlationBunemanBfield(particle) + correlationBunemanNewBfield(particle)) / 2.0;
	} else {
		E = correlationTempEfield(particle);
		//E = correlationNewEfield(particle);
		//E = correlationEfield(particle);
		//B = correlationBfield(particle)*(1-theta) + correlationNewBfield(particle)*theta;
		B = correlationBfield(particle);
	}

	Vector3d velocity = particle->getVelocity();
	Vector3d newVelocity = velocity;
	Vector3d middleVelocity = velocity;


	//see Noguchi
	double gamma = particle->gammaFactor();
	double beta = 0.5 * particle->charge * deltaT / particle->mass;

	int particleIterations = 50 * gamma;
	particleIterations = 3;


	//Particle tempParticle = *particle;
	Vector3d oldMomentum = particle->getMomentum();
	Vector3d newMomentum = oldMomentum;
	Vector3d oldCoordinates = particle->coordinates;


	double etaDeltaT = eta * deltaT;
	double restEtaDeltaT = (1.0 - eta) * deltaT;

	//newMomentum += (E + (velocity.vectorMult(B) / (speed_of_light_normalized))) * particle->charge * deltaT;
	newMomentum += (E + (velocity.vectorMult(B))) * particle->charge * deltaT;

	newVelocity = Particle::evaluateVelocity(newMomentum, particle->mass);

	middleVelocity = velocity * (1 - eta) + newVelocity * eta;


	particle->coordinates.x += middleVelocity.x * etaDeltaT;
	particle->coordinates.y += middleVelocity.y * etaDeltaT;
	particle->coordinates.z += middleVelocity.z * etaDeltaT;


	if ((particle->coordinates.x < xgrid[1 + additionalBinNumber]) && (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT || boundaryConditionTypeX ==
		FREE_MIRROR_BOTH) && (cartCoord[0] == 0)) {
		particle->coordinates.x = 2 * xgrid[1 + additionalBinNumber] - particle->coordinates.x + fabs(
			middleVelocity.x * restEtaDeltaT);
		particle->coordinates.y = particle->coordinates.y + middleVelocity.y * restEtaDeltaT;
		particle->coordinates.z = particle->coordinates.z + middleVelocity.z * restEtaDeltaT;
		newVelocity.x = fabs(newVelocity.x);
		particle->setMomentumByV(newVelocity);
		theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
		sortParticleToEscaped(particle);
		return;
	}

	if ((particle->coordinates.x > xgrid[xnumberAdded - 1 - additionalBinNumber]) && (boundaryConditionTypeX == FREE_MIRROR_BOTH) && (cartCoord[0] ==
		cartDim[0] - 1)) {
		particle->coordinates.x = 2 * xgrid[xnumberAdded - 1 - additionalBinNumber] - particle->coordinates.x - fabs(
			middleVelocity.x * restEtaDeltaT);
		particle->coordinates.y = particle->coordinates.y + middleVelocity.y * restEtaDeltaT;
		particle->coordinates.z = particle->coordinates.z + middleVelocity.z * restEtaDeltaT;
		newVelocity.x = -fabs(newVelocity.x);
		particle->setMomentumByV(newVelocity);
		theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
		sortParticleToEscaped(particle);
		return;
	}


	//Vector3d prevVelocity = velocity;
	int i = 0;
	Vector3d velocityHat = (particle->rotationTensor * particle->gammaFactor(
	) * velocity);


	Vector3d prevMomentum = particle->getMomentum();
	Vector3d momentum = oldMomentum;
	double momentumNorm = momentum.norm();

	double error = (prevMomentum - newMomentum).norm();
	Vector3d prevCoordinates = particle->coordinates;
	Vector3d newCoordinates = particle->coordinates;
	double coordinatesError = deltaX;
	//error = 0;


	while (error > particleVelocityErrorLevel * momentumNorm && coordinatesError > particleCoordinatesErrorLevel * deltaX && i < particleIterations) {
		++i;
		//prevVelocity = newVelocity;
		prevMomentum = newMomentum;
		prevCoordinates = newCoordinates;

		//tempParticle = *particle;
		Vector3d rotatedE = particle->rotationTensor * E;

		double Enorm = E.norm();
		double rotatedEnorm = rotatedE.norm();

		middleVelocity = velocityHat + rotatedE * beta;


		particle->coordinates.x = oldCoordinates.x + (middleVelocity.x * etaDeltaT);
		particle->coordinates.y = oldCoordinates.y + (middleVelocity.y * etaDeltaT);
		particle->coordinates.z = oldCoordinates.z + (middleVelocity.z * etaDeltaT);

		//newCoordinates = tempParticle.coordinates;
		newCoordinates = particle->coordinates;

		//if ((tempParticle.coordinates.x < xgrid[1 + additionalBinNumber]) && (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) && (cartCoord[0] == 0)) {
		if ((particle->coordinates.x < xgrid[1 + additionalBinNumber]) && (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT || boundaryConditionTypeX ==
			FREE_MIRROR_BOTH) && (cartCoord[0] == 0)) {
			//particle->coordinates.x = 2 * xgrid[1 + additionalBinNumber] - tempParticle.coordinates.x + fabs(
			//	middleVelocity.x * restEtaDeltaT);
			//particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * restEtaDeltaT;
			//particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * restEtaDeltaT;

			particle->coordinates.x = 2 * xgrid[1 + additionalBinNumber] - particle->coordinates.x + fabs(
				middleVelocity.x * restEtaDeltaT);
			particle->coordinates.y = particle->coordinates.y + middleVelocity.y * restEtaDeltaT;
			particle->coordinates.z = particle->coordinates.z + middleVelocity.z * restEtaDeltaT;

			newVelocity.x = fabs(newVelocity.x);
			particle->setMomentumByV(newVelocity);
			theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
			sortParticleToEscaped(particle);
			return;

		}

		if ((particle->coordinates.x > xgrid[xnumberAdded - 1 - additionalBinNumber]) && (boundaryConditionTypeX == FREE_MIRROR_BOTH) && (cartCoord[0] ==
			cartDim[0] - 1)) {
			particle->coordinates.x = 2 * xgrid[xnumberAdded - 1 - additionalBinNumber] - particle->coordinates.x - fabs(
				middleVelocity.x * restEtaDeltaT);
			particle->coordinates.y = particle->coordinates.y + middleVelocity.y * restEtaDeltaT;
			particle->coordinates.z = particle->coordinates.z + middleVelocity.z * restEtaDeltaT;

			newVelocity.x = -fabs(newVelocity.x);
			particle->setMomentumByV(newVelocity);
			theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
			sortParticleToEscaped(particle);
			return;
		}


		//correctParticlePosition(tempParticle);
		//updateCorrelationMaps(tempParticle);
		updateCorrelationMaps(particle);

		if (solverType == BUNEMAN) {
			//E = (correlationBunemanEfield(tempParticle) + correlationBunemanNewEfield(tempParticle)) / 2.0;
			//B = (correlationBunemanBfield(tempParticle) + correlationBunemanNewBfield(tempParticle)) / 2.0;

			E = (correlationBunemanEfield(particle) + correlationBunemanNewEfield(particle)) / 2.0;
			B = (correlationBunemanBfield(particle) + correlationBunemanNewBfield(particle)) / 2.0;
		} else {
			//E = correlationTempEfield(tempParticle);
			//B = correlationBfield(tempParticle);

			E = correlationTempEfield(particle);
			B = correlationBfield(particle);
		}
		//E = E0;
		//B = B0;

		//tempParticle.addMomentum((E + (middleVelocity.vectorMult(B) / (speed_of_light_normalized))) * (particle->charge *deltaT));
		//newMomentum = tempParticle.getMomentum();

		//newMomentum = oldMomentum +(E + (middleVelocity.vectorMult(B) / (speed_of_light_normalized))) * (particle->charge * deltaT);
		newMomentum = oldMomentum +(E + middleVelocity.vectorMult(B)) * (particle->charge * deltaT);

		error = (prevMomentum - newMomentum).norm();
		coordinatesError = (prevCoordinates - newCoordinates).norm();
	}

	//particle->copyMomentum(tempParticle);
	particle->setMomentum(newMomentum);

	//particle->coordinates.x += middleVelocity.x * deltaT;
	//particle->coordinates.y += middleVelocity.y * deltaT;
	//particle->coordinates.z += middleVelocity.z * deltaT;

	particle->coordinates.x = oldCoordinates.x + middleVelocity.x * deltaT;
	particle->coordinates.y = oldCoordinates.y + middleVelocity.y * deltaT;
	particle->coordinates.z = oldCoordinates.z + middleVelocity.z * deltaT;

	//double newGamma = particle->gammaFactor(speed_of_light_normalized);

	//double deltaGammaTheor = particle->charge * deltaT * E.scalarMult(middleVelocity) / (particle->mass *speed_of_light_normalized_sqr);
	double deltaGammaTheor = particle->charge * deltaT * E.scalarMult(middleVelocity) / (particle->mass);
	//double theorNewGamma = oldGamma + deltaGammaTheor;
	/*if(theorNewGamma - 1 > relativisticPrecision){
		double newTheorMomentumNorm = particle->mass*speed_of_light_normalized*sqrt(theorNewGamma*theorNewGamma - 1.0);
		Vector3d p = particle->getMomentum();
		double newMomentumNorm = p.norm();
		p = p*(newTheorMomentumNorm/newMomentumNorm);
		particle->setMomentum(p);
	}*/

	/*if(i >= particleIterations){
		printf("i >= particle iterations\n");
		printf("error = %g\n", error);
		printf("relative error = %g\n", error/momentumNorm);
	}*/

	//if (fabs(newGamma - oldGamma) > 0.2) {
	//printf("delta gamma > 0.2\n");
	//printf("oldGamma = %g newGamma = %g delta gamma = %g\n", oldGamma, newGamma, newGamma - oldGamma);
	//printf("theoretical delta gamma = %g\n", deltaGammaTheor);
	/*printf("particle iterations = %d\n", i);
	printf("cartx = %d carty = %d cartz = %d\n", cartCoord[0], cartCoord[1], cartCoord[2]);
	printf("particle number = %d\n", particle->number);
	printf("x = %g y = %g z = %g\n", particle->coordinates.x, particle->coordinates.y, particle->coordinates.z);
	printf("xindex = %d yindex = %d zindex = %d\n", particle->correlationMapCell.xindex[(2 + splineOrder)/2], particle->correlationMapCell.yindex[(2 + splineOrder)/2], particle->correlationMapCell.zindex[(2 + splineOrder)/2]);
	printf("vx = %g vy = %g vz = %g\n", velocity.x, velocity.y, velocity.z);
	printf("newvx = %g newvy = %g newvz = %g\n", newVelocity.x, newVelocity.y, newVelocity.z);
	printf("Ex = %g Ey = %g Ez = %g\n", E.x, E.y, E.z);
	printf("Bx = %g By = %g Bz = %g\n", B.x, B.y, B.z);
	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			printf("alpha[%d][%d] = %g\n", i, j, particle->rotationTensor.matrix[i][j]);
		}
	}*/
	//}

	//correctParticlePosition(particle);
	/*if(particle->coordinates.x > xgrid[xnumberAdded - additionalBinNumber]){
		printf("aaa\n");
	}*/

	if (particle->coordinates.x < xgrid[1 + additionalBinNumber]) {
		if ((boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT || boundaryConditionTypeX == FREE_MIRROR_BOTH) && cartCoord[0] == 0) {
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
		if ((boundaryConditionTypeX == FREE_MIRROR_BOTH) && (cartCoord[0] == cartDim[0] - 1)) {
			particle->coordinates.x = 2 * xgrid[xnumberAdded - 1 - additionalBinNumber] - particle->coordinates.x;
			particle->reflectMomentumX();
			theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
		} else {
			escapedParticlesRight.push_back(particle);
			particle->crossBoundaryCount++;
			particle->escaped = true;
			return;
		}
	}

	if (ynumberGeneral > 1) {
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
	}

	if (znumberGeneral > 1) {
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
}

void Simulation::moveParticle(Particle* particle, int cur, int N) {
	double localTheta = cur * 1.0 / N;
	if (cur == N) {
		sortParticleToEscaped(particle);
		return;
	}
	updateCorrelationMaps(particle);
	particle->prevMomentum = particle->getMomentum();
	Vector3d E;
	Vector3d B;

	Vector3d oldE;
	Vector3d oldB;

	if (solverType == BUNEMAN) {
		E = (correlationBunemanEfield(particle) + correlationBunemanNewEfield(particle)) / 2.0;
		B = (correlationBunemanBfield(particle) + correlationBunemanNewBfield(particle)) / 2.0;
	} else {
		E = correlationEfield(particle) * (1 - localTheta - theta / N) + correlationNewEfield(particle) * (localTheta + theta
			/ N);
		//E = correlationNewEfield(particle);
		//B = correlationBfield(particle)*(1-theta) + correlationNewBfield(particle)*theta;
		B = correlationBfield(particle) * (1 - localTheta) + correlationNewBfield(particle) * localTheta;
	}
	//B = B0;
	//printf("E = %g %g %g\n", E.x, E.y, E.z);
	//printf("B = %g %g %g\n", B.x, B.y, B.z);

	Vector3d velocity = particle->getVelocity();
	Vector3d newVelocity = velocity;
	Vector3d middleVelocity = velocity;

	oldE = correlationEfield(particle) * (1 - localTheta - theta / N);
	oldB = correlationBfield(particle) * (1 - localTheta);

	//see Noguchi
	double gamma = particle->gammaFactor();
	double beta = 0.5 * particle->charge * deltaT / particle->mass;

	particle->rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);

	int particleIterations = 5;


	Particle tempParticle = *particle;

	//tempParticle.addMomentum((E + (velocity.vectorMult(B) / speed_of_light_normalized)) * particle->charge * deltaT);


	//tempParticle.addMomentum((E + (velocity.vectorMult(B) / speed_of_light_normalized)) * particle->charge * deltaT / N);
	tempParticle.addMomentum((E + velocity.vectorMult(B)) * particle->charge * deltaT / N);
	//alertNaNOrInfinity(E.x, "E.x = Nan in move particle\n");

	newVelocity = tempParticle.getVelocity();

	middleVelocity = velocity * (1 - eta) + newVelocity * eta;

	tempParticle.coordinates.x += middleVelocity.x * eta * deltaT;
	tempParticle.coordinates.y += middleVelocity.y * eta * deltaT;
	tempParticle.coordinates.z += middleVelocity.z * eta * deltaT;

	if ((tempParticle.coordinates.x < xgrid[1 + additionalBinNumber]) && (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT || boundaryConditionTypeX ==
			FREE_MIRROR_BOTH)
		&& (cartCoord[0] == 0)) {
		particle->coordinates.x = 2 * xgrid[1 + additionalBinNumber] - tempParticle.coordinates.x + fabs(
			middleVelocity.x * (1 - eta) * deltaT / N);
		particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * (1 - eta) * deltaT / N;
		particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * (1 - eta) * deltaT / N;
		newVelocity.x = fabs(newVelocity.x);
		particle->setMomentumByV(newVelocity);
		theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
		//sortParticleToEscaped(particle);
		moveParticle(particle, cur + 1, N);
		return;
	}

	if ((tempParticle.coordinates.x > xgrid[xnumberAdded - 1 - additionalBinNumber]) && (boundaryConditionTypeX == FREE_MIRROR_BOTH)
		&& (cartCoord[0] == cartDim[0] - 1)) {
		particle->coordinates.x = 2 * xgrid[xnumberAdded - 1 - additionalBinNumber] - tempParticle.coordinates.x - fabs(
			middleVelocity.x * (1 - eta) * deltaT / N);
		particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * (1 - eta) * deltaT / N;
		particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * (1 - eta) * deltaT / N;
		newVelocity.x = -fabs(newVelocity.x);
		particle->setMomentumByV(newVelocity);
		theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
		//sortParticleToEscaped(particle);
		moveParticle(particle, cur + 1, N);
		return;
	}


	Vector3d prevVelocity = velocity;
	int i = 0;
	Vector3d velocityHat = (particle->rotationTensor * particle->gammaFactor(
	) * velocity);
	Vector3d Eperp = E - velocity * (velocity.scalarMult(E) / (velocity.scalarMult(velocity)));
	Vector3d electricVelocityShift = (Eperp * (2 * eta * beta / gamma));
	//velocityHat += electricVelocityShift;

	//if (velocityHat.norm() > speed_of_light_normalized) {
	if (velocityHat.norm() > 1.0) {
		//printf("velocity Hat norm > c\n");
		//MPI_Finalize();
		//exit(0);
	}
	double a = tempParticle.gammaFactor();
	double etaDeltaT = eta * deltaT / N;
	double restEtaDeltaT = (1.0 - eta) * deltaT / N;
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

		if ((tempParticle.coordinates.x < xgrid[1 + additionalBinNumber]) && (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT || boundaryConditionTypeX ==
				FREE_MIRROR_BOTH)
			&& (cartCoord[0] == 0)) {
			particle->coordinates.x = 2 * xgrid[1 + additionalBinNumber] - tempParticle.coordinates.x + fabs(
				middleVelocity.x * restEtaDeltaT);
			particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * restEtaDeltaT;
			particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * restEtaDeltaT;

			newVelocity.x = fabs(newVelocity.x);
			particle->setMomentumByV(newVelocity);
			theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
			//sortParticleToEscaped(particle);
			moveParticle(particle, cur + 1, N);
			return;

		}

		if ((tempParticle.coordinates.x > xgrid[xnumberAdded - 1 - additionalBinNumber]) && (boundaryConditionTypeX == FREE_MIRROR_BOTH)
			&& (cartCoord[0] == cartDim[0] - 1)) {
			particle->coordinates.x = 2 * xgrid[xnumberAdded - 1 - additionalBinNumber] - tempParticle.coordinates.x - fabs(
				middleVelocity.x * restEtaDeltaT);
			particle->coordinates.y = tempParticle.coordinates.y + middleVelocity.y * restEtaDeltaT;
			particle->coordinates.z = tempParticle.coordinates.z + middleVelocity.z * restEtaDeltaT;

			newVelocity.x = -fabs(newVelocity.x);
			particle->setMomentumByV(newVelocity);
			theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
			//sortParticleToEscaped(particle);
			moveParticle(particle, cur + 1, N);
			return;

		}


		//correctParticlePosition(tempParticle);
		updateCorrelationMaps(&tempParticle);

		if (solverType == BUNEMAN) {
			E = correlationBunemanEfield(&tempParticle);
			B = correlationBunemanBfield(&tempParticle);
		} else {
			E = correlationEfield(&tempParticle) * (1 - localTheta - theta / N) + correlationNewEfield(&tempParticle) * (localTheta
				+ theta / N);
			//E = correlationNewEfield(particle);
			//B = correlationBfield(particle)*(1-theta) + correlationNewBfield(particle)*theta;
			B = correlationBfield(&tempParticle) * (1 - localTheta) + correlationNewBfield(&tempParticle) * localTheta;
		}
		//E = E0;
		//B = B0;

		//tempParticle.addMomentum((E + (middleVelocity.vectorMult(B) / speed_of_light_normalized)) * (particle->charge * deltaT / N));
		tempParticle.addMomentum((E + middleVelocity.vectorMult(B)) * (particle->charge * deltaT / N));
		newVelocity = tempParticle.getVelocity();
		newMomentum = tempParticle.getMomentum();
		//error = (prevVelocity - newVelocity).norm();
		error = (prevMomentum - newMomentum).norm();
	}

	particle->copyMomentum(tempParticle);

	particle->coordinates.x += middleVelocity.x * deltaT / N;

	particle->coordinates.y += middleVelocity.y * deltaT / N;
	particle->coordinates.z += middleVelocity.z * deltaT / N;

	//correctParticlePosition(particle);
	/*if(particle->coordinates.x > xgrid[xnumberAdded - additionalBinNumber]){
		printf("aaa\n");
	}*/

	if (particle->coordinates.x < xgrid[1 + additionalBinNumber]) {
		if (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
			particle->coordinates.x = 2 * xgrid[1 + additionalBinNumber] - particle->coordinates.x;
			particle->reflectMomentumX();
			theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
			//return;
		}
	}

	if (particle->coordinates.x < xgrid[xnumberAdded - 1 - additionalBinNumber]) {
		if (boundaryConditionTypeX == FREE_MIRROR_BOTH && cartCoord[0] == cartDim[0] - 1) {
			particle->coordinates.x = 2 * xgrid[1 + additionalBinNumber] - particle->coordinates.x;
			particle->reflectMomentumX();
			theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
			//return;
		}
	}

	moveParticle(particle, cur + 1, N);
}

void Simulation::moveParticleTristan(Particle* particle) {

	Vector3d E;
	Vector3d B;

	//only for Bunman solver
	updateBunemanCorrelationMaps(particle);
	correlationBunemanEBfields(particle, bunemanEx, bunemanEy, bunemanEz, bunemanBx, bunemanBy, bunemanBz, E, B);

	//correlationBunemanEBfieldsWithoutMaps(particle, bunemanEx, bunemanEy, bunemanEz, bunemanBx, bunemanBy, bunemanBz, E, B);

	Vector3d dp = E * (particle->charge * deltaT * 0.5);
	//Vector3d momentum = particle->getMomentum();
	//momentum += dp;
	particle->addMomentum(dp);
	double gamma = particle->gammaFactor();
	//double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	//double gamma = sqrt(p2/(particle->mass*particle->mass) + 1.0);
	//double beta = particle->charge * deltaT / (2.0 * particle->mass * speed_of_light_normalized);
	double betaShift = particle->beta / gamma;
	double f = 2.0 / sqrt(1.0 + betaShift * betaShift * B.norm2());
	Vector3d momentum = particle->getMomentum();
	Vector3d tempMomentum = momentum + momentum.vectorMult(B) * (betaShift);
	momentum += tempMomentum.vectorMult(B) * (f * betaShift) + dp;
	particle->setMomentum(momentum);
	Vector3d velocity = particle->getVelocity();
	particle->coordinates += velocity * deltaT;

	if ((boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT || boundaryConditionTypeX == FREE_MIRROR_BOTH) && (cartCoord[0] == 0) && particle->coordinates.x < xgrid[reflectingWallPoint]) {
		particle->coordinates.x = 2 * xgrid[reflectingWallPoint] - particle->coordinates.x;
		particle->reflectMomentumX();
		theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
		//return;
	} else {
		if (particle->coordinates.x < xgrid[1 + additionalBinNumber]) {
			particle->escaped = true;
			particle->crossBoundaryCount++;
			escapedParticlesLeft.push_back(particle);
			return;
		}
	}
	if (particle->coordinates.x > xgrid[xnumberAdded - 1 - additionalBinNumber]) {
		if ((boundaryConditionTypeX == FREE_MIRROR_BOTH) && cartCoord[0] == cartDim[0] - 1) {
			particle->coordinates.x = 2 * xgrid[xnumberAdded - 1 - additionalBinNumber] - particle->coordinates.x;
			particle->reflectMomentumX();
			theoreticalMomentum.x += particle->getMomentum().x * (2 * particle->weight * scaleFactor / plasma_period);
			//return;
		} else {
			particle->escaped = true;
			particle->crossBoundaryCount++;
			escapedParticlesRight.push_back(particle);
			return;
		}
	}
	/*if (particle->coordinates.x > xgrid[xnumberAdded - 1 - additionalBinNumber]) {
		escapedParticlesRight.push_back(particle);
		particle->crossBoundaryCount++;
		particle->escaped = true;
		return;
	}*/

	if (ynumberGeneral > 1) {
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
	}
	if (znumberGeneral > 1) {
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
}

void Simulation::evaluateParticlesRotationTensor() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	for (int i = 0; i < particles.size(); ++i) {
		Particle* particle = particles[i];
		double beta = 0.5 * particle->charge * deltaT / particle->mass;
		double gamma = particle->gammaFactor();
		Vector3d velocity = particle->getVelocity();

		Vector3d oldE;
		Vector3d oldB;

		if (solverType == BUNEMAN) {
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
	//MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("evaluating ParticlesRotationTensor time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::updateParticlesBeta() {
	for (int i = 0; i < particles.size(); ++i) {
		Particle* particle = particles[i];
		//particle->beta = particle->charge * deltaT / (2.0 * particle->mass * speed_of_light_normalized);
		particle->beta = particle->charge * deltaT / (2.0 * particle->mass);
	}
}

Matrix3d Simulation::evaluateAlphaRotationTensor(const double& beta, Vector3d& velocity, double& gamma, Vector3d& EField,
                                                 Vector3d& BField) {
	Matrix3d result = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);

	//double G = ((beta * (EField.scalarMult(velocity)) / (speed_of_light_normalized_sqr * speed_of_light_correction_sqr)) + gamma);
	//double G = ((beta * (EField.scalarMult(velocity)) / (speed_of_light_normalized_sqr)) + gamma);
	double G = ((beta * EField.scalarMult(velocity)) + gamma);
	double betaShift = beta / G;
	//double beta2c = betaShift * betaShift / (speed_of_light_normalized_sqr * speed_of_light_correction_sqr);
	//double beta2c = betaShift * betaShift / (speed_of_light_normalized_sqr);
	double beta2c = betaShift * betaShift;
	double denominator = G * (1 + beta2c * BField.scalarMult(BField));

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			result.matrix[i][j] = Kronecker.matrix[i][j] + (beta2c * BField[i] * BField[j]);
			//for (int k = 0; k < 3; ++k) {
			for (int l = 0; l < 3; ++l) {
				if (LeviCivita[j][i][l] != 0) {
					//result.matrix[i][j] -= (betaShift * LeviCivita[j][i][l] * BField[l] / (speed_of_light_normalized *speed_of_light_correction));
					//result.matrix[i][j] -= (betaShift * LeviCivita[j][i][l] * BField[l] / (speed_of_light_normalized));
					result.matrix[i][j] -= (betaShift * LeviCivita[j][i][l] * BField[l]);
					//result.matrix[i][j] += (betaShift * LeviCivita[j][k][l] * Kronecker.matrix[i][k] * BField[l] / speed_of_light_normalized);
				}
			}
			//}

			result.matrix[i][j] /= denominator;
			//if (debugMode) alertNaNOrInfinity(result.matrix[i][j], "rotation tensor = NaN");
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

void Simulation::injectNewParticles() {
	for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
		for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
			if (cartJ == cartCoord[1] && cartK == cartCoord[2] && cartCoord[0] == cartDim[0] - 1) {
				for (int typeCounter = 0; typeCounter < typesNumber; ++typeCounter) {
					if (types[typeCounter].particlesPerBin > 0) {
						if (types[typeCounter].particlesPerBin * types[typeCounter].injectionLength >= deltaX) {
							int newParticlesCount = types[typeCounter].injectionLength * types[typeCounter].particlesPerBin / deltaX;
							for (int i = newParticlesCount; i > 0; --i) {
								types[typeCounter].injectionLength -= deltaX / types[typeCounter].particlesPerBin;
								injectNewParticles(1, types[typeCounter], types[typeCounter].injectionLength);
							}
						}
					}
				}
			}
			int tempRank;
			int tempCartCoord[3];
			tempCartCoord[0] = cartDim[0] - 1;
			tempCartCoord[1] = cartJ;
			tempCartCoord[2] = cartK;
			MPI_Cart_rank(cartComm, tempCartCoord, &tempRank);
			MPI_Barrier(cartComm);
			int tempParticleNumber[1];
			tempParticleNumber[0] = particlesNumber;
			MPI_Bcast(tempParticleNumber, 1, MPI_INT, tempRank, cartComm);
			particlesNumber = tempParticleNumber[0];
		}
	}
}

void Simulation::injectNewParticles(int count, ParticleTypeContainer& typeContainer, const double& length) {
	if ((rank == 0) && (verbosity > 0)) printf("inject new particles\n");

	//Particle* tempParticle = particles[0];

	double x = xgrid[xnumberAdded - 1 - additionalBinNumber] - length;
	double tempDeltaY = deltaY * uniformDistribution();
	double tempDeltaZ = deltaZ * uniformDistribution();

	/*if (typeContainer.type == ELECTRON && preserveChargeLocal) {
	    return;
	}*/

	int minJ = 1 + additionalBinNumber;
	int minK = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - 1 - additionalBinNumber;
	int maxK = znumberAdded - 1 - additionalBinNumber;
	if (ynumberGeneral == 1) {
		minJ = 0;
		maxJ = 1;
	}
	if (znumberGeneral == 1) {
		minK = 0;
		maxK = 1;
	}
	for (int j = minJ; j < maxJ; ++j) {
		for (int k = minK; k < maxK; ++k) {
			double weight = (typeContainer.concentration / typeContainer.particlesPerBin) * volumeB();
			for (int l = 0; l < count; ++l) {
				ParticleTypes type = typeContainer.type;
				if (verbosity > 1) {
					printf("inject particle number = %d\n", particlesNumber);
				}
				Particle* particle = createParticle(particlesNumber, xnumberAdded - 1 - additionalBinNumber, j, k, weight, type,
				                                    typeContainer,
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

				particle->addVelocity(V0);
				Vector3d momentum = particle->getMomentum();
				particle->initialMomentum = momentum;
				//particle->prevMomentum = momentum;
				particles.push_back(particle);
				double en = particle->energy() * particle->weight * sqr(
					scaleFactor / plasma_period);
				theoreticalEnergy += particle->energy() * particle->weight * sqr(
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
		if (boundaryConditionTypeX == PERIODIC) {
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
		if (cartCoord[0] == 0 && boundaryConditionTypeX == PERIODIC) {
			for (int i = 0; i < escapedParticlesLeft.size(); ++i) {
				Particle* particle = escapedParticlesLeft[i];
				particle->coordinates.x += xsizeGeneral;
			}
		}
		if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX == PERIODIC) {
			for (int i = 0; i < escapedParticlesRight.size(); ++i) {
				Particle* particle = escapedParticlesRight[i];
				particle->coordinates.x -= xsizeGeneral;
			}
		}
		/*for(int i = 0; i < escapedParticlesLeft.size(); ++i) {
			Particle* particle = escapedParticlesLeft[i];
			reservedParticles.push_back(particle);
		}*/
		if (verbosity > 2) printf("send particles left rank = %d\n", rank);
		sendLeftReceiveRightParticles(escapedParticlesLeft, tempParticles, reservedParticles, types,
		                              typesNumber, boundaryConditionTypeX == PERIODIC, verbosity, cartComm, rank, leftRank, rightRank, deltaT);
		MPI_Barrier(cartComm);
		/*for(int i = 0; i < escapedParticlesRight.size(); ++i) {
			Particle* particle = escapedParticlesRight[i];
			reservedParticles.push_back(particle);
		}*/
		if (verbosity > 2) printf("send particles right rank = %d\n", rank);
		sendRightReceiveLeftParticles(escapedParticlesRight, tempParticles, reservedParticles, types,
		                              typesNumber, boundaryConditionTypeX == PERIODIC, verbosity, cartComm, rank, leftRank, rightRank, deltaT);
	}

	for (int pcount = 0; pcount < tempParticles.size(); ++pcount) {
		Particle* particle = tempParticles[pcount];
		if ((particle->coordinates.y < ygrid[1 + additionalBinNumber]) && (ynumberGeneral > 1)) {
			escapedParticlesFront.push_back(particle);
			particle->escaped = true;
			particle->crossBoundaryCount++;
		} else if ((particle->coordinates.y > ygrid[ynumberAdded - 1 - additionalBinNumber]) && (ynumberGeneral > 1)) {
			escapedParticlesBack.push_back(particle);
			particle->escaped = true;
			particle->crossBoundaryCount++;
		} else if ((particle->coordinates.z < zgrid[1 + additionalBinNumber]) && (znumberGeneral > 1)) {
			escapedParticlesBottom.push_back(particle);
			particle->escaped = true;
			particle->crossBoundaryCount++;
		} else if ((particle->coordinates.z > zgrid[znumberAdded - 1 - additionalBinNumber]) && (znumberGeneral > 1)) {
			escapedParticlesTop.push_back(particle);
			particle->escaped = true;
			particle->crossBoundaryCount++;
		} else {
			theoreticalEnergy += particle->energy() * particle->weight *
				sqr(scaleFactor / plasma_period);
			theoreticalMomentum += particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
			particle->escaped = false;
			particles.push_back(particle);
		}
	}
	tempParticles.clear();

	MPI_Barrier(cartComm);

	if (ynumberGeneral > 1) {
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
			/*for(int i = 0; i < escapedParticlesFront.size(); ++i) {
				Particle* particle = escapedParticlesFront[i];
				reservedParticles.push_back(particle);
			}*/
			if (verbosity > 2) printf("send particles front rank = %d\n", rank);
			sendFrontReceiveBackParticles(escapedParticlesFront, tempParticles, reservedParticles, types,
			                              typesNumber, boundaryConditionTypeY == PERIODIC, verbosity, cartComm, rank, frontRank, backRank, deltaT);
			MPI_Barrier(cartComm);
			/*for(int i = 0; i < escapedParticlesBack.size(); ++i) {
				Particle* particle = escapedParticlesBack[i];
				reservedParticles.push_back(particle);
			}*/
			if (verbosity > 2) printf("send particles back rank = %d\n", rank);
			sendBackReceiveFrontParticles(escapedParticlesBack, tempParticles, reservedParticles, types,
			                              typesNumber, boundaryConditionTypeY == PERIODIC, verbosity, cartComm, rank, frontRank, backRank, deltaT);
		}

		for (int pcount = 0; pcount < tempParticles.size(); ++pcount) {
			Particle* particle = tempParticles[pcount];
			if ((particle->coordinates.z < zgrid[1 + additionalBinNumber]) && (znumberGeneral > 1)) {
				escapedParticlesBottom.push_back(particle);
				particle->escaped = true;
				particle->crossBoundaryCount++;
			} else if ((particle->coordinates.z > zgrid[znumberAdded - 1 - additionalBinNumber]) && (znumberGeneral > 1)) {
				escapedParticlesTop.push_back(particle);
				particle->escaped = true;
				particle->crossBoundaryCount++;
			} else {
				theoreticalEnergy += particle->energy() * particle->weight *
					sqr(scaleFactor / plasma_period);
				theoreticalMomentum += particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
				particle->escaped = false;
				particles.push_back(particle);
			}
		}
		tempParticles.clear();
	}

	MPI_Barrier(cartComm);
	if (znumberGeneral > 1) {
		if (cartDim[2] == 1) {
			for (int i = 0; i < escapedParticlesBottom.size(); ++i) {
				Particle* particle = escapedParticlesBottom[i];
				particle->coordinates.z += zsizeGeneral;
				particle->escaped = false;
				tempParticles.push_back(particle);
			}
			for (int i = 0; i < escapedParticlesTop.size(); ++i) {
				Particle* particle = escapedParticlesTop[i];
				particle->coordinates.z -= zsizeGeneral;
				particle->escaped = false;
				tempParticles.push_back(particle);
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
			/*for(int i = 0; i < escapedParticlesBottom.size(); ++i) {
				Particle* particle = escapedParticlesBottom[i];
				reservedParticles.push_back(particle);
			}*/
			if (verbosity > 2) printf("send particles bottom rank = %d\n", rank);
			sendBottomReceiveTopParticles(escapedParticlesBottom, tempParticles, reservedParticles, types,
			                              typesNumber, boundaryConditionTypeZ == PERIODIC, verbosity, cartComm, rank, bottomRank, topRank, deltaT);
			MPI_Barrier(cartComm);
			/*for(int i = 0; i < escapedParticlesTop.size(); ++i) {
				Particle* particle = escapedParticlesTop[i];
				reservedParticles.push_back(particle);
			}*/
			if (verbosity > 2) printf("send particles top rank = %d\n", rank);
			sendTopReceiveBottomParticles(escapedParticlesTop, tempParticles, reservedParticles, types,
			                              typesNumber, boundaryConditionTypeZ == PERIODIC, verbosity, cartComm, rank, bottomRank, topRank, deltaT);
		}

		for (int pcount = 0; pcount < tempParticles.size(); ++pcount) {
			Particle* particle = tempParticles[pcount];
			theoreticalEnergy += particle->energy() * particle->weight *
				sqr(scaleFactor / plasma_period);
			theoreticalMomentum += particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
			particle->escaped = false;
			particles.push_back(particle);
		}
		tempParticles.clear();
	}

	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("exchange particles time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}
