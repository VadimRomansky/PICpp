#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>

#include "simulation.h"
#include "util.h"
#include "output.h"
#include "constants.h"
#include "matrix3d.h"
#include "random.h"


void Simulation::simulate() {
	if(newlyStarted){
		createArrays();
		initialize();
		createFiles();
		//initializeTwoStream();
		createParticles();
		//collectParticlesIntoBins();
		//initializeExternalFluxInstability();
		initializeAlfvenWave();
		//initializeSimpleElectroMagneticWave();
		//initializeLangmuirWave();
	}
	collectParticlesIntoBins();
	updateParameters();
	updateDeltaT();
	evaluateParticlesRotationTensor();
	updateElectroMagneticParameters();
	updateDensityParameters();

	evaluateExplicitDerivative();
	//cleanupDivergence();
	updateFields();
	updateEnergy();

	while (time * plasma_period < maxTime && currentIteration < maxIteration) {
		printf("iteration number %d time = %15.10g\n", currentIteration, time * plasma_period);
		printf(" dt/plasma_period = %15.10g\n", deltaT);

		if (currentIteration % writeParameter == 0) {
			output();
		}

		updateDeltaT();
		evaluateParticlesRotationTensor();
		updateElectroMagneticParameters();
		evaluateFields();
		evaluateMagneticField();
		moveParticles();
		//cleanupDivergence();
		updateDensityParameters();
		updateFields();
		updateEnergy();
		updateParameters();

		time += deltaT;
		currentIteration++;

		if(currentIteration % writeBackupParameter == 0){
			outputBackup();
		}
	}
}

void Simulation::output() {
	printf("outputing\n");

	if (particles.size() > 0) {
		distributionFile = fopen("./output/distribution_protons.dat", "a");
		outputDistribution(distributionFile, particles, PROTON);
		fclose(distributionFile);
		protonTraectoryFile = fopen("./output/traectory_proton.dat", "a");
		outputTraectory(protonTraectoryFile, getFirstProton(), time, plasma_period, gyroradius);
		fclose(protonTraectoryFile);
		electronTraectoryFile = fopen("./output/traectory_electron.dat", "a");
		outputTraectory(electronTraectoryFile, getFirstElectron(), time, plasma_period, gyroradius);
		fclose(electronTraectoryFile);
	}

	EfieldFile = fopen("./output/Efield.dat", "a");
	BfieldFile = fopen("./output/Bfield.dat", "a");
	outputFields(EfieldFile, BfieldFile, Efield, Bfield, xnumber, plasma_period, gyroradius, fieldScale);
	fclose(EfieldFile);
	fclose(BfieldFile);

	Xfile = fopen("./output/Xfile.dat", "w");
	outputGrid(Xfile, xgrid, xnumber, gyroradius);
	fclose(Xfile);

	densityFile = fopen("./output/concentrations.dat", "a");
	outputConcentrations(densityFile, electronConcentration, protonConcentration, chargeDensity, electricDensity, xnumber, plasma_period, gyroradius, fieldScale);
	fclose(densityFile);

	velocityFile = fopen("./output/velocity.dat", "a");
	velocityElectronFile = fopen("./output/velocity_electron.dat", "a");
	outputVelocity(velocityFile, velocityElectronFile, velocityBulk, velocityBulkElectron, xnumber, plasma_period, gyroradius);
	fclose(velocityFile);
	fclose(velocityElectronFile);

	fluxFile = fopen("./output/fluxFile.dat", "a");
	outputFlux(fluxFile, electricFlux, externalElectricFlux, xnumber, plasma_period, gyroradius, fieldScale);
	fclose(fluxFile);

	divergenceErrorFile = fopen("./output/divergence_error.dat", "a");
	outputDivergenceError(divergenceErrorFile, this, plasma_period, gyroradius, fieldScale);
	fclose(divergenceErrorFile);

	double rotBscale = fieldScale/(plasma_period * plasma_period * sqrt(gyroradius));

	rotBFile = fopen("./output/rotBFile.dat", "a");
	outputVectorArray(rotBFile, rotB, xnumber + 1, rotBscale);
	fclose(rotBFile);

	EderivativeFile = fopen("./output/EderivativeFile.dat", "a");
	outputVectorArray(EderivativeFile, Ederivative, xnumber + 1, rotBscale);
	fclose(EderivativeFile);

	generalFile = fopen("./output/general.dat", "a");
	outputGeneral(generalFile, this);
	fclose(generalFile);
}

void Simulation::outputBackup(){
	printf("writing backup\n");

	backupGeneralFile = fopen("./backup/general.dat", "w");
	backupEfieldFile = fopen("./backup/Efield.dat", "w");
	backupBfieldFile = fopen("./backup/Bfield.dat", "w");
	backupParticlesFile = fopen("./backup/particles.dat", "w");

	outputSimulationBackup(backupGeneralFile, backupEfieldFile, backupBfieldFile, backupParticlesFile, this);

	fclose(backupGeneralFile);
	fclose(backupEfieldFile);
	fclose(backupBfieldFile);
	fclose(backupParticlesFile);
}

void Simulation::updateDeltaT() {
	printf("updating time step\n");
	double timeEpsilon = 0.002;
	double delta = deltaX;
	deltaT = timeEpsilon * delta / speed_of_light_normalized;
	double B = B0.norm();
	double E = E0.norm();
	for (int i = 0; i < xnumber; ++i) {
		if (Bfield[i].norm() > B) {
			B = Bfield[i].norm();
		}
	}
	for (int i = 0; i < xnumber; ++i) {
		if (Efield[i].norm() > E) {
			E = Efield[i].norm();
		}
	}

	B *= fieldScale;
	E *= fieldScale;

	double minFlux = electricFlux[0].norm();
	for(int i = 0; i < xnumber; ++i){
		if(electricFlux[i].norm() < minFlux){
			minFlux = electricFlux[i].norm();
		}
	}

	/*if(E > 0 && minFlux > 0){
		double maxResistance = 0;
		double minResistance = Efield[0].norm()*fieldScale/electricFlux[0].norm();
		for(int i = 0; i < xnumber; ++i){
			if (Efield[i].norm()*fieldScale/electricFlux[i].norm() < minResistance){
				minResistance = Efield[i].norm()*fieldScale/electricFlux[i].norm();
			}
			if (Efield[i].norm()*fieldScale/electricFlux[i].norm() > maxResistance){
				maxResistance = Efield[i].norm()*fieldScale/electricFlux[i].norm();
			}
		}

		deltaT = min2(deltaT, 0.05*(minResistance/(4*pi)));
		//note that omega plasma = 1 in our system
		deltaT = min2(deltaT, 0.05/maxResistance);
	}*/

	double concentration = density / (massProton + massElectron);

	omegaPlasmaProton = sqrt(4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massProton);
	omegaPlasmaElectron = sqrt(4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massElectron);

	deltaT = min2(deltaT, timeEpsilon/omegaPlasmaElectron);

	omegaGyroProton = electron_charge_normalized * B / (massProton * speed_of_light);
	omegaGyroElectron = electron_charge_normalized * B / (massElectron * speed_of_light);

	double thermalMomentum = sqrt(massElectron*kBoltzman_normalized*temperature);

	if (B > 0) {
		deltaT = min2(deltaT, timeEpsilon * massElectron * speed_of_light_normalized / (electron_charge_normalized * B));
	}
	/*if (E > 0) {
		deltaT = min2(deltaT, 0.1 * thermalMomentum / (electron_charge_normalized * E));
	}*/


	double Vthermal = sqrt(2*kBoltzman_normalized*temperature/massElectron);
	double minDeltaT = deltaX/Vthermal;
	if(minDeltaT > deltaT){
		printf("deltaT < dx/Vthermal\n");
	}
}

double Simulation::volume(int i) {
	return deltaX;
}

void Simulation::checkParticleInBox(Particle& particle) {
	if (particle.x < 0) {
		printf("particle.x < 0\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "particle.x = %15.10g < 0\n", particle.x);
		fclose(errorLogFile);
		exit(0);
	}
	if (particle.x > xgrid[xnumber]) {
		printf("particle.x > xsize\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "particle.x = %15.10g > %15.10g\n", particle.x, xgrid[xnumber]);
		fclose(errorLogFile);
		exit(0);
	}
}

void Simulation::checkParticlesInBin() {
	for(int pcount = 0; pcount < particlesInEbin[0].size(); ++pcount) {
		Particle* particle = particlesInEbin[0][pcount];
		for(int pcountRight = 0; pcountRight < particlesInEbin[xnumber].size(); ++pcountRight) {
			if(particle == particlesInEbin[xnumber][pcountRight]) {
				printf("particle is in 0 and xnumber bin\n");
				errorLogFile = fopen("./output/errorLog.dat", "w");
				fprintf(errorLogFile, "particle is in 0 and xnumber bin, particle.x = %15.10g xmax = %15.10g\n", particle->x, xgrid[xnumber]);
				fclose(errorLogFile);
				exit(0);
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		if(particlesInEbin[i].size() > 1){
			for(int pcount = 0; pcount < particlesInEbin[i].size() - 1; ++pcount) {
				Particle* particle = particlesInEbin[i][pcount];
				for(int pcount2 = pcount + 1; pcount2 < particlesInEbin[i].size(); ++pcount2) {
					if(particle == particlesInEbin[i][pcount2]) {
						printf("particle is twice in Ebin number %d\n", i);
						errorLogFile = fopen("./output/errorLog.dat", "w");
						fprintf(errorLogFile, "particle is twice in Ebin number %d\n", i);
						fclose(errorLogFile);
						exit(0);
					}
				}
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		if(particlesInBbin[i].size() > 1){
			for(int pcount = 0; pcount < particlesInBbin[i].size() - 1; ++pcount) {
				Particle* particle = particlesInBbin[i][pcount];
				for(int pcount2 = pcount + 1; pcount2 < particlesInBbin[i].size(); ++pcount2) {
					if(particle == particlesInBbin[i][pcount2]) {
						printf("particle is twice in Bbin number %d\n", i);
						errorLogFile = fopen("./output/errorLog.dat", "w");
						fprintf(errorLogFile, "particle is twice in Bbi number %d\n", i);
						fclose(errorLogFile);
						exit(0);
					}
				}
			}
		}
	}
}

void Simulation::updateElectroMagneticParameters() {
	printf("updating flux, density snd dielectric tensor\n");
	//check particle only in one boundary box
	if(debugMode) {
		//checkParticlesInBin();
	}
	//collectParticlesIntoBins();
	for (int i = 0; i < xnumber + 1; ++i) {
		electricFlux[i] = Vector3d(0, 0, 0);
		dielectricTensor[i] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
		for (int pcount = 0; pcount < particlesInEbin[i].size(); ++pcount) {
			Particle* particle = particlesInEbin[i][pcount];
			double correlation = correlationWithEbin(*particle, i) / volume(i);

			Vector3d velocity = particle->velocity(speed_of_light_normalized);
			double gamma = particle->gammaFactor(speed_of_light_normalized);
			Vector3d rotatedVelocity = particle->rotationTensor * (velocity * gamma);

			if(solverType == IMPLICIT){
				electricFlux[i] += rotatedVelocity * (particle->charge * particle->weight * correlation);
				dielectricTensor[i] = dielectricTensor[i] - particle->rotationTensor * (particle->weight*theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge * correlation / particle->mass);
			}
			if(solverType == EXPLICIT){
				electricFlux[i] += velocity*particle->charge*particle->weight*correlation;
			}
			alertNaNOrInfinity(electricFlux[i].x, "right part x = NaN");
			alertNaNOrInfinity(electricFlux[i].y, "right part y = NaN");
			alertNaNOrInfinity(electricFlux[i].z, "right part z = NaN");
		}
	}

	//for i = 0


	/*for (int pcount = 0; pcount < particlesInEbin[0].size(); ++pcount) {
		Particle* particle = particlesInEbin[0][pcount];
		double correlation = correlationWithEbin(*particle, -1) / volume(0);

		if(correlation > 0){
			Vector3d velocity = particle->velocity(speed_of_light_normalized);

			double gamma = particle->gammaFactor(speed_of_light_normalized);
			Vector3d rotatedVelocity = particle->rotationTensor * (velocity * gamma);

			electricFlux[xnumber - 1] += rotatedVelocity * (particle->charge * particle->weight * correlation);
			dielectricTensor[xnumber - 1] = dielectricTensor[xnumber - 1] - particle->rotationTensor * (theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge * correlation / particle->mass);
		}
		if(correlation < 0){
			printf("aaa\n");
		}
	}*/
	/*for (int pcount = 0; pcount < particlesInEbin[xnumber].size(); ++pcount) {
		Particle* particle = particlesInEbin[xnumber][pcount];
		double correlation = correlationWithEbin(*particle, xnumber + 1) / volume(0);
		if(correlation > 0){

			Vector3d velocity = particle->velocity(speed_of_light_normalized);

			double gamma = particle->gammaFactor(speed_of_light_normalized);
			Vector3d rotatedVelocity = particle->rotationTensor * velocity * gamma;
	
			electricFlux[1] += rotatedVelocity * particle->charge * particle->weight * correlation;
			dielectricTensor[1] = dielectricTensor[1] - particle->rotationTensor * (theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge * correlation / particle->mass);
		}
		if(correlation < 0){
			printf("aaa\n");
		}
	}*/


	//for periodic conditions we must summ sides parameters
	electricFlux[0] = electricFlux[0] + electricFlux[xnumber];
	electricFlux[xnumber] = electricFlux[0];

	dielectricTensor[0] = dielectricTensor[0] + dielectricTensor[xnumber];
	dielectricTensor[xnumber] = dielectricTensor[0];


	for (int i = 0; i < xnumber; ++i) {
		electricDensity[i] = 0;
		pressureTensor[i] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
		if(solverType == IMPLICIT){
			for (int pcount = 0; pcount < particlesInBbin[i].size(); ++pcount) {
				Particle* particle = particlesInBbin[i][pcount];
				double correlation = correlationWithBbin(*particle, i) / volume(i);

				double gamma = particle->gammaFactor(speed_of_light_normalized);
				Vector3d velocity = particle->velocity(speed_of_light_normalized);
				Vector3d rotatedVelocity = particle->rotationTensor * velocity * gamma;

				electricDensity[i] += particle->weight * particle->charge * correlation;

				pressureTensor[i] += rotatedVelocity.tensorMult(rotatedVelocity) * particle->weight * particle->charge * correlation;
			}
		}
	}

	if(solverType == IMPLICIT){
		for (int i = 0; i <= xnumber; ++i) {
			Vector3d divPressureTensor = evaluateDivPressureTensor(i);
			electricFlux[i] = electricFlux[i] - divPressureTensor * deltaT / 2;
		}
	}

	//smoothFlux();
	//smoothFlux();
	//smoothFlux();

	if(solverType == IMPLICIT){
		for (int i = 0; i < xnumber; ++i) {
			double divJ = evaluateDivFlux(i);

			electricDensity[i] -= deltaT * theta * divJ;
		}
	}

	updateExternalFlux();

	for(int i = 0; i < xnumber +1 ; ++i){
		electricFlux[i] = electricFlux[i] + externalElectricFlux[i];
	}
}

void Simulation::updateDensityParameters() {
	double full_density = 0;
	double full_p_concentration = 0;
	double full_e_concentration = 0;
	collectParticlesIntoBins();
	//FILE* debugFile = fopen("./output/particleCorrelations.dat","w");
	for (int i = 0; i < xnumber; ++i) {
		electronConcentration[i] = 0;
		protonConcentration[i] = 0;
		chargeDensity[i] = 0;
		velocityBulk[i] = Vector3d(0, 0, 0);
		velocityBulkElectron[i] = Vector3d(0, 0, 0);
		//fprintf(debugFile, "%d %d %d\n", i, j, k);
		for (int pcount = 0; pcount < particlesInBbin[i].size(); ++pcount) {
			Particle* particle = particlesInBbin[i][pcount];

			double correlation = correlationWithBbin(*particle, i) / volume(i);

			chargeDensity[i] += correlation * particle->charge * particle->weight;
			if (particle->type == ELECTRON) {
				electronConcentration[i] += correlation * particle->weight;
				velocityBulkElectron[i] += particle->momentum * particle->weight * correlation;
			} else if (particle->type == PROTON) {
				protonConcentration[i] += correlation * particle->weight;
			}
			velocityBulk[i] += particle->momentum * particle->weight * correlation;

			if (correlation == 0) {
				//printf("aaa\n");
			}

			//fprintf(debugFile, "%d %15.10g\n", particle->number, correlation*volume(i, j, k));
		}

		//fprintf(debugFile, "charge %15.10g proton %15.10g electron %15.10g\n", chargeDensity[i][j][k], protonConcentration[i][j][k], electronConcentration[i][j][k]);

 		velocityBulk[i] = velocityBulk[i] / (electronConcentration[i] * massElectron + protonConcentration[i] * massProton);
		velocityBulkElectron[i] = velocityBulkElectron[i] / (electronConcentration[i] * massElectron);
		//if(i < xnumber - 1  && i > 0){
		full_density += chargeDensity[i] * volume(i);
		full_p_concentration += protonConcentration[i] * volume(i);
		full_e_concentration += electronConcentration[i] * volume(i);
		//electronConcentration[i][j][k] /= cube(gyroradius);
		//protonConcentration[i][j][k] /= cube(gyroradius);
		//chargeDensity[i][j][k] /= (sqrt(cube(gyroradius))*plasma_period);
		//}
	}
	full_density /= xsize;
	full_p_concentration /= (xsize * gyroradius);
	full_e_concentration /= (xsize * gyroradius);
	//fprintf(debugFile, "charge %15.10g proton %15.10g electron %15.10g\n", full_density, full_p_concentration, full_e_concentration);
	//fclose(debugFile);
}

void Simulation::updateEnergy() {
	particleEnergy = 0;
	electricFieldEnergy = 0;
	magneticFieldEnergy = 0;
	energy = 0;

	momentum = Vector3d(0, 0, 0);
	for (int i = 0; i < xnumber + 1; ++i) {
		double factor = 1;
		if (i == 0 || i == xnumber) {
			factor = factor / 2;
		}
		Vector3d E = Efield[i]*fieldScale;
		electricFieldEnergy += E.scalarMult(E) * volume(i) * factor / (8 * pi);
	}

	for (int i = 0; i < xnumber; ++i) {
		Vector3d B = (Bfield[i] - B0)*fieldScale;
		magneticFieldEnergy += B.scalarMult(B) * volume(i) / (8 * pi);
	}

	for (int i = 0; i < xnumber; ++i) {
		Vector3d E = ((Efield[i] + Efield[i + 1]) / 2)*fieldScale;
		Vector3d B = Bfield[i]*fieldScale;

		momentum += (E.vectorMult(B) / (4 * pi * speed_of_light_normalized)) * volume(i);
	}

	for(int pcount = 0; pcount < particles.size(); ++pcount){
		particleEnergy += particles[pcount]->energy(speed_of_light_normalized)*particles[pcount]->weight;
		momentum += particles[pcount]->momentum*particles[pcount]->weight;
	}


	particleEnergy *= sqr(gyroradius / plasma_period);
	electricFieldEnergy *= sqr(gyroradius / plasma_period);
	magneticFieldEnergy *= sqr(gyroradius / plasma_period);
	momentum = momentum * gyroradius / plasma_period;


	energy = particleEnergy + electricFieldEnergy + magneticFieldEnergy;
}

Vector3d Simulation::getBfield(int i) {
	if (i == -1) {
		i = xnumber - 1;
	} else if (i == xnumber) {
		i = 0;
	} else if (i < -1 || i > xnumber) {
		printf("i < -1 || i > xnumber in getBfied i = %d\n", i);
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "i = %d, xnumber = %d\n", i, xnumber);
		fclose(errorLogFile);
		exit(0);
	}

	return Bfield[i];
}

Vector3d Simulation::getTempEfield(int i) {
	if (i < 0) {
		return tempEfield[xnumber-1];
	} else if (i > xnumber) {
		return tempEfield[1];
	}

	return tempEfield[i];
}

Vector3d Simulation::getNewEfield(int i) {
	if (i < 0) {
		return newEfield[xnumber-1];
	} else if (i > xnumber) {
		return newEfield[1];
	}

	return newEfield[i];
}

Vector3d Simulation::getEfield(int i) {
	if (i < 0) {
		return Efield[xnumber-1];
	} else if (i > xnumber) {
		return Efield[1];
	}

	return Efield[i];
}

Matrix3d Simulation::getPressureTensor(int i) {
	if (i < 0) {
		i = xnumber - 1;
	} else if (i >= xnumber) {
		i = 0;
	}

	return pressureTensor[i];
}

double Simulation::getDensity(int i) {
	if (i < 0) {
		i = xnumber - 1;
	}
	if (i >= xnumber) {
		i = 0;
	}


	return electricDensity[i];
}

void Simulation::smoothFlux(){
	Vector3d* newFlux = new Vector3d[xnumber + 1];
	for(int i = 0; i < xnumber + 1; ++i){
		newFlux[i] = Vector3d(0, 0, 0);
		int prevI = i - 1;
		if(prevI < 0){
			prevI = prevI + xnumber;
		}
		int nextI = i + 1;
		if(nextI >= xnumber + 1){
			nextI = nextI - xnumber;
		}

		newFlux[i] = (electricFlux[prevI] + electricFlux[i]*2 + electricFlux[nextI])/4.0;
	}

	for(int i = 0; i < xnumber + 1; ++i){
		electricFlux[i] = newFlux[i];
	}

	delete[] newFlux;
}

void Simulation::smoothEderivative(){
	Vector3d* newDerivative = new Vector3d[xnumber + 1];
	for(int i = 0; i < xnumber + 1; ++i){
		newDerivative[i] = Vector3d(0, 0, 0);
		int prevI = i - 1;
		if(prevI < 0){
			prevI = prevI + xnumber;
		}
		int nextI = i + 1;
		if(nextI >= xnumber + 1){
			nextI = nextI - xnumber;
		}

		newDerivative[i] = (Ederivative[prevI] + Ederivative[i]*2 + Ederivative[nextI])/4.0;
	}

	for(int i = 0; i < xnumber + 1; ++i){
		Ederivative[i] = newDerivative[i];
	}

	delete[] newDerivative;
}

void Simulation::updateParameters(){
	maxBfield = Bfield[0];
	maxEfield = Efield[0];
	for(int i = 0; i < xnumber ; ++i){
		if(Bfield[i].norm() > maxBfield.norm()){
			maxBfield = Bfield[i];
		}
	}

	for(int i = 0; i < xnumber + 1; ++i){
		if(Efield[i].norm() > maxEfield.norm()){
			maxEfield = Efield[i];
		}
	}
}

void Simulation::updateExternalFlux(){
	double alfvenV = B0.norm()*fieldScale/sqrt(4*pi*density);
	double concentration = density/(massProton + massElectron);
	double phaseV = 2*alfvenV;
	double kw = 2*pi/xsize;
	double omega = kw*phaseV;

	for(int i = 0; i < xnumber +1 ; ++i){
		externalElectricFlux[i] = Vector3d(0, 0, 1.0)*extJ*cos(kw*xgrid[i] - omega*time);
	}
}

void Simulation::resetNewTempFields(){
	for(int i = 0; i < xnumber; ++i){
		newBfield[i] = Bfield[i];
	}

	for(int i = 0; i < xnumber + 1; ++i){
		tempEfield[i] = Efield[i];
		newEfield[i] = Efield[i];
		implicitEfield[i] = Efield[i];
	}
}