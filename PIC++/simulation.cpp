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
		createFiles();
		initialize();
		//initializeTwoStream();
		//createParticles();
		//initializeExternalFluxInstability();
		//initializeAlfvenWave(1, 0.01);
		//initializeRotatedAlfvenWave(1, 0.01);
		//initializeFluxFromRight();
		initializeSimpleElectroMagneticWave();
		//initializeLangmuirWave();
	}
	collectParticlesIntoBins();
	updateParameters();

	//double because dielectric tensor needs deltaT;
	updateDeltaT();
	evaluateParticlesRotationTensor();
	updateElectroMagneticParameters();
	updateDensityParameters();
	updateDeltaT();
	evaluateParticlesRotationTensor();
	updateElectroMagneticParameters();
	updateDensityParameters();

	evaluateExplicitDerivative();
	//cleanupDivergence();
	updateFields();
	updateEnergy();

	double length = deltaX/particlesPerBin - 0.0001*deltaX;

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
		updateParameters();
		moveParticles();

		length += fabs(V0.x*deltaT);
		if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
			if(length > deltaX/particlesPerBin){
				length -= deltaX/particlesPerBin;
				injectNewParticles(1);
			}
		}
		//cleanupDivergence();
		updateDensityParameters();
		updateFields();
		updateEnergy();

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
	outputFields(EfieldFile, BfieldFile, Efield, Bfield, xnumber, ynumber, znumber, plasma_period, gyroradius, fieldScale);
	fclose(EfieldFile);
	fclose(BfieldFile);

	Xfile = fopen("./output/Xfile.dat", "w");
	outputGrid(Xfile, xgrid, xnumber, gyroradius);
	fclose(Xfile);

	Yfile = fopen("./output/Yfile.dat", "w");
	outputGrid(Yfile, ygrid, ynumber, gyroradius);
	fclose(Yfile);

	Zfile = fopen("./output/Zfile.dat", "w");
	outputGrid(Zfile, zgrid, znumber, gyroradius);
	fclose(Zfile);

	densityFile = fopen("./output/concentrations.dat", "a");
	outputConcentrations(densityFile, electronConcentration, protonConcentration, chargeDensity, electricDensity, xnumber, ynumber, znumber, plasma_period, gyroradius, fieldScale);
	fclose(densityFile);

	velocityFile = fopen("./output/velocity.dat", "a");
	velocityElectronFile = fopen("./output/velocity_electron.dat", "a");
	outputVelocity(velocityFile, velocityElectronFile, velocityBulk, velocityBulkElectron, xnumber, ynumber, znumber, plasma_period, gyroradius);
	fclose(velocityFile);
	fclose(velocityElectronFile);

	fluxFile = fopen("./output/fluxFile.dat", "a");
	outputFlux(fluxFile, electricFlux, externalElectricFlux, xnumber, ynumber, znumber, plasma_period, gyroradius, fieldScale);
	fclose(fluxFile);

	divergenceErrorFile = fopen("./output/divergence_error.dat", "a");
	outputDivergenceError(divergenceErrorFile, this);
	fclose(divergenceErrorFile);

	double rotBscale = fieldScale/(plasma_period * plasma_period * sqrt(gyroradius));

	rotBFile = fopen("./output/rotBFile.dat", "a");
	outputVectorArray(rotBFile, rotB, xnumber + 1, ynumber + 1, znumber + 1, rotBscale);
	fclose(rotBFile);

	EderivativeFile = fopen("./output/EderivativeFile.dat", "a");
	outputVectorArray(EderivativeFile, Ederivative, xnumber + 1, ynumber + 1, znumber + 1, rotBscale);
	fclose(EderivativeFile);

	dielectricTensorFile = fopen("./output/dielectricTensorFile.dat", "a");
	outputMatrixArray(dielectricTensorFile, dielectricTensor, xnumber + 1, ynumber + 1, znumber + 1);
	fclose(dielectricTensorFile);

	particleProtonsFile = fopen("./output/protons.dat", "w");
	particleElectronsFile = fopen("./output/electrons.dat", "w");
	outputParticles(particleProtonsFile, particleElectronsFile, this);
	fclose(particleProtonsFile);
	fclose(particleElectronsFile);

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
	double delta = deltaX;
	deltaT = timeEpsilon * delta / speed_of_light_normalized;
	if(particles.size() > 0){
		double B = B0.norm();
		double E = E0.norm();
		for (int i = 0; i < xnumber; ++i) {
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					if (Bfield[i][j][k].norm() > B) {
						B = Bfield[i][j][k].norm();
					}
				}
			}
		}
		for (int i = 0; i < xnumber; ++i) {
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					if (Efield[i][j][k].norm() > E) {
						E = Efield[i][j][k].norm();
					}
				}
			}
		}

		B *= fieldScale;
		E *= fieldScale;

		double minFlux = electricFlux[0][0][0].norm();
		for(int i = 0; i < xnumber; ++i){
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k< znumber; ++k){
					if(electricFlux[i][j][k].norm() < minFlux){
						minFlux = electricFlux[i][j][k].norm();
					}
				}
			}
		}

		/*if(E > 0 && minFlux > 0){
			double maxResistance = 0;
			double minResistance = Efield[0].norm()*fieldScale/electricFlux[0].norm();
			for(int i = 0; i < xnumber; ++i){
				for(int j = 0; j < ynumber; ++j){
					for(int k = 0; k < znumber; ++K){
						if (Efield[i][j][k].norm()*fieldScale/electricFlux[i][j][k].norm() < minResistance){
							minResistance = Efield[i][j][k].norm()*fieldScale/electricFlux[i][j][k].norm();
						}
						if (Efield[i][j][k].norm()*fieldScale/electricFlux[i][j][k].norm() > maxResistance){
							maxResistance = Efield[i][j][k].norm()*fieldScale/electricFlux[i][j][k].norm();
						}
					}
				}
			}

			deltaT = min2(deltaT, timeEpsilon*(minResistance/(4*pi)));
			//note that omega plasma = 1 in our system
			deltaT = min2(deltaT,timeEpsilon/maxResistance);
		}*/

		double concentration = density / (massProton + massElectron);

		omegaPlasmaProton = sqrt(4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massProton);
		omegaPlasmaElectron = sqrt(4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massElectron);

		deltaT = min2(deltaT, timeEpsilon/omegaPlasmaElectron);

		omegaGyroProton = electron_charge_normalized * B / (massProton * speed_of_light);
		omegaGyroElectron = electron_charge_normalized * B / (massElectron * speed_of_light);

		double thermalMomentum = sqrt(massElectron*kBoltzman_normalized*temperature) + massElectron*V0.norm();

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
}

double Simulation::volumeE(int i, int j, int k) {
	double dx = deltaX;
	if((boundaryConditionType == PERIODIC) || ((i > 0) && (i < xnumber))){
		dx = deltaX;
	} else {
		dx = deltaX/2;
	}
	return dx*deltaY*deltaZ;
}

double Simulation::volumeB(int i, int j, int k) {
	return deltaX*deltaY*deltaZ;
}

void Simulation::checkParticleInBox(Particle& particle) {
	if (particle.coordinates.x < 0) {
		printf("particle.x < 0\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "particle.x = %15.10g < 0\n", particle.coordinates.x);
		fclose(errorLogFile);
		exit(0);
	}
	if (particle.coordinates.x > xgrid[xnumber]) {
		printf("particle.x > xsize\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "particle.x = %15.10g > %15.10g\n", particle.coordinates.x, xgrid[xnumber]);
		fclose(errorLogFile);
		exit(0);
	}

	if (particle.coordinates.y < 0) {
		printf("particle.y < 0\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "particle.y = %15.10g < 0\n", particle.coordinates.y);
		fclose(errorLogFile);
		exit(0);
	}
	if (particle.coordinates.y > ygrid[ynumber]) {
		printf("particle.y > ysize\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "particle.y = %15.10g > %15.10g\n", particle.coordinates.y, ygrid[ynumber]);
		fclose(errorLogFile);
		exit(0);
	}

	if (particle.coordinates.z < 0) {
		printf("particle.z < 0\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "particle.z = %15.10g < 0\n", particle.coordinates.z);
		fclose(errorLogFile);
		exit(0);
	}
	if (particle.coordinates.z > zgrid[znumber]) {
		printf("particle.z > zsize\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "particle.z = %15.10g > %15.10g\n", particle.coordinates.z, zgrid[znumber]);
		fclose(errorLogFile);
		exit(0);
	}
}

void Simulation::checkParticlesInBin() {
	for(int j = 0; j <ynumber; ++j){
		for(int k = 0; k < znumber; ++k){
			for(int pcount = 0; pcount < particlesInEbin[0][j][k].size(); ++pcount) {
				Particle* particle = particlesInEbin[0][j][k][pcount];
				for(int pcountRight = 0; pcountRight < particlesInEbin[xnumber][j][k].size(); ++pcountRight) {
					if(particle == particlesInEbin[xnumber][j][k][pcountRight]) {
						printf("particle is in 0 and xnumber bin\n");
						errorLogFile = fopen("./output/errorLog.dat", "w");
						fprintf(errorLogFile, "particle is in 0 and xnumber bin, particle.x = %15.10g xmax = %15.10g\n", particle->coordinates.x, xgrid[xnumber]);
						fclose(errorLogFile);
						exit(0);
					}
				}
			}
		}
	}

	for(int i = 0; i <xnumber; ++i){
		for(int k = 0; k < znumber; ++k){
			for(int pcount = 0; pcount < particlesInEbin[i][0][k].size(); ++pcount) {
				Particle* particle = particlesInEbin[i][0][k][pcount];
				for(int pcountRight = 0; pcountRight < particlesInEbin[i][ynumber][k].size(); ++pcountRight) {
					if(particle == particlesInEbin[i][ynumber][k][pcountRight]) {
						printf("particle is in 0 and ynumber bin\n");
						errorLogFile = fopen("./output/errorLog.dat", "w");
						fprintf(errorLogFile, "particle is in 0 and ynumber bin, particle.y = %15.10g ymax = %15.10g\n", particle->coordinates.y, ygrid[ynumber]);
						fclose(errorLogFile);
						exit(0);
					}
				}
			}
		}
	}

	for(int j = 0; j <ynumber; ++j){
		for(int i = 0; i < xnumber; ++i){
			for(int pcount = 0; pcount < particlesInEbin[i][j][0].size(); ++pcount) {
				Particle* particle = particlesInEbin[i][j][0][pcount];
				for(int pcountRight = 0; pcountRight < particlesInEbin[i][j][znumber].size(); ++pcountRight) {
					if(particle == particlesInEbin[i][j][znumber][pcountRight]) {
						printf("particle is in 0 and znumber bin\n");
						errorLogFile = fopen("./output/errorLog.dat", "w");
						fprintf(errorLogFile, "particle is in 0 and znumber bin, particle.z = %15.10g zmax = %15.10g\n", particle->coordinates.z, zgrid[znumber]);
						fclose(errorLogFile);
						exit(0);
					}
				}
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				if(particlesInEbin[i][j][k].size() > 1){
					for(int pcount = 0; pcount < particlesInEbin[i][j][k].size() - 1; ++pcount) {
						Particle* particle = particlesInEbin[i][j][k][pcount];
						for(int pcount2 = pcount + 1; pcount2 < particlesInEbin[i][j][k].size(); ++pcount2) {
							if(particle == particlesInEbin[i][j][k][pcount2]) {
								printf("particle is twice in Ebin number %d %d %d\n", i, j, k);
								errorLogFile = fopen("./output/errorLog.dat", "w");
								fprintf(errorLogFile, "particle is twice in Ebin number %d %d %d\n", i, j, k);
								fclose(errorLogFile);
								exit(0);
							}
						}
					}
				}
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				if(particlesInBbin[i][j][k].size() > 1){
					for(int pcount = 0; pcount < particlesInBbin[i][j][k].size() - 1; ++pcount) {
						Particle* particle = particlesInBbin[i][j][k][pcount];
						for(int pcount2 = pcount + 1; pcount2 < particlesInBbin[i][j][k].size(); ++pcount2) {
							if(particle == particlesInBbin[i][j][k][pcount2]) {
								printf("particle is twice in Bbin number %d %d %d\n", i, j, k);
								errorLogFile = fopen("./output/errorLog.dat", "w");
								fprintf(errorLogFile, "particle is twice in Bbi number %d %d 5d\n", i, j, k);
								fclose(errorLogFile);
								exit(0);
							}
						}
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
		for(int j = 0; j  < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				electricFlux[i][j][k] = Vector3d(0, 0, 0);
				dielectricTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				divPressureTensor[i][j][k] = Vector3d(0, 0, 0);
				for (int pcount = 0; pcount < particlesInEbin[i][j][k].size(); ++pcount) {
					Particle* particle = particlesInEbin[i][j][k][pcount];
					double correlation = correlationWithEbin(*particle, i, j, k) / volumeE(i, j, k);

					Vector3d velocity = particle->velocity(speed_of_light_normalized);
					double gamma = particle->gammaFactor(speed_of_light_normalized);
					Vector3d rotatedVelocity = particle->rotationTensor * (velocity * gamma);

					if(solverType == IMPLICIT){
						electricFlux[i][j][k] += rotatedVelocity * (particle->charge * particle->weight * correlation);
						//electricFlux[i] += velocity * (particle->charge * particle->weight * correlation);
						dielectricTensor[i][j][k] = dielectricTensor[i][j][k] - particle->rotationTensor * (particle->weight*theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge * correlation / particle->mass);
						//dielectricTensor[i] = dielectricTensor[i] + particle->rotationTensor * (particle->weight*theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge * correlation / particle->mass);
				
						Particle tempParticle = *particle;
						Matrix3d pressureTensorDerX;
						Matrix3d pressureTensorDerY;
						Matrix3d pressureTensorDerZ;
						Matrix3d tensor = rotatedVelocity.tensorMult(rotatedVelocity);

						double shiftX = 0.01*deltaX;
						if(particle->coordinates.x + shiftX >xgrid[xnumber]){
							shiftX = -shiftX;
						}
						tempParticle.coordinates.x = particle->coordinates.x + shiftX;

						double tempCorrelation = correlationWithEbin(tempParticle, i, j, k) / volumeE(i, j, k);

						pressureTensorDerX = tensor*particle->weight*particle->charge*(tempCorrelation - correlation)/shiftX;

						double shiftY = 0.01*deltaY;
						if(particle->coordinates.y + shiftY >ygrid[ynumber]){
							shiftY = -shiftY;
						}
						tempParticle.coordinates.x = particle->coordinates.x;
						tempParticle.coordinates.y = particle->coordinates.y + shiftY;

						tempCorrelation = correlationWithEbin(tempParticle, i, j, k) / volumeE(i, j, k);

						pressureTensorDerY = tensor*particle->weight*particle->charge*(tempCorrelation - correlation)/shiftY;

						double shiftZ = 0.01*deltaZ;
						if(particle->coordinates.z + shiftZ >zgrid[znumber]){
							shiftZ = -shiftZ;
						}
						tempParticle.coordinates.y = particle->coordinates.y;
						tempParticle.coordinates.z = particle->coordinates.z + shiftZ;

						tempCorrelation = correlationWithEbin(tempParticle, i, j, k) / volumeE(i, j, k);

						pressureTensorDerZ = tensor*particle->weight*particle->charge*(tempCorrelation - correlation)/shiftZ;

						divPressureTensor[i][j][k].x += pressureTensorDerX.matrix[0][0] + pressureTensorDerY.matrix[1][0] + pressureTensorDerZ.matrix[2][0];
						divPressureTensor[i][j][k].y += pressureTensorDerX.matrix[0][1] + pressureTensorDerY.matrix[1][1] + pressureTensorDerZ.matrix[2][1];
						divPressureTensor[i][j][k].z += pressureTensorDerX.matrix[0][2] + pressureTensorDerY.matrix[1][2] + pressureTensorDerZ.matrix[2][2];
					}
					if(solverType == EXPLICIT){
						electricFlux[i][j][k] += velocity*particle->charge*particle->weight*correlation;
					}
					alertNaNOrInfinity(electricFlux[i][j][k].x, "right part x = NaN");
					alertNaNOrInfinity(electricFlux[i][j][k].y, "right part y = NaN");
					alertNaNOrInfinity(electricFlux[i][j][k].z, "right part z = NaN");
				}
			}
		}
	}

	//for periodic conditions we must summ sides parameters
	if(boundaryConditionType == PERIODIC){
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				electricFlux[0][j][k] = electricFlux[0][j][k] + electricFlux[xnumber][j][k];
				electricFlux[xnumber][j][k] = electricFlux[0][j][k];

				dielectricTensor[0][j][k] = dielectricTensor[0][j][k] + dielectricTensor[xnumber][j][k];
				dielectricTensor[xnumber][j][k] = dielectricTensor[0][j][k];

				divPressureTensor[0][j][k] = divPressureTensor[0][j][k] + divPressureTensor[xnumber][j][k];
				divPressureTensor[xnumber][j][k] = divPressureTensor[0][j][k];
			}
		}
	}

	//for periodic y
	for(int i = 0; i < xnumber + 1; ++i){
		for(int k = 0; k < znumber + 1; ++k){
			electricFlux[i][0][k] = electricFlux[i][0][k] + electricFlux[i][ynumber][k];
			electricFlux[i][ynumber][k] = electricFlux[i][0][k];

			dielectricTensor[i][0][k] = dielectricTensor[i][0][k] + dielectricTensor[i][ynumber][k];
			dielectricTensor[i][ynumber][k] = dielectricTensor[i][0][k];

			divPressureTensor[i][0][k] = divPressureTensor[i][0][k] + divPressureTensor[i][ynumber][k];
			divPressureTensor[i][ynumber][k] = divPressureTensor[i][0][k];
		}
	}

	//for periodic z
	for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			electricFlux[i][j][0] = electricFlux[i][j][0] + electricFlux[i][j][znumber];
			electricFlux[i][j][znumber] = electricFlux[i][j][0];

			dielectricTensor[i][j][0] = dielectricTensor[i][j][0] + dielectricTensor[i][j][znumber];
			dielectricTensor[i][j][znumber] = dielectricTensor[i][0][0];

			divPressureTensor[i][j][0] = divPressureTensor[i][j][0] + divPressureTensor[i][j][znumber];
			divPressureTensor[i][j][znumber] = divPressureTensor[i][j][0];
		}
	}


	for (int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				electricDensity[i][j][k] = 0;
				pressureTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				if(solverType == IMPLICIT){
					for (int pcount = 0; pcount < particlesInBbin[i][j][k].size(); ++pcount) {
						Particle* particle = particlesInBbin[i][j][k][pcount];
						double correlation = correlationWithBbin(*particle, i, j, k) / volumeB(i, j, k);

						double gamma = particle->gammaFactor(speed_of_light_normalized);
						Vector3d velocity = particle->velocity(speed_of_light_normalized);
						Vector3d rotatedVelocity = particle->rotationTensor * velocity * gamma;

						electricDensity[i][j][k] += particle->weight * particle->charge * correlation;

						pressureTensor[i][j][k] += rotatedVelocity.tensorMult(rotatedVelocity) * particle->weight * particle->charge * correlation;
					}
				}
			}
		}
	}

	if(solverType == IMPLICIT){
		for (int i = 0; i < xnumber; ++i) {
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					Vector3d divPressureTensorEvaluated = evaluateDivPressureTensor(i, j, k);
					electricFlux[i][j][k] = electricFlux[i][j][k] - divPressureTensor[i][j][k] * eta * deltaT;
					//electricFlux[i] = electricFlux[i] - divPressureTensorEvaluated * eta * deltaT;
				}
			}
		}
	}

	//smoothFlux();

	if(solverType == IMPLICIT){
		for (int i = 0; i < xnumber; ++i) {
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					double divJ = evaluateDivFlux(i, j, k);

					electricDensity[i][j][k] -= deltaT * theta * divJ;
				}
			}
		}
	}

	updateExternalFlux();

	for(int i = 0; i < xnumber +1 ; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				electricFlux[i][j][k] = electricFlux[i][j][k] + externalElectricFlux[i][j][k];
			}
		}
	}

	//for debug only
		/*double kw = 2*pi/xsize;
		double concentration = density/(massProton + massElectron);
		for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
		for(int k = 0l k < znumber + 1; ++k){
			electricFlux[i][j][k].y = electron_charge_normalized*concentration*(VyamplitudeProton - VyamplitudeElectron)*sin(kw*xgrid[i] - omega*(time + theta*deltaT));
			electricFlux[i][j][k].z = electron_charge_normalized*concentration*(VzamplitudeProton - VzamplitudeElectron)*cos(kw*xgrid[i] - omega*(time + theta*deltaT));
		}
		}
		}
		*/
	//
}

void Simulation::updateDensityParameters() {
	double full_density = 0;
	double full_p_concentration = 0;
	double full_e_concentration = 0;
	collectParticlesIntoBins();
	//FILE* debugFile = fopen("./output/particleCorrelations.dat","w");
	for (int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				electronConcentration[i][j][k] = 0;
				protonConcentration[i][j][k] = 0;
				chargeDensity[i][j][k] = 0;
				velocityBulk[i][j][k] = Vector3d(0, 0, 0);
				velocityBulkElectron[i][j][k] = Vector3d(0, 0, 0);
				//fprintf(debugFile, "%d %d %d\n", i, j, k);
				for (int pcount = 0; pcount < particlesInBbin[i][j][k].size(); ++pcount) {
					Particle* particle = particlesInBbin[i][j][k][pcount];

					double correlation = correlationWithBbin(*particle, i, j, k) / volumeB(i, j, k);

					chargeDensity[i][j][k] += correlation * particle->charge * particle->weight;
					if (particle->type == ELECTRON) {
						electronConcentration[i][j][k] += correlation * particle->weight;
						velocityBulkElectron[i][j][k] += particle->momentum * particle->weight * correlation;
					} else if (particle->type == PROTON) {
						protonConcentration[i][j][k] += correlation * particle->weight;
					}
					velocityBulk[i][j][k] += particle->momentum * particle->weight * correlation;

					if (correlation == 0) {
						//printf("aaa\n");
					}

					//fprintf(debugFile, "%d %15.10g\n", particle->number, correlation*volume(i, j, k));
				}

				//fprintf(debugFile, "charge %15.10g proton %15.10g electron %15.10g\n", chargeDensity[i][j][k], protonConcentration[i][j][k], electronConcentration[i][j][k]);

 				velocityBulk[i][j][k] = velocityBulk[i][j][k] / (electronConcentration[i][j][k] * massElectron + protonConcentration[i][j][k] * massProton);
				velocityBulkElectron[i][j][k] = velocityBulkElectron[i][j][k] / (electronConcentration[i][j][k] * massElectron);
				full_density += chargeDensity[i][j][k] * volumeB(i, j, k);
				full_p_concentration += protonConcentration[i][j][k] * volumeB(i, j, k);
				full_e_concentration += electronConcentration[i][j][k] * volumeB(i, j, k);
			}
		}
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
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				Vector3d E = Efield[i][j][k]*fieldScale;
				electricFieldEnergy += E.scalarMult(E) * volumeE(i, j, k) / (8 * pi);
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				Vector3d B = (Bfield[i][j][k] - B0)*fieldScale;
				magneticFieldEnergy += B.scalarMult(B) * volumeB(i, j, k) / (8 * pi);
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				Vector3d E = ((Efield[i][j][k] + Efield[i][j][k+1] + Efield[i][j+1][k] + Efield[i][j+1][k+1] + Efield[i+1][j][k] + Efield[i+1][j][k+1] + Efield[i+1][j+1][k] + Efield[i + 1][j+1][k+1]) / 8)*fieldScale;
				Vector3d B = Bfield[i][j][k]*fieldScale;

				momentum += (E.vectorMult(B) / (4 * pi * speed_of_light_normalized)) * volumeB(i, j, k);
			}
		}
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

Vector3d Simulation::getBfield(int i, int j, int k) {
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

	if (j == -1) {
		j = ynumber - 1;
	} else if (j == ynumber) {
		j = 0;
	} else if (j < -1 || j > ynumber) {
		printf("j < -1 || j > ynumber in getBfied i = %d\n", j);
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "j = %d, ynumber = %d\n", j, ynumber);
		fclose(errorLogFile);
		exit(0);
	}

	if (k == -1) {
		k = znumber - 1;
	} else if (k == znumber) {
		k = 0;
	} else if (k < -1 || k > znumber) {
		printf("k < -1 || k > znumber in getBfied i = %d\n", k);
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "k = %d, znumber = %d\n", k, znumber);
		fclose(errorLogFile);
		exit(0);
	}

	return Bfield[i][j][k];
}

Vector3d Simulation::getTempEfield(int i, int j, int k) {
	if (i < 0) {
		i = xnumber - 1;
	} else if (i > xnumber) {
		i = 1;
	}

	if (j < 0) {
		j = ynumber - 1;
	} else if (j > ynumber) {
		j = 1;
	}

	if (k < 0) {
		k = znumber - 1;
	} else if (k > znumber) {
		k = 1;
	}

	return tempEfield[i][j][k];
}

Vector3d Simulation::getNewEfield(int i, int j, int k) {
	if (i < 0) {
		i = xnumber - 1;
	} else if (i > xnumber) {
		i = 1;
	}

	if (j < 0) {
		j = ynumber - 1;
	} else if (j > ynumber) {
		j = 1;
	}

	if (k < 0) {
		k = znumber - 1;
	} else if (k > znumber) {
		k = 1;
	}

	return newEfield[i][j][k];
}

Vector3d Simulation::getEfield(int i, int j, int k) {
	if (i < 0) {
		i = xnumber - 1;
	} else if (i > xnumber) {
		i = 1;
	}

	if (j < 0) {
		j = ynumber - 1;
	} else if (j > ynumber) {
		j = 1;
	}

	if (k < 0) {
		k = znumber - 1;
	} else if (k > znumber) {
		k = 1;
	}

	return Efield[i][j][k];
}

Matrix3d Simulation::getPressureTensor(int i, int j, int k) {
	if (i < 0) {
		i = xnumber - 1;
	} else if (i >= xnumber) {
		i = 0;
	}

	if (j < 0) {
		j = ynumber - 1;
	} else if (j >= ynumber) {
		j = 0;
	}

	if (k < 0) {
		k = znumber - 1;
	} else if (k >= znumber) {
		k = 0;
	}

	return pressureTensor[i][j][k];
}

double Simulation::getDensity(int i, int j, int k) {
	if (i < 0) {
		i = xnumber - 1;
	} else if (i >= xnumber) {
		i = 0;
	}

	if (j < 0) {
		j = ynumber - 1;
	} else if (j >= ynumber) {
		j = 0;
	}

	if (k < 0) {
		k = znumber - 1;
	} else if (k >= znumber) {
		k = 0;
	}


	return electricDensity[i][j][k];
}

void Simulation::updateParameters(){
	maxBfield = Bfield[0][0][0];
	maxEfield = Efield[0][0][0];
	for(int i = 0; i < xnumber ; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				if(Bfield[i][j][k].norm() > maxBfield.norm()){
					maxBfield = Bfield[i][j][k];
				}
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				if(Efield[i][j][k].norm() > maxEfield.norm()){
					maxEfield = Efield[i][j][k];
				}
			}
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
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				externalElectricFlux[i][j][k] = Vector3d(0, 0, 1.0)*extJ*cos(kw*xgrid[i] - omega*time);
			}
		}
	}
}

void Simulation::resetNewTempFields(){
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}
}