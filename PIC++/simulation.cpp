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
	if (newlyStarted) {
		createArrays();
		createFiles();
		initialize();
		//initializeTwoStream();
		//initializeExternalFluxInstability();
		//initializeAlfvenWave(1, 0.01);
		//initializeRotatedAlfvenWave(1, 0.01);
		initializeFluxFromRight();
		//initializeSimpleElectroMagneticWave();
		//initializeRotatedSimpleElectroMagneticWave(1);
		//initializeLangmuirWave();
		//createParticles();
	}
	//printf("E.x = %g E.y = %g E.z = %g\n", Efield[0][0][0].x, Efield[0][0][0].y, Efield[0][0][0].z);
	//printf("tempE.x = %g tempE.y = %g tempE.z = %g\n", tempEfield[0][0][0].x, tempEfield[0][0][0].y, tempEfield[0][0][0].z);
	//outputEverythingFile = fopen("./output/outputEverythingFile.dat", "w");
	collectParticlesIntoBins();
	//fclose(outputEverythingFile);
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
	cleanupDivergence();
	updateFields();
	updateEnergy();
	theoreticalEnergy = energy;
	theoreticalMomentum = momentum;

	//printf("E.x = %g E.y = %g E.z = %g\n", Efield[0][0][0].x, Efield[0][0][0].y, Efield[0][0][0].z);
	//printf("tempE.x = %g tempE.y = %g tempE.z = %g\n", tempEfield[0][0][0].x, tempEfield[0][0][0].y, tempEfield[0][0][0].z);

	for (int i = 0; i < typesNumber; ++i) {
		types[i].injectionLength = types[i].particesDeltaX - 0.0001 * deltaX;
	}

	while (time * plasma_period < maxTime && currentIteration < maxIteration) {
		printf("iteration number %d time = %15.10g\n", currentIteration, time * plasma_period);
		printLog("start iteration\n");
		printf(" dt/plasma_period = %15.10g\n", deltaT);

		if (currentIteration % writeParameter == 0) {
			output();
		}

		updateDeltaT();
		evaluateParticlesRotationTensor();
		updateElectroMagneticParameters();
		evaluateFields();
		evaluateMagneticField();

		//printf("E.x = %g E.y = %g E.z = %g\n", Efield[0][0][0].x, Efield[0][0][0].y, Efield[0][0][0].z);
		//printf("tempE.x = %g tempE.y = %g tempE.z = %g\n", tempEfield[0][0][0].x, tempEfield[0][0][0].y, tempEfield[0][0][0].z);
		moveParticles();

		for (int i = 0; i < typesNumber; ++i) {
			types[i].injectionLength += fabs(V0.x * deltaT);
		}
		if (boundaryConditionType == SUPER_CONDUCTOR_LEFT || boundaryConditionType == FREE_BOTH) {
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
		if(preserveChargeGlobal){
			addToPreserveChargeGlobal();
		}
		updateDensityParameters();
		cleanupDivergence();
		updateFields();
		updateEnergy();
		updateParameters();

		time += deltaT;
		currentIteration++;

		if (currentIteration % writeBackupParameter == 0) {
			outputBackup();
		}
	}
}

void Simulation::output() {
	printf("outputing\n");
	printLog("outputing\n");

	if (particles.size() > 0) {
		distributionFileProton = fopen((outputDir + "distribution_protons.dat").c_str(), "a");
		outputDistribution(distributionFileProton, particles, PROTON, scaleFactor, plasma_period);
		fclose(distributionFileProton);
		distributionFileElectron = fopen((outputDir + "distribution_electrons.dat").c_str(), "a");
		outputDistribution(distributionFileElectron, particles, ELECTRON, scaleFactor, plasma_period);
		fclose(distributionFileElectron);
		distributionFileProtonUpstream = fopen((outputDir + "distribution_protons_upstream.dat").c_str(), "a");
		outputDistributionUpstream(distributionFileProtonUpstream, particles, PROTON, xgrid[shockWavePoint], scaleFactor, plasma_period);
		fclose(distributionFileProtonUpstream);
		distributionFileElectronUpstream = fopen((outputDir + "distribution_electrons_upstream.dat").c_str(), "a");
		outputDistributionUpstream(distributionFileElectronUpstream, particles, ELECTRON, xgrid[shockWavePoint], scaleFactor, plasma_period);
		fclose(distributionFileElectronUpstream);
		protonTraectoryFile = fopen((outputDir + "trajectory_proton.dat").c_str(), "a");
		outputTrajectory(protonTraectoryFile, getFirstProton(), time, plasma_period, scaleFactor);
		fclose(protonTraectoryFile);
		electronTraectoryFile = fopen((outputDir + "trajectory_electron.dat").c_str(), "a");
		outputTrajectory(electronTraectoryFile, getFirstElectron(), time, plasma_period, scaleFactor);
		fclose(electronTraectoryFile);
	}

	EfieldFile = fopen((outputDir + "Efield.dat").c_str(), "a");
	BfieldFile = fopen((outputDir + "Bfield.dat").c_str(), "a");
	outputFields(EfieldFile, BfieldFile, Efield, Bfield, xnumber, ynumber, znumber, plasma_period, scaleFactor, fieldScale);
	fclose(EfieldFile);
	fclose(BfieldFile);

	Xfile = fopen((outputDir + "Xfile.dat").c_str(), "w");
	outputGrid(Xfile, xgrid, xnumber, scaleFactor);
	fclose(Xfile);

	Yfile = fopen((outputDir + "Yfile.dat").c_str(), "w");
	outputGrid(Yfile, ygrid, ynumber, scaleFactor);
	fclose(Yfile);

	Zfile = fopen((outputDir + "Zfile.dat").c_str(), "w");
	outputGrid(Zfile, zgrid, znumber, scaleFactor);
	fclose(Zfile);

	densityFile = fopen((outputDir + "concentrations.dat").c_str(), "a");
	outputConcentrations(densityFile, electronConcentration, protonConcentration, chargeDensity, electricDensity, xnumber, ynumber, znumber, plasma_period, scaleFactor, fieldScale);
	fclose(densityFile);

	velocityFile = fopen((outputDir + "velocity.dat").c_str(), "a");
	velocityElectronFile = fopen((outputDir + "velocity_electron.dat").c_str(), "a");
	outputVelocity(velocityFile, velocityElectronFile, velocityBulkProton, velocityBulkElectron, xnumber, ynumber, znumber, plasma_period, scaleFactor);
	fclose(velocityFile);
	fclose(velocityElectronFile);

	fluxFile = fopen((outputDir + "flux.dat").c_str(), "a");
	outputFlux(fluxFile, electricFlux, externalElectricFlux, xnumber + 1, ynumber + 1, znumber + 1, plasma_period, scaleFactor, fieldScale);
	fclose(fluxFile);

	divergenceErrorFile = fopen((outputDir + "divergence_error.dat").c_str(), "a");
	outputDivergenceError(divergenceErrorFile, this, plasma_period, scaleFactor, fieldScale);
	fclose(divergenceErrorFile);

	double rotBscale = fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor));

	rotBFile = fopen((outputDir + "rotBFile.dat").c_str(), "a");
	outputVectorArray(rotBFile, rotB, xnumber + 1, ynumber + 1, znumber + 1, rotBscale);
	fclose(rotBFile);

	EderivativeFile = fopen((outputDir + "EderivativeFile.dat").c_str(), "a");
	outputVectorArray(EderivativeFile, Ederivative, xnumber + 1, ynumber + 1, znumber + 1, rotBscale);
	fclose(EderivativeFile);

	dielectricTensorFile = fopen((outputDir + "dielectricTensorFile.dat").c_str(), "a");
	outputMatrixArray(dielectricTensorFile, dielectricTensor, xnumber + 1, ynumber + 1, znumber + 1);
	fclose(dielectricTensorFile);

	particleProtonsFile = fopen((outputDir + "protons.dat").c_str(), "w");
	particleElectronsFile = fopen((outputDir + "electrons.dat").c_str(), "w");
	particlePositronsFile = fopen((outputDir + ".positrons.dat").c_str(), "w");
	particleAlphaFile = fopen((outputDir + "alphas.dat").c_str(), "w");
	outputParticles(particleProtonsFile, particleElectronsFile, particlePositronsFile, particleAlphaFile, this);
	fclose(particleProtonsFile);
	fclose(particleElectronsFile);
	fclose(particlePositronsFile);
	fclose(particleAlphaFile);

	/*maxwellMatrixFile = fopen("./output/maxwellMatrixFile.dat", "w");
	outputMaxwellEquationMatrixFull(maxwellMatrixFile, maxwellEquationMatrix, xnumber, ynumber, znumber, maxwellEquationMatrixSize);
	fclose(maxwellMatrixFile);*/

	generalFile = fopen((outputDir + "general.dat").c_str(), "a");
	outputGeneral(generalFile, this);
	fclose(generalFile);
}

void Simulation::outputBackup() {
	printf("writing backup\n");
	printLog("writing backup\n");

	std::string backupDir = backupDirectory;
	FILE* backupGeneralFile = fopen((backupDir + "general.dat").c_str(), "r");
	FILE* backupEfieldFile = fopen((backupDir + "Efield.dat").c_str(), "r");
	FILE* backupBfieldFile = fopen((backupDir + "Bfield.dat").c_str(), "r");
	FILE* backupParticlesFile = fopen((backupDir + "particles.dat").c_str(), "r");

	outputSimulationBackup(backupGeneralFile, backupEfieldFile, backupBfieldFile, backupParticlesFile, this);

	fclose(backupGeneralFile);
	fclose(backupEfieldFile);
	fclose(backupBfieldFile);
	fclose(backupParticlesFile);
}

void Simulation::updateDeltaT() {
	printf("updating time step\n");
	printLog("updating time step\n");
	//double delta = min2(deltaX, min2(deltaY, deltaZ));
	double delta = deltaX;
	deltaT = timeEpsilon * delta / speed_of_light_normalized;
	printLog("dx/c\n");
	if (particles.size() > 0) {
		double B = B0.norm();
		double E = E0.norm();
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					if (Bfield[i][j][k].norm() > B) {
						B = Bfield[i][j][k].norm();
					}
				}
			}
		}
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					if (Efield[i][j][k].norm() > E) {
						E = Efield[i][j][k].norm();
					}
				}
			}
		}

		B *= fieldScale;
		E *= fieldScale;

		double minFlux = electricFlux[0][0][0].norm();
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					if (electricFlux[i][j][k].norm() < minFlux) {
						minFlux = electricFlux[i][j][k].norm();
					}
				}
			}
		}
		printLog("evaluatd E and minFlux\n");

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

		deltaT = min2(deltaT, timeEpsilon / omegaPlasmaElectron);

		omegaGyroProton = electron_charge_normalized * B / (massProton * speed_of_light);
		omegaGyroElectron = electron_charge_normalized * B / (massElectron * speed_of_light);

		double thermalMomentum = sqrt(massElectron * kBoltzman_normalized * temperature) + massElectron * V0.norm();

		printLog("evaluzted omega\n");

		if (B > 0) {
			deltaT = min2(deltaT, timeEpsilon * massElectron * speed_of_light_normalized / (electron_charge_normalized * B));
		}
		if (E > 0) {
			deltaT = min2(deltaT, timeEpsilon * massElectron * speed_of_light_normalized / (electron_charge_normalized * E));
		}
		/*if (E > 0) {
			deltaT = min2(deltaT, 0.1 * thermalMomentum / (electron_charge_normalized * E));
		}*/


		double Vthermal = sqrt(2 * kBoltzman_normalized * temperature / massElectron);
		double minDeltaT = deltaX / Vthermal;
		/*if(minDeltaT > deltaT){
			printf("deltaT < dx/Vthermal\n");
		}*/
		printLog("end updating time step\n");
	}
}

double Simulation::volumeE(int i, int j, int k) {
	double dx = deltaX;
	if ((boundaryConditionType == PERIODIC) || ((i > 0) && (i < xnumber))) {
		dx = deltaX;
	} else {
		dx = deltaX / 2;
	}
	return dx * deltaY * deltaZ;
}

double Simulation::volumeB(int i, int j, int k) {
	return deltaX * deltaY * deltaZ;
}

void Simulation::checkParticleInBox(Particle& particle) {
	if (particle.coordinates.x < xgrid[0]) {
		printf("particle.x < 0\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.x = %15.10g < 0\n", particle.coordinates.x);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x, particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.velocity(speed_of_light_normalized).norm()/speed_of_light_normalized));
		fclose(errorLogFile);
		exit(0);
	}
	if (particle.coordinates.x > xgrid[xnumber]) {
		printf("particle.x > xsize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.x = %15.10g > %15.10g\n", particle.coordinates.x, xgrid[xnumber]);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x, particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.velocity(speed_of_light_normalized).norm()/speed_of_light_normalized));
		fclose(errorLogFile);
		exit(0);
	}

	if (particle.coordinates.y < ygrid[0]) {
		printf("particle.y < ygrid[0]\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.y = %15.10g < ygrid[0]\n", particle.coordinates.y);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x, particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.velocity(speed_of_light_normalized).norm()/speed_of_light_normalized));
		fclose(errorLogFile);
		exit(0);
	}
	if (particle.coordinates.y > ygrid[ynumber]) {
		printf("particle.y > ysize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.y = %15.10g > %15.10g\n", particle.coordinates.y, ygrid[ynumber]);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x, particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.velocity(speed_of_light_normalized).norm()/speed_of_light_normalized));
		fclose(errorLogFile);
		exit(0);
	}

	if (particle.coordinates.z < zgrid[0]) {
		printf("particle.z < zgrid[0]\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.z = %15.10g < zgrid[0]\n", particle.coordinates.z);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x, particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.velocity(speed_of_light_normalized).norm()/speed_of_light_normalized));
		fclose(errorLogFile);
		exit(0);
	}
	if (particle.coordinates.z > zgrid[znumber]) {
		printf("particle.z > zsize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.z = %15.10g > %15.10g\n", particle.coordinates.z, zgrid[znumber]);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x, particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.velocity(speed_of_light_normalized).norm()/speed_of_light_normalized));
		fclose(errorLogFile);
		exit(0);
	}
}

void Simulation::checkParticlesInBin() {
	for (int j = 0; j < ynumber; ++j) {
		for (int k = 0; k < znumber; ++k) {
			for (int pcount = 0; pcount < particlesInEbin[0][j][k].size(); ++pcount) {
				Particle* particle = particlesInEbin[0][j][k][pcount];
				for (int pcountRight = 0; pcountRight < particlesInEbin[xnumber][j][k].size(); ++pcountRight) {
					if (particle == particlesInEbin[xnumber][j][k][pcountRight]) {
						printf("particle is in 0 and xnumber bin\n");
						errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
						fprintf(errorLogFile, "particle is in 0 and xnumber bin, particle.x = %15.10g xmax = %15.10g\n", particle->coordinates.x, xgrid[xnumber]);
						fclose(errorLogFile);
						exit(0);
					}
				}
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int k = 0; k < znumber; ++k) {
			for (int pcount = 0; pcount < particlesInEbin[i][0][k].size(); ++pcount) {
				Particle* particle = particlesInEbin[i][0][k][pcount];
				for (int pcountRight = 0; pcountRight < particlesInEbin[i][ynumber][k].size(); ++pcountRight) {
					if (particle == particlesInEbin[i][ynumber][k][pcountRight]) {
						printf("particle is in 0 and ynumber bin\n");
						errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
						fprintf(errorLogFile, "particle is in 0 and ynumber bin, particle.y = %15.10g ymax = %15.10g\n", particle->coordinates.y, ygrid[ynumber]);
						fclose(errorLogFile);
						exit(0);
					}
				}
			}
		}
	}

	for (int j = 0; j < ynumber; ++j) {
		for (int i = 0; i < xnumber; ++i) {
			for (int pcount = 0; pcount < particlesInEbin[i][j][0].size(); ++pcount) {
				Particle* particle = particlesInEbin[i][j][0][pcount];
				for (int pcountRight = 0; pcountRight < particlesInEbin[i][j][znumber].size(); ++pcountRight) {
					if (particle == particlesInEbin[i][j][znumber][pcountRight]) {
						printf("particle is in 0 and znumber bin\n");
						errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
						fprintf(errorLogFile, "particle is in 0 and znumber bin, particle.z = %15.10g zmax = %15.10g\n", particle->coordinates.z, zgrid[znumber]);
						fclose(errorLogFile);
						exit(0);
					}
				}
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				if (particlesInEbin[i][j][k].size() > 1) {
					for (int pcount = 0; pcount < particlesInEbin[i][j][k].size() - 1; ++pcount) {
						Particle* particle = particlesInEbin[i][j][k][pcount];
						for (int pcount2 = pcount + 1; pcount2 < particlesInEbin[i][j][k].size(); ++pcount2) {
							if (particle == particlesInEbin[i][j][k][pcount2]) {
								printf("particle is twice in Ebin number %d %d %d\n", i, j, k);
								errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
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

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				if (particlesInBbin[i][j][k].size() > 1) {
					for (int pcount = 0; pcount < particlesInBbin[i][j][k].size() - 1; ++pcount) {
						Particle* particle = particlesInBbin[i][j][k][pcount];
						for (int pcount2 = pcount + 1; pcount2 < particlesInBbin[i][j][k].size(); ++pcount2) {
							if (particle == particlesInBbin[i][j][k][pcount2]) {
								printf("particle is twice in Bbin number %d %d %d\n", i, j, k);
								errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
								fprintf(errorLogFile, "particle is twice in Bbi number %d %d %d\n", i, j, k);
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
	printLog("updating flux, density and dielectric tensor\n");
	//check particle only in one boundary box
	if (debugMode) {
		//checkParticlesInBin();
	}
	//collectParticlesIntoBins();
	//fopen("./output/outputEverythingFile.dat","a");
	int particlePartsCount = 0;
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				electricFlux[i][j][k] = Vector3d(0, 0, 0);
				dielectricTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				divPressureTensor[i][j][k] = Vector3d(0, 0, 0);
				for (int pcount = 0; pcount < particlesInEbin[i][j][k].size(); ++pcount) {
					particlePartsCount++;
					Particle* particle = particlesInEbin[i][j][k][pcount];
					double correlation = correlationWithEbin(*particle, i, j, k) / volumeE(i, j, k);

					Vector3d velocity = particle->velocity(speed_of_light_normalized);
					double gamma = particle->gammaFactor(speed_of_light_normalized);
					Vector3d rotatedVelocity = particle->rotationTensor * (velocity * gamma);

					if(i == debugPoint){
						//fprintf(outputEverythingFile, "particle number %d correlation = %28.22g\n", particle->number, correlation);
						//fprintf(outputEverythingFile, "%28.22g %28.22g %28.22g\n", rotatedVelocity.x, rotatedVelocity.y, rotatedVelocity.z);
					}

					if (solverType == IMPLICIT) {
						electricFlux[i][j][k] += rotatedVelocity * (particle->charge * particle->weight * correlation);
						//electricFlux[i][j][k] += velocity * (particle->charge * particle->weight * correlation);
						dielectricTensor[i][j][k] = dielectricTensor[i][j][k] - particle->rotationTensor * (particle->weight * theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge * correlation / particle->mass);
						//dielectricTensor[i] = dielectricTensor[i] + particle->rotationTensor * (particle->weight*theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge * correlation / particle->mass);
						if(i == debugPoint){
							//fprintf(outputEverythingFile, "electricFlux %d = %28.22g %28.22g %28.22g\n", debugPoint, electricFlux[i][j][k].x, electricFlux[i][j][k].y, electricFlux[i][j][k].z);
						}
						Particle tempParticle = *particle;
						Matrix3d pressureTensorDerX;
						Matrix3d pressureTensorDerY;
						Matrix3d pressureTensorDerZ;
						Matrix3d tensor = rotatedVelocity.tensorMult(rotatedVelocity);

						double shiftX = 0.01 * deltaX;
						if (particle->coordinates.x + shiftX > xgrid[xnumber]) {
							shiftX = -shiftX;
						}
						tempParticle.coordinates.x = particle->coordinates.x + shiftX;

						double tempCorrelation = correlationWithEbin(tempParticle, i, j, k) / volumeE(i, j, k);

						pressureTensorDerX = tensor * particle->weight * particle->charge * (tempCorrelation - correlation) / shiftX;

						double shiftY = 0.01 * deltaY;
						if (particle->coordinates.y + shiftY > ygrid[ynumber]) {
							shiftY = -shiftY;
						}
						tempParticle.coordinates.x = particle->coordinates.x;
						tempParticle.coordinates.y = particle->coordinates.y + shiftY;

						tempCorrelation = correlationWithEbin(tempParticle, i, j, k) / volumeE(i, j, k);

						pressureTensorDerY = tensor * particle->weight * particle->charge * (tempCorrelation - correlation) / shiftY;

						double shiftZ = 0.01 * deltaZ;
						if (particle->coordinates.z + shiftZ > zgrid[znumber]) {
							shiftZ = -shiftZ;
						}
						tempParticle.coordinates.y = particle->coordinates.y;
						tempParticle.coordinates.z = particle->coordinates.z + shiftZ;

						tempCorrelation = correlationWithEbin(tempParticle, i, j, k) / volumeE(i, j, k);

						pressureTensorDerZ = tensor * particle->weight * particle->charge * (tempCorrelation - correlation) / shiftZ;

						divPressureTensor[i][j][k].x += pressureTensorDerX.matrix[0][0] + pressureTensorDerY.matrix[1][0] + pressureTensorDerZ.matrix[2][0];
						divPressureTensor[i][j][k].y += pressureTensorDerX.matrix[0][1] + pressureTensorDerY.matrix[1][1] + pressureTensorDerZ.matrix[2][1];
						divPressureTensor[i][j][k].z += pressureTensorDerX.matrix[0][2] + pressureTensorDerY.matrix[1][2] + pressureTensorDerZ.matrix[2][2];
					}
					if (solverType == EXPLICIT) {
						electricFlux[i][j][k] += velocity * particle->charge * particle->weight * correlation;
					}
					if (i == 0 && boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
						if (particle->coordinates.x - particle->dx < 0) {
							addReflectedParticleToElectroMagneticParameters(particle, j, k);
						}
					}
					alertNaNOrInfinity(electricFlux[i][j][k].x, "right part x = NaN");
					alertNaNOrInfinity(electricFlux[i][j][k].y, "right part y = NaN");
					alertNaNOrInfinity(electricFlux[i][j][k].z, "right part z = NaN");
				}
			}
		}
	}

	//for periodic conditions we must summ sides parameters
	if (boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
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
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int k = 0; k < znumber + 1; ++k) {
			electricFlux[i][0][k] = electricFlux[i][0][k] + electricFlux[i][ynumber][k];
			electricFlux[i][ynumber][k] = electricFlux[i][0][k];

			dielectricTensor[i][0][k] = dielectricTensor[i][0][k] + dielectricTensor[i][ynumber][k];
			dielectricTensor[i][ynumber][k] = dielectricTensor[i][0][k];

			divPressureTensor[i][0][k] = divPressureTensor[i][0][k] + divPressureTensor[i][ynumber][k];
			divPressureTensor[i][ynumber][k] = divPressureTensor[i][0][k];
		}
	}

	//for periodic z
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			electricFlux[i][j][0] = electricFlux[i][j][0] + electricFlux[i][j][znumber];
			electricFlux[i][j][znumber] = electricFlux[i][j][0];

			dielectricTensor[i][j][0] = dielectricTensor[i][j][0] + dielectricTensor[i][j][znumber];
			dielectricTensor[i][j][znumber] = dielectricTensor[i][0][0];

			divPressureTensor[i][j][0] = divPressureTensor[i][j][0] + divPressureTensor[i][j][znumber];
			divPressureTensor[i][j][znumber] = divPressureTensor[i][j][0];
		}
	}

	//fprintf(outputEverythingFile, "electricFlux %d after boundaries = %28.22g %28.22g %28.22g\n", debugPoint, electricFlux[debugPoint][0][0].x, electricFlux[debugPoint][0][0].y, electricFlux[debugPoint][0][0].z);


	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				electricDensity[i][j][k] = 0;
				pressureTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				if (solverType == IMPLICIT) {
					if(i == debugPoint){
						//fprintf(outputEverythingFile, "particles in bin %d\n", particlesInBbin[i][j][k].size());
					}
					for (int pcount = 0; pcount < particlesInBbin[i][j][k].size(); ++pcount) {
						Particle* particle = particlesInBbin[i][j][k][pcount];
						double correlation = correlationWithBbin(*particle, i, j, k) / volumeB(i, j, k);
						if(i == debugPoint){
							//fprintf(outputEverythingFile, "particle number %d correlation = %28.22g\n", particle->number, correlation);
							//fprintf(outputEverythingFile, "%28.22g %28.22g %28.22g\n", particle->coordinates.x, particle->coordinates.y, particle->coordinates.z);
						}

						double gamma = particle->gammaFactor(speed_of_light_normalized);
						Vector3d velocity = particle->velocity(speed_of_light_normalized);
						Vector3d rotatedVelocity = particle->rotationTensor * velocity * gamma;

						electricDensity[i][j][k] += particle->weight * particle->charge * correlation;

						if(i == debugPoint){
							//fprintf(outputEverythingFile, "density %d = %28.22g\n", debugPoint, electricDensity[debugPoint][0][0]);
						}

						pressureTensor[i][j][k] += rotatedVelocity.tensorMult(rotatedVelocity) * particle->weight * particle->charge * correlation;
					}
				}
			}
		}
	}
	//fprintf(outputEverythingFile, "density %d = %28.22g\n", debugPoint - 1, electricDensity[debugPoint - 1][0][0]);
	//fprintf(outputEverythingFile, "density %d = %28.22g\n", debugPoint, electricDensity[debugPoint][0][0]);
	//fprintf(outputEverythingFile, "density %d = %28.22g\n", debugPoint + 1, electricDensity[debugPoint + 1][0][0]);

	if (solverType == IMPLICIT) {
		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber + 1; ++j) {
				for (int k = 0; k < znumber + 1; ++k) {
					//Vector3d divPressureTensorEvaluated = evaluateDivPressureTensor(i, j, k);
					electricFlux[i][j][k] = electricFlux[i][j][k] - divPressureTensor[i][j][k] * eta * deltaT;
					//electricFlux[i] = electricFlux[i] - divPressureTensorEvaluated * eta * deltaT;
				}
			}
		}
	}

	//fprintf(outputEverythingFile, "electricFlux %d after div tensor = %28.22g %28.22g %28.22g\n", debugPoint, electricFlux[debugPoint][0][0].x, electricFlux[debugPoint][0][0].y, electricFlux[debugPoint][0][0].z);

	//smoothFlux();

	FILE* tempFluxFile = fopen((outputDir + "tempFluxFile.dat").c_str(), "w");
	for(int i = 0; i < xnumber; ++i){
		fprintf(tempFluxFile, "%28.22g %28.22g %28.22g\n", electricFlux[i][0][0][0], electricFlux[i][0][0][1], electricFlux[i][0][0][2]);
	}
	fclose(tempFluxFile);

	if (solverType == IMPLICIT) {
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					double divJ = evaluateDivFlux(i, j, k);

					electricDensity[i][j][k] -= deltaT * theta * divJ;
					//electricDensity[i][j][k] = 0;
				}
			}
		}
	}
	//fprintf(outputEverythingFile, "density %d afterFlux = %28.22g\n", electricDensity[debugPoint - 1][0][0]);
	//fprintf(outputEverythingFile, "density %d afterFlux = %28.22g\n", electricDensity[debugPoint][0][0]);
	//fprintf(outputEverythingFile, "density %d afterFlux = %28.22g\n", electricDensity[debugPoint + 1][0][0]);

	FILE* tempDensityFile = fopen((outputDir + "tempDensityFile.dat").c_str(), "w");
	for(int i = 0; i < xnumber; ++i){
		fprintf(tempDensityFile, "%28.22g\n", electricDensity[i][0][0]);
	}
	fclose(tempDensityFile);

	updateExternalFlux();

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				electricFlux[i][j][k] = electricFlux[i][j][k] + externalElectricFlux[i][j][k];
			}
		}
	}

	//for debug only
	/*double kw = 2 * pi / xsize;
	double concentration = density / (massProton + massElectron);
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0l k<znumber + 1; ++k) {
				electricFlux[i][j][k].y = electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * sin(kw * xgrid[i] - omega * (time + theta * deltaT));
				electricFlux[i][j][k].z = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * cos(kw * xgrid[i] - omega * (time + theta * deltaT));
			}
		}
	}*/

	/*double kx = 1 * 2 * pi / xsize;
	double ky = 1 * 2 * pi / ysize;
	double kz = 1 * 2 * pi / zsize;
	kz = 0;
	double concentration = density / (massProton + massElectron);

	double kw = sqrt(kx * kx + ky * ky + kz * kz);

	double kxy = sqrt(kx * kx + ky * ky);
	double rotatedZortNorm = sqrt(kx * kx + ky * ky + sqr(kx * kx + ky * ky) / (kz * kz));
	double matrixzz = (kx * kx + ky * ky) / (kz * rotatedZortNorm);

	Matrix3d rotationMatrix = Matrix3d(kx / kw, -ky / kxy, -kx / rotatedZortNorm,
	                                   ky / kw, kx / kxy, -ky / rotatedZortNorm,
	                                   kz / kw, 0, matrixzz);

	if (kz == 0) {
		rotationMatrix = Matrix3d(kx / kw, -ky / kw, 0,
		                          ky / kw, kx / kw, 0,
		                          0, 0, 1);
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Vector3d newFlux;

				newFlux.x = 0;
				newFlux.y = electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * sin(kx * xgrid[i] + ky * ygrid[j] + kz * zgrid[k] - omega * (time + theta * deltaT));
				newFlux.z = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * cos(kx * xgrid[i] + ky * ygrid[j] + kz * zgrid[k] - omega * (time + theta * deltaT));

				electricFlux[i][j][k] = rotationMatrix * newFlux;
			}
		}
	}*/

	//
	//fclose(outputEverythingFile);
}

void Simulation::smoothDensity() {
	for (int j = 0; j < ynumber; ++j) {
		for (int k = 0; k < znumber; ++k) {
			double newLeftDensity = (electricDensity[0][j][k] + electricDensity[1][j][k]) / 2.0;
			double newRightDensity = (electricDensity[xnumber - 1][j][k] + electricDensity[xnumber - 2][j][k]) / 2.0;
			if (boundaryConditionType == PERIODIC) {
				newLeftDensity = (electricDensity[xnumber - 1][j][k] + 2.0 * electricDensity[0][j][k] + electricDensity[1][j][k]) / 4.0;
				newRightDensity = (electricDensity[xnumber - 2][j][k] + 2.0 * electricDensity[xnumber - 1][j][k] + electricDensity[0][j][k]) / 4.0;
			}
			double prevDensity = electricDensity[0][j][k];
			for (int i = 1; i < xnumber - 1; ++i) {
				double tempDensity = electricDensity[i][j][k];

				electricDensity[i][j][k] = (prevDensity + 2.0 * electricDensity[i][j][k] + electricDensity[i + 1][j][k]) / 4.0;

				prevDensity = tempDensity;
			}

			electricDensity[0][j][k] = newLeftDensity;
			electricDensity[xnumber - 1][j][k] = newRightDensity;
		}
	}
}

void Simulation::addReflectedParticleToElectroMagneticParameters(const Particle* particle, int j, int k) {
	Particle tempParticle = *particle;

	tempParticle.momentum.x = - particle->momentum.x;

	double correlation = correlationWithEbin(tempParticle, -1, j, k) / volumeE(0, j, k);

	double beta = 0.5 * particle->charge * deltaT / particle->mass;
	Vector3d velocity = tempParticle.velocity(speed_of_light_normalized);

	Vector3d oldE = correlationEfield(tempParticle) * fieldScale;
	Vector3d oldB = correlationBfield(tempParticle) * fieldScale;

	tempParticle.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, oldE, oldB);
	double gamma = tempParticle.gammaFactor(speed_of_light_normalized);

	Vector3d rotatedVelocity = tempParticle.rotationTensor * (velocity * gamma);

	if (solverType == IMPLICIT) {
		electricFlux[0][j][k] += rotatedVelocity * (particle->charge * particle->weight * correlation);
		//electricFlux[i] += velocity * (particle->charge * particle->weight * correlation);
		dielectricTensor[0][j][k] = dielectricTensor[0][j][k] - tempParticle.rotationTensor * (particle->weight * theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge * correlation / particle->mass);
		//dielectricTensor[0] = dielectricTensor[i] + particle->rotationTensor * (particle->weight*theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge * correlation / particle->mass);

		Particle tempParticle2 = tempParticle;
		double shiftX = 0.01 * deltaX;
		if (tempParticle.coordinates.x + shiftX > xgrid[xnumber]) {
			shiftX = -shiftX;
		}
		tempParticle2.coordinates.x += shiftX;

		double tempCorrelation = correlationWithEbin(tempParticle2, -1, j, k) / volumeE(0, j, k);

		divPressureTensor[0][j][k].x += (rotatedVelocity.tensorMult(rotatedVelocity)).matrix[0][0] * particle->weight * particle->charge * (tempCorrelation - correlation) / shiftX;
		divPressureTensor[0][j][k].y += (rotatedVelocity.tensorMult(rotatedVelocity)).matrix[0][1] * particle->weight * particle->charge * (tempCorrelation - correlation) / shiftX;
		divPressureTensor[0][j][k].z += (rotatedVelocity.tensorMult(rotatedVelocity)).matrix[0][2] * particle->weight * particle->charge * (tempCorrelation - correlation) / shiftX;
	}
	if (solverType == EXPLICIT) {
		electricFlux[0][j][k] += velocity * particle->charge * particle->weight * correlation;
	}
}

void Simulation::updateDensityParameters() {
	double full_density = 0;
	double full_p_concentration = 0;
	double full_e_concentration = 0;
	collectParticlesIntoBins();
	//FILE* debugFile = fopen("./output/particleCorrelations.dat","w");
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				electronConcentration[i][j][k] = 0;
				protonConcentration[i][j][k] = 0;
				chargeDensity[i][j][k] = 0;
				velocityBulkProton[i][j][k] = Vector3d(0, 0, 0);
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
						velocityBulkProton[i][j][k] += particle->momentum * particle->weight * correlation;
					}

					if (correlation == 0) {
						//printf("aaa\n");
					}

					//fprintf(debugFile, "%d %15.10g\n", particle->number, correlation*volume(i, j, k));
				}

				//fprintf(debugFile, "charge %15.10g proton %15.10g electron %15.10g\n", chargeDensity[i][j][k], protonConcentration[i][j][k], electronConcentration[i][j][k]);

				velocityBulkProton[i][j][k] = velocityBulkProton[i][j][k] / (protonConcentration[i][j][k] * massProton);
				double gamma = sqrt((velocityBulkProton[i][j][k].scalarMult(velocityBulkProton[i][j][k]) / speed_of_light_normalized_sqr) + 1);
				velocityBulkProton[i][j][k] = velocityBulkProton[i][j][k] / gamma;
				velocityBulkElectron[i][j][k] = velocityBulkElectron[i][j][k] / (electronConcentration[i][j][k] * massElectron);
				gamma = sqrt((velocityBulkElectron[i][j][k].scalarMult(velocityBulkElectron[i][j][k]) / speed_of_light_normalized_sqr) + 1);
				velocityBulkElectron[i][j][k] = velocityBulkElectron[i][j][k] / gamma;
				full_density += chargeDensity[i][j][k] * volumeB(i, j, k);
				full_p_concentration += protonConcentration[i][j][k] * volumeB(i, j, k);
				full_e_concentration += electronConcentration[i][j][k] * volumeB(i, j, k);
			}
		}
	}
	full_density /= xsize;
	full_p_concentration /= (xsize * scaleFactor);
	full_e_concentration /= (xsize * scaleFactor);
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
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Vector3d E = Efield[i][j][k] * fieldScale;
				electricFieldEnergy += E.scalarMult(E) * volumeE(i, j, k) / (8 * pi);
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Vector3d B = Bfield[i][j][k] * fieldScale;
				magneticFieldEnergy += B.scalarMult(B) * volumeB(i, j, k) / (8 * pi);
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Vector3d E = ((Efield[i][j][k] + Efield[i][j][k + 1] + Efield[i][j + 1][k] + Efield[i][j + 1][k + 1] + Efield[i + 1][j][k] + Efield[i + 1][j][k + 1] + Efield[i + 1][j + 1][k] + Efield[i + 1][j + 1][k + 1]) / 8) * fieldScale;
				Vector3d B = Bfield[i][j][k] * fieldScale;

				momentum += (E.vectorMult(B) / (4 * pi * speed_of_light_normalized)) * volumeB(i, j, k);
			}
		}
	}

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		particleEnergy += particles[pcount]->energy(speed_of_light_normalized) * particles[pcount]->weight;
		momentum += particles[pcount]->momentum * particles[pcount]->weight;
	}


	particleEnergy *= sqr(scaleFactor / plasma_period);
	electricFieldEnergy *= sqr(scaleFactor / plasma_period);
	magneticFieldEnergy *= sqr(scaleFactor / plasma_period);
	momentum = momentum * scaleFactor / plasma_period;


	energy = particleEnergy + electricFieldEnergy + magneticFieldEnergy;

	double concentration = density / (massProton + massElectron);

	if (boundaryConditionType != PERIODIC) {
		//theoreticalEnergy -= (density * V0.scalarMult(V0) / 2.0) * V0.x * deltaT * ysize * zsize * sqr(scaleFactor / plasma_period);
		//theoreticalEnergy -= (2 * (3.0 / 2.0) * concentration * kBoltzman_normalized * temperature) * V0.x * deltaT * ysize * zsize * sqr(scaleFactor / plasma_period);

		//theoreticalMomentum -= V0 * V0.x * density * deltaT * ysize * zsize * scaleFactor / plasma_period;
		//theoreticalMomentum -= Vector3d(1, 0, 0) * (2 * concentration * kBoltzman_normalized * temperature) * deltaT * ysize * zsize * scaleFactor / plasma_period;
		//it is taken in account in inject new particle

		for (int i = 0; i < escapedParticles.size(); ++i) {
			Particle* particle = escapedParticles[i];
			theoreticalEnergy -= particle->energy(speed_of_light_normalized) * particle->weight * sqr(scaleFactor / plasma_period);
			theoreticalMomentum -= particle->momentum * particle->weight * scaleFactor / plasma_period;
		}

		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				theoreticalEnergy -= (Efield[xnumber][j][k].vectorMult(Bfield[xnumber - 1][j][k])).x * deltaT * speed_of_light_normalized * deltaZ * deltaY * sqr(scaleFactor / plasma_period)/ (4 * pi);
				theoreticalEnergy += (Efield[0][j][k].vectorMult(Bfield[0][j][k])).x * deltaT * speed_of_light_normalized * deltaZ * deltaY * sqr(scaleFactor / plasma_period)/ (4 * pi);

				theoreticalMomentum -= (Efield[xnumber][j][k].vectorMult(Bfield[xnumber - 1][j][k])) * deltaT * deltaZ * deltaY * scaleFactor / plasma_period / (4 * pi);
				theoreticalMomentum += (Efield[0][j][k].vectorMult(Bfield[0][j][k])) * deltaT * deltaZ * deltaY * scaleFactor / plasma_period / (4 * pi);
			}
		}
	}
}

Vector3d Simulation::getBfield(int i, int j, int k) {
	if (i == -1) {
		i = xnumber - 1;
	} else if (i == xnumber) {
		i = 0;
	} else if (i < -1 || i > xnumber) {
		printf("i < -1 || i > xnumber in getBfied i = %d\n", i);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
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
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
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
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
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

void Simulation::updateParameters() {
	maxBfield = Bfield[0][0][0];
	maxEfield = Efield[0][0][0];
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				if (Bfield[i][j][k].norm() > maxBfield.norm()) {
					maxBfield = Bfield[i][j][k];
				}
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				if (Efield[i][j][k].norm() > maxEfield.norm()) {
					maxEfield = Efield[i][j][k];
				}
			}
		}
	}
}

void Simulation::updateExternalFlux() {
	double alfvenV = B0.norm() * fieldScale / sqrt(4 * pi * density);
	double concentration = density / (massProton + massElectron);
	double phaseV = 2 * alfvenV;
	double kw = 2 * pi / xsize;
	double omega = kw * phaseV;

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				externalElectricFlux[i][j][k] = Vector3d(0, 0, 1.0) * extJ * cos(kw * xgrid[i] - omega * time);
			}
		}
	}
}

void Simulation::resetNewTempFields() {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}
}
