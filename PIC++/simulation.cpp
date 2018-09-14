#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include <time.h>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "mpi_util.h"
#include "util.h"
#include "output.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "particle.h"
#include "random.h"
#include "simulation.h"
#include "paths.h"

double Simulation::massProton;
double Simulation::massElectron;
double Simulation::massAlpha;
double Simulation::massDeuterium;
double Simulation::massHelium3;
double Simulation::massOxygen;
double Simulation::massSilicon;

//double Simulation::speed_of_light_normalized;
//double Simulation::speed_of_light_normalized_sqr;
double Simulation::kBoltzman_normalized;
double Simulation::electron_charge_normalized;
DimensionType Simulation::dimensionType;

void Simulation::simulate() {
	MPI_Barrier(cartComm);
	double initializationTime = 0;
	MPI_Barrier(cartComm);
	if (rank == 0) {
		initializationTime = clock();
	}
	if (newlyStarted) {
		//if(rank == 0) printf("create arrays\n");
		createArrays();
		initialize();
		updateDeltaT();
		createFiles();
		//initializeTwoStream();
		//initializeExternalFluxInstability();
		//initializeAlfvenWaveX(1, 0.01);
		//initializeAlfvenWaveY(1, 0.01);
		//initializeAlfvenWaveZ(1, 0.01);
		//initializeRotatedAlfvenWave(1, 1, 0, 0.01);
		//initializeAnisotropic();
		//initializeAnisotropicSilicon();
		//initializeWeibel();
		//initializeRingWeibel();
		initializeFluxFromRight();
		//initializeHarris();
		//initializeBell();
		//initializeSimpleElectroMagneticWave();
		//initializeSimpleElectroMagneticWaveY();
		//initializeSimpleElectroMagneticWaveZ();
		//initializeRotatedSimpleElectroMagneticWave(1, 1, 0);
		//initializeHomogenouseFlow();
		//initializeLangmuirWave();
		//createParticles();
		//initializeShockWave();
		//initializeFake();
		//initializeTestOneParticle();
	}
	if (rank == 0) {
		initializationTime = clock() - initializationTime;
		printf("initialization time time = %g sec\n", initializationTime / CLOCKS_PER_SEC);
	}

	outputGeneralInitialParameters((outputDir + "initialParameters.dat").c_str(),
	                               (outputDir + "initialParametersWithText.dat").c_str(), this);


	MPI_Barrier(cartComm);

	if ((rank == 0) && (verbosity > 1)) printf("initial exchanging fields\n");
	if(solverType == BUNEMAN) {
		exchangeBunemanBfield(bunemanBx, bunemanBy, bunemanBz);
		exchangeBunemanEfield(bunemanEx, bunemanEy, bunemanEz);
	} else {
		exchangeEfield();
		exchangeGeneralBfield(Bfield);
		exchangeGeneralBfield(newBfield);
	}

	if ((rank == 0) && (verbosity > 1)) printf("initial collecting particles\n");
	updateParticleCorrelationMaps();
	updateParameters();
	//updateAnisotropy();

	//double because dielectric tensor needs deltaT;
	updateDeltaT();
	//updateElectroMagneticParameters();
	//updateDensityParameters();
	//updateDeltaT();
	if(solverType != BUNEMAN){
		evaluateParticlesRotationTensor();
		updateElectroMagneticParameters();
	}
	updateDensityParameters();
	if (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
		//updateShockWaveX();
	}

	//evaluateExplicitDerivative();
	if (currentIteration % divergenceCleanUpParameter == 0) {
		if (solverType == BUNEMAN) {
			//cleanupDivergenceBuneman();
			//cleanUpDivergenceBunemanMagnetic();
		} else {
			//cleanupDivergence(newEfield, chargeDensity);
			//cleanupDivergenceMagnetic();
		}
	}
	updateFields();
	//exchangeGeneralBfield(Bfield);
	updateEnergy();

	for (int i = 0; i < typesNumber; ++i) {
		types[i].injectionLength = types[i].particesDeltaX - 0.0001 * deltaX;
	}

	while (time * plasma_period < maxTime && currentIteration < maxIteration) {
		if (currentIteration % writeMemoryParameter == 0) {
			outputMemory((outputDir + "memory_using.dat").c_str(), cartComm, cartCoord, cartDim);
		}
		if ((rank == 0) && (verbosity > 0)) {
			FILE* logFile = fopen((outputDir + "log.dat").c_str(), "a");
			fprintf(logFile, "start iteration number = %d time = %15.10g\n", currentIteration, time);
			fflush(logFile);
			fclose(logFile);
		}
		if ((rank == 0)) printf("start iteration number = %d time = %15.10g\n", currentIteration, time);
		if(currentIteration % writeGeneralParameter == 0) {
			updateEnergy();
		}
		if (currentIteration % writeParameter == 0) {
			output();
			currentWriteNumber++;
		}
		if (currentIteration % writeGeneralParameter == 0) {
			if (rank == 0) outputGeneral((outputDir + "general.dat").c_str(), this);
		}

		if (currentIteration % writeTrajectoryNumber == 0) {
			outputTrajectories();
		}

		MPI_Barrier(cartComm);
		double iterationTime = 0;
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			iterationTime = clock();
		}

		//updateDeltaT();
		/////////////////////////////////////////////////////
		double procTime = 0;
		if (solverType == BUNEMAN) {
			tristanEvaluateBhalfStep();
			exchangeBunemanBfield(bunemanBx, bunemanBy, bunemanBz);
			/*if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock();
			}
			for(int i = 0; i < particles.size(); ++i) {
				Particle* particle = particles[i];
				updateBunemanCorrelationMaps(particle);
			}
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock() - procTime;
				printf("updating correlation maps = %g sec\n", procTime / CLOCKS_PER_SEC);
			}*/
			moveParticles();
			tristanEvaluateBhalfStep();
			exchangeBunemanBfield(bunemanBx, bunemanBy, bunemanBz);
			tristanUpdateFlux();
			////////////////////////////////
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock();
			}
			for (int n = 0; n < smoothingCount; ++n) {
				smoothBunemanEfieldGeneral(bunemanJx, bunemanJy, bunemanJz);
			}
			MPI_Barrier(cartComm);
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock() - procTime;
				printf("smoothing flux = %g sec\n", procTime / CLOCKS_PER_SEC);
			}
			//////////////////////////
			tristanEvaluateE();
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock();
			}
			exchangeBunemanEfield(bunemanEx, bunemanEy, bunemanEz);
			MPI_Barrier(cartComm);
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock() - procTime;
				printf("exchange buneman E field = %g sec\n", procTime / CLOCKS_PER_SEC);
			}
		} else {
			evaluateParticlesRotationTensor();

			updateElectroMagneticParameters();
			procTime = 0;
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock();
			}
			for (int n = 0; n < smoothingCount; ++n) {
				smoothChargeDensityHat();
				smoothFlux();
				smoothMatrixNodeParameter(dielectricTensor);
			}
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock() - procTime;
				printf("smoothing chatge density and flux = %g sec\n", procTime / CLOCKS_PER_SEC);
			}

			evaluateElectricField();
			exchangeEfield();

			/*cleanupDivergence(tempEfield, chargeDensityHat);
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						newEfield[i][j][k] = (tempEfield[i][j][k] - Efield[i][j][k] * (1 - theta)) / theta;
					}
				}
			}
			exchangeEfield();*/


			evaluateMagneticField();

			procTime = 0;
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock();
			}
			exchangeGeneralBfield(newBfield);
			exchangeGeneralEfield(tempEfield);
			exchangeGeneralEfield(newEfield);
			//MPI_Barrier(cartComm);
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock() - procTime;
				printf("exchanging all B field time = %g sec\n", procTime / CLOCKS_PER_SEC);
			}

			for (int n = 0; n < smoothingCount; ++n) {
				//smoothTempEfield();
			}

			//cleanupDivergence(tempEfield, chargeDensityHat);
			//cleanupDivergenceMagnetic();
			moveParticles();
		}


		if ((rank == 0) && (verbosity > 1)) printLog("erasing escaped particles\n");
		if ((rank == 0) && (verbosity > 1)) printf("erasing escaped particles\n");

		eraseEscapedPaticles();

		//MPI_Barrier(cartComm);
		if ((rank == 0) && (verbosity > 1)) printLog("start exchange particles\n");
		if ((rank == 0) && (verbosity > 1)) printf("start exchange particles\n");
		//MPI_Barrier(cartComm);
		updateTheoreticalEnergy();
		exchangeParticles();

		if ((rank == 0) && (verbosity > 1)) printLog("start injecting new particles\n");
		if ((rank == 0) && (verbosity > 1)) printf("start injecting new particles\n");


		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock();
		}
		for (int i = 0; i < typesNumber; ++i) {
			types[i].injectionLength += fabs(V0.x * deltaT);
		}
		if (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT || boundaryConditionTypeX == FREE_BOTH) {
			injectNewParticles();
		}
		if (preserveChargeGlobal && (boundaryConditionTypeX != PERIODIC)) {
			addToPreserveChargeGlobal();
		}
		//MPI_Barrier(cartComm);
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("injecting and preserving time = %g sec\n", procTime / CLOCKS_PER_SEC);
		}

		if(solverType != BUNEMAN || ((currentIteration + 1) % writeParameter == 0)){
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock();
			}
			updateParticleCorrelationMaps();
			//MPI_Barrier(cartComm);
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock() - procTime;
				printf("updating correlation maps = %g sec\n", procTime / CLOCKS_PER_SEC);
			}


			updateDensityParameters();
		}

		if (solverType != BUNEMAN) {

			/*for(int n = 0; n < smoothingCount; ++n){
				smoothChargeDensity();
			}*/

			if ((rank == 0) && (verbosity > 0)) {
				printf("finish update density parameters\n");
			}
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock();
			}
			exchangeGeneralEfield(newEfield);
			if ((rank == 0) && (verbosity > 0)) {
				printf("finish exchange new efield\n");
			}
			//MPI_Barrier(cartComm);
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock() - procTime;
				printf("exchanging fields time = %g sec\n", procTime / CLOCKS_PER_SEC);
			}

			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock();
			}

			/*for(int n = 0; n < smoothingCount; ++n){
				smoothNewEfield();
				smoothNewBfield();
			}*/

			//MPI_Barrier(cartComm);
			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock() - procTime;
				printf("smoothing fields time = %g sec\n", procTime / CLOCKS_PER_SEC);
			}
			if (currentIteration % divergenceCleanUpParameter == 0) {
				//cleanupDivergence(newEfield, chargeDensity);
				//cleanupDivergenceMagnetic();

			}

			/*
			for(int n = 0; n < smoothingCount; ++n){
				smoothNewEfield();
				smoothNewBfield();
			}
			*/

			updateFields();
			if ((rank == 0) && (verbosity > 0)) {
				printf("finish update fields\n");
			}

			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock();
			}
			exchangeGeneralEfield(newEfield);
			exchangeGeneralEfield(Efield);
			exchangeGeneralBfield(Bfield);

			if ((rank == 0) && (verbosity > 0)) {
				printf("finish exchange fields\n");
			}

			if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
				procTime = clock() - procTime;
				printf("exchanging fields time = %g sec\n", procTime / CLOCKS_PER_SEC);
			}
		}


		if ((currentIteration + 1) % writeGeneralParameter == 0) {
			updateParameters();
			//updateAnisotropy();
		}

		removeEscapedParticles();
		if ((rank == 0) && (verbosity > 0)) {
			printf("finish remove escaped particles\n");
		}

		/*if(currentIteration % splitParticlesParameter == 0){
		    if ((rank == 0) && (verbosity > 0)) printf("split particles\n");
		    splitParticles();
		}*/

		time += deltaT;
		MPI_Barrier(cartComm);
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			iterationTime = clock() - iterationTime;
			printf("iteration except outputing time = %g sec\n", iterationTime / CLOCKS_PER_SEC);
		}

		currentIteration++;
		if (currentIteration % writeBackupParameter == 0) {
			outputBackup();
		}
	}
}

void Simulation::outputTrajectories() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 1)) printf("outputing trajectories\n");
	outputParticlesTrajectories((outputDir + "particlesTrajectories.dat").c_str(),
	                            (outputDir + "electronsTrajectories.dat").c_str(), particles, trackedParticlesNumbers,
	                            trackedParticlesNumber, time, plasma_period, scaleFactor, this);

	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("outputing trajectories time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::output() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0)) printf("outputing iteration number %d\n", currentIteration);
	if ((rank == 0) && (verbosity > 0)) printLog("collecting most accelerate particles\n");
	collectMostAcceleratedParticles();
	if ((rank == 0) && (verbosity > 0)) printLog("outputing\n");

	std::string fileNumber = "";
	if (multiplyFileOutput) {
		fileNumber = std::string("_") + convertIntToString(currentWriteNumber);
	}

	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution protons\n");
	outputDistribution((outputDir + "distribution_protons" + fileNumber + ".dat").c_str(), particles, PROTON, scaleFactor,
	                   plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution electrons\n");
	outputDistribution((outputDir + "distribution_electrons" + fileNumber + ".dat").c_str(), particles, ELECTRON, scaleFactor,
	                   plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution alphas\n");
	outputDistribution((outputDir + "distribution_alphas" + fileNumber + ".dat").c_str(), particles, ALPHA, scaleFactor,
	                   plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution positrons\n");
	outputDistribution((outputDir + "distribution_positrons" + fileNumber + ".dat").c_str(), particles, POSITRON, scaleFactor,
	                   plasma_period, verbosity, multiplyFileOutput);

	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution protons\n");
	outputDistributionByXgrid((outputDir + "distribution_protons_grid" + fileNumber + ".dat").c_str(), this, particles, PROTON, scaleFactor,
	                   plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution electrons\n");
	outputDistributionByXgrid((outputDir + "distribution_electrons_grid" + fileNumber + ".dat").c_str(), this, particles, ELECTRON, scaleFactor,
	                   plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution alphas\n");
	outputDistributionByXgrid((outputDir + "distribution_alphas_grid" + fileNumber + ".dat").c_str(), this, particles, ALPHA, scaleFactor,
	                   plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution positrons\n");
	outputDistributionByXgrid((outputDir + "distribution_positrons_grid" + fileNumber + ".dat").c_str(), this, particles, POSITRON, scaleFactor,
	                   plasma_period, verbosity, multiplyFileOutput);

	Vector3d shockWaveV = V0 / 3;

	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution protons shock wave\n");
	/*outputDistributionShiftedSystem((outputDir + "distribution_protons_sw" + fileNumber + ".dat").c_str(), particles, shockWaveV,
	                                speed_of_light_normalized, PROTON, scaleFactor,
	                                plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution electrons shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_electrons_sw" + fileNumber + ".dat").c_str(), particles, shockWaveV,
	                                speed_of_light_normalized, ELECTRON, scaleFactor,
	                                plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution alphas shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_alphas_sw" + fileNumber + ".dat").c_str(), particles, shockWaveV,
	                                speed_of_light_normalized, ALPHA, scaleFactor,
	                                plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution positrons shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_positrons_sw" + fileNumber + ".dat").c_str(), particles, shockWaveV,
	                                speed_of_light_normalized, POSITRON, scaleFactor,
	                                plasma_period, verbosity, multiplyFileOutput);*/
	/*outputDistributionShiftedSystem((outputDir + "distribution_protons_sw" + fileNumber + ".dat").c_str(), particles, shockWaveV,
	                                1.0, PROTON, scaleFactor,
	                                plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution electrons shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_electrons_sw" + fileNumber + ".dat").c_str(), particles, shockWaveV,
	                                1.0, ELECTRON, scaleFactor,
	                                plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution alphas shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_alphas_sw" + fileNumber + ".dat").c_str(), particles, shockWaveV,
	                                1.0, ALPHA, scaleFactor,
	                                plasma_period, verbosity, multiplyFileOutput);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution positrons shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_positrons_sw" + fileNumber + ".dat").c_str(), particles, shockWaveV,
	                                1.0, POSITRON, scaleFactor,
	                                plasma_period, verbosity, multiplyFileOutput);*/

	int coordX = getCartCoordWithAbsoluteIndexX(xnumberGeneral / 2);
	int coordY = getCartCoordWithAbsoluteIndexY(ynumberGeneral / 2);
	int coordZ = getCartCoordWithAbsoluteIndexZ(znumberGeneral / 2);

	if ((rank == 0) && (verbosity > 1)) printf("outputing fields\n");

	if (solverType == BUNEMAN) {
		double fieldScale = 1.0 / (plasma_period * sqrt(scaleFactor));

		resetBunemanFieldToCellVectorParameter(bunemanEx, bunemanEy, bunemanEz);

		//outputVectorCellArray((outputDir + "Efield.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim, multiplyFileOutput, fieldScale);

		if (verbosity > 2) printf("get cart coord with absolute index rank = %d\n", rank);

		MPI_Barrier(cartComm);
		if (verbosity > 2) printf("output crossection fields yz\n");
		if (coordX == cartCoord[0]) {
			int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
			if (verbosity > 2) printf("x local index = %d\n", xindex);
			outputVectorCellArrayCrossectionYZ((outputDir + "EfieldYZ.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                                   additionalBinNumber, cartCommYZ, cartCommZ, cartCoord, cartDim, xindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);

		if (verbosity > 2) printf("output crossection fields xz\n");
		if (coordY == cartCoord[1]) {
			int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
			if (verbosity > 2) printf("y local index = %d\n", yindex);
			outputVectorCellArrayCrossectionXZ((outputDir + "EfieldXZ.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                                   additionalBinNumber, cartCommXZ, cartCommZ, cartCoord, cartDim, yindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);

		if (verbosity > 2) printf("output crossection fields xy\n");
		if (coordZ == cartCoord[2]) {
			int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
			if (verbosity > 2) printf("z local index = %d\n", zindex);
			outputVectorCellArrayCrossectionXY((outputDir + "EfieldXY.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                                   additionalBinNumber, cartCommXY, cartCommY, cartCoord, cartDim, zindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);

		if (verbosity > 2) printf("output line field x\n");
		if (coordY == cartCoord[1] && coordZ == cartCoord[2]) {
			int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
			int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
			if (verbosity > 2) printf("y local index = %d\n", yindex);
			if (verbosity > 2) printf("z local index = %d\n", zindex);
			outputVectorCellArrayLineX((outputDir + "EfieldX.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                           additionalBinNumber, cartCommX, cartCoord, cartDim, yindex, zindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);

		if (verbosity > 2) printf("output line field y\n");
		if (coordX == cartCoord[0] && coordZ == cartCoord[2]) {
			int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
			int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
			if (verbosity > 2) printf("x local index = %d\n", xindex);
			if (verbosity > 2) printf("z local index = %d\n", zindex);
			outputVectorCellArrayLineY((outputDir + "EfieldY.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                           additionalBinNumber, cartCommY, cartCoord, cartDim, xindex, zindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);

		if (verbosity > 2) printf("output line field z\n");
		if (coordX == cartCoord[0] && coordY == cartCoord[1]) {
			int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
			int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
			if (verbosity > 2) printf("y local index = %d\n", yindex);
			if (verbosity > 2) printf("x local index = %d\n", xindex);
			outputVectorCellArrayLineZ((outputDir + "EfieldZ.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                           additionalBinNumber, cartCommZ, cartCoord, cartDim, xindex, yindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);

		resetBunemanFieldToCellVectorParameter(bunemanBx, bunemanBy, bunemanBz);

		//outputVectorCellArray((outputDir + "Bfield.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim, multiplyFileOutput, fieldScale);

		MPI_Barrier(cartComm);
		if (verbosity > 2) printf("output crossection fields yz\n");
		if (coordX == cartCoord[0]) {
			int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
			if (verbosity > 2) printf("x local index = %d\n", xindex);
			outputVectorCellArrayCrossectionYZ((outputDir + "BfieldYZ.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                                   additionalBinNumber, cartCommYZ, cartCommZ, cartCoord, cartDim, xindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);

		if (verbosity > 2) printf("output crossection fields xz\n");
		if (coordY == cartCoord[1]) {
			int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
			if (verbosity > 2) printf("y local index = %d\n", yindex);
			outputVectorCellArrayCrossectionXZ((outputDir + "BfieldXZ.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                                   additionalBinNumber, cartCommXZ, cartCommZ, cartCoord, cartDim, yindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);

		if (verbosity > 2) printf("output crossection fields xy\n");
		if (coordZ == cartCoord[2]) {
			int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
			if (verbosity > 2) printf("z local index = %d\n", zindex);
			outputVectorCellArrayCrossectionXY((outputDir + "BfieldXY.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                                   additionalBinNumber, cartCommXY, cartCommY, cartCoord, cartDim, zindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);

		if (verbosity > 2) printf("output line field x\n");
		if (coordY == cartCoord[1] && coordZ == cartCoord[2]) {
			int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
			int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
			if (verbosity > 2) printf("y local index = %d\n", yindex);
			if (verbosity > 2) printf("z local index = %d\n", zindex);
			outputVectorCellArrayLineX((outputDir + "BfieldX.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                           additionalBinNumber, cartCommX, cartCoord, cartDim, yindex, zindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);

		if (verbosity > 2) printf("output line field y\n");
		if (coordX == cartCoord[0] && coordZ == cartCoord[2]) {
			int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
			int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
			if (verbosity > 2) printf("x local index = %d\n", xindex);
			if (verbosity > 2) printf("z local index = %d\n", zindex);
			outputVectorCellArrayLineY((outputDir + "BfieldY.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                           additionalBinNumber, cartCommY, cartCoord, cartDim, xindex, zindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);

		if (verbosity > 2) printf("output line field z\n");
		if (coordX == cartCoord[0] && coordY == cartCoord[1]) {
			int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
			int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
			if (verbosity > 2) printf("y local index = %d\n", yindex);
			if (verbosity > 2) printf("x local index = %d\n", xindex);
			outputVectorCellArrayLineZ((outputDir + "BfieldZ.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded,
			                           additionalBinNumber, cartCommZ, cartCoord, cartDim, xindex, yindex, multiplyFileOutput, fieldScale);
		}
		MPI_Barrier(cartComm);
	} else {
		if (verbosity > 2) printf("get cart coord with absolute index rank = %d\n", rank);
		if (verbosity > 2) printf("x coord with absolute index = %d\n", coordX);
		MPI_Barrier(cartComm);
		if (verbosity > 2) printf("output crossection fields x\n");
		if (coordX == cartCoord[0]) {
			int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
			if (verbosity > 2) printf("x local index = %d\n", xindex);
			outputFieldsCrossectionYZ((outputDir + "EfieldYZ.dat").c_str(), (outputDir + "BfieldYZ" + fileNumber + ".dat").c_str(), Efield, Bfield,
			                          xnumberAdded,
			                          ynumberAdded, znumberAdded, additionalBinNumber, plasma_period, scaleFactor, cartCommYZ,
			                          cartCommZ, cartCoord, cartDim, xindex, multiplyFileOutput);
		}
		MPI_Barrier(cartComm);
		if (verbosity > 2) printf("output crossection fields y\n");
		if (coordY == cartCoord[1]) {
			int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
			if (verbosity > 2) printf("y local index = %d\n", yindex);
			outputFieldsCrossectionXZ((outputDir + "EfieldXZ.dat").c_str(), (outputDir + "BfieldXZ" + fileNumber + ".dat").c_str(), Efield, Bfield,
			                          xnumberAdded,
			                          ynumberAdded, znumberAdded, additionalBinNumber, plasma_period, scaleFactor, cartCommXZ,
			                          cartCommZ, cartCoord, cartDim, yindex, multiplyFileOutput);
		}
		MPI_Barrier(cartComm);
		if (verbosity > 2) printf("output crossection fields z\n");
		if (coordZ == cartCoord[2]) {
			int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
			if (verbosity > 2) printf("z local index = %d\n", zindex);
			outputFieldsCrossectionXY((outputDir + "EfieldXY.dat").c_str(), (outputDir + "BfieldXY" + fileNumber + ".dat").c_str(), Efield, Bfield,
			                          xnumberAdded,
			                          ynumberAdded, znumberAdded, additionalBinNumber, plasma_period, scaleFactor, cartCommXY,
			                          cartCommY, cartCoord, cartDim, zindex, multiplyFileOutput);
		}
		MPI_Barrier(cartComm);
		if (verbosity > 2) printf("output line fields x\n");
		if (coordY == cartCoord[1] && coordZ == cartCoord[2]) {
			int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
			int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
			if (verbosity > 2) printf("y local index = %d\n", yindex);
			if (verbosity > 2) printf("z local index = %d\n", zindex);
			outputFieldsLineX((outputDir + "EfieldX.dat").c_str(), (outputDir + "BfieldX" + fileNumber + ".dat").c_str(), Efield, Bfield,
			                  xnumberAdded,
			                  ynumberAdded, znumberAdded, additionalBinNumber, plasma_period, scaleFactor, cartCommX, cartCoord,
			                  cartDim, yindex, zindex, multiplyFileOutput);
		}
		MPI_Barrier(cartComm);
		if (verbosity > 2) printf("output line fields y\n");
		if (coordX == cartCoord[0] && coordZ == cartCoord[2]) {
			int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
			int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
			if (verbosity > 2) printf("x local index = %d\n", xindex);
			if (verbosity > 2) printf("z local index = %d\n", zindex);
			outputFieldsLineY((outputDir + "EfieldY.dat").c_str(), (outputDir + "BfieldY" + fileNumber + ".dat").c_str(), Efield, Bfield,
			                  xnumberAdded,
			                  ynumberAdded, znumberAdded, additionalBinNumber, plasma_period, scaleFactor, cartCommY, cartCoord,
			                  cartDim, xindex, zindex, multiplyFileOutput);
		}
		MPI_Barrier(cartComm);
		if (verbosity > 2) printf("output line fields z\n");
		if (coordY == cartCoord[1] && coordX == cartCoord[0]) {
			int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
			int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
			if (verbosity > 2) printf("y local index = %d\n", yindex);
			if (verbosity > 2) printf("x local index = %d\n", xindex);
			outputFieldsLineZ((outputDir + "EfieldZ.dat").c_str(), (outputDir + "BfieldZ" + fileNumber + ".dat").c_str(), Efield, Bfield,
			                  xnumberAdded,
			                  ynumberAdded, znumberAdded, additionalBinNumber, plasma_period, scaleFactor, cartCommZ, cartCoord,
			                  cartDim, xindex, yindex, multiplyFileOutput);
		}
	}

	if ((rank == 0) && (verbosity > 1)) printf("outputing grid\n");
	outputGridX((outputDir + "Xfile.dat").c_str(), xgrid, xnumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim,
	            true, scaleFactor);
	/*outputGridReducedX((outputDir + "XfileReduced.dat").c_str(), xgrid, xnumberAdded, additionalBinNumber, reduceStepX,
	                   rank, leftRank, rightRank, cartComm, cartCoord, cartDim, true, scaleFactor);*/

	outputGridY((outputDir + "Yfile.dat").c_str(), ygrid, ynumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim,
	            true, scaleFactor);
	/*outputGridReducedY((outputDir + "YfileReduced.dat").c_str(), ygrid, ynumberAdded, additionalBinNumber, reduceStepY,
	                   rank, frontRank, backRank, cartComm, cartCoord, cartDim, true, scaleFactor);*/

	outputGridZ((outputDir + "Zfile.dat").c_str(), zgrid, znumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim,
	            true, scaleFactor);
	/*outputGridReducedZ((outputDir + "ZfileReduced.dat").c_str(), zgrid, znumberAdded, additionalBinNumber, reduceStepZ,
	                   rank, bottomRank, topRank, cartComm, cartCoord, cartDim, true, scaleFactor);*/


	//if ((rank == 0) && (verbosity > 1)) printf("outputing concentrations\n");
	/*outputConcentrations((outputDir + "concentrations.dat").c_str(), particleConcentrations, chargeDensity,
	                     chargeDensityHat, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, typesNumber,
	                     plasma_period, scaleFactor, cartComm, cartCoord, cartDim);*/

	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output crossection concentrations x\n");
	if (coordX == cartCoord[0]) {
		int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
		if (verbosity > 2) printf("x local index = %d\n", xindex);
		outputConcentrationsCrossectionYZ((outputDir + "concentrationsYZ" + fileNumber + ".dat").c_str(), particleConcentrations, chargeDensity,
		                                  chargeDensityHat, xnumberAdded,
		                                  ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period,
		                                  scaleFactor, cartCommYZ, cartCommZ, cartCoord, cartDim, xindex, multiplyFileOutput);
	}
	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output crossection concentrations y\n");
	if (coordY == cartCoord[1]) {
		int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
		if (verbosity > 2) printf("y local index = %d\n", yindex);
		outputConcentrationsCrossectionXZ((outputDir + "concentrationsXZ" + fileNumber + ".dat").c_str(), particleConcentrations, chargeDensity,
		                                  chargeDensityHat, xnumberAdded,
		                                  ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period,
		                                  scaleFactor, cartCommXZ, cartCommZ, cartCoord, cartDim, yindex, multiplyFileOutput);
	}
	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output crossection concentrations z\n");
	if (coordZ == cartCoord[2]) {
		int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
		if (verbosity > 2) printf("z local index = %d\n", zindex);
		outputConcentrationsCrossectionXY((outputDir + "concentrationsXY" + fileNumber + ".dat").c_str(), particleConcentrations, chargeDensity,
		                                  chargeDensityHat, xnumberAdded,
		                                  ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period,
		                                  scaleFactor, cartCommXY, cartCommY, cartCoord, cartDim, zindex, multiplyFileOutput);
	}
	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output line concentrations x\n");
	if (coordY == cartCoord[1] && coordZ == cartCoord[2]) {
		int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
		int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
		if (verbosity > 2) printf("y local index = %d\n", yindex);
		if (verbosity > 2) printf("z local index = %d\n", zindex);
		outputConcentrationsLineX((outputDir + "concentrationsX" + fileNumber + ".dat").c_str(), particleConcentrations, chargeDensity,
		                          chargeDensityHat, xnumberAdded,
		                          ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period, scaleFactor,
		                          cartCommX, cartCoord, cartDim, yindex, zindex, multiplyFileOutput);
	}
	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output line concentrations y\n");
	if (coordX == cartCoord[0] && coordZ == cartCoord[2]) {
		int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
		int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
		if (verbosity > 2) printf("x local index = %d\n", xindex);
		if (verbosity > 2) printf("z local index = %d\n", zindex);
		outputConcentrationsLineY((outputDir + "concentrationsY" + fileNumber + ".dat").c_str(), particleConcentrations, chargeDensity,
		                          chargeDensityHat, xnumberAdded,
		                          ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period, scaleFactor,
		                          cartCommY, cartCoord, cartDim, xindex, zindex, multiplyFileOutput);
	}
	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output line concentrations z\n");
	if (coordY == cartCoord[1] && coordX == cartCoord[0]) {
		int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
		int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
		if (verbosity > 2) printf("y local index = %d\n", yindex);
		if (verbosity > 2) printf("x local index = %d\n", xindex);
		outputConcentrationsLineZ((outputDir + "concentrationsZ" + fileNumber + ".dat").c_str(), particleConcentrations, chargeDensity,
		                          chargeDensityHat, xnumberAdded,
		                          ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period, scaleFactor,
		                          cartCommZ, cartCoord, cartDim, xindex, yindex, multiplyFileOutput);
	}

	if ((rank == 0) && (verbosity > 1)) printf("outputing velocity\n");
	/*outputVelocity((outputDir + "velocity.dat").c_str(),
	               particleBulkVelocities, types, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber,
	               typesNumber, plasma_period, scaleFactor, cartComm, cartCoord, cartDim);*/

	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output crossection velocity x\n");
	if (coordX == cartCoord[0]) {
		int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
		if (verbosity > 2) printf("x local index = %d\n", xindex);
		outputVelocityCrossectionYZ((outputDir + "velocityYZ" + fileNumber + ".dat").c_str(), particleBulkVelocities, types, xnumberAdded,
		                            ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period, scaleFactor,
		                            cartCommYZ, cartCommZ, cartCoord, cartDim, xindex, multiplyFileOutput);
	}
	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output crossection velocity y\n");
	if (coordY == cartCoord[1]) {
		int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
		if (verbosity > 2) printf("y local index = %d\n", yindex);
		outputVelocityCrossectionXZ((outputDir + "velocityXZ" + fileNumber + ".dat").c_str(), particleBulkVelocities, types, xnumberAdded,
		                            ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period, scaleFactor,
		                            cartCommXZ, cartCommZ, cartCoord, cartDim, yindex, multiplyFileOutput);
	}
	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output crossection velocity z\n");
	if (coordZ == cartCoord[2]) {
		int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
		if (verbosity > 2) printf("z local index = %d\n", zindex);
		outputVelocityCrossectionXY((outputDir + "velocityXY" + fileNumber + ".dat").c_str(), particleBulkVelocities, types, xnumberAdded,
		                            ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period, scaleFactor,
		                            cartCommXY, cartCommY, cartCoord, cartDim, zindex, multiplyFileOutput);
	}
	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output line velocity x\n");
	if (coordY == cartCoord[1] && coordZ == cartCoord[2]) {
		int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
		int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
		if (verbosity > 2) printf("y local index = %d\n", yindex);
		if (verbosity > 2) printf("z local index = %d\n", zindex);
		outputVelocityLineX((outputDir + "velocityX" + fileNumber + ".dat").c_str(), particleBulkVelocities, types, xnumberAdded,
		                    ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period, scaleFactor,
		                    cartCommX, cartCoord, cartDim, yindex, zindex, multiplyFileOutput);
	}
	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output line velocity y\n");
	if (coordX == cartCoord[0] && coordZ == cartCoord[2]) {
		int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
		int zindex = getLocalIndexByAbsoluteZ(znumberGeneral / 2);
		if (verbosity > 2) printf("x local index = %d\n", xindex);
		if (verbosity > 2) printf("z local index = %d\n", zindex);
		outputVelocityLineY((outputDir + "velocityY" + fileNumber + ".dat").c_str(), particleBulkVelocities, types, xnumberAdded,
		                    ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period, scaleFactor,
		                    cartCommY, cartCoord, cartDim, xindex, zindex, multiplyFileOutput);
	}
	MPI_Barrier(cartComm);
	if (verbosity > 2) printf("output line velocity z\n");
	if (coordY == cartCoord[1] && coordX == cartCoord[0]) {
		int yindex = getLocalIndexByAbsoluteY(ynumberGeneral / 2);
		int xindex = getLocalIndexByAbsoluteX(xnumberGeneral / 2);
		if (verbosity > 2) printf("y local index = %d\n", yindex);
		if (verbosity > 2) printf("x local index = %d\n", xindex);
		outputVelocityLineZ((outputDir + "velocityZ" + fileNumber + ".dat").c_str(), particleBulkVelocities, types, xnumberAdded,
		                    ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period, scaleFactor,
		                    cartCommZ, cartCoord, cartDim, xindex, yindex, multiplyFileOutput);
	}

	if ((rank == 0) && (verbosity > 1)) printf("outputing flux\n");
	if(solverType == BUNEMAN){
		double fluxScale = 1.0/(plasma_period*plasma_period*sqrt(scaleFactor));
		resetBunemanFieldToCellVectorParameter(bunemanJx, bunemanJy, bunemanJz);
		outputVectorCellArray((outputDir + "flux.dat").c_str(), tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim, multiplyFileOutput, fluxScale);
	} else {
		outputFlux((outputDir + "flux.dat").c_str(), electricFlux, externalElectricFlux, xnumberAdded, ynumberAdded,
	           znumberAdded, additionalBinNumber, plasma_period, scaleFactor, cartComm, cartCoord, cartDim, multiplyFileOutput);
	}

	if ((rank == 0) && (verbosity > 1)) printf("outputing divergence\n");
	outputDivergenceError((outputDir + "divergence_error.dat").c_str(), this, plasma_period, scaleFactor, multiplyFileOutput);

	double rotBscale = 1.0 / (plasma_period * plasma_period * sqrt(scaleFactor));

	if ((rank == 0) && (verbosity > 1)) printf("outputing rotB\n");
	outputVectorNodeArray((outputDir + "rotBFile.dat").c_str(), rotB, xnumberAdded, ynumberAdded, znumberAdded,
	                      additionalBinNumber, cartComm, cartCoord, cartDim, rotBscale);

	if ((rank == 0) && (verbosity > 1)) printf("outputing Ederivative\n");
	outputVectorNodeArray((outputDir + "EderivativeFile.dat").c_str(), Ederivative, xnumberAdded, ynumberAdded,
	                      znumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim, rotBscale);

	if ((rank == 0) && (verbosity > 1)) printf("outputing rotE\n");
	outputVectorCellArray((outputDir + "rotEFile.dat").c_str(), rotE, xnumberAdded, ynumberAdded, znumberAdded,
	                      additionalBinNumber, cartComm, cartCoord, cartDim, rotBscale);

	/*if (rank == 0) printf("outputing dielectricTensor\n");
	outputMatrixArray((outputDir + "dielectricTensorFile.dat").c_str(), dielectricTensor, xnumber + 1, ynumber + 1,
	                  znumber + 1);*/

	if ((rank == 0) && (verbosity > 1)) printf("outputing particles\n");
	for (int i = 0; i < typesNumber; ++i) {
		outputParticles((outputDir + types[i].typeName + ".dat").c_str(), this, types[i].type);
	}

	if ((rank == 0) && (verbosity > 1)) printf("outputing accelerated particles\n");
	outputAcceleratedParticlesNumbers((outputDir + "acceleratedParticlesNumbers.dat").c_str(), this);

	/*maxwellMatrixFile = fopen("./output/maxwellMatrixFile.dat", "w");
	outputMaxwellEquationMatrixFull(maxwellMatrixFile, maxwellEquationMatrix, xnumber, ynumber, znumber, maxwellEquationMatrixSize);
	fclose(maxwellMatrixFile);*/

	//if (rank == 0) outputGeneral((outputDir + "general.dat").c_str(), this);
	if (rank == 0) outputGeneralAnisotropy((outputDir + "generalAnisotropy.dat").c_str(), this);
	if ((rank == 0) && (verbosity > 0)) printf("finish outputing\n");
	if ((rank == 0) && (verbosity > 0)) printLog("finish outputing\n");
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("outputing time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::resetBunemanFieldToCellVectorParameter(double*** bunemanEx, double*** bunemanEy, double*** bunemanEz) {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				tempCellVectorParameter[i][j][k].x = bunemanEx[i][j][k];
				tempCellVectorParameter[i][j][k].y = bunemanEy[i][j][k];
				tempCellVectorParameter[i][j][k].z = bunemanEz[i][j][k];
			}
		}
	}
}

void Simulation::outputBackup() {
	if ((rank == 0) && (verbosity > 0)) printf("writing backup\n");
	if ((rank == 0) && (verbosity > 0)) printLog("writing backup\n");

	std::string backupDir = backupDirectory;
	//FILE* backupGeneralFile = fopen((backupDir + "general.dat").c_str(), "r");
	//FILE* backupEfieldFile = fopen((backupDir + "Efield.dat").c_str(), "r");
	//FILE* backupBfieldFile = fopen((backupDir + "Bfield.dat").c_str(), "r");
	//FILE* backupParticlesFile = fopen((backupDir + "particles.dat").c_str(), "r");

	outputSimulationBackup((backupDir + "general.dat").c_str(), (backupDir + "Efield.dat").c_str(),
	                       (backupDir + "Bfield.dat").c_str(), (backupDir + "particles.dat").c_str(), this);

	//fclose(backupGeneralFile);
	//fclose(backupEfieldFile);
	//fclose(backupBfieldFile);
	//fclose(backupParticlesFile);
}

void Simulation::collectMostAcceleratedParticles() {
	for (int t = 0; t < typesNumber; ++t) {
		mostAcceleratedParticlesNumbers[t].clear();
	}
	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		int t = getTypeNumber(particle);
		double momentum = particle->getMomentum().norm();
		if (mostAcceleratedParticlesNumbers[t].size() == 0) {
			mostAcceleratedParticlesNumbers[t].push_back(std::pair < int, double >(particle->number, momentum));
		} else {
			std::list < std::pair < int, double > >::iterator it = mostAcceleratedParticlesNumbers[t].begin();
			std::pair < int, double > pair = *it;
			if (momentum > pair.second) {
				bool added = false;
				for (it = mostAcceleratedParticlesNumbers[t].begin(); it != mostAcceleratedParticlesNumbers[t].end(); ++it) {
					pair = *it;
					if (momentum < pair.second) {
						mostAcceleratedParticlesNumbers[t].insert(it, std::pair < int, double >(particle->number, momentum));
						added = true;
						break;
					}
				}
				if (!added) {
					mostAcceleratedParticlesNumbers[t].push_back(std::pair < int, double >(particle->number, momentum));
				}
				if (mostAcceleratedParticlesNumbers[t].size() > mostFastParticlesNumber) {
					mostAcceleratedParticlesNumbers[t].pop_front();
				}
			} else if (mostAcceleratedParticlesNumbers[t].size() < mostFastParticlesNumber) {
				mostAcceleratedParticlesNumbers[t].push_front(std::pair < int, double >(particle->number, momentum));
			}
		}
	}

	for (int t = 0; t < typesNumber; ++t) {
		if (types[t].particlesPerBin > 0) {
			while (mostAcceleratedParticlesNumbers[t].size() < mostFastParticlesNumber) {
				mostAcceleratedParticlesNumbers[t].push_front(std::pair < int, double >(-1, 0));
			}
		}
	}

	for (int t = 0; t < typesNumber; ++t) {
		if (types[t].particlesPerBin > 0) {

			int tempAcceleratedNumbers[mostFastParticlesNumber];
			double tempAcceleratedMomentum[mostFastParticlesNumber];

			std::list < std::pair < int, double > >::iterator it = mostAcceleratedParticlesNumbers[t].begin();
			for (int i = 0; i < mostFastParticlesNumber; ++i) {
				std::pair < int, double > pair = *it;
				tempAcceleratedNumbers[i] = pair.first;
				tempAcceleratedMomentum[i] = pair.second;
				++it;
			}


			int inAcceleratedNumbers[mostFastParticlesNumber];
			double inAcceleratedMomentum[mostFastParticlesNumber];

			if (rank > 0) {
				MPI_Send(tempAcceleratedNumbers, mostFastParticlesNumber, MPI_INT, 0, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm);
				MPI_Send(tempAcceleratedMomentum, mostFastParticlesNumber, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_ALL_TO_FIRST, cartComm);
			} else {
				for (int i = 1; i < nprocs; ++i) {
					MPI_Status status;
					MPI_Recv(inAcceleratedNumbers, mostFastParticlesNumber, MPI_INT, i, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm,
					         &status);
					MPI_Recv(inAcceleratedMomentum, mostFastParticlesNumber, MPI_DOUBLE, i, MPI_SEND_DOUBLE_ALL_TO_FIRST, cartComm,
					         &status);

					for (int j = 0; j < mostFastParticlesNumber; ++j) {
						int curNumber = inAcceleratedNumbers[j];
						int momentum = inAcceleratedMomentum[j];
						if (mostAcceleratedParticlesNumbers[t].size() == 0) {
							mostAcceleratedParticlesNumbers[t].push_back(std::pair < int, double >(curNumber, momentum));
						} else {
							it = mostAcceleratedParticlesNumbers[t].begin();
							std::pair < int, double > pair = *it;
							if (momentum > pair.second) {
								bool added = false;
								for (it = mostAcceleratedParticlesNumbers[t].begin(); it != mostAcceleratedParticlesNumbers[t].end(); ++it) {
									pair = *it;
									if (momentum < pair.second) {
										mostAcceleratedParticlesNumbers[t].insert(it, std::pair < int, double >(curNumber, momentum));
										added = true;
										break;
									}
								}
								if (!added) {
									mostAcceleratedParticlesNumbers[t].push_back(std::pair < int, double >(curNumber, momentum));
								}
								if (mostAcceleratedParticlesNumbers[t].size() > mostFastParticlesNumber) {
									mostAcceleratedParticlesNumbers[t].pop_front();
								}
							} else if (mostAcceleratedParticlesNumbers[t].size() < mostFastParticlesNumber) {
								mostAcceleratedParticlesNumbers[t].push_front(std::pair < int, double >(curNumber, momentum));
							}
						}
					}
				}
			}

			it = mostAcceleratedParticlesNumbers[t].begin();
			for (int i = 0; i < mostFastParticlesNumber; ++i) {
				std::pair < int, double > pair = *it;
				tempAcceleratedNumbers[i] = pair.first;
				tempAcceleratedMomentum[i] = pair.second;
				++it;
			}

			MPI_Bcast(tempAcceleratedNumbers, mostFastParticlesNumber, MPI_INT, 0, cartComm);
			MPI_Bcast(tempAcceleratedMomentum, mostFastParticlesNumber, MPI_DOUBLE, 0, cartComm);

			mostAcceleratedParticlesNumbers[t].clear();

			for (int i = 0; i < mostFastParticlesNumber; ++i) {
				mostAcceleratedParticlesNumbers[t].push_back(
					std::pair < int, double >(tempAcceleratedNumbers[i], tempAcceleratedMomentum[i]));
			}
		}
	}
}

void Simulation::updateDeltaT() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) printf("updating time step\n");
	if ((rank == 0) && (verbosity > 0)) printLog("updating time step\n");
	deltaT = preferedDeltaT / plasma_period;

	double delta = min2(deltaX, min2(deltaY, deltaZ));
	//double delta = deltaX;
	//deltaT = min2(deltaT, timeEpsilonKourant * delta / speed_of_light_normalized);
	deltaT = min2(deltaT, timeEpsilonKourant * delta);
	if ((rank == 0) && (verbosity > 1)) printLog("dx/c\n");
	if (particles.size() > 0) {
		double derEmax = 0;
		double B = B0.norm();
		double E = E0.norm();
		double nmax = types[0].concentration;
		if (solverType == BUNEMAN) {
			for (int i = 0; i < xnumberAdded; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						double Bnorm = sqrt(sqr(bunemanBx[i][j][k]) + sqr(bunemanBy[i][j][k]) + sqr(bunemanBz[i][j][k]));
						if (Bnorm > B) {
							B = Bnorm;
						}
					}
				}
			}
			for (int i = 0; i < xnumberAdded - 1; ++i) {
				for (int j = 0; j < ynumberAdded - 1; ++j) {
					for (int k = 0; k < znumberAdded - 1; ++k) {
						double Enorm = sqrt(sqr(bunemanEx[i][j][k]) + sqr(bunemanEy[i][j][k]) + sqr(bunemanEz[i][j][k]));
						if (Enorm > E) {
							E = Enorm;
						}
						if (fabs(bunemanEx[i + 1][j][k] - bunemanEx[i][j][k]) / deltaX > derEmax) {
							derEmax = fabs(bunemanEx[i + 1][j][k] - bunemanEx[i][j][k]) / deltaX;
						}
						if (fabs(bunemanEy[i][j + 1][k] - bunemanEy[i][j][k]) / deltaY > derEmax) {
							derEmax = fabs(bunemanEy[i][j + 1][k] - bunemanEy[i][j][k]) / deltaY;
						}
						if (fabs(bunemanEz[i][j][k + 1] - bunemanEz[i][j][k]) / deltaZ > derEmax) {
							derEmax = fabs(bunemanEz[i][j][k + 1] - bunemanEz[i][j][k]) / deltaZ;
						}
					}
				}
			}
		} else {
			for (int i = 0; i < xnumberAdded; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						if (Bfield[i][j][k].norm() > B) {
							B = Bfield[i][j][k].norm();
						}
					}
				}
			}
			for (int i = 0; i < xnumberAdded; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						if (Efield[i][j][k].norm() > E) {
							E = Efield[i][j][k].norm();
						}
						if ((Efield[i + 1][j][k] - Efield[i][j][k]).norm() / deltaX > derEmax) {
							derEmax = (Efield[i + 1][j][k] - Efield[i][j][k]).norm() / deltaX;
						}
						if ((Efield[i][j + 1][k] - Efield[i][j][k]).norm() / deltaY > derEmax) {
							derEmax = (Efield[i][j + 1][k] - Efield[i][j][k]).norm() / deltaY;
						}
						if ((Efield[i][j][k + 1] - Efield[i][j][k]).norm() / deltaZ > derEmax) {
							derEmax = (Efield[i][j + 1][k] - Efield[i][j][k]).norm() / deltaZ;
						}
					}
				}
			}
		}

		/*double minFlux = electricFlux[0][0][0].norm();
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					if (electricFlux[i][j][k].norm() < minFlux) {
						minFlux = electricFlux[i][j][k].norm();
					}
				}
			}
		}
		if ((rank == 0) && (verbosity > 1)) printLog("evaluatd E and minFlux\n");*/

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

		//double gamma = 1.0/sqrt(1.0 - V0.norm2()/speed_of_light_normalized_sqr);
		double gamma = 1.0/sqrt(1.0 - V0.norm2());

		omegaPlasmaProton = sqrt(
			4 * pi * types[1].concentration * electron_charge_normalized * electron_charge_normalized / (massProton*gamma));
		omegaPlasmaElectron = sqrt(
			4 * pi *types[0].concentration * electron_charge_normalized * electron_charge_normalized / (massElectron*gamma));
		omegaPlasmaTotal = 0;
		for(int i = 0; i < typesNumber; ++i) {
			omegaPlasmaTotal += 4 * pi *types[i].concentration * electron_charge_normalized * electron_charge_normalized / (types[i].mass*gamma);
		}
		omegaPlasmaTotal = sqrt(omegaPlasmaTotal);
		double omegaCapture = sqrt(electron_charge_normalized * derEmax / massElectron);


		//omegaGyroProton = electron_charge_normalized * B / (massProton * speed_of_light_normalized);
		omegaGyroProton = electron_charge_normalized * B / (massProton);
		//omegaGyroElectron = electron_charge_normalized * B / (massElectron * speed_of_light_normalized);
		omegaGyroElectron = electron_charge_normalized * B / (massElectron);


		if (timeEpsilonPlasma / omegaPlasmaTotal < deltaT) {
			if (rank == 0) {
				printf("plasma frequency time limitation\n");
			}
		}

		deltaT = min2(deltaT, timeEpsilonPlasma / omegaPlasmaTotal);

		double maxOmega = sqrt(4 * pi * nmax * types[0].charge * types[0].charge / types[0].mass);

		deltaT = min2(deltaT, timeEpsilonPlasma / maxOmega);


		double thermalMomentum = sqrt(massElectron * kBoltzman_normalized * temperature) + massElectron * V0.norm();


		//if ((rank == 0) && writeLog)) printLog("evaluzted omega\n");

		if (B > 0) {
			//if (timeEpsilonFields * gamma * massElectron * speed_of_light_normalized / (electron_charge_normalized * B) < deltaT) {
			if (timeEpsilonFields * gamma * massElectron / (electron_charge_normalized * B) < deltaT) {
				if (rank == 0) {
					printf("gyro frequency time limitation\n");
				}
			}
			//deltaT = min2(deltaT,timeEpsilonFields * gamma * massElectron * speed_of_light_normalized / (electron_charge_normalized * B));
			deltaT = min2(deltaT,timeEpsilonFields * gamma * massElectron / (electron_charge_normalized * B));
		}
		if (E > 0) {
			//if (timeEpsilonFields * massElectron * speed_of_light_normalized / (electron_charge_normalized * E) < deltaT) {
			if (timeEpsilonFields * massElectron / (electron_charge_normalized * E) < deltaT) {
				if (rank == 0) {
					printf("electric acceleration time limitation\n");
				}
			}
			//deltaT = min2(deltaT,timeEpsilonFields * massElectron * speed_of_light_normalized / (electron_charge_normalized * E));
			deltaT = min2(deltaT,timeEpsilonFields * massElectron / (electron_charge_normalized * E));
		}

		if (omegaCapture > 0) {
			/*if(timeEpsilonCapture / omegaCapture < deltaT){
				if(rank == 0){
					printf("capture frecuency acceleration time limitation\n");
				}
			}
			deltaT = min2(deltaT, timeEpsilonCapture / omegaCapture);*/
		}
		/*if (E > 0) {
		    deltaT = min2(deltaT, 0.1 * thermalMomentum / (electron_charge_normalized * E));
		}*/


		double Vthermal = sqrt(2 * kBoltzman_normalized * temperature / massElectron);
		double minDeltaT = deltaX / Vthermal;
		/*if(minDeltaT > deltaT){
		    printf("deltaT < dx/Vthermal\n");
		}*/
	}

	if ((rank == 0) && (verbosity > 1)) printLog("exchanging time step\n");
	if (nprocs > 1) {

		MPI_Barrier(cartComm);
		double temp[1];
		temp[0] = deltaT;
		if (rank != 0) {
			MPI_Send(temp, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, cartComm);
		} else {
			for (int i = 1; i < nprocs; ++i) {
				MPI_Status status;
				MPI_Recv(temp, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, cartComm, &status);
				if (temp[0] < deltaT) {
					deltaT = temp[0];
				}
			}
		}
		MPI_Barrier(cartComm);

		if (rank == 0) {
			temp[0] = deltaT;
			for (int i = 1; i < nprocs; ++i) {
				MPI_Send(temp, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, cartComm);
			}
		} else {
			MPI_Status status;
			MPI_Recv(temp, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, cartComm, &status);
			deltaT = temp[0];
		}
		MPI_Barrier(cartComm);
	}
	if ((rank == 0) && (verbosity > 0)) printLog("end updating time step\n");
	//MPI_Barrier(cartComm);
	updateParticlesBeta();
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating deltaT = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

double Simulation::volumeE() {
	return cellVolume;
}

double Simulation::volumeB() {
	return cellVolume;
}

void Simulation::checkParticleInBox(Particle& particle) {
	//alertNaNOrInfinity(particle.coordinates.x, "particle.x = NaN in check particle in box\n");
	if (particle.coordinates.x < xgrid[additionalBinNumber]) {
		printf("particle.x < 0\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.x = %15.10g < 0\n", particle.coordinates.x);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x,
		        particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		//fprintf(errorLogFile, "particle.v/c = %15.10g\n",(particle.getVelocity().norm() / speed_of_light_normalized));
		fprintf(errorLogFile, "particle.v/c = %15.10g\n",(particle.getVelocity().norm()));
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (particle.coordinates.x > xgrid[xnumberAdded - additionalBinNumber]) {
		printf("particle.x > xsize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.x = %15.10g > %15.10g\n", particle.coordinates.x,
		        xgrid[xnumberAdded - additionalBinNumber]);
		printf("particle.x = %15.10g > %15.10g\n", particle.coordinates.x, xgrid[xnumberAdded - additionalBinNumber]);
		printf("2*xsizeGeneral = %15.10g \n", 2 * xsizeGeneral);
		printf("xgrid[0] = %15.10g \n", xgrid[0]);
		printf("xsize = %15.10g \n", xsize);
		fprintf(errorLogFile, "2*xsizeGeneral = %15.10g \n", 2 * xsizeGeneral);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x,
		        particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		//fprintf(errorLogFile, "particle.v/c = %15.10g\n",(particle.getVelocity().norm() / speed_of_light_normalized));
		fprintf(errorLogFile, "particle.v/c = %15.10g\n",(particle.getVelocity().norm()));
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	if (particle.coordinates.y < ygrid[additionalBinNumber]) {
		if (rank == 0) printf("particle.y < ygrid[0]\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.y = %15.10g < %15.10g\n", particle.coordinates.y, ygrid[0]);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x,
		        particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		//fprintf(errorLogFile, "particle.v/c = %15.10g\n",(particle.getVelocity().norm() / speed_of_light_normalized));
		fprintf(errorLogFile, "particle.v/c = %15.10g\n",(particle.getVelocity().norm()));
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (particle.coordinates.y > ygrid[ynumberAdded - additionalBinNumber]) {
		printf("particle.y > ysize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.y = %15.10g > %15.10g\n", particle.coordinates.y,
		        ygrid[ynumberAdded - additionalBinNumber]);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x,
		        particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		//fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.getVelocity().norm() / speed_of_light_normalized));
		fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.getVelocity().norm()));
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	if (particle.coordinates.z < zgrid[additionalBinNumber]) {
		if (rank == 0) printf("particle.z < zgrid[0]\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.z = %15.10g < %15.10g\n", particle.coordinates.z, zgrid[0]);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x,
		        particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		//fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.getVelocity().norm() / speed_of_light_normalized));
		fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.getVelocity().norm()));
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (particle.coordinates.z > zgrid[znumberAdded - additionalBinNumber]) {
		if (rank == 0) printf("particle.z > zsize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.z = %15.10g > %15.10g\n", particle.coordinates.z,
		        zgrid[znumberAdded - additionalBinNumber]);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x,
		        particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		//fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.getVelocity().norm() / speed_of_light_normalized));
		fprintf(errorLogFile, "particle.v/c = %15.10g\n", (particle.getVelocity().norm()));
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
}

void Simulation::updateTheoreticalEnergy() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	int minI = 1 + additionalBinNumber;
	if (cartCoord[0] == 0 && boundaryConditionTypeX != PERIODIC) {
		minI = 0;
	}
	int maxI = xnumberAdded - 1 - additionalBinNumber;
	if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
		maxI = xnumberAdded - 1;
	}
	int minJ = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - 1 - additionalBinNumber;
	int minK = 1 + additionalBinNumber;
	int maxK = znumberAdded - 1 - additionalBinNumber;
	if(ynumberGeneral == 1) {
		minJ = 0;
		maxJ = 1;
	}
	if(znumberGeneral == 1) {
		minK = 0;
		maxK = 1;
	}

	//if (currentIteration > 0) {
	//if (boundaryConditionTypeX != PERIODIC) {

	for (int i = 0; i < escapedParticlesLeft.size(); ++i) {
		Particle* particle = escapedParticlesLeft[i];
		theoreticalEnergy -= particle->energy() * particle->weight *
			sqr(scaleFactor / plasma_period);
		theoreticalMomentum -= particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
	}
	for (int i = 0; i < escapedParticlesRight.size(); ++i) {
		Particle* particle = escapedParticlesRight[i];
		theoreticalEnergy -= particle->energy() * particle->weight *
			sqr(scaleFactor / plasma_period);
		theoreticalMomentum -= particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
	}

	for (int i = 0; i < escapedParticlesFront.size(); ++i) {
		Particle* particle = escapedParticlesFront[i];
		theoreticalEnergy -= particle->energy() * particle->weight *
			sqr(scaleFactor / plasma_period);
		theoreticalMomentum -= particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
	}
	for (int i = 0; i < escapedParticlesBack.size(); ++i) {
		Particle* particle = escapedParticlesBack[i];
		theoreticalEnergy -= particle->energy() * particle->weight *
			sqr(scaleFactor / plasma_period);
		theoreticalMomentum -= particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
	}

	for (int i = 0; i < escapedParticlesTop.size(); ++i) {
		Particle* particle = escapedParticlesTop[i];
		theoreticalEnergy -= particle->energy() * particle->weight *
			sqr(scaleFactor / plasma_period);
		theoreticalMomentum -= particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
	}
	for (int i = 0; i < escapedParticlesBottom.size(); ++i) {
		Particle* particle = escapedParticlesBottom[i];
		theoreticalEnergy -= particle->energy() * particle->weight *
			sqr(scaleFactor / plasma_period);
		theoreticalMomentum -= particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
	}

	for (int j = minJ; j < maxJ; ++j) {
		for (int k = minK; k < maxK; ++k) {
			Vector3d Eleft;
			Vector3d Eright;
			Vector3d Bleft;
			Vector3d Bright;
			if (solverType == BUNEMAN) {
				Eleft = getBunemanElectricField(1 + additionalBinNumber, j, k);
				Eright = getBunemanElectricField(xnumberAdded - additionalBinNumber - 1, j, k);
				Bleft = getBunemanMagneticField(1 + additionalBinNumber, j, k);
				Bright = getBunemanMagneticField(xnumberAdded - additionalBinNumber - 1, j, k);
			} else {
				Eleft = Efield[1 + additionalBinNumber][j][k];
				Eright = Efield[xnumberAdded - additionalBinNumber - 1][j][k];
				Bleft = Bfield[1 + additionalBinNumber][j][k];
				Bright = Bfield[xnumberAdded - additionalBinNumber - 1][j][k];
			}
			//theoreticalEnergy -= (Eright.vectorMult(Bright).x * deltaT * speed_of_light_normalized * deltaZ * deltaY * sqr(scaleFactor / plasma_period)) / (4 * pi);
			theoreticalEnergy -= (Eright.vectorMult(Bright).x * deltaT * deltaZ * deltaY * sqr(scaleFactor / plasma_period)) / (4 * pi);
			//theoreticalEnergy += (Eleft.vectorMult(Bleft).x * deltaT * speed_of_light_normalized * deltaZ * deltaY * sqr(scaleFactor / plasma_period)) / (4 * pi);
			theoreticalEnergy += (Eleft.vectorMult(Bleft).x * deltaT * deltaZ * deltaY * sqr(scaleFactor / plasma_period)) / (4 * pi);

			theoreticalMomentum -=
				(Eright.vectorMult(Bright)) * deltaT * deltaZ * deltaY *
				scaleFactor / plasma_period / (4 * pi);
			theoreticalMomentum +=
				(Eleft.vectorMult(Bleft)) * deltaT * deltaZ * deltaY * scaleFactor /
				plasma_period / (4 * pi);
		}
	}
	//}
	//}
	//MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("ipdate theoretical energy time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::updateEnergy() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	particleEnergy = 0;
	electricFieldEnergy = 0;
	electricFieldEnergyX = 0;
	electricFieldEnergyY = 0;
	electricFieldEnergyZ = 0;
	magneticFieldEnergy = 0;
	magneticFieldEnergyX = 0;
	magneticFieldEnergyY = 0;
	magneticFieldEnergyZ = 0;
	energy = 0;

	int minI = 1 + additionalBinNumber;
	if (cartCoord[0] == 0 && boundaryConditionTypeX != PERIODIC) {
		minI = 0;
	}
	int maxI = xnumberAdded - 1 - additionalBinNumber;
	if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
		maxI = xnumberAdded - 1;
	}
	int minJ = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - 1 - additionalBinNumber;
	int minK = 1 + additionalBinNumber;
	int maxK = znumberAdded - 1 - additionalBinNumber;
	if(ynumberGeneral == 1) {
		minJ = 0;
		maxJ = 1;
	}
	if(znumberGeneral == 1) {
		minK = 0;
		maxK = 1;
	}

	globalMomentum = Vector3d(0, 0, 0);
	electromagneticMomentum = Vector3d(0, 0, 0);
	particleMomentum = Vector3d(0, 0, 0);
	if ((currentIteration) % writeGeneralParameter == 0) {
		//particlesNumber = particles.size();
		for (int i = minI; i < maxI; ++i) {
			for (int j = minJ; j < maxJ; ++j) {
				for (int k = minK; k < maxK; ++k) {
					Vector3d E;
					if (solverType == BUNEMAN) {
						E = getBunemanElectricField(i, j, k);
					} else {
						E = Efield[i][j][k];
					}
					electricFieldEnergy += E.scalarMult(E) * volumeE() / (8 * pi);
					electricFieldEnergyX += E.x * E.x * volumeE() / (8 * pi);
					electricFieldEnergyY += E.y * E.y * volumeE() / (8 * pi);
					electricFieldEnergyZ += E.z * E.z * volumeE() / (8 * pi);
					//if (boundaryConditionTypeX == PERIODIC || rank > 0 || i > 1) {
					//electricFieldEnergy -= E0.scalarMult(E0) * volumeE() / (8 * pi);
					//}
				}
			}
		}

		for (int i = minI; i < maxI; ++i) {
			for (int j = minJ; j < maxJ; ++j) {
				for (int k = minK; k < maxK; ++k) {
					Vector3d B;
					if (solverType == BUNEMAN) {
						B = getBunemanMagneticField(i, j, k);
					} else {
						B = Bfield[i][j][k];
					}
					magneticFieldEnergy += (B.scalarMult(B)) * volumeB() / (8 * pi);
					magneticFieldEnergyX += B.x * B.x * volumeB() / (8 * pi);
					magneticFieldEnergyY += B.y * B.y * volumeB() / (8 * pi);
					magneticFieldEnergyZ += B.z * B.z * volumeB() / (8 * pi);
					//if (boundaryConditionTypeX == PERIODIC || rank > 0 || i > 0) {
					//magneticFieldEnergy -= B0.scalarMult(B0) * volumeB() / (8 * pi);
					//}
				}
			}
		}

		for (int i = minI; i < maxI; ++i) {
			for (int j = minJ; j < maxJ; ++j) {
				for (int k = minK; k < maxK; ++k) {
					Vector3d E;
					if (solverType == BUNEMAN) {
						E = getBunemanElectricField(i, j, k);
					} else {
						E = ((Efield[i][j][k] + Efield[i][j][k + 1] + Efield[i][j + 1][k] + Efield[i][j + 1][k + 1] + Efield[i + 1][j][k]
							+ Efield[i + 1][j][k + 1] + Efield[i + 1][j + 1][k] + Efield[i + 1][j + 1][k + 1]) / 8);
					}
					Vector3d B;
					if (solverType == BUNEMAN) {
						B = getBunemanMagneticField(i, j, k);
					} else {
						B = Bfield[i][j][k];
					}

					//globalMomentum += (E.vectorMult(B) / (4 * pi * speed_of_light_normalized)) * volumeB();
					globalMomentum += (E.vectorMult(B) / (4 * pi)) * volumeB();
					//electromagneticMomentum += (E.vectorMult(B) / (4 * pi * speed_of_light_normalized)) * volumeB();
					electromagneticMomentum += (E.vectorMult(B) / (4 * pi)) * volumeB();
				}
			}
		}

		for (int pcount = 0; pcount < particles.size(); ++pcount) {
			double particleE = particles[pcount]->energy();
			if (particleE < 0) {
				printf("particle energy < o\n");
				particleE = particles[pcount]->energy();
			}
			particleEnergy += particleE * particles[pcount]->weight;
			globalMomentum += particles[pcount]->getMomentum() * particles[pcount]->weight;
			particleMomentum += particles[pcount]->getMomentum() * particles[pcount]->weight;
		}


		particleEnergy *= sqr(scaleFactor / plasma_period);
		electricFieldEnergy *= sqr(scaleFactor / plasma_period);
		electricFieldEnergyX *= sqr(scaleFactor / plasma_period);
		electricFieldEnergyY *= sqr(scaleFactor / plasma_period);
		electricFieldEnergyZ *= sqr(scaleFactor / plasma_period);
		magneticFieldEnergy *= sqr(scaleFactor / plasma_period);
		magneticFieldEnergyX *= sqr(scaleFactor / plasma_period);
		magneticFieldEnergyY *= sqr(scaleFactor / plasma_period);
		magneticFieldEnergyZ *= sqr(scaleFactor / plasma_period);
		globalMomentum = globalMomentum * scaleFactor / plasma_period;
		particleMomentum = particleMomentum * scaleFactor / plasma_period;
		electromagneticMomentum = electromagneticMomentum * scaleFactor / plasma_period;


		energy = particleEnergy + electricFieldEnergy + magneticFieldEnergy;
	}

	if (currentIteration == 0) {
		theoreticalEnergy = energy;
		theoreticalMomentum = globalMomentum;
	}

	if ((currentIteration) % writeGeneralParameter == 0) {
		double buffer[23];
		if (rank != 0) {
			buffer[0] = particleEnergy;
			buffer[1] = electricFieldEnergy;
			buffer[2] = electricFieldEnergyX;
			buffer[3] = electricFieldEnergyY;
			buffer[4] = electricFieldEnergyZ;
			buffer[5] = magneticFieldEnergy;
			buffer[6] = magneticFieldEnergyX;
			buffer[7] = magneticFieldEnergyY;
			buffer[8] = magneticFieldEnergyZ;
			buffer[9] = energy;
			buffer[10] = globalMomentum.x;
			buffer[11] = globalMomentum.y;
			buffer[12] = globalMomentum.z;
			buffer[13] = theoreticalEnergy;
			buffer[14] = theoreticalMomentum.x;
			buffer[15] = theoreticalMomentum.y;
			buffer[16] = theoreticalMomentum.z;
			buffer[17] = electromagneticMomentum.x;
			buffer[18] = electromagneticMomentum.y;
			buffer[19] = electromagneticMomentum.z;
			buffer[20] = particleMomentum.x;
			buffer[21] = particleMomentum.y;
			buffer[22] = particleMomentum.z;
			//buffer[11] = (double)particlesNumber;

			MPI_Send(buffer, 23, MPI_DOUBLE, 0, MPI_SEND_GENERAL_PARAMETERS, cartComm);
		} else {
			generalTheoreticalEnergy = theoreticalEnergy;
			generalTheoreticalMomentum = theoreticalMomentum;
			for (int i = 1; i < nprocs; ++i) {
				MPI_Status status;
				MPI_Recv(buffer, 23, MPI_DOUBLE, i, MPI_SEND_GENERAL_PARAMETERS, cartComm, &status);
				particleEnergy += buffer[0];
				electricFieldEnergy += buffer[1];
				electricFieldEnergyX += buffer[2];
				electricFieldEnergyY += buffer[3];
				electricFieldEnergyZ += buffer[4];
				magneticFieldEnergy += buffer[5];
				magneticFieldEnergyX += buffer[6];
				magneticFieldEnergyY += buffer[7];
				magneticFieldEnergyZ += buffer[8];
				energy += buffer[9];
				globalMomentum.x += buffer[10];
				globalMomentum.y += buffer[11];
				globalMomentum.z += buffer[12];
				generalTheoreticalEnergy += buffer[13];
				generalTheoreticalMomentum.x += buffer[14];
				generalTheoreticalMomentum.y += buffer[15];
				generalTheoreticalMomentum.z += buffer[16];
				electromagneticMomentum.x += buffer[17];
				electromagneticMomentum.y += buffer[18];
				electromagneticMomentum.z += buffer[19];
				particleMomentum.x += buffer[20];
				particleMomentum.y += buffer[21];
				particleMomentum.z += buffer[22];
				//particlesNumber += buffer[11];
			}
		}
	}

	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating energy time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::updateParameters() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if (solverType == BUNEMAN) {
		maxBfield = Vector3d(bunemanBx[0][0][0], bunemanBy[0][0][0], bunemanBz[0][0][0]) - B0;
		maxEfield = Vector3d(bunemanEx[0][0][0], bunemanEy[0][0][0], bunemanEz[0][0][0]);

		int minI = 1 + additionalBinNumber;
		if (cartCoord[0] == 0 && boundaryConditionTypeX != PERIODIC) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1 - additionalBinNumber;
		if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
			maxI = xnumberAdded - 1;
		}
		int minJ = 1 + additionalBinNumber;
		int maxJ = ynumberAdded - 1 - additionalBinNumber;
		if(ynumberGeneral == 1) {
			minJ = 0;
			maxJ = 1;
		}
		int minK = 1 + additionalBinNumber;
		int maxK = znumberAdded - 1 - additionalBinNumber;
		if(znumberGeneral == 1) {
			minK = 0;
			maxK = 1;
		}
		for (int i = minI; i < maxI; ++i) {
			for (int j = minJ; j < maxJ; ++j) {
				for (int k = minK; k < maxK; ++k) {
					Vector3d normalField = Vector3d(bunemanBx[i][j][k], bunemanBy[i][j][k], bunemanBz[i][j][k]);
					//if (Bfield[i][j][k].norm() > maxBfield.norm()) {
					if (normalField.norm() > maxBfield.norm()) {
						//maxBfield = Bfield[i][j][k];
						maxBfield = normalField;
					}
				}
			}
		}

		meanSquaredEfield[0] = 0;
		meanSquaredEfield[1] = 0;
		meanSquaredEfield[2] = 0;
		for (int i = minI; i < maxI; ++i) {
			for (int j = minJ; j < maxJ; ++j) {
				for (int k = minK; k < maxK; ++k) {
					Vector3d tempE = Vector3d(bunemanEx[i][j][k], bunemanEy[i][j][k], bunemanEz[i][j][k]);
					if (tempE.norm() > maxEfield.norm()) {
						maxEfield = tempE;
					}
					//if (i > 1 + additionalBinNumber && i < xnumberAdded - additionalBinNumber - 1) {
					meanSquaredEfield[0] += bunemanEx[i][j][k] * bunemanEx[i][j][k];
					meanSquaredEfield[1] += bunemanEy[i][j][k] * bunemanEy[i][j][k];
					meanSquaredEfield[2] += bunemanEz[i][j][k] * bunemanEz[i][j][k];
					//}
				}
			}
		}
	} else {
		maxBfield = Bfield[0][0][0] - B0;
		maxEfield = Efield[0][0][0];

		int minI = 1 + additionalBinNumber;
		if (cartCoord[0] == 0 && boundaryConditionTypeX != PERIODIC) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1 - additionalBinNumber;
		if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
			maxI = xnumberAdded - 1;
		}
		int minJ = 1 + additionalBinNumber;
		int maxJ = ynumberAdded - 1 - additionalBinNumber;
		if(ynumberGeneral == 1) {
			minJ = 0;
			maxJ = 1;
		}
		int minK = 1 + additionalBinNumber;
		int maxK = znumberAdded - 1 - additionalBinNumber;
		if(znumberGeneral == 1) {
			minK = 0;
			maxK = 1;
		}

		for (int i = minI; i < maxI; ++i) {
			for (int j = minJ; j < maxJ; ++j) {
				for (int k = minK; k < maxK; ++k) {
					Vector3d normalField = Bfield[i][j][k];
					//if (Bfield[i][j][k].norm() > maxBfield.norm()) {
					if (normalField.norm() > maxBfield.norm()) {
						//maxBfield = Bfield[i][j][k];
						maxBfield = normalField;
					}
				}
			}
		}

		meanSquaredEfield[0] = 0;
		meanSquaredEfield[1] = 0;
		meanSquaredEfield[2] = 0;
		for (int i = minI; i < maxI; ++i) {
			for (int j = minJ; j < maxJ; ++j) {
				for (int k = minK; k < maxK; ++k) {
					if (Efield[i][j][k].norm() > maxEfield.norm()) {
						maxEfield = Efield[i][j][k];
					}
					//if (i > 1 + additionalBinNumber && i < xnumberAdded - additionalBinNumber - 1) {
					meanSquaredEfield[0] += Efield[i][j][k].x * Efield[i][j][k].x;
					meanSquaredEfield[1] += Efield[i][j][k].y * Efield[i][j][k].y;
					meanSquaredEfield[2] += Efield[i][j][k].z * Efield[i][j][k].z;
					//}
				}
			}
		}
	}
	double tempEsquared[3];

	MPI_Allreduce(meanSquaredEfield, tempEsquared, 3, MPI_DOUBLE, MPI_SUM, cartComm);

	meanSquaredEfield[0] = sqrt(tempEsquared[0] / (xnumberGeneral * ynumberGeneral * znumberGeneral));
	meanSquaredEfield[1] = sqrt(tempEsquared[1] / (xnumberGeneral * ynumberGeneral * znumberGeneral));
	meanSquaredEfield[2] = sqrt(tempEsquared[2] / (xnumberGeneral * ynumberGeneral * znumberGeneral));
	MPI_Barrier(cartComm);
	double temp[7];
	temp[0] = maxBfield.x;
	temp[1] = maxBfield.y;
	temp[2] = maxBfield.z;
	temp[3] = maxEfield.x;
	temp[4] = maxEfield.y;
	temp[5] = maxEfield.z;
	if (rank != 0) {
		MPI_Send(temp, 6, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, cartComm);
	} else {
		for (int i = 1; i < nprocs; ++i) {
			MPI_Status status;
			MPI_Recv(temp, 6, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, cartComm, &status);
			Vector3d tempB = Vector3d(temp[0], temp[1], temp[2]);
			Vector3d tempE = Vector3d(temp[3], temp[4], temp[5]);
			if (tempB.norm() > maxBfield.norm()) {
				maxBfield = tempB;
			}
			if (tempE.norm() > maxEfield.norm()) {
				maxEfield = tempE;
			}
		}
	}
	MPI_Barrier(cartComm);

	if (rank == 0) {
		temp[0] = maxBfield.x;
		temp[1] = maxBfield.y;
		temp[2] = maxBfield.z;
		temp[3] = maxEfield.x;
		temp[4] = maxEfield.y;
		temp[5] = maxEfield.z;
		for (int i = 1; i < nprocs; ++i) {
			MPI_Send(temp, 6, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, cartComm);
		}
	} else {
		MPI_Status status;
		MPI_Recv(temp, 6, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, cartComm, &status);
		maxBfield = Vector3d(temp[0], temp[1], temp[2]);
		maxEfield = Vector3d(temp[3], temp[4], temp[5]);
	}
	//MPI_Barrier(cartComm);

	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating parameters time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::updateAnisotropy() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	double* parallelV2 = new double[typesNumber];
	double* normalV2 = new double[typesNumber];
	double* inBuffer = new double[3 * typesNumber];
	double* outBuffer = new double[3 * typesNumber];
	for (int t = 0; t < typesNumber; ++t) {
		parallelV2[t] = 0;
		normalV2[t] = 0;
		types[t].generalWeight = 0;
	}
	for (int p = 0; p < particles.size(); ++p) {
		Particle* particle = particles[p];
		int t = getTypeNumber(particle);
		types[t].generalWeight += particle->weight;
		parallelV2[t] = parallelV2[t] + sqr(particle->velocityX()) * particle->weight;
		normalV2[t] = normalV2[t] + (sqr(particle->velocityY())
			+ sqr(particle->velocityZ())) * particle->weight;
	}

	MPI_Barrier(cartComm);

	if (nprocs > 1) {
		if (rank != 0) {
			for (int i = 0; i < typesNumber; ++i) {
				outBuffer[3 * i] = parallelV2[i];
				outBuffer[3 * i + 1] = normalV2[i];
				outBuffer[3 * i + 2] = types[i].generalWeight;
			}
			MPI_Send(outBuffer, 3 * typesNumber, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_ALL_TO_FIRST, cartComm);
		} else {
			for (int i = 1; i < nprocs; ++i) {
				MPI_Status status;
				MPI_Recv(inBuffer, 3 * typesNumber, MPI_DOUBLE, i, MPI_SEND_DOUBLE_ALL_TO_FIRST, cartComm, &status);
				for (int t = 0; t < typesNumber; ++t) {
					parallelV2[t] += inBuffer[3 * t];
					normalV2[t] += inBuffer[3 * t + 1];
					types[t].generalWeight += inBuffer[3 * t + 2];
				}
			}
		}
	}

	MPI_Barrier(cartComm);

	for (int t = 0; t < typesNumber; ++t) {
		if (types[t].generalWeight > 0) {
			types[t].anisotropy = (0.5 * normalV2[t] / parallelV2[t]) - (2 * parallelV2[t] / normalV2[t]);
			types[t].parallelTemperatureEvaluated = types[t].mass * parallelV2[t] / (kBoltzman_normalized * types[t].
				generalWeight);
			types[t].normalTemperatureEvaluated = 0.5 * types[t].mass * normalV2[t] / (kBoltzman_normalized * types[t].
				generalWeight);
		} else {
			types[t].anisotropy = 0;
			types[t].parallelTemperatureEvaluated = 0;
			types[t].normalTemperatureEvaluated = 0;
		}
	}

	delete[] parallelV2;
	delete[] normalV2;
	delete[] inBuffer;
	delete[] outBuffer;
	//MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating anisotropy time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::resetNewTempFields() {
	if(solverType == BUNEMAN) {
		for(int i = 0; i < xnumberAdded; ++i) {
			for(int j = 0; j <= ynumberAdded; ++j) {
				for(int k = 0; k <= znumberAdded; ++k) {
					bunemanNewEx[i][j][k] = bunemanEx[i][j][k];
				}
			}
		}

		for(int i = 0; i <= xnumberAdded; ++i) {
			for(int j = 0; j < ynumberAdded; ++j) {
				for(int k = 0; k <= znumberAdded; ++k) {
					bunemanNewEy[i][j][k] = bunemanEy[i][j][k];
				}
			}
		}

		for(int i = 0; i <= xnumberAdded; ++i) {
			for(int j = 0; j <= ynumberAdded; ++j) {
				for(int k = 0; k <= znumberAdded; ++k) {
					bunemanNewEz[i][j][k] = bunemanEz[i][j][k];
				}
			}
		}

		for(int i = 0; i <= xnumberAdded; ++i) {
			for(int j = 0; j < ynumberAdded; ++j) {
				for(int k = 0; k < znumberAdded; ++k) {
					bunemanNewBx[i][j][k] = bunemanBx[i][j][k];
				}
			}
		}

		for(int i = 0; i < xnumberAdded; ++i) {
			for(int j = 0; j <= ynumberAdded; ++j) {
				for(int k = 0; k < znumberAdded; ++k) {
					bunemanNewBy[i][j][k] = bunemanBy[i][j][k];
				}
			}
		}

		for(int i = 0; i < xnumberAdded; ++i) {
			for(int j = 0; j < ynumberAdded; ++j) {
				for(int k = 0; k <= znumberAdded; ++k) {
					bunemanNewBz[i][j][k] = bunemanBz[i][j][k];
				}
			}
		}
		
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					newBfield[i][j][k] = Bfield[i][j][k];
				}
			}
		}

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					tempEfield[i][j][k] = Efield[i][j][k];
					newEfield[i][j][k] = Efield[i][j][k];
					explicitEfield[i][j][k] = Efield[i][j][k];
				}
			}
		}
	}
}

void Simulation::updateShockWaveX() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}

	double tempShockWavex[1];
	tempShockWavex[0] = -1.0;
	if (rank < nprocs - 1) {
		MPI_Status status;
		MPI_Recv(tempShockWavex, 1, MPI_DOUBLE, rank + 1, MPI_SEND_DOUBLE_NUMBER_LEFT, cartComm, &status);
		if (tempShockWavex[0] < 0) {
			for (int i = xnumberAdded - 1; i > 0; i--) {
				if (particleBulkVelocities[1][i][0][0].x > V0.x / 2) {
					tempShockWavex[0] = xgrid[i + 1];
					break;
				}
			}
		}
		if (rank > 0) {
			MPI_Send(tempShockWavex, 1, MPI_DOUBLE, rank - 1, MPI_SEND_DOUBLE_NUMBER_LEFT, cartComm);
		}
	} else {
		for (int i = xnumberAdded - 1; i > 0; i--) {
			if (particleBulkVelocities[1][i][0][0].x > V0.x / 2) {
				tempShockWavex[0] = xgrid[i + 1];
				break;
			}
		}
		if (rank > 0) {
			MPI_Send(tempShockWavex, 1, MPI_DOUBLE, rank - 1, MPI_SEND_DOUBLE_NUMBER_LEFT, cartComm);
		}
	}

	if (tempShockWavex[0] < shockWaveX) {
		tempShockWavex[0] = shockWaveX;
	}

	MPI_Barrier(cartComm);

	MPI_Bcast(tempShockWavex, 1, MPI_DOUBLE, 0, cartComm);

	shockWaveX = tempShockWavex[0];

	//MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating shock wave = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::splitParticles() {
	int tempParticlesNumber[1];
	for (int i = 0; i < nprocs; ++i) {
		if (i == rank) {
			if (rank > 0) {
				MPI_Status status;
				MPI_Recv(tempParticlesNumber, 1, MPI_INT, rank - 1, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
				particlesNumber = tempParticlesNumber[0];
			}
			int n = particles.size();
			for (int p = 0; p < n; ++p) {
				Particle* particle = particles[p];
				/*if (particle->getMomentum().norm() > particleSplitLevel * particle->prevMomentum.norm()) {
					splitParticle(particle);
					particlesNumber++;
				}*/
			}
			if (rank < nprocs - 1) {
				tempParticlesNumber[0] = particlesNumber;
				MPI_Send(tempParticlesNumber, 1, MPI_INT, rank + 1, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm);
			}
		}
	}
	tempParticlesNumber[0] = particlesNumber;
	MPI_Bcast(tempParticlesNumber, 1, MPI_INT, nprocs - 1, cartComm);
	particlesNumber = tempParticlesNumber[0];
}

void Simulation::splitParticle(Particle* particle) {
	Vector3d momentum = particle->getMomentum();
	//particle->prevMomentum = momentum;
	particle->weight /= 2.0;
	Particle* newParticle = new Particle(*particle);
	newParticle->number = particlesNumber;

	Vector3d xort = Vector3d(1, 0, 0);
	Vector3d yort = Vector3d(0, 1, 0);
	Vector3d zort = Vector3d(0, 0, 1);

	double norm = momentum.norm();

	if (norm > 1E-100) {
		Vector3d tempXort = momentum / norm;
		Vector3d tempYort = momentum.vectorMult(xort);
		double normY = tempYort.norm();
		if (normY < 1E-100) {
			tempYort = momentum.vectorMult(yort);
			normY = tempYort.norm();
		}
		if (normY < 1E-100) {
			tempYort = momentum.vectorMult(zort);
			normY = tempYort.norm();
		}
		if (normY < 1E-100) {
			printf("something wrong with normY =  %g norm p = %g\n", normY, norm);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "somtehing wrong with normY = %g norm p = %g\n", normY, norm);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}
		tempYort = tempYort / normY;
		Vector3d tempZort = tempXort.vectorMult(tempYort);
		double normZ = tempZort.norm();
		if (normZ < 1E-100) {
			printf("something wrong with normZ =  %g norm p = %g\n", normZ, norm);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "somtehing wrong with normZ = %g norm p = %g\n", normZ, norm);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}
		tempZort = tempZort / normZ;

		double phi = 2 * pi * uniformDistribution();

		double x = norm * cos(particleSplitAngle);
		double y = norm * sin(particleSplitAngle) * cos(phi);
		double z = norm * sin(particleSplitAngle) * sin(phi);

		Vector3d newMomentum1 = (tempXort * x) + (tempYort * y) + (tempZort * z);
		Vector3d newMomentum2 = (tempXort * x) + (tempYort * y) + (tempZort * (-z));

		particle->setMomentum(newMomentum1);
		newParticle->setMomentum(newMomentum2);
	}

	particles.push_back(newParticle);
}

