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


void Simulation::simulate() {
	MPI_Barrier(cartComm);
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
	}

	outputGeneralInitialParameters((outputDir + "initialParameters.dat").c_str(), (outputDir + "initialParametersWithText.dat").c_str(), this);


	MPI_Barrier(cartComm);

	if ((rank == 0) && (verbosity > 1)) printf("initial exchanging fields\n");
	exchangeEfield();
	exchangeGeneralBfield(Bfield);
	exchangeGeneralBfield(newBfield);

	if ((rank == 0) && (verbosity > 1)) printf("initial collecting particles\n");
	updateParticleCorrelationMaps();
	updateParameters();
	//updateAnisotropy();

	//double because dielectric tensor needs deltaT;
	updateDeltaT();
	evaluateParticlesRotationTensor();
	updateElectroMagneticParameters();
	updateDensityParameters();
	updateDeltaT();
	evaluateParticlesRotationTensor();
	updateElectroMagneticParameters();
	updateDensityParameters();
	if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
		//updateShockWaveX();
	}

	evaluateExplicitDerivative();
	if (currentIteration % divergenceCleanUpParameter == 0){
		if(solverType == BUNEMAN){
			//cleanupDivergenceBuneman();
			//cleanUpDivergenceBunemanMagnetic();
		} else {
			//cleanupDivergence(newEfield, chargeDensity);
			//cleanupDivergenceMagnetic();
		}
	}
	updateFields();
	exchangeGeneralBfield(Bfield);
	//updateEnergy();

	for (int i = 0; i < typesNumber; ++i) {
		types[i].injectionLength = types[i].particesDeltaX - 0.0001 * deltaX;
	}

	while (time * plasma_period < maxTime && currentIteration < maxIteration) {
		if(currentIteration % writeMemoryParameter == 0){
			outputMemory((outputDir + "memory_using.dat").c_str(), cartComm, cartCoord, cartDim);
		}
		if ((rank == 0) && (verbosity > 0)) {
			FILE* logFile = fopen((outputDir + "log.dat").c_str(), "a");
			fprintf(logFile, "start iteration number = %d time = %15.10g\n", currentIteration, time);
			fflush(logFile);
			fclose(logFile);
		}
		if ((rank == 0)) printf("start iteration number = %d time = %15.10g\n", currentIteration, time);
		updateEnergy();
		if (currentIteration % writeParameter == 0) {
			output();
		}
		if(currentIteration % writeGeneralParameter == 0) {
			if (rank == 0) outputGeneral((outputDir + "general.dat").c_str(), this);
		}

		if (currentIteration % writeTrajectoryNumber == 0) {
			outputTrajectories();
		}

		double iterationTime = 0;
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			iterationTime = clock();
		}

		//updateDeltaT();

 		evaluateParticlesRotationTensor();

		updateElectroMagneticParameters();
		//smoothChargeDensityHat();

		evaluateElectricField();
		exchangeEfield();
		if(solverType == BUNEMAN){
		exchangeBunemanEfield(bunemanEx, bunemanEy, bunemanEz);
		exchangeBunemanEfield(bunemanNewEx, bunemanNewEy, bunemanNewEz);
		}

		if(solverType == BUNEMAN){
		} else if(solverType == IMPLICIT || solverType == IMPLICIT_EC){
			/*cleanupDivergence(tempEfield, chargeDensityHat);
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						newEfield[i][j][k] = (tempEfield[i][j][k] - Efield[i][j][k] * (1 - theta)) / theta;
					}
				}
			}
			exchangeEfield();*/
		} else {
		}

		evaluateMagneticField();
		if(solverType == BUNEMAN){
		exchangeBunemanBfield(bunemanBx, bunemanBy, bunemanBz);
		exchangeBunemanBfield(bunemanNewBx, bunemanNewBy, bunemanNewBz);
		}
		double procTime = 0;
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock();
		}
		exchangeGeneralBfield(newBfield);
		exchangeGeneralEfield(tempEfield);
		exchangeGeneralEfield(newEfield);
		MPI_Barrier(cartComm);
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("exchanging field time = %g sec\n", procTime / CLOCKS_PER_SEC);
		}
		

		moveParticles();

		if ((rank == 0) && (verbosity > 1)) printLog("erasing escaped particles\n");
		if ((rank == 0) && (verbosity > 1)) printf("erasing escaped particles\n");

		eraseEscapedPaticles();

		//MPI_Barrier(cartComm);
		if ((rank == 0) && (verbosity > 1)) printLog("start exchange particles\n");
		if ((rank == 0) && (verbosity > 1)) printf("start exchange particles\n");
		MPI_Barrier(cartComm);

		exchangeParticles();

		if ((rank == 0) && (verbosity > 1)) printLog("start injecting new particles\n");
		if ((rank == 0) && (verbosity > 1)) printf("start injecting new particles\n");


		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock();
		}
		for (int i = 0; i < typesNumber; ++i) {
			types[i].injectionLength += fabs(V0.x * deltaT);
		}
		if (boundaryConditionType == SUPER_CONDUCTOR_LEFT || boundaryConditionType == FREE_BOTH) {
			for(int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
				for(int cartK = 0; cartK < cartDim[2]; ++cartK){
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
		if (preserveChargeGlobal && (boundaryConditionType != PERIODIC)) {
			addToPreserveChargeGlobal();
		}
		MPI_Barrier(cartComm);
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("injecting and preserving time = %g sec\n", procTime / CLOCKS_PER_SEC);
		}

		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock();
		}
		updateParticleCorrelationMaps();
		MPI_Barrier(cartComm);
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("updating correlation maps = %g sec\n", procTime / CLOCKS_PER_SEC);
		}

		updateDensityParameters();
		if((rank == 0) && (verbosity > 0)) {
			printf("finish update density parameters\n");
		}
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock();
		}
		exchangeGeneralEfield(newEfield);
		if((rank == 0) && (verbosity > 0)) {
			printf("finish exchange new efield\n");
		}
		MPI_Barrier(cartComm);
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("exchanging fields time = %g sec\n", procTime / CLOCKS_PER_SEC);
		}

		if(solverType == BUNEMAN){
			//smoothBunemanEfieldGeneral(bunemanNewEx, bunemanNewEy, bunemanNewEz);
			//smoothBunemanBfieldGeneral(bunemanNewBx, bunemanNewBy, bunemanNewBz);
		} else {
			//smoothNewEfield();
			//smoothNewBfield();
		}
		if (currentIteration % divergenceCleanUpParameter == 0) {
			if(solverType == BUNEMAN){
				//cleanupDivergenceBuneman();
				//cleanupDivergenceBunemanMagnetic();
			} else {
				//cleanupDivergence(newEfield, chargeDensity);
				//cleanupDivergenceMagnetic();
			}
		}

		filterFields(5);
		//filterFieldsLocal(5);

		updateFields();
		if((rank == 0) && (verbosity > 0)) {
			printf("finish update fields\n");
		}

		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock();
		}
		exchangeGeneralEfield(newEfield);
		exchangeGeneralEfield(Efield);
		exchangeGeneralBfield(Bfield);
		if(solverType == BUNEMAN){
		exchangeBunemanEfield(bunemanEx, bunemanEy, bunemanEz);
		exchangeBunemanEfield(bunemanNewEx, bunemanNewEy, bunemanNewEz);
		exchangeBunemanBfield(bunemanBx, bunemanBy, bunemanBz);
		exchangeBunemanBfield(bunemanNewBx, bunemanNewBy, bunemanNewBz);
		}

		/*if(solverType == BUNEMAN){
			smoothBunemanEfieldGeneral(bunemanEx, bunemanEy, bunemanEz);
			smoothBunemanBfieldGeneral(bunemanBx, bunemanBy, bunemanBz);
		}*/
		if((rank == 0) && (verbosity > 0)) {
			printf("finish exchange fields\n");
		}
		MPI_Barrier(cartComm);
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("exchanging fields time = %g sec\n", procTime / CLOCKS_PER_SEC);
		}

		//smoothNewEfield();

		//updateEnergy();

		if ((currentIteration + 1) % writeGeneralParameter == 0) {
			updateParameters();
			updateAnisotropy();
		}

		if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
			//todo
			//updateShockWaveX();
		}

		removeEscapedParticles();
		if((rank == 0) && (verbosity > 0)) {
			printf("finish remove escaped particles\n");
		}

		/*if(currentIteration % splitParticlesParameter == 0){
		    if ((rank == 0) && (verbosity > 0)) printf("split particles\n");
		    splitParticles();
		}*/

		time += deltaT;

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
	outputParticlesTrajectories((outputDir + "particlesTrajectories.dat").c_str(), (outputDir + "electronsTrajectories.dat").c_str(), particles, trackedParticlesNumbers, trackedParticlesNumber, time, plasma_period, scaleFactor, this);

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
	if ((rank == 0) && (verbosity > 0)) printf("outputing iteration number %d\n", currentIteration);
	if ((rank == 0) && (verbosity > 0)) printLog("collecting most accelerate particles\n");
	collectMostAcceleratedParticles();
	if ((rank == 0) && (verbosity > 0)) printLog("outputing\n");

	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution protons\n");
	outputDistribution((outputDir + "distribution_protons.dat").c_str(), particles, PROTON, scaleFactor,
	                   plasma_period, verbosity);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution electrons\n");
	outputDistribution((outputDir + "distribution_electrons.dat").c_str(), particles, ELECTRON, scaleFactor,
	                   plasma_period, verbosity);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution alphas\n");
	outputDistribution((outputDir + "distribution_alphas.dat").c_str(), particles, ALPHA, scaleFactor,
	                   plasma_period, verbosity);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution positrons\n");
	outputDistribution((outputDir + "distribution_positrons.dat").c_str(), particles, POSITRON, scaleFactor,
	                   plasma_period, verbosity);

	Vector3d shockWaveV = V0 / 3;

	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution protons shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_protons_sw.dat").c_str(), particles, shockWaveV, speed_of_light_normalized, PROTON, scaleFactor,
	                                plasma_period, verbosity);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution electrons shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_electrons_sw.dat").c_str(), particles, shockWaveV, speed_of_light_normalized, ELECTRON, scaleFactor,
	                                plasma_period, verbosity);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution alphas shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_alphas_sw.dat").c_str(), particles, shockWaveV, speed_of_light_normalized, ALPHA, scaleFactor,
	                                plasma_period, verbosity);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution positrons shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_positrons_sw.dat").c_str(), particles, shockWaveV, speed_of_light_normalized, POSITRON, scaleFactor,
	                                plasma_period, verbosity);

	//if ((rank == 0) && (verbosity > 1)) printf("outputing anisotropy electrons\n");
	//outputAnisotropy((outputDir + "anisotropy_electrons.dat").c_str(), this, ELECTRON, scaleFactor, plasma_period);

	//anisotropyFileProton = fopen((outputDir + "anisotropy_protons.dat").c_str(), "a");
	//if ((rank == 0) && (verbosity > 1)) printf("outputing anisotropy protons\n");

	//if ((rank == 0) && (verbosity > 1)) printf("outputing anisotropy alphas\n");
	//outputAnisotropy((outputDir + "anisotropy_alphas.dat").c_str(), this, ALPHA, scaleFactor, plasma_period);

	//if ((rank == 0) && (verbosity > 1)) printf("outputing anisotropy positrons\n");
	//outputAnisotropy((outputDir + "anisotropy_positrons.dat").c_str(), this, POSITRON, scaleFactor, plasma_period);

	//if ((rank == 0) && (verbosity > 1)) printf("outputing anisotropy deuterium\n");
	//outputAnisotropy((outputDir + "anisotropy_deuterium.dat").c_str(), this, DEUTERIUM, scaleFactor, plasma_period);
	//if ((rank == 0) && (verbosity > 1)) printf("outputing anisotropy helium3\n");
	//outputAnisotropy((outputDir + "anisotropy_helium3.dat").c_str(), this, HELIUM3, scaleFactor, plasma_period);
	//if ((rank == 0) && (verbosity > 1)) printf("outputing anisotropy oxygen+3\n");
	//outputAnisotropy((outputDir + "anisotropy_oxygen+3.dat").c_str(), this, OXYGEN_PLUS3, scaleFactor, plasma_period);
	//if ((rank == 0) && (verbosity > 1)) printf("outputing anisotropy silicon\n");
	//outputAnisotropy((outputDir + "anisotropy_silicon.dat").c_str(), this, SILICON_PLUS1, scaleFactor, plasma_period);
	//}

	if ((rank == 0) && (verbosity > 1)) printf("outputing fields\n");
	outputFields((outputDir + "Efield.dat").c_str(), (outputDir + "Bfield.dat").c_str(), Efield, Bfield, xnumberAdded,
	             ynumberAdded, znumberAdded, additionalBinNumber, plasma_period, scaleFactor, cartComm, cartCoord, cartDim);
	/*outputFieldsReduced((outputDir + "EfieldReduced.dat").c_str(), (outputDir + "BfieldReduced.dat").c_str(), Efield, Bfield, xnumberAdded,
	                    ynumberAdded, znumberAdded, additionalBinNumber, reduceStepX, reduceStepY, reduceStepZ, plasma_period, scaleFactor);*/
	/*for(int curI = 0; curI < cartDim[0]; ++curI){
		for(int curJ = 0; curJ < cartDim[1]; ++curJ){
			for(int curK = 0; curK < cartDim[2]; ++curK){
				if(curI == cartCoord[0] && curJ == cartCoord[1] && curK == cartCoord[2]){
					std::string fileName = outputDir + "Xfile_" + convertIntToString(curI) +"_" + convertIntToString(curJ) + "_" + convertIntToString(curK) + ".dat";
					outputGridSimple(fileName.c_str(), xgrid, xnumberAdded, scaleFactor);
					fileName = outputDir + "Yfile_" + convertIntToString(curI) +"_" + convertIntToString(curJ) + "_" + convertIntToString(curK) + ".dat";
					outputGridSimple(fileName.c_str(), ygrid, ynumberAdded, scaleFactor);
					fileName = outputDir + "Zfile_" + convertIntToString(curI) +"_" + convertIntToString(curJ) + "_" + convertIntToString(curK) + ".dat";
					outputGridSimple(fileName.c_str(), zgrid, znumberAdded, scaleFactor);
				}
			}
		}
	}*/

	int dimsYZ[3];
	dimsYZ[0] = 0;
	dimsYZ[1] = 1;
	dimsYZ[2] = 1;
	MPI_Comm subCommYZ;
	MPI_Cart_sub(cartComm, dimsYZ, &subCommYZ);
	int dimsZ[3];
	dimsZ[0] = 0;
	dimsZ[1] = 0;
	dimsZ[2] = 1;
	MPI_Comm subCommZ;
	MPI_Cart_sub(cartComm, dimsZ, &subCommZ);
	int dimsXY[3];
	dimsXY[0] = 1;
	dimsXY[1] = 1;
	dimsXY[2] = 0;
	MPI_Comm subCommXY;
	MPI_Cart_sub(cartComm, dimsXY, &subCommXY);
	int dimsY[3];
	dimsY[0] = 0;
	dimsY[1] = 1;
	dimsY[2] = 0;
	MPI_Comm subCommY;
	MPI_Cart_sub(cartComm, dimsY, &subCommY);
	int dimsXZ[3];
	dimsXZ[0] = 1;
	dimsXZ[1] = 0;
	dimsXZ[2] = 1;
	MPI_Comm subCommXZ;
	MPI_Cart_sub(cartComm, dimsXZ, &subCommXZ);
	int dimsX[3];
	dimsX[0] = 1;
	dimsX[1] = 0;
	dimsX[2] = 0;
	MPI_Comm subCommX;
	MPI_Cart_sub(cartComm, dimsX, &subCommX);

	if(verbosity > 2) printf("get cart coord with absolute index rank = %d\n", rank);
	int coordX = getCartCoordWithAbsoluteIndexX(xnumberGeneral/2);
	int coordY = getCartCoordWithAbsoluteIndexY(ynumberGeneral/2);
	int coordZ = getCartCoordWithAbsoluteIndexZ(znumberGeneral/2);
	if(verbosity > 2) printf("x coord with absolute index = %d\n", coordX);
	MPI_Barrier(cartComm);
	if(verbosity > 2) printf("output crossection fields x\n");
	if(coordX == cartCoord[0]){
		int xindex = getLocalIndexByAbsoluteX(xnumberGeneral/2);
		if(verbosity > 2) printf("x local index = %d\n", xindex);
		outputFieldsCrossectionYZ((outputDir + "EfieldYZ.dat").c_str(), (outputDir + "BfieldYZ.dat").c_str(), Efield, Bfield, xnumberAdded,
	             ynumberAdded, znumberAdded, additionalBinNumber, plasma_period, scaleFactor, subCommYZ, subCommZ, cartCoord, cartDim, xindex);
	}
	MPI_Barrier(cartComm);
	if(verbosity > 2) printf("output crossection fields y\n");
	if(coordY == cartCoord[1]){
		int yindex = getLocalIndexByAbsoluteY(ynumberGeneral/2);
		if(verbosity > 2) printf("y local index = %d\n", yindex);
		outputFieldsCrossectionXZ((outputDir + "EfieldXZ.dat").c_str(), (outputDir + "BfieldXZ.dat").c_str(), Efield, Bfield, xnumberAdded,
	             ynumberAdded, znumberAdded, additionalBinNumber, plasma_period, scaleFactor, subCommXZ, subCommZ, cartCoord, cartDim, yindex);
	}
	MPI_Barrier(cartComm);
	if(verbosity > 2) printf("output crossection fields z\n");
	if(coordZ == cartCoord[2]){
		int zindex = getLocalIndexByAbsoluteZ(znumberGeneral/2);
		if(verbosity > 2) printf("z local index = %d\n", zindex);
		outputFieldsCrossectionXY((outputDir + "EfieldXY.dat").c_str(), (outputDir + "BfieldXY.dat").c_str(), Efield, Bfield, xnumberAdded,
	             ynumberAdded, znumberAdded, additionalBinNumber, plasma_period, scaleFactor, subCommXY, subCommY, cartCoord, cartDim, zindex);
	}

	if ((rank == 0) && (verbosity > 1)) printf("outputing grid\n");
	outputGridX((outputDir + "Xfile.dat").c_str(), xgrid, xnumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim, scaleFactor);
	outputGridReducedX((outputDir + "XfileReduced.dat").c_str(), xgrid, xnumberAdded, additionalBinNumber, reduceStepX, rank, leftRank, rightRank, cartComm, cartCoord, cartDim, scaleFactor);

	outputGridY((outputDir + "Yfile.dat").c_str(), ygrid, ynumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim, scaleFactor);
	outputGridReducedY((outputDir + "YfileReduced.dat").c_str(), ygrid, ynumberAdded, additionalBinNumber, reduceStepY, rank, frontRank, backRank, cartComm, cartCoord, cartDim, scaleFactor);

	outputGridZ((outputDir + "Zfile.dat").c_str(), zgrid, znumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim, scaleFactor);
	outputGridReducedZ((outputDir + "ZfileReduced.dat").c_str(), zgrid, znumberAdded, additionalBinNumber, reduceStepZ, rank, bottomRank, topRank, cartComm, cartCoord, cartDim, scaleFactor);


	//if ((rank == 0) && (verbosity > 1)) printf("outputing concentrations\n");
	outputConcentrations((outputDir + "concentrations.dat").c_str(), particleConcentrations, chargeDensity,
	                     chargeDensityHat, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period, scaleFactor, cartComm, cartCoord, cartDim);


	if ((rank == 0) && (verbosity > 1)) printf("outputing velocity\n");
	outputVelocity((outputDir + "velocity.dat").c_str(),
	               particleBulkVelocities, types, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, typesNumber, plasma_period, scaleFactor, cartComm, cartCoord, cartDim);

	if ((rank == 0) && (verbosity > 1)) printf("outputing flux\n");
	outputFlux((outputDir + "flux.dat").c_str(), electricFlux, externalElectricFlux, xnumberAdded, ynumberAdded,
	           znumberAdded, additionalBinNumber, plasma_period, scaleFactor, cartComm, cartCoord, cartDim);

	if ((rank == 0) && (verbosity > 1)) printf("outputing divergence\n");
	outputDivergenceError((outputDir + "divergence_error.dat").c_str(), this, plasma_period, scaleFactor);

	double rotBscale = 1.0 / (plasma_period * plasma_period * sqrt(scaleFactor));

	if ((rank == 0) && (verbosity > 1)) printf("outputing rotB\n");
	outputVectorNodeArray((outputDir + "rotBFile.dat").c_str(), rotB, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim, rotBscale);

	if ((rank == 0) && (verbosity > 1)) printf("outputing Ederivative\n");
	outputVectorNodeArray((outputDir + "EderivativeFile.dat").c_str(), Ederivative, xnumberAdded, ynumberAdded,
	                      znumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim, rotBscale);

	if ((rank == 0) && (verbosity > 1)) printf("outputing rotE\n");
	outputVectorCellArray((outputDir + "rotEFile.dat").c_str(), rotE, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim, rotBscale);

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
			mostAcceleratedParticlesNumbers[t].push_back(std::pair<int, double>(particle->number, momentum));
		} else {
			std::list<std::pair<int, double> >::iterator it = mostAcceleratedParticlesNumbers[t].begin();
			std::pair<int, double> pair = *it;
			if (momentum > pair.second) {
				bool added = false;
				for (it = mostAcceleratedParticlesNumbers[t].begin(); it != mostAcceleratedParticlesNumbers[t].end(); ++it) {
					pair = *it;
					if (momentum < pair.second) {
						mostAcceleratedParticlesNumbers[t].insert(it, std::pair<int, double>(particle->number, momentum));
						added = true;
						break;
					}
				}
				if (!added) {
					mostAcceleratedParticlesNumbers[t].push_back(std::pair<int, double>(particle->number, momentum));
				}
				if (mostAcceleratedParticlesNumbers[t].size() > mostFastParticlesNumber) {
					mostAcceleratedParticlesNumbers[t].pop_front();
				}
			} else if (mostAcceleratedParticlesNumbers[t].size() < mostFastParticlesNumber) {
				mostAcceleratedParticlesNumbers[t].push_front(std::pair<int, double>(particle->number, momentum));
			}
		}
	}

	for (int t = 0; t < typesNumber; ++t) {
		if (types[t].particlesPerBin > 0) {
			while (mostAcceleratedParticlesNumbers[t].size() < mostFastParticlesNumber) {
				mostAcceleratedParticlesNumbers[t].push_front(std::pair<int, double>(-1, 0));
			}
		}
	}

	for (int t = 0; t < typesNumber; ++t) {
		if (types[t].particlesPerBin > 0) {

			int tempAcceleratedNumbers[mostFastParticlesNumber];
			double tempAcceleratedMomentum[mostFastParticlesNumber];

			std::list<std::pair<int, double> >::iterator it = mostAcceleratedParticlesNumbers[t].begin();
			for (int i = 0; i < mostFastParticlesNumber; ++i) {
				std::pair<int, double> pair = *it;
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
					MPI_Recv(inAcceleratedNumbers, mostFastParticlesNumber, MPI_INT, i, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm, &status);
					MPI_Recv(inAcceleratedMomentum, mostFastParticlesNumber, MPI_DOUBLE, i, MPI_SEND_DOUBLE_ALL_TO_FIRST, cartComm, &status);

					for (int j = 0; j < mostFastParticlesNumber; ++j) {
						int curNumber = inAcceleratedNumbers[j];
						int momentum = inAcceleratedMomentum[j];
						if (mostAcceleratedParticlesNumbers[t].size() == 0) {
							mostAcceleratedParticlesNumbers[t].push_back(std::pair<int, double>(curNumber, momentum));
						} else {
							it = mostAcceleratedParticlesNumbers[t].begin();
							std::pair<int, double> pair = *it;
							if (momentum > pair.second) {
								bool added = false;
								for (it = mostAcceleratedParticlesNumbers[t].begin(); it != mostAcceleratedParticlesNumbers[t].end(); ++it) {
									pair = *it;
									if (momentum < pair.second) {
										mostAcceleratedParticlesNumbers[t].insert(it, std::pair<int, double>(curNumber, momentum));
										added = true;
										break;
									}
								}
								if (!added) {
									mostAcceleratedParticlesNumbers[t].push_back(std::pair<int, double>(curNumber, momentum));
								}
								if (mostAcceleratedParticlesNumbers[t].size() > mostFastParticlesNumber) {
									mostAcceleratedParticlesNumbers[t].pop_front();
								}
							} else if (mostAcceleratedParticlesNumbers[t].size() < mostFastParticlesNumber) {
								mostAcceleratedParticlesNumbers[t].push_front(std::pair<int, double>(curNumber, momentum));
							}
						}
					}
				}
			}

			it = mostAcceleratedParticlesNumbers[t].begin();
			for (int i = 0; i < mostFastParticlesNumber; ++i) {
				std::pair<int, double> pair = *it;
				tempAcceleratedNumbers[i] = pair.first;
				tempAcceleratedMomentum[i] = pair.second;
				++it;
			}

			MPI_Bcast(tempAcceleratedNumbers, mostFastParticlesNumber, MPI_INT, 0, cartComm);
			MPI_Bcast(tempAcceleratedMomentum, mostFastParticlesNumber, MPI_DOUBLE, 0, cartComm);

			mostAcceleratedParticlesNumbers[t].clear();

			for (int i = 0; i < mostFastParticlesNumber; ++i) {
				mostAcceleratedParticlesNumbers[t].push_back(std::pair<int, double>(tempAcceleratedNumbers[i], tempAcceleratedMomentum[i]));
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
	double delta = min2(deltaX, min2(deltaY, deltaZ));
	//double delta = deltaX;
	deltaT = timeEpsilonKourant * delta / speed_of_light_normalized;
	if ((rank == 0) && (verbosity > 1)) printLog("dx/c\n");
	if (particles.size() > 0) {
		double derEmax = 0;
		double B = B0.norm();
		double E = E0.norm();
		if(solverType == BUNEMAN){
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
						if (fabs(bunemanEx[i + 1][j][k] - bunemanEx[i][j][k])/ deltaX > derEmax) {
							derEmax = fabs(bunemanEx[i + 1][j][k] - bunemanEx[i][j][k]) / deltaX;
						}
						if (fabs(bunemanEy[i][j+1][k] - bunemanEy[i][j][k]) / deltaY > derEmax) {
							derEmax = fabs(bunemanEy[i][j+1][k] - bunemanEy[i][j][k]) / deltaY;
						}
						if (fabs(bunemanEz[i][j][k+1] - bunemanEz[i][j][k]) / deltaZ > derEmax) {
							derEmax = fabs(bunemanEz[i][j][k+1] - bunemanEz[i][j][k]) / deltaZ;
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
						if ((Efield[i][j+1][k] - Efield[i][j][k]).norm() / deltaY > derEmax) {
							derEmax = (Efield[i][j+1][k] - Efield[i][j][k]).norm() / deltaY;
						}
						if ((Efield[i][j][k+1] - Efield[i][j][k]).norm() / deltaZ > derEmax) {
							derEmax = (Efield[i][j+1][k] - Efield[i][j][k]).norm() / deltaZ;
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



		omegaPlasmaProton = sqrt(
			4 * pi * types[1].concentration * types[1].charge * types[1].charge / types[1].mass);
		omegaPlasmaElectron = sqrt(
			4 * pi * types[0].concentration * types[0].charge * types[0].charge / types[0].mass);
		double omegaCapture = sqrt(electron_charge_normalized * derEmax / massElectron);


		omegaGyroProton = electron_charge_normalized * B / (massProton * speed_of_light_normalized);
		omegaGyroElectron = electron_charge_normalized * B / (massElectron * speed_of_light_normalized);

		double omega2 = 0;
		for (int i = 0; i < typesNumber; ++i) {
			omega2 += 4 * pi * types[i].concentration * types[i].charge * types[i].charge / types[i].mass;
		}

		double gamma = 1.0/sqrt(1 - V0.scalarMult(V0)/speed_of_light_normalized_sqr);
		double omega = sqrt(omega2/gamma);

		if(timeEpsilonPlasma/omega < deltaT){
			if(rank == 0){
				printf("plasma frequency time limitation\n");
			}
		}

		deltaT = min2(deltaT, timeEpsilonPlasma / omega);


		double thermalMomentum = sqrt(massElectron * kBoltzman_normalized * temperature) + massElectron * V0.norm();

		

		//if ((rank == 0) && writeLog)) printLog("evaluzted omega\n");

		if (B > 0) {
			if(timeEpsilonFields * gamma * massElectron * speed_of_light_normalized / (electron_charge_normalized * B) < deltaT){
				if(rank == 0){
					printf("gyro frequency time limitation\n");
				}
			}
			deltaT = min2(deltaT,
			              timeEpsilonFields * gamma * massElectron * speed_of_light_normalized / (electron_charge_normalized * B));
		}
		if (E > 0) {
			/*if(timeEpsilonFields * massElectron * speed_of_light_normalized / (electron_charge_normalized * E) < deltaT){
				if(rank == 0){
					printf("electric acceleration time limitation\n");
				}
			}
			deltaT = min2(deltaT,
			              timeEpsilonFields * massElectron * speed_of_light_normalized / (electron_charge_normalized * E));*/
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
	MPI_Barrier(cartComm);
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
		fprintf(errorLogFile, "particle.v/c = %15.10g\n",
		        (particle.getVelocity(speed_of_light_normalized).norm() / speed_of_light_normalized));
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (particle.coordinates.x > xgrid[xnumberAdded - additionalBinNumber]) {
		printf("particle.x > xsize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.x = %15.10g > %15.10g\n", particle.coordinates.x, xgrid[xnumberAdded - additionalBinNumber]);
		printf("particle.x = %15.10g > %15.10g\n", particle.coordinates.x, xgrid[xnumberAdded - additionalBinNumber]);
		printf("2*xsizeGeneral = %15.10g \n", 2 * xsizeGeneral);
		printf("xgrid[0] = %15.10g \n", xgrid[0]);
		printf("xsize = %15.10g \n", xsize);
		fprintf(errorLogFile, "2*xsizeGeneral = %15.10g \n", 2 * xsizeGeneral);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x,
		        particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		fprintf(errorLogFile, "particle.v/c = %15.10g\n",
		        (particle.getVelocity(speed_of_light_normalized).norm() / speed_of_light_normalized));
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
		fprintf(errorLogFile, "particle.v/c = %15.10g\n",
		        (particle.getVelocity(speed_of_light_normalized).norm() / speed_of_light_normalized));
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (particle.coordinates.y > ygrid[ynumberAdded - additionalBinNumber]) {
		printf("particle.y > ysize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.y = %15.10g > %15.10g\n", particle.coordinates.y, ygrid[ynumberAdded - additionalBinNumber]);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x,
		        particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		fprintf(errorLogFile, "particle.v/c = %15.10g\n",
		        (particle.getVelocity(speed_of_light_normalized).norm() / speed_of_light_normalized));
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
		fprintf(errorLogFile, "particle.v/c = %15.10g\n",
		        (particle.getVelocity(speed_of_light_normalized).norm() / speed_of_light_normalized));
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (particle.coordinates.z > zgrid[znumberAdded - additionalBinNumber]) {
		if (rank == 0) printf("particle.z > zsize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.z = %15.10g > %15.10g\n", particle.coordinates.z, zgrid[znumberAdded - additionalBinNumber]);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x,
		        particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		fprintf(errorLogFile, "particle.v/c = %15.10g\n",
		        (particle.getVelocity(speed_of_light_normalized).norm() / speed_of_light_normalized));
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
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
	if (cartCoord[0] == 0 && boundaryConditionType != PERIODIC) {
		minI = 0;
	}
	int maxI = xnumberAdded - 1 - additionalBinNumber;
	if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionType != PERIODIC) {
		maxI = xnumberAdded - 1;
	}
	int minJ = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - 1 - additionalBinNumber;
	int minK = 1 + additionalBinNumber;
	int maxK = znumberAdded - 1 - additionalBinNumber;

	globalMomentum = Vector3d(0, 0, 0);
	electromagneticMomentum = Vector3d(0, 0, 0);
	particleMomentum = Vector3d(0, 0, 0);
	if (currentIteration % writeGeneralParameter == 0) {
		//particlesNumber = particles.size();
		for (int i = minI; i < maxI; ++i) {
			for (int j = minJ; j < maxJ; ++j) {
				for (int k = minK; k < maxK; ++k) {
					Vector3d E;
					if(solverType == BUNEMAN){
						E = getBunemanElectricField(i, j, k);
					} else {
						E = Efield[i][j][k];
					}
					electricFieldEnergy += E.scalarMult(E) * volumeE() / (8 * pi);
					electricFieldEnergyX += E.x*E.x* volumeE() / (8 * pi);
					electricFieldEnergyY += E.y*E.y* volumeE() / (8 * pi);
					electricFieldEnergyZ += E.z*E.z* volumeE() / (8 * pi);
					//if (boundaryConditionType == PERIODIC || rank > 0 || i > 1) {
						//electricFieldEnergy -= E0.scalarMult(E0) * volumeE() / (8 * pi);
					//}
				}
			}
		}

		for (int i = minI; i < maxI; ++i) {
			for (int j = minJ; j < maxJ; ++j) {
				for (int k = minK; k < maxK; ++k) {
					Vector3d B;
					if(solverType == BUNEMAN){
						B = getBunemanMagneticField(i, j, k);
					} else {
						B = Bfield[i][j][k];
					}
					magneticFieldEnergy += (B.scalarMult(B)) * volumeB() / (8 * pi);
					magneticFieldEnergyX += B.x*B.x* volumeB() / (8 * pi);
					magneticFieldEnergyY += B.y*B.y* volumeB() / (8 * pi);
					magneticFieldEnergyZ += B.z*B.z* volumeB() / (8 * pi);
					//if (boundaryConditionType == PERIODIC || rank > 0 || i > 0) {
						//magneticFieldEnergy -= B0.scalarMult(B0) * volumeB() / (8 * pi);
					//}
				}
			}
		}

		for (int i = minI; i < maxI; ++i) {
			for (int j = minJ; j < maxJ; ++j) {
				for (int k = minK; k < maxK; ++k) {
					Vector3d E;
					if(solverType == BUNEMAN){
						E = getBunemanElectricField(i, j, k);
					} else {
						E = ((Efield[i][j][k] + Efield[i][j][k + 1] + Efield[i][j + 1][k] + Efield[i][j + 1][k + 1] + Efield[i + 1][j][k] + Efield[i + 1][j][k + 1] + Efield[i + 1][j + 1][k] + Efield[i + 1][j + 1][k + 1]) / 8);
					}
					Vector3d B;
					if(solverType == BUNEMAN){
						B = getBunemanMagneticField(i, j, k);
					} else {
						B = Bfield[i][j][k];
					}

					globalMomentum += (E.vectorMult(B) / (4 * pi * speed_of_light_normalized)) * volumeB();
					electromagneticMomentum += (E.vectorMult(B) / (4 * pi * speed_of_light_normalized)) * volumeB();
				}
			}
		}

		for (int pcount = 0; pcount < particles.size(); ++pcount) {
			double particleE = particles[pcount]->energy(speed_of_light_normalized);
			if (particleE < 0) {
				printf("particle energy < o\n");
				particleE = particles[pcount]->energy(speed_of_light_normalized);
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

	if (currentIteration > 0) {
		if (boundaryConditionType != PERIODIC) {

			for (int i = 0; i < escapedParticlesLeft.size(); ++i) {
				Particle* particle = escapedParticlesLeft[i];
				theoreticalEnergy -= particle->energy(speed_of_light_normalized) * particle->weight *
					sqr(scaleFactor / plasma_period);
				theoreticalMomentum -= particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
			}
			for (int i = 0; i < escapedParticlesRight.size(); ++i) {
				Particle* particle = escapedParticlesRight[i];
				theoreticalEnergy -= particle->energy(speed_of_light_normalized) * particle->weight *
					sqr(scaleFactor / plasma_period);
				theoreticalMomentum -= particle->getMomentum() * particle->weight * scaleFactor / plasma_period;
			}

			for (int j = minJ; j < maxJ; ++j) {
				for (int k = minK; k < maxK; ++k) {
					Vector3d Eleft;
					Vector3d Eright;
					Vector3d Bleft;
					Vector3d Bright;
					if(solverType == BUNEMAN){
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
					theoreticalEnergy -= (Eright.vectorMult(Bright).x * deltaT *speed_of_light_normalized * deltaZ * deltaY *sqr(scaleFactor / plasma_period)) / (4 * pi);
					theoreticalEnergy += (Eleft.vectorMult(Bleft).x * deltaT * speed_of_light_normalized * deltaZ * deltaY * sqr(scaleFactor / plasma_period)) / (4 * pi);

					theoreticalMomentum -=
						(Eright.vectorMult(Bright)) * deltaT * deltaZ * deltaY *
						scaleFactor / plasma_period / (4 * pi);
					theoreticalMomentum +=
						(Eleft.vectorMult(Bleft)) * deltaT * deltaZ * deltaY * scaleFactor /
						plasma_period / (4 * pi);
				}
			}
		}
	}

	if (currentIteration % writeGeneralParameter == 0) {
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

void Simulation:: updateParameters() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if(solverType == BUNEMAN){
		maxBfield = Vector3d(bunemanBx[0][0][0], bunemanBy[0][0][0], bunemanBz[0][0][0]) - B0;
		maxEfield = Vector3d(bunemanEx[0][0][0], bunemanEy[0][0][0], bunemanEz[0][0][0]);

		int minI = 1 + additionalBinNumber;
		if (cartCoord[0] == 0 && boundaryConditionType != PERIODIC) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1 - additionalBinNumber;
		if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionType != PERIODIC) {
			maxI = xnumberAdded - 1;
		}
		int minJ = 1 + additionalBinNumber;
		int maxJ = ynumberAdded - 1 - additionalBinNumber;
		int minK = 1 + additionalBinNumber;
		int maxK = znumberAdded - 1 - additionalBinNumber;

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
		if (cartCoord[0] == 0 && boundaryConditionType != PERIODIC) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1 - additionalBinNumber;
		if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionType != PERIODIC) {
			maxI = xnumberAdded - 1;
		}
		int minJ = 1 + additionalBinNumber;
		int maxJ = ynumberAdded - 1 - additionalBinNumber;
		int minK = 1 + additionalBinNumber;
		int maxK = znumberAdded - 1 - additionalBinNumber;

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
	MPI_Barrier(cartComm);

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
		parallelV2[t] = parallelV2[t] + sqr(particle->velocityX(speed_of_light_normalized)) * particle->weight;
		normalV2[t] = normalV2[t] + (sqr(particle->velocityY(speed_of_light_normalized))
			+ sqr(particle->velocityZ(speed_of_light_normalized))) * particle->weight;
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
			types[t].parallelTemperatureEvaluated = types[t].mass * parallelV2[t] / (kBoltzman_normalized * types[t].generalWeight);
			types[t].normalTemperatureEvaluated = 0.5 * types[t].mass * normalV2[t] / (kBoltzman_normalized * types[t].generalWeight);
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
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating anisotropy time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::resetNewTempFields() {
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

	MPI_Barrier(cartComm);
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

