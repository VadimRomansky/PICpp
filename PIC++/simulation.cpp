#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>
#include <time.h>
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "util.h"
#include "output.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "particle.h"
#include "random.h"
#include "simulation.h"
#include "mpi_util.h"


void Simulation::simulate() {
	MPI_Barrier(MPI_COMM_WORLD);
	if (newlyStarted) {
		//if(rank == 0) printf("create arrays\n");
		createArrays();
		initialize();
		createFiles();
		//initializeTwoStream();
		//initializeExternalFluxInstability();
		//initializeAlfvenWaveX(1, 0.01);
		//initializeRotatedAlfvenWave(1, 0.01);
		//initializeAnisotropic();
		//initializeAnisotropicSilicon();
		//initializeWeibel();
		//initializeRingWeibel();
		initializeFluxFromRight();
		//initializeSimpleElectroMagneticWave();
		//initializeHomogenouseFlow();
		//initializeRotatedSimpleElectroMagneticWave(1);
		//initializeLangmuirWave();
		//createParticles();
	}

	outputGeneralInitialParameters((outputDir + "initialParameters.dat").c_str(), (outputDir + "initialParametersWithText.dat").c_str(), this);


	MPI_Barrier(MPI_COMM_WORLD);

	if ((rank == 0) && (verbosity > 1)) printf("initial exchanging fields\n");
	exchangeEfield();
	exchangeGeneralBfield(Bfield, additionalBfieldLeft, additionalBfieldRight);
	exchangeGeneralBfield(newBfield, additionalNewBfieldLeft, additionalNewBfieldRight);

	//printf("E.x = %g E.y = %g E.z = %g\n", Efield[0][0][0].x, Efield[0][0][0].y, Efield[0][0][0].z);
	//printf("tempE.x = %g tempE.y = %g tempE.z = %g\n", tempEfield[0][0][0].x, tempEfield[0][0][0].y, tempEfield[0][0][0].z);
	//outputEverythingFile = fopen("./output/outputEverythingFile.dat", "w");
	if ((rank == 0) && (verbosity > 1)) printf("initial collecting particles\n");
	updateParticleCorrelationMaps();
	//fclose(outputEverythingFile);
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
		updateShockWaveX();
	}

	evaluateExplicitDerivative();
	if (currentIteration % divergenceCleanUpParameter == 0) {
		cleanupDivergence();
	}
	updateFields();
	exchangeGeneralBfield(Bfield, additionalBfieldLeft, additionalBfieldRight);
	//updateEnergy();

	//printf("E.x = %g E.y = %g E.z = %g\n", Efield[0][0][0].x, Efield[0][0][0].y, Efield[0][0][0].z);
	//printf("tempE.x = %g tempE.y = %g tempE.z = %g\n", tempEfield[0][0][0].x, tempEfield[0][0][0].y, tempEfield[0][0][0].z);

	for (int i = 0; i < typesNumber; ++i) {
		types[i].injectionLength = types[i].particesDeltaX - 0.0001 * deltaX;
	}

	while (time * plasma_period < maxTime && currentIteration < maxIteration) {
		if ((rank == 0) && (verbosity > 0)) {
			FILE* logFile = fopen((outputDir + "log.dat").c_str(), "a");
			fprintf(logFile, "start iteration number = %d time = %15.10g\n", currentIteration, time);
			fflush(logFile);
			fclose(logFile);
		}
		//if ((rank == 0)) printLog("start iteration\n");
		//if ((rank == 0)) printf(" dt/plasma_period = %15.10g\n", deltaT);
		if ((rank == 0)) printf("start iteration number = %d time = %15.10g\n", currentIteration, time);
		updateEnergy();
		if (currentIteration % writeParameter == 0) {
			output();
		}

		if (currentIteration % writeTrajectoryNumber == 0) {
			outputTrajectories();
		}

		double iterationTime = 0;
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			iterationTime = clock();
		}

		updateDeltaT();

		evaluateParticlesRotationTensor();

		updateElectroMagneticParameters();
		//smoothChargeDensityHat();

		evaluateElectricField();

		evaluateMagneticField();
		double procTime = 0;
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock();
		}
		exchangeGeneralBfield(newBfield, additionalNewBfieldLeft, additionalNewBfieldRight);
		exchangeGeneralEfield(tempEfield, additionalTempEfieldLeft, additionalTempEfieldRight);
		exchangeGeneralEfield(newEfield, additionalNewEfieldLeft, additionalNewEfieldRight);
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("exchanging field time = %g sec\n", procTime / CLOCKS_PER_SEC);
		}
		//smoothTempEfield();
		//smoothBfield();


		//MPI_Barrier(MPI_COMM_WORLD);
		/*if(rank == 0){
		    printf("Efield[-1].y = %15.10g\n", Efield[0][0][0].y);
		    printf("Efield[0].y = %15.10g\n", Efield[1][0][0].y);
		    printf("Efield[n/2 - 1].y left = %15.10g\n", Efield[xnumber - 1][0][0].y);
		    printf("Efield[n/2].y left = %15.10g\n", Efield[xnumber][0][0].y);

		} else {
		    printf("Efield[n/2 - 1].y right = %15.10g\n", Efield[0][0][0].y);
		    printf("Efield[n/2].y right = %15.10g\n", Efield[1][0][0].y);
		    printf("Efield[n-1].y right = %15.10g\n", Efield[xnumber - 1][0][0].y);
		    printf("Efield[n].y right = %15.10g\n", Efield[xnumber][0][0].y);
		}
		MPI_Barrier(MPI_COMM_WORLD);*/


		//printf("E.x = %g E.y = %g E.z = %g\n", Efield[0][0][0].x, Efield[0][0][0].y, Efield[0][0][0].z);
		//printf("tempE.x = %g tempE.y = %g tempE.z = %g\n", tempEfield[0][0][0].x, tempEfield[0][0][0].y, tempEfield[0][0][0].z);

		moveParticles();

		if ((rank == 0) && (verbosity > 1)) printLog("erasing escaped particles\n");
		if ((rank == 0) && (verbosity > 1)) printf("erasing escaped particles\n");

		eraseEscapedPaticles();

		//MPI_Barrier(MPI_COMM_WORLD);
		if ((rank == 0) && (verbosity > 1)) printLog("start exchange particles\n");
		if ((rank == 0) && (verbosity > 1)) printf("start exchange particles\n");
		MPI_Barrier(MPI_COMM_WORLD);

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
			if (rank == nprocs - 1) {
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
			if ((rank == 0) && (verbosity > 1)) printf("exchanging particles number\n");
			MPI_Barrier(MPI_COMM_WORLD);
			int tempParticleNumber[1];
			tempParticleNumber[0] = particlesNumber;
			MPI_Bcast(tempParticleNumber, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD);
			particlesNumber = tempParticleNumber[0];
		}
		if (preserveChargeGlobal && (boundaryConditionType != PERIODIC)) {
			addToPreserveChargeGlobal();
		}
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("injecting and preserving time = %g sec\n", procTime / CLOCKS_PER_SEC);
		}

		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock();
		}
		updateParticleCorrelationMaps();
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("updating correlation maps = %g sec\n", procTime / CLOCKS_PER_SEC);
		}

		updateDensityParameters();
		exchangeGeneralEfield(newEfield, additionalNewEfieldLeft, additionalNewEfieldRight);
		//smoothChargeDensity();
		smoothNewEfield();
		if (currentIteration % divergenceCleanUpParameter == 0) {
			cleanupDivergence();
		}
		//smoothNewBfield();

		updateFields();

		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock();
		}
		exchangeGeneralEfield(newEfield, additionalNewEfieldLeft, additionalNewEfieldRight);
		exchangeGeneralBfield(Bfield, additionalBfieldLeft, additionalBfieldRight);
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("exchanging fields time = %g sec\n", procTime / CLOCKS_PER_SEC);
		}

		//smoothNewEfield();

		//updateEnergy();

		if ((currentIteration + 1) % writeParameter == 0) {
			updateParameters();
			updateAnisotropy();
		}

		if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
			updateShockWaveX();
		}

		removeEscapedParticles();

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
	if ((rank == 0) && (verbosity > 1)) printf("outputing trajectory 1\n");
	outputTrajectoryByNumber((outputDir + "trajectory_proton.dat").c_str(), protonNumber, this);
	if ((rank == 0) && (verbosity > 1)) printf("outputing trajectory 2\n");
	outputTrajectoryByNumber((outputDir + "trajectory_electron.dat").c_str(), electronNumber, this);
	if ((rank == 0) && (verbosity > 1)) printf("outputing trajectory 3\n");
	outputTrajectoryByNumber((outputDir + "trajectory_proton_1.dat").c_str(), protonNumber1, this);
	if ((rank == 0) && (verbosity > 1)) printf("outputing trajectory 4\n");
	outputTrajectoryByNumber((outputDir + "trajectory_electron_1.dat").c_str(), electronNumber1, this);
	if ((rank == 0) && (verbosity > 1)) printf("outputing trajectory 5\n");
	outputTrajectoryByNumber((outputDir + "trajectory_proton_2.dat").c_str(), protonNumber2, this);
	if ((rank == 0) && (verbosity > 1)) printf("outputing trajectory 6\n");
	outputTrajectoryByNumber((outputDir + "trajectory_electron_2.dat").c_str(), electronNumber2, this);

	outputTrajectoryByNumber((outputDir + "trajectory_proton_3.dat").c_str(), protonNumber3, this);
	outputTrajectoryByNumber((outputDir + "trajectory_proton_4.dat").c_str(), protonNumber4, this);
	outputTrajectoryByNumber((outputDir + "trajectory_proton_5.dat").c_str(), protonNumber5, this);
	outputTrajectoryByNumber((outputDir + "trajectory_proton_6.dat").c_str(), protonNumber6, this);
	outputTrajectoryByNumber((outputDir + "trajectory_proton_7.dat").c_str(), protonNumber7, this);
	outputTrajectoryByNumber((outputDir + "trajectory_proton_8.dat").c_str(), protonNumber8, this);
	outputTrajectoryByNumber((outputDir + "trajectory_proton_9.dat").c_str(), protonNumber9, this);

	outputTrajectoryByNumber((outputDir + "trajectory_electron_3.dat").c_str(), electronNumber3, this);
	outputTrajectoryByNumber((outputDir + "trajectory_electron_4.dat").c_str(), electronNumber4, this);
	outputTrajectoryByNumber((outputDir + "trajectory_electron_5.dat").c_str(), electronNumber5, this);
	outputTrajectoryByNumber((outputDir + "trajectory_electron_6.dat").c_str(), electronNumber6, this);
	outputTrajectoryByNumber((outputDir + "trajectory_electron_7.dat").c_str(), electronNumber7, this);
	outputTrajectoryByNumber((outputDir + "trajectory_electron_8.dat").c_str(), electronNumber8, this);
	outputTrajectoryByNumber((outputDir + "trajectory_electron_9.dat").c_str(), electronNumber9, this);
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

	Vector3d shockWaveV = V0/3;

	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution protons shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_protons_sw.dat").c_str(), particles, shockWaveV, speed_of_light_normalized, PROTON, scaleFactor,
	                   plasma_period, verbosity);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution electrons shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_electrons_sw.dat").c_str(), particles, shockWaveV, speed_of_light_normalized, ELECTRON, scaleFactor,
	                   plasma_period, verbosity);
	if ((rank == 0) && (verbosity > 1)) printf("outputing distribution alphas shock wave\n");
	outputDistributionShiftedSystem((outputDir + "distribution_alphas_sw.dat").c_str(), particles, shockWaveV, speed_of_light_normalized, ALPHA, scaleFactor,
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
	outputFields((outputDir + "Efield.dat").c_str(), (outputDir + "Bfield.dat").c_str(), Efield, Bfield, xnumber,
	             ynumber, znumber, plasma_period, scaleFactor);

	if ((rank == 0) && (verbosity > 1)) printf("outputing grid\n");
	outputGrid((outputDir + "Xfile.dat").c_str(), xgrid, xnumber, scaleFactor);

	if (rank == 0) outputGridSimple((outputDir + "Yfile.dat").c_str(), ygrid, ynumber, scaleFactor);

	if (rank == 0) outputGridSimple((outputDir + "Zfile.dat").c_str(), zgrid, znumber, scaleFactor);

	//todo!!! for debug only!
	switch (rank) {
	case 0:
		outputGridSimple((outputDir + "Xfile0.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield0.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	case 1:
		outputGridSimple((outputDir + "Xfile1.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield1.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	case 2:
		outputGridSimple((outputDir + "Xfile2.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield2.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	case 3:
		outputGridSimple((outputDir + "Xfile3.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield3.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	case 4:
		outputGridSimple((outputDir + "Xfile4.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield4.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	case 5:
		outputGridSimple((outputDir + "Xfile5.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield5.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	case 6:
		outputGridSimple((outputDir + "Xfile6.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield6.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	case 7:
		outputGridSimple((outputDir + "Xfile7.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield7.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	case 8:
		outputGridSimple((outputDir + "Xfile8.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield8.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	case 9:
		outputGridSimple((outputDir + "Xfile9.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield9.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	case 10:
		outputGridSimple((outputDir + "Xfile10.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield10.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	case 11:
		outputGridSimple((outputDir + "Xfile11.dat").c_str(), xgrid, xnumber, scaleFactor);
		outputVectorNodeArraySimple((outputDir + "Efield11.dat").c_str(), Efield, xnumber, ynumber, znumber, 1.0 / (plasma_period * sqrt(scaleFactor)));
		break;
	}

	if ((rank == 0) && (verbosity > 1)) printf("outputing concentrations\n");
	outputConcentrations((outputDir + "concentrations.dat").c_str(), particleConcentrations, chargeDensity,
	                     chargeDensityHat, xnumber, ynumber, znumber, typesNumber, plasma_period, scaleFactor);



	if ((rank == 0) && (verbosity > 1)) printf("outputing velocity\n");
	outputVelocity((outputDir + "velocity.dat").c_str(),
	               particleBulkVelocities, types, xnumber, ynumber, znumber, typesNumber, plasma_period, scaleFactor);

	/*if ((rank == 0) && (verbosity > 1)) printf("outputing flux\n");
	outputFlux((outputDir + "flux.dat").c_str(), electricFlux, externalElectricFlux, xnumber + 1, ynumber + 1,
	           znumber + 1, plasma_period, scaleFactor);*/

	/*if ((rank == 0) && (verbosity > 1)) printf("outputing divergence\n");
	outputDivergenceError((outputDir + "divergence_error.dat").c_str(), this, plasma_period, scaleFactor);*/

	double rotBscale = 1.0 / (plasma_period * plasma_period * sqrt(scaleFactor));

	/*if ((rank == 0) && (verbosity > 1)) printf("outputing rotB\n");
	outputVectorNodeArray((outputDir + "rotBFile.dat").c_str(), rotB, xnumber + 1, ynumber + 1, znumber + 1, rotBscale);*/

	/*if ((rank == 0) && (verbosity > 1)) printf("outputing Ederivative\n");
	outputVectorNodeArray((outputDir + "EderivativeFile.dat").c_str(), Ederivative, xnumber + 1, ynumber + 1,
	                      znumber + 1, rotBscale);*/

	/*if ((rank == 0) && (verbosity > 1)) printf("outputing rotE\n");
	outputVectorCellArray((outputDir + "rotEFile.dat").c_str(), rotE, xnumber, ynumber, znumber, rotBscale);*/

	/*if (rank == 0) printf("outputing dielectricTensor\n");
	outputMatrixArray((outputDir + "dielectricTensorFile.dat").c_str(), dielectricTensor, xnumber + 1, ynumber + 1,
	                  znumber + 1);*/

	if ((rank == 0) && (verbosity > 1)) printf("outputing particles\n");
	for (int i = 0; i < typesNumber; ++i) {
		outputParticles((outputDir + types[i].typeName + ".dat").c_str(), this, types[i].type);
	}


	/*maxwellMatrixFile = fopen("./output/maxwellMatrixFile.dat", "w");
	outputMaxwellEquationMatrixFull(maxwellMatrixFile, maxwellEquationMatrix, xnumber, ynumber, znumber, maxwellEquationMatrixSize);
	fclose(maxwellMatrixFile);*/

	if (rank == 0) outputGeneral((outputDir + "general.dat").c_str(), this);
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

void Simulation::updateDeltaT() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) printf("updating time step\n");
	if ((rank == 0) && (verbosity > 0)) printLog("updating time step\n");
	//double delta = min2(deltaX, min2(deltaY, deltaZ));
	double delta = deltaX;
	deltaT = timeEpsilonKourant * delta / speed_of_light_normalized;
	if ((rank == 0) && (verbosity > 1)) printLog("dx/c\n");
	if (particles.size() > 0) {
		double derEmax = 0;
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
					if ((Efield[i + 1][j][k] - Efield[i][j][k]).norm() / deltaX > derEmax) {
						derEmax = (Efield[i + 1][j][k] - Efield[i][j][k]).norm() / deltaX;
					}
				}
			}
		}


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
		if ((rank == 0) && (verbosity > 1)) printLog("evaluatd E and minFlux\n");

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


		omegaGyroProton = electron_charge_normalized * B / (massProton * speed_of_light);
		omegaGyroElectron = electron_charge_normalized * B / (massElectron * speed_of_light);

		double omega2 = 0;
		for (int i = 0; i < typesNumber; ++i) {
			omega2 += 4 * pi * types[i].concentration * types[i].charge * types[i].charge / types[i].mass;
		}

		double omega = sqrt(omega2);

		deltaT = min2(deltaT, timeEpsilonPlasma / omega);


		double thermalMomentum = sqrt(massElectron * kBoltzman_normalized * temperature) + massElectron * V0.norm();

		//if ((rank == 0) && writeLog)) printLog("evaluzted omega\n");

		if (B > 0) {
			deltaT = min2(deltaT,
			              timeEpsilonFields * massElectron * speed_of_light_normalized / (electron_charge_normalized * B));
		}
		if (E > 0) {
			deltaT = min2(deltaT,
			              timeEpsilonFields * massElectron * speed_of_light_normalized / (electron_charge_normalized * E));
		}

		if (omegaCapture > 0) {
			deltaT = min2(deltaT, timeEpsilonCapture / omegaCapture);
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

		MPI_Barrier(MPI_COMM_WORLD);
		double temp[1];
		temp[0] = deltaT;
		if (rank != 0) {
			MPI_Send(temp, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD);
		} else {
			for (int i = 1; i < nprocs; ++i) {
				MPI_Status status;
				MPI_Recv(temp, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
				if (temp[0] < deltaT) {
					deltaT = temp[0];
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0) {
			temp[0] = deltaT;
			for (int i = 1; i < nprocs; ++i) {
				MPI_Send(temp, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD);
			}
		} else {
			MPI_Status status;
			MPI_Recv(temp, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
			deltaT = temp[0];
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if ((rank == 0) && (verbosity > 0)) printLog("end updating time step\n");
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating deltaT = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

double Simulation::volumeE(int i, int j, int k) {
	return deltaX * deltaY * deltaZ;
}

double Simulation::volumeB(int i, int j, int k) {
	return deltaX * deltaY * deltaZ;
}

void Simulation::checkParticleInBox(Particle& particle) {
	alertNaNOrInfinity(particle.coordinates.x, "particle.x = NaN in check particle in box\n");
	if (particle.coordinates.x < xgrid[0]) {
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
	if (particle.coordinates.x > xgrid[xnumber]) {
		printf("particle.x > xsize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.x = %15.10g > %15.10g\n", particle.coordinates.x, xgrid[xnumber]);
		printf("particle.x = %15.10g > %15.10g\n", particle.coordinates.x, xgrid[xnumber]);
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

	if (particle.coordinates.y < ygrid[0]) {
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
	if (particle.coordinates.y > ygrid[ynumber]) {
		printf("particle.y > ysize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.y = %15.10g > %15.10g\n", particle.coordinates.y, ygrid[ynumber]);
		fprintf(errorLogFile, "particle.coordinates = %15.10g %15.10g %15.10g\n", particle.coordinates.x,
		        particle.coordinates.y, particle.coordinates.z);
		fprintf(errorLogFile, "particle.n = %d\n", particle.number);
		fprintf(errorLogFile, "particle.v/c = %15.10g\n",
		        (particle.getVelocity(speed_of_light_normalized).norm() / speed_of_light_normalized));
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	if (particle.coordinates.z < zgrid[0]) {
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
	if (particle.coordinates.z > zgrid[znumber]) {
		if (rank == 0) printf("particle.z > zsize\n");
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "particle.z = %15.10g > %15.10g\n", particle.coordinates.z, zgrid[znumber]);
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

void Simulation::updateElectroMagneticParameters() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	if ((rank == 0) && (verbosity > 0)) printf("updating flux, density snd dielectric tensor\n");
	if ((rank == 0) && (verbosity > 0)) printLog("updating flux, density and dielectric tensor\n");
	//check particle only in one boundary box
	if (debugMode) {
		//checkParticlesInB
	}
	//collectParticlesIntoBins();
	//fopen("./output/outputEverythingFile.dat","a");
	int particlePartsCount = 0;
	for (int i = 0; i < xnumber + 2; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				electricFlux[i][j][k] = Vector3d(0, 0, 0);
				electricFluxMinus[i][j][k] = Vector3d(0, 0, 0);
				dielectricTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				divPressureTensor[i][j][k] = Vector3d(0, 0, 0);
			}
		}
	}
	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				additionalElectricFluxLeft[i][j][k] = Vector3d(0, 0, 0);
				additionalElectricFluxMinusLeft[i][j][k] = Vector3d(0, 0, 0);
				additionalElectricFluxRight[i][j][k] = Vector3d(0, 0, 0);
				additionalElectricFluxMinusRight[i][j][k] = Vector3d(0, 0, 0);
				additionalDielectricTensorLeft[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				additionalDielectricTensorRight[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				additionalDivPressureTensorLeft[i][j][k] = Vector3d(0, 0, 0);
				additionalDivPressureTensorRight[i][j][k] = Vector3d(0, 0, 0);
			}
		}
	}

	//todo charges of different sign
	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];

		Particle tempParticleShiftX = *particle;
		Particle tempParticleShiftY = *particle;
		Particle tempParticleShiftZ = *particle;
		Vector3d tempRotatedVelocityShiftX;
		Vector3d tempRotatedVelocityShiftY;
		Vector3d tempRotatedVelocityShiftZ;

		Vector3d velocity = particle->getVelocity(speed_of_light_normalized);
		double gamma = particle->gammaFactor(speed_of_light_normalized);
		Vector3d rotatedVelocity = particle->rotationTensor * (velocity * gamma);
		double beta = 0.5 * particle->charge * deltaT / particle->mass;
		Matrix3d tensor = rotatedVelocity.tensorMult(rotatedVelocity);
		double particleOmega = particle->weight * theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge / particle->mass;

		double shiftX = velocity.x * deltaT * 0.1;
		if (fabs(shiftX) < particle->dx * 0.01) {
			shiftX = particle->dx * 0.01;
		}

		double shiftY = velocity.y * deltaT * 0.1;
		if (fabs(shiftY) < particle->dy * 0.01) {
			shiftY = particle->dy * 0.01;
		}

		double shiftZ = velocity.z * deltaT * 0.1;
		if (fabs(shiftZ) < particle->dz * 0.01) {
			shiftZ = particle->dz * 0.01;
		}

		if (solverType == IMPLICIT) {
			tempParticleShiftX.coordinates.x = particle->coordinates.x + shiftX;
			updateCorrelationMaps(tempParticleShiftX);
			Vector3d oldE = correlationEfield(tempParticleShiftX);
			Vector3d oldB = correlationBfield(tempParticleShiftX);
			tempParticleShiftX.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);
			tempRotatedVelocityShiftX = tempParticleShiftX.rotationTensor * (velocity * gamma);
			if (ynumber > 1) {
				tempParticleShiftY.coordinates.y = particle->coordinates.y + shiftY;
				updateCorrelationMaps(tempParticleShiftY);
				oldE = correlationEfield(tempParticleShiftY);
				oldB = correlationBfield(tempParticleShiftY);
				tempParticleShiftY.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);
				tempRotatedVelocityShiftY = tempParticleShiftY.rotationTensor * (velocity * gamma);
			}
			if (znumber > 1) {
				tempParticleShiftZ.coordinates.z = particle->coordinates.z + shiftZ;
				updateCorrelationMaps(tempParticleShiftZ);
				oldE = correlationEfield(tempParticleShiftZ);
				oldB = correlationBfield(tempParticleShiftZ);
				tempParticleShiftZ.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);
				tempRotatedVelocityShiftZ = tempParticleShiftZ.rotationTensor * (velocity * gamma);
			}
		}

		for (int i = 0; i < splineOrder + 2; ++i) {
			for (int j = 0; j < min2(ynumber, splineOrder + 2); ++j) {
				for (int k = 0; k < min2(znumber, splineOrder + 2); ++k) {
					int curI = particle->correlationMapNode.xindex[i];
					int curJ = particle->correlationMapNode.yindex[j];
					int curK = particle->correlationMapNode.zindex[k];
					while (curJ > ynumber) {
						curJ = curJ - ynumber;
					}
					while (curJ < 0) {
						curJ = curJ + ynumber;
					}
					while (curK > znumber) {
						curK = curK - znumber;
					}
					while (curK < 0) {
						curK = curK + znumber;
					}
					double correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.ycorrelation[j] * particle->correlationMapNode.zcorrelation[k] / volumeE(curI, curJ, curK);
					double particleCharge = particle->charge * particle->weight;
					if (curI <= 0 && boundaryConditionType == SUPER_CONDUCTOR_LEFT && rank == 0) {
						Particle tempParticle = *particle;
						updateCorrelationMaps(tempParticle);
						Vector3d oldE = correlationEfield(tempParticleShiftX);
						Vector3d oldB = correlationBfield(tempParticleShiftX);
						tempParticle.reflectMomentumX();
						tempParticleShiftX.reflectMomentumX();
						tempParticleShiftY.reflectMomentumX();
						tempParticleShiftZ.reflectMomentumX();
						velocity.x = -velocity.x;
						tempParticle.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);
						rotatedVelocity = tempParticle.rotationTensor * (velocity * gamma);
						tensor = rotatedVelocity.tensorMult(rotatedVelocity);
						double tempCorrelation = correlationWithEbin(tempParticleShiftX, curI, curJ, curK) / volumeE(curI, curJ, curK);
						if (solverType == IMPLICIT) {
							tempParticleShiftX.coordinates.x = particle->coordinates.x + shiftX;
							updateCorrelationMaps(tempParticleShiftX);
							oldE = correlationEfield(tempParticleShiftX);
							oldB = correlationBfield(tempParticleShiftX);
							tempParticleShiftX.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);
							tempRotatedVelocityShiftX = tempParticleShiftX.rotationTensor * (velocity * gamma);
							if (ynumber > 1) {
								tempParticleShiftY.coordinates.y = particle->coordinates.y + shiftY;
								updateCorrelationMaps(tempParticleShiftY);
								oldE = correlationEfield(tempParticleShiftY);
								oldB = correlationBfield(tempParticleShiftY);
								tempParticleShiftY.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);
								tempRotatedVelocityShiftY = tempParticleShiftY.rotationTensor * (velocity * gamma);
							}
							if (znumber > 1) {
								tempParticleShiftZ.coordinates.z = particle->coordinates.z + shiftZ;
								updateCorrelationMaps(tempParticleShiftZ);
								oldE = correlationEfield(tempParticleShiftZ);
								oldB = correlationBfield(tempParticleShiftZ);
								tempParticleShiftZ.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);
								tempRotatedVelocityShiftZ = tempParticleShiftZ.rotationTensor * (velocity * gamma);
							}
							double tensorDer00X = ((tempRotatedVelocityShiftX.x * tempRotatedVelocityShiftX.x) * tempCorrelation - tensor.matrix[0][0] * correlation) / shiftX;
							double tensorDer01X = ((tempRotatedVelocityShiftX.x * tempRotatedVelocityShiftX.y) * tempCorrelation - tensor.matrix[0][1] * correlation) / shiftX;
							double tensorDer02X = ((tempRotatedVelocityShiftX.x * tempRotatedVelocityShiftX.z) * tempCorrelation - tensor.matrix[0][2] * correlation) / shiftX;

							double tensorDer10Y = 0;
							double tensorDer11Y = 0;
							double tensorDer12Y = 0;

							double tensorDer20Z = 0;
							double tensorDer21Z = 0;
							double tensorDer22Z = 0;


							if (ynumber > 1) {
								tempCorrelation = correlationWithEbin(tempParticleShiftY, curI, j, k) / volumeE(curI, j, k);

								tensorDer10Y = ((tempRotatedVelocityShiftY.y * tempRotatedVelocityShiftY.x) * tempCorrelation - tensor.matrix[1][0] * correlation) / shiftY;
								tensorDer11Y = ((tempRotatedVelocityShiftY.y * tempRotatedVelocityShiftY.y) * tempCorrelation - tensor.matrix[1][1] * correlation) / shiftY;
								tensorDer12Y = ((tempRotatedVelocityShiftY.y * tempRotatedVelocityShiftY.z) * tempCorrelation - tensor.matrix[1][2] * correlation) / shiftY;
							}

							if (znumber > 1) {
								tempCorrelation = correlationWithEbin(tempParticleShiftZ, curI, j, k) / volumeE(curI, j, k);

								tensorDer20Z = ((tempRotatedVelocityShiftZ.z * tempRotatedVelocityShiftZ.x) * tempCorrelation - tensor.matrix[0][0] * correlation) / shiftZ;
								tensorDer21Z = ((tempRotatedVelocityShiftZ.z * tempRotatedVelocityShiftZ.y) * tempCorrelation - tensor.matrix[0][1] * correlation) / shiftZ;
								tensorDer22Z = ((tempRotatedVelocityShiftZ.z * tempRotatedVelocityShiftZ.z) * tempCorrelation - tensor.matrix[0][2] * correlation) / shiftZ;
							}
							if (particleCharge > 0) {
								electricFlux[2 - curI][curJ][curK] += rotatedVelocity * (particleCharge * correlation);
								//electricFlux[2 - curI][curJ][curK] += velocity * (particleCharge * correlation);
							} else {
								electricFluxMinus[2 - curI][curJ][curK] += rotatedVelocity * (particleCharge * correlation);
								//electricFluxMinus[2 - curI][curJ][curK] += velocity * (particleCharge * correlation);						
							}
							dielectricTensor[2 - curI][curJ][curK] = dielectricTensor[curI][curJ][curK] - particle->rotationTensor * (particleOmega * correlation);

							divPressureTensor[2 - curI][curJ][curK].x += (tensorDer00X + tensorDer10Y + tensorDer20Z) * particleCharge;
							divPressureTensor[2 - curI][curJ][curK].y += (tensorDer01X + tensorDer11Y + tensorDer21Z) * particleCharge;
							divPressureTensor[2 - curI][curJ][curK].z += (tensorDer02X + tensorDer12Y + tensorDer22Z) * particleCharge;
						} else {
							if (particleCharge > 0) {
								electricFlux[2 - curI][j][k] += velocity * particleCharge * correlation;
							} else {
								electricFluxMinus[2 - curI][j][k] += velocity * particleCharge * correlation;
							}
						}
					} else {
						double tempCorrelation = correlationWithEbin(tempParticleShiftX, curI, curJ, curK) / volumeE(curI, curJ, curK);

						//double tensorDer00X = 0;
						//double tensorDer01X = 0;
						//double tensorDer02X = 0;

						double tensorDer00X = ((tempRotatedVelocityShiftX.x * tempRotatedVelocityShiftX.x) * tempCorrelation - tensor.matrix[0][0] * correlation) / shiftX;
						double tensorDer01X = ((tempRotatedVelocityShiftX.x * tempRotatedVelocityShiftX.y) * tempCorrelation - tensor.matrix[0][1] * correlation) / shiftX;
						double tensorDer02X = ((tempRotatedVelocityShiftX.x * tempRotatedVelocityShiftX.z) * tempCorrelation - tensor.matrix[0][2] * correlation) / shiftX;

						double tensorDer10Y = 0;
						double tensorDer11Y = 0;
						double tensorDer12Y = 0;

						double tensorDer20Z = 0;
						double tensorDer21Z = 0;
						double tensorDer22Z = 0;


						if (ynumber > 1) {
							tempCorrelation = correlationWithEbin(tempParticleShiftY, curI, j, k) / volumeE(curI, j, k);

							tensorDer10Y = ((tempRotatedVelocityShiftY.y * tempRotatedVelocityShiftY.x) * tempCorrelation - tensor.matrix[1][0] * correlation) / shiftY;
							tensorDer11Y = ((tempRotatedVelocityShiftY.y * tempRotatedVelocityShiftY.y) * tempCorrelation - tensor.matrix[1][1] * correlation) / shiftY;
							tensorDer12Y = ((tempRotatedVelocityShiftY.y * tempRotatedVelocityShiftY.z) * tempCorrelation - tensor.matrix[1][2] * correlation) / shiftY;
						}

						if (znumber > 1) {
							tempCorrelation = correlationWithEbin(tempParticleShiftZ, curI, j, k) / volumeE(curI, j, k);

							tensorDer20Z = ((tempRotatedVelocityShiftZ.z * tempRotatedVelocityShiftZ.x) * tempCorrelation - tensor.matrix[0][0] * correlation) / shiftZ;
							tensorDer21Z = ((tempRotatedVelocityShiftZ.z * tempRotatedVelocityShiftZ.y) * tempCorrelation - tensor.matrix[0][1] * correlation) / shiftZ;
							tensorDer22Z = ((tempRotatedVelocityShiftZ.z * tempRotatedVelocityShiftZ.z) * tempCorrelation - tensor.matrix[0][2] * correlation) / shiftZ;
						}

						if (curI < 0) {
							if (particleCharge > 0) {
								additionalElectricFluxLeft[-curI - 1][curJ][curK] += rotatedVelocity * (particleCharge * correlation);
								//additionalElectricFluxLeft[-curI - 1][curJ][curK] += velocity * (particleCharge * correlation);
							} else {
								//additionalElectricFluxMinusLeft[-curI - 1][curJ][curK] += rotatedVelocity * (particleCharge * correlation);
								additionalElectricFluxMinusLeft[-curI - 1][curJ][curK] += velocity * (particleCharge * correlation);
							}
							additionalDielectricTensorLeft[-curI - 1][curJ][curK] = additionalDielectricTensorLeft[-curI - 1][curJ][curK] - particle->rotationTensor * (particleOmega * correlation);

							additionalDivPressureTensorLeft[-curI - 1][curJ][curK].x += (tensorDer00X + tensorDer10Y + tensorDer20Z) * particleCharge;
							additionalDivPressureTensorLeft[-curI - 1][curJ][curK].y += (tensorDer01X + tensorDer11Y + tensorDer21Z) * particleCharge;
							additionalDivPressureTensorLeft[-curI - 1][curJ][curK].z += (tensorDer02X + tensorDer12Y + tensorDer22Z) * particleCharge;
						} else if (curI > xnumber + 1) {
							if (particleCharge > 0) {
								additionalElectricFluxRight[curI - xnumber - 2][curJ][curK] += rotatedVelocity * (particleCharge * correlation);
								//additionalElectricFluxRight[curI - xnumber - 2][curJ][curK] += velocity * (particleCharge * correlation);
							} else {
								additionalElectricFluxMinusRight[curI - xnumber - 2][curJ][curK] += rotatedVelocity * (particleCharge * correlation);
								//additionalElectricFluxMinusRight[curI - xnumber - 2][curJ][curK] += velocity * (particleCharge * correlation);
							}
							additionalDielectricTensorRight[curI - xnumber - 2][curJ][curK] = additionalDielectricTensorRight[curI - xnumber - 2][curJ][curK] - particle->rotationTensor * (particleOmega * correlation);

							additionalDivPressureTensorRight[curI - xnumber - 2][curJ][curK].x += (tensorDer00X + tensorDer10Y + tensorDer20Z) * particleCharge;
							additionalDivPressureTensorRight[curI - xnumber - 2][curJ][curK].y += (tensorDer01X + tensorDer11Y + tensorDer21Z) * particleCharge;
							additionalDivPressureTensorRight[curI - xnumber - 2][curJ][curK].z += (tensorDer02X + tensorDer12Y + tensorDer22Z) * particleCharge;
						} else {
							if (solverType == IMPLICIT) {
								if (particleCharge > 0) {
									electricFlux[curI][curJ][curK] += rotatedVelocity * (particleCharge * correlation);
									//electricFlux[curI][curJ][curK] += velocity * (particleCharge * correlation);
								} else {
									electricFluxMinus[curI][curJ][curK] += rotatedVelocity * (particleCharge * correlation);
									//electricFluxMinus[curI][curJ][curK] += velocity * (particleCharge * correlation);
								}
								dielectricTensor[curI][curJ][curK] = dielectricTensor[curI][curJ][curK] - particle->rotationTensor * (particleOmega * correlation);

								divPressureTensor[curI][curJ][curK].x += (tensorDer00X + tensorDer10Y + tensorDer20Z) * particleCharge;
								divPressureTensor[curI][curJ][curK].y += (tensorDer01X + tensorDer11Y + tensorDer21Z) * particleCharge;
								divPressureTensor[curI][curJ][curK].z += (tensorDer02X + tensorDer12Y + tensorDer22Z) * particleCharge;
							}
							if (solverType == EXPLICIT) {
								if (particleCharge > 0) {
									electricFlux[curI][curJ][curK] += velocity * particleCharge * correlation;
								} else {
									electricFluxMinus[curI][curJ][curK] += velocity * particleCharge * correlation;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < xnumber + 2; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				electricFlux[i][j][k] += electricFluxMinus[i][j][k];
			}
		}
	}
	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				additionalElectricFluxLeft[i][j][k] += additionalElectricFluxMinusLeft[i][j][k];
				additionalElectricFluxRight[i][j][k] += additionalElectricFluxMinusRight[i][j][k];
			}
		}
	}

	if (rank == 0 && boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				electricFlux[1][j][k] = Vector3d(0, 0, 0);
				electricFluxMinus[1][j][k] = Vector3d(0, 0, 0);
				divPressureTensor[1][j][k] = Vector3d(0, 0, 0);
				dielectricTensor[1][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
			}
		}
	}

	//for periodic y
	for (int i = 0; i < xnumber + 2; ++i) {
		for (int k = 0; k < znumber + 1; ++k) {
			electricFlux[i][0][k] = electricFlux[i][0][k] + electricFlux[i][ynumber][k];
			electricFlux[i][ynumber][k] = electricFlux[i][0][k];

			dielectricTensor[i][0][k] = dielectricTensor[i][0][k] + dielectricTensor[i][ynumber][k];
			dielectricTensor[i][ynumber][k] = dielectricTensor[i][0][k];

			divPressureTensor[i][0][k] = divPressureTensor[i][0][k] + divPressureTensor[i][ynumber][k];
			divPressureTensor[i][ynumber][k] = divPressureTensor[i][0][k];
		}
	}

	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int k = 0; k < znumber + 1; ++k) {
			additionalElectricFluxLeft[i][0][k] = additionalElectricFluxLeft[i][0][k] + additionalElectricFluxLeft[i][ynumber][k];
			additionalElectricFluxLeft[i][ynumber][k] = additionalElectricFluxLeft[i][0][k];
			additionalElectricFluxRight[i][0][k] = additionalElectricFluxRight[i][0][k] + additionalElectricFluxRight[i][ynumber][k];
			additionalElectricFluxRight[i][ynumber][k] = additionalElectricFluxRight[i][0][k];

			additionalDielectricTensorLeft[i][0][k] = additionalDielectricTensorLeft[i][0][k] + additionalDielectricTensorLeft[i][ynumber][k];
			additionalDielectricTensorLeft[i][ynumber][k] = additionalDielectricTensorLeft[i][0][k];
			additionalDielectricTensorRight[i][0][k] = additionalDielectricTensorRight[i][0][k] + additionalDielectricTensorRight[i][ynumber][k];
			additionalDielectricTensorRight[i][ynumber][k] = additionalDielectricTensorRight[i][0][k];

			additionalDivPressureTensorLeft[i][0][k] = additionalDivPressureTensorLeft[i][0][k] + additionalDivPressureTensorLeft[i][ynumber][k];
			additionalDivPressureTensorLeft[i][ynumber][k] = additionalDivPressureTensorLeft[i][0][k];
			additionalDivPressureTensorRight[i][0][k] = additionalDivPressureTensorRight[i][0][k] + additionalDivPressureTensorRight[i][ynumber][k];
			additionalDivPressureTensorRight[i][ynumber][k] = additionalDivPressureTensorRight[i][0][k];
		}
	}

	//for periodic z
	for (int i = 0; i < xnumber + 2; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			electricFlux[i][j][0] = electricFlux[i][j][0] + electricFlux[i][j][znumber];
			electricFlux[i][j][znumber] = electricFlux[i][j][0];

			dielectricTensor[i][j][0] = dielectricTensor[i][j][0] + dielectricTensor[i][j][znumber];
			dielectricTensor[i][j][znumber] = dielectricTensor[i][0][0];

			divPressureTensor[i][j][0] = divPressureTensor[i][j][0] + divPressureTensor[i][j][znumber];
			divPressureTensor[i][j][znumber] = divPressureTensor[i][j][0];
		}
	}

	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			additionalElectricFluxLeft[i][j][0] = additionalElectricFluxLeft[i][j][0] + additionalElectricFluxLeft[i][j][znumber];
			additionalElectricFluxLeft[i][j][znumber] = additionalElectricFluxLeft[i][j][0];
			additionalElectricFluxRight[i][j][0] = additionalElectricFluxRight[i][j][0] + additionalElectricFluxRight[i][j][znumber];
			additionalElectricFluxRight[i][j][znumber] = additionalElectricFluxRight[i][j][0];

			additionalDielectricTensorLeft[i][j][0] = additionalDielectricTensorLeft[i][j][0] + additionalDielectricTensorLeft[i][j][znumber];
			additionalDielectricTensorLeft[i][j][znumber] = additionalDielectricTensorLeft[i][j][0];
			additionalDielectricTensorRight[i][j][0] = additionalDielectricTensorRight[i][j][0] + additionalDielectricTensorRight[i][j][znumber];
			additionalDielectricTensorRight[i][j][znumber] = additionalDielectricTensorRight[i][j][0];

			additionalDivPressureTensorLeft[i][j][0] = additionalDivPressureTensorLeft[i][j][0] + additionalDivPressureTensorLeft[i][j][znumber];
			additionalDivPressureTensorLeft[i][j][znumber] = additionalDivPressureTensorLeft[i][j][0];
			additionalDivPressureTensorRight[i][j][0] = additionalDivPressureTensorRight[i][j][0] + additionalDivPressureTensorRight[i][j][znumber];
			additionalDivPressureTensorRight[i][j][znumber] = additionalDivPressureTensorRight[i][j][0];
		}
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	//fprintf(outputEverythingFile, "electricFlux %d after boundaries = %28.22g %28.22g %28.22g\n", debugPoint, electricFlux[debugPoint][0][0].x, electricFlux[debugPoint][0][0].y, electricFlux[debugPoint][0][0].z);
	if ((rank == 0) && (verbosity > 1)) printf("updating electricDensityHat and pressure tensor\n");

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				chargeDensityHat[i][j][k] = 0;
				chargeDensity[i][j][k] = 0;
				chargeDensityMinus[i][j][k] = 0;
				pressureTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
			}
		}
	}
	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				additionalChargeDensityHatLeft[i][j][k] = 0;
				additionalChargeDensityMinusLeft[i][j][k] = 0;
				additionalChargeDensityHatRight[i][j][k] = 0;
				additionalChargeDensityMinusRight[i][j][k] = 0;
				additionalPressureTensorLeft[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				additionalPressureTensorRight[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
			}
		}
	}

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		for (int i = 0; i < splineOrder + 2; ++i) {
			for (int j = 0; j < min2(ynumber, splineOrder + 2); ++j) {
				for (int k = 0; k < min2(znumber, splineOrder + 2); ++k) {
					int curI = particle->correlationMapCell.xindex[i];
					int curJ = particle->correlationMapCell.yindex[j];
					int curK = particle->correlationMapCell.zindex[k];
					while (curJ >= ynumber) {
						curJ = curJ - ynumber;
					}
					while (curJ < 0) {
						curJ = curJ + ynumber;
					}
					while (curK >= znumber) {
						curK = curK - znumber;
					}
					while (curK < 0) {
						curK = curK + znumber;
					}
					double correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.ycorrelation[j] * particle->correlationMapCell.zcorrelation[k] / volumeB(curI, curJ, curK);
					double particleCharge = particle->charge * particle->weight;


					if (curI <= 0 && boundaryConditionType == SUPER_CONDUCTOR_LEFT && rank == 0) {
						//Particle tempParticle = *particle;
						//tempParticle.reflectMomentumX();
						//double beta = 0.5 * particle->charge * deltaT / particle->mass;
						//double gamma = particle->gammaFactor(speed_of_light_normalized);
						//Vector3d oldE = correlationEfield(particle);
						//Vector3d oldB = correlationBfield(particle);
						//Vector3d tempVelocity = tempParticle.getVelocity(speed_of_light_normalized);
						//tempParticle.rotationTensor = evaluateAlphaRotationTensor(beta, tempVelocity, gamma, oldE, oldB);
						//Vector3d rotatedVelocity = particle->rotationTensor * tempVelocity * gamma;
						if (particleCharge > 0) {
							chargeDensityHat[1 - curI][curJ][curK] += particleCharge * correlation;
						} else {
							chargeDensityMinus[1 - curI][curJ][curK] += particleCharge * correlation;
						}
						//pressureTensor[1 - curI][curJ][curK] += rotatedVelocity.tensorMult(rotatedVelocity) * (particleCharge * correlation);
					} else {
						//Vector3d velocity = particle->getVelocity(speed_of_light_normalized);
						//double gamma = particle->gammaFactor(speed_of_light_normalized);
						//Vector3d rotatedVelocity = particle->rotationTensor * velocity * gamma;
						if (curI < 0) {
							if (particleCharge > 0) {
								additionalChargeDensityHatLeft[-curI - 1][curJ][curK] += particleCharge * correlation;
							} else {
								additionalChargeDensityMinusLeft[-curI - 1][curJ][curK] += particleCharge * correlation;
							}
							//pressureTensor[-curI-1][curJ][curK] += rotatedVelocity.tensorMult(rotatedVelocity) * (particleCharge * correlation);
						} else if (curI >= xnumber + 1) {
							if (particleCharge > 0) {
								additionalChargeDensityHatRight[curI - xnumber - 1][curJ][curK] += particleCharge * correlation;
							} else {
								additionalChargeDensityMinusRight[curI - xnumber - 1][curJ][curK] += particleCharge * correlation;
							}
							//pressureTensor[curI - xnumber - 1][curJ][curK] += rotatedVelocity.tensorMult(rotatedVelocity) * (particleCharge * correlation);
						} else {
							if (particleCharge > 0) {
								chargeDensityHat[curI][curJ][curK] += particleCharge * correlation;
							} else {
								chargeDensityMinus[curI][curJ][curK] += particleCharge * correlation;
							}
							//pressureTensor[curI][curJ][curK] += rotatedVelocity.tensorMult(rotatedVelocity) * (particleCharge * correlation);
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				chargeDensityHat[i][j][k] += chargeDensityMinus[i][j][k];
			}
		}
	}

	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				additionalChargeDensityHatLeft[i][j][k] += additionalChargeDensityMinusLeft[i][j][k];
				additionalChargeDensityHatRight[i][j][k] += additionalChargeDensityMinusRight[i][j][k];
			}
		}
	}

	if (solverType == IMPLICIT) {
		for (int i = 0; i < xnumber + 2; ++i) {
			for (int j = 0; j < ynumber + 1; ++j) {
				for (int k = 0; k < znumber + 1; ++k) {
					//Vector3d divPressureTensorEvaluated = evaluateDivPressureTensor(i, j, k);
					//electricFlux[i] = electricFlux[i] - divPressureTensorEvaluated * eta * deltaT;
					electricFlux[i][j][k] = electricFlux[i][j][k] - divPressureTensor[i][j][k] * eta * deltaT;
					if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].x, "electricFlux[i][j][k].x = NaN");
					if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].y, "electricFlux[i][j][k].y = NaN");
					if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].z, "electricFlux[i][j][k].z = NaN");
				}
			}
		}

		for (int i = 0; i < additionalBinNumber; ++i) {
			for (int j = 0; j < ynumber + 1; ++j) {
				for (int k = 0; k < znumber + 1; ++k) {
					additionalElectricFluxLeft[i][j][k] = additionalElectricFluxLeft[i][j][k] - additionalDivPressureTensorLeft[i][j][k] * eta * deltaT;
					additionalElectricFluxRight[i][j][k] = additionalElectricFluxRight[i][j][k] - additionalDivPressureTensorRight[i][j][k] * eta * deltaT;
				}
			}
		}
	}

	//todo realy?
	if (boundaryConditionType == SUPER_CONDUCTOR_LEFT && rank == 0) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				electricFlux[1][j][k] = Vector3d(0, 0, 0);
				electricFlux[0][j][k] = Vector3d(0, 0, 0);
			}
		}
	}

	//here for evaluating div J

	if (solverType == IMPLICIT) {
		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					double divJ = evaluateDivFlux(i, j, k);

					chargeDensityHat[i][j][k] -= deltaT * theta * divJ;
					//chargeDensityHat[i][j][k] = 0;
				}
			}
		}
		for (int i = 0; i < additionalBinNumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {

					double divJ = evaluateDivFlux(-i - 1, j, k);

					additionalChargeDensityHatLeft[i][j][k] -= deltaT * theta * divJ;

					divJ = evaluateDivFlux(xnumber + 1 + i, j, k);

					additionalChargeDensityHatRight[i][j][k] -= deltaT * theta * divJ;
				}
			}
		}
	}
	//fprintf(outputEverythingFile, "density %d afterFlux = %28.22g\n", chargeDensityHat[debugPoint - 1][0][0]);
	//fprintf(outputEverythingFile, "density %d afterFlux = %28.22g\n", chargeDensityHat[debugPoint][0][0]);
	//fprintf(outputEverythingFile, "density %d afterFlux = %28.22g\n", chargeDensityHat[debugPoint + 1][0][0]);

	/*FILE* tempDensityFile = fopen((outputDir + "tempDensityFile.dat").c_str(), "w");
	for(int i = 0; i < xnumber; ++i){
	    fprintf(tempDensityFile, "%28.22g\n", chargeDensityHat[i][0][0]);
	}
	fclose(tempDensityFile);*/
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating electromagnetic parameters time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	if ((rank == 0) && (verbosity > 0)) printf("start exchange parameters\n");
	if ((rank == 0) && (verbosity > 1)) printf("sum charge density hat\n");
	//MPI_Barrier(MPI_COMM_WORLD);
	sumChargeDensityHat();
	if ((rank == 0) && (verbosity > 1)) printf("sumcell matrix parameters\n");
	//MPI_Barrier(MPI_COMM_WORLD);
	sumCellMatrixParameters();
	if ((rank == 0) && (verbosity > 1)) printf("sum node vector parameters\n");
	//MPI_Barrier(MPI_COMM_WORLD);
	sumNodeVectorParameters();
	if ((rank == 0) && (verbosity > 1)) printf("sum node matrix parameters\n");
	//MPI_Barrier(MPI_COMM_WORLD);
	sumNodeMatrixParameters();
	if ((rank == 0) && (verbosity > 1)) printf("update external flux\n");
	updateExternalFlux();

	for (int i = 0; i < xnumber + 2; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				electricFlux[i][j][k] = electricFlux[i][j][k] + externalElectricFlux[i][j][k];
				//electricFlux[i][j][k].x = 0;
				//electricFlux[i][j][k].y = 0;
				//electricFlux[i][j][k].z = 0;
				if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].x, "electricFlux[i][j][k].x = NaN");
				if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].y, "electricFlux[i][j][k].y = NaN");
				if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].z, "electricFlux[i][j][k].z = NaN");
			}
		}
	}

	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("exchanging electromagnetic parameters time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::sumNodeVectorParameters() {
	if (nprocs > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum node vector parameters\n");
		double* inBufferRight = new double[(2 + additionalBinNumber) * (ynumber + 1) * (znumber + 1) * 3];
		double* outBufferRight = new double[(2 + additionalBinNumber) * (ynumber + 1) * (znumber + 1) * 3];
		double* inBufferLeft = new double[(2 + additionalBinNumber) * (ynumber + 1) * (znumber + 1) * 3];
		double* outBufferLeft = new double[(2 + additionalBinNumber) * (ynumber + 1) * (znumber + 1) * 3];


		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("sending left flux sum node vector parameters rank = %d\n", rank);
		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				sendNodeVectorParametersLeft(electricFlux, additionalElectricFluxLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendNodeVectorParametersLeft(electricFlux, additionalElectricFluxLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if ((verbosity > 2)) printf("receiving right flux sum node vector parameters rank = %d\n", rank);
		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				receiveNodeVectorParametersRight(tempNodeVectorParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveNodeVectorParametersRight(tempNodeVectorParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if ((verbosity > 2)) printf("sending right flux sum node vector parameters rank = %d\n", rank);
		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				sendNodeVectorParametersRight(electricFlux, additionalElectricFluxRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendNodeVectorParametersRight(electricFlux, additionalElectricFluxRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("receiving left flux sum node vector parameters rank = %d\n", rank);
		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				receiveNodeVectorParametersLeft(tempNodeVectorParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveNodeVectorParametersLeft(tempNodeVectorParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}

		sumTempNodeVectorParameters(electricFlux);

		//MPI_Barrier(MPI_COMM_WORLD);
		//printf("send divPressureTensor\n");
		if ((verbosity > 2)) printf("sending left div pressure tensor sum node vector parameters rank = %d\n", rank);
		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				sendNodeVectorParametersLeft(divPressureTensor, additionalDivPressureTensorLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendNodeVectorParametersLeft(divPressureTensor, additionalDivPressureTensorLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if ((verbosity > 2)) printf("receiving right divpressure tensor sum node vector parameters\n");
		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				receiveNodeVectorParametersRight(tempNodeVectorParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveNodeVectorParametersRight(tempNodeVectorParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if ((verbosity > 2)) printf("sending right div pressure tensor sum node vector parameters rank = %d\n", rank);
		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				sendNodeVectorParametersRight(divPressureTensor, additionalDivPressureTensorRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendNodeVectorParametersRight(divPressureTensor, additionalDivPressureTensorRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("receiving left divpressure tensor sum node vector parameters\n");
		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				receiveNodeVectorParametersLeft(tempNodeVectorParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveNodeVectorParametersLeft(tempNodeVectorParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}

		sumTempNodeVectorParameters(divPressureTensor);
		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters rank = %d\n", rank);
		delete[] inBufferLeft;
		delete[] inBufferRight;
		delete[] outBufferLeft;
		delete[] outBufferRight;
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j <= ynumber; ++j) {
				for (int k = 0; k <= znumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						electricFlux[xnumber][j][k][l] += electricFlux[1][j][k][l];
						electricFlux[1][j][k][l] = electricFlux[xnumber][j][k][l];
						divPressureTensor[xnumber][j][k][l] += divPressureTensor[1][j][k][l];
						divPressureTensor[1][j][k][l] = divPressureTensor[xnumber][j][k][l];

						electricFlux[xnumber - 1][j][k][l] += electricFlux[0][j][k][l];
						electricFlux[0][j][k][l] = electricFlux[xnumber - 1][j][k][l];
						divPressureTensor[xnumber - 1][j][k][l] += divPressureTensor[0][j][k][l];
						divPressureTensor[0][j][k][l] = divPressureTensor[xnumber - 1][j][k][l];

						electricFlux[2][j][k][l] += electricFlux[xnumber + 1][j][k][l];
						electricFlux[xnumber + 1][j][k][l] = electricFlux[2][j][k][l];
						divPressureTensor[2][j][k][l] += divPressureTensor[xnumber + 1][j][k][l];
						divPressureTensor[xnumber + 1][j][k][l] = divPressureTensor[2][j][k][l];

						for (int i = 0; i < additionalBinNumber; ++i) {
							electricFlux[2 + i][j][k][l] += additionalElectricFluxRight[i][j][k][l];
							electricFlux[xnumber - 2 - i][j][k][l] += additionalElectricFluxLeft[i][j][k][l];

							divPressureTensor[2 + i][j][k][l] += additionalDivPressureTensorRight[i][j][k][l];
							divPressureTensor[xnumber - 2 - i][j][k][l] += additionalDivPressureTensorLeft[i][j][k][l];
						}
					}
				}
			}
		}
	}
}

void Simulation::sumTempNodeVectorParameters(Vector3d*** array) {
	if (rank > 0 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				array[1][j][k] = array[1][j][k] + tempNodeVectorParameterLeft[0][j][k];
				array[2][j][k] = array[2][j][k] + tempNodeVectorParameterLeft[1][j][k];
				for (int i = 0; i < additionalBinNumber; ++i) {
					array[3 + i][j][k] = array[3 + i][j][k] + tempNodeVectorParameterLeft[2 + i][j][k];
				}
			}
		}
	}

	if (rank < nprocs - 1 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				array[xnumber - 1][j][k] = array[xnumber - 1][j][k] + tempNodeVectorParameterRight[0][j][k];
				array[xnumber][j][k] = array[xnumber][j][k] + tempNodeVectorParameterRight[1][j][k];

				for (int i = 0; i < additionalBinNumber; ++i) {
					array[xnumber - 2 - i][j][k] = array[xnumber - 2 - i][j][k] + tempNodeVectorParameterRight[2 + i][j][k];
				}
			}
		}
	}
}

void Simulation::sumNodeMatrixParameters() {
	if (nprocs > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum node matrix parameters rank = %d\n", rank);
		double* inBufferRight = new double[(2 + additionalBinNumber) * (ynumber + 1) * (znumber + 1) * 9];
		double* outBufferRight = new double[(2 + additionalBinNumber) * (ynumber + 1) * (znumber + 1) * 9];
		double* inBufferLeft = new double[(2 + additionalBinNumber) * (ynumber + 1) * (znumber + 1) * 9];
		double* outBufferLeft = new double[(2 + additionalBinNumber) * (ynumber + 1) * (znumber + 1) * 9];

		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("send left dielectric tensor sum node matrix parameters rank = %d\n", rank);
		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				sendNodeMatrixParametersLeft(dielectricTensor, additionalDielectricTensorLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendNodeMatrixParametersLeft(dielectricTensor, additionalDielectricTensorLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if ((verbosity > 2)) printf("receiving right dielectric tensor sum node matrix parameters rank = %d\n", rank);
		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				receiveNodeMatrixParametersRight(tempNodeMatrixParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveNodeMatrixParametersRight(tempNodeMatrixParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if ((verbosity > 2)) printf("send right dielectric tensor sum node matrix parameters rank = %d\n", rank);
		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				sendNodeMatrixParametersRight(dielectricTensor, additionalDielectricTensorRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendNodeMatrixParametersRight(dielectricTensor, additionalDielectricTensorRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("receiving left dielectric tensor sum node matrix parameters rank = %d\n", rank);
		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				receiveNodeMatrixParametersLeft(tempNodeMatrixParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveNodeMatrixParametersLeft(tempNodeMatrixParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}

		sumTempNodeMatrixParameters(dielectricTensor);
		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("deleting buffer in sum node matrix parameters rank = %d\n", rank);
		delete[] inBufferLeft;
		delete[] inBufferRight;
		delete[] outBufferLeft;
		delete[] outBufferRight;
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j <= ynumber; ++j) {
				for (int k = 0; k <= znumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						for (int m = 0; m < 3; ++m) {
							dielectricTensor[xnumber][j][k].matrix[l][m] += dielectricTensor[1][j][k].matrix[l][m];
							dielectricTensor[1][j][k].matrix[l][m] = dielectricTensor[xnumber][j][k].matrix[l][m];

							dielectricTensor[xnumber - 1][j][k].matrix[l][m] += dielectricTensor[0][j][k].matrix[l][m];
							dielectricTensor[0][j][k].matrix[l][m] = dielectricTensor[xnumber - 1][j][k].matrix[l][m];

							dielectricTensor[xnumber + 1][j][k].matrix[l][m] += dielectricTensor[2][j][k].matrix[l][m];
							dielectricTensor[2][j][k].matrix[l][m] = dielectricTensor[xnumber + 1][j][k].matrix[l][m];

							for (int i = 0; i < additionalBinNumber; ++i) {
								dielectricTensor[2 + i][j][k].matrix[l][m] += additionalDielectricTensorRight[i][j][k].matrix[l][m];
								dielectricTensor[xnumber - 2 - i][j][k].matrix[l][m] += additionalDielectricTensorLeft[i][j][k].matrix[l][m];
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumTempNodeMatrixParameters(Matrix3d*** array) {
	if (rank > 0 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				array[1][j][k] = array[1][j][k] + tempNodeMatrixParameterLeft[0][j][k];
				array[2][j][k] = array[2][j][k] + tempNodeMatrixParameterLeft[1][j][k];
				for (int i = 0; i < additionalBinNumber; ++i) {
					array[3 + i][j][k] = array[3 + i][j][k] + tempNodeMatrixParameterLeft[2 + i][j][k];
				}
			}
		}
	}

	if (rank < nprocs - 1 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				array[xnumber - 1][j][k] = array[xnumber - 1][j][k] + tempNodeMatrixParameterRight[0][j][k];
				array[xnumber][j][k] = array[xnumber][j][k] + tempNodeMatrixParameterRight[1][j][k];

				for (int i = 0; i < additionalBinNumber; ++i) {
					array[xnumber - 2 - i][j][k] = array[xnumber - 2 - i][j][k] + tempNodeMatrixParameterRight[2 + i][j][k];
				}
			}
		}
	}
}

void Simulation::sumTempNodeParameters(double*** array) {
	if (rank > 0 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				array[1][j][k] = array[1][j][k] + tempNodeParameterLeft[0][j][k];
				array[2][j][k] = array[2][j][k] + tempNodeParameterLeft[1][j][k];
				for (int i = 0; i < additionalBinNumber; ++i) {
					array[3 + i][j][k] = array[3 + i][j][k] + tempNodeParameterLeft[2 + i][j][k];
				}
			}
		}
	}

	if (rank < nprocs - 1 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				array[xnumber - 1][j][k] = array[xnumber - 1][j][k] + tempNodeParameterRight[0][j][k];
				array[xnumber][j][k] = array[xnumber][j][k] + tempNodeParameterRight[1][j][k];

				for (int i = 0; i < additionalBinNumber; ++i) {
					array[xnumber - 2 - i][j][k] = array[xnumber - 2 - i][j][k] + tempNodeParameterRight[2 + i][j][k];
				}
			}
		}
	}
}

void Simulation::smoothDensity() {
	for (int j = 0; j < ynumber; ++j) {
		for (int k = 0; k < znumber; ++k) {
			double newLeftDensity = (chargeDensityHat[0][j][k] + chargeDensityHat[1][j][k]) / 2.0;
			double newRightDensity = (chargeDensityHat[xnumber - 1][j][k] + chargeDensityHat[xnumber - 2][j][k]) / 2.0;
			if (boundaryConditionType == PERIODIC) {
				newLeftDensity = (chargeDensityHat[xnumber - 1][j][k] + 2.0 * chargeDensityHat[0][j][k] + chargeDensityHat[1][j][k]) / 4.0;
				newRightDensity = (chargeDensityHat[xnumber - 2][j][k] + 2.0 * chargeDensityHat[xnumber - 1][j][k] + chargeDensityHat[0][j][k]) / 4.0;
			}
			double prevDensity = chargeDensityHat[0][j][k];
			for (int i = 1; i < xnumber - 1; ++i) {
				double tempDensity = chargeDensityHat[i][j][k];

				chargeDensityHat[i][j][k] = (prevDensity + 2.0 * chargeDensityHat[i][j][k] + chargeDensityHat[i + 1][j][k]) / 4.0;

				prevDensity = tempDensity;
			}

			chargeDensityHat[0][j][k] = newLeftDensity;
			chargeDensityHat[xnumber - 1][j][k] = newRightDensity;
		}
	}
}

void Simulation::updateDensityParameters() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) {
		printf("updating densityParameters\n");
	}
	//FILE* debugFile = fopen("./output/particleCorrelations.dat","w");
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int t = 0; t < typesNumber; ++t) {
					particleConcentrations[t][i][j][k] = 0;
					particleBulkVelocities[t][i][j][k] = Vector3d(0, 0, 0);
				}
				chargeDensity[i][j][k] = 0;
				chargeDensityMinus[i][j][k] = 0;
			}
		}
	}

	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int t = 0; t < typesNumber; ++t) {
					additionalParticleConcentrationsLeft[t][i][j][k] = 0;
					additionalParticleConcentrationsRight[t][i][j][k] = 0;
					additionalParticleBulkVelocitiesLeft[t][i][j][k] = Vector3d(0, 0, 0);
					additionalParticleBulkVelocitiesRight[t][i][j][k] = Vector3d(0, 0, 0);
				}
				additionalChargeDensityLeft[i][j][k] = 0;
				additionalChargeDensityMinusLeft[i][j][k] = 0;
				additionalChargeDensityRight[i][j][k] = 0;
				additionalChargeDensityMinusRight[i][j][k] = 0;
			}
		}
	}

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		for (int i = 0; i < splineOrder + 2; ++i) {
			for (int j = 0; j < min2(ynumber, splineOrder + 2); ++j) {
				for (int k = 0; k < min2(znumber, splineOrder + 2); ++k) {
					int curI = particle->correlationMapCell.xindex[i];
					int curJ = particle->correlationMapCell.yindex[j];
					int curK = particle->correlationMapCell.zindex[k];
					while (curJ >= ynumber) {
						curJ = curJ - ynumber;
					}
					while (curJ < 0) {
						curJ = curJ + ynumber;
					}
					while (curK >= znumber) {
						curK = curK - znumber;
					}
					while (curK < 0) {
						curK = curK + znumber;
					}
					double correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.ycorrelation[j] * particle->correlationMapCell.zcorrelation[k] / volumeB(curI, curJ, curK);
					double particleCharge = particle->charge * particle->weight;
					int typeN = getTypeNumber(particle);

					if (curI <= 0 && boundaryConditionType == SUPER_CONDUCTOR_LEFT && rank == 0) {
						Vector3d reflectedMomentum = particle->getMomentum();
						reflectedMomentum.x = -reflectedMomentum.x;
						particleBulkVelocities[typeN][1 - curI][curJ][curK] += reflectedMomentum * particle->weight * correlation;
						if (particleCharge > 0) {
							chargeDensity[1 - curI][curJ][curK] += particleCharge * correlation;
						} else {
							chargeDensityMinus[1 - curI][curJ][curK] += particleCharge * correlation;
						}
						particleConcentrations[typeN][1 - curI][curJ][curK] += particle->weight * correlation;
					} else {
						if (curI < 0) {
							additionalParticleBulkVelocitiesLeft[typeN][-curI - 1][curJ][curK] += particle->getMomentum() * particle->weight * correlation;
							if (particleCharge > 0) {
								additionalChargeDensityLeft[-curI - 1][curJ][curK] += particleCharge * correlation;
							} else {
								additionalChargeDensityMinusLeft[-curI - 1][curJ][curK] += particleCharge * correlation;
							}
							additionalParticleConcentrationsLeft[typeN][- curI - 1][curJ][curK] += particle->weight * correlation;
						} else if (curI >= xnumber + 1) {
							additionalParticleBulkVelocitiesRight[typeN][curI - xnumber - 1][curJ][curK] += particle->getMomentum() * particle->weight * correlation;
							if (particleCharge > 0) {
								additionalChargeDensityRight[curI - xnumber - 1][curJ][curK] += particleCharge * correlation;
							} else {
								additionalChargeDensityMinusRight[curI - xnumber - 1][curJ][curK] += particleCharge * correlation;
							}
							additionalParticleConcentrationsRight[typeN][curI - xnumber - 1][curJ][curK] += particle->weight * correlation;
						} else {
							particleBulkVelocities[typeN][curI][curJ][curK] += particle->getMomentum() * particle->weight * correlation;
							if (particleCharge > 0) {
								chargeDensity[curI][curJ][curK] += particleCharge * correlation;
							} else {
								chargeDensityMinus[curI][curJ][curK] += particleCharge * correlation;
							}
							particleConcentrations[typeN][curI][curJ][curK] += particle->weight * correlation;
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				chargeDensity[i][j][k] += chargeDensityMinus[i][j][k];
			}
		}
	}

	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				additionalChargeDensityLeft[i][j][k] += additionalChargeDensityMinusLeft[i][j][k];
				additionalChargeDensityRight[i][j][k] += additionalChargeDensityMinusRight[i][j][k];
			}
		}
	}

	if ((rank == 0) && (verbosity > 1)) {
		printf("sum densityParameters\n");
	}
	sumCellParameters();
	sumCellVectorParameters();

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				//fprintf(debugFile, "charge %15.10g proton %15.10g electron %15.10g\n", chargeDensity[i][j][k], protonConcentration[i][j][k], electronConcentration[i][j][k]);
				for (int t = 0; t < typesNumber; ++t) {
					if (particleConcentrations[t][i][j][k] > 0) {
						particleBulkVelocities[t][i][j][k] = particleBulkVelocities[t][i][j][k] / (particleConcentrations[t][i][j][k] * types[t].mass);
						double gamma = sqrt((particleBulkVelocities[t][i][j][k].scalarMult(
							particleBulkVelocities[t][i][j][k]) / speed_of_light_normalized_sqr) + 1);
						particleBulkVelocities[t][i][j][k] = particleBulkVelocities[t][i][j][k] / gamma;
					}
				}
			}
		}
	}
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating density parameters time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::sumChargeDensityHat() {
	if (nprocs > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum charge density hat rank = %d\n", rank);
		double* inBufferRight = new double[(2 + additionalBinNumber) * ynumber * znumber];
		double* outBufferRight = new double[(2 + additionalBinNumber) * ynumber * znumber];
		double* inBufferLeft = new double[(2 + additionalBinNumber) * ynumber * znumber];
		double* outBufferLeft = new double[(2 + additionalBinNumber) * ynumber * znumber];

		//MPI_Barrier(MPI_COMM_WORLD);
		if ((rank == 0) && (verbosity > 2)) printf("sending left sum charge density hat rank = %d\n", rank);
		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				sendCellParametersLeft(chargeDensityHat, additionalChargeDensityHatLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendCellParametersLeft(chargeDensityHat, additionalChargeDensityHatLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if ((verbosity > 2)) printf("receiving right in sum charge density hat rank = %d\n", rank);
		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				receiveCellParametersRight(tempCellParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveCellParametersRight(tempCellParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if ((rank == 0) && (verbosity > 2)) printf("sending right sum charge density hat rank = %d\n", rank);
		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				sendCellParametersRight(chargeDensityHat, additionalChargeDensityHatRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendCellParametersRight(chargeDensityHat, additionalChargeDensityHatRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("receiving left in sum charge density hat rank = %d\n", rank);
		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				receiveCellParametersLeft(tempCellParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveCellParametersLeft(tempCellParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}
		sumCellTempParameters(chargeDensityHat);

		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("deleting buffer in sum charge density hat rank = %d\n", rank);
		delete[] inBufferRight;
		delete[] outBufferRight;
		delete[] inBufferLeft;
		delete[] outBufferLeft;
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					chargeDensityHat[xnumber][j][k] += chargeDensityHat[1][j][k];
					chargeDensityHat[1][j][k] = chargeDensityHat[xnumber][j][k];

					chargeDensityHat[xnumber - 1][j][k] += chargeDensityHat[0][j][k];
					chargeDensityHat[0][j][k] = chargeDensityHat[xnumber - 1][j][k];

					for (int i = 0; i < additionalBinNumber; ++i) {
						chargeDensityHat[xnumber - 2 - i][j][k] += additionalChargeDensityHatLeft[i][j][k];
						chargeDensityHat[2 + i][j][k] += additionalChargeDensityHatRight[i][j][k];
					}
				}
			}
		}
	}
}

void Simulation::sumCellParameters() {
	if (nprocs > 1) {
		double* inBufferRight = new double[(2 + additionalBinNumber) * ynumber * znumber];
		double* outBufferRight = new double[(2 + additionalBinNumber) * ynumber * znumber];
		double* inBufferLeft = new double[(2 + additionalBinNumber) * ynumber * znumber];
		double* outBufferLeft = new double[(2 + additionalBinNumber) * ynumber * znumber];

		for (int t = 0; t < typesNumber; ++t) {
			//MPI_Barrier(MPI_COMM_WORLD);
			if (types[t].particlesPerBin > 0) {
				if (rank == 0) {
					if (boundaryConditionType == PERIODIC) {
						sendCellParametersLeft(particleConcentrations[t], additionalParticleConcentrationsLeft[t], outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
					}
				} else {
					sendCellParametersLeft(particleConcentrations[t], additionalParticleConcentrationsLeft[t], outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
				}

				if (rank == nprocs - 1) {
					if (boundaryConditionType == PERIODIC) {
						receiveCellParametersRight(tempCellParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
					}
				} else {
					receiveCellParametersRight(tempCellParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
				}


				if (rank == nprocs - 1) {
					if (boundaryConditionType == PERIODIC) {
						sendCellParametersRight(particleConcentrations[t], additionalParticleConcentrationsRight[t], outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
					}
				} else {
					sendCellParametersRight(particleConcentrations[t], additionalParticleConcentrationsRight[t], outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
				}


				//MPI_Barrier(MPI_COMM_WORLD);
				if (rank == 0) {
					if (boundaryConditionType == PERIODIC) {
						receiveCellParametersLeft(tempCellParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
					}
				} else {
					receiveCellParametersLeft(tempCellParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
				}

				sumCellTempParameters(particleConcentrations[t]);
			}
			//MPI_Barrier(MPI_COMM_WORLD);
		}


		//MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				sendCellParametersLeft(chargeDensity, additionalChargeDensityLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendCellParametersLeft(chargeDensity, additionalChargeDensityLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				receiveCellParametersRight(tempCellParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveCellParametersRight(tempCellParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				sendCellParametersRight(chargeDensity, additionalChargeDensityRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendCellParametersRight(chargeDensity, additionalChargeDensityRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		//MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				receiveCellParametersLeft(tempCellParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveCellParametersLeft(tempCellParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}

		sumCellTempParameters(chargeDensity);
		//MPI_Barrier(MPI_COMM_WORLD);

		if (rank == nprocs - 1 && boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				chargeDensity[xnumber - 1][j][k] = chargeDensity[xnumber - 1][j][k] + chargeDensity[xnumber][j][k];
				chargeDensity[xnumber][j][k] = 0;
				for (int i = 0; i < additionalBinNumber; ++i) {
					chargeDensity[xnumber - 1][j][k] = chargeDensity[xnumber - 1][j][k] + additionalChargeDensityRight[i][j][k];
				}
				for(int t = 0; t < typesNumber; ++t) {
					particleConcentrations[t][xnumber - 1][j][k] = particleConcentrations[t][xnumber - 1][j][k] + particleConcentrations[t][xnumber][j][k];
					particleConcentrations[t][xnumber][j][k] = 0;
					for (int i = 0; i < additionalBinNumber; ++i) {
						particleConcentrations[t][xnumber - 1][j][k] = particleConcentrations[t][xnumber - 1][j][k] + additionalParticleConcentrationsRight[t][i][j][k];
					}
				}
			}
		}
	}

		delete[] inBufferRight;
		delete[] outBufferRight;
		delete[] inBufferLeft;
		delete[] outBufferLeft;
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int t = 0; t < typesNumber; ++t) {
						particleConcentrations[t][xnumber][j][k] += particleConcentrations[t][1][j][k];
						particleConcentrations[t][1][j][k] = particleConcentrations[t][xnumber][j][k];
						particleConcentrations[t][xnumber - 1][j][k] += particleConcentrations[t][0][j][k];
						particleConcentrations[t][0][j][k] = particleConcentrations[t][xnumber - 1][j][k];
					}
					chargeDensity[xnumber][j][k] += chargeDensity[1][j][k];
					chargeDensity[1][j][k] = chargeDensity[xnumber][j][k];


					chargeDensity[xnumber - 1][j][k] += chargeDensity[0][j][k];
					chargeDensity[0][j][k] = chargeDensity[xnumber - 1][j][k];

					for (int i = 0; i < additionalBinNumber; ++i) {

						chargeDensity[xnumber - 2 - i][j][k] += additionalChargeDensityLeft[i][j][k];
						chargeDensity[2 + i][j][k] += additionalChargeDensityRight[i][j][k];

						for (int t = 0; t < typesNumber; ++t) {
							particleConcentrations[t][xnumber - 2 - i][j][k] += additionalParticleConcentrationsLeft[t][i][j][k];
							particleConcentrations[t][2 + i][j][k] += additionalParticleConcentrationsRight[t][i][j][k];
						}
					}
				}
			}
		} else {
			//manualy set zero to the last bin
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int t = 0; t < typesNumber; ++t) {
						particleConcentrations[t][xnumber - 1][j][k] += particleConcentrations[t][xnumber][j][k];
						particleConcentrations[t][xnumber][j][k] = 0;
					}
					chargeDensity[xnumber - 1][j][k] += chargeDensity[xnumber][j][k];
					chargeDensity[xnumber][j][k] = 0;
					for (int i = 0; i < additionalBinNumber; ++i) {
						chargeDensity[xnumber - 1][j][k] += additionalChargeDensityRight[i][j][k];
						for (int t = 0; t < typesNumber; ++t) {
							particleConcentrations[t][xnumber - 1][j][k] += additionalParticleConcentrationsRight[t][i][j][k];
						}
					}
				}
			}
		}
	}

}

void Simulation::sumCellTempParameters(double*** array) {
	if (rank > 0 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				array[0][j][k] = array[0][j][k] + tempCellParameterLeft[0][j][k];
				array[1][j][k] = array[1][j][k] + tempCellParameterLeft[1][j][k];
				for (int i = 0; i < additionalBinNumber; ++i) {
					array[2 + i][j][k] = array[2 + i][j][k] + tempCellParameterLeft[2 + i][j][k];
				}
			}
		}
	}

	if (rank < nprocs - 1 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				array[xnumber - 1][j][k] = array[xnumber - 1][j][k] + tempCellParameterRight[0][j][k];
				array[xnumber][j][k] = array[xnumber][j][k] + tempCellParameterRight[1][j][k];

				for (int i = 0; i < additionalBinNumber; ++i) {
					array[xnumber - 2 - i][j][k] = array[xnumber - 2 - i][j][k] + tempCellParameterRight[2 + i][j][k];
				}
			}
		}
	}
}

void Simulation::sumCellVectorParameters() {
	if (nprocs > 1) {
		double* inBufferRight = new double[(2 + additionalBinNumber) * 3 * ynumber * znumber];
		double* outBufferRight = new double[(2 + additionalBinNumber) * 3 * ynumber * znumber];
		double* inBufferLeft = new double[(2 + additionalBinNumber) * 3 * ynumber * znumber];
		double* outBufferLeft = new double[(2 + additionalBinNumber) * 3 * ynumber * znumber];

		for (int t = 0; t < typesNumber; ++t) {
			//MPI_Barrier(MPI_COMM_WORLD);
			if (types[t].particlesPerBin > 0) {
				if (rank == 0) {
					if (boundaryConditionType == PERIODIC) {
						sendCellVectorParametersLeft(particleBulkVelocities[t], additionalParticleBulkVelocitiesLeft[t], outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
					}
				} else {
					sendCellVectorParametersLeft(particleBulkVelocities[t], additionalParticleBulkVelocitiesLeft[t], outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
				}


				if (rank == nprocs - 1) {
					if (boundaryConditionType == PERIODIC) {
						receiveCellVectorParametersRight(tempCellVectorParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
					}
				} else {
					receiveCellVectorParametersRight(tempCellVectorParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
				}

				if (rank == nprocs - 1) {
					if (boundaryConditionType == PERIODIC) {
						sendCellVectorParametersRight(particleBulkVelocities[t], additionalParticleBulkVelocitiesRight[t], outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
					}
				} else {
					sendCellVectorParametersRight(particleBulkVelocities[t], additionalParticleBulkVelocitiesRight[t], outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
				}

				//MPI_Barrier(MPI_COMM_WORLD);

				if (rank == 0) {
					if (boundaryConditionType == PERIODIC) {
						receiveCellVectorParametersLeft(tempCellVectorParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
					}
				} else {
					receiveCellVectorParametersLeft(tempCellVectorParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
				}

				sumCellTempVectorParameters(particleBulkVelocities[t]);
			}
			//MPI_Barrier(MPI_COMM_WORLD);
		}

		delete[] inBufferRight;
		delete[] outBufferRight;
		delete[] inBufferLeft;
		delete[] outBufferLeft;
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int t = 0; t < typesNumber; ++t) {
						particleBulkVelocities[t][xnumber][j][k] += particleBulkVelocities[t][1][j][k];
						particleBulkVelocities[t][1][j][k] = particleBulkVelocities[t][xnumber][j][k];
						particleBulkVelocities[t][xnumber - 1][j][k] += particleBulkVelocities[t][0][j][k];
						particleBulkVelocities[t][0][j][k] = particleBulkVelocities[t][xnumber - 1][j][k];
					}

					for (int i = 0; i < additionalBinNumber; ++i) {
						for (int t = 0; t < typesNumber; ++t) {
							particleBulkVelocities[t][xnumber - 2 - i][j][k] += additionalParticleBulkVelocitiesLeft[t][i][j][k];
							particleBulkVelocities[t][2 + i][j][k] += additionalParticleBulkVelocitiesRight[t][i][j][k];
						}
					}
				}
			}
		}
	}
}

void Simulation::sumCellTempVectorParameters(Vector3d*** array) {
	if (rank > 0 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				array[0][j][k] = array[0][j][k] + tempCellVectorParameterLeft[0][j][k];
				array[1][j][k] = array[1][j][k] + tempCellVectorParameterLeft[1][j][k];
				for (int i = 0; i < additionalBinNumber; ++i) {
					array[2 + i][j][k] = array[2 + i][j][k] + tempCellVectorParameterLeft[2 + i][j][k];
				}
			}
		}
	}

	if (rank < nprocs - 1 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				array[xnumber - 1][j][k] = array[xnumber - 1][j][k] + tempCellVectorParameterRight[0][j][k];
				array[xnumber][j][k] = array[xnumber][j][k] + tempCellVectorParameterRight[1][j][k];

				for (int i = 0; i < additionalBinNumber; ++i) {
					array[xnumber - 2 - i][j][k] = array[xnumber - 2 - i][j][k] + tempCellVectorParameterRight[2 + i][j][k];
				}
			}
		}
	}
}

void Simulation::sumCellMatrixParameters() {
	if (nprocs > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum cell matrix parameters rank = %d\n", rank);
		double* inBufferRight = new double[(2 + additionalBinNumber) * 9 * ynumber * znumber];
		double* outBufferRight = new double[(2 + additionalBinNumber) * 9 * ynumber * znumber];
		double* inBufferLeft = new double[(2 + additionalBinNumber) * 9 * ynumber * znumber];
		double* outBufferLeft = new double[(2 + additionalBinNumber) * 9 * ynumber * znumber];

		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("sending left  pressure tensor sum cell matrix parameters rank = %d\n", rank);
		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				sendCellMatrixParametersLeft(pressureTensor, additionalPressureTensorLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendCellMatrixParametersLeft(pressureTensor, additionalPressureTensorLeft, outBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if ((verbosity > 2)) printf("receiving right pressure tensor in sum cell matrix parameters rank = %d\n", rank);
		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				receiveCellMatrixParametersRight(tempCellMatrixParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveCellMatrixParametersRight(tempCellMatrixParameterRight, inBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		if ((verbosity > 2)) printf("sending right pressure tensor sum cell matrix parameters rank = %d\n", rank);
		if (rank == nprocs - 1) {
			if (boundaryConditionType == PERIODIC) {
				sendCellMatrixParametersRight(pressureTensor, additionalPressureTensorRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			sendCellMatrixParametersRight(pressureTensor, additionalPressureTensorRight, outBufferRight, xnumber, ynumber, znumber, additionalBinNumber);
		}

		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("receiving left pressure tensor in sum cell matrix parameters rank = %d\n", rank);
		if (rank == 0) {
			if (boundaryConditionType == PERIODIC) {
				receiveCellMatrixParametersLeft(tempCellMatrixParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
			}
		} else {
			receiveCellMatrixParametersLeft(tempCellMatrixParameterLeft, inBufferLeft, xnumber, ynumber, znumber, additionalBinNumber);
		}

		sumCellTempMatrixParameters(pressureTensor);
		//MPI_Barrier(MPI_COMM_WORLD);
		if ((verbosity > 2)) printf("deleting buffer in sum cell matrix parameters rank = %d\n", rank);
		delete[] inBufferRight;
		delete[] outBufferRight;
		delete[] inBufferLeft;
		delete[] outBufferLeft;
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						for (int m = 0; m < 3; ++m) {
							pressureTensor[xnumber][j][k].matrix[l][m] += pressureTensor[1][j][k].matrix[l][m];
							pressureTensor[1][j][k].matrix[l][m] = pressureTensor[xnumber][j][k].matrix[l][m];

							pressureTensor[xnumber - 1][j][k].matrix[l][m] += pressureTensor[0][j][k].matrix[l][m];
							pressureTensor[0][j][k].matrix[l][m] = pressureTensor[xnumber - 1][j][k].matrix[l][m];

							for (int i = 0; i < additionalBinNumber; ++i) {
								pressureTensor[xnumber - 2 - i][j][k].matrix[l][m] += additionalPressureTensorLeft[i][j][k].matrix[l][m];
								pressureTensor[2 + i][j][k].matrix[l][m] += additionalPressureTensorRight[i][j][k].matrix[l][m];
							}
						}
					}
				}
			}
		}
	}

}

void Simulation::sumCellTempMatrixParameters(Matrix3d*** array) {
	if (rank > 0 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				array[0][j][k] = array[0][j][k] + tempCellMatrixParameterLeft[0][j][k];
				array[1][j][k] = array[1][j][k] + tempCellMatrixParameterLeft[1][j][k];
				for (int i = 0; i < additionalBinNumber; ++i) {
					array[2 + i][j][k] = array[2 + i][j][k] + tempCellMatrixParameterLeft[2 + i][j][k];
				}
			}
		}
	}

	if (rank < nprocs - 1 || boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				array[xnumber - 1][j][k] = array[xnumber - 1][j][k] + tempCellMatrixParameterRight[0][j][k];
				array[xnumber][j][k] = array[xnumber][j][k] + tempCellMatrixParameterRight[1][j][k];

				for (int i = 0; i < additionalBinNumber; ++i) {
					array[xnumber - 2 - i][j][k] = array[xnumber - 2 - i][j][k] + tempCellMatrixParameterRight[2 + i][j][k];
				}
			}
		}
	}
}

void Simulation::updateEnergy() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	particleEnergy = 0;
	electricFieldEnergy = 0;
	magneticFieldEnergy = 0;
	energy = 0;

	globalMomentum = Vector3d(0, 0, 0);
	if (currentIteration % writeParameter == 0) {
		//particlesNumber = particles.size();
		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					Vector3d E = Efield[i][j][k];
					electricFieldEnergy += E.scalarMult(E) * volumeE(i, j, k) / (8 * pi);
					if (boundaryConditionType == PERIODIC || rank > 0 || i > 1) {
						electricFieldEnergy -= E0.scalarMult(E0) * volumeB(i, j, k) / (8 * pi);
					}
				}
			}
		}

		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					Vector3d B = Bfield[i][j][k];
					magneticFieldEnergy += (B.scalarMult(B)) * volumeB(i, j, k) / (8 * pi);
					if (boundaryConditionType == PERIODIC || rank > 0 || i > 0) {
						magneticFieldEnergy -= B0.scalarMult(B0) * volumeB(i, j, k) / (8 * pi);
					}
				}
			}
		}

		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					Vector3d E = ((Efield[i][j][k] + Efield[i][j][k + 1] + Efield[i][j + 1][k] + Efield[i][j + 1][k + 1] + Efield[i + 1][j][k] + Efield[i + 1][j][k + 1] + Efield[i + 1][j + 1][k] + Efield[i + 1][j + 1][k + 1]) / 8);
					Vector3d B = Bfield[i][j][k];

					globalMomentum += (E.vectorMult(B) / (4 * pi * speed_of_light_normalized)) * volumeB(i, j, k);
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
		}


		particleEnergy *= sqr(scaleFactor / plasma_period);
		electricFieldEnergy *= sqr(scaleFactor / plasma_period);
		magneticFieldEnergy *= sqr(scaleFactor / plasma_period);
		globalMomentum = globalMomentum * scaleFactor / plasma_period;


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

			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					theoreticalEnergy -= (Efield[xnumber][j][k].vectorMult(Bfield[xnumber][j][k])).x * deltaT *
						speed_of_light_normalized * deltaZ * deltaY *
						sqr(scaleFactor / plasma_period) / (4 * pi);
					theoreticalEnergy +=
						(Efield[1][j][k].vectorMult(Bfield[1][j][k])).x * deltaT * speed_of_light_normalized *
						deltaZ * deltaY * sqr(scaleFactor / plasma_period) / (4 * pi);

					theoreticalMomentum -=
						(Efield[xnumber][j][k].vectorMult(Bfield[xnumber][j][k])) * deltaT * deltaZ * deltaY *
						scaleFactor / plasma_period / (4 * pi);
					theoreticalMomentum +=
						(Efield[1][j][k].vectorMult(Bfield[1][j][k])) * deltaT * deltaZ * deltaY * scaleFactor /
						plasma_period / (4 * pi);
				}
			}
		}
	}

	if (currentIteration % writeParameter == 0) {
		double buffer[12];
		if (rank != 0) {
			buffer[0] = particleEnergy;
			buffer[1] = electricFieldEnergy;
			buffer[2] = magneticFieldEnergy;
			buffer[3] = energy;
			buffer[4] = globalMomentum.x;
			buffer[5] = globalMomentum.y;
			buffer[6] = globalMomentum.z;
			buffer[7] = theoreticalEnergy;
			buffer[8] = theoreticalMomentum.x;
			buffer[9] = theoreticalMomentum.y;
			buffer[10] = theoreticalMomentum.z;
			//buffer[11] = (double)particlesNumber;

			MPI_Send(buffer, 11, MPI_DOUBLE, 0, MPI_SEND_GENERAL_PARAMETERS, MPI_COMM_WORLD);
		} else {
			generalTheoreticalEnergy = theoreticalEnergy;
			generalTheoreticalMomentum = theoreticalMomentum;
			for (int i = 1; i < nprocs; ++i) {
				MPI_Status status;
				MPI_Recv(buffer, 11, MPI_DOUBLE, i, MPI_SEND_GENERAL_PARAMETERS, MPI_COMM_WORLD, &status);
				particleEnergy += buffer[0];
				electricFieldEnergy += buffer[1];
				magneticFieldEnergy += buffer[2];
				energy += buffer[3];
				globalMomentum.x += buffer[4];
				globalMomentum.y += buffer[5];
				globalMomentum.z += buffer[6];
				generalTheoreticalEnergy += buffer[7];
				generalTheoreticalMomentum.x += buffer[8];
				generalTheoreticalMomentum.y += buffer[9];
				generalTheoreticalMomentum.z += buffer[10];
				//particlesNumber += buffer[11];
			}
		}
	}

	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating energy time = %g sec\n", procTime / CLOCKS_PER_SEC);
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
		MPI_Finalize();
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
		MPI_Finalize();
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
		MPI_Finalize();
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


	return chargeDensityHat[i][j][k];
}

void Simulation::updateParameters() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	maxBfield = Bfield[0][0][0] - B0;
	maxEfield = Efield[0][0][0];
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Vector3d normalField = Bfield[i][j][k] - B0;
				//if (Bfield[i][j][k].norm() > maxBfield.norm()) {
				if (normalField.norm() > maxBfield.norm()) {
					//maxBfield = Bfield[i][j][k];
					maxBfield = normalField;
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

	MPI_Barrier(MPI_COMM_WORLD);
	double temp[7];
	temp[0] = maxBfield.x;
	temp[1] = maxBfield.y;
	temp[2] = maxBfield.z;
	temp[3] = maxEfield.x;
	temp[4] = maxEfield.y;
	temp[5] = maxEfield.z;
	if (rank != 0) {
		MPI_Send(temp, 6, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD);
	} else {
		for (int i = 1; i < nprocs; ++i) {
			MPI_Status status;
			MPI_Recv(temp, 6, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
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
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {
		temp[0] = maxBfield.x;
		temp[1] = maxBfield.y;
		temp[2] = maxBfield.z;
		temp[3] = maxEfield.x;
		temp[4] = maxEfield.y;
		temp[5] = maxEfield.z;
		for (int i = 1; i < nprocs; ++i) {
			MPI_Send(temp, 6, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD);
		}
	} else {
		MPI_Status status;
		MPI_Recv(temp, 6, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
		maxBfield = Vector3d(temp[0], temp[1], temp[2]);
		maxEfield = Vector3d(temp[3], temp[4], temp[5]);
	}
	MPI_Barrier(MPI_COMM_WORLD);

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

	MPI_Barrier(MPI_COMM_WORLD);

	if (nprocs > 1) {
		if (rank != 0) {
			for (int i = 0; i < typesNumber; ++i) {
				outBuffer[3 * i] = parallelV2[i];
				outBuffer[3 * i + 1] = normalV2[i];
				outBuffer[3 * i + 2] = types[i].generalWeight;
			}
			MPI_Send(outBuffer, 3 * typesNumber, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD);
		} else {
			for (int i = 1; i < nprocs; ++i) {
				MPI_Status status;
				MPI_Recv(inBuffer, 3 * typesNumber, MPI_DOUBLE, i, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD, &status);
				for (int t = 0; t < typesNumber; ++t) {
					parallelV2[t] += inBuffer[3 * t];
					normalV2[t] += inBuffer[3 * t + 1];
					types[t].generalWeight += inBuffer[3 * t + 2];
				}
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

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

	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating anisotropy time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::updateExternalFlux() {
	double alfvenV;
	if (density <= 1E-100) {
		alfvenV = speed_of_light_normalized;
	} else {
		alfvenV = B0.norm() / sqrt(4 * pi * density);
	}
	double concentration = density / (massProton + massElectron);
	double phaseV = 2 * alfvenV;
	double kw = 2 * pi / xsize;
	double omega = kw * phaseV;

	for (int i = 0; i < xnumber + 2; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				externalElectricFlux[i][j][k] = Vector3d(0, 0, 1.0) * extJ * cos(kw * xgrid[i] - omega * time);
				alertNaNOrInfinity(externalElectricFlux[i][j][k].x, "externalFlux.x = NaN\n");
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

void Simulation::updateShockWaveX() {
	double tempShockWavex[1];
	tempShockWavex[0] = -1.0;
	if (rank < nprocs - 1) {
		MPI_Status status;
		MPI_Recv(tempShockWavex, 1, MPI_DOUBLE, rank + 1, MPI_SEND_DOUBLE_NUMBER_LEFT, MPI_COMM_WORLD, &status);
		if (tempShockWavex[0] < 0) {
			for (int i = xnumber - 1; i > 0; i--) {
				if (particleBulkVelocities[1][i][0][0].x > V0.x / 2) {
					tempShockWavex[0] = xgrid[i + 1];
					break;
				}
			}
		}
		if (rank > 0) {
			MPI_Send(tempShockWavex, 1, MPI_DOUBLE, rank - 1, MPI_SEND_DOUBLE_NUMBER_LEFT, MPI_COMM_WORLD);
		}
	} else {
		for (int i = xnumber - 1; i > 0; i--) {
			if (particleBulkVelocities[1][i][0][0].x > V0.x / 2) {
				tempShockWavex[0] = xgrid[i + 1];
				break;
			}
		}
		if (rank > 0) {
			MPI_Send(tempShockWavex, 1, MPI_DOUBLE, rank - 1, MPI_SEND_DOUBLE_NUMBER_LEFT, MPI_COMM_WORLD);
		}
	}

	if (tempShockWavex[0] < shockWaveX) {
		tempShockWavex[0] = shockWaveX;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(tempShockWavex, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	shockWaveX = tempShockWavex[0];
}

void Simulation::splitParticles() {
	int tempParticlesNumber[1];
	for (int i = 0; i < nprocs; ++i) {
		if (i == rank) {
			if (rank > 0) {
				MPI_Status status;
				MPI_Recv(tempParticlesNumber, 1, MPI_INT, rank - 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD, &status);
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
				MPI_Send(tempParticlesNumber, 1, MPI_INT, rank + 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD);
			}
		}
	}
	tempParticlesNumber[0] = particlesNumber;
	MPI_Bcast(tempParticlesNumber, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD);
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

