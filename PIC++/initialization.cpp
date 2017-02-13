#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "dichotomousSolver.h"
#include "specialmath.h"
#include "util.h"
#include "output.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "matrixElement.h"
#include "particle.h"
#include "random.h"
#include "simulation.h"

Simulation::Simulation() {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	newlyStarted = false;
	preserveChargeGlobal = true;
	arrayCreated = false;
	outputDir = std::string(outputDirectory);

	maxEfield = Vector3d(0, 0, 0);
	maxBfield = Vector3d(0, 0, 0);

	shockWavePoint = 0;
	chargeBalance = 0;

	Kronecker = Matrix3d(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				LeviCivita[i][j][k] = 0;
			}
		}
	}
	LeviCivita[0][1][2] = 1.0;
	LeviCivita[0][2][1] = -1.0;
	LeviCivita[1][0][2] = -1.0;
	LeviCivita[1][2][0] = 1.0;
	LeviCivita[2][0][1] = 1.0;
	LeviCivita[2][1][0] = -1.0;
	typesNumber = 8;
	additionalBinNumber = (splineOrder + 1) / 2;
	verbosity = 0;
	//types = new ParticleTypeContainer[typesNumber];
	//concentrations = new double[typesNumber];
	//particlesPerBin = new int[typesNumber];
	shockWaveX = -1.0;
	protonNumber = 0;
	protonNumber1 = 0;
	protonNumber2 = 0;
	protonNumber3 = 0;
	protonNumber4 = 0;
	protonNumber5 = 0;
	protonNumber6 = 0;
	protonNumber7 = 0;
	protonNumber8 = 0;
	protonNumber9 = 0;
	electronNumber = 0;
	electronNumber1 = 0;
	electronNumber2 = 0;
	electronNumber3 = 0;
	electronNumber4 = 0;
	electronNumber5 = 0;
	electronNumber6 = 0;
	electronNumber7 = 0;
	electronNumber8 = 0;
	electronNumber9 = 0;
}

void Simulation::setSpaceForProc() {
	int tempXnumber = ((xnumberGeneral) / nprocs) + 1;
	int modXnumber = (xnumberGeneral) % nprocs;

	if (rank >= nprocs - modXnumber) {
		xnumber = tempXnumber + 1;
		firstAbsoluteXindex = xnumberGeneral - (xnumber - 1) * (nprocs - rank);
	}
	else {
		xnumber = tempXnumber;
		firstAbsoluteXindex = (xnumber - 1) * rank;

	}

	//todo boundary conditiontype
	if (nprocs > 5) {
		int firstSmallRegions = 2 * nprocs / 3;
		int firstSmallXnumber = xnumberGeneral / 3;
		int lastLargeRegions = nprocs - firstSmallRegions;
		int lastLargeXnumber = xnumberGeneral - firstSmallXnumber;
		if (rank < firstSmallRegions) {
			tempXnumber = firstSmallXnumber / firstSmallRegions + 1;
			modXnumber = firstSmallXnumber % firstSmallRegions;
			if (rank >= firstSmallRegions - modXnumber) {
				xnumber = tempXnumber + 1;
				firstAbsoluteXindex = firstSmallXnumber - (xnumber - 1) * (firstSmallRegions - rank);
			}
			else {
				xnumber = tempXnumber;
				firstAbsoluteXindex = (xnumber - 1) * rank;
			}
		}
		else {
			tempXnumber = lastLargeXnumber / lastLargeRegions + 1;
			modXnumber = lastLargeXnumber % lastLargeRegions;
			if (rank >= nprocs - modXnumber) {
				xnumber = tempXnumber + 1;
				firstAbsoluteXindex = xnumberGeneral - (xnumber - 1) * (nprocs - rank);
			}
			else {
				xnumber = tempXnumber;
				firstAbsoluteXindex = firstSmallXnumber + (xnumber - 1) * (rank - firstSmallRegions);
			}
		}

	}
	//printf("xnumber = %d\n", xnumber);


	xsize = xnumber * xsizeGeneral / xnumberGeneral;
	ysize = ysizeGeneral;
	zsize = zsizeGeneral;

	deltaX = xsize / (xnumber);
	deltaY = ysize / (ynumber);
	deltaZ = zsize / (znumber);

	deltaX2 = deltaX * deltaX;
	deltaY2 = deltaY * deltaY;
	deltaZ2 = deltaZ * deltaZ;

	/*if (rank >= nprocs - modXnumber) {
		leftX = 2 * xsizeGeneral - (xsize - deltaX) * (nprocs - rank) - deltaX;
	} else {
		leftX = xsizeGeneral + (xsize - deltaX) * rank - deltaX;
	}*/
	leftX = xsizeGeneral + (firstAbsoluteXindex - 1) * deltaX;
	rightX = leftX + xsize;
}

Simulation::Simulation(int xn, int yn, int zn, double xsizev, double ysizev, double zsizev, double temp, double Vx,
                       double Vy, double Vz, double Ex, double Ey, double Ez, double Bx, double By, double Bz,
                       int maxIterations, double maxTimeV, int typesNumberV, int* particlesPerBinV,
                       double* concentrationsV, int inType, int nprocsV, int verbosityV) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	nprocs = nprocsV;
	arrayCreated = false;
	timing = true;
	outputDir = std::string(outputDirectory);
	if (inType == 0) {
		inputType = CGS;
	}
	else if (inType == 1) {
		inputType = Theoretical;
	}
	else {
		printf("input type must be 1 or 0\n");
		fflush(stdout);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "input type must be 1 or 0\n");
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	debugMode = false;
	newlyStarted = true;
	preserveChargeGlobal = true;
	solverType = IMPLICIT; //не явный
	//solverType = EXPLICIT; //явный
	boundaryConditionType = PERIODIC;
	//boundaryConditionType = SUPER_CONDUCTOR_LEFT;
	maxwellEquationMatrixSize = 3;
	verbosity = verbosityV;

	typesNumber = typesNumberV;
	if (typesNumber != 8) {
		printf(
			"PIC++ support only 8 types of ions, typesNumber must be = 8. if you need less, ypu can initialize them with 0 concentration\n");
		fflush(stdout);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile,
		        "PIC++ support only 8 types of ions, typesNumber must be = 8. if you need less, ypu can initialize them with 0 concentration\n");
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	concentrations = concentrationsV;
	particlesPerBin = particlesPerBinV;

	currentIteration = 0;
	time = 0;
	particlesNumber = 0;

	particleEnergy = 0;
	electricFieldEnergy = 0;
	magneticFieldEnergy = 0;
	chargeBalance = 0;

	globalMomentum = Vector3d(0, 0, 0);


	theta = initialTheta;
	//eta = theta;
	eta = 0.5;

	xnumberGeneral = xn;
	ynumberGeneral = yn;
	znumberGeneral = zn;

	ynumber = yn;
	znumber = zn;

	xsizeGeneral = xsizev;
	ysizeGeneral = ysizev;
	zsizeGeneral = zsizev;

	setSpaceForProc();

	temperature = temp;

	maxIteration = maxIterations;
	maxTime = maxTimeV;

	extJ = 0;

	V0 = Vector3d(Vx, Vy, Vz);

	B0 = Vector3d(Bx, By, Bz);
	E0 = Vector3d(Ex, Ey, Ez);

	additionalBinNumber = (splineOrder + 1) / 2;

	if (inputType == CGS) {
		massProton = massProtonReal;
		massElectron = massElectronFactor * massElectronReal;
		massAlpha = massAlphaReal;
		massDeuterium = massDeuteriumReal;
		massHelium3 = massHelium3Real;

		double omega2 = 0;

		omega2 += 4 * pi * concentrationsV[0] * electron_charge * electron_charge / massElectron;
		omega2 += 4 * pi * concentrationsV[1] * electron_charge * electron_charge / massProton;
		omega2 += 4 * pi * concentrationsV[2] * electron_charge * electron_charge / massElectron;
		omega2 += 4 * pi * concentrationsV[3] * 4 * electron_charge * electron_charge / massAlpha;
		omega2 += 4 * pi * concentrationsV[4] * electron_charge * electron_charge / massDeuterium;
		omega2 += 4 * pi * concentrationsV[5] * 4 * electron_charge * electron_charge / massHelium3;

		double gamma = sqrt(1.0 - V0.scalarMult(V0) / (speed_of_light * speed_of_light));

		plasma_period = sqrt(1 / omega2) * (2 * pi) * gamma * sqrt(gamma);
		//plasma_period = sqrt(1 / omega2) * (2 * pi) * gamma * sqrt(gamma)/(2*pi);
		double thermal_momentum;
		if (kBoltzman * temperature > massElectron * speed_of_light * speed_of_light) {
			thermal_momentum = kBoltzman * temperature / speed_of_light;
		}
		else {
			thermal_momentum = sqrt(2 * massElectron * kBoltzman * temperature);
		}
		thermal_momentum += V0.norm() * massElectron;
		double gyro_radius = thermal_momentum * speed_of_light / (electron_charge * B0.norm());
		//if (B0.norm() <= 0) {
		//scaleFactor = 1.0;
		//}
		scaleFactor = speed_of_light * plasma_period;

		plasma_period = 1.0;
		scaleFactor = 1.0;

		//scaleFactor = xsize;


		E0 = E0 * (plasma_period * sqrt(scaleFactor));
		B0 = B0 * (plasma_period * sqrt(scaleFactor));
		V0 = V0 * plasma_period / scaleFactor;

		rescaleConstants();

		density = density * cube(scaleFactor);
		for (int i = 0; i < typesNumber; ++i) {
			concentrations[i] = concentrations[i] * cube(scaleFactor);
		}

		printf("scaleFactor = %lf\n", scaleFactor);

		deltaX /= scaleFactor;
		deltaY /= scaleFactor;
		deltaZ /= scaleFactor;
		deltaX2 /= scaleFactor * scaleFactor;
		deltaY2 /= scaleFactor * scaleFactor;
		deltaZ2 /= scaleFactor * scaleFactor;

		leftX /= scaleFactor;
		rightX /= scaleFactor;
		xsize /= scaleFactor;
		xsizeGeneral /= scaleFactor;
		ysize /= scaleFactor;
		ysizeGeneral /= scaleFactor;
		zsize /= scaleFactor;
		zsizeGeneral /= scaleFactor;
		if (rank == 0) printf("xsize/scaleFactor = %lf\n", xsize);
		//fflush(stdout);
	}
	else {
		massProton = 1.0;
		massElectron = massProton / 256;
		massAlpha = massProton * massAlphaReal / massProtonReal;
		massDeuterium = massProton * massDeuteriumReal / massProtonReal;
		massHelium3 = massProton * massHelium3Real / massProtonReal;
		double protonScale = massProton / massProtonReal;
		plasma_period = cube(speed_of_light) / (protonScale * electron_charge);
		scaleFactor = speed_of_light * plasma_period;
		rescaleConstantsToTheoretical();
		double densityForUnits = massElectron * (massProton + massElectron) / (4 * pi * electron_charge_normalized * electron_charge_normalized);
		if (fabs(density - densityForUnits) / (densityForUnits + densityForUnits) > 1E-3) {
			if (rank == 0) printf("density must be changed\n");
			if (rank == 0) printf("density = %g, must be = %g\n", density, densityForUnits);
			fflush(stdout);
		}
	}

	maxEfield = Vector3d(0, 0, 0);
	maxBfield = Vector3d(0, 0, 0);
	shockWavePoint = 0;

	Kronecker = Matrix3d(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				LeviCivita[i][j][k] = 0;
			}
		}
	}
	LeviCivita[0][1][2] = 1.0;
	LeviCivita[0][2][1] = -1.0;
	LeviCivita[1][0][2] = -1.0;
	LeviCivita[1][2][0] = 1.0;
	LeviCivita[2][0][1] = 1.0;
	LeviCivita[2][1][0] = -1.0;
	shockWaveX = -1.0;
	//if(rank == 0) printf("end constructor\n");
	//fflush(stdout);
	protonNumber = 0;
	protonNumber1 = 0;
	protonNumber2 = 0;
	protonNumber3 = 0;
	protonNumber4 = 0;
	protonNumber5 = 0;
	protonNumber6 = 0;
	protonNumber7 = 0;
	protonNumber8 = 0;
	protonNumber9 = 0;
	electronNumber = 0;
	electronNumber1 = 0;
	electronNumber2 = 0;
	electronNumber3 = 0;
	electronNumber4 = 0;
	electronNumber5 = 0;
	electronNumber6 = 0;
	electronNumber7 = 0;
	electronNumber8 = 0;
	electronNumber9 = 0;
}

Simulation::Simulation(int xn, int yn, int zn, double xsizev, double ysizev, double zsizev, double temp, double Vx,
                       double Vy, double Vz, double Ex, double Ey, double Ez, double Bx, double By, double Bz,
                       int maxIterations, double maxTimeV, int typesNumberV, int* particlesPerBinV,
                       double* concentrationsV, int inType, int nprocsV, int verbosityV, double plasmaPeriodV,
                       double scaleFactorV) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	nprocs = nprocsV;
	arrayCreated = false;
	timing = true;
	outputDir = std::string(outputDirectory);
	if (inType == 0) {
		inputType = CGS;
	}
	else if (inType == 1) {
		inputType = Theoretical;
	}
	else {
		printf("input type must be 1 or 0\n");
		fflush(stdout);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "input type must be 1 or 0\n");
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	debugMode = false;
	newlyStarted = true;
	preserveChargeGlobal = true;
	solverType = IMPLICIT; //не явный
	//solverType = EXPLICIT; //явный
	boundaryConditionType = PERIODIC;
	//boundaryConditionType = SUPER_CONDUCTOR_LEFT;
	maxwellEquationMatrixSize = 3;
	verbosity = verbosityV;

	typesNumber = typesNumberV;
	if (typesNumber != 8) {
		printf(
			"PIC++ support only 8 types of ions, typesNumber must be = 8. if you need less, ypu can initialize them with 0 concentration\n");
		fflush(stdout);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile,
		        "PIC++ support only 8 types of ions, typesNumber must be = 8. if you need less, ypu can initialize them with 0 concentration\n");
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	concentrations = concentrationsV;
	particlesPerBin = particlesPerBinV;

	currentIteration = 0;
	time = 0;
	particlesNumber = 0;

	particleEnergy = 0;
	electricFieldEnergy = 0;
	magneticFieldEnergy = 0;
	chargeBalance = 0;

	globalMomentum = Vector3d(0, 0, 0);


	theta = initialTheta;
	//eta = theta;
	eta = 0.5;

	xnumberGeneral = xn;
	ynumberGeneral = yn;
	znumberGeneral = zn;

	ynumber = yn;
	znumber = zn;

	xsizeGeneral = xsizev;
	ysizeGeneral = ysizev;
	zsizeGeneral = zsizev;

	setSpaceForProc();

	temperature = temp;

	maxIteration = maxIterations;
	maxTime = maxTimeV;

	extJ = 0;

	V0 = Vector3d(Vx, Vy, Vz);

	B0 = Vector3d(Bx, By, Bz);
	E0 = Vector3d(Ex, Ey, Ez);

	additionalBinNumber = (splineOrder + 1) / 2;

	if (inputType == CGS) {
		massProton = massProtonReal;
		massElectron = massElectronFactor * massElectronReal;
		massAlpha = massAlphaReal;
		massDeuterium = massDeuteriumReal;
		massHelium3 = massHelium3Real;


		plasma_period = plasmaPeriodV;

		scaleFactor = scaleFactorV;


		rescaleConstants();


		printf("scaleFactor = %lf\n", scaleFactor);

		if (rank == 0) printf("xsize/scaleFactor = %lf\n", xsize);
		//fflush(stdout);
	}
	else {
		massProton = 1.0;
		massElectron = massProton / 256;
		massAlpha = massProton * massAlphaReal / massProtonReal;
		massDeuterium = massProton * massDeuteriumReal / massProtonReal;
		massHelium3 = massProton * massHelium3Real / massProtonReal;
		double protonScale = massProton / massProtonReal;
		plasma_period = plasmaPeriodV;
		scaleFactor = scaleFactorV;
		rescaleConstantsToTheoretical();
		double densityForUnits = massElectron * (massProton + massElectron) / (4 * pi * electron_charge_normalized * electron_charge_normalized);
		if (fabs(density - densityForUnits) / (densityForUnits + densityForUnits) > 1E-3) {
			if (rank == 0) printf("density must be changed\n");
			if (rank == 0) printf("density = %g, must be = %g\n", density, densityForUnits);
			fflush(stdout);
		};
	}

	maxEfield = Vector3d(0, 0, 0);
	maxBfield = Vector3d(0, 0, 0);
	shockWavePoint = 0;

	Kronecker = Matrix3d(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				LeviCivita[i][j][k] = 0;
			}
		}
	}
	LeviCivita[0][1][2] = 1.0;
	LeviCivita[0][2][1] = -1.0;
	LeviCivita[1][0][2] = -1.0;
	LeviCivita[1][2][0] = 1.0;
	LeviCivita[2][0][1] = 1.0;
	LeviCivita[2][1][0] = -1.0;
	shockWaveX = -1.0;
	if (rank == 0) printf("end constructor\n");
	fflush(stdout);
	protonNumber = 0;
	protonNumber1 = 0;
	protonNumber2 = 0;
	protonNumber3 = 0;
	protonNumber4 = 0;
	protonNumber5 = 0;
	protonNumber6 = 0;
	protonNumber7 = 0;
	protonNumber8 = 0;
	protonNumber9 = 0;
	electronNumber = 0;
	electronNumber1 = 0;
	electronNumber2 = 0;
	electronNumber3 = 0;
	electronNumber4 = 0;
	electronNumber5 = 0;
	electronNumber6 = 0;
	electronNumber7 = 0;
	electronNumber8 = 0;
	electronNumber9 = 0;
}

Simulation::~Simulation() {
	if (arrayCreated) {
		delete[] types;
		delete[] concentrations;
		delete[] particlesPerBin;

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					delete[] gmresOutput[i][j][k];
				}
				delete[] gmresOutput[i][j];
			}
			delete[] gmresOutput[i];
		}
		delete[] gmresOutput;

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					delete[] maxwellEquationMatrix[i][j][k];
					delete[] maxwellEquationRightPart[i][j][k];
				}
				delete[] maxwellEquationMatrix[i][j];
				delete[] maxwellEquationRightPart[i][j];
			}
			delete[] maxwellEquationMatrix[i];
			delete[] maxwellEquationRightPart[i];
		}
		delete[] maxwellEquationMatrix;
		delete[] maxwellEquationRightPart;

		for (int i = 0; i < xnumber + 2; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					delete[] divergenceCleaningField[i][j][k];
					delete[] divergenceCleanUpMatrix[i][j][k];
					delete[] divergenceCleanUpRightPart[i][j][k];
				}
				delete[] divergenceCleaningField[i][j];
				delete[] divergenceCleanUpMatrix[i][j];
				delete[] divergenceCleanUpRightPart[i][j];
			}
			delete[] divergenceCleaningField[i];
			delete[] divergenceCleanUpMatrix[i];
			delete[] divergenceCleanUpRightPart[i];
		}
		delete[] divergenceCleaningField;
		delete[] divergenceCleanUpMatrix;
		delete[] divergenceCleanUpRightPart;

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					delete[] divergenceCleaningPotential[i][j][k];
					delete[] tempDivergenceCleaningPotential[i][j][k];
				}
				delete[] divergenceCleaningPotential[i][j];
				delete[] tempDivergenceCleaningPotential[i][j];
				delete[] divergenceCleaningPotentialFourier[i][j];
			}
			delete[] divergenceCleaningPotential[i];
			delete[] tempDivergenceCleaningPotential[i];
			delete[] divergenceCleaningPotentialFourier[i];
		}
		delete[] divergenceCleaningPotential;
		delete[] tempDivergenceCleaningPotential;
		delete[] divergenceCleaningPotentialFourier;

		for (int i = 0; i < xnumber + 2; ++i) {
			for (int j = 0; j < ynumber + 1; ++j) {
				delete[] Efield[i][j];
				delete[] newEfield[i][j];
				delete[] tempEfield[i][j];
				delete[] smoothingEfield[i][j];
				delete[] explicitEfield[i][j];
				delete[] rotB[i][j];
				delete[] Ederivative[i][j];
			}
			delete[] Efield[i];
			delete[] newEfield[i];
			delete[] tempEfield[i];
			delete[] smoothingEfield[i];
			delete[] explicitEfield[i];
			delete[] rotB[i];
			delete[] Ederivative[i];
		}
		for (int i = 0; i < xnumber + 2; ++i) {
			for (int j = 0; j < ynumber + 1; ++j) {
				delete[] electricFlux[i][j];
				delete[] electricFluxMinus[i][j];
				delete[] externalElectricFlux[i][j];
				delete[] divPressureTensor[i][j];
				delete[] dielectricTensor[i][j];
			}

			delete[] electricFlux[i];
			delete[] electricFluxMinus[i];
			delete[] externalElectricFlux[i];
			delete[] divPressureTensor[i];
			delete[] dielectricTensor[i];
		}

		delete[] Efield;
		delete[] newEfield;
		delete[] tempEfield;
		delete[] smoothingEfield;
		delete[] explicitEfield;
		delete[] rotB;
		delete[] Ederivative;

		delete[] electricFlux;
		delete[] electricFluxMinus;
		delete[] externalElectricFlux;
		delete[] divPressureTensor;
		delete[] dielectricTensor;

		for (int t = 0; t < typesNumber; ++t) {
			for (int i = 0; i < xnumber + 1; ++i) {
				for (int j = 0; j < ynumber; ++j) {
					delete[] particleConcentrations[t][i][j];
					delete[] particleBulkVelocities[t][i][j];
				}
				delete[] particleConcentrations[t][i];
				delete[] particleBulkVelocities[t][i];
			}
			delete[] particleConcentrations[t];
			delete[] particleBulkVelocities[t];
		}
		delete[] particleConcentrations;
		delete[] particleBulkVelocities;

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				delete[] Bfield[i][j];
				delete[] newBfield[i][j];
				delete[] smoothingBfield[i][j];
				delete[] chargeDensity[i][j];
				delete[] chargeDensityMinus[i][j];
				delete[] pressureTensor[i][j];
				delete[] chargeDensityHat[i][j];
				delete[] tempCellParameter[i][j];
			}
			delete[] Bfield[i];
			delete[] newBfield[i];
			delete[] smoothingBfield[i];
			delete[] chargeDensity[i];
			delete[] chargeDensityMinus[i];
			delete[] pressureTensor[i];
			delete[] chargeDensityHat[i];
			delete[] tempCellParameter[i];
		}

		delete[] Bfield;
		delete[] newBfield;
		delete[] smoothingBfield;

		delete[] chargeDensity;
		delete[] chargeDensityMinus;
		delete[] pressureTensor;
		delete[] chargeDensityHat;
		delete[] tempCellParameter;

		delete[] xgrid;
		delete[] middleXgrid;
		delete[] ygrid;
		delete[] middleYgrid;
		delete[] zgrid;
		delete[] middleZgrid;

		delete[] rightEbuffer;
		delete[] rightEinBuffer;
		delete[] leftEbuffer;
		delete[] leftEoutBuffer;

		if (additionalBinNumber > 0) {
			for (int i = 0; i < additionalBinNumber; ++i) {
				for (int j = 0; j < ynumber + 1; ++j) {
					delete[] additionalEfieldLeft[i][j];
					delete[] additionalEfieldRight[i][j];
					delete[] additionalTempEfieldLeft[i][j];
					delete[] additionalTempEfieldRight[i][j];
					delete[] additionalNewEfieldLeft[i][j];
					delete[] additionalNewEfieldRight[i][j];
					delete[] additionalElectricFluxLeft[i][j];
					delete[] additionalElectricFluxMinusLeft[i][j];
					delete[] additionalElectricFluxRight[i][j];
					delete[] additionalElectricFluxMinusRight[i][j];
					delete[] additionalDielectricTensorLeft[i][j];
					delete[] additionalDielectricTensorRight[i][j];
					delete[] additionalDivPressureTensorLeft[i][j];
					delete[] additionalDivPressureTensorRight[i][j];
				}
				delete[] additionalEfieldLeft[i];
				delete[] additionalEfieldRight[i];
				delete[] additionalTempEfieldLeft[i];
				delete[] additionalTempEfieldRight[i];
				delete[] additionalNewEfieldLeft[i];
				delete[] additionalNewEfieldRight[i];
				delete[] additionalElectricFluxLeft[i];
				delete[] additionalElectricFluxMinusLeft[i];
				delete[] additionalElectricFluxRight[i];
				delete[] additionalElectricFluxMinusRight[i];
				delete[] additionalDielectricTensorLeft[i];
				delete[] additionalDielectricTensorRight[i];
				delete[] additionalDivPressureTensorLeft[i];
				delete[] additionalDivPressureTensorRight[i];

				for (int j = 0; j < ynumber; ++j) {
					delete[] additionalBfieldLeft[i][j];
					delete[] additionalBfieldRight[i][j];
					delete[] additionalNewBfieldLeft[i][j];
					delete[] additionalNewBfieldRight[i][j];
					delete[] additionalChargeDensityHatLeft[i][j];
					delete[] additionalChargeDensityHatRight[i][j];
					delete[] additionalChargeDensityLeft[i][j];
					delete[] additionalChargeDensityMinusLeft[i][j];
					delete[] additionalChargeDensityRight[i][j];
					delete[] additionalChargeDensityMinusRight[i][j];
					delete[] additionalPressureTensorLeft[i][j];
					delete[] additionalPressureTensorRight[i][j];
				}
				delete[] additionalBfieldLeft[i];
				delete[] additionalBfieldRight[i];
				delete[] additionalNewBfieldLeft[i];
				delete[] additionalNewBfieldRight[i];
				delete[] additionalChargeDensityHatLeft[i];
				delete[] additionalChargeDensityHatRight[i];
				delete[] additionalChargeDensityLeft[i];
				delete[] additionalChargeDensityMinusLeft[i];
				delete[] additionalChargeDensityRight[i];
				delete[] additionalChargeDensityMinusRight[i];
				delete[] additionalPressureTensorLeft[i];
				delete[] additionalPressureTensorRight[i];
			}
			delete[] additionalEfieldLeft;
			delete[] additionalEfieldRight;
			delete[] additionalTempEfieldLeft;
			delete[] additionalTempEfieldRight;
			delete[] additionalNewEfieldLeft;
			delete[] additionalNewEfieldRight;
			delete[] additionalElectricFluxLeft;
			delete[] additionalElectricFluxMinusLeft;
			delete[] additionalElectricFluxRight;
			delete[] additionalElectricFluxMinusRight;
			delete[] additionalDielectricTensorLeft;
			delete[] additionalDielectricTensorRight;
			delete[] additionalDivPressureTensorLeft;
			delete[] additionalDivPressureTensorRight;

			delete[] additionalBfieldLeft;
			delete[] additionalBfieldRight;
			delete[] additionalNewBfieldLeft;
			delete[] additionalNewBfieldRight;
			delete[] additionalChargeDensityHatLeft;
			delete[] additionalChargeDensityHatRight;
			delete[] additionalChargeDensityLeft;
			delete[] additionalChargeDensityMinusLeft;
			delete[] additionalChargeDensityRight;
			delete[] additionalChargeDensityMinusRight;
			delete[] additionalPressureTensorLeft;
			delete[] additionalPressureTensorRight;

			for (int t = 0; t < typesNumber; ++t) {
				for (int i = 0; i < additionalBinNumber; ++i) {
					for (int j = 0; j < ynumber; ++j) {
						delete[] additionalParticleConcentrationsLeft[t][i][j];
						delete[] additionalParticleConcentrationsRight[t][i][j];
						delete[] additionalParticleBulkVelocitiesLeft[t][i][j];
						delete[] additionalParticleBulkVelocitiesRight[t][i][j];
					}
					delete[] additionalParticleConcentrationsLeft[t][i];
					delete[] additionalParticleConcentrationsRight[t][i];
					delete[] additionalParticleBulkVelocitiesLeft[t][i];
					delete[] additionalParticleBulkVelocitiesRight[t][i];
				}
				delete[] additionalParticleConcentrationsLeft[t];
				delete[] additionalParticleConcentrationsRight[t];
				delete[] additionalParticleBulkVelocitiesLeft[t];
				delete[] additionalParticleBulkVelocitiesRight[t];
			}
		}
		delete[] additionalParticleConcentrationsLeft;
		delete[] additionalParticleConcentrationsRight;
		delete[] additionalParticleBulkVelocitiesLeft;
		delete[] additionalParticleBulkVelocitiesRight;

		for (int i = 0; i < particles.size(); ++i) {
			Particle* particle = particles[i];
			delete particle;
		}
		particles.clear();

		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				delete[] tempCellParameterLeft[i][j];
				delete[] tempCellParameterRight[i][j];
				delete[] tempCellVectorParameterLeft[i][j];
				delete[] tempCellVectorParameterRight[i][j];
				delete[] tempCellMatrixParameterLeft[i][j];
				delete[] tempCellMatrixParameterRight[i][j];
			}
			delete[] tempCellParameterLeft[i];
			delete[] tempCellParameterRight[i];
			delete[] tempCellVectorParameterLeft[i];
			delete[] tempCellVectorParameterRight[i];
			delete[] tempCellMatrixParameterLeft[i];
			delete[] tempCellMatrixParameterRight[i];
		}
		delete[] tempCellParameterLeft;
		delete[] tempCellParameterRight;
		delete[] tempCellVectorParameterLeft;
		delete[] tempCellVectorParameterRight;
		delete[] tempCellMatrixParameterLeft;
		delete[] tempCellMatrixParameterRight;

		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumber + 1; ++j) {
				delete[] tempNodeParameterLeft[i][j];
				delete[] tempNodeParameterRight[i][j];
				delete[] tempNodeVectorParameterLeft[i][j];
				delete[] tempNodeVectorParameterRight[i][j];
				delete[] tempNodeMatrixParameterLeft[i][j];
				delete[] tempNodeMatrixParameterRight[i][j];
			}
			delete[] tempNodeParameterLeft[i];
			delete[] tempNodeParameterRight[i];
			delete[] tempNodeVectorParameterLeft[i];
			delete[] tempNodeVectorParameterRight[i];
			delete[] tempNodeMatrixParameterLeft[i];
			delete[] tempNodeMatrixParameterRight[i];
		}
		delete[] tempNodeParameterLeft;
		delete[] tempNodeParameterRight;
		delete[] tempNodeVectorParameterLeft;
		delete[] tempNodeVectorParameterRight;
		delete[] tempNodeMatrixParameterLeft;
		delete[] tempNodeMatrixParameterRight;
	}
}

void Simulation::rescaleConstants() {
	kBoltzman_normalized = kBoltzman * plasma_period * plasma_period / (scaleFactor * scaleFactor);
	speed_of_light_normalized = speed_of_light * plasma_period / scaleFactor;
	speed_of_light_normalized_sqr = speed_of_light_normalized * speed_of_light_normalized;
	electron_charge_normalized = electron_charge * (plasma_period / sqrt(cube(scaleFactor)));
}

void Simulation::rescaleConstantsToTheoretical() {
	//kBoltzman_normalized = kBoltzman * plasma_period * plasma_period / (scaleFactor * scaleFactor);
	kBoltzman_normalized = kBoltzman * (plasma_period * plasma_period / (scaleFactor * scaleFactor)) * massProton / massProtonReal;
	speed_of_light_normalized = 1.0;
	speed_of_light_normalized_sqr = 1.0;
	electron_charge_normalized = 1.0;
}

void Simulation::initialize() {
	if (rank == 0) printf("initialization\n");
	fflush(stdout);
	if (rank == 0) printLog("initialization\n");

	xgrid[0] = leftX;

	for (int i = 1; i <= xnumber + 1; ++i) {
		xgrid[i] = xgrid[0] + i * deltaX;
	}

	xgrid[xnumber] = rightX;
	//xgrid[xnumber + 1] = rightX + deltaX;
	xgrid[xnumber + 1] = xgrid[0]+(xnumber+1)*deltaX;

	//printf("xgrid[0] = %lf xgrid[xnumber] = %lf\n", xgrid[0], xgrid[xnumber]);

	ygrid[0] = ysize;
	for (int j = 0; j <= ynumber; ++j) {
		ygrid[j] = ygrid[0] + j * deltaY;
	}
	ygrid[ynumber] = 2 * ysize;

	zgrid[0] = zsize;
	for (int k = 0; k <= znumber; ++k) {
		zgrid[k] = zgrid[0] + k * deltaZ;
	}
	zgrid[znumber] = 2 * zsize;

	for (int i = 0; i < xnumber + 1; ++i) {
		middleXgrid[i] = (xgrid[i] + xgrid[i + 1]) / 2;
	}

	for (int j = 0; j < ynumber; ++j) {
		middleYgrid[j] = (ygrid[j] + ygrid[j + 1]) / 2;
	}

	for (int k = 0; k < znumber; ++k) {
		middleZgrid[k] = (zgrid[k] + zgrid[k + 1]) / 2;
	}


	for (int i = 0; i < xnumber + 2; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k] = E0;
				newEfield[i][j][k] = Efield[i][j][k];
				tempEfield[i][j][k] = Efield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
				rotB[i][j][k] = Vector3d(0, 0, 0);
				Ederivative[i][j][k] = Vector3d(0, 0, 0);
				externalElectricFlux[i][j][k] = Vector3d(0, 0, 0);
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Bfield[i][j][k] = B0;
				//Bfield[i][j][k].x = B0.x;
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}

	createParticleTypes(concentrations, particlesPerBin);

	density = 0;
	for (int i = 0; i < typesNumber; ++i) {
		density += types[i].concentration * types[i].mass;
	}

	checkDebyeParameter();

	double concentration = density / (massProton + massElectron);

	omegaPlasmaProton = sqrt(
		4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massProton);
	omegaPlasmaElectron = sqrt(
		4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massElectron);

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");
	if (omegaPlasmaElectron * xsize / speed_of_light_normalized < 5) {
		if (rank == 0) printf("omegaPlasmaElectron*xsize/speed_of_light_normalized < 5\n");
		if (rank == 0) fprintf(informationFile, "omegaPlasmaElectron*xsize/speed_of_light_normalized < 5\n");
	}
	if (rank == 0)
		printf("omegaPlasmaElectron*xsize/speed_of_light_normalized = %g\n",
		       omegaPlasmaElectron * xsize / speed_of_light_normalized);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "omegaPlasmaElectron*xsize/speed_of_light_normalized = %g\n",
		        omegaPlasmaElectron * xsize / speed_of_light_normalized);

	if (omegaPlasmaElectron * deltaX / speed_of_light_normalized > 1) {
		if (rank == 0) printf("omegaPlasmaElectron*deltaX/speed_of_light_normalized > 1\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omegaPlasmaElectron*deltaX/speed_of_light_normalized > 1\n");
	}
	if (rank == 0)
		printf("omegaPlasmaElectron*deltaX/speed_of_light_normalized = %g\n",
		       omegaPlasmaElectron * deltaX / speed_of_light_normalized);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "omegaPlasmaElectron*deltaX/speed_of_light_normalized = %g\n",
		        omegaPlasmaElectron * deltaX / speed_of_light_normalized);
	if (rank == 0) fclose(informationFile);

	//checkDebyeParameter();
	checkGyroRadius();

	omegaGyroProton = electron_charge * B0.norm() / (massProton * speed_of_light);
	omegaGyroProton = electron_charge * B0.norm() / (massElectron * speed_of_light);

}

void Simulation::initializeSimpleElectroMagneticWave() {
	boundaryConditionType = PERIODIC;
	E0 = Vector3d(0, 0, 0);
	B0 = Vector3d(0, 0, 0);
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Bfield[i][j][k] = Vector3d(0, 0, 0);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
	double kw = (2 * pi / xsizeGeneral);
	printf("kw = %15.10g\n", kw);
	fflush(stdout);
	double E = 1E-5;

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Efield[i][j][k].x = 0;
				Efield[i][j][k].y = E * sin(kw * xgrid[i]);
				Efield[i][j][k].z = 0;
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = tempEfield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}

	if (nprocs == 1) {
		for (int k = 0; k < znumber; ++k) {
			for (int j = 0; j < ynumber; ++j) {
				Efield[xnumber][j][k] = Efield[1][j][k];
				tempEfield[xnumber][j][k] = Efield[1][j][k];
				newEfield[xnumber][j][k] = Efield[1][j][k];
				explicitEfield[xnumber][j][k] = explicitEfield[1][j][k];
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			Efield[i][j][znumber] = Efield[i][j][0];
			tempEfield[i][j][znumber] = Efield[i][j][0];
			newEfield[i][j][znumber] = Efield[i][j][0];
			explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
		}
	}


	for (int k = 0; k < znumber + 1; ++k) {
		for (int i = 0; i < xnumber + 1; ++i) {
			Efield[i][ynumber][k] = Efield[i][0][k];
			tempEfield[i][ynumber][k] = Efield[i][0][k];
			newEfield[i][ynumber][k] = Efield[i][0][k];
			explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Bfield[i][j][k].z = E * sin(kw * middleXgrid[i]);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}


	double t = 2 * pi / (kw * speed_of_light_normalized);
}

void Simulation::initializeRotatedSimpleElectroMagneticWave(int wavesCount) {
	boundaryConditionType = PERIODIC;

	Eyamplitude = 1;
	Ezamplitude = 0;
	Bzamplitude = Eyamplitude;
	Byamplitude = Ezamplitude;

	double kx = wavesCount * 2 * pi / xsizeGeneral;
	double ky = wavesCount * 2 * pi / ysizeGeneral;
	double kz = wavesCount * 2 * pi / zsizeGeneral;
	kz = 0;

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
				Efield[i][j][k].x = 0;
				Efield[i][j][k].y = Eyamplitude * cos(kx * xgrid[i] + ky * ygrid[j] + kz * zgrid[k]);
				Efield[i][j][k].z = Ezamplitude * sin(kx * xgrid[i] + ky * ygrid[j] + kz * zgrid[k]);
				Efield[i][j][k] = rotationMatrix * Efield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}

	for (int k = 0; k < znumber; ++k) {
		for (int j = 0; j < ynumber; ++j) {
			Efield[xnumber][j][k] = Efield[0][j][k];
			tempEfield[xnumber][j][k] = Efield[0][j][k];
			newEfield[xnumber][j][k] = Efield[0][j][k];
			explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			Efield[i][j][znumber] = Efield[i][j][0];
			tempEfield[i][j][znumber] = Efield[i][j][0];
			newEfield[i][j][znumber] = Efield[i][j][0];
			explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
		}
	}

	for (int k = 0; k < znumber + 1; ++k) {
		for (int i = 0; i < xnumber + 1; ++i) {
			Efield[i][ynumber][k] = Efield[i][0][k];
			tempEfield[i][ynumber][k] = Efield[i][0][k];
			newEfield[i][ynumber][k] = Efield[i][0][k];
			explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Bfield[i][j][k].x = 0;
				Bfield[i][j][k].y = Byamplitude * sin(kx * middleXgrid[i] + ky * middleYgrid[j] + kz * middleZgrid[k]);
				Bfield[i][j][k].z = Bzamplitude * cos(kx * middleXgrid[i] + ky * middleYgrid[j] + kz * middleZgrid[k]);
				Bfield[i][j][k] = rotationMatrix * Bfield[i][j][k];
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
}

void Simulation:: initializeAlfvenWaveX(int wavesCount, double amplitudeRelation) {
	boundaryConditionType = PERIODIC;
	if (rank == 0) printf("initialization alfven wave\n");
	fflush(stdout);


	double concentration = density / (massProton + massElectron);
	types[1].particesDeltaX = types[0].particesDeltaX;
	types[1].particlesPerBin = types[0].particlesPerBin;
	types[0].concentration = concentration;
	types[1].concentration = concentration;
	for (int i = 2; i < typesNumber; ++i) {
		types[i].particlesPerBin = 0;
		types[i].concentration = 0;
		types[i].particesDeltaX = xsize;
	}
	createParticles();
	E0 = Vector3d(0, 0, 0);

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");

	double alfvenV = B0.norm() / sqrt(4 * pi * density);
	if (alfvenV > speed_of_light_normalized) {
		printf("alfven velocity > c\n");
		if (rank == 0) fprintf(informationFile, "alfven velocity > c\n");
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "alfvenV/c = %15.10g > 1\n", alfvenV / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (rank == 0) fprintf(informationFile, "alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
	if (rank == 0) fprintf(informationFile, "alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);
	if (rank == 0) printf("alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
	if (rank == 0) printf("alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);
	fflush(stdout);

	double kw = (wavesCount * 2 * pi / xsizeGeneral);

	omegaPlasmaProton = sqrt(
		4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massProton);
	omegaPlasmaElectron = sqrt(
		4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
	omegaGyroProton = B0.norm() * electron_charge_normalized / (massProton * speed_of_light_normalized);
	omegaGyroElectron = B0.norm() * electron_charge_normalized / (massElectron * speed_of_light_normalized);

	if (omegaGyroProton < 5 * speed_of_light_normalized * kw) {
		if (rank == 0) printf("omegaGyroProton < 5*k*c\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omegaGyroProton < 5*k*c\n");
		//fclose(informationFile);
		//exit(0);
	}
	if (rank == 0) printf("omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));

	if (omegaPlasmaProton < 5 * omegaGyroProton) {
		if (rank == 0) printf("omegaPlasmaProton < 5*omegaGyroProton\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omegaPlasmaProton < 5*omegaGyroProton\n");
		//fclose(informationFile);
		//exit(0);
	}
	if (rank == 0) printf("omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);

	//w = q*kw*B/mP * 0.5*(sqrt(d)+-b)/a
	double b = speed_of_light_normalized * kw * (massProton - massElectron) / massProton;
	double discriminant = speed_of_light_normalized_sqr * kw * kw * sqr(
		massProton + massElectron) + 16 * pi * concentration * sqr(
		electron_charge_normalized) * (massProton + massElectron) / sqr(massProton);
	double a = (kw * kw * speed_of_light_normalized_sqr * massProton * massElectron + 4 * pi * concentration * sqr(
		electron_charge_normalized) * (massProton + massElectron)) / sqr(massProton);

	if (discriminant < 0) {
		printf("discriminant < 0\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "discriminant < 0\n");
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "discriminant = %15.10g\n", discriminant);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	double fakeOmega = (kw * electron_charge_normalized * B0.norm() / massProton) * (sqrt(
		discriminant) - b) / (2.0 * a);

	//a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0

	double a4 = sqr(speed_of_light_normalized_sqr * massProton * massElectron);
	a4 = a4 * sqr(sqr(sqr(fakeOmega)));
	double a3 = -2 * cube(speed_of_light_normalized_sqr) * sqr(kw * massElectron * massProton)
		- 8 * pi * concentration * sqr(
			speed_of_light_normalized_sqr * electron_charge_normalized) * massElectron * massProton * (massElectron + massProton)
		- sqr(B0.norm() * speed_of_light_normalized * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron));
	a3 = a3 * cube(sqr(fakeOmega));
	double a2 = sqr(sqr(speed_of_light_normalized_sqr * kw) * massProton * massElectron)
		+ 8 * pi * cube(speed_of_light_normalized_sqr) * concentration * sqr(
			kw * electron_charge_normalized) * massProton * massElectron * (massProton + massElectron)
		+ 16 * sqr(pi * speed_of_light_normalized_sqr * concentration * sqr(
			electron_charge_normalized) * (massProton + massElectron))
		+ 2 * sqr(
			B0.norm() * speed_of_light_normalized_sqr * kw * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron))
		+ 8 * pi * concentration * sqr(B0.norm() * speed_of_light_normalized * sqr(
			electron_charge_normalized)) * (massProton + massElectron)
		+ sqr(sqr(B0.norm() * electron_charge_normalized));
	a2 = a2 * sqr(sqr(fakeOmega));
	double a1 = -sqr(
			B0.norm() * cube(speed_of_light_normalized) * kw * kw * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron))
		- 8 * pi * concentration * sqr(B0.norm() * speed_of_light_normalized_sqr * kw * sqr(
			electron_charge_normalized)) * (massProton + massElectron)
		- 2 * sqr(speed_of_light_normalized * kw * sqr(B0.norm() * electron_charge_normalized));
	a1 = a1 * sqr(fakeOmega);
	double a0 = sqr(sqr(B0.norm() * speed_of_light_normalized * kw * electron_charge_normalized));

	a4 = a4 / a0;
	a3 = a3 / a0;
	a2 = a2 / a0;
	a1 = a1 / a0;
	a0 = 1.0;

	if (rank == 0) printf("a4 = %g\n", a4);
	if (rank == 0) fprintf(informationFile, "a4 = %g\n", a4);
	if (rank == 0) printf("a3 = %g\n", a3);
	if (rank == 0) fprintf(informationFile, "a3 = %g\n", a3);
	if (rank == 0) printf("a2 = %g\n", a2);
	if (rank == 0) fprintf(informationFile, "a2 = %g\n", a2);
	if (rank == 0) printf("a1 = %g\n", a1);
	if (rank == 0) fprintf(informationFile, "a1 = %g\n", a1);
	if (rank == 0) printf("a0 = %g\n", a0);
	if (rank == 0) fprintf(informationFile, "a0 = %g\n", a0);
	fflush(stdout);

	double fakeOmega1 = kw * alfvenV;
	if (rank == 0) printf("fakeOmega = %g\n", fakeOmega1 / plasma_period);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "fakeOmega = %g\n", fakeOmega1 / plasma_period);
	double realOmega2 = solve4orderEquation(a4, a3, a2, a1, a0, 1.0);
	if (realOmega2 < 0) {
		printf("omega^2 < 0\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega^2 < 0\n");
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "omega^2 = %15.10g > 1\n", realOmega2);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	double error = (((a4 * realOmega2 + a3) * realOmega2 + a2) * realOmega2 + a1) * realOmega2 + a0;
	if (rank == 0) printf("error = %15.10g\n", error);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "error = %15.10g\n", error);
	//double
	omega = sqrt(realOmega2) * fakeOmega;
	if (omega < 0) {
		omega = -omega;
	}

	if (omega > speed_of_light_normalized * kw / 5.0) {
		if (rank == 0) printf("omega > k*c/5\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega > k*c/5\n");
		if (rank == 0) printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
		//fclose(informationFile);
		if (rank == 0) errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		if (rank == 0) fprintf(errorLogFile, "omega/kc = %15.10g > 0.2\n", omega / (kw * speed_of_light_normalized));
		if (rank == 0) fclose(errorLogFile);
		//exit(0);
	}
	if (rank == 0) printf("omega = %g\n", omega / plasma_period);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "omega = %g\n", omega / plasma_period);

	if (rank == 0) printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));

	if (fabs(omega) > omegaGyroProton / 2) {
		if (rank == 0) printf("omega > omegaGyroProton/2\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega > omegaGyroProton/2\n");
	}
	if (rank == 0) printf("omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
	if (rank == 0) fclose(informationFile);

	checkFrequency(omega);

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");
	//checkCollisionTime(omega);
	//checkMagneticReynolds(alfvenV);
	//checkDissipation(kw, alfvenV);

	double epsilonAmplitude = amplitudeRelation;

	double alfvenVReal = omega / kw;

	//double
	Bzamplitude = B0.norm() * epsilonAmplitude;

	double Omegae = omegaGyroElectron;
	double Omegae2 = Omegae * Omegae;
	double Omegap = omegaGyroProton;
	double Omegap2 = Omegap * Omegap;
	double omegae = omegaPlasmaElectron;
	double omegae2 = omegae * omegae;
	double omegap = omegaPlasmaProton;
	double omegap2 = omegap * omegap;

	double kc = kw * speed_of_light_normalized;
	double kc2 = kc * kc;

	double denominator = omega * omega - Omegae2 - (omegae2 * omega * omega / (omega * omega - kc2));

	//double
	VzamplitudeProton = -((1.0 / (4 * pi * concentration * electron_charge_normalized)) * (kc + ((omegae2 + omegap2 - omega * omega) / kc) + (omegae2 * Omegae2 / (kc * denominator))) / ((Omegae * omegae2 * omega / ((kc2 - omega * omega) * denominator)) + (Omegap / omega))) * Bzamplitude;
	//double
	VzamplitudeElectron = (((electron_charge_normalized * omega * Omegae) / (massElectron * kc)) * Bzamplitude + (omegae2 * omega * omega / (kc2 - omega * omega)) * VzamplitudeProton) / denominator;

	//double
	Byamplitude = (4 * pi * concentration * electron_charge_normalized / ((omega * omega / kc) - kc)) * (VzamplitudeElectron - VzamplitudeProton);

	//double
	VyamplitudeProton = -(Omegap / omega) * VzamplitudeProton - (electron_charge_normalized / (massProton * kc)) * Bzamplitude;
	//double
	VyamplitudeElectron = (Omegae / omega) * VzamplitudeElectron + (electron_charge_normalized / (massElectron * kc)) * Bzamplitude;

	//double
	Eyamplitude = (omega / kc) * Bzamplitude;
	//double
	Ezamplitude = -(omega / kc) * Byamplitude;

	double xshift = 0;
	double phase = 0;
	//double xshift = xsize/4;

	//Eyamplitude = 0.0;
	//VzamplitudeElectron = 0.0;
	//VzamplitudeProton = 0.0;
	//Byamplitude = 0.0;

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k].x = 0;
				Efield[i][j][k].y = Eyamplitude * cos(kw * xgrid[i] - phase);
				Efield[i][j][k].z = Ezamplitude * sin(kw * xgrid[i] - phase);
				explicitEfield[i][j][k] = Efield[i][j][k];
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}

	if (nprocs == 1) {
		for (int k = 0; k < znumber; ++k) {
			for (int j = 0; j < ynumber; ++j) {
				Efield[xnumber][j][k] = Efield[1][j][k];
				tempEfield[xnumber][j][k] = Efield[1][j][k];
				newEfield[xnumber][j][k] = Efield[1][j][k];
				explicitEfield[xnumber][j][k] = explicitEfield[1][j][k];

				Efield[xnumber - 1][j][k] = Efield[0][j][k];
				tempEfield[xnumber - 1][j][k] = Efield[0][j][k];
				newEfield[xnumber - 1][j][k] = Efield[0][j][k];
				explicitEfield[xnumber - 1][j][k] = explicitEfield[0][j][k];
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			Efield[i][j][znumber] = Efield[i][j][0];
			tempEfield[i][j][znumber] = Efield[i][j][0];
			newEfield[i][j][znumber] = Efield[i][j][0];
			explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
		}
	}

	for (int k = 0; k < znumber + 1; ++k) {
		for (int i = 0; i < xnumber + 1; ++i) {
			Efield[i][ynumber][k] = Efield[i][0][k];
			tempEfield[i][ynumber][k] = Efield[i][0][k];
			newEfield[i][ynumber][k] = Efield[i][0][k];
			explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Bfield[i][j][k].x = B0.norm();
				Bfield[i][j][k].y = Byamplitude * sin(kw * middleXgrid[i] - phase);
				Bfield[i][j][k].z = Bzamplitude * cos(kw * middleXgrid[i] - phase);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}

	if (fabs(VzamplitudeProton) > speed_of_light_normalized) {
		printf("VzamplitudeProton > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VzamplitudeProton > speed_of_light_normalized\n");
		printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VzamplitudeProton/c = %15.10g > 1\n", VzamplitudeProton / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (rank == 0) printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);

	if (fabs(VzamplitudeElectron) > speed_of_light_normalized) {
		printf("VzamplitudeElectron > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VzamplitudeElectron > speed_of_light_normalized\n");
		printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VzamplitudeElectron/c = %15.10g > 1\n", VzamplitudeElectron / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (rank == 0) printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);

	if (fabs(VyamplitudeProton) > speed_of_light_normalized) {
		printf("VyamplitudeProton > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VyamplitudeProton > speed_of_light_normalized\n");
		printf("VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VyamplitudeProton/c = %15.10g > 1\n", VyamplitudeProton / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	if (fabs(VyamplitudeElectron) > speed_of_light_normalized) {
		printf("VyamplitudeElectron > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VyamplitudeElectron > speed_of_light_normalized\n");
		printf("VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VyamplitudeElectron/c = %15.10g > 1\n", VyamplitudeElectron / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	//k > 0, w > 0, kx-wt, Bz > 0, Vz < 0, Vyp > 0, Vye < 0

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		Vector3d velocity = particle->getVelocity(speed_of_light_normalized);
		int xn = (particle->coordinates.x - xgrid[0]) / deltaX;
		double rightWeight = (particle->coordinates.x - xgrid[xn]) / deltaX;
		double leftWeight = (xgrid[xn + 1] - particle->coordinates.x) / deltaX;
		if (particle->type == PROTON) {
			velocity = (Vector3d(0, 0, 1) * (VzamplitudeProton) * cos(kw * particle->coordinates.x - phase) + Vector3d(
				0, 1, 0) * VyamplitudeProton * sin(kw * particle->coordinates.x - phase));
			//velocity = Vector3d(0, 0, 1) * (VzamplitudeProton) * (leftWeight*(cos(kw * (xgrid[xn] - xshift)) + rightWeight*cos(kw*(xgrid[xn+1] - xshift)))) + Vector3d(0, 1, 0) * VyamplitudeProton * (leftWeight*(sin(kw * (xgrid[xn] - xshift)) + rightWeight*sin(kw*(xgrid[xn+1] - xshift))));

		}
		if (particle->type == ELECTRON) {
			velocity = (Vector3d(0, 0, 1) * (VzamplitudeElectron) * cos(
				kw * particle->coordinates.x - phase) + Vector3d(0, 1, 0) * VyamplitudeElectron * sin(
				kw * particle->coordinates.x - phase));
			//velocity = Vector3d(0, 0, 1) * (VzamplitudeElectron) * (leftWeight*(cos(kw * (xgrid[xn] - xshift)) + rightWeight*cos(kw*(xgrid[xn+1] - xshift)))) + Vector3d(0, 1, 0) * VyamplitudeElectron * (leftWeight*(sin(kw * (xgrid[xn] - xshift)) + rightWeight*sin(kw*(xgrid[xn+1] - xshift))));
		}
		double beta = velocity.norm() / speed_of_light_normalized;
		particle->addVelocity(velocity, speed_of_light_normalized);
		Vector3d momentum = particle->getMomentum();
		particle->initialMomentum = momentum;
		//particle->prevMomentum = momentum;
	}

	updateDeltaT();

	if (rank == 0) printf("dt/Talfven = %g\n", deltaT * omega / (2 * pi));
	if (rank == 0) printf("dt = %g\n", deltaT * plasma_period);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "dt/Talfven = %g\n", deltaT * omega / (2 * pi));
	if (rank == 0) fprintf(informationFile, "dt = %g\n", deltaT * plasma_period);

	double Vthermal = sqrt(2 * kBoltzman_normalized * temperature / massElectron);
	double thermalFlux = Vthermal * concentration * electron_charge_normalized / sqrt(1.0 * types[0].particlesPerBin);
	double alfvenFlux = (VyamplitudeProton - VyamplitudeElectron) * concentration * electron_charge_normalized;
	if (thermalFlux > alfvenFlux / 2) {
		if (rank == 0) printf("thermalFlux > alfvenFlux/2\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "thermalFlux > alfvenFlux/2\n");
	}
	if (rank == 0) printf("alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
	double minDeltaT = deltaX / Vthermal;
	if (minDeltaT > deltaT) {
		if (rank == 0) printf("deltaT < dx/Vthermal\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "deltaT < dx/Vthermal\n");

		//printf("deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);
		//fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);

		//fclose(informationFile);
		//exit(0);
	}
	if (rank == 0) {
		printf("deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
		fflush(stdout);
		fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
		fprintf(informationFile, "\n");

		fprintf(informationFile, "Bz amplitude = %g\n", Bzamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "By amplitude = %g\n", Byamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "Vz amplitude p = %g\n", VzamplitudeProton * scaleFactor / plasma_period);
		fprintf(informationFile, "Vz amplitude e = %g\n", VzamplitudeElectron * scaleFactor / plasma_period);
		fprintf(informationFile, "Vy amplitude p = %g\n", VyamplitudeProton * scaleFactor / plasma_period);
		fprintf(informationFile, "Vy amplitude e = %g\n", VyamplitudeElectron * scaleFactor / plasma_period);
		fprintf(informationFile, "Ey amplitude = %g\n", Eyamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "Ez amplitude = %g\n", Ezamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "By/Ez = %g\n", Byamplitude / Ezamplitude);
		fprintf(informationFile, "Bz/Ey = %g\n", Bzamplitude / Eyamplitude);
		fprintf(informationFile, "4*pi*Jy amplitude = %g\n",
		        4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "c*rotBy amplitude = %g\n",
		        speed_of_light_normalized * kw * Bzamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "4*pi*Jz amplitude = %g\n",
		        4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "c*rotBz amplitude = %g\n",
		        speed_of_light_normalized * kw * Byamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative By amplitude = %g\n",
		        -omega * Byamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "-c*rotEy = %g\n",
		        speed_of_light_normalized * kw * Ezamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Bz amplitude = %g\n",
		        omega * Bzamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "-c*rotEz = %g\n",
		        speed_of_light_normalized * kw * Eyamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Ey amplitude = %g\n",
		        omega * Eyamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "c*rotBy - 4*pi*Jy = %g\n",
		        (speed_of_light_normalized * kw * Bzamplitude - 4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Ez amplitude = %g\n",
		        -omega * Ezamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "c*rotBz - 4*pi*Jz = %g\n",
		        (speed_of_light_normalized * kw * Byamplitude - 4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativJy = -electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jy amplitude = %g\n",
		        derivativJy / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesY = electron_charge_normalized * ((Eyamplitude * ((1.0 / massProton) + (1.0 / massElectron))) + B0.norm() * ((VzamplitudeProton / massProton) + (VzamplitudeElectron / massElectron)) / speed_of_light_normalized);
		fprintf(informationFile, "dJy/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesY / (plasma_period * plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		double derivativJz = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jz amplitude = %g\n",
		        derivativJz / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesZ = electron_charge_normalized * ((Ezamplitude * ((1.0 / massProton) + (1.0 / massElectron))) - B0.norm() * ((VyamplitudeProton / massProton) + (VyamplitudeElectron / massElectron)) / speed_of_light_normalized);
		fprintf(informationFile, "dJz/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesZ / (plasma_period * plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativVyp = -omega * VyamplitudeProton;
		fprintf(informationFile, "-w*Vyp amplitude = %g\n", derivativVyp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonY = electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeProton / speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVyp/dt amplitude = %g\n",
		        derivativeVelocityProtonY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVzp = omega * VzamplitudeProton;
		fprintf(informationFile, "w*Vzp amplitude = %g\n", derivativVzp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonZ = electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeProton / speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVzp/dt amplitude = %g\n",
		        derivativeVelocityProtonZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVye = -omega * VyamplitudeElectron;
		fprintf(informationFile, "-w*Vye amplitude = %g\n", derivativVye * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronY = -electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeElectron / speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVye/dt amplitude = %g\n",
		        derivativeVelocityElectronY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVze = omega * VzamplitudeElectron;
		fprintf(informationFile, "w*Vze amplitude = %g\n", derivativVze * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronZ = -electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeElectron / speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVze/dt amplitude = %g\n",
		        derivativeVelocityElectronZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");
	}

	if (rank == 0) fclose(informationFile);
}


void Simulation::initializeAlfvenWaveY(int wavesCount, double amplitudeRelation) {
	boundaryConditionType = PERIODIC;
	if (rank == 0)printf("initialization alfven wave\n");
	fflush(stdout);

	double concentration = density / (massProton + massElectron);
	types[1].particesDeltaX = types[0].particesDeltaX;
	types[1].particlesPerBin = types[0].particlesPerBin;
	types[0].concentration = concentration;
	types[1].concentration = concentration;
	for (int i = 2; i < typesNumber; ++i) {
		types[i].particlesPerBin = 0;
		types[i].concentration = 0;
		types[i].particesDeltaX = xsize;
	}
	createParticles();
	E0 = Vector3d(0, 0, 0);

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");

	double alfvenV = B0.norm() / sqrt(4 * pi * density);
	if (alfvenV > speed_of_light_normalized) {
		printf("alfven velocity > c\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "alfven velocity > c\n");
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "alfvenV/c = %15.10g > 1\n", alfvenV / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (rank == 0) fprintf(informationFile, "alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
	if (rank == 0) fprintf(informationFile, "alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);
	if (rank == 0) printf("alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
	if (rank == 0) printf("alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);
	fflush(stdout);

	//double kw = wavesCount * 2 * pi / xsize;
	double kw = wavesCount * 2 * pi / ysizeGeneral;

	omegaPlasmaProton = sqrt(
		4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massProton);
	omegaPlasmaElectron = sqrt(
		4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
	omegaGyroProton = B0.norm() * electron_charge_normalized / (massProton * speed_of_light_normalized);
	omegaGyroElectron = B0.norm() * electron_charge_normalized / (massElectron * speed_of_light_normalized);

	if (omegaGyroProton < 5 * speed_of_light_normalized * kw) {
		if (rank == 0) printf("omegaGyroProton < 5*k*c\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omegaGyroProton < 5*k*c\n");
		//fclose(informationFile);
		//exit(0);
	}
	if (rank == 0) printf("omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));

	if (omegaPlasmaProton < 5 * omegaGyroProton) {
		if (rank == 0) printf("omegaPlasmaProton < 5*omegaGyroProton\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omegaPlasmaProton < 5*omegaGyroProton\n");
		//fclose(informationFile);
		//exit(0);
	}
	if (rank == 0) printf("omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);

	//w = q*kw*B/mP * 0.5*(sqrt(d)+-b)/a
	double b = speed_of_light_normalized * kw * (massProton - massElectron) / massProton;
	double discriminant = speed_of_light_normalized_sqr * kw * kw * sqr(
		massProton + massElectron) + 16 * pi * concentration * sqr(
		electron_charge_normalized) * (massProton + massElectron) / sqr(massProton);
	double a = (kw * kw * speed_of_light_normalized_sqr * massProton * massElectron + 4 * pi * concentration * sqr(
		electron_charge_normalized) * (massProton + massElectron)) / sqr(massProton);

	if (discriminant < 0) {
		printf("discriminant < 0\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "discriminant < 0\n");
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "discriminant = %15.10g\n", discriminant);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	double fakeOmega = (kw * electron_charge_normalized * B0.norm() / massProton) * (sqrt(
		discriminant) - b) / (2.0 * a);

	//a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0

	double a4 = sqr(speed_of_light_normalized_sqr * massProton * massElectron);
	a4 = a4 * sqr(sqr(sqr(fakeOmega)));
	double a3 = -2 * cube(speed_of_light_normalized_sqr) * sqr(kw * massElectron * massProton)
		- 8 * pi * concentration * sqr(
			speed_of_light_normalized_sqr * electron_charge_normalized) * massElectron * massProton * (massElectron + massProton)
		- sqr(B0.norm() * speed_of_light_normalized * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron));
	a3 = a3 * cube(sqr(fakeOmega));
	double a2 = sqr(sqr(speed_of_light_normalized_sqr * kw) * massProton * massElectron)
		+ 8 * pi * cube(speed_of_light_normalized_sqr) * concentration * sqr(
			kw * electron_charge_normalized) * massProton * massElectron * (massProton + massElectron)
		+ 16 * sqr(pi * speed_of_light_normalized_sqr * concentration * sqr(
			electron_charge_normalized) * (massProton + massElectron))
		+ 2 * sqr(
			B0.norm() * speed_of_light_normalized_sqr * kw * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron))
		+ 8 * pi * concentration * sqr(B0.norm() * speed_of_light_normalized * sqr(
			electron_charge_normalized)) * (massProton + massElectron)
		+ sqr(sqr(B0.norm() * electron_charge_normalized));
	a2 = a2 * sqr(sqr(fakeOmega));
	double a1 = -sqr(
			B0.norm() * cube(speed_of_light_normalized) * kw * kw * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron))
		- 8 * pi * concentration * sqr(B0.norm() * speed_of_light_normalized_sqr * kw * sqr(
			electron_charge_normalized)) * (massProton + massElectron)
		- 2 * sqr(speed_of_light_normalized * kw * sqr(B0.norm() * electron_charge_normalized));
	a1 = a1 * sqr(fakeOmega);
	double a0 = sqr(sqr(B0.norm() * speed_of_light_normalized * kw * electron_charge_normalized));

	a4 = a4 / a0;
	a3 = a3 / a0;
	a2 = a2 / a0;
	a1 = a1 / a0;
	a0 = 1.0;

	if (rank == 0) printf("a4 = %g\n", a4);
	if (rank == 0) fprintf(informationFile, "a4 = %g\n", a4);
	if (rank == 0) printf("a3 = %g\n", a3);
	if (rank == 0) fprintf(informationFile, "a3 = %g\n", a3);
	if (rank == 0) printf("a2 = %g\n", a2);
	if (rank == 0) fprintf(informationFile, "a2 = %g\n", a2);
	if (rank == 0) printf("a1 = %g\n", a1);
	if (rank == 0) fprintf(informationFile, "a1 = %g\n", a1);
	if (rank == 0) printf("a0 = %g\n", a0);
	if (rank == 0) fprintf(informationFile, "a0 = %g\n", a0);
	fflush(stdout);

	double fakeOmega1 = kw * alfvenV;
	if (rank == 0) printf("fakeOmega = %g\n", fakeOmega1 / plasma_period);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "fakeOmega = %g\n", fakeOmega1 / plasma_period);
	double realOmega2 = solve4orderEquation(a4, a3, a2, a1, a0, 1.0);
	if (realOmega2 < 0) {
		printf("omega^2 < 0\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega^2 < 0\n");
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "omega^2 = %15.10g > 1\n", realOmega2);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	double error = (((a4 * realOmega2 + a3) * realOmega2 + a2) * realOmega2 + a1) * realOmega2 + a0;
	if (rank == 0) printf("error = %15.10g\n", error);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "error = %15.10g\n", error);
	//double
	omega = sqrt(realOmega2) * fakeOmega;
	if (omega < 0) {
		omega = -omega;
	}

	if (omega > speed_of_light_normalized * kw / 5.0) {
		printf("omega > k*c/5\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega > k*c/5\n");
		if (rank == 0) printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
		//fclose(informationFile);
		if (rank == 0) errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		if (rank == 0) fprintf(errorLogFile, "omega/kc = %15.10g > 0.2\n", omega / (kw * speed_of_light_normalized));
		if (rank == 0) fclose(errorLogFile);
		//exit(0);
	}
	if (rank == 0) printf("omega = %g\n", omega / plasma_period);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "omega = %g\n", omega / plasma_period);

	if (rank == 0) printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));

	if (fabs(omega) > omegaGyroProton / 2) {
		if (rank == 0) printf("omega > omegaGyroProton/2\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega > omegaGyroProton/2\n");
	}
	if (rank == 0) printf("omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
	if (rank == 0) fclose(informationFile);

	checkFrequency(omega);

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");
	//checkCollisionTime(omega);
	//checkMagneticReynolds(alfvenV);
	//checkDissipation(kw, alfvenV);

	double epsilonAmplitude = amplitudeRelation;

	double alfvenVReal = omega / kw;

	//double
	Bzamplitude = B0.norm() * epsilonAmplitude;

	double Omegae = omegaGyroElectron;
	double Omegae2 = Omegae * Omegae;
	double Omegap = omegaGyroProton;
	double Omegap2 = Omegap * Omegap;
	double omegae = omegaPlasmaElectron;
	double omegae2 = omegae * omegae;
	double omegap = omegaPlasmaProton;
	double omegap2 = omegap * omegap;

	double kc = kw * speed_of_light_normalized;
	double kc2 = kc * kc;

	double denominator = omega * omega - Omegae2 - (omegae2 * omega * omega / (omega * omega - kc2));

	//double
	VzamplitudeProton = -((1.0 / (4 * pi * concentration * electron_charge_normalized)) * (kc + ((omegae2 + omegap2 - omega * omega) / kc) + (omegae2 * Omegae2 / (kc * denominator))) / ((Omegae * omegae2 * omega / ((kc2 - omega * omega) * denominator)) + (Omegap / omega))) * Bzamplitude;
	//double
	VzamplitudeElectron = (((electron_charge_normalized * omega * Omegae) / (massElectron * kc)) * Bzamplitude + (omegae2 * omega * omega / (kc2 - omega * omega)) * VzamplitudeProton) / denominator;

	//double
	Byamplitude = (4 * pi * concentration * electron_charge_normalized / ((omega * omega / kc) - kc)) * (VzamplitudeElectron - VzamplitudeProton);

	//double
	VyamplitudeProton = -(Omegap / omega) * VzamplitudeProton - (electron_charge_normalized / (massProton * kc)) * Bzamplitude;
	//double
	VyamplitudeElectron = (Omegae / omega) * VzamplitudeElectron + (electron_charge_normalized / (massElectron * kc)) * Bzamplitude;

	//double
	Eyamplitude = (omega / kc) * Bzamplitude;
	//double
	Ezamplitude = -(omega / kc) * Byamplitude;

	double xshift = 0;
	//double xshift = xsize/4;

	//Eyamplitude = 0.0;
	//VzamplitudeElectron = 0.0;
	//VzamplitudeProton = 0.0;
	//Byamplitude = 0.0;

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				/*Efield[i][j][k].x = 0;
				Efield[i][j][k].y = Eyamplitude * cos(kw * xgrid[i] - kw * xshift);
				Efield[i][j][k].z = Ezamplitude * sin(kw * xgrid[i] - kw * xshift);*/
				Efield[i][j][k].x = -Eyamplitude * cos(kw * ygrid[j] - kw * xshift);
				Efield[i][j][k].y = 0;
				Efield[i][j][k].z = Ezamplitude * sin(kw * ygrid[j] - kw * xshift);
				explicitEfield[i][j][k] = Efield[i][j][k];
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}

	if (nprocs == 1) {
		for (int k = 0; k < znumber; ++k) {
			for (int j = 0; j < ynumber; ++j) {
				Efield[xnumber][j][k] = Efield[0][j][k];
				tempEfield[xnumber][j][k] = Efield[0][j][k];
				newEfield[xnumber][j][k] = Efield[0][j][k];
				explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			Efield[i][j][znumber] = Efield[i][j][0];
			tempEfield[i][j][znumber] = Efield[i][j][0];
			newEfield[i][j][znumber] = Efield[i][j][0];
			explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
		}
	}

	for (int k = 0; k < znumber + 1; ++k) {
		for (int i = 0; i < xnumber + 1; ++i) {
			Efield[i][ynumber][k] = Efield[i][0][k];
			tempEfield[i][ynumber][k] = Efield[i][0][k];
			newEfield[i][ynumber][k] = Efield[i][0][k];
			explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				/*Bfield[i][j][k].x = B0.norm();
				Bfield[i][j][k].y = Byamplitude * sin(kw * middleXgrid[i] - kw * xshift);
				Bfield[i][j][k].z = Bzamplitude * cos(kw * middleXgrid[i] - kw * xshift);*/
				Bfield[i][j][k].x = -Byamplitude * sin(kw * middleYgrid[j] - kw * xshift);
				Bfield[i][j][k].y = B0.norm();
				Bfield[i][j][k].z = Bzamplitude * cos(kw * middleYgrid[j] - kw * xshift);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}

	if (fabs(VzamplitudeProton) > speed_of_light_normalized) {
		printf("VzamplitudeProton > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VzamplitudeProton > speed_of_light_normalized\n");
		printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VzamplitudeProton/c = %15.10g > 1\n", VzamplitudeProton / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (rank == 0) printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);

	if (fabs(VzamplitudeElectron) > speed_of_light_normalized) {
		printf("VzamplitudeElectron > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VzamplitudeElectron > speed_of_light_normalized\n");
		printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VzamplitudeElectron/c = %15.10g > 1\n", VzamplitudeElectron / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (rank == 0) printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);

	if (fabs(VyamplitudeProton) > speed_of_light_normalized) {
		printf("VyamplitudeProton > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VyamplitudeProton > speed_of_light_normalized\n");
		printf("VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VyamplitudeProton/c = %15.10g > 1\n", VyamplitudeProton / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	if (fabs(VyamplitudeElectron) > speed_of_light_normalized) {
		printf("VyamplitudeElectron > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VyamplitudeElectron > speed_of_light_normalized\n");
		printf("VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VyamplitudeElectron/c = %15.10g > 1\n", VyamplitudeElectron / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	//k > 0, w > 0, kx-wt, Bz > 0, Vz < 0, Vyp > 0, Vye < 0

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		Vector3d velocity = particle->getVelocity(speed_of_light_normalized);
		int yn = (particle->coordinates.y - ygrid[0]) / deltaY;
		double rightWeight = (particle->coordinates.y - ygrid[yn]) / deltaY;
		double leftWeight = (ygrid[yn + 1] - particle->coordinates.y) / deltaY;
		if (particle->type == PROTON) {
			//velocity = Vector3d(0, 0, 1) * (VzamplitudeProton) * (leftWeight*cos(kw * (ygrid[yn] - xshift)) + rightWeight*cos(kw*(ygrid[yn+1] - xshift))) + Vector3d(0, 1, 0) * VyamplitudeProton * (leftWeight*sin(kw * (ygrid[yn] - xshift)) + rightWeight*sin(kw*(ygrid[yn+1] - xshift)));
			velocity = (Vector3d(0, 0, 1) * (VzamplitudeProton) * cos(
				kw * particle->coordinates.y - kw * xshift) - Vector3d(1, 0, 0) * VyamplitudeProton * sin(
				kw * particle->coordinates.y - kw * xshift));

		}
		if (particle->type == ELECTRON) {

			//velocity = Vector3d(0, 0, 1) * (VzamplitudeElectron) * (leftWeight*cos(kw * (ygrid[yn] - xshift)) + rightWeight*cos(kw*(ygrid[yn+1] - xshift))) + Vector3d(0, 1, 0) * VyamplitudeElectron * (leftWeight*sin(kw * (ygrid[yn] - xshift)) + rightWeight*sin(kw*(ygrid[yn+1] - xshift)));
			velocity = (Vector3d(0, 0, 1) * (VzamplitudeElectron) * cos(
				kw * particle->coordinates.y - kw * xshift) - Vector3d(1, 0, 0) * VyamplitudeElectron * sin(
				kw * particle->coordinates.y - kw * xshift));
		}
		double beta = velocity.norm() / speed_of_light_normalized;
		particle->addVelocity(velocity, speed_of_light_normalized);
		Vector3d momentum = particle->getMomentum();
		particle->initialMomentum = momentum;
		//particle->prevMomentum = momentum;
	}

	updateDeltaT();

	if (rank == 0) printf("dt/Talfven = %g\n", deltaT * omega / (2 * pi));
	if (rank == 0) printf("dt = %g\n", deltaT * plasma_period);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "dt/Talfven = %g\n", deltaT * omega / (2 * pi));
	if (rank == 0) fprintf(informationFile, "dt = %g\n", deltaT * plasma_period);

	double Vthermal = sqrt(2 * kBoltzman_normalized * temperature / massElectron);
	double thermalFlux = Vthermal * concentration * electron_charge_normalized / sqrt(1.0 * types[0].particlesPerBin);
	double alfvenFlux = (VyamplitudeProton - VyamplitudeElectron) * concentration * electron_charge_normalized;
	if (thermalFlux > alfvenFlux / 2) {
		if (rank == 0) printf("thermalFlux > alfvenFlux/2\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "thermalFlux > alfvenFlux/2\n");
	}
	if (rank == 0) printf("alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
	double minDeltaT = deltaX / Vthermal;
	if (minDeltaT > deltaT) {
		if (rank == 0) printf("deltaT < dx/Vthermal\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "deltaT < dx/Vthermal\n");

		//printf("deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);
		//fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);

		//fclose(informationFile);
		//exit(0);
	}
	if (rank == 0) {
		printf("deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
		fflush(stdout);
		fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
		fprintf(informationFile, "\n");

		fprintf(informationFile, "Bz amplitude = %g\n", Bzamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "By amplitude = %g\n", Byamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "Vz amplitude p = %g\n", VzamplitudeProton * scaleFactor / plasma_period);
		fprintf(informationFile, "Vz amplitude e = %g\n", VzamplitudeElectron * scaleFactor / plasma_period);
		fprintf(informationFile, "Vy amplitude p = %g\n", VyamplitudeProton * scaleFactor / plasma_period);
		fprintf(informationFile, "Vy amplitude e = %g\n", VyamplitudeElectron * scaleFactor / plasma_period);
		fprintf(informationFile, "Ey amplitude = %g\n", Eyamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "Ez amplitude = %g\n", Ezamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "By/Ez = %g\n", Byamplitude / Ezamplitude);
		fprintf(informationFile, "Bz/Ey = %g\n", Bzamplitude / Eyamplitude);
		fprintf(informationFile, "4*pi*Jy amplitude = %g\n",
		        4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "c*rotBy amplitude = %g\n",
		        speed_of_light_normalized * kw * Bzamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "4*pi*Jz amplitude = %g\n",
		        4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "c*rotBz amplitude = %g\n",
		        speed_of_light_normalized * kw * Byamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative By amplitude = %g\n",
		        -omega * Byamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "-c*rotEy = %g\n",
		        speed_of_light_normalized * kw * Ezamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Bz amplitude = %g\n",
		        omega * Bzamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "-c*rotEz = %g\n",
		        speed_of_light_normalized * kw * Eyamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Ey amplitude = %g\n",
		        omega * Eyamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "c*rotBy - 4*pi*Jy = %g\n",
		        (speed_of_light_normalized * kw * Bzamplitude - 4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Ez amplitude = %g\n",
		        -omega * Ezamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "c*rotBz - 4*pi*Jz = %g\n",
		        (speed_of_light_normalized * kw * Byamplitude - 4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativJy = -electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jy amplitude = %g\n",
		        derivativJy / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesY = electron_charge_normalized * ((Eyamplitude * ((1.0 / massProton) + (1.0 / massElectron))) + B0.norm() * ((VzamplitudeProton / massProton) + (VzamplitudeElectron / massElectron)) / speed_of_light_normalized);
		fprintf(informationFile, "dJy/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesY / (plasma_period * plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		double derivativJz = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jz amplitude = %g\n",
		        derivativJz / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesZ = electron_charge_normalized * ((Ezamplitude * ((1.0 / massProton) + (1.0 / massElectron))) - B0.norm() * ((VyamplitudeProton / massProton) + (VyamplitudeElectron / massElectron)) / speed_of_light_normalized);
		fprintf(informationFile, "dJz/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesZ / (plasma_period * plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativVyp = -omega * VyamplitudeProton;
		fprintf(informationFile, "-w*Vyp amplitude = %g\n", derivativVyp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonY = electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeProton / speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVyp/dt amplitude = %g\n",
		        derivativeVelocityProtonY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVzp = omega * VzamplitudeProton;
		fprintf(informationFile, "w*Vzp amplitude = %g\n", derivativVzp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonZ = electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeProton / speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVzp/dt amplitude = %g\n",
		        derivativeVelocityProtonZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVye = -omega * VyamplitudeElectron;
		fprintf(informationFile, "-w*Vye amplitude = %g\n", derivativVye * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronY = -electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeElectron / speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVye/dt amplitude = %g\n",
		        derivativeVelocityElectronY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVze = omega * VzamplitudeElectron;
		fprintf(informationFile, "w*Vze amplitude = %g\n", derivativVze * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronZ = -electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeElectron / speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVze/dt amplitude = %g\n",
		        derivativeVelocityElectronZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");
	}

	if (rank == 0) fclose(informationFile);
}

void Simulation::initializeRotatedAlfvenWave(int wavesCount, double amplitudeRelation) {
	boundaryConditionType = PERIODIC;
	if (rank == 0) printf("initialization alfven wave\n");
	fflush(stdout);

	double concentration = density / (massProton + massElectron);
	types[1].particesDeltaX = types[0].particesDeltaX;
	types[1].particlesPerBin = types[0].particlesPerBin;
	types[0].concentration = concentration;
	types[1].concentration = concentration;
	for (int i = 2; i < typesNumber; ++i) {
		types[i].particlesPerBin = 0;
		types[i].concentration = 0;
		types[i].particesDeltaX = xsize;
	}
	createParticles();
	E0 = Vector3d(0, 0, 0);

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");

	double alfvenV = B0.norm() / sqrt(4 * pi * density);
	if (alfvenV > speed_of_light_normalized) {
		printf("alfven velocity > c\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "alfven velocity > c\n");
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "alfvenV/c = %15.10g > 1\n", alfvenV / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (rank == 0) fprintf(informationFile, "alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
	if (rank == 0) fprintf(informationFile, "alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);
	if (rank == 0) printf("alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
	if (rank == 0) printf("alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);
	fflush(stdout);

	double kx = wavesCount * 2 * pi / xsize;
	double ky = wavesCount * 2 * pi / ysize;
	double kz = wavesCount * 2 * pi / zsize;
	kz = 0;

	double kw = sqrt(kx * kx + ky * ky + kz * kz);

	double weight = concentration * volumeB(0, 0, 0) / types[0].particlesPerBin;

	omegaPlasmaProton = sqrt(
		4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massProton);
	omegaPlasmaElectron = sqrt(
		4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
	omegaGyroProton = B0.norm() * electron_charge_normalized / (massProton * speed_of_light_normalized);
	omegaGyroElectron = B0.norm() * electron_charge_normalized / (massElectron * speed_of_light_normalized);

	if (omegaGyroProton < 5 * speed_of_light_normalized * kw) {
		if (rank == 0) printf("omegaGyroProton < 5*k*c\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omegaGyroProton < 5*k*c\n");
		//fclose(informationFile);
		//exit(0);
	}
	if (rank == 0) printf("omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));

	if (omegaPlasmaProton < 5 * omegaGyroProton) {
		if (rank == 0) printf("omegaPlasmaProton < 5*omegaGyroProton\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omegaPlasmaProton < 5*omegaGyroProton\n");
		//fclose(informationFile);
		//exit(0);
	}
	if (rank == 0) printf("omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);

	//w = q*kw*B/mP * 0.5*(sqrt(d)+-b)/a
	double b = speed_of_light_normalized * kw * (massProton - massElectron) / massProton;
	double discriminant = speed_of_light_normalized_sqr * kw * kw * sqr(
		massProton + massElectron) + 16 * pi * concentration * sqr(
		electron_charge_normalized) * (massProton + massElectron) / sqr(massProton);
	double a = (kw * kw * speed_of_light_normalized_sqr * massProton * massElectron + 4 * pi * concentration * sqr(
		electron_charge_normalized) * (massProton + massElectron)) / sqr(massProton);

	if (discriminant < 0) {
		printf("discriminant < 0\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "discriminant < 0\n");
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "discriminant = %15.10g\n", discriminant);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	double fakeOmega = (kw * electron_charge_normalized * B0.norm() / massProton) * (sqrt(
		discriminant) - b) / (2.0 * a);

	//a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0

	double a4 = sqr(speed_of_light_normalized_sqr * massProton * massElectron);
	a4 = a4 * sqr(sqr(sqr(fakeOmega)));
	double a3 = -2 * cube(speed_of_light_normalized_sqr) * sqr(kw * massElectron * massProton)
		- 8 * pi * concentration * sqr(
			speed_of_light_normalized_sqr * electron_charge_normalized) * massElectron * massProton * (massElectron + massProton)
		- sqr(B0.norm() * speed_of_light_normalized * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron));
	a3 = a3 * cube(sqr(fakeOmega));
	double a2 = sqr(sqr(speed_of_light_normalized_sqr * kw) * massProton * massElectron)
		+ 8 * pi * cube(speed_of_light_normalized_sqr) * concentration * sqr(
			kw * electron_charge_normalized) * massProton * massElectron * (massProton + massElectron)
		+ 16 * sqr(pi * speed_of_light_normalized_sqr * concentration * sqr(
			electron_charge_normalized) * (massProton + massElectron))
		+ 2 * sqr(
			B0.norm() * speed_of_light_normalized_sqr * kw * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron))
		+ 8 * pi * concentration * sqr(B0.norm() * speed_of_light_normalized * sqr(
			electron_charge_normalized)) * (massProton + massElectron)
		+ sqr(sqr(B0.norm() * electron_charge_normalized));
	a2 = a2 * sqr(sqr(fakeOmega));
	double a1 = -sqr(
			B0.norm() * cube(speed_of_light_normalized) * kw * kw * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron))
		- 8 * pi * concentration * sqr(B0.norm() * speed_of_light_normalized_sqr * kw * sqr(
			electron_charge_normalized)) * (massProton + massElectron)
		- 2 * sqr(speed_of_light_normalized * kw * sqr(B0.norm() * electron_charge_normalized));
	a1 = a1 * sqr(fakeOmega);
	double a0 = sqr(sqr(B0.norm() * speed_of_light_normalized * kw * electron_charge_normalized));

	a4 = a4 / a0;
	a3 = a3 / a0;
	a2 = a2 / a0;
	a1 = a1 / a0;
	a0 = 1.0;

	if (rank == 0) printf("a4 = %g\n", a4);
	if (rank == 0) fprintf(informationFile, "a4 = %g\n", a4);
	if (rank == 0) printf("a3 = %g\n", a3);
	if (rank == 0) fprintf(informationFile, "a3 = %g\n", a3);
	if (rank == 0) printf("a2 = %g\n", a2);
	if (rank == 0) fprintf(informationFile, "a2 = %g\n", a2);
	if (rank == 0) printf("a1 = %g\n", a1);
	if (rank == 0) fprintf(informationFile, "a1 = %g\n", a1);
	if (rank == 0) printf("a0 = %g\n", a0);
	if (rank == 0) fprintf(informationFile, "a0 = %g\n", a0);
	fflush(stdout);

	double fakeOmega1 = kw * alfvenV;
	if (rank == 0) printf("fakeOmega = %g\n", fakeOmega1 / plasma_period);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "fakeOmega = %g\n", fakeOmega1 / plasma_period);
	double realOmega2 = solve4orderEquation(a4, a3, a2, a1, a0, 1.0);
	if (realOmega2 < 0) {
		printf("omega^2 < 0\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega^2 < 0\n");
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "omega^2 = %15.10g > 1\n", realOmega2);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	double error = (((a4 * realOmega2 + a3) * realOmega2 + a2) * realOmega2 + a1) * realOmega2 + a0;
	if (rank == 0) printf("error = %15.10g\n", error);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "error = %15.10g\n", error);
	//double
	omega = sqrt(realOmega2) * fakeOmega;
	if (omega < 0) {
		omega = -omega;
	}

	if (omega > speed_of_light_normalized * kw / 5.0) {
		if (rank == 0) printf("omega > k*c/5\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega > k*c/5\n");
		if (rank == 0) printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
		//fclose(informationFile);
		if (rank == 0) errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		if (rank == 0) fprintf(errorLogFile, "omega/kc = %15.10g > 0.2\n", omega / (kw * speed_of_light_normalized));
		if (rank == 0) fclose(errorLogFile);
		//exit(0);
	}
	if (rank == 0) printf("omega = %g\n", omega / plasma_period);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "omega = %g\n", omega / plasma_period);

	if (rank == 0) printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));

	if (fabs(omega) > omegaGyroProton / 2) {
		if (rank == 0) printf("omega > omegaGyroProton/2\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega > omegaGyroProton/2\n");
	}
	if (rank == 0) printf("omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
	if (rank == 0) fclose(informationFile);

	checkFrequency(omega);

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");
	//checkCollisionTime(omega);
	//checkMagneticReynolds(alfvenV);
	//checkDissipation(kw, alfvenV);

	double epsilonAmplitude = amplitudeRelation;

	double alfvenVReal = omega / kw;

	if (rank == 0) fprintf(informationFile, "alfven V real = %15.10g\n", alfvenVReal * scaleFactor / plasma_period);
	if (rank == 0)
		fprintf(informationFile, "alfven V real x = %15.10g\n", alfvenVReal * (kx / kw) * scaleFactor / plasma_period);

	//double
	Bzamplitude = B0.norm() * epsilonAmplitude;

	double Omegae = omegaGyroElectron;
	double Omegae2 = Omegae * Omegae;
	double Omegap = omegaGyroProton;
	double Omegap2 = Omegap * Omegap;
	double omegae = omegaPlasmaElectron;
	double omegae2 = omegae * omegae;
	double omegap = omegaPlasmaProton;
	double omegap2 = omegap * omegap;

	double kc = kw * speed_of_light_normalized;
	double kc2 = kc * kc;

	double denominator = omega * omega - Omegae2 - (omegae2 * omega * omega / (omega * omega - kc2));

	//double
	VzamplitudeProton = -((1.0 / (4 * pi * concentration * electron_charge_normalized)) * (kc + ((omegae2 + omegap2 - omega * omega) / kc) + (omegae2 * Omegae2 / (kc * denominator))) / ((Omegae * omegae2 * omega / ((kc2 - omega * omega) * denominator)) + (Omegap / omega))) * Bzamplitude;
	//double
	VzamplitudeElectron = (((electron_charge_normalized * omega * Omegae) / (massElectron * kc)) * Bzamplitude + (omegae2 * omega * omega / (kc2 - omega * omega)) * VzamplitudeProton) / denominator;

	//double
	Byamplitude = (4 * pi * concentration * electron_charge_normalized / ((omega * omega / kc) - kc)) * (VzamplitudeElectron - VzamplitudeProton);

	//double
	VyamplitudeProton = -(Omegap / omega) * VzamplitudeProton - (electron_charge_normalized / (massProton * kc)) * Bzamplitude;
	//double
	VyamplitudeElectron = (Omegae / omega) * VzamplitudeElectron + (electron_charge_normalized / (massElectron * kc)) * Bzamplitude;

	//double
	Eyamplitude = (omega / kc) * Bzamplitude;
	//double
	Ezamplitude = -(omega / kc) * Byamplitude;

	double xshift = 0.0;

	//Eyamplitude = 0.0;
	//VzamplitudeElectron = 0.0;
	//VzamplitudeProton = 0.0;
	//Byamplitude = 0.0;


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
	//Matrix3d inverse = *(rotationMatrix.Inverse());

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k].x = 0;
				Efield[i][j][k].y = Eyamplitude * cos(kx * xgrid[i] + ky * ygrid[j] + kz * zgrid[k]);
				Efield[i][j][k].z = Ezamplitude * sin(kx * xgrid[i] + ky * ygrid[j] + kz * zgrid[k]);
				Efield[i][j][k] = rotationMatrix * Efield[i][j][k];
				//Efield[i][j][k] = inverse * Efield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}


	if (nprocs == 1) {
		for (int k = 0; k < znumber; ++k) {
			for (int j = 0; j < ynumber; ++j) {
				Efield[xnumber][j][k] = Efield[1][j][k];
				tempEfield[xnumber][j][k] = Efield[1][j][k];
				newEfield[xnumber][j][k] = Efield[1][j][k];
				explicitEfield[xnumber][j][k] = explicitEfield[1][j][k];
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			Efield[i][j][znumber] = Efield[i][j][0];
			tempEfield[i][j][znumber] = Efield[i][j][0];
			newEfield[i][j][znumber] = Efield[i][j][0];
			explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
		}
	}

	for (int k = 0; k < znumber + 1; ++k) {
		for (int i = 0; i < xnumber + 1; ++i) {
			Efield[i][ynumber][k] = Efield[i][0][k];
			tempEfield[i][ynumber][k] = Efield[i][0][k];
			newEfield[i][ynumber][k] = Efield[i][0][k];
			explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Bfield[i][j][k].x = B0.norm();
				Bfield[i][j][k].y = Byamplitude * sin(kx * middleXgrid[i] + ky * middleYgrid[j] + kz * middleZgrid[k]);
				Bfield[i][j][k].z = Bzamplitude * cos(kx * middleXgrid[i] + ky * middleYgrid[j] + kz * middleZgrid[k]);
				Bfield[i][j][k] = rotationMatrix * Bfield[i][j][k];
				//Bfield[i][j][k] = inverse * Bfield[i][j][k];
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}

	if (fabs(VzamplitudeProton) > speed_of_light_normalized) {
		printf("VzamplitudeProton > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VzamplitudeProton > speed_of_light_normalized\n");
		printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VzamplitudeProton/c = %15.10g > 1\n", VzamplitudeProton / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (rank == 0) printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);

	if (fabs(VzamplitudeElectron) > speed_of_light_normalized) {
		printf("VzamplitudeElectron > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VzamplitudeElectron > speed_of_light_normalized\n");
		printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VzamplitudeElectron/c = %15.10g > 1\n", VzamplitudeElectron / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (rank == 0) printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);

	if (fabs(VyamplitudeProton) > speed_of_light_normalized) {
		printf("VyamplitudeProton > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VyamplitudeProton > speed_of_light_normalized\n");
		printf("VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VyamplitudeProton/c = %15.10g > 1\n", VyamplitudeProton / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	if (fabs(VyamplitudeElectron) > speed_of_light_normalized) {
		printf("VyamplitudeElectron > speed_of_light_normalized\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "VyamplitudeElectron > speed_of_light_normalized\n");
		printf("VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
		fflush(stdout);
		if (rank == 0)
			fprintf(informationFile, "VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
		if (rank == 0) fclose(informationFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "VyamplitudeElectron/c = %15.10g > 1\n", VyamplitudeElectron / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	//k > 0, w > 0, kx-wt, Bz > 0, Vz < 0, Vyp > 0, Vye < 0

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		Vector3d velocity = particle->getVelocity(speed_of_light_normalized);
		int xn = (particle->coordinates.x - xgrid[0]) / deltaX;
		int yn = (particle->coordinates.y - ygrid[0]) / deltaY;
		double rightWeightX = (particle->coordinates.x - xgrid[xn]) / deltaX;
		double leftWeightX = (xgrid[xn + 1] - particle->coordinates.x) / deltaX;
		double rightWeightY = (particle->coordinates.y - ygrid[yn]) / deltaY;
		double leftWeightY = (ygrid[yn + 1] - particle->coordinates.y) / deltaY;
		if (particle->type == PROTON) {
			velocity = (Vector3d(0, 0, 1) * (VzamplitudeProton) * cos(
				kx * particle->coordinates.x + ky * particle->coordinates.y + kz * particle->coordinates.z) + Vector3d(
				0, 1, 0) * VyamplitudeProton * sin(
				kx * particle->coordinates.x + ky * particle->coordinates.y + kz * particle->coordinates.z));
			/*velocity = Vector3d(0, 0, 1) * (VzamplitudeProton) * (leftWeightX * (leftWeightY * cos(
				kx * (xgrid[xn] - xshift) + ky * ygrid[yn]) + rightWeightY * cos(
				kx * (xgrid[xn] - xshift) + ky * ygrid[yn + 1])) + rightWeightX * (leftWeightY * cos(
				kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn]) + rightWeightY * cos(
				kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn + 1])))
					 + Vector3d(0, 1, 0) * VyamplitudeProton * (leftWeightX * (leftWeightY * sin(
				kx * (xgrid[xn] - xshift) + ky * ygrid[yn]) + rightWeightY * sin(
				kx * (xgrid[xn] - xshift) + ky * ygrid[yn + 1])) + rightWeightX * (leftWeightY * sin(
				kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn]) + rightWeightY * sin(
				kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn + 1])));*/
		}
		if (particle->type == ELECTRON) {
			velocity = (Vector3d(0, 0, 1) * (VzamplitudeElectron) * cos(
				kx * particle->coordinates.x + ky * particle->coordinates.y + kz * particle->coordinates.z) + Vector3d(
				0, 1, 0) * VyamplitudeElectron * sin(
				kx * particle->coordinates.x + ky * particle->coordinates.y + kz * particle->coordinates.z));
			/*velocity = Vector3d(0, 0, 1) * (VzamplitudeElectron) * (leftWeightX * (leftWeightY * cos(
				kx * (xgrid[xn] - xshift) + ky * ygrid[yn]) + rightWeightY * cos(
				kx * (xgrid[xn] - xshift) + ky * ygrid[yn + 1])) + rightWeightX * (leftWeightY * cos(
				kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn]) + rightWeightY * cos(
				kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn + 1])))
					 + Vector3d(0, 1, 0) * VyamplitudeElectron * (leftWeightX * (leftWeightY * sin(
				kx * (xgrid[xn] - xshift) + ky * ygrid[yn]) + rightWeightY * sin(
				kx * (xgrid[xn] - xshift) + ky * ygrid[yn + 1])) + rightWeightX * (leftWeightY * sin(
				kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn]) + rightWeightY * sin(
				kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn + 1])));*/
		}
		velocity = rotationMatrix * velocity;
		//velocity = inverse * velocity;
		double beta = velocity.norm() / speed_of_light_normalized;
		particle->addVelocity(velocity, speed_of_light_normalized);
	}

	updateDeltaT();

	if (rank == 0) printf("dt/Talfven = %g\n", deltaT * omega / (2 * pi));
	if (rank == 0) printf("dt = %g\n", deltaT * plasma_period);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "dt/Talfven = %g\n", deltaT * omega / (2 * pi));
	if (rank == 0) fprintf(informationFile, "dt = %g\n", deltaT * plasma_period);

	double Vthermal = sqrt(2 * kBoltzman_normalized * temperature / massElectron);
	double thermalFlux = Vthermal * concentration * electron_charge_normalized / sqrt(1.0 * types[0].particlesPerBin);
	double alfvenFlux = (VyamplitudeProton - VyamplitudeElectron) * concentration * electron_charge_normalized;
	if (thermalFlux > alfvenFlux / 2) {
		if (rank == 0) printf("thermalFlux > alfvenFlux/2\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "thermalFlux > alfvenFlux/2\n");
	}
	if (rank == 0) printf("alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
	double minDeltaT = deltaX / Vthermal;
	if (minDeltaT > deltaT) {
		if (rank == 0) printf("deltaT < dx/Vthermal\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "deltaT < dx/Vthermal\n");

		//printf("deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);
		//fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);

		//fclose(informationFile);
		//exit(0);
	}
	if (rank == 0) {
		printf("deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
		fflush(stdout);
		fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
		fprintf(informationFile, "\n");

		fprintf(informationFile, "Bz amplitude = %g\n", Bzamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "By amplitude = %g\n", Byamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "Vz amplitude p = %g\n", VzamplitudeProton * scaleFactor / plasma_period);
		fprintf(informationFile, "Vz amplitude e = %g\n", VzamplitudeElectron * scaleFactor / plasma_period);
		fprintf(informationFile, "Vy amplitude p = %g\n", VyamplitudeProton * scaleFactor / plasma_period);
		fprintf(informationFile, "Vy amplitude e = %g\n", VyamplitudeElectron * scaleFactor / plasma_period);
		fprintf(informationFile, "Ey amplitude = %g\n", Eyamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "Ez amplitude = %g\n", Ezamplitude / (plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "By/Ez = %g\n", Byamplitude / Ezamplitude);
		fprintf(informationFile, "Bz/Ey = %g\n", Bzamplitude / Eyamplitude);
		fprintf(informationFile, "4*pi*Jy amplitude = %g\n",
		        4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "c*rotBy amplitude = %g\n",
		        speed_of_light_normalized * kw * Bzamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "4*pi*Jz amplitude = %g\n",
		        4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "c*rotBz amplitude = %g\n",
		        speed_of_light_normalized * kw * Byamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative By amplitude = %g\n",
		        -omega * Byamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "-c*rotEy = %g\n",
		        speed_of_light_normalized * kw * Ezamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Bz amplitude = %g\n",
		        omega * Bzamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "-c*rotEz = %g\n",
		        speed_of_light_normalized * kw * Eyamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Ey amplitude = %g\n",
		        omega * Eyamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "c*rotBy - 4*pi*Jy = %g\n",
		        (speed_of_light_normalized * kw * Bzamplitude - 4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Ez amplitude = %g\n",
		        -omega * Ezamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "c*rotBz - 4*pi*Jz = %g\n",
		        (speed_of_light_normalized * kw * Byamplitude - 4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativJy = -electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jy amplitude = %g\n",
		        derivativJy / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesY = electron_charge_normalized * ((Eyamplitude * ((1.0 / massProton) + (1.0 / massElectron))) + B0.norm() * ((VzamplitudeProton / massProton) + (VzamplitudeElectron / massElectron)) / speed_of_light_normalized);
		fprintf(informationFile, "dJy/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesY / (plasma_period * plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		double derivativJz = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jz amplitude = %g\n",
		        derivativJz / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesZ = electron_charge_normalized * ((Ezamplitude * ((1.0 / massProton) + (1.0 / massElectron))) - B0.norm() * ((VyamplitudeProton / massProton) + (VyamplitudeElectron / massElectron)) / speed_of_light_normalized);
		fprintf(informationFile, "dJz/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesZ / (plasma_period * plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativVyp = -omega * VyamplitudeProton;
		fprintf(informationFile, "-w*Vyp amplitude = %g\n", derivativVyp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonY = electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeProton / speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVyp/dt amplitude = %g\n",
		        derivativeVelocityProtonY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVzp = omega * VzamplitudeProton;
		fprintf(informationFile, "w*Vzp amplitude = %g\n", derivativVzp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonZ = electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeProton / speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVzp/dt amplitude = %g\n",
		        derivativeVelocityProtonZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVye = -omega * VyamplitudeElectron;
		fprintf(informationFile, "-w*Vye amplitude = %g\n", derivativVye * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronY = -electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeElectron / speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVye/dt amplitude = %g\n",
		        derivativeVelocityElectronY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVze = omega * VzamplitudeElectron;
		fprintf(informationFile, "w*Vze amplitude = %g\n", derivativVze * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronZ = -electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeElectron / speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVze/dt amplitude = %g\n",
		        derivativeVelocityElectronZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");
	}

	if (rank == 0) fclose(informationFile);
}

void Simulation::initializeLangmuirWave() {
	boundaryConditionType = PERIODIC;
	double concentration = density / (massProton + massElectron);
	types[1].particesDeltaX = types[0].particesDeltaX;
	types[1].particlesPerBin = types[0].particlesPerBin;
	types[0].concentration = concentration;
	types[1].concentration = concentration;
	for (int i = 2; i < typesNumber; ++i) {
		types[i].particlesPerBin = 0;
		types[i].concentration = 0;
		types[i].particesDeltaX = xsize;
	}
	double epsilon = 0.1;
	double kw = 2 * 2 * pi / xsize;
	double omega = 2 * pi;
	double langmuirV = omega / kw;

	checkDebyeParameter();
	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");
	if (rank == 0) fprintf(informationFile, "lengmuir V = %lf\n", langmuirV * scaleFactor / plasma_period);
	if (rank == 0) fclose(informationFile);
	if (langmuirV > speed_of_light_normalized) {
		printf("langmuirV > c\n");
		fflush(stdout);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "langmuireV/c = %15.10g > 1\n", langmuirV / speed_of_light_normalized);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	if (rank == 0) printf("creating particles\n");
	fflush(stdout);
	int nproton = 0;
	int nelectron = 0;
	double weight = (concentration / types[0].particlesPerBin) * volumeB(0, 0, 0);
	/*for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
		    for (int k = 0; k < znumber; ++k) {
		        double x;
		        for (int l = 0; l < particlesPerBin; ++l) {
		            ParticleTypes type;
		            type = PROTON;
		            Particle* particle = createParticle(particlesNumber, i, j, k, weight, type, temperature);
		            nproton++;
		            particles.push_back(particle);
		            particlesNumber++;
		            if (particlesNumber % 1000 == 0) {
		                printf("create particle number %d\n", particlesNumber);
		            }
		        }
		        for (int l = 0; l < particlesPerBin * (1 + epsilon * cos(kw * middleXgrid[i])); ++l) {
		            ParticleTypes type;
		            type = ELECTRON;
		            Particle* particle = createParticle(particlesNumber, i, j, k, weight, type, temperature);
		            nelectron++;
		            particles.push_back(particle);
		            particlesNumber++;
		            if (particlesNumber % 1000 == 0) {
		                printf("create particle number %d\n", particlesNumber);
		            }
		        }
		    }
		}
	}
	if (nproton != nelectron) {
		printf("nproton != nelectron\n");
		int n;
		ParticleTypes type;
		if (nproton > nelectron) {
		    n = nproton - nelectron;
		    type = ELECTRON;
		}
		else {
		    n = nelectron - nproton;
		    type = PROTON;
		}
		int i = 0;
		while (n > 0) {
		    Particle* particle = createParticle(particlesNumber, i, 0, 0, weight, type, temperature);
		    particles.push_back(particle);
		    particlesNumber++;
		    ++i;
		    n--;
		    if (i >= xnumber) {
		        i = 0;
		    }
		    if (type == PROTON) {
		        nproton++;
		    }
		    else {
		        nelectron++;
		    }
		}
	}
	if (nproton != nelectron) {
		printf("nproton != nelectron\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "nproton = %d nelectron = %d\n", nproton, nelectron);
		fclose(errorLogFile);
		exit(0);
	}*/
	createParticles();

	double chargeDensityAmplitude = epsilon * concentration * electron_charge_normalized;
	double Eamplitude = -4 * pi * chargeDensityAmplitude / (kw);

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k].x = Eamplitude * sin(kw * xgrid[i]);
				Efield[i][j][k].y = 0;
				Efield[i][j][k].z = 0;

				tempEfield[i] = Efield[i];
				explicitEfield[i] = Efield[i];
			}
		}
	}

	for (int j = 0; j < ynumber + 1; ++j) {
		for (int k = 0; k < znumber + 1; ++k) {
			Efield[xnumber][j][k] = Efield[0][j][k];
			tempEfield[xnumber][j][k] = tempEfield[0][j][k];
			explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
		}
	}

	double Vamplitude = electron_charge_normalized * Eamplitude / (massElectron * omega);

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		if (particle->type == ELECTRON) {
			Vector3d velocity = Vector3d(1, 0, 0) * Vamplitude * cos(kw * particle->coordinates.x);
			particle->addVelocity(velocity, speed_of_light_normalized);
		}
	}

	if (rank == 0) {
		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", 0.0, 0.0, 1.0, 1.0);
		fclose(incrementFile);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void Simulation::initializeFluxFromRight() {
	boundaryConditionType = SUPER_CONDUCTOR_LEFT;
	//boundaryConditionType = PERIODIC;
	createParticles();
	E0 = E0 - V0.vectorMult(B0) / (speed_of_light_normalized);
	//initializeAlfvenWaveY(10, 1.0E-4);
	for (int i = 0; i < xnumber + 2; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k] = E0;
				if (rank == 0) {
					if (i == 0 || i == 1) {
						Efield[i][j][k].y = 0;
						Efield[i][j][k].z = 0;
					}
				}
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}

	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				if (rank == 0) {
					additionalEfieldLeft[i][j][k] = Vector3d(0, 0, 0);
				}
				else {
					additionalEfieldLeft[i][j][k] = E0;
				}
				additionalEfieldRight[i][j][k] = E0;
			}
		}
	}

	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				additionalBfieldLeft[i][j][k] = B0;
				additionalBfieldRight[i][j][k] = B0;
			}
		}
	}

	//fieldsLorentzTransitionX(V0.x);

	/*for (int i = 0; i < particles.size(); ++i) {
		Particle* particle = particles[i];
		particle->addVelocity(V0, speed_of_light_normalized);
	}*/

	double magneticEnergy = B0.scalarMult(B0) / (8 * pi);
	double kineticEnergy = density * V0.scalarMult(V0) / 2;

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");
	if (rank == 0) fprintf(informationFile, "magneticEnergy/kineticEnergy = %15.10g\n", magneticEnergy / kineticEnergy);
	if (rank == 0) printf("magneticEnergy/kinetikEnergy = %15.10g\n", magneticEnergy / kineticEnergy);
	fflush(stdout);
	if (rank == 0) fclose(informationFile);

	checkDebyeParameter();

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", 0.0, 0.0, 1.0, 1.0);
		fclose(incrementFile);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void Simulation::fieldsLorentzTransitionX(const double& v) {
	double gamma = 1.0 / sqrt(1 - v * v / speed_of_light_normalized_sqr);
	for (int i = 0; i < xnumber; ++i) {
		int prevI = i - 1;
		if (prevI < 0) {
			if (boundaryConditionType == PERIODIC) {
				prevI = xnumber - 1;
			}
			else {
				prevI = 0;
			}
		}
		for (int j = 0; j < ynumber; ++j) {
			int prevJ = j - 1;
			if (prevJ < 0) {
				prevJ = ynumber - 1;
			}
			for (int k = 0; k < znumber; ++k) {
				int prevK = k - 1;
				if (prevK < 0) {
					prevK = znumber - 1;
				}
				Vector3d middleB = (Bfield[prevI][prevJ][prevK] + Bfield[prevI][j][prevK] + Bfield[prevI][prevJ][k] + Bfield[prevI][j][k]
					+ Bfield[i][prevJ][prevK] + Bfield[i][j][prevK] + Bfield[i][prevJ][k] + Bfield[i][j][k]) * 0.125;
				newEfield[i][j][k].y = gamma * (Efield[i][j][k].y - v * middleB.z / speed_of_light_normalized);
				newEfield[i][j][k].z = gamma * (Efield[i][j][k].z + v * middleB.y / speed_of_light_normalized);
			}
		}
	}
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Vector3d middleE = (Efield[i][j][k] + Efield[i][j + 1][k] + Efield[i][j][k + 1] + Efield[i][j + 1][k + 1]
					+ Efield[i + 1][j][k] + Efield[i + 1][j + 1][k] + Efield[i + 1][j][k + 1] + Efield[i + 1][j + 1][k + 1]) * 0.125;
				newBfield[i][j][k].y = gamma * (Bfield[i][j][k].y + v * middleE.z / speed_of_light_normalized);
				newBfield[i][j][k].z = gamma * (Bfield[i][j][k].z - v * middleE.y / speed_of_light_normalized);
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int k = 0; k < znumber; ++k) {
			newEfield[i][ynumber][k] = newEfield[i][0][k];
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			newEfield[i][j][znumber] = newEfield[i][j][0];
		}
	}

	if (boundaryConditionType == PERIODIC) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				newEfield[xnumber][j][k] = newEfield[0][j][k];
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k] = newEfield[i][j][k];
				tempEfield[i][j][k] = newEfield[i][j][k];
				explicitEfield[i][j][k] = newEfield[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Bfield[i][j][k] = newBfield[i][j][k];
			}
		}
	}
}

void Simulation::initializeShockWave() {
	boundaryConditionType = FREE_BOTH;

	E0 = Vector3d(0, 0, 0);
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k] = Vector3d(0, 0, 0);
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Bfield[i][j][k] = Vector3d(B0.x, 0, 0);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}

	if (rank == 0) printf("creating particles\n");
	fflush(stdout);
	double concentration = density / (massProton + massElectron);
	double downstreamTemperature = 1000 * temperature;
	double upstreamTemperature = temperature;
	Vector3d upstreamVelocity = V0;
	Vector3d downstreamVelocity = Vector3d(V0.x / 4, 0, 0);
	double alfvenV = B0.norm() / sqrt(4 * pi * density);
	double soundVelectron = sqrt(5 * kBoltzman_normalized * downstreamTemperature / (3 * massElectron));

	if (alfvenV > V0.norm()) {
		if (rank == 0) printf("alfvenV > V0\n");
		fflush(stdout);
	}

	if (rank == 0) printf("alfvenV/V0 = %15.10g\n", alfvenV / V0.norm());
	fflush(stdout);

	if (soundVelectron > V0.norm()) {
		if (rank == 0) printf("soundV > V0\n");
		fflush(stdout);
	}
	if (rank == 0) printf("soundV/V0 = %15.10g\n", soundVelectron / V0.norm());
	fflush(stdout);
	//Vector3d downstreamVelocity = Vector3d(0, 0, 0);
	shockWavePoint = xnumber / 2;
	int n = 0;
	//for (int i = 0; i < xnumber; ++i) {
	for (int i = 1; i < xnumber; ++i) {
		for (int typeCounter = 0; typeCounter < typesNumber; ++typeCounter) {
			double x = xgrid[i] + 0.0001 * deltaX;
			int localParticlesPerBin = types[typeCounter].particlesPerBin;
			double localTemperature = upstreamTemperature;
			if (i < shockWavePoint) {
				localParticlesPerBin = localParticlesPerBin * 4;
				localTemperature = upstreamTemperature;
			}
			double deltaXParticles = deltaX / localParticlesPerBin;
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					double weight = (types[typeCounter].concentration / types[typeCounter].particlesPerBin) * volumeB(i,
					                                                                                                  j,
					                                                                                                  k);
					for (int l = 0; l < localParticlesPerBin; ++l) {
						ParticleTypes type = types[typeCounter].type;
						Particle* particle = createParticle(n, i, j, k, weight, type, types[typeCounter],
						                                    localTemperature, localTemperature, localTemperature);
						//particle->x = middleXgrid[i];
						n++;
						/*if (l % 2 == 0) {
							x = particle->x;
						} else {
							particle->x= x;
						}*/
						if (i >= shockWavePoint) {
							particle->addVelocity(upstreamVelocity, speed_of_light_normalized);
						}
						else {
							particle->addVelocity(downstreamVelocity, speed_of_light_normalized);
						}
						particle->coordinates.x = x + deltaXParticles * l;
						Vector3d momentum = particle->getMomentum();
						particle->initialMomentum = momentum;
						//particle->prevMomentum = momentum;
						particles.push_back(particle);
						particlesNumber++;
						if (particlesNumber % 1000 == 0) {
							if (rank == 0) printf("create particle number %d\n", particlesNumber);
							fflush(stdout);
						}
					}
				}
			}
		}
	}
	synchronizeParticleNumber();

	/*for(int i = 0; i < shockWavePoint; ++i){
		//Bfield[i].y = B0.x*sin(2*20*pi*middleXgrid[i]/xsize);
		double amplitude = 0.1*B0.x;
		Bfield[i].y = amplitude*(uniformDistribution() - 0.5);
		Bfield[i].z = amplitude*(uniformDistribution() - 0.5);
		newBfield[i] = Bfield[i];
	}*/

	initializeKolmogorovSpectrum(1, 100, 0.1);

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				double v = upstreamVelocity.x;
				if (i < shockWavePoint) {
					//v = V0.x/4;
					v = downstreamVelocity.x;
				}
				double gamma = 1.0 / sqrt(1 - v * v / speed_of_light_normalized_sqr);
				Vector3d middleE = (Efield[i][j][k] + Efield[i + 1][j][k]) * 0.5;
				newBfield[i][j][k].y = gamma * (Bfield[i][j][k].y + v * middleE.z / speed_of_light_normalized);
				newBfield[i][j][k].z = gamma * (Bfield[i][j][k].z - v * middleE.y / speed_of_light_normalized);
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				double v = upstreamVelocity.x;
				if (i < shockWavePoint) {
					//v = V0.x/4;
					v = downstreamVelocity.x;
				}
				Vector3d middleE = (Efield[i][j][k] + Efield[i + 1][j][k]) * 0.5;
				double gamma = 1.0 / sqrt(1 - v * v / speed_of_light_normalized_sqr);
				newBfield[i][j][k].y = gamma * (Bfield[i][j][k].y + v * middleE.z / speed_of_light_normalized);
				newBfield[i][j][k].z = gamma * (Bfield[i][j][k].z - v * middleE.y / speed_of_light_normalized);
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			newEfield[i][j][znumber] = newEfield[i][j][0];
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int k = 0; k < znumber + 1; ++k) {
			newEfield[i][ynumber][k] = newEfield[i][0][k];
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k] = newEfield[i][j][k];
				tempEfield[i][j][k] = newEfield[i][j][k];
				explicitEfield[i][j][k] = newEfield[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Bfield[i] = newBfield[i];
			}
		}
	}
}

void Simulation::initializeKolmogorovSpectrum(int first, int last, double turbulenceFraction) {
	//use if defined shockWavePoint
	double length = xsizeGeneral;

	if (last > xnumberGeneral - 1) {
		last = xnumberGeneral - 1;
	}

	double minK = 2 * pi * first / length;
	double maxK = 2 * pi * last / length;
	double k0 = 2 * pi / length;

	double energy = turbulenceFraction * B0.scalarMult(B0);
	double amplitude = sqrt(1.5 * energy * k0 / (power(minK, -2.0 / 3.0) - power(maxK, -2.0 / 3.0)));

	double* phases = new double[2 * (last - first + 1)];

	if (rank == 0) {
		for (int i = 0; i < 2 * (last - first + 1); ++i) {
			phases[i] = 2 * pi * uniformDistribution();
		}
	}


	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(phases, 2 * (last - first + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	for (int harmCounter = first; harmCounter <= last; ++harmCounter) {
		double kw = 2 * pi * harmCounter / length;
		double Bamplitude = amplitude * power(kw, -5.0 / 6.0);
		///double phiY = 2 * pi * uniformDistribution();
		//double phiZ = 2 * pi * uniformDistribution();

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					Bfield[i][j][k].y += Bamplitude * sin(kw * middleXgrid[i] + phases[2 * harmCounter]);
					Bfield[i][j][k].z += Bamplitude * cos(kw * middleXgrid[i] + phases[2 * harmCounter + 1]);
					newBfield[i][j][k] = Bfield[i][j][k];
				}
			}
		}
	}

	delete[] phases;
}

void Simulation::initializeTwoStream() {
	boundaryConditionType = PERIODIC;
	createParticles();
	double u = 0.99 * speed_of_light_normalized;
	double gamma = 1 / sqrt(1 - u * u / speed_of_light_normalized_sqr);
	Vector3d electronsVelocityPlus = Vector3d(0, u, 0);
	Vector3d electronsVelocityMinus = Vector3d(0, -u, 0);
	//B0 = Vector3d(0, 0, 0);

	double Bamplitude = 1E-12 * (plasma_period * sqrt(scaleFactor));
	//Bamplitude = 0;
	double Eamplitude = 0;

	checkDebyeParameter();

	omegaGyroProton = B0.norm() * electron_charge_normalized / (massProton * speed_of_light_normalized);
	omegaGyroElectron = B0.norm() * electron_charge_normalized / (massElectron * speed_of_light_normalized);

	double kw = 2 * pi / xsize;

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");

	if (xsize * omegaPlasmaElectron / speed_of_light_normalized < 5) {
		if (rank == 0) printf("xsize*omegaPlasmaElectron/speed_of_light_normalized < 5\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "xsize*omegaPlasmaElectron/speed_of_light_normalized < 5\n");
	}
	if (rank == 0)
		printf("xsize*omegaPlasmaElectron/speed_of_light_normalized = %g\n",
		       xsize * omegaPlasmaElectron / speed_of_light_normalized);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "xsize*omegaPlasmaElectron/speed_of_light_normalized = %g\n",
		        xsize * omegaPlasmaElectron / speed_of_light_normalized);

	if (deltaX * omegaPlasmaElectron / speed_of_light_normalized > 0.2) {
		if (rank == 0) printf("deltaX*omegaPlasmaElectron/speed_of_light_normalized > 0.2\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "deltaX*omegaPlasmaElectron/speed_of_light_normalized > 0.2\n");
	}
	if (rank == 0)
		printf("deltaX*omegaPlasmaElectron/speed_of_light_normalized = %g\n",
		       xsize * omegaPlasmaElectron / speed_of_light_normalized);
	fflush(stdout);
	if (rank == 0)
		fprintf(informationFile, "deltaX*omegaPlasmaElectron/speed_of_light_normalized = %g\n",
		        xsize * omegaPlasmaElectron / speed_of_light_normalized);

	if (kw > omegaPlasmaElectron / u) {
		if (rank == 0) printf("k > omegaPlasmaElectron/u\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "k > omegaPlasmaElectron/u\n");
	}

	if (rank == 0) printf("k u/omegaPlasmaElectron = %g\n", kw * u / omegaPlasmaElectron);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "k u/omegaPlasmaElectron = %g\n", kw * u / omegaPlasmaElectron);
	if (rank == 0) fclose(informationFile);

	double omegaGyroHelium = B0.norm() * electron_charge_normalized / (massHelium3 * speed_of_light_normalized);

	if (rank == 0) {
		double increment = (u * omegaPlasmaElectron / (speed_of_light_normalized)) / sqrt(gamma);

		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", increment, increment / plasma_period, 1.0, 1.0);
		fclose(incrementFile);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	/*for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
		    for (int k = 0; k < znumber; ++k) {
		        Bfield[i][j][k] += Vector3d(0, 1, 0) * Bamplitude * cos(kw * middleXgrid[i]);
		        newBfield[i][j][k] = Bfield[i][j][k];
		    }
		}
	}*/

	/*for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
		    for (int k = 0; k < znumber + 1; ++k) {
		        Efield[i][j][k] = Vector3d(0, 0, 1) * Eamplitude * cos(kw * xgrid[i]);
		        tempEfield[i][j][k] = Efield[i][j][k];
		        newEfield[i][j][k] = Efield[i][j][k];
		        explicitEfield[i][j][k] = Efield[i][j][k];
		    }
		}
	}*/
	int electronCount = 0;
	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		if (particle->type == ELECTRON) {
			if (electronCount % 2 == 0) {
				particle->addVelocity(electronsVelocityPlus, speed_of_light_normalized);
			}
			else {
				particle->addVelocity(electronsVelocityMinus, speed_of_light_normalized);
			}
			electronCount++;
		}
	}
}

void Simulation::initializeExternalFluxInstability() {
	boundaryConditionType = PERIODIC;
	createParticles();
	double alfvenV = B0.norm() / sqrt(4 * pi * density);
	double concentration = density / (massProton + massElectron);
	double phaseV = 2 * alfvenV;
	double kw = 2 * pi / xsize;
	double omega = kw * phaseV;

	extJ = 0.001 * electron_charge_normalized * alfvenV * concentration;
	checkDebyeParameter();

	double Byamplitude = 4 * pi * sqr(alfvenV) * extJ / (speed_of_light_normalized * kw * (sqr(phaseV) - sqr(
		alfvenV))) / (plasma_period * sqrt(scaleFactor));
	double Uyamplitude = B0.norm() * phaseV * extJ / (kw * density * speed_of_light_normalized * (sqr(
		phaseV) - sqr(alfvenV))) * (scaleFactor / plasma_period);
	double cyclothronOmegaElectron = electron_charge_normalized * B0.norm() / (massElectron * speed_of_light_normalized);
	double cyclothronOmegaProton = electron_charge_normalized * B0.norm() / (massProton * speed_of_light_normalized);


	checkGyroRadius();
	if (rank == 0) {
		informationFile = fopen((outputDir + "information.dat").c_str(), "a");
		fprintf(informationFile, "alfven V = %g\n", alfvenV * scaleFactor / plasma_period);
		fprintf(informationFile, "phase V = %g\n", phaseV * scaleFactor / plasma_period);
		fprintf(informationFile, "alfven V/c = %g\n", alfvenV / speed_of_light_normalized);
		fprintf(informationFile, "phase V/c = %g\n", phaseV / speed_of_light_normalized);
		fprintf(informationFile, "external flux = %g\n", extJ * plasma_period * plasma_period * sqrt(scaleFactor));
		fprintf(informationFile, "By max amplitude = %g\n", Byamplitude);
		fprintf(informationFile, "Uy max amplitude = %g\n", Uyamplitude);
		fprintf(informationFile, "omega/omega_plasma = %g\n", omega);
		if (omega > cyclothronOmegaProton) {
			printf("omega > cyclothron Omega Proton\n");
			fflush(stdout);
			fprintf(informationFile, "omega > cyclothron Omega Proton\n");
		}
		else if (omega > cyclothronOmegaProton / 100.0) {
			printf("omega > cyclothrone Omega Proton/100\n");
			fflush(stdout);
			fprintf(informationFile, "omega > cyclothron Omega Proton/100\n");
		}
		printf("omega/cyclothronOmega = %g\n", omega / cyclothronOmegaProton);
		fprintf(informationFile, "omega/cyclothronOmega = %g\n", omega / cyclothronOmegaProton);

		fclose(informationFile);
	}
}

void Simulation::initializeAnisotropic() {
	boundaryConditionType = PERIODIC;
	/*double Tx = temperature*10;
	double Ty = temperature*1000;
	double Tz = temperature*1000;*/

	double temperatureIon = 6.9E7;
	double temperatureElectron = 6.9E7;
	//double temperatureParallelMinority = 4.63E8;
	//double temperatureParallelMinority = 115.9E8;
	double temperatureParallelMinority = 6.9E7;
	double temperatureNormalMinority = 115.9E8;
	//double temperatureNormalMinority = 6.9E7;

	for (int i = 0; i < typesNumber; ++i) {
		types[i].temperatureX = temperatureIon;
		types[i].temperatureY = temperatureIon;
		types[i].temperatureZ = temperatureIon;
	}

	types[0].temperatureX = temperatureElectron;
	types[0].temperatureY = temperatureElectron;
	types[0].temperatureZ = temperatureElectron;

	types[5].temperatureX = temperatureParallelMinority;
	types[5].temperatureY = temperatureNormalMinority;
	types[5].temperatureZ = temperatureNormalMinority;

	checkDebyeParameter();

	createParticles();

	omegaPlasmaProton = sqrt(
		4 * pi * types[1].concentration * electron_charge_normalized * electron_charge_normalized / massProton);
	omegaPlasmaElectron = sqrt(
		4 * pi * types[0].concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
	double omegaPlasmaAlpha = sqrt(
		4 * pi * types[3].concentration * 2 * electron_charge_normalized * 2 * electron_charge_normalized / massAlpha);
	omegaGyroProton = B0.norm() * electron_charge_normalized / (massProton * speed_of_light_normalized);
	omegaGyroElectron = B0.norm() * electron_charge_normalized / (massElectron * speed_of_light_normalized);
	double omegaGyroAlpha = B0.norm() * 2 * electron_charge_normalized / (massAlpha * speed_of_light_normalized);

	double kmin = 2 * pi / xsizeGeneral;
	double kmax = pi / deltaX;

	double vthermalProton = sqrt(kBoltzman_normalized * types[1].temperatureX / massProton);

	if (rank == 0) printf("delta Omega/kmax vth = %g\n", (omegaGyroAlpha - omegaGyroProton) / (kmax * vthermalProton));
	if (rank == 0) printf("delta Omega/kmin vth = %g\n", (omegaGyroAlpha - omegaGyroProton) / (kmin * vthermalProton));

	if (rank == 0) printf("omega plasma/gyro omega protons = %g\n", omegaPlasmaProton / omegaGyroProton);
	if (rank == 0) printf("omega plasma/gyro omega alphas = %g\n", omegaPlasmaAlpha / omegaGyroAlpha);
	if (rank == 0) printf("omega plasma/gyro omega electrons = %g\n", omegaPlasmaElectron / omegaGyroElectron);
	fflush(stdout);

	initializeKolmogorovSpectrum(1, 100, 0.0000001);

	double omegaGyroHelium = B0.norm() * electron_charge_normalized / (massHelium3 * speed_of_light_normalized);

	if (rank == 0) {
		double increment = 0.01 * omegaGyroHelium;

		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		//todo length
		fprintf(incrementFile, "%g %g %g %g\n", increment, increment / plasma_period, 0.0, 0.0);
		fclose(incrementFile);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void Simulation::initializeAnisotropicSilicon() {
	boundaryConditionType = PERIODIC;
	/*double Tx = temperature*10;
	double Ty = temperature*1000;
	double Tz = temperature*1000;*/

	double temperatureIon = 3.67E7;
	double temperatureElectron = 3.69E7;
	//double temperatureParallelMinority = 4.63E8;
	//double temperatureParallelMinority = 115.9E8;
	double temperatureParallelMinority = 8.61E7;
	double temperatureNormalMinority = 1.31E8;
	//double temperatureNormalMinority = 6.9E7;

	for (int i = 0; i < typesNumber; ++i) {
		types[i].temperatureX = temperatureIon;
		types[i].temperatureY = temperatureIon;
		types[i].temperatureZ = temperatureIon;
	}

	types[0].temperatureX = temperatureElectron;
	types[0].temperatureY = temperatureElectron;
	types[0].temperatureZ = temperatureElectron;

	types[7].temperatureX = temperatureParallelMinority;
	types[7].temperatureY = temperatureNormalMinority;
	types[7].temperatureZ = temperatureNormalMinority;

	checkDebyeParameter();

	createParticles();

	omegaPlasmaProton = sqrt(
		4 * pi * types[1].concentration * electron_charge_normalized * electron_charge_normalized / massProton);
	omegaPlasmaElectron = sqrt(
		4 * pi * types[0].concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
	double omegaPlasmaOxygen = sqrt(
		4 * pi * types[6].concentration * 3 * electron_charge_normalized * 3 * electron_charge_normalized / massOxygen);
	omegaGyroProton = B0.norm() * electron_charge_normalized / (massProton * speed_of_light_normalized);
	omegaGyroElectron = B0.norm() * electron_charge_normalized / (massElectron * speed_of_light_normalized);
	double omegaGyroSilicon = B0.norm() * electron_charge_normalized / (massSilicon * speed_of_light_normalized);

	double kmin = 2 * pi / xsizeGeneral;
	double kmax = pi / deltaX;

	double vthermalProton = sqrt(kBoltzman_normalized * types[1].temperatureX / massProton);


	initializeKolmogorovSpectrum(1, 100, 0.0000001);

	if (rank == 0) {
		double increment = 0.011 * omegaGyroSilicon;
		double length = 180 * speed_of_light_normalized / omegaPlasmaOxygen;

		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		//todo length
		fprintf(incrementFile, "%g %g %g %g\n", increment, increment / plasma_period, length, length * scaleFactor);
		fclose(incrementFile);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}


void Simulation::initializeWeibel() {
	boundaryConditionType = PERIODIC;


	types[0].alphaNormal = 10;
	//types[0].alphaNormal = 1;
	//types[0].alphaParallel = 1.3*types[0].alphaNormal;
	types[0].alphaParallel = 1.1 * types[0].alphaNormal;

	types[0].temperatureX = ((types[0].mass * speed_of_light_normalized_sqr / types[0].alphaParallel) / (1 + (types[0].alphaParallel / types[0].alphaNormal - 1) * McDonaldFunction(
		types[0].alphaParallel, 1) / (types[0].alphaParallel * McDonaldFunction(types[0].alphaParallel,
		                                                                        2)))) / kBoltzman_normalized;
	types[0].temperatureY = (types[0].mass * speed_of_light_normalized_sqr / types[0].alphaNormal) / kBoltzman_normalized;
	types[0].temperatureZ = types[0].temperatureY;

	double temperatureIon = types[0].temperatureX / 1000;

	for (int i = 1; i < typesNumber; ++i) {
		types[i].temperatureX = temperatureIon;
		types[i].temperatureY = temperatureIon;
		types[i].temperatureZ = temperatureIon;
	}

	checkDebyeParameter();

	createParticles();

	omegaPlasmaProton = sqrt(
		4 * pi * types[1].concentration * electron_charge_normalized * electron_charge_normalized / massProton);
	omegaPlasmaElectron = sqrt(
		4 * pi * types[0].concentration * electron_charge_normalized * electron_charge_normalized / massElectron);

	omegaGyroProton = B0.norm() * electron_charge_normalized / (massProton * speed_of_light_normalized);
	omegaGyroElectron = B0.norm() * electron_charge_normalized / (massElectron * speed_of_light_normalized);

	if (rank == 0) printf("omega plasma/gyro omega protons = %g\n", omegaPlasmaProton / omegaGyroProton);
	if (rank == 0) printf("omega plasma/gyro omega electrons = %g\n", omegaPlasmaElectron / omegaGyroElectron);

	initializeKolmogorovSpectrum(1, 100, 0.0000001);

	if (rank == 0) printf("evaluating increment\n");
	if (rank == 0) {
		double alphaNormal = types[0].alphaNormal;
		double alphaParallel = types[0].alphaParallel;

		//todo check indexes of alpha and mcDonald ... somehow

		double k02 = (omegaPlasmaElectron * omegaPlasmaElectron / speed_of_light_normalized_sqr) * (alphaParallel / alphaNormal - 1)
			* (McDonaldFunction(alphaParallel, 1) / McDonaldFunction(alphaParallel, 2) + McDonaldFunction(
				alphaParallel, 0) / (alphaParallel * McDonaldFunction(alphaParallel, 2)))
			/ (1 + (alphaParallel / alphaNormal - 1) * McDonaldFunction(alphaParallel,
			                                                            1) / (alphaParallel * McDonaldFunction(
				alphaParallel, 2)));
		double kmaxIncrement = sqrt(k02 / 3);
		double length = 2 * pi / kmaxIncrement;
		double integral = 0;
		double maxTau = 20 / alphaParallel;
		double dtau = maxTau / 5000;
		double tempTau = 0;
		for (int i = 0; i < 5000; ++i) {
			double tau = tempTau + dtau / 2;
			double dIntegral = dtau * McDonaldFunction(alphaParallel * sqrt(1.0 + tau * tau), 2) / (1.0 + tau * tau);
			alertNaNOrInfinity(dIntegral, "temp = NaN in increment weibel\n");
			integral = integral + dIntegral;
			tempTau += dtau;
		}

		double temp = ((kmaxIncrement * kmaxIncrement / k02) * sqr(
			1 - kmaxIncrement * kmaxIncrement / k02) * (kBoltzman_normalized * types[0].temperatureX / (types[0].mass * speed_of_light_normalized)) * cube(
			types[0].temperatureY / types[0].temperatureX - 1) * cube(
			alphaNormal * McDonaldFunction(alphaParallel, 2) + McDonaldFunction(alphaParallel, 0)) / cube(
			alphaNormal * McDonaldFunction(alphaParallel, 2) + McDonaldFunction(alphaParallel, 1))) * sqr(
			McDonaldFunction(alphaParallel, 2)) / sqr(integral);
		double increment = omegaPlasmaElectron * sqrt(temp);
		alertNaNOrInfinity(increment, "increment = NaN in weibel\n");

		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", increment, increment / plasma_period, length, length / scaleFactor);
		fclose(incrementFile);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) printf("finish initialize weibel\n");
}


void Simulation::initializeRingWeibel() {
	boundaryConditionType = PERIODIC;

	/*double temperatureIon = types[1].temperatureX;

	for(int i = 0; i < typesNumber; ++i){
		types[i].temperatureX = temperatureIon;
		types[i].temperatureY = temperatureIon;
		types[i].temperatureZ = temperatureIon;
	}*/

	checkDebyeParameter();

	createParticles();

	omegaPlasmaProton = sqrt(
		4 * pi * types[1].concentration * electron_charge_normalized * electron_charge_normalized / massProton);
	omegaPlasmaElectron = sqrt(
		4 * pi * types[0].concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
	double omegaPlasmaAlpha = sqrt(
		4 * pi * types[3].concentration * 2 * electron_charge_normalized * 2 * electron_charge_normalized / massAlpha);
	omegaGyroProton = B0.norm() * electron_charge_normalized / (massProton * speed_of_light_normalized);
	omegaGyroElectron = B0.norm() * electron_charge_normalized / (massElectron * speed_of_light_normalized);

	if (rank == 0) printf("omega plasma/gyro omega protons = %g\n", omegaPlasmaProton / omegaGyroProton);
	if (rank == 0) printf("omega plasma/gyro omega electrons = %g\n", omegaPlasmaElectron / omegaGyroElectron);
	fflush(stdout);

	initializeKolmogorovSpectrum(1, 100, 0.0000001);

	double pNormal = 20 * massElectron * speed_of_light_normalized;
	double pParallel = 2 * massElectron * speed_of_light_normalized;

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		if (particle->type == ELECTRON) {
			Vector3d momentum;
			momentum.x = pParallel * (2 * uniformDistribution() - 1.0);
			double phi = 2 * pi * uniformDistribution();
			momentum.y = pNormal * cos(phi);
			momentum.z = pNormal * sin(phi);

			particle->setMomentum(momentum);
		}
	}

	double m2c2 = massElectron * massElectron * speed_of_light_normalized_sqr;

	double gamma = sqrt(1 + pNormal * pNormal / m2c2 + pParallel * pParallel / m2c2);
	double betaNormal = pNormal / (gamma * massElectron * speed_of_light_normalized);
	double betaParallel = pParallel / (gamma * massElectron * speed_of_light_normalized);

	double G = 0.5 * (log((1 + betaParallel) / (1 - betaParallel))) / betaParallel;


	if (rank == 0) {

		double increment;
		double length;


		double k02 = (omegaPlasmaElectron * omegaPlasmaElectron / (gamma * speed_of_light_normalized_sqr)) * ((betaNormal * betaNormal / (2 * betaParallel * betaParallel * (1 - betaParallel * betaParallel))) - G);
		if (k02 < 0) {
			increment = 0;
			length = 0;
		}
		else {

			double a = omegaPlasmaElectron * omegaPlasmaElectron * betaNormal * betaNormal / (2 * gamma * betaParallel * betaParallel);

			double bParallel2 = betaParallel * betaParallel;

			double kmaxIncrement2 = ((speed_of_light_normalized_sqr * k02 - a - bParallel2 * (speed_of_light_normalized_sqr * k02 + a)) +
					sqrt(a * sqr(
						bParallel2 + 1) * (bParallel2 * speed_of_light_normalized_sqr * k02 - speed_of_light_normalized_sqr * k02 + a))) /
				sqr(speed_of_light_normalized * (betaParallel + 1) * (betaParallel - 1));

			if (kmaxIncrement2 < 0) {
				increment = 0;
				length = 0;
			}
			else {

				double kmaxIncrement = sqrt(kmaxIncrement2);
				length = 2 * pi / kmaxIncrement;


				increment = (1 / sqrt(2.0)) * sqrt(
					sqrt(
						sqr(speed_of_light_normalized_sqr * kmaxIncrement2 * bParallel2 + a - speed_of_light_normalized_sqr * (k02 - kmaxIncrement2)) + 4 * speed_of_light_normalized_sqr * speed_of_light_normalized_sqr * kmaxIncrement2 * bParallel2 * (k02 - kmaxIncrement2)) -
					(speed_of_light_normalized_sqr * kmaxIncrement2 * bParallel2 + a - speed_of_light_normalized_sqr * (k02 - kmaxIncrement2))
				);
			}
		}

		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", increment, increment / plasma_period, length, length / scaleFactor);
		fclose(incrementFile);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void Simulation::initializeHomogenouseFlow() {
	createParticles();
	E0 = E0 - V0.vectorMult(B0) / (speed_of_light_normalized);
	//initializeAlfvenWaveY(10, 1.0E-4);
	for (int i = 0; i < xnumber + 2; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k] = E0;
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}

	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				additionalEfieldLeft[i][j][k] = E0;
				additionalEfieldRight[i][j][k] = E0;
			}
		}
	}

	for (int i = 0; i < additionalBinNumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				additionalBfieldLeft[i][j][k] = B0;
				additionalBfieldRight[i][j][k] = B0;
			}
		}
	}
	for (int p = 0; p < particles.size(); ++p) {
		Particle* particle = particles[p];
		Vector3d momentum = V0 * particle->mass / sqrt(1 - V0.scalarMult(V0) / speed_of_light_normalized_sqr);
		particle->setMomentum(momentum);
	}
}

void Simulation::createArrays() {
	if (rank == 0) printf("creating arrays\n");
	if (rank == 0) fflush(stdout);
	//if(rank == 0) printLog("creating arrays\n");

	if (rank == 0) printf("creating grid arrays\n");
	if (rank == 0) fflush(stdout);
	// if(rank == 0) printLog("creating grid arrays\n");
	xgrid = new double[xnumber + 2];
	ygrid = new double[ynumber + 1];
	zgrid = new double[znumber + 1];

	middleXgrid = new double[xnumber + 1];
	middleYgrid = new double[ynumber];
	middleZgrid = new double[znumber];

	if (rank == 0) printf("creating gmresOutput arrays\n");
	if (rank == 0) fflush(stdout);

	gmresOutput = new double ***[xnumber + 1];
	for (int i = 0; i < xnumber + 1; ++i) {
		gmresOutput[i] = new double **[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			gmresOutput[i][j] = new double *[znumber];
			for (int k = 0; k < znumber; ++k) {
				gmresOutput[i][j][k] = new double[maxwellEquationMatrixSize];
			}
		}
	}

	if (rank == 0) printf("creating fields arrays\n");
	if (rank == 0) fflush(stdout);
	// if(rank == 0) printLog("creating fields arrays\n");

	Efield = new Vector3d **[xnumber + 2];
	newEfield = new Vector3d **[xnumber + 2];
	tempEfield = new Vector3d **[xnumber + 2];
	smoothingEfield = new Vector3d **[xnumber + 2];
	explicitEfield = new Vector3d **[xnumber + 2];
	rotB = new Vector3d **[xnumber + 2];
	Ederivative = new Vector3d **[xnumber + 2];
	Bfield = new Vector3d **[xnumber + 1];
	newBfield = new Vector3d **[xnumber + 1];
	smoothingBfield = new Vector3d **[xnumber + 1];
	rotE = new Vector3d **[xnumber + 1];
	Bderivative = new Vector3d **[xnumber + 1];

	for (int i = 0; i < xnumber + 1; ++i) {
		Bfield[i] = new Vector3d *[ynumber];
		newBfield[i] = new Vector3d *[ynumber];
		smoothingBfield[i] = new Vector3d *[ynumber];
		rotE[i] = new Vector3d *[ynumber];
		Bderivative[i] = new Vector3d *[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			//if(rank == 0) printf("%d %d\n", i, j);
			//if(rank == 0) fflush((stdout));
			Bfield[i][j] = new Vector3d[znumber];
			newBfield[i][j] = new Vector3d[znumber];
			smoothingBfield[i][j] = new Vector3d[znumber];
			rotE[i][j] = new Vector3d[znumber];
			Bderivative[i][j] = new Vector3d[znumber];
			for (int k = 0; k < znumber; ++k) {
				Bfield[i][j][k] = Vector3d(0, 0, 0);
				newBfield[i][j][k] = Vector3d(0, 0, 0);
				smoothingBfield[i][j][k] = Vector3d(0, 0, 0);
				rotE[i][j][k] = Vector3d(0, 0, 0);
				Bderivative[i][j][k] = Vector3d(0, 0, 0);
			}
		}
	}

	for (int i = 0; i < xnumber + 2; ++i) {
		Efield[i] = new Vector3d *[ynumber + 1];
		newEfield[i] = new Vector3d *[ynumber + 1];
		tempEfield[i] = new Vector3d *[ynumber + 1];
		smoothingEfield[i] = new Vector3d *[ynumber + 1];
		explicitEfield[i] = new Vector3d *[ynumber + 1];
		rotB[i] = new Vector3d *[ynumber + 1];
		Ederivative[i] = new Vector3d *[ynumber + 1];
		for (int j = 0; j < ynumber + 1; ++j) {
			//if(rank == 0) printf("%d %d\n", i, j);
			//if(rank == 0) fflush((stdout));
			Efield[i][j] = new Vector3d[znumber + 1];
			newEfield[i][j] = new Vector3d[znumber + 1];
			tempEfield[i][j] = new Vector3d[znumber + 1];
			smoothingEfield[i][j] = new Vector3d[znumber + 1];
			explicitEfield[i][j] = new Vector3d[znumber + 1];
			rotB[i][j] = new Vector3d[znumber + 1];
			Ederivative[i][j] = new Vector3d[znumber + 1];
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k] = Vector3d(0, 0, 0);
				newEfield[i][j][k] = Vector3d(0, 0, 0);
				tempEfield[i][j][k] = Vector3d(0, 0, 0);
				smoothingEfield[i][j][k] = Vector3d(0, 0, 0);
				explicitEfield[i][j][k] = Vector3d(0, 0, 0);
				rotB[i][j][k] = Vector3d(0, 0, 0);
				Ederivative[i][j][k] = Vector3d(0, 0, 0);
			}
		}
	}

	if (rank == 0) printf("creating maxwellequation matrix arrays\n");
	fflush(stdout);
	//if(rank == 0) printLog("creating maxwell equation matrix arrays\n");

	maxwellEquationMatrix = new std::vector<MatrixElement> ***[xnumber + 1];
	maxwellEquationRightPart = new double ***[xnumber + 1];
	for (int i = 0; i < xnumber + 1; ++i) {
		maxwellEquationMatrix[i] = new std::vector<MatrixElement> **[ynumber];
		maxwellEquationRightPart[i] = new double **[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			maxwellEquationMatrix[i][j] = new std::vector<MatrixElement> *[znumber];
			maxwellEquationRightPart[i][j] = new double *[znumber];
			for (int k = 0; k < znumber; ++k) {
				//if(rank == 0) printf("%d %d %d\n", i, j, k);
				//if(rank == 0) fflush((stdout));
				maxwellEquationMatrix[i][j][k] = new std::vector<MatrixElement>[maxwellEquationMatrixSize];
				maxwellEquationRightPart[i][j][k] = new double[maxwellEquationMatrixSize];
			}
		}
	}

	if (rank == 0) printf("creating arrays for divergence\n");
	if (rank == 0) fflush(stdout);
	//if(rank == 0) printLog("creating arrays for divergence\n");

	divergenceCleanUpMatrix = new std::vector<MatrixElement> ***[xnumber + 2];
	divergenceCleanUpRightPart = new double ***[xnumber + 2];

	divergenceCleaningField = new double ***[xnumber + 2];
	divergenceCleaningPotential = new double ***[xnumber + 2];
	tempDivergenceCleaningPotential = new double ***[xnumber + 2];
	divergenceCleaningPotentialFourier = new double **[xnumber + 1];

	for (int i = 0; i < xnumber + 2; ++i) {
		divergenceCleanUpMatrix[i] = new std::vector<MatrixElement> **[ynumber];
		divergenceCleanUpRightPart[i] = new double **[ynumber];
		divergenceCleaningField[i] = new double **[ynumber];
		divergenceCleaningPotential[i] = new double **[ynumber];
		tempDivergenceCleaningPotential[i] = new double **[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			divergenceCleanUpMatrix[i][j] = new std::vector<MatrixElement> *[znumber];
			divergenceCleanUpRightPart[i][j] = new double *[znumber];
			divergenceCleaningField[i][j] = new double *[znumber];
			divergenceCleaningPotential[i][j] = new double *[znumber];
			tempDivergenceCleaningPotential[i][j] = new double *[znumber];
			for (int k = 0; k < znumber; ++k) {
				//if(rank == 0) printf("%d %d %d\n", i, j, k);
				//if(rank == 0) fflush((stdout));
				divergenceCleaningField[i][j][k] = new double[3];
				divergenceCleanUpRightPart[i][j][k] = new double[3];
				divergenceCleanUpMatrix[i][j][k] = new std::vector<MatrixElement>[3];
				divergenceCleaningPotential[i][j][k] = new double[1];
				tempDivergenceCleaningPotential[i][j][k] = new double[1];
				divergenceCleaningPotential[i][j][k][0] = 0;
				tempDivergenceCleaningPotential[i][j][k][0] = 0;
				for (int l = 0; l < 3; ++l) {
					divergenceCleaningField[i][j][k][l] = 0;
					divergenceCleanUpRightPart[i][j][k][l] = 0;
				}
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		divergenceCleaningPotentialFourier[i] = new double *[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			divergenceCleaningPotentialFourier[i][j] = new double[znumber];
			for (int k = 0; k < znumber; ++k) {
				divergenceCleaningPotentialFourier[i][j][k] = 0;
			}
		}
	}

	if (rank == 0) printf("creating arrays for particlesInBins\n");
	if (rank == 0) fflush(stdout);

	if (rank == 0) printf("creating arrays for parameters\n");
	fflush(stdout);
	//if(rank == 0) printLog("creating arrays for parameters\n");

	particleConcentrations = new double ***[typesNumber];
	particleBulkVelocities = new Vector3d ***[typesNumber];
	for (int t = 0; t < typesNumber; ++t) {
		particleConcentrations[t] = new double **[xnumber + 1];
		particleBulkVelocities[t] = new Vector3d **[xnumber + 1];
		for (int i = 0; i < xnumber + 1; ++i) {
			particleConcentrations[t][i] = new double *[ynumber];
			particleBulkVelocities[t][i] = new Vector3d *[ynumber];
			for (int j = 0; j < ynumber; ++j) {
				particleConcentrations[t][i][j] = new double[znumber];
				particleBulkVelocities[t][i][j] = new Vector3d[znumber];
				for (int k = 0; k < znumber; ++k) {
					particleConcentrations[t][i][j][k] = 0;
					particleBulkVelocities[t][i][j][k] = Vector3d(0, 0, 0);
				}
			}
		}
	}

	additionalParticleConcentrationsLeft = new double ***[typesNumber];
	additionalParticleConcentrationsRight = new double ***[typesNumber];
	additionalParticleBulkVelocitiesLeft = new Vector3d ***[typesNumber];
	additionalParticleBulkVelocitiesRight = new Vector3d ***[typesNumber];
	if (additionalBinNumber > 0) {
		for (int t = 0; t < typesNumber; ++t) {
			additionalParticleConcentrationsLeft[t] = new double **[additionalBinNumber];
			additionalParticleConcentrationsRight[t] = new double **[additionalBinNumber];
			additionalParticleBulkVelocitiesLeft[t] = new Vector3d **[additionalBinNumber];
			additionalParticleBulkVelocitiesRight[t] = new Vector3d **[additionalBinNumber];
			for (int i = 0; i < additionalBinNumber; ++i) {
				additionalParticleConcentrationsLeft[t][i] = new double *[ynumber];
				additionalParticleConcentrationsRight[t][i] = new double *[ynumber];
				additionalParticleBulkVelocitiesLeft[t][i] = new Vector3d *[ynumber];
				additionalParticleBulkVelocitiesRight[t][i] = new Vector3d *[ynumber];
				for (int j = 0; j < ynumber; ++j) {
					additionalParticleConcentrationsLeft[t][i][j] = new double[znumber];
					additionalParticleConcentrationsRight[t][i][j] = new double[znumber];
					additionalParticleBulkVelocitiesLeft[t][i][j] = new Vector3d[znumber];
					additionalParticleBulkVelocitiesRight[t][i][j] = new Vector3d[znumber];
					for (int k = 0; k < znumber; ++k) {
						additionalParticleConcentrationsLeft[t][i][j][k] = 0;
						additionalParticleConcentrationsRight[t][i][j][k] = 0;
						additionalParticleBulkVelocitiesLeft[t][i][j][k] = Vector3d(0, 0, 0);
						additionalParticleBulkVelocitiesRight[t][i][j][k] = Vector3d(0, 0, 0);
					}
				}
			}
		}
	}

	chargeDensity = new double **[xnumber + 1];
	chargeDensityMinus = new double **[xnumber + 1];
	chargeDensityHat = new double **[xnumber + 1];
	pressureTensor = new Matrix3d **[xnumber + 1];
	tempCellParameter = new double **[xnumber + 1];

	for (int i = 0; i < xnumber + 1; ++i) {
		chargeDensity[i] = new double *[ynumber];
		chargeDensityMinus[i] = new double *[ynumber];
		chargeDensityHat[i] = new double *[ynumber];
		pressureTensor[i] = new Matrix3d *[ynumber];
		tempCellParameter[i] = new double *[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			chargeDensity[i][j] = new double[znumber];
			chargeDensityMinus[i][j] = new double[znumber];
			chargeDensityHat[i][j] = new double[znumber];
			pressureTensor[i][j] = new Matrix3d[znumber];
			tempCellParameter[i][j] = new double[znumber];
			for (int k = 0; k < znumber; ++k) {
				//if(rank == 0) printf("%d %d %d\n", i, j, k);
				//if(rank == 0) fflush((stdout));
				chargeDensity[i][j][k] = 0;
				chargeDensityMinus[i][j][k] = 0;
				chargeDensityHat[i][j][k] = 0;
				pressureTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				tempCellParameter[i][j][k] = 0;
			}
		}
	}

	if (additionalBinNumber > 0) {
		additionalBfieldLeft = new Vector3d **[additionalBinNumber];
		additionalBfieldRight = new Vector3d **[additionalBinNumber];
		additionalNewBfieldLeft = new Vector3d **[additionalBinNumber];
		additionalNewBfieldRight = new Vector3d **[additionalBinNumber];
		additionalChargeDensityHatLeft = new double **[additionalBinNumber];
		additionalChargeDensityHatRight = new double **[additionalBinNumber];
		additionalChargeDensityLeft = new double **[additionalBinNumber];
		additionalChargeDensityMinusLeft = new double **[additionalBinNumber];
		additionalChargeDensityRight = new double **[additionalBinNumber];
		additionalChargeDensityMinusRight = new double **[additionalBinNumber];
		additionalChargeDensityLeft = new double **[additionalBinNumber];
		additionalChargeDensityRight = new double **[additionalBinNumber];
		additionalPressureTensorLeft = new Matrix3d **[additionalBinNumber];
		additionalPressureTensorRight = new Matrix3d **[additionalBinNumber];
		for (int i = 0; i < additionalBinNumber; ++i) {
			additionalBfieldLeft[i] = new Vector3d *[ynumber];
			additionalBfieldRight[i] = new Vector3d *[ynumber];
			additionalNewBfieldLeft[i] = new Vector3d *[ynumber];
			additionalNewBfieldRight[i] = new Vector3d *[ynumber];
			additionalChargeDensityHatLeft[i] = new double *[ynumber];
			additionalChargeDensityHatRight[i] = new double *[ynumber];
			additionalChargeDensityLeft[i] = new double *[ynumber];
			additionalChargeDensityMinusLeft[i] = new double *[ynumber];
			additionalChargeDensityRight[i] = new double *[ynumber];
			additionalChargeDensityMinusRight[i] = new double *[ynumber];
			additionalChargeDensityLeft[i] = new double *[ynumber];
			additionalChargeDensityRight[i] = new double *[ynumber];
			additionalPressureTensorLeft[i] = new Matrix3d *[ynumber];
			additionalPressureTensorRight[i] = new Matrix3d *[ynumber];
			for (int j = 0; j < ynumber; ++j) {
				additionalBfieldLeft[i][j] = new Vector3d[znumber];
				additionalBfieldRight[i][j] = new Vector3d[znumber];
				additionalNewBfieldLeft[i][j] = new Vector3d[znumber];
				additionalNewBfieldRight[i][j] = new Vector3d[znumber];
				additionalChargeDensityHatLeft[i][j] = new double[znumber];
				additionalChargeDensityHatRight[i][j] = new double[znumber];
				additionalChargeDensityLeft[i][j] = new double[znumber];
				additionalChargeDensityMinusLeft[i][j] = new double[znumber];
				additionalChargeDensityRight[i][j] = new double[znumber];
				additionalChargeDensityMinusRight[i][j] = new double[znumber];
				additionalChargeDensityLeft[i][j] = new double[znumber];
				additionalChargeDensityRight[i][j] = new double[znumber];
				additionalPressureTensorLeft[i][j] = new Matrix3d[znumber];
				additionalPressureTensorRight[i][j] = new Matrix3d[znumber];
				for (int k = 0; k < znumber; ++k) {
					additionalBfieldLeft[i][j][k] = Vector3d(0, 0, 0);
					additionalBfieldRight[i][j][k] = Vector3d(0, 0, 0);
					additionalNewBfieldLeft[i][j][k] = Vector3d(0, 0, 0);
					additionalNewBfieldRight[i][j][k] = Vector3d(0, 0, 0);
					additionalChargeDensityHatLeft[i][j][k] = 0;
					additionalChargeDensityHatRight[i][j][k] = 0;
					additionalChargeDensityLeft[i][j][k] = 0;
					additionalChargeDensityMinusLeft[i][j][k] = 0;
					additionalChargeDensityRight[i][j][k] = 0;
					additionalChargeDensityMinusRight[i][j][k] = 0;
					additionalChargeDensityLeft[i][j][k] = 0;
					additionalChargeDensityRight[i][j][k] = 0;
					additionalPressureTensorLeft[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
					additionalPressureTensorRight[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				}
			}
		}
	}

	if (rank == 0) printf("creating arrays for fluxes\n");
	fflush(stdout);
	//if(rank == 0) printLog("creating arrays for fluxes\n");

	electricFlux = new Vector3d **[xnumber + 2];
	electricFluxMinus = new Vector3d **[xnumber + 2];
	dielectricTensor = new Matrix3d **[xnumber + 2];
	externalElectricFlux = new Vector3d **[xnumber + 2];
	divPressureTensor = new Vector3d **[xnumber + 2];

	for (int i = 0; i < xnumber + 2; ++i) {
		electricFlux[i] = new Vector3d *[ynumber + 1];
		electricFluxMinus[i] = new Vector3d *[ynumber + 1];
		dielectricTensor[i] = new Matrix3d *[ynumber + 1];
		externalElectricFlux[i] = new Vector3d *[ynumber + 1];
		divPressureTensor[i] = new Vector3d *[ynumber + 1];
		for (int j = 0; j < ynumber + 1; ++j) {
			electricFlux[i][j] = new Vector3d[znumber + 1];
			electricFluxMinus[i][j] = new Vector3d[znumber + 1];
			dielectricTensor[i][j] = new Matrix3d[znumber + 1];
			externalElectricFlux[i][j] = new Vector3d[znumber + 1];
			divPressureTensor[i][j] = new Vector3d[znumber + 1];
			for (int k = 0; k < znumber + 1; ++k) {
				//if(rank == 0) printf("%d %d %d\n", i, j, k);
				//if(rank == 0) fflush((stdout));
				electricFlux[i][j][k] = Vector3d(0, 0, 0);
				electricFluxMinus[i][j][k] = Vector3d(0, 0, 0);
				externalElectricFlux[i][j][k] = Vector3d(0, 0, 0);
				divPressureTensor[i][j][k] = Vector3d(0, 0, 0);
				dielectricTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
			}
		}
	}

	if (additionalBinNumber > 0) {
		additionalEfieldLeft = new Vector3d **[additionalBinNumber];
		additionalEfieldRight = new Vector3d **[additionalBinNumber];
		additionalTempEfieldLeft = new Vector3d **[additionalBinNumber];
		additionalTempEfieldRight = new Vector3d **[additionalBinNumber];
		additionalNewEfieldLeft = new Vector3d **[additionalBinNumber];
		additionalNewEfieldRight = new Vector3d **[additionalBinNumber];
		additionalElectricFluxLeft = new Vector3d **[additionalBinNumber];
		additionalElectricFluxMinusLeft = new Vector3d **[additionalBinNumber];
		additionalElectricFluxRight = new Vector3d **[additionalBinNumber];
		additionalElectricFluxMinusRight = new Vector3d **[additionalBinNumber];
		additionalDielectricTensorLeft = new Matrix3d **[additionalBinNumber];
		additionalDielectricTensorRight = new Matrix3d **[additionalBinNumber];
		additionalDivPressureTensorLeft = new Vector3d **[additionalBinNumber];
		additionalDivPressureTensorRight = new Vector3d **[additionalBinNumber];
		for (int i = 0; i < additionalBinNumber; ++i) {
			additionalEfieldLeft[i] = new Vector3d *[ynumber + 1];
			additionalEfieldRight[i] = new Vector3d *[ynumber + 1];
			additionalTempEfieldLeft[i] = new Vector3d *[ynumber + 1];
			additionalTempEfieldRight[i] = new Vector3d *[ynumber + 1];
			additionalNewEfieldLeft[i] = new Vector3d *[ynumber + 1];
			additionalNewEfieldRight[i] = new Vector3d *[ynumber + 1];
			additionalElectricFluxLeft[i] = new Vector3d *[ynumber + 1];
			additionalElectricFluxMinusLeft[i] = new Vector3d *[ynumber + 1];
			additionalElectricFluxRight[i] = new Vector3d *[ynumber + 1];
			additionalElectricFluxMinusRight[i] = new Vector3d *[ynumber + 1];
			additionalDielectricTensorLeft[i] = new Matrix3d *[ynumber + 1];
			additionalDielectricTensorRight[i] = new Matrix3d *[ynumber + 1];
			additionalDivPressureTensorLeft[i] = new Vector3d *[ynumber + 1];
			additionalDivPressureTensorRight[i] = new Vector3d *[ynumber + 1];
			for (int j = 0; j < ynumber + 1; ++j) {
				additionalEfieldLeft[i][j] = new Vector3d[znumber + 1];
				additionalEfieldRight[i][j] = new Vector3d[znumber + 1];
				additionalTempEfieldLeft[i][j] = new Vector3d[znumber + 1];
				additionalTempEfieldRight[i][j] = new Vector3d[znumber + 1];
				additionalNewEfieldLeft[i][j] = new Vector3d[znumber + 1];
				additionalNewEfieldRight[i][j] = new Vector3d[znumber + 1];
				additionalElectricFluxLeft[i][j] = new Vector3d[znumber + 1];
				additionalElectricFluxMinusLeft[i][j] = new Vector3d[znumber + 1];
				additionalElectricFluxRight[i][j] = new Vector3d[znumber + 1];
				additionalElectricFluxMinusRight[i][j] = new Vector3d[znumber + 1];
				additionalDielectricTensorLeft[i][j] = new Matrix3d[znumber + 1];
				additionalDielectricTensorRight[i][j] = new Matrix3d[znumber + 1];
				additionalDivPressureTensorLeft[i][j] = new Vector3d[znumber + 1];
				additionalDivPressureTensorRight[i][j] = new Vector3d[znumber + 1];
				for (int k = 0; k < znumber + 1; ++k) {
					additionalEfieldLeft[i][j][k] = Vector3d(0, 0, 0);
					additionalEfieldRight[i][j][k] = Vector3d(0, 0, 0);
					additionalTempEfieldLeft[i][j][k] = Vector3d(0, 0, 0);
					additionalTempEfieldRight[i][j][k] = Vector3d(0, 0, 0);
					additionalNewEfieldLeft[i][j][k] = Vector3d(0, 0, 0);
					additionalNewEfieldRight[i][j][k] = Vector3d(0, 0, 0);
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
	}

	if (rank == 0) printf("creating arrays for buffers\n");
	fflush(stdout);
	//if(rank == 0) printLog("creating arrays for buffers\n");

	rightEbuffer = new double[(ynumber + 1) * (znumber + 1) * 3];
	leftEbuffer = new double[(ynumber + 1) * (znumber + 1) * 3];
	leftEoutBuffer = new double[(ynumber + 1) * (znumber + 1) * 3];
	rightEinBuffer = new double[(ynumber + 1) * (znumber + 1) * 3];

	tempCellParameterLeft = new double **[2 + additionalBinNumber];
	tempCellParameterRight = new double **[2 + additionalBinNumber];
	tempNodeParameterLeft = new double **[2 + additionalBinNumber];
	tempNodeParameterRight = new double **[2 + additionalBinNumber];

	tempCellVectorParameterLeft = new Vector3d **[2 + additionalBinNumber];
	tempCellVectorParameterRight = new Vector3d **[2 + additionalBinNumber];
	tempNodeVectorParameterLeft = new Vector3d **[2 + additionalBinNumber];
	tempNodeVectorParameterRight = new Vector3d **[2 + additionalBinNumber];

	tempCellMatrixParameterLeft = new Matrix3d **[2 + additionalBinNumber];
	tempCellMatrixParameterRight = new Matrix3d **[2 + additionalBinNumber];
	tempNodeMatrixParameterLeft = new Matrix3d **[2 + additionalBinNumber];
	tempNodeMatrixParameterRight = new Matrix3d **[2 + additionalBinNumber];

	for (int i = 0; i < 2 + additionalBinNumber; ++i) {
		tempCellParameterLeft[i] = new double *[ynumber];
		tempCellParameterRight[i] = new double *[ynumber];
		tempCellVectorParameterLeft[i] = new Vector3d *[ynumber];
		tempCellVectorParameterRight[i] = new Vector3d *[ynumber];
		tempCellMatrixParameterLeft[i] = new Matrix3d *[ynumber];
		tempCellMatrixParameterRight[i] = new Matrix3d *[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			tempCellParameterLeft[i][j] = new double[znumber];
			tempCellParameterRight[i][j] = new double[znumber];
			tempCellVectorParameterLeft[i][j] = new Vector3d[znumber];
			tempCellVectorParameterRight[i][j] = new Vector3d[znumber];
			tempCellMatrixParameterLeft[i][j] = new Matrix3d[znumber];
			tempCellMatrixParameterRight[i][j] = new Matrix3d[znumber];
			for (int k = 0; k < znumber; ++k) {
				tempCellParameterLeft[i][j][k] = 0;
				tempCellParameterRight[i][j][k] = 0;
			}
		}
	}

	for (int i = 0; i < 2 + additionalBinNumber; ++i) {
		tempNodeParameterLeft[i] = new double *[ynumber + 1];
		tempNodeParameterRight[i] = new double *[ynumber + 1];
		tempNodeVectorParameterLeft[i] = new Vector3d *[ynumber + 1];
		tempNodeVectorParameterRight[i] = new Vector3d *[ynumber + 1];
		tempNodeMatrixParameterLeft[i] = new Matrix3d *[ynumber + 1];
		tempNodeMatrixParameterRight[i] = new Matrix3d *[ynumber + 1];
		for (int j = 0; j < ynumber + 1; ++j) {
			tempNodeParameterLeft[i][j] = new double[znumber + 1];
			tempNodeParameterRight[i][j] = new double[znumber + 1];
			tempNodeVectorParameterLeft[i][j] = new Vector3d[znumber + 1];
			tempNodeVectorParameterRight[i][j] = new Vector3d[znumber + 1];
			tempNodeMatrixParameterLeft[i][j] = new Matrix3d[znumber + 1];
			tempNodeMatrixParameterRight[i][j] = new Matrix3d[znumber + 1];
			for (int k = 0; k < znumber + 1; ++k) {
				tempNodeParameterLeft[i][j][k] = 0;
				tempNodeParameterRight[i][j][k] = 0;
			}
		}
	}


	arrayCreated = true;

	if (rank == 0) printf("finish creating arrays\n");
	fflush(stdout);
	// if(rank == 0) printLog("finish creating arrays\n");
}

//add new particle types here
void Simulation::createParticleTypes(double* concentrations, int* particlesPerBin) {
	types = new ParticleTypeContainer[typesNumber];
	ParticleTypeContainer type;

	type.type = ELECTRON;
	type.number = 0;
	type.mass = massElectron;
	type.chargeCount = -1;
	type.charge = -electron_charge_normalized;
	type.particlesPerBin = particlesPerBin[0];
	type.concentration = concentrations[0];
	type.typeName = "electron";

	types[0] = type;

	type.type = PROTON;
	type.number = 1;
	type.mass = massProton;
	type.chargeCount = 1;
	type.charge = electron_charge_normalized;
	type.particlesPerBin = particlesPerBin[1];
	type.concentration = concentrations[1];
	type.typeName = "proton";

	types[1] = type;

	type.type = POSITRON;
	type.number = 2;
	type.mass = massElectron;
	type.chargeCount = 1;
	type.charge = electron_charge_normalized;
	type.particlesPerBin = particlesPerBin[2];
	type.concentration = concentrations[2];
	type.typeName = "positron";

	types[2] = type;

	type.type = ALPHA;
	type.number = 3;
	type.mass = massAlpha;
	type.chargeCount = 2;
	type.charge = 2 * electron_charge_normalized;
	type.particlesPerBin = particlesPerBin[3];
	type.concentration = concentrations[3];
	type.typeName = "alpha";

	types[3] = type;

	type.type = DEUTERIUM;
	type.number = 4;
	type.mass = massDeuterium;
	type.chargeCount = 1;
	type.charge = electron_charge_normalized;
	type.particlesPerBin = particlesPerBin[4];
	type.concentration = concentrations[4];
	type.typeName = "deuterium";

	types[4] = type;

	type.type = HELIUM3;
	type.number = 5;
	type.mass = massHelium3;
	type.chargeCount = 2;
	type.charge = 2 * electron_charge_normalized;
	type.particlesPerBin = particlesPerBin[5];
	type.concentration = concentrations[5];
	type.typeName = "helium3";

	types[5] = type;

	type.type = OXYGEN_PLUS3;
	type.number = 6;
	type.mass = massOxygen;
	type.chargeCount = 3;
	type.charge = 3 * electron_charge_normalized;
	type.particlesPerBin = particlesPerBin[6];
	type.concentration = concentrations[6];
	type.typeName = "oxygen_plus_3";

	types[6] = type;

	type.type = SILICON_PLUS1;
	type.number = 7;
	type.mass = massSilicon;
	type.chargeCount = 1;
	type.charge = electron_charge_normalized;
	type.particlesPerBin = particlesPerBin[7];
	type.concentration = concentrations[7];
	type.typeName = "silicon_plus_1";

	types[7] = type;

	for (int i = 0; i < typesNumber; ++i) {
		if (types[i].particlesPerBin > 0) {
			types[i].particesDeltaX = deltaX / types[i].particlesPerBin;
			types[i].particesDeltaY = deltaY / types[i].particlesPerBin;
			types[i].particesDeltaZ = deltaZ / types[i].particlesPerBin;
		}
		else {
			types[i].particesDeltaX = xsize;
			types[i].particesDeltaY = ysize;
			types[i].particesDeltaZ = zsize;
		}
		types[i].temperatureX = temperature;
		types[i].temperatureY = temperature;
		types[i].temperatureZ = temperature;
		types[i].anisotropy = 0;
		types[i].generalWeight = 0;
	}

	if (rank == 0) {
		particleTypesFile = fopen((outputDir + "particleTypes.dat").c_str(), "w");
		for (int t = 0; t < typesNumber; ++t) {
			fprintf(particleTypesFile, "%d\n", types[t].particlesPerBin);
		}
		fclose(particleTypesFile);
	}
}

void Simulation::createFiles() {
	if (rank == 0) {
		printf("creating files\n");
		fflush(stdout);
		FILE* logFile = fopen((outputDir + "log.dat").c_str(), "w");
		fflush(logFile);
		fclose(logFile);
		printLog("creatingFiles\n");
		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", 0.0, 0.0, 0.0, 0.0);
		fclose(incrementFile);
		electronTraectoryFile = fopen((outputDir + "trajectory_electron.dat").c_str(), "w");
		fclose(electronTraectoryFile);
		electronTraectoryFile = fopen((outputDir + "trajectory_electron_1.dat").c_str(), "w");
		fclose(electronTraectoryFile);
		electronTraectoryFile = fopen((outputDir + "trajectory_electron_2.dat").c_str(), "w");
		fclose(electronTraectoryFile);
		electronTraectoryFile = fopen((outputDir + "trajectory_electron_3.dat").c_str(), "w");
		fclose(electronTraectoryFile);
		electronTraectoryFile = fopen((outputDir + "trajectory_electron_4.dat").c_str(), "w");
		fclose(electronTraectoryFile);
		electronTraectoryFile = fopen((outputDir + "trajectory_electron_5.dat").c_str(), "w");
		fclose(electronTraectoryFile);
		electronTraectoryFile = fopen((outputDir + "trajectory_electron_6.dat").c_str(), "w");
		fclose(electronTraectoryFile);
		electronTraectoryFile = fopen((outputDir + "trajectory_electron_7.dat").c_str(), "w");
		fclose(electronTraectoryFile);
		electronTraectoryFile = fopen((outputDir + "trajectory_electron_8.dat").c_str(), "w");
		fclose(electronTraectoryFile);
		electronTraectoryFile = fopen((outputDir + "trajectory_electron_9.dat").c_str(), "w");
		fclose(electronTraectoryFile);
		protonTraectoryFile = fopen((outputDir + "trajectory_proton.dat").c_str(), "w");
		fclose(protonTraectoryFile);
		protonTraectoryFile = fopen((outputDir + "trajectory_proton_1.dat").c_str(), "w");
		fclose(protonTraectoryFile);
		protonTraectoryFile = fopen((outputDir + "trajectory_proton_2.dat").c_str(), "w");
		fclose(protonTraectoryFile);
		protonTraectoryFile = fopen((outputDir + "trajectory_proton_3.dat").c_str(), "w");
		fclose(protonTraectoryFile);
		protonTraectoryFile = fopen((outputDir + "trajectory_proton_4.dat").c_str(), "w");
		fclose(protonTraectoryFile);
		protonTraectoryFile = fopen((outputDir + "trajectory_proton_5.dat").c_str(), "w");
		fclose(protonTraectoryFile);
		protonTraectoryFile = fopen((outputDir + "trajectory_proton_6.dat").c_str(), "w");
		fclose(protonTraectoryFile);
		protonTraectoryFile = fopen((outputDir + "trajectory_proton_7.dat").c_str(), "w");
		fclose(protonTraectoryFile);
		protonTraectoryFile = fopen((outputDir + "trajectory_proton_8.dat").c_str(), "w");
		fclose(protonTraectoryFile);
		protonTraectoryFile = fopen((outputDir + "trajectory_proton_9.dat").c_str(), "w");
		fclose(protonTraectoryFile);

		distributionFileProton = fopen((outputDir + "distribution_protons.dat").c_str(), "w");
		fclose(distributionFileProton);
		distributionFileElectron = fopen((outputDir + "distribution_electrons.dat").c_str(), "w");
		fclose(distributionFileElectron);
		distributionFileAlpha = fopen((outputDir + "distribution_alphas.dat").c_str(), "w");
		fclose(distributionFileAlpha);
		distributionFilePositron = fopen((outputDir + "distribution_positrons.dat").c_str(), "w");
		fclose(distributionFilePositron);
		distributionFileHelium3 = fopen((outputDir + "distribution_helium3.dat").c_str(), "w");
		fclose(distributionFileHelium3);
		distributionFileDeuterium = fopen((outputDir + "distribution_deuterium.dat").c_str(), "w");
		fclose(distributionFileDeuterium);
		distributionFileOxygen = fopen((outputDir + "distribution_oxygen.dat").c_str(), "w");
		fclose(distributionFileOxygen);
		distributionFileSilicon = fopen((outputDir + "distribution_silicon.dat").c_str(), "w");
		fclose(distributionFileSilicon);

		distributionFileProton = fopen((outputDir + "distribution_protons_sw.dat").c_str(), "w");
		fclose(distributionFileProton);
		distributionFileElectron = fopen((outputDir + "distribution_electrons_sw.dat").c_str(), "w");
		fclose(distributionFileElectron);
		distributionFileAlpha = fopen((outputDir + "distribution_alphas_sw.dat").c_str(), "w");
		fclose(distributionFileAlpha);
		distributionFilePositron = fopen((outputDir + "distribution_positrons_sw.dat").c_str(), "w");
		fclose(distributionFilePositron);
		distributionFileHelium3 = fopen((outputDir + "distribution_helium3_sw.dat").c_str(), "w");
		fclose(distributionFileHelium3);
		distributionFileDeuterium = fopen((outputDir + "distribution_deuterium_sw.dat").c_str(), "w");
		fclose(distributionFileDeuterium);
		distributionFileOxygen = fopen((outputDir + "distribution_oxygen_sw.dat").c_str(), "w");
		fclose(distributionFileOxygen);
		distributionFileSilicon = fopen((outputDir + "distribution_silicon_sw.dat").c_str(), "w");
		fclose(distributionFileSilicon);
		
		anisotropyFileElectron = fopen((outputDir + "anisotropy_electrons.dat").c_str(), "w");
		fclose(anisotropyFileElectron);
		anisotropyFileProton = fopen((outputDir + "anisotropy_protons.dat").c_str(), "w");
		fclose(anisotropyFileProton);
		anisotropyFileAlpha = fopen((outputDir + "anisotropy_alphas.dat").c_str(), "w");
		fclose(anisotropyFileAlpha);
		anisotropyFilePositron = fopen((outputDir + "anisotropy_positrons.dat").c_str(), "w");
		fclose(anisotropyFilePositron);
		anisotropyFileDeuterium = fopen((outputDir + "anisotropy_deuterium.dat").c_str(), "w");
		fclose(anisotropyFileDeuterium);
		anisotropyFileHelium3 = fopen((outputDir + "anisotropy_helium3.dat").c_str(), "w");
		fclose(anisotropyFileHelium3);
		anisotropyFileHelium3 = fopen((outputDir + "anisotropy_oxygen+3.dat").c_str(), "w");
		fclose(anisotropyFileHelium3);
		anisotropyFileHelium3 = fopen((outputDir + "anisotropy_silicon.dat").c_str(), "w");
		fclose(anisotropyFileHelium3);
		EfieldFile = fopen((outputDir + "Efield.dat").c_str(), "w");
		fclose(EfieldFile);
		BfieldFile = fopen((outputDir + "Bfield.dat").c_str(), "w");
		fclose(BfieldFile);
		velocityFile = fopen((outputDir + "velocity.dat").c_str(), "w");
		fclose(velocityFile);
		Xfile = fopen((outputDir + "Xfile.dat").c_str(), "w");
		fclose(Xfile);
		Yfile = fopen((outputDir + "Yfile.dat").c_str(), "w");
		fclose(Yfile);
		Zfile = fopen((outputDir + "Zfile.dat").c_str(), "w");
		fclose(Zfile);
		generalFile = fopen((outputDir + "general.dat").c_str(), "w");
		fclose(generalFile);
		generalAnisotropyFile = fopen((outputDir + "generalAnisotropy.dat").c_str(), "w");
		fclose(generalAnisotropyFile);
		densityFile = fopen((outputDir + "concentrations.dat").c_str(), "w");
		fclose(densityFile);
		divergenceErrorFile = fopen((outputDir + "divergence_error.dat").c_str(), "w");
		fclose(divergenceErrorFile);
		informationFile = fopen((outputDir + "information.dat").c_str(), "w");
		fclose(informationFile);
		fluxFile = fopen((outputDir + "flux.dat").c_str(), "w");
		fclose(fluxFile);
		rotBFile = fopen((outputDir + "rotBFile.dat").c_str(), "w");
		fclose(rotBFile);
		rotEFile = fopen((outputDir + "rotEFile.dat").c_str(), "w");
		fclose(rotBFile);
		EderivativeFile = fopen((outputDir + "EderivativeFile.dat").c_str(), "w");
		fclose(EderivativeFile);
		dielectricTensorFile = fopen((outputDir + "dielectricTensorFile.dat").c_str(), "w");
		fclose(dielectricTensorFile);
		maxwellMatrixFile = fopen((outputDir + "maxwellMatrixFile.dat").c_str(), "w");
		fclose(maxwellMatrixFile);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fclose(errorLogFile);
		particleProtonsFile = fopen((outputDir + "proton.dat").c_str(), "w");
		fclose(particleProtonsFile);
		particleElectronsFile = fopen((outputDir + "electron.dat").c_str(), "w");
		fclose(particleElectronsFile);
		particlePositronsFile = fopen((outputDir + "positron.dat").c_str(), "w");
		fclose(particlePositronsFile);
		particleAlphaFile = fopen((outputDir + "alpha.dat").c_str(), "w");
		fclose(particleAlphaFile);
		particleDeuteriumFile = fopen((outputDir + "deuterium.dat").c_str(), "w");
		fclose(particleDeuteriumFile);
		particleHelium3File = fopen((outputDir + "helium3.dat").c_str(), "w");
		fclose(particleHelium3File);
		//outputEverythingFile = fopen("./output/everything.dat","w");
		//fclose(outputEverythingFile);
	}
}

void Simulation::checkFrequency(double omega) {
	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");
	double cyclothronOmegaElectron = electron_charge_normalized * B0.norm() / (massElectron * speed_of_light_normalized);
	double cyclothronOmegaProton = electron_charge_normalized * B0.norm() / (massProton * speed_of_light_normalized);
	if (omega > cyclothronOmegaProton) {
		if (rank == 0) printf("omega > cyclothron Omega Proton\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega > cyclothron Omega Proton\n");
	}
	else if (omega > cyclothronOmegaProton / 100.0) {
		if (rank == 0) printf("omega > cyclothrone Omega Proton/100\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega > cyclothron Omega Proton/100\n");
	}
	if (rank == 0) printf("omega/cyclothronOmega = %g\n", omega / cyclothronOmegaProton);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "omega/cyclothronOmega = %g\n", omega / cyclothronOmegaProton);

	if (omega > 1.0) {
		if (rank == 0) printf("omega > omega plasma\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega > omega plasma\n");
	}
	else if (omega > 0.01) {
		if (rank == 0) printf("omega > omega plasma/100\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "omega > omega plasma/100\n");
	}
	if (rank == 0) printf("omega/omega plasma = %g\n", omega);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "omega/omega plasma = %g\n", omega);

	if (rank == 0) fclose(informationFile);
}

void Simulation::checkDebyeParameter() {
	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");


	double debyeLength2 = 0;
	for (int i = 0; i < typesNumber; ++i) {
		if (types[i].concentration > 0) {
			debyeLength2 += (4 * pi * types[i].concentration * sqr(types[i].charge)) / (kBoltzman_normalized * types[i].temperatureX);
		}
	}

	double debyeLength = 1 / sqrt(debyeLength2);

	double debyeNumber = 0;
	for (int i = 0; i < typesNumber; ++i) {
		debyeNumber += 4 * pi * cube(debyeLength) * types[i].concentration / 3;
	}

	if (debyeLength < deltaX) {
		if (rank == 0) printf("debye length < deltaX\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "debye length < deltaX\n");
	}
	if (rank == 0) printf("debye length/deltaX = %g\n", debyeLength / deltaX);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "debye length/deltaX = %g\n", debyeLength / deltaX);

	if (debyeNumber < 1.0) {
		if (rank == 0) printf("debye number < 1\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "debye number < 1\n");
	}
	else if (debyeNumber < 100.0) {
		if (rank == 0) printf("debye number < 100\n");
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "debye number < 100\n");
	}
	if (rank == 0) printf("debye number = %g\n", debyeNumber);
	fflush(stdout);
	if (rank == 0) fprintf(informationFile, "debye number = %g\n", debyeNumber);

	/*double superParticleDebyeLength = 1 / sqrt(
		4 * pi * superParticleCharge * superParticleCharge * superParticleConcentration / (kBoltzman_normalized * superParticleTemperature));
	double superParticleDebyeNumber = 4 * pi * cube(superParticleDebyeLength) * superParticleConcentration / 3;

	if (superParticleDebyeLength > deltaX) {
		if(rank == 0) printf("super particle debye length > deltaX\n");
		fflush(stdout);
		if(rank == 0) fprintf(informationFile, "super particle debye length > deltaX\n");
	}
	if(rank == 0) printf("super particle debye length/deltaX = %g\n", superParticleDebyeLength / deltaX);
	fflush(stdout);
	if(rank == 0) fprintf(informationFile, "super particle debye length/deltaX = %g\n", superParticleDebyeLength / deltaX);

	if (superParticleDebyeNumber < 1.0) {
		if(rank == 0) printf("superparticle debye number < 1\n");
		fflush(stdout);
		if(rank == 0) fprintf(informationFile, "superparticle debye number < 1\n");
	} else if (superParticleDebyeNumber < 100.0) {
		if(rank == 0) printf("superparticle debye number < 100\n");
		fflush(stdout);
		if(rank == 0) fprintf(informationFile, "superparticle debye number < 100\n");
	}
	if(rank == 0) printf("superparticle debye number = %g\n", superParticleDebyeNumber);
	fflush(stdout);
	if(rank == 0) fprintf(informationFile, "superparticle debye number = %g\n", superParticleDebyeNumber);*/
	if (rank == 0) fclose(informationFile);
}

void Simulation::checkGyroRadius() {
	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");

	if (B0.norm() > 0) {
		double thermalMomentumElectron = sqrt(
			massElectron * kBoltzman_normalized * temperature) + massElectron * V0.norm();
		double gyroRadiusElectron = thermalMomentumElectron * speed_of_light_normalized / (electron_charge_normalized * B0.norm());
		double thermalMomentumProton = sqrt(massProton * kBoltzman_normalized * temperature) + massProton * V0.norm();
		double gyroRadiusProton = thermalMomentumProton * speed_of_light_normalized / (electron_charge_normalized * B0.norm());
		if (deltaX > 0.5 * gyroRadiusElectron) {
			if (rank == 0) printf("deltaX > 0.5*gyroRadiusElectron\n");
			fflush(stdout);
			if (rank == 0) fprintf(informationFile, "deltaX > 0.5*gyroRadiusElectron\n");
		}

		if (rank == 0) printf("deltaX/gyroRadiusElectron = %g\n", deltaX / gyroRadiusElectron);
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "deltaX/gyroRadiusElectron = %g\n", deltaX / gyroRadiusElectron);

		if (xsize < 2 * gyroRadiusProton) {
			if (rank == 0) printf("xsize < 2*gyroRadiusProton\n");
			fflush(stdout);
			if (rank == 0) fprintf(informationFile, "xsize < 2*gyroRadiusProton\n");
		}

		if (rank == 0) printf("xsize/gyroRadiusProton= %g\n", xsize / gyroRadiusProton);
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "xsize/gyroRadiusProton = %g\n", xsize / gyroRadiusProton);
	}

	if (rank == 0) fclose(informationFile);
}

void Simulation::createParticles() {
	evaluateParticleTypesAlpha();
	if (rank == 0) printf("creating particles\n");
	fflush(stdout);
	if (rank == 0) printLog("creating particles\n");
	int n = 0;
	//for (int i = 0; i < xnumber; ++i) {
	for (int i = 1; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				//int maxParticlesPerBin = types[0].particlesPerBin;
				double x = xgrid[i] + 0.0001 * deltaX;
				double y = ygrid[j] + 0.0001 * deltaY;
				double z = zgrid[k] + 0.0001 * deltaZ;
				//for (int l = 0; l < maxParticlesPerBin; ++l) {
				for (int typeCounter = 0; typeCounter < typesNumber; ++typeCounter) {
					double weight = (types[typeCounter].concentration / types[typeCounter].particlesPerBin) * volumeB(
						i, j, k);
					double deltaXParticles = types[typeCounter].particesDeltaX;
					double deltaYParticles = types[typeCounter].particesDeltaY;
					double deltaZParticles = types[typeCounter].particesDeltaZ;
					//if (l < types[typeCounter].particlesPerBin) {
					for (int l = 0; l < types[typeCounter].particlesPerBin; ++l) {
						ParticleTypes type = types[typeCounter].type;
						Particle* particle = createParticle(n, i, j, k, weight, type, types[typeCounter],
						                                    types[typeCounter].temperatureX,
						                                    types[typeCounter].temperatureY,
						                                    types[typeCounter].temperatureZ);
						n++;
						particle->coordinates.x = x + deltaXParticles * l;
						particle->coordinates.y = y + deltaYParticles * l;
						particle->coordinates.z = z + deltaZParticles * l;
						particle->initialCoordinates = particle->coordinates;
						if (V0.norm() > 0) {
							particle->addVelocity(V0, speed_of_light_normalized);
						}
						Vector3d momentum = particle->getMomentum();
						particle->initialMomentum = momentum;
						//particle->prevMomentum = momentum;
						particles.push_back(particle);
						particlesNumber++;
						if (particlesNumber % 1000 == 0) {
							if ((rank == 0) && (verbosity > 0))printf("create particle number %d\n", particlesNumber);
						}
						alertNaNOrInfinity(particle->coordinates.x,"particle.x = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.y,"particle.y = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.z,"particle.z = NaN in createParticles\n");
					}
				}
			}
		}
	}


	synchronizeParticleNumber();

	protonNumber = getParticleNumber(0, PROTON);
	protonNumber1 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[1].particlesPerBin * 0.05, PROTON);
	protonNumber2 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[1].particlesPerBin * 0.10, PROTON);
	protonNumber3 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[1].particlesPerBin * 0.15, PROTON);
	protonNumber4 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[1].particlesPerBin * 0.20, PROTON);
	protonNumber5 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[1].particlesPerBin * 0.25, PROTON);
	protonNumber6 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[1].particlesPerBin * 0.30, PROTON);
	protonNumber7 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[1].particlesPerBin * 0.35, PROTON);
	protonNumber8 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[1].particlesPerBin * 0.40, PROTON);
	protonNumber9 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[1].particlesPerBin * 0.45, PROTON);

	electronNumber = getParticleNumber(0, ELECTRON);
	electronNumber1 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[0].particlesPerBin * 0.05, ELECTRON);
	electronNumber2 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[0].particlesPerBin * 0.10, ELECTRON);
	electronNumber3 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[0].particlesPerBin * 0.15, ELECTRON);
	electronNumber4 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[0].particlesPerBin * 0.20, ELECTRON);
	electronNumber5 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[0].particlesPerBin * 0.25, ELECTRON);
	electronNumber6 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[0].particlesPerBin * 0.30, ELECTRON);
	electronNumber7 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[0].particlesPerBin * 0.35, ELECTRON);
	electronNumber8 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[0].particlesPerBin * 0.40, ELECTRON);
	electronNumber9 = getParticleNumber(xnumberGeneral * ynumberGeneral * znumberGeneral * types[0].particlesPerBin * 0.45, ELECTRON);

	printf("rank = %d p0 = %d p1 = %d p2 = %d e0 = %d e1 = %d e2 = %d Np = %d\n", rank, protonNumber, protonNumber1, protonNumber2, electronNumber, electronNumber1, electronNumber2, particlesNumber);

	/*if (preserveChargeLocal) {
		moveToPreserveChargeLocal();
	}*/
}

void Simulation::moveToPreserveChargeLocal() {
	Particle* electron = NULL;
	Particle* notElectron = NULL;
	Particle* particle;
	int electronCount = 0;
	int notElectronCount = 0;
	int particleCount = 0;

	while (particleCount < particles.size()) {
		particle = particles[particleCount];
		if (particle->type != ELECTRON) {
			notElectron = particle;
			notElectronCount = particleCount;

			int necessaryElectrons = 1;

			if (particle->type == ALPHA) {
				necessaryElectrons = 2;
			}

			while (necessaryElectrons > 0) {
				electron = NULL;
				while (electron == NULL) {
					if (electronCount >= particles.size()) {
						errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
						printf("error in preserving charge\n");
						fflush(stdout);
						fprintf(errorLogFile, "error in preserving charge\n");
						fclose(errorLogFile);
						MPI_Finalize();
						exit(0);
					}
					if (particles[electronCount]->type == ELECTRON) {
						electron = particles[electronCount];
						electron->coordinates = notElectron->coordinates;
					}
					electronCount++;
				}
				necessaryElectrons--;
			}
		}
		particleCount++;
	}
}

void Simulation::addToPreserveChargeGlobal() {
	if (rank == nprocs - 1) {
		for (int i = 0; i < escapedParticlesRight.size(); ++i) {
			Particle* particle = escapedParticlesRight[i];
			int typeNumber = getTypeNumber(particle);
			ParticleTypeContainer type = types[typeNumber];

			Particle* newParticle = createParticle(particlesNumber, xnumber, 0, 0, particle->weight, type.type, type, type.temperatureX, type.temperatureY, type.temperatureZ);
			newParticle->coordinates.x = xgrid[xnumber] - 0.5 * deltaX * splineOrder;
			newParticle->coordinates.y = particle->coordinates.y;
			newParticle->coordinates.z = particle->coordinates.z;
			newParticle->initialCoordinates = newParticle->coordinates;
			newParticle->addVelocity(V0, speed_of_light_normalized);
			Vector3d momentum = particle->getMomentum();
			newParticle->initialMomentum = momentum;
			//newParticle->prevMomentum = momentum;
			particles.push_back(newParticle);
			theoreticalEnergy += newParticle->energy(speed_of_light_normalized) * newParticle->weight * 
				sqr(scaleFactor / plasma_period);
			theoreticalMomentum += newParticle->getMomentum() * newParticle->weight * scaleFactor / plasma_period;
			particlesNumber++;
		}
	}
	if ((rank == 0) && (verbosity > 1)) printf("exchanging particles number\n");
	MPI_Barrier(MPI_COMM_WORLD);
	int tempParticleNumber[1];
	tempParticleNumber[0] = particlesNumber;
	MPI_Bcast(tempParticleNumber, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD);
	particlesNumber = tempParticleNumber[0];
}

Particle* Simulation::getFirstProton() {
	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		if (particle->type == PROTON) {
			return particle;
		}
	}
	return NULL;
}

Particle* Simulation::getFirstElectron() {
	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		if (particle->type == ELECTRON) {
			return particle;
		}
	}
	return NULL;
}

Particle* Simulation::getLastProton() {
	for (int pcount = particles.size() - 1; pcount >= 0; --pcount) {
		Particle* particle = particles[pcount];
		if (particle->type == PROTON) {
			return particle;
		}
	}
	return NULL;
}

Particle* Simulation::getLastElectron() {
	for (int pcount = particles.size() - 1; pcount >= 0; --pcount) {
		Particle* particle = particles[pcount];
		if (particle->type == ELECTRON) {
			return particle;
		}
	}
	return NULL;
}

Particle* Simulation::getProton(int n) {
	int count = 0;
	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		if (particle->type == PROTON) {
			count++;
			if (count == n) {
				return particle;
			}
		}
	}
	errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
	fprintf(errorLogFile, "can not find proton number %d\n", n);
	printf("can not find proton number %d\n", n);
	fflush(stdout);
	fclose(errorLogFile);
	MPI_Finalize();
	exit(0);
	return NULL;
}

Particle* Simulation::getElectron(int n) {
	int count = 0;
	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		if (particle->type == ELECTRON) {
			count++;
			if (count == n) {
				return particle;
			}
		}
	}
	errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
	fprintf(errorLogFile, "can not find electron number %d\n", n);
	printf("can not find electron number %d\n", n);
	fflush(stdout);
	fclose(errorLogFile);
	MPI_Finalize();
	exit(0);
	return NULL;
}

int Simulation::getParticleNumber(int n, const ParticleTypes& type) {
	int tempParticleNumber[1];
	int curParticleNumber[1];
	int particleNumber = -1;
	tempParticleNumber[0] = -1;
	curParticleNumber[0] = 0;

	if (rank > 0) {
		MPI_Status status;
		MPI_Recv(curParticleNumber, 1, MPI_INT, rank - 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD, &status);
	}

	if (curParticleNumber[0] > -1) {
		for (int p = 0; p < particles.size(); ++p) {
			Particle* particle = particles[p];
			if (particle->type == type) {
				if (curParticleNumber[0] == n) {
					tempParticleNumber[0] = particle->number;
					curParticleNumber[0] = -1;
					break;
				}
				curParticleNumber[0] = curParticleNumber[0] + 1;

			}
		}
	}

	if (rank < nprocs - 1) {
		MPI_Send(curParticleNumber, 1, MPI_INT, rank + 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (rank > 0) {
		MPI_Send(tempParticleNumber, 1, MPI_INT, 0, MPI_SEND_INTEGER_ALL_TO_FIRST, MPI_COMM_WORLD);
	}
	else {
		particleNumber = tempParticleNumber[0];
		for (int i = 1; i < nprocs; ++i) {
			MPI_Status status;
			MPI_Recv(tempParticleNumber, 1, MPI_INT, i, MPI_SEND_INTEGER_ALL_TO_FIRST, MPI_COMM_WORLD, &status);
			if (tempParticleNumber[0] >= 0) {
				particleNumber = tempParticleNumber[0];
			}
		}
	}

	tempParticleNumber[0] = particleNumber;

	MPI_Bcast(tempParticleNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);

	particleNumber = tempParticleNumber[0];

	return particleNumber;
}


Particle* Simulation::createParticle(int n, int i, int j, int k, const double& weight,
                                     const ParticleTypes& type, const ParticleTypeContainer& typeContainer, const double& localTemperatureX,
                                     const double& localTemperatureY, const double& localTemperatureZ) {
	int chargeCount = typeContainer.chargeCount;
	double charge = typeContainer.charge;
	double mass = typeContainer.mass;

	double x = xgrid[i] + deltaX * uniformDistribution();
	double y = ygrid[j] + deltaY * uniformDistribution();
	double z = zgrid[k] + deltaZ * uniformDistribution();

	double dx = (deltaX / 2) * (splineOrder + 1);
	double dy = (deltaY / 2) * (splineOrder + 1);
	double dz = (deltaZ / 2) * (splineOrder + 1);

	/*double dx = (deltaX / 2);
	double dy = (deltaY / 2);
	double dz = (deltaZ / 2);*/

	double energy = mass * speed_of_light_normalized_sqr;
	double p = 0;
	double px, py, pz = 0;

	double thetaParamter = kBoltzman_normalized * (localTemperatureX + localTemperatureY + localTemperatureZ) / (3 * mass * speed_of_light_normalized_sqr);

	if (thetaParamter < 0.01) {
		//energy = maxwellDistribution(localTemparature, kBoltzman_normalized);
		px = sqrt(mass * kBoltzman_normalized * localTemperatureX) * normalDistribution();
		py = sqrt(mass * kBoltzman_normalized * localTemperatureY) * normalDistribution();
		pz = sqrt(mass * kBoltzman_normalized * localTemperatureZ) * normalDistribution();
		p = sqrt(px * px + py * py + pz * pz);
	}
	else if (localTemperatureX == localTemperatureY && localTemperatureX == localTemperatureZ) {
		energy = maxwellJuttnerDistribution(localTemperatureX, mass, speed_of_light_normalized, kBoltzman_normalized);
		p = sqrt(energy * energy - sqr(mass * speed_of_light_normalized_sqr)) / speed_of_light_normalized;


		//p = 0;

		pz = p * (2 * uniformDistribution() - 1);
		double phi = 2 * pi * uniformDistribution();
		double pnormal = sqrt(p * p - pz * pz);
		px = pnormal * cos(phi);
		py = pnormal * sin(phi);
	}
	else if (localTemperatureY == localTemperatureZ) {
		double momentumParallel;
		double momentumNormal;
		anisotropicMaxwellJuttnerDistribution(momentumNormal, momentumParallel, types[typeContainer.number].alphaNormal,
		                                      types[typeContainer.number].alphaParallel,
		                                      mass * speed_of_light_normalized);


		double sign = uniformDistribution() - 0.5;
		double phi = 2 * pi * uniformDistribution();
		if (sign > 0) {
			px = momentumParallel;
		}
		else {
			px = -momentumParallel;
		}
		py = momentumNormal * cos(phi);
		pz = momentumNormal * sin(phi);
		alertNaNOrInfinity(px, "px = NaN\n");
		alertNaNOrInfinity(py, "py = NaN\n");
		alertNaNOrInfinity(pz, "pz = NaN\n");
	}
	else {
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "can not find electron number %d\n", n);
		printf("can not find electron number %d\n", n);
		fflush(stdout);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}


	Particle* particle = new Particle(n, mass, chargeCount, charge, weight, type, x, y, z, px, py, pz, dx, dy, dz);

	return particle;
}

int Simulation::getTypeNumber(Particle* particle) {
	for (int t = 0; t < typesNumber; ++t) {
		if (types[t].type == particle->type) {
			return t;
		}
	}
	errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
	fprintf(errorLogFile, "particle has no type\n");
	printf("particle has no type\n");
	fflush(stdout);
	fclose(errorLogFile);
	MPI_Finalize();
	exit(0);

	return -1;
}

void Simulation::evaluateParticleTypesAlpha() {
	double* alphas = new double[2 * typesNumber];
	if (rank == 0) {
		for (int i = 0; i < typesNumber; ++i) {
			if (types[i].particlesPerBin > 0) {
				double alphaNormal = types[i].mass * speed_of_light_normalized_sqr / (kBoltzman_normalized * types[i].temperatureY);
				TemperatureRelativisticMaxwellSolver solver = TemperatureRelativisticMaxwellSolver(alphaNormal,
				                                                                                   kBoltzman_normalized * types[i].temperatureX / (types[i].mass * speed_of_light_normalized_sqr));
				double tempValue = solver.function(alphaNormal);
				double minAlpha;
				double maxAlpha;
				if (tempValue > 0) {
					minAlpha = alphaNormal;
					maxAlpha = 2 * alphaNormal;
					while (solver.function(maxAlpha) > 0) {
						minAlpha = maxAlpha;
						maxAlpha = 2 * maxAlpha;
					}
				}
				else {
					maxAlpha = alphaNormal;
					minAlpha = 0.1 * types[i].mass * speed_of_light_normalized_sqr / (kBoltzman_normalized * types[i].temperatureX);
				}
				double simpleAlphaParallel = types[i].mass * speed_of_light_normalized_sqr / (kBoltzman_normalized * types[i].temperatureX);
				double alphaParallel = solver.solve(minAlpha, maxAlpha);
				types[i].alphaNormal = alphaNormal;
				types[i].alphaParallel = alphaParallel;
			}
			else {
				types[i].alphaNormal = 1.0;
				types[i].alphaParallel = 1.0;
			}
			alphas[2 * i] = types[i].alphaNormal;
			alphas[2 * i + 1] = types[i].alphaParallel;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(alphas, 2 * typesNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < typesNumber; ++i) {
		types[i].alphaNormal = alphas[2 * i];
		types[i].alphaParallel = alphas[2 * i + 1];
	}
	MPI_Barrier(MPI_COMM_WORLD);

	delete[] alphas;
}

void Simulation::synchronizeParticleNumber() {
	if (nprocs > 1) {
		int particleCount[1];
		particleCount[0] = 0;
		if (rank > 0) {
			particleCount[0] = particles.size();
			MPI_Send(particleCount, 1, MPI_INT, 0, MPI_SEND_INTEGER_ALL_TO_FIRST, MPI_COMM_WORLD);
		}
		else {
			particlesNumber = particles.size();
			for (int i = 1; i < nprocs; ++i) {
				MPI_Status status;
				MPI_Recv(particleCount, 1, MPI_INT, i, MPI_SEND_INTEGER_ALL_TO_FIRST, MPI_COMM_WORLD, &status);
				particlesNumber += particleCount[0];
			}
		}
		particleCount[0] = particlesNumber;
		MPI_Bcast(particleCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
		particlesNumber = particleCount[0];

		particleCount[0] = 0;

		if (rank > 0) {
			MPI_Status status;
			MPI_Recv(particleCount, 1, MPI_INT, rank - 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD, &status);
		}
		for (int p = 0; p < particles.size(); ++p) {
			Particle* particle = particles[p];
			particle->number = particleCount[0] + p;
		}
		particleCount[0] = particleCount[0] + particles.size();
		if (rank < nprocs - 1) {
			MPI_Send(particleCount, 1, MPI_INT, rank + 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD);
		}
	}
	else {
		particlesNumber = particles.size();
		for (int p = 0; p < particles.size(); ++p) {
			Particle* particle = particles[p];
			particle->number = p;
		}
	}
}

