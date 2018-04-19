#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>
#include <mpi.h>
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
#include "mpi_util.h"
#include "input.h"
#include "complex.h"
#include "paths.h"
#include "creating_arrays.h"

Simulation::Simulation() {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	newlyStarted = false;
	preserveChargeGlobal = true;
	arrayCreated = false;
	outputDir = std::string(outputDirectory);
	inputDir = std::string(inputDirectory);
	reducedOutputDir = std::string(reducedOutputDirectory);

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
	verbosity = 0;
	//types = new ParticleTypeContainer[typesNumber];
	//concentrations = new double[typesNumber];
	//particlesPerBin = new int[typesNumber];
	shockWaveX = -1.0;
	derExPoint = 0;
	derConcentrationPoint = 0;
	constMeanElevelPoint = 0;
}

void Simulation::setSpaceForProc() {
	int tempXnumber = ((xnumberGeneral) / cartDim[0]);
	int modXnumber = (xnumberGeneral) % cartDim[0];

	if (cartCoord[0] >= cartDim[0] - modXnumber) {
		xnumber = tempXnumber + 1;
		firstAbsoluteXindex = xnumberGeneral - (xnumber) * (cartDim[0] - cartCoord[0]) - additionalBinNumber - 1;
	} else {
		xnumber = tempXnumber;
		firstAbsoluteXindex = (xnumber) * cartCoord[0] - additionalBinNumber - 1;
	}

	int tempYnumber = ((ynumberGeneral) / cartDim[1]);
	int modYnumber = (ynumberGeneral) % cartDim[1];

	if (cartCoord[1] >= cartDim[1] - modYnumber) {
		ynumber = tempYnumber + 1;
		firstAbsoluteYindex = ynumberGeneral - (ynumber) * (cartDim[1] - cartCoord[1]) - additionalBinNumber - 1;
	} else {
		ynumber = tempYnumber;
		firstAbsoluteYindex = (ynumber) * cartCoord[1] - additionalBinNumber - 1;
	}

	int tempZnumber = ((znumberGeneral) / cartDim[2]);
	int modZnumber = (znumberGeneral) % cartDim[2];

	if (cartCoord[2] >= cartDim[2] - modZnumber) {
		znumber = tempZnumber + 1;
		firstAbsoluteZindex = znumberGeneral - (znumber) * (cartDim[2] - cartCoord[2]) - additionalBinNumber - 1;
	} else {
		znumber = tempZnumber;
		firstAbsoluteZindex = (znumber) * cartCoord[2] - additionalBinNumber - 1;
	}

	//todo boundary conditiontype
	/*if (cartDim[0] > 5) {
		int firstSmallRegions = 2 * cartDim[0] / 3;
		int firstSmallXnumber = xnumberGeneral / 3;
		int lastLargeRegions = cartDim[0] - firstSmallRegions;
		int lastLargeXnumber = xnumberGeneral - firstSmallXnumber;
		if (cartCoord[0] < firstSmallRegions) {
			tempXnumber = firstSmallXnumber / firstSmallRegions;
			modXnumber = firstSmallXnumber % firstSmallRegions;
			if (cartCoord[0] >= firstSmallRegions - modXnumber) {
				xnumber = tempXnumber + 1;
				firstAbsoluteXindex = firstSmallXnumber - (xnumber) * (firstSmallRegions - cartCoord[0]) - additionalBinNumber - 1;
			} else {
				xnumber = tempXnumber;
				firstAbsoluteXindex = (xnumber) * cartCoord[0] - additionalBinNumber - 1;
			}
		} else {
			tempXnumber = lastLargeXnumber / lastLargeRegions;
			modXnumber = lastLargeXnumber % lastLargeRegions;
			if (cartCoord[0] >= cartDim[0] - modXnumber) {
				xnumber = tempXnumber + 1;
				firstAbsoluteXindex = xnumberGeneral - (xnumber) * (cartDim[0] - cartCoord[0]) - additionalBinNumber - 1;
			} else {
				xnumber = tempXnumber;
				firstAbsoluteXindex = firstSmallXnumber + (xnumber) * (cartCoord[0] - firstSmallRegions) - additionalBinNumber - 1;
			}
		}

	}*/
	//printf("xnumber = %d\n", xnumber);
	xnumberAdded = xnumber + 2 + 2 * additionalBinNumber;
	ynumberAdded = ynumber + 2 + 2 * additionalBinNumber;
	znumberAdded = znumber + 2 + 2 * additionalBinNumber;

	xsize = xnumberAdded * xsizeGeneral / xnumberGeneral;
	ysize = ynumberAdded * ysizeGeneral / ynumberGeneral;
	zsize = znumberAdded * zsizeGeneral / znumberGeneral;

	//deltaX = xsize / (xnumber);
	//deltaY = ysize / (ynumber);
	//deltaZ = zsize / (znumber);

	deltaX = xsizeGeneral / (xnumberGeneral);
	deltaY = ysizeGeneral / (ynumberGeneral);
	deltaZ = zsizeGeneral / (znumberGeneral);

	deltaX2 = deltaX * deltaX;
	deltaY2 = deltaY * deltaY;
	deltaZ2 = deltaZ * deltaZ;

	cellVolume = deltaX * deltaY * deltaZ;

	leftX = xsizeGeneral + (firstAbsoluteXindex) * deltaX;
	rightX = leftX + xsize;
	leftY = ysizeGeneral + (firstAbsoluteYindex) * deltaY;
	rightY = leftY + ysize;
	leftZ = zsizeGeneral + (firstAbsoluteZindex) * deltaZ;
	rightZ = leftZ + zsize;
	derExPoint = 0;
	constMeanElevelPoint = 0;
}

Simulation::Simulation(int xn, int yn, int zn, double xsizev, double ysizev, double zsizev, double temp, double Vx,
                       double Vy, double Vz, double Ex, double Ey, double Ez, double Bx, double By, double Bz,
                       int maxIterations, double maxTimeV, int typesNumberV, int* particlesPerBinV,
                       double* concentrationsV, int inType, int nprocsV, int verbosityV, double preferedTimeStepV,
                       double massElectronInputV, MPI_Comm& comm) {
	nprocs = nprocsV;
	cartComm = comm;
	int periods[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	MPI_Comm_rank(cartComm, &rank);
	int tempCoord[3];
	for (int i = 0; i < 3; ++i) {
		tempCoord[i] = cartCoord[i];
	}
	tempCoord[0] -= 1;
	if (tempCoord[0] < 0) {
		tempCoord[0] = cartDim[0] - 1;
	}
	MPI_Cart_rank(cartComm, tempCoord, &leftRank);
	tempCoord[0] = cartCoord[0] + 1;
	if (tempCoord[0] >= cartDim[0]) {
		tempCoord[0] = 0;
	}
	MPI_Cart_rank(cartComm, tempCoord, &rightRank);
	tempCoord[0] = cartCoord[0];
	tempCoord[1] -= 1;
	if (tempCoord[1] < 0) {
		tempCoord[1] = cartDim[1] - 1;
	}
	MPI_Cart_rank(cartComm, tempCoord, &frontRank);
	tempCoord[1] = cartCoord[1] + 1;
	if (tempCoord[1] >= cartDim[1]) {
		tempCoord[1] = 0;
	}
	MPI_Cart_rank(cartComm, tempCoord, &backRank);
	tempCoord[1] = cartCoord[1];
	tempCoord[2] -= 1;
	if (tempCoord[2] < 0) {
		tempCoord[2] = cartDim[2] - 1;
	}
	MPI_Cart_rank(cartComm, tempCoord, &bottomRank);
	tempCoord[2] = cartCoord[2] + 1;
	if (tempCoord[2] >= cartDim[2]) {
		tempCoord[2] = 0;
	}
	MPI_Cart_rank(cartComm, tempCoord, &topRank);

	arrayCreated = false;
	timing = true;
	outputDir = std::string(outputDirectory);
	inputDir = std::string(inputDirectory);
	reducedOutputDir = std::string(reducedOutputDirectory);
	if (inType == 0) {
		inputType = CGS;
	} else if (inType == 1) {
		inputType = Theoretical;
	} else {
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
	//solverType = IMPLICIT_EC; //не явный с сохранением энергии
	//solverType = EXPLICIT; //явный
	//solverType = BUNEMAN;
	boundaryConditionTypeX = PERIODIC;
	//boundaryConditionTypeX = SUPER_CONDUCTOR_LEFT;
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
	//eta = 1.0;

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

	electronMassInput = massElectronInputV;
	preferedDeltaT = preferedTimeStepV;


	if (inputType == CGS) {
		massProton = massProtonReal;
		//massElectron = massElectronFactor * massElectronReal;
		massElectron = electronMassInput;
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

		double gamma = 1.0 / sqrt(1.0 - V0.scalarMult(V0) / (speed_of_light * speed_of_light));

		//plasma_period = sqrt(1 / omega2) * (2 * pi) * gamma * sqrt(gamma);
		plasma_period = sqrt(1 / omega2) * sqrt(gamma);
		//plasma_period = sqrt(1 / omega2) * (2 * pi) * gamma * sqrt(gamma)/(2*pi);
		double thermal_momentum;
		if (kBoltzman * temperature > massElectron * speed_of_light * speed_of_light) {
			thermal_momentum = kBoltzman * temperature / speed_of_light;
		} else {
			thermal_momentum = sqrt(2 * massElectron * kBoltzman * temperature);
		}
		thermal_momentum += V0.norm() * massElectron;
		double gyro_radius = thermal_momentum * speed_of_light / (electron_charge * B0.norm());
		//if (B0.norm() <= 0) {
		//scaleFactor = 1.0;
		//}
		scaleFactor = speed_of_light * plasma_period;

		//plasma_period = 1.0;
		//scaleFactor = 1.0;

		//scaleFactor = xsize;
		//scaleFactor = xsizeGeneral;
		//plasma_period = xsizeGeneral/(xnumberGeneral*speed_of_light);

		E0 = E0 * (plasma_period * sqrt(scaleFactor));
		B0 = B0 * (plasma_period * sqrt(scaleFactor));
		V0 = V0 * plasma_period / scaleFactor;

		rescaleConstants();

		density = density * cube(scaleFactor);
		for (int i = 0; i < typesNumber; ++i) {
			concentrations[i] = concentrations[i] * cube(scaleFactor);
		}

		//printf("scaleFactor = %lf\n", scaleFactor);

		deltaX /= scaleFactor;
		deltaY /= scaleFactor;
		deltaZ /= scaleFactor;
		deltaX2 /= scaleFactor * scaleFactor;
		deltaY2 /= scaleFactor * scaleFactor;
		deltaZ2 /= scaleFactor * scaleFactor;
        cellVolume /= scaleFactor * scaleFactor * scaleFactor;

		leftX /= scaleFactor;
		rightX /= scaleFactor;
		xsize /= scaleFactor;
		xsizeGeneral /= scaleFactor;
		leftY /= scaleFactor;
		rightY /= scaleFactor;
		ysize /= scaleFactor;
		ysizeGeneral /= scaleFactor;
		leftZ /= scaleFactor;
		rightZ /= scaleFactor;
		zsize /= scaleFactor;
		zsizeGeneral /= scaleFactor;
		if (rank == 0) printf("xsize/scaleFactor = %lf\n", xsize);
		//fflush(stdout);
	} else {
		massProton = 1.0;
		massElectron = massProton / 256;
		massAlpha = massProton * massAlphaReal / massProtonReal;
		massDeuterium = massProton * massDeuteriumReal / massProtonReal;
		massHelium3 = massProton * massHelium3Real / massProtonReal;
		double protonScale = massProton / massProtonReal;
		plasma_period = cube(speed_of_light) / (protonScale * electron_charge);
		scaleFactor = speed_of_light * plasma_period;
		rescaleConstantsToTheoretical();
		double densityForUnits = massElectron * (massProton + massElectron) / (4 * pi * electron_charge_normalized *
			electron_charge_normalized);
		if (fabs(density - densityForUnits) / (densityForUnits + densityForUnits) > 1E-3) {
			if (rank == 0) printf("density must be changed\n");
			if (rank == 0) printf("density = %g, must be = %g\n", density, densityForUnits);
			fflush(stdout);
		}
	}

	resistiveLayerWidth = defaulResistiveLayerWidth;
	fakeCondactivity = 2 * speed_of_light_normalized / (deltaX * (resistiveLayerWidth + 1));

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
	derExPoint = 0;
	constMeanElevelPoint = 0;
	derConcentrationPoint = 0;
	//if(rank == 0) printf("end constructor\n");
	//fflush(stdout);
}

Simulation::Simulation(int xn, int yn, int zn, double xsizev, double ysizev, double zsizev, double temp, double Vx,
                       double Vy, double Vz, double Ex, double Ey, double Ez, double Bx, double By, double Bz,
                       int maxIterations, double maxTimeV, int typesNumberV, int* particlesPerBinV,
                       double* concentrationsV, int inType, int nprocsV, int verbosityV, double preferedTimeStepV,
                       double massElectronInputV, double plasmaPeriodV,
                       double scaleFactorV, SolverType solverTypev, MPI_Comm& comm) {
	nprocs = nprocsV;
	cartComm = comm;
	int periods[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	MPI_Comm_rank(cartComm, &rank);
	int tempCoord[3];
	for (int i = 0; i < 3; ++i) {
		tempCoord[i] = cartCoord[i];
	}
	tempCoord[0] -= 1;
	if (tempCoord[0] < 0) {
		tempCoord[0] = cartDim[0] - 1;
	}
	MPI_Cart_rank(cartComm, tempCoord, &leftRank);
	tempCoord[0] = cartCoord[0] + 1;
	if (tempCoord[0] >= cartDim[0]) {
		tempCoord[0] = 0;
	}
	MPI_Cart_rank(cartComm, tempCoord, &rightRank);
	tempCoord[0] = cartCoord[0];
	tempCoord[1] -= 1;
	if (tempCoord[1] < 0) {
		tempCoord[1] = cartDim[1] - 1;
	}
	MPI_Cart_rank(cartComm, tempCoord, &frontRank);
	tempCoord[1] = cartCoord[1] + 1;
	if (tempCoord[1] >= cartDim[1]) {
		tempCoord[1] = 0;
	}
	MPI_Cart_rank(cartComm, tempCoord, &backRank);
	tempCoord[1] = cartCoord[1];
	tempCoord[2] -= 1;
	if (tempCoord[2] < 0) {
		tempCoord[2] = cartDim[2] - 1;
	}
	MPI_Cart_rank(cartComm, tempCoord, &bottomRank);
	tempCoord[2] = cartCoord[2] + 1;
	if (tempCoord[2] >= cartDim[2]) {
		tempCoord[2] = 0;
	}
	MPI_Cart_rank(cartComm, tempCoord, &topRank);

	arrayCreated = false;
	timing = true;
	outputDir = std::string(outputDirectory);
	inputDir = std::string(inputDirectory);
	reducedOutputDir = std::string(reducedOutputDirectory);
	if (inType == 0) {
		inputType = CGS;
	} else if (inType == 1) {
		inputType = Theoretical;
	} else {
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
	solverType = solverTypev;
	boundaryConditionTypeX = PERIODIC;
	//boundaryConditionTypeX = SUPER_CONDUCTOR_LEFT;
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
	electronMassInput = massElectronInputV;
	preferedDeltaT = preferedTimeStepV;

	if (inputType == CGS) {
		massProton = massProtonReal;
		//massElectron = massElectronFactor * massElectronReal;
		massElectron = electronMassInput;
		massAlpha = massAlphaReal;
		massDeuterium = massDeuteriumReal;
		massHelium3 = massHelium3Real;


		plasma_period = plasmaPeriodV;

		scaleFactor = scaleFactorV;


		//todo?

		rescaleConstants();

		/*E0 = E0 * (plasma_period * sqrt(scaleFactor));
		B0 = B0 * (plasma_period * sqrt(scaleFactor));
		V0 = V0 * plasma_period / scaleFactor;

		density = density * cube(scaleFactor);
		for (int i = 0; i < typesNumber; ++i) {
			concentrations[i] = concentrations[i] * cube(scaleFactor);
		}

		//printf("scaleFactor = %lf\n", scaleFactor);

		deltaX /= scaleFactor;
		deltaY /= scaleFactor;
		deltaZ /= scaleFactor;
		deltaX2 /= scaleFactor * scaleFactor;
		deltaY2 /= scaleFactor * scaleFactor;
		deltaZ2 /= scaleFactor * scaleFactor;
        cellVolume /= scaleFactor * scaleFactor * scaleFactor;

		leftX /= scaleFactor;
		rightX /= scaleFactor;
		xsize /= scaleFactor;
		xsizeGeneral /= scaleFactor;
		leftY /= scaleFactor;
		rightY /= scaleFactor;
		ysize /= scaleFactor;
		ysizeGeneral /= scaleFactor;
		leftZ /= scaleFactor;
		rightZ /= scaleFactor;
		zsize /= scaleFactor;
		zsizeGeneral /= scaleFactor;*/


		//printf("scaleFactor = %lf\n", scaleFactor);

		if (rank == 0) printf("xsize/scaleFactor = %lf\n", xsize);
		//fflush(stdout);
	} else {
		massProton = 1.0;
		massElectron = massProton / 256;
		massAlpha = massProton * massAlphaReal / massProtonReal;
		massDeuterium = massProton * massDeuteriumReal / massProtonReal;
		massHelium3 = massProton * massHelium3Real / massProtonReal;
		double protonScale = massProton / massProtonReal;
		plasma_period = plasmaPeriodV;
		scaleFactor = scaleFactorV;
		rescaleConstantsToTheoretical();
		double densityForUnits = massElectron * (massProton + massElectron) / (4 * pi * electron_charge_normalized *
			electron_charge_normalized);
		if (fabs(density - densityForUnits) / (densityForUnits + densityForUnits) > 1E-3) {
			if (rank == 0) printf("density must be changed\n");
			if (rank == 0) printf("density = %g, must be = %g\n", density, densityForUnits);
			fflush(stdout);
		}
	}

	resistiveLayerWidth = defaulResistiveLayerWidth;
	fakeCondactivity = 2 * speed_of_light_normalized / (deltaX * (resistiveLayerWidth + 1));

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
	shockWaveX = -1.0;
	derExPoint = 0;
	constMeanElevelPoint = 0;
	derConcentrationPoint = 0;
}

Simulation::~Simulation() {
	if (arrayCreated) {
		//if(false){
		delete[] mostAcceleratedParticlesNumbers;
		for (int i = 0; i < trackedParticlesNumber; ++i) {
			delete[] trackedParticlesNumbers[i];
		}
		delete[] trackedParticlesNumbers;
		delete gmresMaxwellBasis;
		delete gmresCleanupBasis;
		delete[] types;
		delete[] concentrations;
		delete[] particlesPerBin;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					delete[] gmresOutput[i][j][k];
				}
				delete[] gmresOutput[i][j];
			}
			delete[] gmresOutput[i];
		}
		delete[] gmresOutput;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
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

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
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

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					delete[] divergenceCleaningPotential[i][j][k];
					delete[] tempDivergenceCleaningPotential[i][j][k];
				}
				delete[] divergenceCleaningPotential[i][j];
				delete[] tempDivergenceCleaningPotential[i][j];
			}
			delete[] divergenceCleaningPotential[i];
			delete[] tempDivergenceCleaningPotential[i];
		}
		delete[] divergenceCleaningPotential;
		delete[] tempDivergenceCleaningPotential;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				delete[] divergenceCleaningPotentialFourier[i][j];
			}
			delete[] divergenceCleaningPotentialFourier[i];
		}
		delete[] divergenceCleaningPotentialFourier;

		delete3complexArray(fourierInput, xnumberAdded, ynumberAdded, znumberAdded);
		delete3complexArray(fourierImage, xnumberAdded, ynumberAdded, znumberAdded);
		delete3complexArray(fourierOutput, xnumberAdded, ynumberAdded, znumberAdded);

		delete3complexArray(fourierScalarInput, xnumberAdded, ynumberAdded, znumberAdded);
		delete3complexArray(fourierScalarOutput, xnumberAdded, ynumberAdded, znumberAdded);
		delete3complexArray(fourierScalarTempOutput, xnumberAdded, ynumberAdded, znumberAdded);
		delete3complexArray(fourierScalarTempOutput1, xnumberAdded, ynumberAdded, znumberAdded);

		delete3complexArray(fourierScalarMirrorInput, 2*xnumberAdded, ynumberAdded, znumberAdded);
		delete3complexArray(fourierScalarMirrorOutput, 2*xnumberAdded, ynumberAdded, znumberAdded);
		delete3complexArray(fourierScalarMirrorTempOutput, 2*xnumberAdded, ynumberAdded, znumberAdded);
		delete3complexArray(fourierScalarMirrorTempOutput1, 2*xnumberAdded, ynumberAdded, znumberAdded);

		delete[] localFactorX;
		delete[] localFactorY;
		delete[] localFactorZ;

		delete[] rightMeanElevel;

		delete3vectorArray(Efield, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(newEfield, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(tempEfield, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(explicitEfield, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(rotB, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(Ederivative, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3array(tempNodeParameter, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(tempNodeVectorParameter, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3matrixArray(tempNodeMatrixParameter, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				delete[] massMatrix[i][j];
				delete[] tempMassMatrix[i][j];
			}
			delete[] massMatrix[i];
			delete[] tempMassMatrix[i];
		}

		delete[] massMatrix;
		delete[] tempMassMatrix;

		delete3vectorArray(electricFlux, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(electricFluxMinus, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(externalElectricFlux, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(divPressureTensor, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(divPressureTensorMinus, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		delete3matrixArray(dielectricTensor, xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);

		for (int j = 0; j < ynumberAdded + 1; ++j) {
			delete[] leftElevel[j];
			delete[] rightElevel[j];
		}
		delete[] leftElevel;
		delete[] rightElevel;

		for (int t = 0; t < typesNumber; ++t) {
			delete3array(particleConcentrations[t], xnumberAdded, ynumberAdded, znumberAdded);
			delete3array(particleEnergies[t], xnumberAdded, ynumberAdded, znumberAdded);
			delete3vectorArray(particleBulkVelocities[t], xnumberAdded, ynumberAdded, znumberAdded);
		}
		delete[] particleConcentrations;
		delete[] particleEnergies;
		delete[] particleBulkVelocities;

		delete3vectorArray(Bfield, xnumberAdded, ynumberAdded, znumberAdded);
		delete3vectorArray(newBfield, xnumberAdded, ynumberAdded, znumberAdded);
		delete3array(chargeDensity, xnumberAdded, ynumberAdded, znumberAdded);
		delete3array(chargeDensityMinus, xnumberAdded, ynumberAdded, znumberAdded);
		delete3matrixArray(pressureTensor, xnumberAdded, ynumberAdded, znumberAdded);
		delete3array(chargeDensityHat, xnumberAdded, ynumberAdded, znumberAdded);
		delete3array(tempCellParameter, xnumberAdded, ynumberAdded, znumberAdded);
		delete3vectorArray(tempCellVectorParameter, xnumberAdded, ynumberAdded, znumberAdded);
		delete3matrixArray(tempCellMatrixParameter, xnumberAdded, ynumberAdded, znumberAdded);
		delete3vectorArray(rotE, xnumberAdded, ynumberAdded, znumberAdded);
		delete3vectorArray(Bderivative, xnumberAdded, ynumberAdded, znumberAdded);

		delete[] xgrid;
		delete[] middleXgrid;
		delete[] ygrid;
		delete[] middleYgrid;
		delete[] zgrid;
		delete[] middleZgrid;

		delete[] rightOutNodeBuffer;
		delete[] rightInNodeBuffer;
		delete[] leftOutNodeBuffer;
		delete[] leftInNodeBuffer;

		delete[] rightOutVectorNodeBuffer;
		delete[] rightInVectorNodeBuffer;
		delete[] leftOutVectorNodeBuffer;
		delete[] leftInVectorNodeBuffer;

		delete[] rightOutMatrixNodeBuffer;
		delete[] rightInMatrixNodeBuffer;
		delete[] leftOutMatrixNodeBuffer;
		delete[] leftInMatrixNodeBuffer;

		delete[] rightOutVectorCellBuffer;
		delete[] rightInVectorCellBuffer;
		delete[] leftOutVectorCellBuffer;
		delete[] leftInVectorCellBuffer;

		delete[] rightOutCellBuffer;
		delete[] rightInCellBuffer;
		delete[] leftOutCellBuffer;
		delete[] leftInCellBuffer;

		delete[] backOutNodeBuffer;
		delete[] backInNodeBuffer;
		delete[] frontOutNodeBuffer;
		delete[] frontInNodeBuffer;

		delete[] backOutVectorNodeBuffer;
		delete[] backInVectorNodeBuffer;
		delete[] frontOutVectorNodeBuffer;
		delete[] frontInVectorNodeBuffer;

		delete[] backOutMatrixNodeBuffer;
		delete[] backInMatrixNodeBuffer;
		delete[] frontOutMatrixNodeBuffer;
		delete[] frontInMatrixNodeBuffer;

		delete[] backOutVectorCellBuffer;
		delete[] backInVectorCellBuffer;
		delete[] frontOutVectorCellBuffer;
		delete[] frontInVectorCellBuffer;

		delete[] backOutCellBuffer;
		delete[] backInCellBuffer;
		delete[] frontOutCellBuffer;
		delete[] frontInCellBuffer;

		delete[] topOutNodeBuffer;
		delete[] topInNodeBuffer;
		delete[] bottomOutNodeBuffer;
		delete[] bottomInNodeBuffer;

		delete[] topOutVectorNodeBuffer;
		delete[] topInVectorNodeBuffer;
		delete[] bottomOutVectorNodeBuffer;
		delete[] bottomInVectorNodeBuffer;

		delete[] topOutMatrixNodeBuffer;
		delete[] topInMatrixNodeBuffer;
		delete[] bottomOutMatrixNodeBuffer;
		delete[] bottomInMatrixNodeBuffer;

		delete[] topOutVectorCellBuffer;
		delete[] topInVectorCellBuffer;
		delete[] bottomOutVectorCellBuffer;
		delete[] bottomInVectorCellBuffer;

		delete[] topOutCellBuffer;
		delete[] topInCellBuffer;
		delete[] bottomOutCellBuffer;
		delete[] bottomInCellBuffer;

		/// Masha's buffers

		delete[] rightOutNodeBufferMasha;
		delete[] rightInNodeBufferMasha;
		delete[] leftOutNodeBufferMasha;
		delete[] leftInNodeBufferMasha;

		delete[] rightOutVectorNodeBufferMasha;
		delete[] rightInVectorNodeBufferMasha;
		delete[] leftOutVectorNodeBufferMasha;
		delete[] leftInVectorNodeBufferMasha;

		delete[] rightOutMatrixNodeBufferMasha;
		delete[] rightInMatrixNodeBufferMasha;
		delete[] leftOutMatrixNodeBufferMasha;
		delete[] leftInMatrixNodeBufferMasha;

		delete[] rightOutVectorCellBufferMasha;
		delete[] rightInVectorCellBufferMasha;
		delete[] leftOutVectorCellBufferMasha;
		delete[] leftInVectorCellBufferMasha;

		delete[] rightOutCellBufferMasha;
		delete[] rightInCellBufferMasha;
		delete[] leftOutCellBufferMasha;
		delete[] leftInCellBufferMasha;

		delete[] backOutNodeBufferMasha;
		delete[] backInNodeBufferMasha;
		delete[] frontOutNodeBufferMasha;
		delete[] frontInNodeBufferMasha;

		delete[] backOutVectorNodeBufferMasha;
		delete[] backInVectorNodeBufferMasha;
		delete[] frontOutVectorNodeBufferMasha;
		delete[] frontInVectorNodeBufferMasha;

		delete[] backOutMatrixNodeBufferMasha;
		delete[] backInMatrixNodeBufferMasha;
		delete[] frontOutMatrixNodeBufferMasha;
		delete[] frontInMatrixNodeBufferMasha;

		delete[] backOutVectorCellBufferMasha;
		delete[] backInVectorCellBufferMasha;
		delete[] frontOutVectorCellBufferMasha;
		delete[] frontInVectorCellBufferMasha;

		delete[] backOutCellBufferMasha;
		delete[] backInCellBufferMasha;
		delete[] frontOutCellBufferMasha;
		delete[] frontInCellBufferMasha;

		delete[] topOutNodeBufferMasha;
		delete[] topInNodeBufferMasha;
		delete[] bottomOutNodeBufferMasha;
		delete[] bottomInNodeBufferMasha;

		delete[] topOutVectorNodeBufferMasha;
		delete[] topInVectorNodeBufferMasha;
		delete[] bottomOutVectorNodeBufferMasha;
		delete[] bottomInVectorNodeBufferMasha;

		delete[] topOutMatrixNodeBufferMasha;
		delete[] topInMatrixNodeBufferMasha;
		delete[] bottomOutMatrixNodeBufferMasha;
		delete[] bottomInMatrixNodeBufferMasha;

		delete[] topOutVectorCellBufferMasha;
		delete[] topInVectorCellBufferMasha;
		delete[] bottomOutVectorCellBufferMasha;
		delete[] bottomInVectorCellBufferMasha;

		delete[] topOutCellBufferMasha;
		delete[] topInCellBufferMasha;
		delete[] bottomOutCellBufferMasha;
		delete[] bottomInCellBufferMasha;

		////buneman E
		delete[] leftOutBunemanExBuffer;
		delete[] rightOutBunemanExBuffer;
		delete[] leftInBunemanExBuffer;
		delete[] rightInBunemanExBuffer;

		delete[] frontOutBunemanExBuffer;
		delete[] backOutBunemanExBuffer;
		delete[] frontInBunemanExBuffer;
		delete[] backInBunemanExBuffer;

		delete[] bottomOutBunemanExBuffer;
		delete[] topOutBunemanExBuffer;
		delete[] bottomInBunemanExBuffer;
		delete[] topInBunemanExBuffer;

		delete[] leftOutBunemanEyBuffer;
		delete[] rightOutBunemanEyBuffer;
		delete[] leftInBunemanEyBuffer;
		delete[] rightInBunemanEyBuffer;

		delete[] frontOutBunemanEyBuffer;
		delete[] backOutBunemanEyBuffer;
		delete[] frontInBunemanEyBuffer;
		delete[] backInBunemanEyBuffer;

		delete[] bottomOutBunemanEyBuffer;
		delete[] topOutBunemanEyBuffer;
		delete[] bottomInBunemanEyBuffer;
		delete[] topInBunemanEyBuffer;

		delete[] leftOutBunemanEzBuffer;
		delete[] rightOutBunemanEzBuffer;
		delete[] leftInBunemanEzBuffer;
		delete[] rightInBunemanEzBuffer;

		delete[] frontOutBunemanEzBuffer;
		delete[] backOutBunemanEzBuffer;
		delete[] frontInBunemanEzBuffer;
		delete[] backInBunemanEzBuffer;

		delete[] bottomOutBunemanEzBuffer;
		delete[] topOutBunemanEzBuffer;
		delete[] bottomInBunemanEzBuffer;
		delete[] topInBunemanEzBuffer;

		///buneman B
		delete[] leftOutBunemanBxBuffer;
		delete[] rightOutBunemanBxBuffer;
		delete[] leftInBunemanBxBuffer;
		delete[] rightInBunemanBxBuffer;

		delete[] frontOutBunemanBxBuffer;
		delete[] backOutBunemanBxBuffer;
		delete[] frontInBunemanBxBuffer;
		delete[] backInBunemanBxBuffer;

		delete[] bottomOutBunemanBxBuffer;
		delete[] topOutBunemanBxBuffer;
		delete[] bottomInBunemanBxBuffer;
		delete[] topInBunemanBxBuffer;

		delete[] leftOutBunemanByBuffer;
		delete[] rightOutBunemanByBuffer;
		delete[] leftInBunemanByBuffer;
		delete[] rightInBunemanByBuffer;

		delete[] frontOutBunemanByBuffer;
		delete[] backOutBunemanByBuffer;
		delete[] frontInBunemanByBuffer;
		delete[] backInBunemanByBuffer;

		delete[] bottomOutBunemanByBuffer;
		delete[] topOutBunemanByBuffer;
		delete[] bottomInBunemanByBuffer;
		delete[] topInBunemanByBuffer;

		delete[] leftOutBunemanBzBuffer;
		delete[] rightOutBunemanBzBuffer;
		delete[] leftInBunemanBzBuffer;
		delete[] rightInBunemanBzBuffer;

		delete[] frontOutBunemanBzBuffer;
		delete[] backOutBunemanBzBuffer;
		delete[] frontInBunemanBzBuffer;
		delete[] backInBunemanBzBuffer;

		delete[] bottomOutBunemanBzBuffer;
		delete[] topOutBunemanBzBuffer;
		delete[] bottomInBunemanBzBuffer;
		delete[] topInBunemanBzBuffer;

		delete[] rightOutGmresBuffer;
		delete[] rightInGmresBuffer;
		delete[] leftOutGmresBuffer;
		delete[] leftInGmresBuffer;

		delete[] backOutGmresBuffer;
		delete[] backInGmresBuffer;
		delete[] frontOutGmresBuffer;
		delete[] frontInGmresBuffer;

		delete[] topOutGmresBuffer;
		delete[] topInGmresBuffer;
		delete[] bottomOutGmresBuffer;
		delete[] bottomInGmresBuffer;

		delete[] rightOutDivergenceBuffer;
		delete[] leftOutDivergenceBuffer;
		delete[] leftInDivergenceBuffer;
		delete[] rightInDivergenceBuffer;

		delete[] frontOutDivergenceBuffer;
		delete[] backOutDivergenceBuffer;
		delete[] frontInDivergenceBuffer;
		delete[] backInDivergenceBuffer;

		delete[] topOutDivergenceBuffer;
		delete[] topInDivergenceBuffer;
		delete[] bottomOutDivergenceBuffer;
		delete[] bottomInDivergenceBuffer;

		for (int i = 0; i < particles.size(); ++i) {
			Particle* particle = particles[i];
			delete particle;
		}
		particles.clear();
		tempParticles.clear();

		delete3array(tempCellParameterLeft, 2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);
		delete3array(tempCellParameterRight, 2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);
		delete3vectorArray(tempCellVectorParameterLeft, 2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);
		delete3vectorArray(tempCellVectorParameterRight, 2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);
		delete3matrixArray(tempCellMatrixParameterLeft, 2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);
		delete3matrixArray(tempCellMatrixParameterRight, 2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);

		delete3array(tempNodeParameterLeft, 3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);
		delete3array(tempNodeParameterRight, 3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(tempNodeVectorParameterLeft, 3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);
		delete3vectorArray(tempNodeVectorParameterRight, 3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);
		delete3matrixArray(tempNodeMatrixParameterLeft, 3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);
		delete3matrixArray(tempNodeMatrixParameterRight, 3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);

		delete3array(tempCellParameterFront, xnumberAdded, 2 + 2 * additionalBinNumber, znumberAdded);
		delete3array(tempCellParameterBack, xnumberAdded, 2 + 2 * additionalBinNumber, znumberAdded);
		delete3vectorArray(tempCellVectorParameterFront, xnumberAdded, 2 + 2 * additionalBinNumber, znumberAdded);
		delete3vectorArray(tempCellVectorParameterBack, xnumberAdded, 2 + 2 * additionalBinNumber, znumberAdded);
		delete3matrixArray(tempCellMatrixParameterFront, xnumberAdded, 2 + 2 * additionalBinNumber, znumberAdded);
		delete3matrixArray(tempCellMatrixParameterBack, xnumberAdded, 2 + 2 * additionalBinNumber, znumberAdded);

		delete3array(tempNodeParameterFront, xnumberAdded + 1, 3 + 2 * additionalBinNumber, znumberAdded + 1);
		delete3array(tempNodeParameterBack, xnumberAdded + 1, 3 + 2 * additionalBinNumber, znumberAdded + 1);
		delete3vectorArray(tempNodeVectorParameterFront, xnumberAdded + 1, 3 + 2 * additionalBinNumber, znumberAdded + 1);
		delete3vectorArray(tempNodeVectorParameterBack, xnumberAdded + 1, 3 + 2 * additionalBinNumber, znumberAdded + 1);
		delete3matrixArray(tempNodeMatrixParameterFront, xnumberAdded + 1, 3 + 2 * additionalBinNumber, znumberAdded + 1);
		delete3matrixArray(tempNodeMatrixParameterBack, xnumberAdded + 1, 3 + 2 * additionalBinNumber, znumberAdded + 1);

		delete3array(tempCellParameterBottom, xnumberAdded, ynumberAdded, 2 + 2 * additionalBinNumber);
		delete3array(tempCellParameterTop, xnumberAdded, ynumberAdded, 2 + 2 * additionalBinNumber);
		delete3vectorArray(tempCellVectorParameterBottom, xnumberAdded, ynumberAdded, 2 + 2 * additionalBinNumber);
		delete3vectorArray(tempCellVectorParameterTop, xnumberAdded, ynumberAdded, 2 + 2 * additionalBinNumber);
		delete3matrixArray(tempCellMatrixParameterBottom, xnumberAdded, ynumberAdded, 2 + 2 * additionalBinNumber);
		delete3matrixArray(tempCellMatrixParameterTop, xnumberAdded, ynumberAdded, 2 + 2 * additionalBinNumber);

		delete3array(tempNodeParameterBottom, xnumberAdded + 1, ynumberAdded + 1, 3 + 2 * additionalBinNumber);
		delete3array(tempNodeParameterTop, xnumberAdded + 1, ynumberAdded + 1, 3 + 2 * additionalBinNumber);
		delete3vectorArray(tempNodeVectorParameterBottom, xnumberAdded + 1, ynumberAdded + 1, 3 + 2 * additionalBinNumber);
		delete3vectorArray(tempNodeVectorParameterTop, xnumberAdded + 1, ynumberAdded + 1, 3 + 2 * additionalBinNumber);
		delete3matrixArray(tempNodeMatrixParameterBottom, xnumberAdded + 1, ynumberAdded + 1, 3 + 2 * additionalBinNumber);
		delete3matrixArray(tempNodeMatrixParameterTop, xnumberAdded + 1, ynumberAdded + 1, 3 + 2 * additionalBinNumber);

		delete4array(residualBiconjugateDivE, xnumberAdded, ynumberAdded, znumberAdded, 1);
		delete4array(firstResidualBiconjugateDivE, xnumberAdded, ynumberAdded, znumberAdded, 1);
		delete4array(vBiconjugateDivE, xnumberAdded, ynumberAdded, znumberAdded, 1);
		delete4array(pBiconjugateDivE, xnumberAdded, ynumberAdded, znumberAdded, 1);
		delete4array(sBiconjugateDivE, xnumberAdded, ynumberAdded, znumberAdded, 1);
		delete4array(tBiconjugateDivE, xnumberAdded, ynumberAdded, znumberAdded, 1);

		delete4array(residualBiconjugateMaxwell, xnumberAdded, ynumberAdded, znumberAdded, 3);
		delete4array(firstResidualBiconjugateMaxwell, xnumberAdded, ynumberAdded, znumberAdded, 3);
		delete4array(vBiconjugateMaxwell, xnumberAdded, ynumberAdded, znumberAdded, 3);
		delete4array(pBiconjugateMaxwell, xnumberAdded, ynumberAdded, znumberAdded, 3);
		delete4array(sBiconjugateMaxwell, xnumberAdded, ynumberAdded, znumberAdded, 3);
		delete4array(tBiconjugateMaxwell, xnumberAdded, ynumberAdded, znumberAdded, 3);

		if (solverType == BUNEMAN) {

			delete3array(bunemanJx, xnumberAdded, ynumberAdded + 1, znumberAdded + 1);
			delete3array(bunemanEx, xnumberAdded, ynumberAdded + 1, znumberAdded + 1);
			delete3array(bunemanNewEx, xnumberAdded, ynumberAdded + 1, znumberAdded + 1);
			delete3array(tempBunemanExParameter, xnumberAdded, ynumberAdded + 1, znumberAdded + 1);
			delete3array(bunemanDivCleaningEx, xnumberAdded, ynumberAdded + 1, znumberAdded + 1);

			delete3array(bunemanJy, xnumberAdded + 1, ynumberAdded, znumberAdded + 1);
			delete3array(bunemanEy, xnumberAdded + 1, ynumberAdded, znumberAdded + 1);
			delete3array(bunemanNewEy, xnumberAdded + 1, ynumberAdded, znumberAdded + 1);
			delete3array(tempBunemanEyParameter, xnumberAdded + 1, ynumberAdded, znumberAdded + 1);
			delete3array(bunemanDivCleaningEy, xnumberAdded + 1, ynumberAdded, znumberAdded + 1);

			delete3array(bunemanJz, xnumberAdded + 1, ynumberAdded + 1, znumberAdded);
			delete3array(bunemanEz, xnumberAdded + 1, ynumberAdded + 1, znumberAdded);
			delete3array(bunemanNewEz, xnumberAdded + 1, ynumberAdded + 1, znumberAdded);
			delete3array(tempBunemanEzParameter, xnumberAdded + 1, ynumberAdded + 1, znumberAdded);
			delete3array(bunemanDivCleaningEz, xnumberAdded + 1, ynumberAdded + 1, znumberAdded);

			delete3array(bunemanBx, xnumberAdded + 1, ynumberAdded, znumberAdded);
			delete3array(bunemanNewBx, xnumberAdded + 1, ynumberAdded, znumberAdded);
			delete3array(tempBunemanBxParameter, xnumberAdded + 1, ynumberAdded, znumberAdded);
			delete3array(bunemanDivCleaningBx, xnumberAdded + 1, ynumberAdded, znumberAdded);

			delete3array(bunemanBy, xnumberAdded, ynumberAdded + 1, znumberAdded);
			delete3array(bunemanNewBy, xnumberAdded, ynumberAdded + 1, znumberAdded);
			delete3array(tempBunemanByParameter, xnumberAdded, ynumberAdded + 1, znumberAdded);
			delete3array(bunemanDivCleaningBy, xnumberAdded, ynumberAdded + 1, znumberAdded);

			delete3array(bunemanBz, xnumberAdded, ynumberAdded, znumberAdded + 1);
			delete3array(bunemanNewBz, xnumberAdded, ynumberAdded, znumberAdded + 1);
			delete3array(tempBunemanBzParameter, xnumberAdded, ynumberAdded, znumberAdded + 1);
			delete3array(bunemanDivCleaningBz, xnumberAdded, ynumberAdded, znumberAdded + 1);

			//// temp buneman j
			// left right
			delete3array(tempBunemanJxLeft, 2 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);
			delete3array(tempBunemanJxRight, 2 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);

			delete3array(tempBunemanJyLeft, 3 + 2*additionalBinNumber, ynumberAdded, znumberAdded + 1);
			delete3array(tempBunemanJyRight, 3 + 2*additionalBinNumber, ynumberAdded, znumberAdded + 1);

			delete3array(tempBunemanJzLeft, 3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded);
			delete3array(tempBunemanJzRight, 3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded);

			///front back
			delete3array(tempBunemanJxFront, xnumberAdded, 3 + 2*additionalBinNumber, znumberAdded + 1);
			delete3array(tempBunemanJxBack, xnumberAdded, 3 + 2*additionalBinNumber, znumberAdded + 1);
			
			delete3array(tempBunemanJyFront, xnumberAdded + 1, 2 + 2*additionalBinNumber, znumberAdded + 1);
			delete3array(tempBunemanJyBack, xnumberAdded + 1, 2 + 2*additionalBinNumber, znumberAdded + 1);
	
			delete3array(tempBunemanJzFront, xnumberAdded + 1, 3 + 2*additionalBinNumber, znumberAdded);
			delete3array(tempBunemanJzBack, xnumberAdded + 1, 3 + 2*additionalBinNumber, znumberAdded);
			
			// bottom top
			delete3array(tempBunemanJxBottom, xnumberAdded, ynumberAdded + 1, 3 + 2*additionalBinNumber);
			delete3array(tempBunemanJxTop, xnumberAdded, ynumberAdded + 1, 3 + 2*additionalBinNumber);

			delete3array(tempBunemanJyBottom, xnumberAdded + 1, ynumberAdded, 3 + 2*additionalBinNumber);
			delete3array(tempBunemanJyTop, xnumberAdded + 1, ynumberAdded, 3 + 2*additionalBinNumber);

			delete3array(tempBunemanJzBottom, xnumberAdded + 1, ynumberAdded + 1, 2 + 2*additionalBinNumber);
			delete3array(tempBunemanJzTop, xnumberAdded + 1, ynumberAdded + 1, 2 + 2*additionalBinNumber);
		}
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
	kBoltzman_normalized = kBoltzman * (plasma_period * plasma_period / (scaleFactor * scaleFactor)) * massProton /
		massProtonReal;
	speed_of_light_normalized = 1.0;
	speed_of_light_normalized_sqr = 1.0;
	electron_charge_normalized = 1.0;
}

void Simulation::initialize() {
	if (rank == 0) printf("initialization\n");
	fflush(stdout);

	boundaryConditionTypeX = PERIODIC;
	boundaryConditionTypeY = PERIODIC;
	boundaryConditionTypeZ = PERIODIC;
	//if (rank == 0) printLog("initialization\n");


	if (verbosity > 2) printf("creating sub cart rank = %d\n", rank);
	int dimsYZ[3];
	dimsYZ[0] = 0;
	dimsYZ[1] = 1;
	dimsYZ[2] = 1;
	MPI_Cart_sub(cartComm, dimsYZ, &cartCommYZ);
	int dimsZ[3];
	dimsZ[0] = 0;
	dimsZ[1] = 0;
	dimsZ[2] = 1;
	MPI_Cart_sub(cartComm, dimsZ, &cartCommZ);
	int dimsXY[3];
	dimsXY[0] = 1;
	dimsXY[1] = 1;
	dimsXY[2] = 0;
	MPI_Cart_sub(cartComm, dimsXY, &cartCommXY);
	int dimsY[3];
	dimsY[0] = 0;
	dimsY[1] = 1;
	dimsY[2] = 0;
	MPI_Cart_sub(cartComm, dimsY, &cartCommY);
	int dimsXZ[3];
	dimsXZ[0] = 1;
	dimsXZ[1] = 0;
	dimsXZ[2] = 1;
	MPI_Cart_sub(cartComm, dimsXZ, &cartCommXZ);
	int dimsX[3];
	dimsX[0] = 1;
	dimsX[1] = 0;
	dimsX[2] = 0;
	MPI_Cart_sub(cartComm, dimsX, &cartCommX);
	if (verbosity > 2) printf("finish creating sub cart rank = %d\n", rank);

	bool readTrackedParticles = false;

	if (readTrackedParticles) {
		if (rank == 0) {
			trackedParticlesNumbers = readTrackedParticlesNumbers((inputDir + "acceleratedParticlesNumbers.dat").c_str(),
			                                                      trackedParticlesNumber);
		}

		if (rank == 0) {
			for (int i = 1; i < nprocs; ++i) {
				int trackedN[1];
				trackedN[0] = trackedParticlesNumber;
				MPI_Send(trackedN, 1, MPI_INT, i, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm);
			}
		} else {
			int trackedN[1];
			MPI_Status status;
			MPI_Recv(trackedN, 1, MPI_INT, 0, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm, &status);
			trackedParticlesNumber = trackedN[0];
		}

		if (rank == 0) {
			int* tempNumbers = new int[trackedParticlesNumber];
			int* tempTypes = new int[trackedParticlesNumber];
			for (int i = 0; i < trackedParticlesNumber; ++i) {
				tempNumbers[i] = trackedParticlesNumbers[i][0];
				tempTypes[i] = trackedParticlesNumbers[i][1];
			}
			for (int i = 1; i < nprocs; ++i) {
				MPI_Send(tempNumbers, trackedParticlesNumber, MPI_INT, i, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm);
				MPI_Send(tempTypes, trackedParticlesNumber, MPI_INT, i, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm);
			}

			delete[] tempNumbers;
			delete[] tempTypes;
		} else {
			MPI_Status status;
			int* tempNumbers = new int[trackedParticlesNumber];
			int* tempTypes = new int[trackedParticlesNumber];
			trackedParticlesNumbers = new int*[trackedParticlesNumber];
			MPI_Recv(tempNumbers, trackedParticlesNumber, MPI_INT, 0, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm, &status);
			MPI_Recv(tempTypes, trackedParticlesNumber, MPI_INT, 0, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm, &status);
			for (int i = 0; i < trackedParticlesNumber; ++i) {
				trackedParticlesNumbers[i] = new int[2];
				trackedParticlesNumbers[i][0] = tempNumbers[i];
				trackedParticlesNumbers[i][1] = tempTypes[i];
			}
			delete[] tempNumbers;
			delete[] tempTypes;
		}
	} else {
		trackedParticlesNumber = 20;
		trackedParticlesNumbers = new int*[trackedParticlesNumber];
		for (int i = 0; i < trackedParticlesNumber; ++i) {
			trackedParticlesNumbers[i] = new int[2];
			trackedParticlesNumbers[i][0] = 2 * i;
			trackedParticlesNumbers[i][1] = 0;
		}
	}

	if (verbosity > 2) printf("initialize grid rank = %d\n", rank);

	xgrid[0] = leftX;

	for (int i = 1; i <= xnumberAdded; ++i) {
		xgrid[i] = xgrid[0] + i * deltaX;
	}

	//xgrid[xnumberAdded-1] = rightX;
	//xgrid[xnumber + 1] = rightX + deltaX;
	//xgrid[xnumberAdded] = xgrid[0] + (xnumberAdded) * deltaX;

	//printf("xgrid[0] = %lf xgrid[xnumber] = %lf\n", xgrid[0], xgrid[xnumber]);

	ygrid[0] = leftY;
	for (int j = 0; j <= ynumberAdded; ++j) {
		ygrid[j] = ygrid[0] + j * deltaY;
	}
	//ygrid[ynumber] = 2 * ysize;

	zgrid[0] = leftZ;
	for (int k = 0; k <= znumberAdded; ++k) {
		zgrid[k] = zgrid[0] + k * deltaZ;
	}
	//zgrid[znumber] = 2 * zsize;

	for (int i = 0; i < xnumberAdded; ++i) {
		middleXgrid[i] = (xgrid[i] + xgrid[i + 1]) / 2;
	}

	for (int j = 0; j < ynumberAdded; ++j) {
		middleYgrid[j] = (ygrid[j] + ygrid[j + 1]) / 2;
	}

	for (int k = 0; k < znumberAdded; ++k) {
		middleZgrid[k] = (zgrid[k] + zgrid[k + 1]) / 2;
	}

	if (verbosity > 2) printf("initialize fields rank = %d\n", rank);

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
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

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k] = B0;
				//Bfield[i][j][k].x = B0.x;
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}

	rightBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(E0, B0);
	leftBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(E0, B0);

	if (verbosity > 2) printf("creating particle types rank = %d\n", rank);

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

	if (rank == 0) printf("finish initialization\n");
	fflush(stdout);
	if (rank == 0) printLog("finish initialization\n");
}

void Simulation::initializeSimpleElectroMagneticWave() {
	boundaryConditionTypeX = PERIODIC;
	//E0 = Vector3d(0, 0, 0);
	//B0 = Vector3d(0, 0, 0);
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k] = Vector3d(0, 0, 0);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
	double kw = (2 * pi / xsizeGeneral);
	printf("kw = %15.10g\n", kw);
	fflush(stdout);
	//double E = 1E-5;
	double E = B0.norm();
	if (solverType == BUNEMAN) {
		double omega = speed_of_light_normalized * kw;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEx[i][j][k] = 0;
					bunemanNewEx[i][j][k] = 0;
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEy[i][j][k] = E * sin(kw * xgrid[i]);
					bunemanNewEy[i][j][k] = bunemanEy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanEz[i][j][k] = 0;
					bunemanNewEz[i][j][k] = bunemanEz[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBx[i][j][k] = 0;
					bunemanNewBx[i][j][k] = 0;
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBy[i][j][k] = 0;
					bunemanNewBy[i][j][k] = 0;
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanBz[i][j][k] = E * sin(kw * middleXgrid[i] - omega * deltaT / 2);
					//bunemanBz[i][j][k] = E * sin(kw * middleXgrid[i]);
					bunemanNewBz[i][j][k] = bunemanBz[i][j][k];
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					Efield[i][j][k].x = 0;
					Efield[i][j][k].y = E * sin(kw * xgrid[i]);
					Efield[i][j][k].z = 0;
					tempEfield[i][j][k] = Efield[i][j][k];
					newEfield[i][j][k] = tempEfield[i][j][k];
					explicitEfield[i][j][k] = Efield[i][j][k];
				}
			}
		}

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					Bfield[i][j][k].z = E * sin(kw * middleXgrid[i]);
					newBfield[i][j][k] = Bfield[i][j][k];
				}
			}
		}
	}
}

void Simulation::initializeSimpleElectroMagneticWaveY() {
	boundaryConditionTypeX = PERIODIC;
	//E0 = Vector3d(0, 0, 0);
	//B0 = Vector3d(0, 0, 0);
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k] = Vector3d(0, 0, 0);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
	double kw = (2 * pi / ysizeGeneral);
	printf("kw = %15.10g\n", kw);
	fflush(stdout);
	//double E = 1E-5;
	double E = B0.norm();

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				Efield[i][j][k].x = 0;
				Efield[i][j][k].z = E * sin(kw * ygrid[j]);
				Efield[i][j][k].y = 0;
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = tempEfield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k].x = E * sin(kw * middleYgrid[j]);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
}

void Simulation::initializeSimpleElectroMagneticWaveZ() {
	boundaryConditionTypeX = PERIODIC;
	//E0 = Vector3d(0, 0, 0);
	//B0 = Vector3d(0, 0, 0);
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k] = Vector3d(0, 0, 0);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
	double kw = (2 * pi / zsizeGeneral);
	printf("kw = %15.10g\n", kw);
	fflush(stdout);
	//double E = 1E-5;
	double E = B0.norm();

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				Efield[i][j][k].y = 0;
				Efield[i][j][k].x = E * sin(kw * zgrid[k]);
				Efield[i][j][k].z = 0;
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = tempEfield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k].y = E * sin(kw * middleZgrid[k]);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
}

void Simulation::initializeRotatedSimpleElectroMagneticWave(int waveCountX, int waveCountY, int waveCountZ) {
	boundaryConditionTypeX = PERIODIC;

	Eyamplitude = 1;
	Ezamplitude = 0;
	Bzamplitude = Eyamplitude;
	Byamplitude = Ezamplitude;

	double kx = waveCountX * 2 * pi / xsizeGeneral;
	double ky = waveCountY * 2 * pi / ysizeGeneral;
	double kz = waveCountZ * 2 * pi / zsizeGeneral;
	//kz = 0;

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


	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
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

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k].x = 0;
				Bfield[i][j][k].y = Byamplitude * sin(kx * middleXgrid[i] + ky * middleYgrid[j] + kz * middleZgrid[k]);
				Bfield[i][j][k].z = Bzamplitude * cos(kx * middleXgrid[i] + ky * middleYgrid[j] + kz * middleZgrid[k]);
				Bfield[i][j][k] = rotationMatrix * Bfield[i][j][k];
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
}

void Simulation::initializeAlfvenWaveX(int wavesCount, double amplitudeRelation) {
	boundaryConditionTypeX = PERIODIC;
	if (rank == 0) printf("initiali   zation alfven wave\n");
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
	VzamplitudeProton = -((1.0 / (4 * pi * concentration * electron_charge_normalized)) * (kc + ((omegae2 + omegap2 - omega
		* omega) / kc) + (omegae2 * Omegae2 / (kc * denominator))) / ((Omegae * omegae2 * omega / ((kc2 - omega * omega) *
		denominator)) + (Omegap / omega))) * Bzamplitude;
	//double
	VzamplitudeElectron = (((electron_charge_normalized * omega * Omegae) / (massElectron * kc)) * Bzamplitude + (omegae2 *
		omega * omega / (kc2 - omega * omega)) * VzamplitudeProton) / denominator;

	//double
	Byamplitude = (4 * pi * concentration * electron_charge_normalized / ((omega * omega / kc) - kc)) * (
		VzamplitudeElectron - VzamplitudeProton);

	//double
	VyamplitudeProton = -(Omegap / omega) * VzamplitudeProton - (electron_charge_normalized / (massProton * kc)) *
		Bzamplitude;
	//double
	VyamplitudeElectron = (Omegae / omega) * VzamplitudeElectron + (electron_charge_normalized / (massElectron * kc)) *
		Bzamplitude;

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

	B0 = Vector3d(B0.norm(), 0, 0);
	if (solverType == BUNEMAN) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEx[i][j][k] = 0;
					bunemanNewEx[i][j][k] = 0;
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEy[i][j][k] = Eyamplitude * cos(kw * xgrid[i] - phase);
					bunemanNewEy[i][j][k] = bunemanEy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanEz[i][j][k] = Ezamplitude * sin(kw * xgrid[i] - phase);
					bunemanNewEz[i][j][k] = bunemanEz[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBx[i][j][k] = B0.norm();
					bunemanNewBx[i][j][k] = bunemanBx[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBy[i][j][k] = Byamplitude * sin(kw * middleXgrid[i] - phase);
					bunemanNewBy[i][j][k] = bunemanBy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanBz[i][j][k] = Bzamplitude * cos(kw * middleXgrid[i] - phase);
					bunemanNewBz[i][j][k] = bunemanBz[i][j][k];
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					Efield[i][j][k].x = 0;
					Efield[i][j][k].y = Eyamplitude * cos(kw * xgrid[i] - phase);
					Efield[i][j][k].z = Ezamplitude * sin(kw * xgrid[i] - phase);
					explicitEfield[i][j][k] = Efield[i][j][k];
					tempEfield[i][j][k] = Efield[i][j][k];
					newEfield[i][j][k] = Efield[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					Bfield[i][j][k].x = B0.norm();
					Bfield[i][j][k].y = Byamplitude * sin(kw * middleXgrid[i] - phase);
					Bfield[i][j][k].z = Bzamplitude * cos(kw * middleXgrid[i] - phase);
					newBfield[i][j][k] = Bfield[i][j][k];
				}
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
		particle->prevMomentum = momentum;
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
		        4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron) / (
			        plasma_period * plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "c*rotBy amplitude = %g\n",
		        speed_of_light_normalized * kw * Bzamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "4*pi*Jz amplitude = %g\n",
		        4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron) / (
			        plasma_period * plasma_period * sqrt(
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
		        (speed_of_light_normalized * kw * Bzamplitude - 4 * pi * concentration * electron_charge_normalized * (
			        VyamplitudeProton - VyamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Ez amplitude = %g\n",
		        -omega * Ezamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "c*rotBz - 4*pi*Jz = %g\n",
		        (speed_of_light_normalized * kw * Byamplitude - 4 * pi * concentration * electron_charge_normalized * (
			        VzamplitudeProton - VzamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativJy = -electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jy amplitude = %g\n",
		        derivativJy / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesY = electron_charge_normalized * ((Eyamplitude * ((1.0 / massProton) + (1.0 / massElectron)
			)) + B0.norm() * ((VzamplitudeProton / massProton) + (VzamplitudeElectron / massElectron)) /
			speed_of_light_normalized
		);
		fprintf(informationFile, "dJy/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesY / (plasma_period * plasma_period *
			        plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "\n");
		double derivativJz = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jz amplitude = %g\n",
		        derivativJz / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesZ = electron_charge_normalized * ((Ezamplitude * ((1.0 / massProton) + (1.0 / massElectron)
			)) - B0.norm() * ((VyamplitudeProton / massProton) + (VyamplitudeElectron / massElectron)) /
			speed_of_light_normalized
		);
		fprintf(informationFile, "dJz/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesZ / (plasma_period * plasma_period *
			        plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativVyp = -omega * VyamplitudeProton;
		fprintf(informationFile, "-w*Vyp amplitude = %g\n", derivativVyp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonY = electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeProton /
			speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVyp/dt amplitude = %g\n",
		        derivativeVelocityProtonY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVzp = omega * VzamplitudeProton;
		fprintf(informationFile, "w*Vzp amplitude = %g\n", derivativVzp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonZ = electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeProton /
			speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVzp/dt amplitude = %g\n",
		        derivativeVelocityProtonZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVye = -omega * VyamplitudeElectron;
		fprintf(informationFile, "-w*Vye amplitude = %g\n", derivativVye * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronY = -electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeElectron /
			speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVye/dt amplitude = %g\n",
		        derivativeVelocityElectronY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVze = omega * VzamplitudeElectron;
		fprintf(informationFile, "w*Vze amplitude = %g\n", derivativVze * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronZ = -electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeElectron /
			speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVze/dt amplitude = %g\n",
		        derivativeVelocityElectronZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");
	}

	if (rank == 0) fclose(informationFile);
}


void Simulation::initializeAlfvenWaveY(int wavesCount, double amplitudeRelation) {
	boundaryConditionTypeX = PERIODIC;
	if (rank == 0)printf("initialization alfven wave\n");
	fflush(stdout);

	double concentration = density / (massProton + massElectron);
	types[1].particesDeltaX = types[0].particesDeltaX;
	types[1].particesDeltaY = types[0].particesDeltaY;
	types[1].particesDeltaZ = types[0].particesDeltaZ;
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
	VzamplitudeProton = -((1.0 / (4 * pi * concentration * electron_charge_normalized)) * (kc + ((omegae2 + omegap2 - omega
		* omega) / kc) + (omegae2 * Omegae2 / (kc * denominator))) / ((Omegae * omegae2 * omega / ((kc2 - omega * omega) *
		denominator)) + (Omegap / omega))) * Bzamplitude;
	//double
	VzamplitudeElectron = (((electron_charge_normalized * omega * Omegae) / (massElectron * kc)) * Bzamplitude + (omegae2 *
		omega * omega / (kc2 - omega * omega)) * VzamplitudeProton) / denominator;

	//double
	Byamplitude = (4 * pi * concentration * electron_charge_normalized / ((omega * omega / kc) - kc)) * (
		VzamplitudeElectron - VzamplitudeProton);

	//double
	VyamplitudeProton = -(Omegap / omega) * VzamplitudeProton - (electron_charge_normalized / (massProton * kc)) *
		Bzamplitude;
	//double
	VyamplitudeElectron = (Omegae / omega) * VzamplitudeElectron + (electron_charge_normalized / (massElectron * kc)) *
		Bzamplitude;

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

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				/*Efield[i][j][k].x = 0;
				Efield[i][j][k].y = Eyamplitude * cos(kw * xgrid[i] - kw * xshift);
				Efield[i][j][k].z = Ezamplitude * sin(kw * xgrid[i] - kw * xshift);*/
				Efield[i][j][k].x = Ezamplitude * sin(kw * ygrid[j] - kw * xshift);
				Efield[i][j][k].y = 0;
				Efield[i][j][k].z = Eyamplitude * cos(kw * ygrid[j] - kw * xshift);
				explicitEfield[i][j][k] = Efield[i][j][k];
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}

	/*if (nprocs == 1) {
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
	}*/

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				/*Bfield[i][j][k].x = B0.norm();
				Bfield[i][j][k].y = Byamplitude * sin(kw * middleXgrid[i] - kw * xshift);
				Bfield[i][j][k].z = Bzamplitude * cos(kw * middleXgrid[i] - kw * xshift);*/
				Bfield[i][j][k].x = Bzamplitude * cos(kw * middleYgrid[j] - kw * xshift);
				Bfield[i][j][k].y = B0.norm();
				Bfield[i][j][k].z = Byamplitude * sin(kw * middleYgrid[j] - kw * xshift);
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
			velocity = (Vector3d(1, 0, 0) * (VzamplitudeProton) * cos(
				kw * particle->coordinates.y - kw * xshift) + Vector3d(0, 0, 1) * VyamplitudeProton * sin(
				kw * particle->coordinates.y - kw * xshift));

		}
		if (particle->type == ELECTRON) {

			//velocity = Vector3d(0, 0, 1) * (VzamplitudeElectron) * (leftWeight*cos(kw * (ygrid[yn] - xshift)) + rightWeight*cos(kw*(ygrid[yn+1] - xshift))) + Vector3d(0, 1, 0) * VyamplitudeElectron * (leftWeight*sin(kw * (ygrid[yn] - xshift)) + rightWeight*sin(kw*(ygrid[yn+1] - xshift)));
			velocity = (Vector3d(1, 0, 0) * (VzamplitudeElectron) * cos(
				kw * particle->coordinates.y - kw * xshift) + Vector3d(0, 0, 1) * VyamplitudeElectron * sin(
				kw * particle->coordinates.y - kw * xshift));
		}
		double beta = velocity.norm() / speed_of_light_normalized;
		particle->addVelocity(velocity, speed_of_light_normalized);
		Vector3d momentum = particle->getMomentum();
		particle->initialMomentum = momentum;
		particle->prevMomentum = momentum;
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
		        4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron) / (
			        plasma_period * plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "c*rotBy amplitude = %g\n",
		        speed_of_light_normalized * kw * Bzamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "4*pi*Jz amplitude = %g\n",
		        4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron) / (
			        plasma_period * plasma_period * sqrt(
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
		        (speed_of_light_normalized * kw * Bzamplitude - 4 * pi * concentration * electron_charge_normalized * (
			        VyamplitudeProton - VyamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Ez amplitude = %g\n",
		        -omega * Ezamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "c*rotBz - 4*pi*Jz = %g\n",
		        (speed_of_light_normalized * kw * Byamplitude - 4 * pi * concentration * electron_charge_normalized * (
			        VzamplitudeProton - VzamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativJy = -electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jy amplitude = %g\n",
		        derivativJy / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesY = electron_charge_normalized * ((Eyamplitude * ((1.0 / massProton) + (1.0 / massElectron)
			)) + B0.norm() * ((VzamplitudeProton / massProton) + (VzamplitudeElectron / massElectron)) /
			speed_of_light_normalized
		);
		fprintf(informationFile, "dJy/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesY / (plasma_period * plasma_period *
			        plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "\n");
		double derivativJz = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jz amplitude = %g\n",
		        derivativJz / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesZ = electron_charge_normalized * ((Ezamplitude * ((1.0 / massProton) + (1.0 / massElectron)
			)) - B0.norm() * ((VyamplitudeProton / massProton) + (VyamplitudeElectron / massElectron)) /
			speed_of_light_normalized
		);
		fprintf(informationFile, "dJz/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesZ / (plasma_period * plasma_period *
			        plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativVyp = -omega * VyamplitudeProton;
		fprintf(informationFile, "-w*Vyp amplitude = %g\n", derivativVyp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonY = electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeProton /
			speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVyp/dt amplitude = %g\n",
		        derivativeVelocityProtonY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVzp = omega * VzamplitudeProton;
		fprintf(informationFile, "w*Vzp amplitude = %g\n", derivativVzp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonZ = electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeProton /
			speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVzp/dt amplitude = %g\n",
		        derivativeVelocityProtonZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVye = -omega * VyamplitudeElectron;
		fprintf(informationFile, "-w*Vye amplitude = %g\n", derivativVye * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronY = -electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeElectron /
			speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVye/dt amplitude = %g\n",
		        derivativeVelocityElectronY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVze = omega * VzamplitudeElectron;
		fprintf(informationFile, "w*Vze amplitude = %g\n", derivativVze * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronZ = -electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeElectron /
			speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVze/dt amplitude = %g\n",
		        derivativeVelocityElectronZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");
	}

	if (rank == 0) fclose(informationFile);
}

void Simulation::initializeAlfvenWaveZ(int wavesCount, double amplitudeRelation) {
	boundaryConditionTypeX = PERIODIC;
	if (rank == 0)printf("initialization alfven wave\n");
	fflush(stdout);

	double concentration = density / (massProton + massElectron);
	types[1].particesDeltaX = types[0].particesDeltaX;
	types[1].particesDeltaY = types[0].particesDeltaY;
	types[1].particesDeltaZ = types[0].particesDeltaZ;
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
	double kw = wavesCount * 2 * pi / zsizeGeneral;

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
	VzamplitudeProton = -((1.0 / (4 * pi * concentration * electron_charge_normalized)) * (kc + ((omegae2 + omegap2 - omega
		* omega) / kc) + (omegae2 * Omegae2 / (kc * denominator))) / ((Omegae * omegae2 * omega / ((kc2 - omega * omega) *
		denominator)) + (Omegap / omega))) * Bzamplitude;
	//double
	VzamplitudeElectron = (((electron_charge_normalized * omega * Omegae) / (massElectron * kc)) * Bzamplitude + (omegae2 *
		omega * omega / (kc2 - omega * omega)) * VzamplitudeProton) / denominator;

	//double
	Byamplitude = (4 * pi * concentration * electron_charge_normalized / ((omega * omega / kc) - kc)) * (
		VzamplitudeElectron - VzamplitudeProton);

	//double
	VyamplitudeProton = -(Omegap / omega) * VzamplitudeProton - (electron_charge_normalized / (massProton * kc)) *
		Bzamplitude;
	//double
	VyamplitudeElectron = (Omegae / omega) * VzamplitudeElectron + (electron_charge_normalized / (massElectron * kc)) *
		Bzamplitude;

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

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				/*Efield[i][j][k].x = 0;
				Efield[i][j][k].y = Eyamplitude * cos(kw * xgrid[i] - kw * xshift);
				Efield[i][j][k].z = Ezamplitude * sin(kw * xgrid[i] - kw * xshift);*/
				Efield[i][j][k].x = Eyamplitude * cos(kw * zgrid[k] - kw * xshift);
				Efield[i][j][k].z = 0;
				Efield[i][j][k].y = Ezamplitude * sin(kw * zgrid[k] - kw * xshift);
				explicitEfield[i][j][k] = Efield[i][j][k];
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}
	B0 = Vector3d(0, 0, B0.norm());
	/*if (nprocs == 1) {
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
	}*/

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				/*Bfield[i][j][k].x = B0.norm();
				Bfield[i][j][k].y = Byamplitude * sin(kw * middleXgrid[i] - kw * xshift);
				Bfield[i][j][k].z = Bzamplitude * cos(kw * middleXgrid[i] - kw * xshift);*/
				Bfield[i][j][k].x = Byamplitude * sin(kw * middleZgrid[k] - kw * xshift);
				Bfield[i][j][k].z = B0.norm();
				Bfield[i][j][k].y = Bzamplitude * cos(kw * middleZgrid[k] - kw * xshift);
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
			velocity = (Vector3d(0, 1, 0) * (VzamplitudeProton) * cos(
				kw * particle->coordinates.z - kw * xshift) + Vector3d(1, 0, 0) * VyamplitudeProton * sin(
				kw * particle->coordinates.z - kw * xshift));

		}
		if (particle->type == ELECTRON) {

			//velocity = Vector3d(0, 0, 1) * (VzamplitudeElectron) * (leftWeight*cos(kw * (ygrid[yn] - xshift)) + rightWeight*cos(kw*(ygrid[yn+1] - xshift))) + Vector3d(0, 1, 0) * VyamplitudeElectron * (leftWeight*sin(kw * (ygrid[yn] - xshift)) + rightWeight*sin(kw*(ygrid[yn+1] - xshift)));
			velocity = (Vector3d(0, 1, 0) * (VzamplitudeElectron) * cos(
				kw * particle->coordinates.z - kw * xshift) + Vector3d(1, 0, 0) * VyamplitudeElectron * sin(
				kw * particle->coordinates.z - kw * xshift));
		}
		double beta = velocity.norm() / speed_of_light_normalized;
		particle->addVelocity(velocity, speed_of_light_normalized);
		Vector3d momentum = particle->getMomentum();
		particle->initialMomentum = momentum;
		particle->prevMomentum = momentum;
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
		        4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron) / (
			        plasma_period * plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "c*rotBy amplitude = %g\n",
		        speed_of_light_normalized * kw * Bzamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "4*pi*Jz amplitude = %g\n",
		        4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron) / (
			        plasma_period * plasma_period * sqrt(
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
		        (speed_of_light_normalized * kw * Bzamplitude - 4 * pi * concentration * electron_charge_normalized * (
			        VyamplitudeProton - VyamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Ez amplitude = %g\n",
		        -omega * Ezamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "c*rotBz - 4*pi*Jz = %g\n",
		        (speed_of_light_normalized * kw * Byamplitude - 4 * pi * concentration * electron_charge_normalized * (
			        VzamplitudeProton - VzamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativJy = -electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jy amplitude = %g\n",
		        derivativJy / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesY = electron_charge_normalized * ((Eyamplitude * ((1.0 / massProton) + (1.0 / massElectron)
			)) + B0.norm() * ((VzamplitudeProton / massProton) + (VzamplitudeElectron / massElectron)) /
			speed_of_light_normalized
		);
		fprintf(informationFile, "dJy/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesY / (plasma_period * plasma_period *
			        plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "\n");
		double derivativJz = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jz amplitude = %g\n",
		        derivativJz / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesZ = electron_charge_normalized * ((Ezamplitude * ((1.0 / massProton) + (1.0 / massElectron)
			)) - B0.norm() * ((VyamplitudeProton / massProton) + (VyamplitudeElectron / massElectron)) /
			speed_of_light_normalized
		);
		fprintf(informationFile, "dJz/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesZ / (plasma_period * plasma_period *
			        plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativVyp = -omega * VyamplitudeProton;
		fprintf(informationFile, "-w*Vyp amplitude = %g\n", derivativVyp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonY = electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeProton /
			speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVyp/dt amplitude = %g\n",
		        derivativeVelocityProtonY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVzp = omega * VzamplitudeProton;
		fprintf(informationFile, "w*Vzp amplitude = %g\n", derivativVzp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonZ = electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeProton /
			speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVzp/dt amplitude = %g\n",
		        derivativeVelocityProtonZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVye = -omega * VyamplitudeElectron;
		fprintf(informationFile, "-w*Vye amplitude = %g\n", derivativVye * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronY = -electron_charge_normalized * (Eyamplitude + B0.norm() * VzamplitudeElectron /
			speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVye/dt amplitude = %g\n",
		        derivativeVelocityElectronY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVze = omega * VzamplitudeElectron;
		fprintf(informationFile, "w*Vze amplitude = %g\n", derivativVze * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronZ = -electron_charge_normalized * (Ezamplitude - B0.norm() * VyamplitudeElectron /
			speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVze/dt amplitude = %g\n",
		        derivativeVelocityElectronZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");
	}

	if (rank == 0) fclose(informationFile);
}

void Simulation::initializeRotatedAlfvenWave(int waveCountX, int waveCountY, int waveCountZ, double amplitudeRelation) {
	boundaryConditionTypeX = PERIODIC;
	if (rank == 0) printf("initialization alfven wave\n");
	fflush(stdout);

	double concentration = density / (massProton + massElectron);
	types[1].particesDeltaX = types[0].particesDeltaX;
	types[1].particesDeltaY = types[0].particesDeltaY;
	types[1].particesDeltaZ = types[0].particesDeltaZ;
	types[1].particlesPerBin = types[0].particlesPerBin;
	types[0].concentration = concentration;
	types[1].concentration = concentration;
	for (int i = 2; i < typesNumber; ++i) {
		types[i].particlesPerBin = 0;
		types[i].concentration = 0;
		types[i].particesDeltaX = xsize;
		types[i].particesDeltaY = ysize;
		types[i].particesDeltaZ = zsize;
	}
	createParticles();
	E0 = Vector3d(0, 0, 0);

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");

	double B0norm = B0.norm();
	double alfvenV = B0norm / sqrt(4 * pi * density);
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

	double kx = waveCountX * 2 * pi / xsizeGeneral;
	double ky = waveCountY * 2 * pi / ysizeGeneral;
	double kz = waveCountZ * 2 * pi / zsizeGeneral;
	//kz = 0;

	double kw = sqrt(kx * kx + ky * ky + kz * kz);

	double weight = concentration * volumeB() / types[0].particlesPerBin;

	omegaPlasmaProton = sqrt(
		4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massProton);
	omegaPlasmaElectron = sqrt(
		4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
	omegaGyroProton = B0norm * electron_charge_normalized / (massProton * speed_of_light_normalized);
	omegaGyroElectron = B0norm * electron_charge_normalized / (massElectron * speed_of_light_normalized);

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

	double fakeOmega = (kw * electron_charge_normalized * B0norm / massProton) * (sqrt(
		discriminant) - b) / (2.0 * a);

	//a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0

	double a4 = sqr(speed_of_light_normalized_sqr * massProton * massElectron);
	a4 = a4 * sqr(sqr(sqr(fakeOmega)));
	double a3 = -2 * cube(speed_of_light_normalized_sqr) * sqr(kw * massElectron * massProton)
		- 8 * pi * concentration * sqr(
			speed_of_light_normalized_sqr * electron_charge_normalized) * massElectron * massProton * (massElectron + massProton)
		- sqr(B0norm * speed_of_light_normalized * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron));
	a3 = a3 * cube(sqr(fakeOmega));
	double a2 = sqr(sqr(speed_of_light_normalized_sqr * kw) * massProton * massElectron)
		+ 8 * pi * cube(speed_of_light_normalized_sqr) * concentration * sqr(
			kw * electron_charge_normalized) * massProton * massElectron * (massProton + massElectron)
		+ 16 * sqr(pi * speed_of_light_normalized_sqr * concentration * sqr(
			electron_charge_normalized) * (massProton + massElectron))
		+ 2 * sqr(
			B0norm * speed_of_light_normalized_sqr * kw * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron))
		+ 8 * pi * concentration * sqr(B0norm * speed_of_light_normalized * sqr(
			electron_charge_normalized)) * (massProton + massElectron)
		+ sqr(sqr(B0norm * electron_charge_normalized));
	a2 = a2 * sqr(sqr(fakeOmega));
	double a1 = -sqr(
			B0norm * cube(speed_of_light_normalized) * kw * kw * electron_charge_normalized) * (sqr(
			massProton) + sqr(massElectron))
		- 8 * pi * concentration * sqr(B0norm * speed_of_light_normalized_sqr * kw * sqr(
			electron_charge_normalized)) * (massProton + massElectron)
		- 2 * sqr(speed_of_light_normalized * kw * sqr(B0norm * electron_charge_normalized));
	a1 = a1 * sqr(fakeOmega);
	double a0 = sqr(sqr(B0norm * speed_of_light_normalized * kw * electron_charge_normalized));

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
	Bzamplitude = B0norm * epsilonAmplitude;

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
	VzamplitudeProton = -((1.0 / (4 * pi * concentration * electron_charge_normalized)) * (kc + ((omegae2 + omegap2 - omega
		* omega) / kc) + (omegae2 * Omegae2 / (kc * denominator))) / ((Omegae * omegae2 * omega / ((kc2 - omega * omega) *
		denominator)) + (Omegap / omega))) * Bzamplitude;
	//double
	VzamplitudeElectron = (((electron_charge_normalized * omega * Omegae) / (massElectron * kc)) * Bzamplitude + (omegae2 *
		omega * omega / (kc2 - omega * omega)) * VzamplitudeProton) / denominator;

	//double
	Byamplitude = (4 * pi * concentration * electron_charge_normalized / ((omega * omega / kc) - kc)) * (
		VzamplitudeElectron - VzamplitudeProton);

	//double
	VyamplitudeProton = -(Omegap / omega) * VzamplitudeProton - (electron_charge_normalized / (massProton * kc)) *
		Bzamplitude;
	//double
	VyamplitudeElectron = (Omegae / omega) * VzamplitudeElectron + (electron_charge_normalized / (massElectron * kc)) *
		Bzamplitude;

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

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
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


	/*if (nprocs == 1) {
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
	}*/

	B0 = Vector3d(B0norm, 0, 0);
	B0 = rotationMatrix * B0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k].x = B0norm;
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
		        4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron) / (
			        plasma_period * plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "c*rotBy amplitude = %g\n",
		        speed_of_light_normalized * kw * Bzamplitude / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "4*pi*Jz amplitude = %g\n",
		        4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron) / (
			        plasma_period * plasma_period * sqrt(
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
		        (speed_of_light_normalized * kw * Bzamplitude - 4 * pi * concentration * electron_charge_normalized * (
			        VyamplitudeProton - VyamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");
		fprintf(informationFile, "derivative Ez amplitude = %g\n",
		        -omega * Ezamplitude / (plasma_period * plasma_period * sqrt(scaleFactor)));
		fprintf(informationFile, "c*rotBz - 4*pi*Jz = %g\n",
		        (speed_of_light_normalized * kw * Byamplitude - 4 * pi * concentration * electron_charge_normalized * (
			        VzamplitudeProton - VzamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
			        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativJy = -electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jy amplitude = %g\n",
		        derivativJy / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesY = electron_charge_normalized * ((Eyamplitude * ((1.0 / massProton) + (1.0 / massElectron)
		)) + B0norm * ((VzamplitudeProton / massProton) + (VzamplitudeElectron / massElectron)) / speed_of_light_normalized);
		fprintf(informationFile, "dJy/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesY / (plasma_period * plasma_period *
			        plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "\n");
		double derivativJz = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * omega;
		fprintf(informationFile, "w*Jz amplitude = %g\n",
		        derivativJz / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

		double derivativeVelocitiesZ = electron_charge_normalized * ((Ezamplitude * ((1.0 / massProton) + (1.0 / massElectron)
		)) - B0norm * ((VyamplitudeProton / massProton) + (VyamplitudeElectron / massElectron)) / speed_of_light_normalized);
		fprintf(informationFile, "dJz/dt amplitude = %g\n",
		        electron_charge_normalized * concentration * derivativeVelocitiesZ / (plasma_period * plasma_period *
			        plasma_period * sqrt(
				        scaleFactor)));
		fprintf(informationFile, "\n");

		double derivativVyp = -omega * VyamplitudeProton;
		fprintf(informationFile, "-w*Vyp amplitude = %g\n", derivativVyp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonY = electron_charge_normalized * (Eyamplitude + B0norm * VzamplitudeProton /
			speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVyp/dt amplitude = %g\n",
		        derivativeVelocityProtonY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVzp = omega * VzamplitudeProton;
		fprintf(informationFile, "w*Vzp amplitude = %g\n", derivativVzp * scaleFactor / sqr(plasma_period));

		double derivativeVelocityProtonZ = electron_charge_normalized * (Ezamplitude - B0norm * VyamplitudeProton /
			speed_of_light_normalized) / massProton;
		fprintf(informationFile, "dVzp/dt amplitude = %g\n",
		        derivativeVelocityProtonZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVye = -omega * VyamplitudeElectron;
		fprintf(informationFile, "-w*Vye amplitude = %g\n", derivativVye * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronY = -electron_charge_normalized * (Eyamplitude + B0norm * VzamplitudeElectron /
			speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVye/dt amplitude = %g\n",
		        derivativeVelocityElectronY * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");

		double derivativVze = omega * VzamplitudeElectron;
		fprintf(informationFile, "w*Vze amplitude = %g\n", derivativVze * scaleFactor / sqr(plasma_period));

		double derivativeVelocityElectronZ = -electron_charge_normalized * (Ezamplitude - B0norm * VyamplitudeElectron /
			speed_of_light_normalized) / massElectron;
		fprintf(informationFile, "dVze/dt amplitude = %g\n",
		        derivativeVelocityElectronZ * scaleFactor / sqr(plasma_period));
		fprintf(informationFile, "\n");
	}

	if (rank == 0) fclose(informationFile);
}

void Simulation::initializeLangmuirWave() {
	boundaryConditionTypeX = PERIODIC;
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
	double weight = (concentration / types[0].particlesPerBin) * volumeB();
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

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				Efield[i][j][k].x = Eamplitude * sin(kw * xgrid[i]);
				Efield[i][j][k].y = 0;
				Efield[i][j][k].z = 0;

				tempEfield[i] = Efield[i];
				explicitEfield[i] = Efield[i];
			}
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
	MPI_Barrier(cartComm);
}

void Simulation::initializeFluxFromRight() {
	//boundaryConditionTypeX = SUPER_CONDUCTOR_LEFT;
	boundaryConditionTypeX = PERIODIC;
	boundaryConditionTypeY = PERIODIC;
	boundaryConditionTypeZ = PERIODIC;
	createParticles();
	//E0 = E0 - V0.vectorMult(B0) / (speed_of_light_normalized * speed_of_light_correction);
	E0 = E0 - V0.vectorMult(B0) / (speed_of_light_normalized);
	rightBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(E0, B0);
	//initializeAlfvenWaveY(10, 1.0E-4);
	if (solverType == BUNEMAN) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEx[i][j][k] = E0.x;
					bunemanNewEx[i][j][k] = bunemanEx[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEy[i][j][k] = E0.y;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i < 1 + additionalBinNumber) {
							bunemanEy[i][j][k] = 0;
						}
					}
					bunemanNewEy[i][j][k] = bunemanEy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanEz[i][j][k] = E0.z;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i < 1 + additionalBinNumber) {
							bunemanEz[i][j][k] = 0;
						}
					}
					bunemanNewEz[i][j][k] = bunemanEz[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBx[i][j][k] = B0.x;
					bunemanNewBx[i][j][k] = bunemanBx[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBy[i][j][k] = B0.y;
					bunemanNewBy[i][j][k] = bunemanBy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanBz[i][j][k] = B0.z;
					bunemanNewBz[i][j][k] = bunemanBz[i][j][k];
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					Efield[i][j][k] = E0;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i <= 1 + additionalBinNumber) {
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
	}

	double gamma = 1.0/sqrt(1 - V0.scalarMult(V0)/speed_of_light_normalized_sqr);
	double p0 = gamma*massProton*V0.norm();
	double protonGyroRadius = p0*speed_of_light_normalized/(fabs(electron_charge_normalized)*B0.norm());
	int countGyroRadius = xsizeGeneral/protonGyroRadius;
	double deltaK = 2*pi/xsizeGeneral;
	double minK = deltaK;
	double maxK = 2*pi/deltaX;
	int maxCount = min2(2*countGyroRadius, xnumberGeneral);
	int minCount = max2(1, countGyroRadius/2);
	int count = maxCount - minCount + 1;

	//initializeRandomModes(count, minCount, 0.5);

	double magneticEnergy = B0.scalarMult(B0) / (8 * pi);
	double kineticEnergy = density * V0.scalarMult(V0) / 2;

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");
	if (rank == 0) fprintf(informationFile, "magneticEnergy/kineticEnergy = %15.10g\n", magneticEnergy / kineticEnergy);
	if (rank == 0) fprintf(informationFile, "magneticEnergy/kinetikEnergy = %15.10g\n", magneticEnergy / kineticEnergy);
	if (rank == 0) fprintf(informationFile, "protonGyroRadius = %15.10g\n", protonGyroRadius);
	if (rank == 0) fprintf(informationFile, "countGyroRadius = %d\n", countGyroRadius);
	if (rank == 0) fprintf(informationFile, "minCount = %d\n", minCount);
	if (rank == 0) fprintf(informationFile, "maxCount = %d\n", maxCount);
	fflush(stdout);
	if (rank == 0) fclose(informationFile);

	checkDebyeParameter();

	MPI_Barrier(cartComm);
	if (rank == 0) {
		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", 0.0, 0.0, 1.0, 1.0);
		fclose(incrementFile);
	}
	MPI_Barrier(cartComm);
}

void Simulation::initializeHarris() {
	boundaryConditionTypeX = FREE_BOTH;
	//boundaryConditionTypeX = PERIODIC;
	boundaryConditionTypeY = PERIODIC;
	boundaryConditionTypeZ = PERIODIC;
	
	E0 = Vector3d(0, 0, 0);
	double harrisWidth = 10*deltaX;
	createParticlesHarris(harrisWidth);

	if(verbosity > 2) printf("initialize harris boundary evaluator rank = %d\n", rank);

	rightBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(E0, B0);
	Vector3d leftB = B0*(-1);
	leftBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(E0, leftB);

	if (solverType == BUNEMAN) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEx[i][j][k] = E0.x;
					bunemanNewEx[i][j][k] = bunemanEx[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEy[i][j][k] = E0.y;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i < 1 + additionalBinNumber) {
							bunemanEy[i][j][k] = 0;
						}
					}
					bunemanNewEy[i][j][k] = bunemanEy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanEz[i][j][k] = E0.z;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i < 1 + additionalBinNumber) {
							bunemanEz[i][j][k] = 0;
						}
					}
					bunemanNewEz[i][j][k] = bunemanEz[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBx[i][j][k] = B0.x;
					bunemanNewBx[i][j][k] = bunemanBx[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBy[i][j][k] = B0.y;
					bunemanNewBy[i][j][k] = bunemanBy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanBz[i][j][k] = B0.z;
					bunemanNewBz[i][j][k] = bunemanBz[i][j][k];
				}
			}
		}
	} else {
		if(verbosity > 2) printf("initialize harris E field rank = %d\n", rank);
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					Efield[i][j][k] = E0;
					tempEfield[i][j][k] = Efield[i][j][k];
					newEfield[i][j][k] = Efield[i][j][k];
					explicitEfield[i][j][k] = Efield[i][j][k];
				}
			}
		}

		if(verbosity > 2) printf("initialize harris B field rank = %d\n", rank);

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					Bfield[i][j][k] = B0*tanh((xgrid[i] - 1.5*xsizeGeneral)/harrisWidth);
					newBfield[i][j][k] = Bfield[i][j][k];
				}
			}
		}
	}

	fflush(stdout);

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");
	if (rank == 0) fclose(informationFile);

	checkDebyeParameter();

	MPI_Barrier(cartComm);
	if (rank == 0) {
		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", 0.0, 0.0, 1.0, 1.0);
		fclose(incrementFile);
	}
	MPI_Barrier(cartComm);
}

void Simulation::fieldsLorentzTransitionX(const double& v) {
	double gamma = 1.0 / sqrt(1 - v * v / speed_of_light_normalized_sqr);
	for (int i = 0; i < xnumberAdded; ++i) {
		int prevI = i - 1;
		if (prevI < 0) {
			if (boundaryConditionTypeX == PERIODIC) {
				prevI = xnumberAdded - 1;
			} else {
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
				Vector3d middleB = (Bfield[prevI][prevJ][prevK] + Bfield[prevI][j][prevK] + Bfield[prevI][prevJ][k] + Bfield[prevI][
						j][k]
					+ Bfield[i][prevJ][prevK] + Bfield[i][j][prevK] + Bfield[i][prevJ][k] + Bfield[i][j][k]) * 0.125;
				newEfield[i][j][k].y = gamma * (Efield[i][j][k].y - v * middleB.z / speed_of_light_normalized);
				newEfield[i][j][k].z = gamma * (Efield[i][j][k].z + v * middleB.y / speed_of_light_normalized);
			}
		}
	}
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Vector3d middleE = (Efield[i][j][k] + Efield[i][j + 1][k] + Efield[i][j][k + 1] + Efield[i][j + 1][k + 1]
					+ Efield[i + 1][j][k] + Efield[i + 1][j + 1][k] + Efield[i + 1][j][k + 1] + Efield[i + 1][j + 1][k + 1]) * 0.125;
				newBfield[i][j][k].y = gamma * (Bfield[i][j][k].y + v * middleE.z / speed_of_light_normalized);
				newBfield[i][j][k].z = gamma * (Bfield[i][j][k].z - v * middleE.y / speed_of_light_normalized);
			}
		}
	}

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				Efield[i][j][k] = newEfield[i][j][k];
				tempEfield[i][j][k] = newEfield[i][j][k];
				explicitEfield[i][j][k] = newEfield[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k] = newBfield[i][j][k];
			}
		}
	}
}

void Simulation::initializeShockWave() {
	boundaryConditionTypeX = FREE_BOTH;

	E0 = E0 - V0.vectorMult(B0) / (speed_of_light_normalized);
	double frontWidth = 500 * deltaX;
	Vector3d downstreamB = Vector3d(0, 0, 0);
	Vector3d downstreamE = Vector3d(0, 0, 0);
	double downstreamDensity;
	Vector3d downstreamVelocity = Vector3d(0, 0, 0);
	double compressionRatio = 3.9;
	double downstreamPressure;
	double upstreamPressure = evaluatePressureByTemperature(temperature);
	solveRankineHugoniot(density, V0, upstreamPressure, B0, E0, downstreamDensity, downstreamVelocity, downstreamPressure,
	                     downstreamB, downstreamE, 5.0 / 3.0, compressionRatio);
	double downstreamTemperature = evaluateTemperatureByPressure(downstreamPressure) / compressionRatio;
	double velocityForRotB = -(downstreamB.z - B0.z) * speed_of_light_normalized / (8 * pi * frontWidth *
		electron_charge_normalized * types[0].concentration);
	int localI = -1;
	int shockWaveRank[1];
	int tempShockWaveRank[1];
	double shockWaveX[1];
	shockWaveX[0] = 0;
	shockWaveRank[0] = -1;
	shockWavePoint = xnumberGeneral / 3;
	for (int i = 0; i < xnumberAdded - 1; ++i) {
		if (firstAbsoluteXindex + i == shockWavePoint) {
			if (i < 3) {
				shockWavePoint += 3;
				localI = i + 3;
			} else {
				localI = i;
			}
			if (localI < 10) {
				localI += 10;
				shockWavePoint += 10;
			}
			if (localI > xnumberAdded - 10) {
				localI -= 10;
				shockWavePoint -= 10;
			}
			shockWaveRank[0] = rank;
			break;
		}
	}
	if (rank > 0) {
		MPI_Send(shockWaveRank, 1, MPI_INT, 0, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm);
	} else {
		tempShockWaveRank[0] = shockWaveRank[0];
		for (int i = 1; i < nprocs; ++i) {
			MPI_Status status;
			MPI_Recv(tempShockWaveRank, 1, MPI_INT, i, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm, &status);
			if (tempShockWaveRank[0] > -1) {
				shockWaveRank[0] = tempShockWaveRank[0];
			}
		}
	}

	MPI_Bcast(shockWaveRank, 1, MPI_INT, 0, cartComm);
	if (rank == shockWaveRank[0]) {
		shockWaveX[0] = xgrid[localI];
	}
	MPI_Bcast(shockWaveX, 1, MPI_DOUBLE, shockWaveRank[0], cartComm);

	int lnumber = 3 * typesNumber + 3;

	double**** outVector = new double***[xnumberAdded];
	double**** tempVector = new double***[xnumberAdded];
	for (int i = 0; i < xnumberAdded; ++i) {
		outVector[i] = new double**[ynumberAdded];
		tempVector[i] = new double**[ynumberAdded];
		for (int j = 0; j < ynumberAdded; ++j) {
			outVector[i][j] = new double* [znumberAdded];
			tempVector[i][j] = new double* [znumberAdded];
			for (int k = 0; k < znumberAdded; ++k) {
				outVector[i][j][k] = new double[lnumber];
				tempVector[i][j][k] = new double[lnumber];
				for (int l = 0; l < lnumber; ++l) {
					outVector[i][j][k][l] = 0;
					tempVector[i][j][k][l] = 0;
				}
			}
		}
	}

	double density = 0;
	for (int i = 0; i < typesNumber; ++i) {
		density += types[i].concentration * types[i].mass;
	}
	double bnorm = B0.norm();
	double alfvenV = B0.norm() / sqrt(4 * pi * density);

	RightPartIterativeEvaluator* evaluator = new InitializeShockRightPartEvaluator(
		downstreamVelocity, V0, downstreamB, B0, E0, deltaX, speed_of_light_normalized, bnorm, alfvenV, types, typesNumber,
		xnumberAdded, 1, 1, lnumber, rank, nprocs, cartComm, cartCoord, cartDim);

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int t = 0; t < typesNumber; ++t) {
					outVector[i][0][0][3 * t] = (0.5 * (downstreamVelocity.x + V0.x) - 0.5 * (downstreamVelocity.x - V0.x) * tanh(
						(xgrid[i] - shockWaveX[0]) / frontWidth)) / alfvenV;
					outVector[i][0][0][3 * t + 1] = (0.5 * (downstreamVelocity.y + V0.y) - 0.5 * (downstreamVelocity.y - V0.y) * tanh(
						(xgrid[i] - shockWaveX[0]) / frontWidth)) / alfvenV;
					outVector[i][0][0][3 * t + 2] = (0.5 * (downstreamVelocity.z + V0.z) - 0.5 * (downstreamVelocity.z - V0.z) * tanh(
						(xgrid[i] - shockWaveX[0]) / frontWidth)) / alfvenV;
				}
				outVector[i][0][0][lnumber - 3] = E0.x / bnorm;
				outVector[i][0][0][lnumber - 2] = (0.5 * (downstreamB.y + B0.y) - 0.5 * (downstreamB.y - B0.y) * tanh(
					(xgrid[i] - shockWaveX[0]) / frontWidth)) / bnorm;
				outVector[i][0][0][lnumber - 1] = (0.5 * (downstreamB.z + B0.z) - 0.5 * (downstreamB.z - B0.z) * tanh(
					(xgrid[i] - shockWaveX[0]) / frontWidth)) / bnorm;
			}
		}
	}

	simpleIterationSolver(outVector, tempVector, xnumberAdded, 1, 1, additionalBinNumber, lnumber, rank, nprocs, xnumberGeneral,
	                      ynumberGeneral, znumberGeneral, maxErrorLevel, maxSimpleIterationSolverIterations,
	                      boundaryConditionTypeX == PERIODIC, boundaryConditionTypeY == PERIODIC,
	                      boundaryConditionTypeZ == PERIODIC, verbosity, leftOutGmresBuffer, rightOutGmresBuffer,
	                      leftInGmresBuffer, rightInGmresBuffer, frontOutGmresBuffer, backOutGmresBuffer,
	                      frontInGmresBuffer, backInGmresBuffer, bottomOutGmresBuffer, topOutGmresBuffer,
	                      bottomInGmresBuffer, topInGmresBuffer, evaluator, cartComm, cartCoord, cartDim);


	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				Efield[i][j][k] = E0;
				Efield[i][j][k].x = outVector[i][0][0][lnumber - 3] * bnorm;
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}
	exchangeGeneralEfield(Efield);
	exchangeGeneralEfield(tempEfield);
	exchangeGeneralEfield(newEfield);
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k] = B0;
				Bfield[i][j][k].y = 0.5 * (outVector[i][0][0][lnumber - 2] + outVector[i + 1][0][0][lnumber - 2]) * bnorm;
				Bfield[i][j][k].z = 0.5 * (outVector[i][0][0][lnumber - 1] + outVector[i + 1][0][0][lnumber - 1]) * bnorm;
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
	exchangeGeneralBfield(Bfield);
	exchangeGeneralBfield(newBfield);

	rightBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(Efield[xnumberAdded - additionalBinNumber][1+additionalBinNumber][1+additionalBinNumber], Bfield[xnumberAdded = additionalBinNumber][1+additionalBinNumber][1+additionalBinNumber]);
	leftBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(Efield[1+additionalBinNumber][1+additionalBinNumber][1+additionalBinNumber], Bfield[1+additionalBinNumber][1+additionalBinNumber][1+additionalBinNumber]);

	if (rank == 0) printf("creating particles\n");

	int n = 0;
	for (int i = 1; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				//int maxParticlesPerBin = types[0].particlesPerBin;
				double x = xgrid[i] + 0.0001 * deltaX;
				double y = ygrid[j] + 0.0001 * deltaY;
				double z = zgrid[k] + 0.0001 * deltaZ;
				//for (int l = 0; l < maxParticlesPerBin; ++l) {
				for (int typeCounter = 0; typeCounter < typesNumber; ++typeCounter) {
					double localConcentration;
					double localWeight;
					double localTemperature;
					int localParticlesPerBin;
					double localParticleDeltaX;
					double localParticleDeltaY;
					double localParticleDeltaZ;
					Vector3d localVelocity;
					if (types[typeCounter].concentration > 0) {
						localConcentration = types[typeCounter].concentration * V0.x / (0.5 * alfvenV * (outVector[i][0][0][3 *
							typeCounter] + outVector[i + 1][0][0][3 * typeCounter]));
						localParticlesPerBin = types[typeCounter].particlesPerBin * V0.x / (0.5 * alfvenV * (outVector[i][0][0][3 *
							typeCounter] + outVector[i + 1][0][0][3 * typeCounter]));
					} else {
						localConcentration = 0;
						localParticlesPerBin = 0;
					}
					localWeight = (localConcentration / localParticlesPerBin) * volumeB();
					if (types[typeCounter].particlesPerBin > 0) {
						localParticleDeltaX = deltaX / localParticlesPerBin;
						localParticleDeltaY = deltaY / localParticlesPerBin;
						localParticleDeltaZ = deltaZ / localParticlesPerBin;
					} else {
						localParticleDeltaX = xsize;
						localParticleDeltaY = ysize;
						localParticleDeltaZ = zsize;
					}
					localTemperature = (downstreamTemperature + types[typeCounter].temperatureX) * 0.5 - (downstreamTemperature - types
						[typeCounter].temperatureX) * 0.5 * tanh((middleXgrid[i] - shockWaveX[0]) / frontWidth);
					localVelocity.x = 0.5 * (outVector[i][0][0][3 * typeCounter] + outVector[i + 1][0][0][3 * typeCounter]) * alfvenV;
					localVelocity.y = 0.5 * (outVector[i][0][0][3 * typeCounter + 1] + outVector[i + 1][0][0][3 * typeCounter + 1]) *
						alfvenV;
					localVelocity.z = 0.5 * (outVector[i][0][0][3 * typeCounter + 2] + outVector[i + 1][0][0][3 * typeCounter + 2]) *
						alfvenV;
					//if (l < types[typeCounter].particlesPerBin) {
					for (int l = 0; l < localParticlesPerBin; ++l) {
						ParticleTypes type = types[typeCounter].type;
						Particle* particle = createParticle(n, i, j, k, localWeight, type, types[typeCounter],
						                                    localTemperature,
						                                    localTemperature,
						                                    localTemperature);
						n++;
						particle->coordinates.x = x + localParticleDeltaX * l;
						particle->coordinates.y = y + localParticleDeltaY * l;
						particle->coordinates.z = z + localParticleDeltaZ * l;
						particle->initialCoordinates = particle->coordinates;
						Vector3d particleVelocity = localVelocity;
						if (typeCounter == 0) {
							double denominator = (sqr(cosh((particle->coordinates.x - shockWaveX[0]) / frontWidth)) * ((compressionRatio + 1)
								* 0.5 - (compressionRatio - 1) * 0.5 * tanh(particle->coordinates.x - shockWaveX[0]) / frontWidth));
							if (fabs(denominator) < 1.0) {
								printf("aaa\n");
							}
							particleVelocity.y += velocityForRotB / denominator;
							if (particleVelocity.norm() > speed_of_light_normalized) {
								printf("aaa\n");
							}
						}
						if (localVelocity.norm() > 0) {
							particle->addVelocity(particleVelocity, speed_of_light_normalized);
						}
						Vector3d momentum = particle->getMomentum();
						particle->initialMomentum = momentum;
						//particle->prevMomentum = momentum;

						particles.push_back(particle);
						particlesNumber++;
						if (particlesNumber % 1000 == 0) {
							if ((rank == 0) && (verbosity > 0))printf("create particle number %d\n", particlesNumber);
						}
						alertNaNOrInfinity(particle->coordinates.x, "particle.x = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.y, "particle.y = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.z, "particle.z = NaN in createParticles\n");
					}
				}
			}
		}
	}


	synchronizeParticleNumber();

	delete evaluator;
	for (int i = 0; i < xnumberAdded; ++i) {;
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				delete[] outVector[i][j][k];
				delete[] tempVector[i][j][k];
			}
			delete[] outVector[i][j];
			delete[] tempVector[i][j];
		}
		delete[] outVector[i];
		delete[] tempVector[i];
	}
	delete[] outVector;
	delete[] tempVector;
}

void Simulation::solveRankineHugoniot(double upstreamDensity, Vector3d upstreamVelocity, double upstreamPressure,
                                      Vector3d upstreamB, Vector3d upstreamE, double& downstreamDensity,
                                      Vector3d& downstreamVelocity, double& downstreamPressure, Vector3d& downstreamB,
                                      Vector3d& downstreamE, double adiabaticParameter, double compressionRatio) {
	double j = upstreamDensity * upstreamVelocity.x;
	double magneticE = upstreamB.x * upstreamB.x / (4 * pi);
	downstreamDensity = compressionRatio * upstreamDensity;
	downstreamVelocity.x = upstreamVelocity.x / compressionRatio;
	downstreamB.x = upstreamB.x;
	double Brelation = ((j * j / upstreamDensity) - magneticE) / ((j * j / downstreamDensity) - magneticE);
	downstreamB.y = upstreamB.y * Brelation;
	downstreamB.z = upstreamB.z * Brelation;

	downstreamVelocity.y = upstreamVelocity.y + upstreamB.x * (downstreamB.y - upstreamB.y) / (4 * pi * j);
	downstreamVelocity.z = upstreamVelocity.z + upstreamB.x * (downstreamB.z - upstreamB.z) / (4 * pi * j);

	downstreamE = upstreamE;

	double a = (1.0 / ((adiabaticParameter - 1) * downstreamDensity)) + (0.5 / downstreamDensity) - (0.5 / upstreamDensity
	);
	double b = (1.0 / ((adiabaticParameter - 1) * upstreamDensity)) + (0.5 / upstreamDensity) - (0.5 / downstreamDensity);
	double c = ((1 / downstreamDensity) - (1 / upstreamDensity)) * (sqr(downstreamB.y) + sqr(downstreamB.z) -
		sqr(upstreamB.y) - sqr(upstreamB.z)) / (16 * pi);

	downstreamPressure = (b * upstreamPressure - c) / a;
}

double Simulation::evaluateTemperatureByPressure(double& pressure) {
	double totalConcentration = 0;
	for (int t = 0; t < typesNumber; ++t) {
		totalConcentration += types[t].concentration;
	}
	double result = pressure / (kBoltzman_normalized * totalConcentration);
	return result;
}

double Simulation::evaluatePressureByTemperature(double& temperature) {
	double pressure = 0;
	for (int t = 0; t < typesNumber; ++t) {
		pressure += types[t].concentration * kBoltzman_normalized * temperature;
	}
	return pressure;
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
	double* amplitudes = new double[2*(last - first + 1)];
	double* knumbers = new double[(last - first + 1)];
	double* omega = new double[(last - first + 1)];

	if (rank == 0) {
		for (int i = 0; i < 2 * (last - first + 1); ++i) {
			phases[i] = 2 * pi * uniformDistribution();
		}
	}


	MPI_Barrier(cartComm);
	MPI_Bcast(phases, 2 * (last - first + 1), MPI_DOUBLE, 0, cartComm);
	MPI_Barrier(cartComm);

	for (int harmCounter = first; harmCounter <= last; ++harmCounter) {
		double kw = 2 * pi * harmCounter / length;
		int l = harmCounter - first;
		double Bamplitude = amplitude * power(kw, -5.0 / 6.0);
		///double phiY = 2 * pi * uniformDistribution();
		//double phiZ = 2 * pi * uniformDistribution();
		knumbers[l] = kw;
		omega[l] = 0;
		amplitudes[2*l] = Bamplitude;
		amplitudes[2*l+1] = Bamplitude;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					Bfield[i][j][k].y += Bamplitude * sin(kw * middleXgrid[i] + phases[2 * harmCounter]);
					Bfield[i][j][k].z += Bamplitude * cos(kw * middleXgrid[i] + phases[2 * harmCounter + 1]);
					newBfield[i][j][k] = Bfield[i][j][k];
				}
			}
		}
	}

	rightBoundaryFieldEvaluator = new TurbulenceBoundaryFieldEvaluator(E0, B0, V0, last-first+1, amplitudes, phases, knumbers, omega, xgrid[xnumberAdded - additionalBinNumber], speed_of_light_normalized);
	leftBoundaryFieldEvaluator = new TurbulenceBoundaryFieldEvaluator(E0, B0, V0, last-first+1, amplitudes, phases, knumbers, omega, xgrid[1 + additionalBinNumber], speed_of_light_normalized);

	delete[] phases;
	delete[] amplitudes;
	delete[] knumbers;
	delete[] omega;
}


void Simulation::initializeRandomModes(int number, int minNumber, double energyFraction) {
	//use if defined shockWavePoint

	double deltaK = 2*pi/xsizeGeneral;
	double minK = minNumber*deltaK;
	double maxK = minK*number;

	double constFieldFraction = sqrt(1.0 - energyFraction);


	double energy = energyFraction * B0.scalarMult(B0);
	double amplitude = sqrt(2*energy/number);

	double* phases = new double[2 * number];
	double* amplitudes = new double[2*number];
	double* knumbers = new double[number];
	double* omega = new double[number];

	if (rank == 0) {
		for (int i = 0; i < 2*number; ++i) {
			phases[i] = 2 * pi * uniformDistribution();
			//phases[i] = 0;
		}
	}


	MPI_Barrier(cartComm);
	MPI_Bcast(phases, number, MPI_DOUBLE, 0, cartComm);
	MPI_Barrier(cartComm);
	
	B0 = B0*constFieldFraction;
	E0 = E0*constFieldFraction;

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k] = B0;
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}

	for (int l = 0; l < number; ++l) {
		double kw = minK + l*deltaK;
		knumbers[l] = kw;
		omega[l] = 0;
		double Bamplitude = amplitude;
		amplitudes[2*l] = 0;
		amplitudes[2*l + 1] = amplitude;
		///double phiY = 2 * pi * uniformDistribution();
		//double phiZ = 2 * pi * uniformDistribution();

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					//Bfield[i][j][k].y += Bamplitude * sin(kw * middleXgrid[i] + phases[2 * harmCounter]);
					Bfield[i][j][k].z += Bamplitude * sin(kw * middleXgrid[i] + phases[2 * l + 1]);
					newBfield[i][j][k] = Bfield[i][j][k];
				}
			}
		}
	}

	for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					Efield[i][j][k] = E0;
					if(i > 0 && i < xnumberAdded && j < ynumberAdded && k < znumberAdded) {
						Vector3d B = (Bfield[i-1][j][k] + Bfield[i][j][k])*0.5;
						Efield[i][j][k] = V0.vectorMult(B)/(speed_of_light_normalized * speed_of_light_correction*(-1.0));
					}
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i <= 1 + additionalBinNumber) {
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

	//rightBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(E0, B0);
	//leftBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(E0, B0);
	rightBoundaryFieldEvaluator = new TurbulenceBoundaryFieldEvaluator(E0, B0, V0, number, amplitudes, phases, knumbers, omega, middleXgrid[xnumberAdded - additionalBinNumber], speed_of_light_normalized);
	leftBoundaryFieldEvaluator = new TurbulenceBoundaryFieldEvaluator(E0, B0, V0, number, amplitudes, phases, knumbers, omega, xgrid[1 + additionalBinNumber], speed_of_light_normalized);

	delete[] phases;
	delete[] amplitudes;
	delete[] knumbers;
	delete[] omega;
}

void Simulation::initializeTwoStream() {
	boundaryConditionTypeX = PERIODIC;
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
	MPI_Barrier(cartComm);

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
			} else {
				particle->addVelocity(electronsVelocityMinus, speed_of_light_normalized);
			}
			electronCount++;
		}
	}
}

void Simulation::initializeExternalFluxInstability() {
	boundaryConditionTypeX = PERIODIC;
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
		} else if (omega > cyclothronOmegaProton / 100.0) {
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
	boundaryConditionTypeX = PERIODIC;
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
	MPI_Barrier(cartComm);
}

void Simulation::initializeAnisotropicSilicon() {
	boundaryConditionTypeX = PERIODIC;
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
	MPI_Barrier(cartComm);
}


void Simulation::initializeWeibel() {
	boundaryConditionTypeX = PERIODIC;


	types[0].alphaNormal = 10;
	//types[0].alphaNormal = 1;
	//types[0].alphaParallel = 1.3*types[0].alphaNormal;
	types[0].alphaParallel = 1.1 * types[0].alphaNormal;

	types[0].temperatureX = ((types[0].mass * speed_of_light_normalized_sqr / types[0].alphaParallel) / (1 + (types[0].
		alphaParallel / types[0].alphaNormal - 1) * McDonaldFunction(
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

		double k02 = (omegaPlasmaElectron * omegaPlasmaElectron / speed_of_light_normalized_sqr) * (alphaParallel /
				alphaNormal - 1)
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
			1 - kmaxIncrement * kmaxIncrement / k02) * (kBoltzman_normalized * types[0].temperatureX / (types[0].mass *
			speed_of_light_normalized)) * cube(
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
	MPI_Barrier(cartComm);
	if (rank == 0) printf("finish initialize weibel\n");
}


void Simulation::initializeRingWeibel() {
	boundaryConditionTypeX = PERIODIC;

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


		double k02 = (omegaPlasmaElectron * omegaPlasmaElectron / (gamma * speed_of_light_normalized_sqr)) * ((betaNormal *
			betaNormal / (2 * betaParallel * betaParallel * (1 - betaParallel * betaParallel))) - G);
		if (k02 < 0) {
			increment = 0;
			length = 0;
		} else {

			double a = omegaPlasmaElectron * omegaPlasmaElectron * betaNormal * betaNormal / (2 * gamma * betaParallel *
				betaParallel);

			double bParallel2 = betaParallel * betaParallel;

			double kmaxIncrement2 = ((speed_of_light_normalized_sqr * k02 - a - bParallel2 * (speed_of_light_normalized_sqr * k02
						+ a)) +
					sqrt(a * sqr(
						bParallel2 + 1) * (bParallel2 * speed_of_light_normalized_sqr * k02 - speed_of_light_normalized_sqr * k02 + a))) /
				sqr(speed_of_light_normalized * (betaParallel + 1) * (betaParallel - 1));

			if (kmaxIncrement2 < 0) {
				increment = 0;
				length = 0;
			} else {

				double kmaxIncrement = sqrt(kmaxIncrement2);
				length = 2 * pi / kmaxIncrement;


				increment = (1 / sqrt(2.0)) * sqrt(
					sqrt(
						sqr(speed_of_light_normalized_sqr * kmaxIncrement2 * bParallel2 + a - speed_of_light_normalized_sqr * (k02 -
							kmaxIncrement2)) + 4 * speed_of_light_normalized_sqr * speed_of_light_normalized_sqr * kmaxIncrement2 *
						bParallel2 * (k02 - kmaxIncrement2)) -
					(speed_of_light_normalized_sqr * kmaxIncrement2 * bParallel2 + a - speed_of_light_normalized_sqr * (k02 -
						kmaxIncrement2))
				);
			}
		}

		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", increment, increment / plasma_period, length, length / scaleFactor);
		fclose(incrementFile);
	}
	MPI_Barrier(cartComm);
}

void Simulation::initializeHomogenouseFlow() {
	boundaryConditionTypeX = FREE_BOTH;
	boundaryConditionTypeX = SUPER_CONDUCTOR_LEFT;
	boundaryConditionTypeX = PERIODIC;
	createParticles();
	E0 = E0 - V0.vectorMult(B0) / (speed_of_light_normalized);
	//initializeAlfvenWaveY(10, 1.0E-4);
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				Efield[i][j][k] = E0;
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}

	/*for (int p = 0; p < particles.size(); ++p) {
		Particle* particle = particles[p];
		Vector3d momentum = V0 * particle->mass / sqrt(1 - V0.scalarMult(V0) / speed_of_light_normalized_sqr);
		particle->setMomentum(momentum);
	}*/
}

void Simulation::initializeFake() {
	boundaryConditionTypeX = PERIODIC;
	createParticles();
	//E0.y = -B0.z;
	E0.y = B0.x;
	E0.x = 0;
	E0.z = 0;
	double kw = 2 * pi / ysizeGeneral;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				Efield[i][j][k] = E0 * cos(kw * ygrid[j]);
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = Efield[i][j][k];
				explicitEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k] = B0;
				//Bfield[i][j][k] = B0*(middleYgrid[j] - ysizeGeneral)/ysizeGeneral;
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
}

void Simulation::createArrays() {
	if (rank == 0) printf("creating arrays\n");
	if (rank == 0) fflush(stdout);
	//if(rank == 0) printLog("creating arrays\n");

	if (rank == 0) printf("creating grid arrays\n");
	if (rank == 0) fflush(stdout);
	// if(rank == 0) printLog("creating grid arrays\n");
	xgrid = new double[xnumberAdded + 1];
	ygrid = new double[ynumberAdded + 1];
	zgrid = new double[znumberAdded + 1];

	middleXgrid = new double[xnumberAdded];
	middleYgrid = new double[ynumberAdded];
	middleZgrid = new double[znumberAdded];

	if (rank == 0) printf("creating gmresOutput arrays\n");
	if (rank == 0) fflush(stdout);

	gmresOutput = create4array(xnumberAdded, ynumberAdded, znumberAdded, maxwellEquationMatrixSize);

	gmresMaxwellBasis = new LargeVectorBasis(20, xnumberAdded, ynumberAdded, znumberAdded, 3);
	gmresCleanupBasis = new LargeVectorBasis(20, xnumberAdded, ynumberAdded, znumberAdded, 1);

	if (rank == 0) printf("creating fields arrays\n");
	if (rank == 0) fflush(stdout);
	// if(rank == 0) printLog("creating fields arrays\n");

	Efield = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	newEfield = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	tempEfield = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	explicitEfield = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	tempNodeParameter = create3array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	tempNodeVectorParameter = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	tempNodeMatrixParameter = create3matrixArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	rotB = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	Ederivative = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);

	Bfield = create3vectorArray(xnumberAdded, ynumberAdded, znumberAdded);
	newBfield = create3vectorArray(xnumberAdded, ynumberAdded, znumberAdded);
	rotE = create3vectorArray(xnumberAdded, ynumberAdded, znumberAdded);
	Bderivative = create3vectorArray(xnumberAdded, ynumberAdded, znumberAdded);
	tempCellParameter = create3array(xnumberAdded, ynumberAdded, znumberAdded);
	tempCellVectorParameter = create3vectorArray(xnumberAdded, ynumberAdded, znumberAdded);
	tempCellMatrixParameter = create3matrixArray(xnumberAdded, ynumberAdded, znumberAdded);

	massMatrix = new MassMatrix**[xnumberAdded + 1];
	tempMassMatrix = new MassMatrix**[xnumberAdded + 1];

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		massMatrix[i] = new MassMatrix*[ynumberAdded + 1];
		tempMassMatrix[i] = new MassMatrix*[ynumberAdded + 1];
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			massMatrix[i][j] = new MassMatrix[znumberAdded + 1];
			tempMassMatrix[i][j] = new MassMatrix[znumberAdded + 1];
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int curI = 0; curI < 3; ++curI) {
								for (int curJ = 0; curJ < 3; ++curJ) {
									massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[curI][curJ] = 0;
									tempMassMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[curI][curJ] = 0;
								}
							}
							massMatrix[i][j][k].xindex[tempI] = i + tempI - splineOrder - 1;
							massMatrix[i][j][k].yindex[tempJ] = j + tempJ - splineOrder - 1;
							massMatrix[i][j][k].zindex[tempK] = k + tempK - splineOrder - 1;

							tempMassMatrix[i][j][k].xindex[tempI] = i + tempI - splineOrder - 1;
							tempMassMatrix[i][j][k].yindex[tempJ] = j + tempJ - splineOrder - 1;
							tempMassMatrix[i][j][k].zindex[tempK] = k + tempK - splineOrder - 1;
						}
					}
				}
			}
		}
	}

	leftElevel = new Vector3d*[ynumberAdded + 1];
	rightElevel = new Vector3d*[ynumberAdded + 1];
	for (int j = 0; j < ynumberAdded + 1; ++j) {
		leftElevel[j] = new Vector3d[znumberAdded + 1];
		rightElevel[j] = new Vector3d[znumberAdded + 1];
		for (int k = 0; k < znumberAdded + 1; ++k) {
			leftElevel[j][k] = Vector3d(0, 0, 0);
			rightElevel[j][k] = Vector3d(0, 0, 0);
		}
	}

	if (rank == 0) printf("creating maxwellequation matrix arrays\n");
	fflush(stdout);
	//if(rank == 0) printLog("creating maxwell equation matrix arrays\n");

	maxwellEquationMatrix = new std::vector < MatrixElement > ***[xnumberAdded];
	maxwellEquationRightPart = new double ***[xnumberAdded];
	for (int i = 0; i < xnumberAdded; ++i) {
		maxwellEquationMatrix[i] = new std::vector < MatrixElement > **[ynumberAdded];
		maxwellEquationRightPart[i] = new double **[ynumberAdded];
		for (int j = 0; j < ynumberAdded; ++j) {
			maxwellEquationMatrix[i][j] = new std::vector < MatrixElement > *[znumberAdded];
			maxwellEquationRightPart[i][j] = new double *[znumberAdded];
			for (int k = 0; k < znumberAdded; ++k) {
				//if(rank == 0) printf("%d %d %d\n", i, j, k);
				//if(rank == 0) fflush((stdout));
				maxwellEquationMatrix[i][j][k] = new std::vector < MatrixElement >[maxwellEquationMatrixSize];
				maxwellEquationRightPart[i][j][k] = new double[maxwellEquationMatrixSize];
			}
		}
	}

	if (rank == 0) printf("creating arrays for divergence\n");
	if (rank == 0) fflush(stdout);
	//if(rank == 0) printLog("creating arrays for divergence\n");

	divergenceCleanUpMatrix = new std::vector < MatrixElement > ***[xnumberAdded + 1];

	divergenceCleanUpRightPart = create4array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1, 3);
	divergenceCleaningField = create4array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1, 3);
	divergenceCleaningPotential = create4array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1, 1);
	tempDivergenceCleaningPotential = create4array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1, 1);
	divergenceCleaningPotentialFourier = create3array(xnumberAdded, ynumberAdded, znumberAdded);

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		divergenceCleanUpMatrix[i] = new std::vector < MatrixElement > **[ynumberAdded + 1];
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			divergenceCleanUpMatrix[i][j] = new std::vector < MatrixElement > *[znumberAdded + 1];
			for (int k = 0; k < znumberAdded + 1; ++k) {
				divergenceCleanUpMatrix[i][j][k] = new std::vector < MatrixElement >[3];
			}
		}
	}

	fourierInput = create3complexArray(xnumberAdded, ynumberAdded, znumberAdded);
	fourierImage = create3complexArray(xnumberAdded, ynumberAdded, znumberAdded);
	fourierOutput = create3complexArray(xnumberAdded, ynumberAdded, znumberAdded);

	fourierScalarInput = create3complexArray(xnumberAdded, ynumberAdded, znumberAdded);
	fourierScalarOutput = create3complexArray(xnumberAdded, ynumberAdded, znumberAdded);
	fourierScalarTempOutput = create3complexArray(xnumberAdded, ynumberAdded, znumberAdded);
	fourierScalarTempOutput1 = create3complexArray(xnumberAdded, ynumberAdded, znumberAdded);

	fourierScalarMirrorInput = create3complexArray(2*xnumberAdded, ynumberAdded, znumberAdded);
	fourierScalarMirrorOutput = create3complexArray(2*xnumberAdded, ynumberAdded, znumberAdded);
	fourierScalarMirrorTempOutput = create3complexArray(2*xnumberAdded, ynumberAdded, znumberAdded);
	fourierScalarMirrorTempOutput1 = create3complexArray(2*xnumberAdded, ynumberAdded, znumberAdded);

	localFactorX = new Complex[2 * (xnumberAdded + 1)];
	localFactorY = new Complex[ynumberAdded + 1];
	localFactorZ = new Complex[znumberAdded + 1];

	rightMeanElevel = new Vector3d[xnumberAdded + 1];

	if (rank == 0) printf("creating arrays for particlesInBins\n");
	if (rank == 0) fflush(stdout);

	if (rank == 0) printf("creating arrays for parameters\n");
	fflush(stdout);
	//if(rank == 0) printLog("creating arrays for parameters\n");

	mostAcceleratedParticlesNumbers = new std::list < std::pair < int, double > >[typesNumber];

	particleConcentrations = new double ***[typesNumber];
	particleEnergies = new double ***[typesNumber];
	particleBulkVelocities = new Vector3d ***[typesNumber];
	for (int t = 0; t < typesNumber; ++t) {
		particleConcentrations[t] = create3array(xnumberAdded, ynumberAdded, znumberAdded);
		particleEnergies[t] = create3array(xnumberAdded, ynumberAdded, znumberAdded);
		particleBulkVelocities[t] = create3vectorArray(xnumberAdded, ynumberAdded, znumberAdded);
	}

	chargeDensity = create3array(xnumberAdded, ynumberAdded, znumberAdded);
	chargeDensityMinus = create3array(xnumberAdded, ynumberAdded, znumberAdded);
	chargeDensityHat = create3array(xnumberAdded, ynumberAdded, znumberAdded);
	pressureTensor = create3matrixArray(xnumberAdded, ynumberAdded, znumberAdded);


	if (rank == 0) printf("creating arrays for fluxes\n");
	fflush(stdout);
	//if(rank == 0) printLog("creating arrays for fluxes\n");

	electricFlux = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	electricFluxMinus = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	dielectricTensor = create3matrixArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	externalElectricFlux = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	divPressureTensor = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
	divPressureTensorMinus = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);

	if (rank == 0) printf("creating arrays for buffers\n");
	fflush(stdout);
	//if(rank == 0) printLog("creating arrays for buffers\n");

	rightOutNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	leftOutNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	leftInNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	rightInNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (2 + additionalBinNumber)];

	rightOutVectorNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 3 * (2 + additionalBinNumber)];
	leftOutVectorNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 3 * (2 + additionalBinNumber)];
	leftInVectorNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 3 * (2 + additionalBinNumber)];
	rightInVectorNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 3 * (2 + additionalBinNumber)];

	rightOutVectorCellBuffer = new double[(ynumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	leftOutVectorCellBuffer = new double[(ynumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	leftInVectorCellBuffer = new double[(ynumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	rightInVectorCellBuffer = new double[(ynumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];

	rightOutMatrixNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 9 * (2 + additionalBinNumber)];
	leftOutMatrixNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 9 * (2 + additionalBinNumber)];
	leftInMatrixNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 9 * (2 + additionalBinNumber)];
	rightInMatrixNodeBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 9 * (2 + additionalBinNumber)];

	rightOutCellBuffer = new double[(ynumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	leftOutCellBuffer = new double[(ynumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	leftInCellBuffer = new double[(ynumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	rightInCellBuffer = new double[(ynumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];

	backOutVectorNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 3 * (2 + additionalBinNumber)];
	frontOutVectorNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 3 * (2 + additionalBinNumber)];
	backInVectorNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 3 * (2 + additionalBinNumber)];
	frontInVectorNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 3 * (2 + additionalBinNumber)];

	backOutNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	frontOutNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	backInNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	frontInNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (2 + additionalBinNumber)];

	backOutMatrixNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 9 * (2 + additionalBinNumber)];
	frontOutMatrixNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 9 * (2 + additionalBinNumber)];
	backInMatrixNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 9 * (2 + additionalBinNumber)];
	frontInMatrixNodeBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 9 * (2 + additionalBinNumber)];

	backOutVectorCellBuffer = new double[(xnumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	frontOutVectorCellBuffer = new double[(xnumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	backInVectorCellBuffer = new double[(xnumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	frontInVectorCellBuffer = new double[(xnumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];

	backOutCellBuffer = new double[(xnumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	frontOutCellBuffer = new double[(xnumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	backInCellBuffer = new double[(xnumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	frontInCellBuffer = new double[(xnumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];

	topOutNodeBuffer = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (2 + additionalBinNumber)];
	bottomOutNodeBuffer = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (2 + additionalBinNumber)];
	topInNodeBuffer = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (2 + additionalBinNumber)];
	bottomInNodeBuffer = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (2 + additionalBinNumber)];

	topOutVectorNodeBuffer = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 3 * (2 + additionalBinNumber)];
	bottomOutVectorNodeBuffer = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 3 * (2 + additionalBinNumber)];
	topInVectorNodeBuffer = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 3 * (2 + additionalBinNumber)];
	bottomInVectorNodeBuffer = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 3 * (2 + additionalBinNumber)];

	topOutMatrixNodeBuffer = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 9 * (2 + additionalBinNumber)];
	bottomOutMatrixNodeBuffer = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 9 * (2 + additionalBinNumber)];
	topInMatrixNodeBuffer = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 9 * (2 + additionalBinNumber)];
	bottomInMatrixNodeBuffer = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 9 * (2 + additionalBinNumber)];

	topOutVectorCellBuffer = new double[(ynumberAdded) * (xnumberAdded) * 3 * (1 + additionalBinNumber)];
	bottomOutVectorCellBuffer = new double[(ynumberAdded) * (xnumberAdded) * 3 * (1 + additionalBinNumber)];
	topInVectorCellBuffer = new double[(ynumberAdded) * (xnumberAdded) * 3 * (1 + additionalBinNumber)];
	bottomInVectorCellBuffer = new double[(ynumberAdded) * (xnumberAdded) * 3 * (1 + additionalBinNumber)];

	topOutCellBuffer = new double[(ynumberAdded) * (xnumberAdded) * (1 + additionalBinNumber)];
	bottomOutCellBuffer = new double[(ynumberAdded) * (xnumberAdded) * (1 + additionalBinNumber)];
	topInCellBuffer = new double[(ynumberAdded) * (xnumberAdded) * (1 + additionalBinNumber)];
	bottomInCellBuffer = new double[(ynumberAdded) * (xnumberAdded) * (1 + additionalBinNumber)];

	////// Masha's buffers
	rightOutNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (3 + 2*additionalBinNumber)];
	leftOutNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (3 + 2*additionalBinNumber)];
	leftInNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (3 + 2*additionalBinNumber)];
	rightInNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (3 + 2*additionalBinNumber)];

	rightOutVectorNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];
	leftOutVectorNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];
	leftInVectorNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];
	rightInVectorNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];

	rightOutMatrixNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];
	leftOutMatrixNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];
	leftInMatrixNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];
	rightInMatrixNodeBufferMasha = new double[(ynumberAdded + 1) * (znumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];

	rightOutCellBufferMasha = new double[(ynumberAdded) * (znumberAdded) * (2 + 2*additionalBinNumber)];
	leftOutCellBufferMasha = new double[(ynumberAdded) * (znumberAdded) * (2 + 2*additionalBinNumber)];
	leftInCellBufferMasha = new double[(ynumberAdded) * (znumberAdded) * (2 + 2*additionalBinNumber)];
	rightInCellBufferMasha = new double[(ynumberAdded) * (znumberAdded) * (2 + 2*additionalBinNumber)];

	rightOutVectorCellBufferMasha = new double[(ynumberAdded) * (znumberAdded) * 3 * (2 + 2*additionalBinNumber)];
	leftOutVectorCellBufferMasha = new double[(ynumberAdded) * (znumberAdded) * 3 * (2 + 2*additionalBinNumber)];
	leftInVectorCellBufferMasha = new double[(ynumberAdded) * (znumberAdded) * 3 * (2 + 2*additionalBinNumber)];
	rightInVectorCellBufferMasha = new double[(ynumberAdded) * (znumberAdded) * 3 * (2 + 2*additionalBinNumber)];

	rightOutMatrixCellBufferMasha = new double[ynumberAdded * znumberAdded * 9 * (2 + 2*additionalBinNumber)];
	leftOutMatrixCellBufferMasha = new double[ynumberAdded * znumberAdded * 9 * (2 + 2*additionalBinNumber)];
	leftInMatrixCellBufferMasha = new double[ynumberAdded * znumberAdded * 9 * (2 + 2*additionalBinNumber)];
	rightInMatrixCellBufferMasha = new double[ynumberAdded * znumberAdded * 9 * (2 + 2*additionalBinNumber)];

	backOutVectorNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];
	frontOutVectorNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];
	backInVectorNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];
	frontInVectorNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];

	backOutNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (3 + 2*additionalBinNumber)];
	frontOutNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (3 + 2*additionalBinNumber)];
	backInNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (3 + 2*additionalBinNumber)];
	frontInNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (3 + 2*additionalBinNumber)];

	backOutMatrixNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];
	frontOutMatrixNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];
	backInMatrixNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];
	frontInMatrixNodeBufferMasha = new double[(xnumberAdded + 1) * (znumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];

	backOutVectorCellBufferMasha = new double[(xnumberAdded) * (znumberAdded) * 3 * (2 + 2*additionalBinNumber)];
	frontOutVectorCellBufferMasha = new double[(xnumberAdded) * (znumberAdded) * 3 * (2 + 2*additionalBinNumber)];
	backInVectorCellBufferMasha = new double[(xnumberAdded) * (znumberAdded) * 3 * (2 + 2*additionalBinNumber)];
	frontInVectorCellBufferMasha = new double[(xnumberAdded) * (znumberAdded) * 3 * (2 + 2*additionalBinNumber)];

	backOutCellBufferMasha = new double[(xnumberAdded) * (znumberAdded) * (2 + 2*additionalBinNumber)];
	frontOutCellBufferMasha = new double[(xnumberAdded) * (znumberAdded) * (2 + 2*additionalBinNumber)];
	backInCellBufferMasha = new double[(xnumberAdded) * (znumberAdded) * (2 + 2*additionalBinNumber)];
	frontInCellBufferMasha = new double[(xnumberAdded) * (znumberAdded) * (2 + 2*additionalBinNumber)];

	backOutMatrixCellBufferMasha = new double[xnumberAdded * znumberAdded * 9 * (2 + 2*additionalBinNumber)];
	frontOutMatrixCellBufferMasha = new double[xnumberAdded * znumberAdded * 9 * (2 + 2*additionalBinNumber)];
	backInMatrixCellBufferMasha = new double[xnumberAdded * znumberAdded * 9 * (2 + 2*additionalBinNumber)];
	frontInMatrixCellBufferMasha = new double[xnumberAdded * znumberAdded * 9 * (2 + 2*additionalBinNumber)];

	topOutNodeBufferMasha = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (3 + 2*additionalBinNumber)];
	bottomOutNodeBufferMasha = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (3 + 2*additionalBinNumber)];
	topInNodeBufferMasha = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (3 + 2*additionalBinNumber)];
	bottomInNodeBufferMasha = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (3 + 2*additionalBinNumber)];

	topOutVectorNodeBufferMasha = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];
	bottomOutVectorNodeBufferMasha = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];
	topInVectorNodeBufferMasha = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];
	bottomInVectorNodeBufferMasha = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 3 * (3 + 2*additionalBinNumber)];

	topOutMatrixNodeBufferMasha = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];
	bottomOutMatrixNodeBufferMasha = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];
	topInMatrixNodeBufferMasha = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];
	bottomInMatrixNodeBufferMasha = new double[(xnumberAdded + 1) * (ynumberAdded + 1) * 9 * (3 + 2*additionalBinNumber)];

	topOutVectorCellBufferMasha = new double[(ynumberAdded) * (xnumberAdded) * 3 * (2 + 2*additionalBinNumber)];
	bottomOutVectorCellBufferMasha = new double[(ynumberAdded) * (xnumberAdded) * 3 * (2 + 2*additionalBinNumber)];
	topInVectorCellBufferMasha = new double[(ynumberAdded) * (xnumberAdded) * 3 * (2 + 2*additionalBinNumber)];
	bottomInVectorCellBufferMasha = new double[(ynumberAdded) * (xnumberAdded) * 3 * (2 + 2*additionalBinNumber)];

	topOutCellBufferMasha = new double[(ynumberAdded) * (xnumberAdded) * (2 + 2*additionalBinNumber)];
	bottomOutCellBufferMasha = new double[(ynumberAdded) * (xnumberAdded) * (2 + 2*additionalBinNumber)];
	topInCellBufferMasha = new double[(ynumberAdded) * (xnumberAdded) * (2 + 2*additionalBinNumber)];
	bottomInCellBufferMasha = new double[(ynumberAdded) * (xnumberAdded) * (2 + 2*additionalBinNumber)];

	topOutMatrixCellBufferMasha = new double[xnumberAdded  * ynumberAdded * 9 * (2 + 2*additionalBinNumber)];
	bottomOutMatrixCellBufferMasha = new double[xnumberAdded * ynumberAdded * 9 * (2 + 2*additionalBinNumber)];
	topInMatrixCellBufferMasha = new double[xnumberAdded * ynumberAdded * 9 * (2 + 2*additionalBinNumber)];
	bottomInMatrixCellBufferMasha = new double[xnumberAdded * ynumberAdded * 9 * (2 + 2*additionalBinNumber)];

	////buneman E
	leftOutBunemanExBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	rightOutBunemanExBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	leftInBunemanExBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	rightInBunemanExBuffer = new double[(ynumberAdded + 1) * (znumberAdded + 1) * (1 + additionalBinNumber)];

	frontOutBunemanExBuffer = new double[(xnumberAdded) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	backOutBunemanExBuffer = new double[(xnumberAdded) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	frontInBunemanExBuffer = new double[(xnumberAdded) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	backInBunemanExBuffer = new double[(xnumberAdded) * (znumberAdded + 1) * (2 + additionalBinNumber)];

	bottomOutBunemanExBuffer = new double[(ynumberAdded + 1) * (xnumberAdded) * (2 + additionalBinNumber)];
	topOutBunemanExBuffer = new double[(ynumberAdded + 1) * (xnumberAdded) * (2 + additionalBinNumber)];
	bottomInBunemanExBuffer = new double[(ynumberAdded + 1) * (xnumberAdded) * (2 + additionalBinNumber)];
	topInBunemanExBuffer = new double[(ynumberAdded + 1) * (xnumberAdded) * (2 + additionalBinNumber)];

	leftOutBunemanEyBuffer = new double[(ynumberAdded) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	rightOutBunemanEyBuffer = new double[(ynumberAdded) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	leftInBunemanEyBuffer = new double[(ynumberAdded) * (znumberAdded + 1) * (2 + additionalBinNumber)];
	rightInBunemanEyBuffer = new double[(ynumberAdded) * (znumberAdded + 1) * (2 + additionalBinNumber)];

	frontOutBunemanEyBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	backOutBunemanEyBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	frontInBunemanEyBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	backInBunemanEyBuffer = new double[(xnumberAdded + 1) * (znumberAdded + 1) * (1 + additionalBinNumber)];

	bottomOutBunemanEyBuffer = new double[(ynumberAdded) * (xnumberAdded + 1) * (2 + additionalBinNumber)];
	topOutBunemanEyBuffer = new double[(ynumberAdded) * (xnumberAdded + 1) * (2 + additionalBinNumber)];
	bottomInBunemanEyBuffer = new double[(ynumberAdded) * (xnumberAdded + 1) * (2 + additionalBinNumber)];
	topInBunemanEyBuffer = new double[(ynumberAdded) * (xnumberAdded + 1) * (2 + additionalBinNumber)];

	leftOutBunemanEzBuffer = new double[(ynumberAdded + 1) * (znumberAdded) * (2 + additionalBinNumber)];
	rightOutBunemanEzBuffer = new double[(ynumberAdded + 1) * (znumberAdded) * (2 + additionalBinNumber)];
	leftInBunemanEzBuffer = new double[(ynumberAdded + 1) * (znumberAdded) * (2 + additionalBinNumber)];
	rightInBunemanEzBuffer = new double[(ynumberAdded + 1) * (znumberAdded) * (2 + additionalBinNumber)];

	frontOutBunemanEzBuffer = new double[(xnumberAdded + 1) * (znumberAdded) * (2 + additionalBinNumber)];
	backOutBunemanEzBuffer = new double[(xnumberAdded + 1) * (znumberAdded) * (2 + additionalBinNumber)];
	frontInBunemanEzBuffer = new double[(xnumberAdded + 1) * (znumberAdded) * (2 + additionalBinNumber)];
	backInBunemanEzBuffer = new double[(xnumberAdded + 1) * (znumberAdded) * (2 + additionalBinNumber)];

	bottomOutBunemanEzBuffer = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (1 + additionalBinNumber)];
	topOutBunemanEzBuffer = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (1 + additionalBinNumber)];
	bottomInBunemanEzBuffer = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (1 + additionalBinNumber)];
	topInBunemanEzBuffer = new double[(ynumberAdded + 1) * (xnumberAdded + 1) * (1 + additionalBinNumber)];

	///buneman B
	leftOutBunemanBxBuffer = new double[(ynumberAdded) * (znumberAdded) * (2 + additionalBinNumber)];
	rightOutBunemanBxBuffer = new double[(ynumberAdded) * (znumberAdded) * (2 + additionalBinNumber)];
	leftInBunemanBxBuffer = new double[(ynumberAdded) * (znumberAdded) * (2 + additionalBinNumber)];
	rightInBunemanBxBuffer = new double[(ynumberAdded) * (znumberAdded) * (2 + additionalBinNumber)];

	frontOutBunemanBxBuffer = new double[(xnumberAdded + 1) * (znumberAdded) * (1 + additionalBinNumber)];
	backOutBunemanBxBuffer = new double[(xnumberAdded + 1) * (znumberAdded) * (1 + additionalBinNumber)];
	frontInBunemanBxBuffer = new double[(xnumberAdded + 1) * (znumberAdded) * (1 + additionalBinNumber)];
	backInBunemanBxBuffer = new double[(xnumberAdded + 1) * (znumberAdded) * (1 + additionalBinNumber)];

	bottomOutBunemanBxBuffer = new double[(ynumberAdded) * (xnumberAdded + 1) * (1 + additionalBinNumber)];
	topOutBunemanBxBuffer = new double[(ynumberAdded) * (xnumberAdded + 1) * (1 + additionalBinNumber)];
	bottomInBunemanBxBuffer = new double[(ynumberAdded) * (xnumberAdded + 1) * (1 + additionalBinNumber)];
	topInBunemanBxBuffer = new double[(ynumberAdded) * (xnumberAdded + 1) * (1 + additionalBinNumber)];

	leftOutBunemanByBuffer = new double[(ynumberAdded + 1) * (znumberAdded) * (1 + additionalBinNumber)];
	rightOutBunemanByBuffer = new double[(ynumberAdded + 1) * (znumberAdded) * (1 + additionalBinNumber)];
	leftInBunemanByBuffer = new double[(ynumberAdded + 1) * (znumberAdded) * (1 + additionalBinNumber)];
	rightInBunemanByBuffer = new double[(ynumberAdded + 1) * (znumberAdded) * (1 + additionalBinNumber)];

	frontOutBunemanByBuffer = new double[(xnumberAdded) * (znumberAdded) * (2 + additionalBinNumber)];
	backOutBunemanByBuffer = new double[(xnumberAdded) * (znumberAdded) * (2 + additionalBinNumber)];
	frontInBunemanByBuffer = new double[(xnumberAdded) * (znumberAdded) * (2 + additionalBinNumber)];
	backInBunemanByBuffer = new double[(xnumberAdded) * (znumberAdded) * (2 + additionalBinNumber)];

	bottomOutBunemanByBuffer = new double[(ynumberAdded + 1) * (xnumberAdded) * (1 + additionalBinNumber)];
	topOutBunemanByBuffer = new double[(ynumberAdded + 1) * (xnumberAdded) * (1 + additionalBinNumber)];
	bottomInBunemanByBuffer = new double[(ynumberAdded + 1) * (xnumberAdded) * (1 + additionalBinNumber)];
	topInBunemanByBuffer = new double[(ynumberAdded + 1) * (xnumberAdded) * (1 + additionalBinNumber)];

	leftOutBunemanBzBuffer = new double[(ynumberAdded) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	rightOutBunemanBzBuffer = new double[(ynumberAdded) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	leftInBunemanBzBuffer = new double[(ynumberAdded) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	rightInBunemanBzBuffer = new double[(ynumberAdded) * (znumberAdded + 1) * (1 + additionalBinNumber)];

	frontOutBunemanBzBuffer = new double[(xnumberAdded) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	backOutBunemanBzBuffer = new double[(xnumberAdded) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	frontInBunemanBzBuffer = new double[(xnumberAdded) * (znumberAdded + 1) * (1 + additionalBinNumber)];
	backInBunemanBzBuffer = new double[(xnumberAdded) * (znumberAdded + 1) * (1 + additionalBinNumber)];

	bottomOutBunemanBzBuffer = new double[(ynumberAdded) * (xnumberAdded) * (2 + additionalBinNumber)];
	topOutBunemanBzBuffer = new double[(ynumberAdded) * (xnumberAdded) * (2 + additionalBinNumber)];
	bottomInBunemanBzBuffer = new double[(ynumberAdded) * (xnumberAdded) * (2 + additionalBinNumber)];
	topInBunemanBzBuffer = new double[(ynumberAdded) * (xnumberAdded) * (2 + additionalBinNumber)];


	rightOutGmresBuffer = new double[(ynumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	leftOutGmresBuffer = new double[(ynumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	leftInGmresBuffer = new double[(ynumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	rightInGmresBuffer = new double[(ynumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];

	frontOutGmresBuffer = new double[(xnumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	backOutGmresBuffer = new double[(xnumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	frontInGmresBuffer = new double[(xnumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];
	backInGmresBuffer = new double[(xnumberAdded) * (znumberAdded) * 3 * (1 + additionalBinNumber)];

	topOutGmresBuffer = new double[(ynumberAdded) * (xnumberAdded) * 3 * (1 + additionalBinNumber)];
	topInGmresBuffer = new double[(ynumberAdded) * (xnumberAdded) * 3 * (1 + additionalBinNumber)];
	bottomOutGmresBuffer = new double[(ynumberAdded) * (xnumberAdded) * 3 * (1 + additionalBinNumber)];
	bottomInGmresBuffer = new double[(ynumberAdded) * (xnumberAdded) * 3 * (1 + additionalBinNumber)];

	rightOutDivergenceBuffer = new double[(ynumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	leftOutDivergenceBuffer = new double[(ynumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	leftInDivergenceBuffer = new double[(ynumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	rightInDivergenceBuffer = new double[(ynumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];

	frontOutDivergenceBuffer = new double[(xnumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	backOutDivergenceBuffer = new double[(xnumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	frontInDivergenceBuffer = new double[(xnumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];
	backInDivergenceBuffer = new double[(xnumberAdded) * (znumberAdded) * (1 + additionalBinNumber)];

	topOutDivergenceBuffer = new double[(ynumberAdded) * (xnumberAdded) * (1 + additionalBinNumber)];
	topInDivergenceBuffer = new double[(ynumberAdded) * (xnumberAdded) * (1 + additionalBinNumber)];
	bottomOutDivergenceBuffer = new double[(ynumberAdded) * (xnumberAdded) * (1 + additionalBinNumber)];
	bottomInDivergenceBuffer = new double[(ynumberAdded) * (xnumberAdded) * (1 + additionalBinNumber)];

	tempCellParameterLeft = create3array(2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);
	tempCellParameterRight = create3array(2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);
	tempNodeParameterLeft = create3array(3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);
	tempNodeParameterRight = create3array(3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);

	tempCellVectorParameterLeft = create3vectorArray(2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);
	tempCellVectorParameterRight = create3vectorArray(2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);
	tempNodeVectorParameterLeft = create3vectorArray(3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);
	tempNodeVectorParameterRight = create3vectorArray(3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);

	tempCellMatrixParameterLeft = create3matrixArray(2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);
	tempCellMatrixParameterRight = create3matrixArray(2 + 2*additionalBinNumber, ynumberAdded, znumberAdded);
	tempNodeMatrixParameterLeft = create3matrixArray(3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);
	tempNodeMatrixParameterRight = create3matrixArray(3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);

	tempNodeMassMatrixParameterLeft = new MassMatrix **[3 + 2 * additionalBinNumber];
	tempNodeMassMatrixParameterRight = new MassMatrix **[3 + 2 * additionalBinNumber];

	for (int i = 0; i < 3 + 2 * additionalBinNumber; ++i) {
		tempNodeMassMatrixParameterLeft[i] = new MassMatrix *[ynumberAdded + 1];
		tempNodeMassMatrixParameterRight[i] = new MassMatrix *[ynumberAdded + 1];
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			tempNodeMassMatrixParameterLeft[i][j] = new MassMatrix[znumberAdded + 1];
			tempNodeMassMatrixParameterRight[i][j] = new MassMatrix[znumberAdded + 1];
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							tempNodeMassMatrixParameterLeft[i][j][k].matrix[tempI][tempJ][tempK] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
							tempNodeMassMatrixParameterRight[i][j][k].matrix[tempI][tempJ][tempK] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
						}
					}
				}
			}
		}
	}

	tempCellParameterFront = create3array(xnumberAdded, 2 + 2*additionalBinNumber, znumberAdded);
	tempCellParameterBack = create3array(xnumberAdded, 2 + 2*additionalBinNumber, znumberAdded);
	tempNodeParameterFront = create3array(xnumberAdded + 1, 3 + 2*additionalBinNumber, znumberAdded + 1);
	tempNodeParameterBack = create3array(xnumberAdded + 1, 3 + 2*additionalBinNumber, znumberAdded + 1);

	tempCellVectorParameterFront = create3vectorArray(xnumberAdded, 2 + 2*additionalBinNumber, znumberAdded);
	tempCellVectorParameterBack = create3vectorArray(xnumberAdded, 2 + 2*additionalBinNumber, znumberAdded);
	tempNodeVectorParameterFront = create3vectorArray(xnumberAdded + 1, 3 + 2*additionalBinNumber, znumberAdded + 1);
	tempNodeVectorParameterBack = create3vectorArray(xnumberAdded + 1, 3 + 2*additionalBinNumber, znumberAdded + 1);

	tempCellMatrixParameterFront = create3matrixArray(xnumberAdded, 2 + 2*additionalBinNumber, znumberAdded);
	tempCellMatrixParameterBack = create3matrixArray(xnumberAdded, 2 + 2*additionalBinNumber, znumberAdded);
	tempNodeMatrixParameterFront = create3matrixArray(xnumberAdded + 1, 3 + 2*additionalBinNumber, znumberAdded + 1);
	tempNodeMatrixParameterBack = create3matrixArray(xnumberAdded + 1, 3 + 2*additionalBinNumber, znumberAdded + 1);

	tempNodeMassMatrixParameterFront = new MassMatrix **[xnumberAdded + 1];
	tempNodeMassMatrixParameterBack = new MassMatrix **[xnumberAdded + 1];

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		tempNodeMassMatrixParameterFront[i] = new MassMatrix *[3 + 2 * additionalBinNumber];
		tempNodeMassMatrixParameterBack[i] = new MassMatrix *[3 + 2 * additionalBinNumber];
		for (int j = 0; j < 3 + 2 * additionalBinNumber; ++j) {
			tempNodeMassMatrixParameterFront[i][j] = new MassMatrix[znumberAdded + 1];
			tempNodeMassMatrixParameterBack[i][j] = new MassMatrix[znumberAdded + 1];
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							tempNodeMassMatrixParameterFront[i][j][k].matrix[tempI][tempJ][tempK] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
							tempNodeMassMatrixParameterBack[i][j][k].matrix[tempI][tempJ][tempK] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
						}
					}
				}
			}
		}
	}

	tempCellParameterBottom = create3array(xnumberAdded, ynumberAdded, 2 + 2*additionalBinNumber);
	tempCellParameterTop = create3array(xnumberAdded, ynumberAdded, 2 + 2*additionalBinNumber);
	tempNodeParameterBottom = create3array(xnumberAdded + 1, ynumberAdded + 1, 3 + 2*additionalBinNumber);
	tempNodeParameterTop = create3array(xnumberAdded + 1, ynumberAdded + 1, 3 + 2*additionalBinNumber);

	tempCellVectorParameterBottom = create3vectorArray(xnumberAdded, ynumberAdded, 2 + 2*additionalBinNumber);
	tempCellVectorParameterTop = create3vectorArray(xnumberAdded, ynumberAdded, 2 + 2*additionalBinNumber);
	tempNodeVectorParameterBottom = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, 3 + 2*additionalBinNumber);
	tempNodeVectorParameterTop = create3vectorArray(xnumberAdded + 1, ynumberAdded + 1, 3 + 2*additionalBinNumber);

	tempCellMatrixParameterBottom = create3matrixArray(xnumberAdded, ynumberAdded, 2 + 2*additionalBinNumber);
	tempCellMatrixParameterTop = create3matrixArray(xnumberAdded, ynumberAdded, 2 + 2*additionalBinNumber);
	tempNodeMatrixParameterBottom = create3matrixArray(xnumberAdded + 1, ynumberAdded + 1, 3 + 2*additionalBinNumber);
	tempNodeMatrixParameterTop = create3matrixArray(xnumberAdded + 1, ynumberAdded + 1, 3 + 2*additionalBinNumber);

	tempNodeMassMatrixParameterBottom = new MassMatrix **[xnumberAdded + 1];
	tempNodeMassMatrixParameterTop = new MassMatrix **[xnumberAdded + 1];

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		tempNodeMassMatrixParameterBottom[i] = new MassMatrix *[ynumberAdded + 1];
		tempNodeMassMatrixParameterTop[i] = new MassMatrix *[ynumberAdded + 1];
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			tempNodeMassMatrixParameterBottom[i][j] = new MassMatrix[3 + 2 * additionalBinNumber];
			tempNodeMassMatrixParameterTop[i][j] = new MassMatrix[3 + 2 * additionalBinNumber];
			for (int k = 0; k < 3 + 2 * additionalBinNumber; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							tempNodeMassMatrixParameterBottom[i][j][k].matrix[tempI][tempJ][tempK] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
							tempNodeMassMatrixParameterTop[i][j][k].matrix[tempI][tempJ][tempK] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
						}
					}
				}
			}
		}
	}

	residualBiconjugateDivE = create4array(xnumberAdded, ynumberAdded, znumberAdded, 1);
	firstResidualBiconjugateDivE = create4array(xnumberAdded, ynumberAdded, znumberAdded, 1);
	pBiconjugateDivE = create4array(xnumberAdded, ynumberAdded, znumberAdded, 1);
	vBiconjugateDivE = create4array(xnumberAdded, ynumberAdded, znumberAdded, 1);
	sBiconjugateDivE = create4array(xnumberAdded, ynumberAdded, znumberAdded, 1);
	tBiconjugateDivE = create4array(xnumberAdded, ynumberAdded, znumberAdded, 1);


	residualBiconjugateMaxwell = create4array(xnumberAdded, ynumberAdded, znumberAdded, 3);
	firstResidualBiconjugateMaxwell = create4array(xnumberAdded, ynumberAdded, znumberAdded, 3);
	pBiconjugateMaxwell = create4array(xnumberAdded, ynumberAdded, znumberAdded, 3);
	vBiconjugateMaxwell = create4array(xnumberAdded, ynumberAdded, znumberAdded, 3);
	sBiconjugateMaxwell = create4array(xnumberAdded, ynumberAdded, znumberAdded, 3);
	tBiconjugateMaxwell = create4array(xnumberAdded, ynumberAdded, znumberAdded, 3);

	///buneman
	if (solverType == BUNEMAN) {
		bunemanJx = create3array(xnumberAdded, ynumberAdded + 1, znumberAdded + 1);
		bunemanEx = create3array(xnumberAdded, ynumberAdded + 1, znumberAdded + 1);
		bunemanNewEx = create3array(xnumberAdded, ynumberAdded + 1, znumberAdded + 1);
		tempBunemanExParameter = create3array(xnumberAdded, ynumberAdded + 1, znumberAdded + 1);
		bunemanDivCleaningEx = create3array(xnumberAdded, ynumberAdded + 1, znumberAdded + 1);

		bunemanJy = create3array(xnumberAdded + 1, ynumberAdded, znumberAdded + 1);
		bunemanEy = create3array(xnumberAdded + 1, ynumberAdded, znumberAdded + 1);
		bunemanNewEy = create3array(xnumberAdded + 1, ynumberAdded, znumberAdded + 1);
		tempBunemanEyParameter = create3array(xnumberAdded + 1, ynumberAdded, znumberAdded + 1);
		bunemanDivCleaningEy = create3array(xnumberAdded + 1, ynumberAdded, znumberAdded + 1);

		bunemanJz = create3array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded);
		bunemanEz = create3array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded);
		bunemanNewEz = create3array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded);
		tempBunemanEzParameter = create3array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded);
		bunemanDivCleaningEz = create3array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded);

		bunemanBx = create3array(xnumberAdded + 1, ynumberAdded, znumberAdded);
		bunemanNewBx = create3array(xnumberAdded + 1, ynumberAdded, znumberAdded);
		tempBunemanBxParameter = create3array(xnumberAdded + 1, ynumberAdded, znumberAdded);
		bunemanDivCleaningBx = create3array(xnumberAdded + 1, ynumberAdded, znumberAdded);

		bunemanBy = create3array(xnumberAdded, ynumberAdded + 1, znumberAdded);
		bunemanNewBy = create3array(xnumberAdded, ynumberAdded + 1, znumberAdded);
		tempBunemanByParameter = create3array(xnumberAdded, ynumberAdded + 1, znumberAdded);
		bunemanDivCleaningBy = create3array(xnumberAdded, ynumberAdded + 1, znumberAdded);

		bunemanBz = create3array(xnumberAdded, ynumberAdded, znumberAdded + 1);
		bunemanNewBz = create3array(xnumberAdded, ynumberAdded, znumberAdded + 1);
		tempBunemanBzParameter = create3array(xnumberAdded, ynumberAdded, znumberAdded + 1);
		bunemanDivCleaningBz = create3array(xnumberAdded, ynumberAdded, znumberAdded + 1);
		

		bunemanChargeDensity = create3array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1);
		bunemanDivergenceCleaningPotential = create4array(xnumberAdded + 1, ynumberAdded + 1, znumberAdded + 1, 1);

		/// left right
		tempBunemanJxLeft = create3array(2 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);
		tempBunemanJxRight = create3array(2 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded + 1);
		
		tempBunemanJyLeft = create3array(3 + 2*additionalBinNumber, ynumberAdded, znumberAdded + 1);
		tempBunemanJyRight = create3array(3 + 2*additionalBinNumber, ynumberAdded, znumberAdded + 1);
		
		tempBunemanJzLeft = create3array(3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded);
		tempBunemanJzRight = create3array(3 + 2*additionalBinNumber, ynumberAdded + 1, znumberAdded);

		///front back
		tempBunemanJxFront = create3array(xnumberAdded, 3 + 2*additionalBinNumber, znumberAdded + 1);
		tempBunemanJxBack = create3array(xnumberAdded, 3 + 2*additionalBinNumber, znumberAdded + 1);

		tempBunemanJyFront = create3array(xnumberAdded + 1, 2 + 2*additionalBinNumber, znumberAdded + 1);
		tempBunemanJyBack = create3array(xnumberAdded + 1, 2 + 2*additionalBinNumber, znumberAdded + 1);

		tempBunemanJzFront = create3array(xnumberAdded + 1, 3 + 2*additionalBinNumber, znumberAdded);
		tempBunemanJzBack = create3array(xnumberAdded + 1, 3 + 2*additionalBinNumber, znumberAdded);
		// bottom top

		tempBunemanJxBottom = create3array(xnumberAdded, ynumberAdded + 1, 3 + 2*additionalBinNumber);
		tempBunemanJxTop = create3array(xnumberAdded, ynumberAdded + 1, 3 + 2*additionalBinNumber);

		tempBunemanJyBottom = create3array(xnumberAdded + 1, ynumberAdded, 3 + 2*additionalBinNumber);
		tempBunemanJyTop = create3array(xnumberAdded + 1, ynumberAdded, 3 + 2*additionalBinNumber);

		tempBunemanJzBottom = create3array(xnumberAdded + 1, ynumberAdded + 1, 2 + 2*additionalBinNumber);
		tempBunemanJzTop = create3array(xnumberAdded + 1, ynumberAdded + 1, 2 + 2*additionalBinNumber);
	}

	arrayCreated = true;

	MPI_Barrier(MPI_COMM_WORLD);
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
		} else {
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

		particlesTrajectoryFile = fopen((outputDir + "particlesTrajectories.dat").c_str(), "w");
		fclose(particlesTrajectoryFile);
		particlesTrajectoryFile = fopen((outputDir + "electronsTrajectories.dat").c_str(), "w");
		fclose(particlesTrajectoryFile);
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
		EfieldFile = fopen((outputDir + "EfieldXY.dat").c_str(), "w");
		fclose(EfieldFile);
		EfieldFile = fopen((outputDir + "EfieldYZ.dat").c_str(), "w");
		fclose(EfieldFile);
		EfieldFile = fopen((outputDir + "EfieldXZ.dat").c_str(), "w");
		fclose(EfieldFile);
		EfieldFile = fopen((outputDir + "EfieldX.dat").c_str(), "w");
		fclose(EfieldFile);
		EfieldFile = fopen((outputDir + "EfieldY.dat").c_str(), "w");
		fclose(EfieldFile);
		EfieldFile = fopen((outputDir + "EfieldZ.dat").c_str(), "w");
		fclose(EfieldFile);
		EfieldFile = fopen((reducedOutputDir + "EfieldReduced.dat").c_str(), "w");
		fclose(EfieldFile);
		BfieldFile = fopen((outputDir + "Bfield.dat").c_str(), "w");
		fclose(BfieldFile);
		BfieldFile = fopen((outputDir + "BfieldXY.dat").c_str(), "w");
		fclose(BfieldFile);
		BfieldFile = fopen((outputDir + "BfieldYZ.dat").c_str(), "w");
		fclose(BfieldFile);
		BfieldFile = fopen((outputDir + "BfieldXZ.dat").c_str(), "w");
		fclose(BfieldFile);
		BfieldFile = fopen((outputDir + "BfieldX.dat").c_str(), "w");
		fclose(BfieldFile);
		BfieldFile = fopen((outputDir + "BfieldY.dat").c_str(), "w");
		fclose(BfieldFile);
		BfieldFile = fopen((outputDir + "BfieldZ.dat").c_str(), "w");
		fclose(BfieldFile);
		BfieldFile = fopen((reducedOutputDir + "BfieldReduced.dat").c_str(), "w");
		fclose(BfieldFile);
		velocityFile = fopen((outputDir + "velocity.dat").c_str(), "w");
		fclose(velocityFile);
		velocityFile = fopen((outputDir + "velocityXY.dat").c_str(), "w");
		fclose(velocityFile);
		velocityFile = fopen((outputDir + "velocityXZ.dat").c_str(), "w");
		fclose(velocityFile);
		velocityFile = fopen((outputDir + "velocityYZ.dat").c_str(), "w");
		fclose(velocityFile);
		velocityFile = fopen((outputDir + "velocityX.dat").c_str(), "w");
		fclose(velocityFile);
		velocityFile = fopen((outputDir + "velocityY.dat").c_str(), "w");
		fclose(velocityFile);
		velocityFile = fopen((outputDir + "velocityZ.dat").c_str(), "w");
		fclose(velocityFile);
		Xfile = fopen((outputDir + "Xfile.dat").c_str(), "w");
		fclose(Xfile);
		Xfile = fopen((reducedOutputDir + "XfileReduced.dat").c_str(), "w");
		fclose(Xfile);
		Yfile = fopen((outputDir + "Yfile.dat").c_str(), "w");
		fclose(Yfile);
		Yfile = fopen((reducedOutputDir + "YfileReduced.dat").c_str(), "w");
		fclose(Yfile);
		Zfile = fopen((outputDir + "Zfile.dat").c_str(), "w");
		fclose(Zfile);
		Zfile = fopen((reducedOutputDir + "ZfileReduced.dat").c_str(), "w");
		fclose(Zfile);
		generalFile = fopen((outputDir + "general.dat").c_str(), "w");
		fclose(generalFile);
		generalAnisotropyFile = fopen((outputDir + "generalAnisotropy.dat").c_str(), "w");
		fclose(generalAnisotropyFile);
		densityFile = fopen((outputDir + "concentrations.dat").c_str(), "w");
		fclose(densityFile);
		densityFile = fopen((outputDir + "concentrationsXY.dat").c_str(), "w");
		fclose(densityFile);
		densityFile = fopen((outputDir + "concentrationsXZ.dat").c_str(), "w");
		fclose(densityFile);
		densityFile = fopen((outputDir + "concentrationsYZ.dat").c_str(), "w");
		fclose(densityFile);
		densityFile = fopen((outputDir + "concentrationsX.dat").c_str(), "w");
		fclose(densityFile);
		densityFile = fopen((outputDir + "concentrationsY.dat").c_str(), "w");
		fclose(densityFile);
		densityFile = fopen((outputDir + "concentrationsZ.dat").c_str(), "w");
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

		FILE* tempFile = fopen((outputDir + "rightPart.dat").c_str(), "w");
		fclose(tempFile);
		printf("finish creating files\n");
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
	} else if (omega > cyclothronOmegaProton / 100.0) {
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
	} else if (omega > 0.01) {
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
			debyeLength2 += (4 * pi * types[i].concentration * sqr(types[i].charge)) / (kBoltzman_normalized * types[i].
				temperatureX);
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
	} else if (debyeNumber < 100.0) {
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
			massElectron * kBoltzman_normalized * temperature);
		double movingMomentumElectron = massElectron * V0.norm();
		double momentumElectron = thermalMomentumElectron + movingMomentumElectron;
		double gyroRadiusElectron = momentumElectron * speed_of_light_normalized / (electron_charge_normalized * B0.norm());
		double thermalMomentumProton = sqrt(massProton * kBoltzman_normalized * temperature);
		double movingMomentumProton = massProton * V0.norm();
		double momentumProton = thermalMomentumProton + movingMomentumProton;
		double gyroRadiusProton = momentumProton * speed_of_light_normalized / (electron_charge_normalized * B0.norm());
		if (deltaX > 0.5 * gyroRadiusElectron) {
			if (rank == 0) printf("deltaX > 0.5*gyroRadiusElectron\n");
			fflush(stdout);
			if (rank == 0) fprintf(informationFile, "deltaX > 0.5*gyroRadiusElectron\n");
		}

		if (rank == 0) printf("deltaX/gyroRadiusElectron = %g\n", deltaX / gyroRadiusElectron);
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "deltaX/gyroRadiusElectron = %g\n", deltaX / gyroRadiusElectron);

		if (xsizeGeneral < 2 * gyroRadiusProton) {
			if (rank == 0) printf("xsize < 2*gyroRadiusProton\n");
			fflush(stdout);
			if (rank == 0) fprintf(informationFile, "xsize < 2*gyroRadiusProton\n");
		}

		if (rank == 0) printf("xsize/gyroRadiusProton= %g\n", xsizeGeneral / gyroRadiusProton);
		fflush(stdout);
		if (rank == 0) fprintf(informationFile, "xsize/gyroRadiusProton = %g\n", xsizeGeneral / gyroRadiusProton);
	}

	if (rank == 0) fclose(informationFile);
}

void Simulation::initializeBell() {
	const double protonCRfraction = 1.0 / 3.0;
	const double gammaCR = 2;
	double protonCRconcentration = types[1].concentration * protonCRfraction;
	int protonCRperbin = types[1].particlesPerBin * protonCRfraction;
	double electronConcentration = types[1].concentration + protonCRconcentration;
	types[0].concentration = electronConcentration;
	types[0].particlesPerBin += protonCRperbin;
	types[0].particesDeltaX = deltaX / types[0].particlesPerBin;
	types[0].particesDeltaY = deltaY / types[0].particlesPerBin;
	types[0].particesDeltaZ = deltaZ / types[0].particlesPerBin;
	double CRdeltaX = deltaX / protonCRperbin;
	double CRdeltaY = deltaY / protonCRperbin;
	double CRdeltaZ = deltaZ / protonCRperbin;
	Vector3d Vsh = V0;
	V0 = Vector3d(0, 0, 0);
	Vector3d Vd = Vsh * (protonCRconcentration / electronConcentration);
	//boundaryConditionTypeX = SUPER_CONDUCTOR_LEFT;
	boundaryConditionTypeX = PERIODIC;
	//createParticles();
	//E0 = E0 - Vsh.vectorMult(B0) / (speed_of_light_normalized);
	//initializeAlfvenWaveY(10, 1.0E-4);
	double alfvenV = sqrt(
		B0.scalarMult(B0) / (4 * pi * (massElectron * types[0].concentration + massProton * types[1].concentration)));
	double alfvenMach = Vsh.norm() / alfvenV;
	double dzeta = protonCRconcentration * speed_of_light_normalized * sqrt(gammaCR * gammaCR - 1) / (types[1].
		concentration * Vsh.norm());
	double rCR = massProton * speed_of_light_normalized_sqr * sqrt(gammaCR * gammaCR - 1) / (electron_charge_normalized *
		B0.norm());
	double kmax = dzeta * 0.5 * sqr(Vsh.norm() / alfvenV) / rCR;
	double lambda = 2 * pi / kmax;
	int wavesCount = xsizeGeneral / lambda;
	double kw = 2 * pi * wavesCount / xsizeGeneral;

	double E = B0.x * 0.1;
	E = 0;
	if (solverType == BUNEMAN) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEx[i][j][k] = E0.x;
					bunemanNewEx[i][j][k] = bunemanEx[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEy[i][j][k] = E0.y + E * cos(kw * xgrid[i]);
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i < 1 + additionalBinNumber) {
							bunemanEy[i][j][k] = 0;
						}
					}
					bunemanNewEy[i][j][k] = bunemanEy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanEz[i][j][k] = E0.z;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i < 1 + additionalBinNumber) {
							bunemanEz[i][j][k] = 0;
						}
					}
					bunemanNewEz[i][j][k] = bunemanEz[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBx[i][j][k] = B0.x;
					bunemanNewBx[i][j][k] = bunemanBx[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBy[i][j][k] = B0.y;
					bunemanNewBy[i][j][k] = bunemanBy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanBz[i][j][k] = B0.z + E * cos(kw * xgrid[i]);
					bunemanNewBz[i][j][k] = bunemanBz[i][j][k];
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					Efield[i][j][k] = E0;
					Efield[i][j][k].y += E * cos(kw * xgrid[i]);
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i <= 1 + additionalBinNumber) {
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

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					Bfield[i][j][k] = B0;
					Bfield[i][j][k].z += E * cos(kw * middleXgrid[i]);
					newBfield[i][j][k] = Bfield[i][j][k];
				}
			}
		}
	}


	evaluateParticleTypesAlpha();
	if (rank == 0) printf("creating particles\n");
	fflush(stdout);
	if (rank == 0) printLog("creating particles\n");
	int n = 0;
	//for (int i = 0; i < xnumber; ++i) {
	for (int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				//int maxParticlesPerBin = types[0].particlesPerBin;
				double x = xgrid[i] + 0.0001 * deltaX;
				double y = ygrid[j] + 0.0001 * deltaY;
				double z = zgrid[k] + 0.0001 * deltaZ;
				//for (int l = 0; l < maxParticlesPerBin; ++l) {
				for (int typeCounter = 0; typeCounter < typesNumber; ++typeCounter) {
					double weight = (types[typeCounter].concentration / types[typeCounter].particlesPerBin) * volumeB(
					);
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
						//particle->coordinates.x = x + deltaXParticles * l;
						//particle->coordinates.y = y + deltaYParticles * l;
						//particle->coordinates.z = z + deltaZParticles * l;
						//particle->coordinates.x = middleXgrid[i];
						//particle->coordinates.y = middleYgrid[j];
						//particle->coordinates.z = middleZgrid[k];

						particle->coordinates.x = xgrid[i] + uniformDistribution() * deltaX;
						particle->coordinates.y = ygrid[j] + uniformDistribution() * deltaY;
						particle->coordinates.z = zgrid[k] + uniformDistribution() * deltaZ;
						particle->initialCoordinates = particle->coordinates;
						/*if (V0.norm() > 0) {
							particle->addVelocity(V0, speed_of_light_normalized);
						}*/
						if (type == 0) {
							particle->addVelocity(Vd, speed_of_light_normalized);
						}
						Vector3d momentum = particle->getMomentum();
						particle->initialMomentum = momentum;
						particle->prevMomentum = momentum;
						particles.push_back(particle);
						particlesNumber++;
						if (particlesNumber % 1000 == 0) {
							if ((rank == 0) && (verbosity > 0))printf("create particle number %d\n", particlesNumber);
						}
						alertNaNOrInfinity(particle->coordinates.x, "particle.x = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.y, "particle.y = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.z, "particle.z = NaN in createParticles\n");
					}
				}
				//create proton CRs
				double weight = (protonCRconcentration / protonCRperbin) * volumeB();
				for (int l = 0; l < protonCRperbin; ++l) {
					ParticleTypes type = types[1].type;
					Particle* particle = createParticle(n, i, j, k, weight, type, types[1],
					                                    types[1].temperatureX,
					                                    types[1].temperatureY,
					                                    types[1].temperatureZ);
					double phi = 2 * pi * uniformDistribution();
					double cosTheta = 2.0 * (uniformDistribution() - 0.5);
					double momentum = particle->mass * speed_of_light_normalized * sqrt(gammaCR * gammaCR - 1);
					double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
					particle->setMomentumX(momentum * cosTheta);
					particle->setMomentumY(momentum * sinTheta * cos(phi));
					particle->setMomentumZ(momentum * sinTheta * sin(phi));
					n++;
					//particle->coordinates.x = x + CRdeltaX * l;
					//particle->coordinates.y = y + CRdeltaY * l;
					//particle->coordinates.z = z + CRdeltaZ * l;
					//particle->coordinates.x = middleXgrid[i];
					//particle->coordinates.y = middleYgrid[j];
					//particle->coordinates.z = middleZgrid[k];

					particle->coordinates.x = xgrid[i] + uniformDistribution() * deltaX;
					particle->coordinates.y = ygrid[j] + uniformDistribution() * deltaY;
					particle->coordinates.z = zgrid[k] + uniformDistribution() * deltaZ;

					particle->initialCoordinates = particle->coordinates;
					particle->addVelocity(Vsh, speed_of_light_normalized);
					Vector3d particleMomentum = particle->getMomentum();
					particle->initialMomentum = particleMomentum;
					particle->prevMomentum = particleMomentum;
					particles.push_back(particle);
					particlesNumber++;
					if (particlesNumber % 1000 == 0) {
						if ((rank == 0) && (verbosity > 0))printf("create particle number %d\n", particlesNumber);
					}
					alertNaNOrInfinity(particle->coordinates.x, "particle.x = NaN in createParticles\n");
					alertNaNOrInfinity(particle->coordinates.y, "particle.y = NaN in createParticles\n");
					alertNaNOrInfinity(particle->coordinates.z, "particle.z = NaN in createParticles\n");
				}
			}
		}
	}


	synchronizeParticleNumber();


	if (rank == 0) {
		informationFile = fopen((outputDir + "information.dat").c_str(), "a");
		fprintf(informationFile, "CR gamma = %g\n", gammaCR);
		fprintf(informationFile, "V shock = %g\n", Vsh.x);
		fprintf(informationFile, "alfvenV = %g\n", alfvenV);
		fprintf(informationFile, "alfven Mach = %g\n", alfvenMach);
		fprintf(informationFile, "dzeta = %g\n", dzeta);
		fprintf(informationFile, "gyro radius CR = %g\n", rCR);
		fprintf(informationFile, "kmax = %g\n", kmax);
		double lmax = 2 * pi / kmax;
		fprintf(informationFile, "lmax = %g\n", lmax);
		double omega_pe = sqrt(
			4 * pi * types[0].concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
		fprintf(informationFile, "omega plasma electron = %g\n", omega_pe);
		double omega_pp = sqrt(
			4 * pi * types[1].concentration * electron_charge_normalized * electron_charge_normalized / massProton);
		fprintf(informationFile, "omega plasma electron = %g\n", omega_pp);
		double electronSkinDepth = speed_of_light_normalized / omega_pe;
		fprintf(informationFile, "electron skin depth = %g\n", electronSkinDepth);
		double protonSkinDepth = speed_of_light_normalized / omega_pp;
		fprintf(informationFile, "proton skin depth = %g\n", protonSkinDepth);
		double debyeLength = 1 / sqrt(
			4 * pi * electron_charge_normalized * electron_charge_normalized * (types[0].concentration + types[1].concentration)
			/ (kBoltzman_normalized * temperature));
		fprintf(informationFile, "debye length = %g\n", debyeLength);
		fprintf(informationFile, "dx/rCR = %g\n", deltaX / rCR);
		fprintf(informationFile, "rCR/Lx = %g\n", rCR / xsizeGeneral);
		fprintf(informationFile, "dx/lmax = %g\n", deltaX / lmax);
		fprintf(informationFile, "lmax/Lx = %g\n", lmax / xsizeGeneral);
		fprintf(informationFile, "dx/electronSkinDepth = %g\n", deltaX / electronSkinDepth);
		fprintf(informationFile, "electronSkinDepth/Lx = %g\n", electronSkinDepth / xsizeGeneral);
		fprintf(informationFile, "dx/protonSkinDepth = %g\n", deltaX / protonSkinDepth);
		fprintf(informationFile, "protonSkinDepth/Lx = %g\n", protonSkinDepth / xsizeGeneral);
		fprintf(informationFile, "dx/debyeLentgth = %g\n", deltaX / debyeLength);
		fprintf(informationFile, "debyeLentgth/Lx = %g\n", debyeLength / xsizeGeneral);
		fclose(informationFile);

		double increment = sqrt((dzeta * Vsh.scalarMult(Vsh) * kmax / rCR) - alfvenV * alfvenV * kmax * kmax);

		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", increment, increment / plasma_period, lmax, lmax * scaleFactor);
		fclose(incrementFile);
	}
}

void Simulation::initializeTestOneParticle() {
	//boundaryConditionTypeX = SUPER_CONDUCTOR_LEFT;
	boundaryConditionTypeX = PERIODIC;

	//initializeAlfvenWaveY(10, 1.0E-4);
	if (solverType == BUNEMAN) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEx[i][j][k] = E0.x;
					bunemanNewEx[i][j][k] = bunemanEx[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEy[i][j][k] = E0.y;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i < 1 + additionalBinNumber) {
							bunemanEy[i][j][k] = 0;
						}
					}
					bunemanNewEy[i][j][k] = bunemanEy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanEz[i][j][k] = E0.z;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i < 1 + additionalBinNumber) {
							bunemanEz[i][j][k] = 0;
						}
					}
					bunemanNewEz[i][j][k] = bunemanEz[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBx[i][j][k] = B0.x;
					bunemanNewBx[i][j][k] = bunemanBx[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBy[i][j][k] = B0.y;
					bunemanNewBy[i][j][k] = bunemanBy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanBz[i][j][k] = B0.z;
					bunemanNewBz[i][j][k] = bunemanBz[i][j][k];
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					Efield[i][j][k] = E0;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i <= 1 + additionalBinNumber) {
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
	}

	if (rank == 0) {
		double weight = (types[0].concentration / types[0].particlesPerBin) * volumeB();
		ParticleTypes type = types[0].type;
		Particle* particle = createParticle(0, 1 + additionalBinNumber, 1 + additionalBinNumber, 1 + additionalBinNumber,
		                                    weight, type, types[0], 0, 0, 0);

		particle->coordinates.x = middleXgrid[xnumberAdded / 2];
		particle->coordinates.y = middleYgrid[ynumberAdded / 2];
		particle->coordinates.z = middleZgrid[znumberAdded / 2];

		//particle->coordinates.x = xgrid[i] + uniformDistribution()*deltaX;
		//particle->coordinates.y = ygrid[j] + uniformDistribution()*deltaY;
		//particle->coordinates.z = zgrid[k] + uniformDistribution()*deltaZ;
		particle->initialCoordinates = particle->coordinates;
		if (V0.norm() > 0) {
			particle->addVelocity(V0, speed_of_light_normalized);
		}
		Vector3d momentum = particle->getMomentum();
		particle->initialMomentum = momentum;
		particle->prevMomentum = momentum;
		particles.push_back(particle);
		particlesNumber++;
	}
	synchronizeParticleNumber();

	MPI_Barrier(cartComm);
	if (rank == 0) {
		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", 0.0, 0.0, 1.0, 1.0);
		fclose(incrementFile);
	}
	MPI_Barrier(cartComm);
}

void Simulation::createParticles() {
	evaluateParticleTypesAlpha();
	if (rank == 0) printf("creating particles\n");
	fflush(stdout);
	if (rank == 0) printLog("creating particles\n");
	int n = 0;
	//for (int i = 0; i < xnumber; ++i) {
	for (int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				//int maxParticlesPerBin = types[0].particlesPerBin;
				double x = xgrid[i] + 0.0001 * deltaX;
				double y = ygrid[j] + 0.0001 * deltaY;
				double z = zgrid[k] + 0.0001 * deltaZ;
				//for (int l = 0; l < maxParticlesPerBin; ++l) {
				for (int typeCounter = 0; typeCounter < typesNumber; ++typeCounter) {
					double weight = (types[typeCounter].concentration / types[typeCounter].particlesPerBin) * volumeB(
					);
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
						//particle->coordinates.x = middleXgrid[i];
						//particle->coordinates.y = middleYgrid[j];
						//particle->coordinates.z = middleZgrid[k];

						//particle->coordinates.x = xgrid[i] + uniformDistribution()*deltaX;
						//particle->coordinates.y = ygrid[j] + uniformDistribution()*deltaY;
						//particle->coordinates.z = zgrid[k] + uniformDistribution()*deltaZ;
						particle->initialCoordinates = particle->coordinates;
						if (V0.norm() > 0) {
							particle->addVelocity(V0, speed_of_light_normalized);
						}
						Vector3d momentum = particle->getMomentum();
						particle->initialMomentum = momentum;
						particle->prevMomentum = momentum;
						particles.push_back(particle);
						particlesNumber++;
						if (particlesNumber % 1000 == 0) {
							if ((rank == 0) && (verbosity > 0))printf("create particle number %d\n", particlesNumber);
						}
						alertNaNOrInfinity(particle->coordinates.x, "particle.x = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.y, "particle.y = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.z, "particle.z = NaN in createParticles\n");
					}
				}
			}
		}
	}


	synchronizeParticleNumber();

	//printf("rank = %d p0 = %d p1 = %d p2 = %d e0 = %d e1 = %d e2 = %d Np = %d\n", rank, protonNumber, protonNumber1, protonNumber2, electronNumber, electronNumber1, electronNumber2, particlesNumber);

	/*if (preserveChargeLocal) {
		moveToPreserveChargeLocal();
	}*/
}

void Simulation::createParticlesHarris(double harrisWidth) {
	evaluateParticleTypesAlpha();
	double concentration = types[0].concentration;
	types[1].particesDeltaX = types[0].particesDeltaX;
	types[1].particlesPerBin = types[0].particlesPerBin;
	types[0].concentration = concentration;
	types[1].concentration = concentration;
	for (int i = 2; i < typesNumber; ++i) {
		types[i].particlesPerBin = 0;
		types[i].concentration = 0;
		types[i].particesDeltaX = xsize;
	}

	double temperature = B0.scalarMult(B0)/(16*pi*concentration*kBoltzman_normalized);
	if (rank == 0) printf("creating particles\n");
	if(verbosity > 2) printf("creating particles rank = %d\n", rank);
	fflush(stdout);
	if (rank == 0) printLog("creating particles\n");
	int n = 0;
	//for (int i = 0; i < xnumber; ++i) {
	for (int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				//int maxParticlesPerBin = types[0].particlesPerBin;
				double x = xgrid[i] + 0.0001 * deltaX;
				double y = ygrid[j] + 0.0001 * deltaY;
				double z = zgrid[k] + 0.0001 * deltaZ;
				//for (int l = 0; l < maxParticlesPerBin; ++l) {
				for (int typeCounter = 0; typeCounter < typesNumber; ++typeCounter) {
					double localConcentration = concentration/sqr(cosh((middleXgrid[i] - xsizeGeneral*1.5)/harrisWidth));
					double weight = (localConcentration / types[typeCounter].particlesPerBin) * volumeB(
					);
					double deltaXParticles = types[typeCounter].particesDeltaX;
					double deltaYParticles = types[typeCounter].particesDeltaY;
					double deltaZParticles = types[typeCounter].particesDeltaZ;
					//if (l < types[typeCounter].particlesPerBin) {
					for (int l = 0; l < types[typeCounter].particlesPerBin; ++l) {
						ParticleTypes type = types[typeCounter].type;
						Particle* particle = createParticle(n, i, j, k, weight, type, types[typeCounter],
						                                    temperature,
						                                    temperature,
						                                    temperature);
						n++;
						particle->coordinates.x = x + deltaXParticles * l;
						particle->coordinates.y = y + deltaYParticles * l;
						particle->coordinates.z = z + deltaZParticles * l;
						//particle->coordinates.x = middleXgrid[i];
						//particle->coordinates.y = middleYgrid[j];
						//particle->coordinates.z = middleZgrid[k];

						//particle->coordinates.x = xgrid[i] + uniformDistribution()*deltaX;
						//particle->coordinates.y = ygrid[j] + uniformDistribution()*deltaY;
						//particle->coordinates.z = zgrid[k] + uniformDistribution()*deltaZ;
						particle->initialCoordinates = particle->coordinates;

						Vector3d localV = Vector3d(0, 0, 1)*B0.norm()*speed_of_light_normalized/(8*pi*concentration*harrisWidth*types[typeCounter].charge);
						
						particle->addVelocity(localV, speed_of_light_normalized);
						
						Vector3d momentum = particle->getMomentum();
						particle->initialMomentum = momentum;
						particle->prevMomentum = momentum;
						particles.push_back(particle);
						particlesNumber++;
						if (particlesNumber % 1000 == 0) {
							if ((rank == 0) && (verbosity > 0))printf("create particle number %d\n", particlesNumber);
						}
						alertNaNOrInfinity(particle->coordinates.x, "particle.x = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.y, "particle.y = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.z, "particle.z = NaN in createParticles\n");
					}
				}
			}
		}
	}


	synchronizeParticleNumber();
	if(verbosity > 2) printf("finish creating particles rank = %d\n", rank);
	//printf("rank = %d p0 = %d p1 = %d p2 = %d e0 = %d e1 = %d e2 = %d Np = %d\n", rank, protonNumber, protonNumber1, protonNumber2, electronNumber, electronNumber1, electronNumber2, particlesNumber);

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

	for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
		for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
			if (cartJ == cartCoord[1] && cartK == cartCoord[2] && cartCoord[0] == cartDim[0] - 1) {
				for (int i = 0; i < escapedParticlesRight.size(); ++i) {
					Particle* particle = escapedParticlesRight[i];
					int typeNumber = getTypeNumber(particle);
					ParticleTypeContainer type = types[typeNumber];

					Particle* newParticle = createParticle(particlesNumber, xnumberAdded - additionalBinNumber - 2, 0, 0,
					                                       particle->weight, type.type, type, type.temperatureX, type.temperatureY,
					                                       type.temperatureZ);
					newParticle->coordinates.x = xgrid[xnumberAdded - 1 - additionalBinNumber];
					newParticle->coordinates.y = particle->coordinates.y;
					if (newParticle->coordinates.y < ygrid[1 + additionalBinNumber]) {
						newParticle->coordinates.y = ygrid[1 + additionalBinNumber] + deltaY * uniformDistribution();
					}
					if (newParticle->coordinates.y > ygrid[ynumberAdded - 1 - additionalBinNumber]) {
						newParticle->coordinates.y = ygrid[ynumberAdded - 1 - additionalBinNumber] - deltaY * uniformDistribution();
					}
					newParticle->coordinates.z = particle->coordinates.z;
					if (newParticle->coordinates.z < zgrid[1 + additionalBinNumber]) {
						newParticle->coordinates.z = zgrid[1 + additionalBinNumber] + deltaZ * uniformDistribution();
					}
					if (newParticle->coordinates.z > zgrid[znumberAdded - 1 - additionalBinNumber]) {
						newParticle->coordinates.z = zgrid[znumberAdded - 1 - additionalBinNumber] - deltaZ * uniformDistribution();
					}
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
		MPI_Recv(curParticleNumber, 1, MPI_INT, rank - 1, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
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
		MPI_Send(curParticleNumber, 1, MPI_INT, rank + 1, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm);
	}

	MPI_Barrier(cartComm);

	if (rank > 0) {
		MPI_Send(tempParticleNumber, 1, MPI_INT, 0, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm);
	} else {
		particleNumber = tempParticleNumber[0];
		for (int i = 1; i < nprocs; ++i) {
			MPI_Status status;
			MPI_Recv(tempParticleNumber, 1, MPI_INT, i, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm, &status);
			if (tempParticleNumber[0] >= 0) {
				particleNumber = tempParticleNumber[0];
			}
		}
	}

	tempParticleNumber[0] = particleNumber;

	MPI_Bcast(tempParticleNumber, 1, MPI_INT, 0, cartComm);

	particleNumber = tempParticleNumber[0];

	return particleNumber;
}


Particle* Simulation::createParticle(int n, int i, int j, int k, const double& weight,
                                     const ParticleTypes& type, const ParticleTypeContainer& typeContainer,
                                     const double& localTemperatureX,
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

	double thetaParamter = kBoltzman_normalized * (localTemperatureX + localTemperatureY + localTemperatureZ) / (3 * mass *
		speed_of_light_normalized_sqr);

	if (thetaParamter < 0.01) {
		//energy = maxwellDistribution(localTemparature, kBoltzman_normalized);
		px = sqrt(mass * kBoltzman_normalized * localTemperatureX) * normalDistribution();
		py = sqrt(mass * kBoltzman_normalized * localTemperatureY) * normalDistribution();
		pz = sqrt(mass * kBoltzman_normalized * localTemperatureZ) * normalDistribution();
		p = sqrt(px * px + py * py + pz * pz);
	} else if (localTemperatureX == localTemperatureY && localTemperatureX == localTemperatureZ) {
		energy = maxwellJuttnerDistribution(localTemperatureX, mass, speed_of_light_normalized, kBoltzman_normalized);
		p = sqrt(energy * energy - sqr(mass * speed_of_light_normalized_sqr)) / speed_of_light_normalized;


		//p = 0;

		pz = p * (2 * uniformDistribution() - 1);
		double phi = 2 * pi * uniformDistribution();
		double pnormal = sqrt(p * p - pz * pz);
		px = pnormal * cos(phi);
		py = pnormal * sin(phi);
	} else if (localTemperatureY == localTemperatureZ) {
		double momentumParallel;
		double momentumNormal;
		anisotropicMaxwellJuttnerDistribution(momentumNormal, momentumParallel, types[typeContainer.number].alphaNormal,
		                                      types[typeContainer.number].alphaParallel,
		                                      mass * speed_of_light_normalized);


		double sign = uniformDistribution() - 0.5;
		double phi = 2 * pi * uniformDistribution();
		if (sign > 0) {
			px = momentumParallel;
		} else {
			px = -momentumParallel;
		}
		py = momentumNormal * cos(phi);
		pz = momentumNormal * sin(phi);
		alertNaNOrInfinity(px, "px = NaN\n");
		alertNaNOrInfinity(py, "py = NaN\n");
		alertNaNOrInfinity(pz, "pz = NaN\n");
	} else {
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
				                                                                                   kBoltzman_normalized * types[i].
				                                                                                   temperatureX / (types[i].mass *
					                                                                                   speed_of_light_normalized_sqr));
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
				} else {
					maxAlpha = alphaNormal;
					minAlpha = 0.1 * types[i].mass * speed_of_light_normalized_sqr / (kBoltzman_normalized * types[i].temperatureX);
				}
				double simpleAlphaParallel = types[i].mass * speed_of_light_normalized_sqr / (kBoltzman_normalized * types[i].
					temperatureX);
				double alphaParallel = solver.solve(minAlpha, maxAlpha);
				types[i].alphaNormal = alphaNormal;
				types[i].alphaParallel = alphaParallel;
			} else {
				types[i].alphaNormal = 1.0;
				types[i].alphaParallel = 1.0;
			}
			alphas[2 * i] = types[i].alphaNormal;
			alphas[2 * i + 1] = types[i].alphaParallel;
		}
	}
	MPI_Barrier(cartComm);
	MPI_Bcast(alphas, 2 * typesNumber, MPI_DOUBLE, 0, cartComm);

	MPI_Barrier(cartComm);
	for (int i = 0; i < typesNumber; ++i) {
		types[i].alphaNormal = alphas[2 * i];
		types[i].alphaParallel = alphas[2 * i + 1];
	}
	MPI_Barrier(cartComm);

	delete[] alphas;
}

void Simulation::synchronizeParticleNumber() {
	if (nprocs > 1) {
		int particleCount[1];
		particleCount[0] = 0;
		if (rank > 0) {
			particleCount[0] = particles.size();
			MPI_Send(particleCount, 1, MPI_INT, 0, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm);
		} else {
			particlesNumber = particles.size();
			for (int i = 1; i < nprocs; ++i) {
				MPI_Status status;
				MPI_Recv(particleCount, 1, MPI_INT, i, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm, &status);
				particlesNumber += particleCount[0];
			}
		}
		particleCount[0] = particlesNumber;
		MPI_Bcast(particleCount, 1, MPI_INT, 0, cartComm);
		particlesNumber = particleCount[0];

		particleCount[0] = 0;

		if (rank > 0) {
			MPI_Status status;
			MPI_Recv(particleCount, 1, MPI_INT, rank - 1, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
		}
		for (int p = 0; p < particles.size(); ++p) {
			Particle* particle = particles[p];
			particle->number = particleCount[0] + p;
		}
		particleCount[0] = particleCount[0] + particles.size();
		if (rank < nprocs - 1) {
			MPI_Send(particleCount, 1, MPI_INT, rank + 1, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm);
		}
	} else {
		particlesNumber = particles.size();
		for (int p = 0; p < particles.size(); ++p) {
			Particle* particle = particles[p];
			particle->number = p;
		}
	}
}

int Simulation::getCartCoordWithAbsoluteIndexX(int i) {
	int isInThread[1];
	isInThread[0] = 0;
	int coordIn[1];
	coordIn[0] = -1;

	if (cartDim[0] == 1) {
		return 0;
	}

	if (i >= firstAbsoluteXindex + additionalBinNumber && i < firstAbsoluteXindex + xnumberAdded - 2 * additionalBinNumber
		- 1) {
		isInThread[0] = 1;
	} else {
		if ((cartCoord[0] == 0) && (i < 0)) {
			isInThread[0] = 1;
		} else if ((cartCoord[0] == cartDim[0] - 1) && (i >= firstAbsoluteXindex + xnumberAdded - 2 * additionalBinNumber - 1)
		) {
			isInThread[0] = 1;
		}
	}

	if (cartCoord[0] == 0) {
		if (isInThread[0] == 1) {
			coordIn[0] = cartCoord[0];
		}
		for (int tempI = 1; tempI < cartDim[0]; ++tempI) {
			int curRank = 0;
			int tempCoord[3];
			tempCoord[0] = tempI;
			tempCoord[1] = cartCoord[1];
			tempCoord[2] = cartCoord[2];
			MPI_Cart_rank(cartComm, tempCoord, &curRank);
			MPI_Status status;
			MPI_Recv(isInThread, 1, MPI_INT, curRank, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm, &status);
			if (isInThread[0] == 1) {
				coordIn[0] = tempI;
			}
		}
	} else {
		int firstRank = 0;
		int tempCoord[3];
		tempCoord[0] = 0;
		tempCoord[1] = cartCoord[1];
		tempCoord[2] = cartCoord[2];
		MPI_Cart_rank(cartComm, tempCoord, &firstRank);
		MPI_Send(isInThread, 1, MPI_INT, firstRank, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm);
	}

	if ((cartCoord[0] == 0) && (coordIn[0] == -1)) {
		printf("wrong cartCoord in getCoordWithAbsoluteX\n");
		MPI_Finalize();
		exit(0);
	}

	if (cartCoord[0] == 0) {
		for (int tempI = 1; tempI < cartDim[0]; ++tempI) {
			int curRank = 0;
			int tempCoord[3];
			tempCoord[0] = tempI;
			tempCoord[1] = cartCoord[1];
			tempCoord[2] = cartCoord[2];
			MPI_Cart_rank(cartComm, tempCoord, &curRank);
			MPI_Send(coordIn, 1, MPI_INT, curRank, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm);
		}
	} else {
		int firstRank = 0;
		int tempCoord[3];
		tempCoord[0] = 0;
		tempCoord[1] = cartCoord[1];
		tempCoord[2] = cartCoord[2];
		MPI_Cart_rank(cartComm, tempCoord, &firstRank);
		MPI_Status status;
		MPI_Recv(coordIn, 1, MPI_INT, firstRank, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm, &status);
	}

	return coordIn[0];
}

int Simulation::getCartCoordWithAbsoluteIndexY(int j) {
	int isInThread[1];
	isInThread[0] = 0;
	int coordIn[1];
	coordIn[0] = -1;

	if (cartDim[1] == 1) {
		return 0;
	}

	if (j >= firstAbsoluteYindex + additionalBinNumber && j < firstAbsoluteYindex + ynumberAdded - 2 * additionalBinNumber
		- 1) {
		isInThread[0] = 1;
	} else {
		if ((cartCoord[1] == 0) && (j < 0)) {
			isInThread[0] = 1;
		} else if ((cartCoord[1] == cartDim[1] - 1) && (j >= firstAbsoluteYindex + ynumberAdded - 2 * additionalBinNumber - 1)
		) {
			isInThread[0] = 1;
		}
	}

	if (cartCoord[1] == 0) {
		if (isInThread[0] == 1) {
			coordIn[0] = cartCoord[1];
		}
		for (int tempJ = 1; tempJ < cartDim[1]; ++tempJ) {
			int curRank = 0;
			int tempCoord[3];
			tempCoord[1] = tempJ;
			tempCoord[0] = cartCoord[0];
			tempCoord[2] = cartCoord[2];
			MPI_Cart_rank(cartComm, tempCoord, &curRank);
			MPI_Status status;
			MPI_Recv(isInThread, 1, MPI_INT, curRank, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm, &status);
			if (isInThread[0] == 1) {
				coordIn[0] = tempJ;
			}
		}
	} else {
		int firstRank = 0;
		int tempCoord[3];
		tempCoord[1] = 0;
		tempCoord[0] = cartCoord[0];
		tempCoord[2] = cartCoord[2];
		MPI_Cart_rank(cartComm, tempCoord, &firstRank);
		MPI_Send(isInThread, 1, MPI_INT, firstRank, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm);
	}

	if ((cartCoord[1] == 0) && (coordIn[0] == -1)) {
		printf("wrong cartCoord in getCoordWithAbsoluteY\n");
		MPI_Finalize();
		exit(0);
	}

	if (cartCoord[1] == 0) {
		for (int tempJ = 1; tempJ < cartDim[1]; ++tempJ) {
			int curRank = 0;
			int tempCoord[3];
			tempCoord[1] = tempJ;
			tempCoord[0] = cartCoord[0];
			tempCoord[2] = cartCoord[2];
			MPI_Cart_rank(cartComm, tempCoord, &curRank);
			MPI_Send(coordIn, 1, MPI_INT, curRank, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm);
		}
	} else {
		int firstRank = 0;
		int tempCoord[3];
		tempCoord[1] = 0;
		tempCoord[0] = cartCoord[0];
		tempCoord[2] = cartCoord[2];
		MPI_Cart_rank(cartComm, tempCoord, &firstRank);
		MPI_Status status;
		MPI_Recv(coordIn, 1, MPI_INT, firstRank, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm, &status);
	}

	return coordIn[0];
}

int Simulation::getCartCoordWithAbsoluteIndexZ(int k) {
	int isInThread[1];
	isInThread[0] = 0;
	int coordIn[1];
	coordIn[0] = -1;

	if (cartDim[2] == 1) {
		return 0;
	}

	if (k >= firstAbsoluteZindex + additionalBinNumber && k < firstAbsoluteZindex + znumberAdded - 2 * additionalBinNumber
		- 1) {
		isInThread[0] = 1;
	} else {
		if ((cartCoord[2] == 0) && (k < 0)) {
			isInThread[0] = 1;
		} else if ((cartCoord[2] == cartDim[2] - 1) && (k >= firstAbsoluteZindex + znumberAdded - 2 * additionalBinNumber - 1)
		) {
			isInThread[0] = 1;
		}
	}

	if (cartCoord[2] == 0) {
		if (isInThread[0] == 1) {
			coordIn[0] = cartCoord[2];
		}
		for (int tempI = 1; tempI < cartDim[2]; ++tempI) {
			int curRank = 0;
			int tempCoord[3];
			tempCoord[2] = tempI;
			tempCoord[1] = cartCoord[1];
			tempCoord[0] = cartCoord[0];
			MPI_Cart_rank(cartComm, tempCoord, &curRank);
			MPI_Status status;
			MPI_Recv(isInThread, 1, MPI_INT, curRank, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm, &status);
			if (isInThread[0] == 1) {
				coordIn[0] = tempI;
			}
		}
	} else {
		int firstRank = 0;
		int tempCoord[3];
		tempCoord[2] = 0;
		tempCoord[1] = cartCoord[1];
		tempCoord[0] = cartCoord[0];
		MPI_Cart_rank(cartComm, tempCoord, &firstRank);
		MPI_Send(isInThread, 1, MPI_INT, firstRank, MPI_SEND_INTEGER_ALL_TO_FIRST, cartComm);
	}

	if ((cartCoord[2] == 0) && (coordIn[0] == -1)) {
		printf("wrong cartCoord in getCoordWithAbsoluteZ\n");
		MPI_Finalize();
		exit(0);
	}

	if (cartCoord[2] == 0) {
		for (int tempI = 1; tempI < cartDim[2]; ++tempI) {
			int curRank = 0;
			int tempCoord[3];
			tempCoord[2] = tempI;
			tempCoord[1] = cartCoord[1];
			tempCoord[0] = cartCoord[0];
			MPI_Cart_rank(cartComm, tempCoord, &curRank);
			MPI_Send(coordIn, 1, MPI_INT, curRank, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm);
		}
	} else {
		int firstRank = 0;
		int tempCoord[3];
		tempCoord[2] = 0;
		tempCoord[1] = cartCoord[1];
		tempCoord[0] = cartCoord[0];
		MPI_Cart_rank(cartComm, tempCoord, &firstRank);
		MPI_Status status;
		MPI_Recv(coordIn, 1, MPI_INT, firstRank, MPI_SEND_INTEGER_FIRST_TO_ALL, cartComm, &status);
	}

	return coordIn[0];
}

int Simulation::getLocalIndexByAbsoluteX(int i) {
	return i - firstAbsoluteXindex;
}

int Simulation::getLocalIndexByAbsoluteY(int j) {
	return j - firstAbsoluteYindex;
}

int Simulation::getLocalIndexByAbsoluteZ(int k) {
	return k - firstAbsoluteZindex;
}

