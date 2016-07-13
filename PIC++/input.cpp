#include "stdio.h"
#include "simulation.h"

#include "input.h"
#include "constants.h"

Simulation readInput(FILE* inputFile) {
	std::string outputDir = outputDirectory;
	int inputType;
	char ch = ' ';

	fscanf(inputFile, "%d", &inputType);

	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	int xnumber;
	fscanf(inputFile, "%d", &xnumber);

	if (xnumber < 0) {
		printf("xnumber must be > 0\n");
		FILE *errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "xnumber must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	int ynumber;
	fscanf(inputFile, "%d", &ynumber);

	if (ynumber < 0) {
		printf("ynumber must be > 0\n");
		FILE *errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "ynumber must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	int znumber;
	fscanf(inputFile, "%d", &znumber);
	if (znumber < 0) {
		printf("znumber must be > 0\n");
		FILE *errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "znumber must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double xsize;
	fscanf(inputFile, "%lf", &xsize);
	if (xsize < 0) {
		printf("xsize must be > 0\n");
		FILE *errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "xsize must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	double ysize;
	fscanf(inputFile, "%lf", &ysize);
	if (ysize < 0) {
		printf("ysize must be > 0\n");
		FILE *errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "ysize must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	double zsize;
	fscanf(inputFile, "%lf", &zsize);
	if (zsize < 0) {
		printf("zsize must be > 0\n");
		FILE *errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "zsize must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double temperature;
	fscanf(inputFile, "%lf", &temperature);

	if (temperature < 0) {
		printf("temperature must be > 0\n");
		FILE *errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "v/c > 1 in setMomentumByV\n");
		fclose(errorLogFile);
		exit(0);
	}

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double Vx;
	fscanf(inputFile, "%lf", &Vx);

	double Vy;
	fscanf(inputFile, "%lf", &Vy);

	double Vz;
	fscanf(inputFile, "%lf", &Vz);

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double Ex;
	fscanf(inputFile, "%lf", &Ex);

	double Ey;
	fscanf(inputFile, "%lf", &Ey);

	double Ez;
	fscanf(inputFile, "%lf", &Ez);

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double Bx;
	fscanf(inputFile, "%lf", &Bx);

	double By;
	fscanf(inputFile, "%lf", &By);

	double Bz;
	fscanf(inputFile, "%lf", &Bz);

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	int maxIterations;
	fscanf(inputFile, "%d", &maxIterations);

	if (maxIterations < 0) {
		printf("max iterations must be > 0\n");
		FILE *errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "maxIterations mast be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double maxTime;
	fscanf(inputFile, "%lf", &maxTime);

	if (maxTime < 0) {
		printf("max time must be > 0\n");
		FILE *errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "maxTime must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	int typesNumber;
	fscanf(inputFile, "%d", &typesNumber);

	if (typesNumber < 1) {
		printf("typesNumber > 1\n");
		FILE *errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "typesNumber > 1\n");
		fclose(errorLogFile);
		exit(0);
	}

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double *concentrations = new double[typesNumber];
	int *particlesPerBin = new int[typesNumber];

	for (int i = 0; i < typesNumber; ++i) {
		int particlePerBin;
		double concentration;
		fscanf(inputFile, "%lf %d", &concentration, &particlePerBin);
		concentrations[i] = concentration;
		particlesPerBin[i] = particlePerBin;

		if (particlesPerBin[i] < 0) {
			printf("particlesPerBin[i] must be > 0\n");
			FILE *errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "particlesPerBin[i] must be > 0\n");
			fclose(errorLogFile);
			exit(0);
		}

		ch = ' ';
		while (ch != '\n') {
			fscanf(inputFile, "%c", &ch);
		}
	}

	return Simulation(xnumber, ynumber, znumber, xsize, ysize, zsize, temperature, Vx, Vy, Vz, Ex, Ey, Ez, Bx,
					  By, Bz, maxIterations, maxTime, typesNumber, particlesPerBin, concentrations, inputType);
}

Simulation readBackup(FILE* generalFile, FILE* Efile, FILE* Bfile, FILE* particlesFile) {
	Simulation simulation = Simulation();
	int inputType = 0;
	fscanf(generalFile, "%d", &inputType);
	if(inputType == 0){
		simulation.inputType = CGS;
	} else if(inputType == 1){
		simulation.inputType = Theoretical;
	} else {
		printf("input type must be 1 or 0\n");
		FILE* errorLogFile = fopen((simulation.outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "input type must be 1 or 0\n");
		fclose(errorLogFile);
		exit(0);
	}
	fscanf(generalFile, "%d", &simulation.xnumber);
	fscanf(generalFile, "%d", &simulation.ynumber);
	fscanf(generalFile, "%d", &simulation.znumber);
	fscanf(generalFile, "%d", &simulation.particlesNumber);
	fscanf(generalFile, "%d", &simulation.types[0].particlesPerBin);
	fscanf(generalFile, "%d", &simulation.types[1].particlesPerBin);
	fscanf(generalFile, "%d", &simulation.types[2].particlesPerBin);
	fscanf(generalFile, "%d", &simulation.types[3].particlesPerBin);
	fscanf(generalFile, "%d", &simulation.types[4].particlesPerBin);
	fscanf(generalFile, "%d", &simulation.types[5].particlesPerBin);

	fscanf(generalFile, "%lf", &simulation.types[0].concentration);
	fscanf(generalFile, "%lf", &simulation.types[1].concentration);
	fscanf(generalFile, "%lf", &simulation.types[2].concentration);
	fscanf(generalFile, "%lf", &simulation.types[3].concentration);
	fscanf(generalFile, "%lf", &simulation.types[4].concentration);
	fscanf(generalFile, "%lf", &simulation.types[5].concentration);


	fscanf(generalFile, "%lf", &simulation.temperature);
	fscanf(generalFile, "%lf", &simulation.plasma_period);
	fscanf(generalFile, "%lf", &simulation.scaleFactor);
	fscanf(generalFile, "%lf", &simulation.fieldScale);

	fscanf(generalFile, "%lf", &simulation.time);
	fscanf(generalFile, "%lf", &simulation.maxTime);

	fscanf(generalFile, "%d", &simulation.currentIteration);
	fscanf(generalFile, "%d", &simulation.maxIteration);

	fscanf(generalFile, "%lf", &simulation.xsize);
	fscanf(generalFile, "%lf", &simulation.ysize);
	fscanf(generalFile, "%lf", &simulation.zsize);
	fscanf(generalFile, "%lf", &simulation.theta);
	fscanf(generalFile, "%lf", &simulation.eta);

	int debugMode = 0;

	fscanf(generalFile, "%d", &debugMode);

	simulation.debugMode = (debugMode == 1);

	int preserveChargeLocal = 0;

	fscanf(generalFile, "%d", &preserveChargeLocal);

	simulation.preserveChargeLocal = (preserveChargeLocal == 1);

	int preserveChargeGlobal = 0;

	fscanf(generalFile, "%d", &preserveChargeGlobal);

	simulation.preserveChargeGlobal = (preserveChargeGlobal == 1);

	int solverType = 0;

	fscanf(generalFile, "%d", &solverType);

	if (solverType == 1) {
		simulation.solverType = IMPLICIT;
	} else {
		simulation.solverType = EXPLICIT;
	}

	int boundaryConditionType = 0;

	fscanf(generalFile, "%d", &boundaryConditionType);

	if (boundaryConditionType == 1) {
		simulation.boundaryConditionType = PERIODIC;
	} else {
		simulation.boundaryConditionType = SUPER_CONDUCTOR_LEFT;
	}

	fscanf(generalFile, "%d", &simulation.maxwellEquationMatrixSize);

	fscanf(generalFile, "%lf", &simulation.extJ);

	fscanf(generalFile, "%lf", &simulation.V0.x);
	fscanf(generalFile, "%lf", &simulation.V0.y);
	fscanf(generalFile, "%lf", &simulation.V0.z);
	fscanf(generalFile, "%lf", &simulation.E0.x);
	fscanf(generalFile, "%lf", &simulation.E0.y);
	fscanf(generalFile, "%lf", &simulation.E0.z);
	fscanf(generalFile, "%lf", &simulation.B0.x);
	fscanf(generalFile, "%lf", &simulation.B0.y);
	fscanf(generalFile, "%lf", &simulation.B0.z);

	simulation.deltaX = simulation.xsize / (simulation.xnumber);
	simulation.deltaY = simulation.ysize / (simulation.ynumber);
	simulation.deltaZ = simulation.zsize / (simulation.znumber);

	simulation.deltaX2 = simulation.deltaX * simulation.deltaX;
	simulation.deltaY2 = simulation.deltaY * simulation.deltaY;
	simulation.deltaZ2 = simulation.deltaZ * simulation.deltaZ;

	simulation.createArrays();

	simulation.rescaleConstants();

	for (int i = 0; i <= simulation.xnumber; ++i) {
		simulation.xgrid[i] = i * simulation.deltaX;
	}

	for (int i = 0; i <= simulation.ynumber; ++i) {
		simulation.ygrid[i] = i * simulation.deltaY;
	}

	for (int i = 0; i <= simulation.znumber; ++i) {
		simulation.zgrid[i] = i * simulation.deltaZ;
	}

	for (int i = 0; i < simulation.xnumber; ++i) {
		simulation.middleXgrid[i] = (simulation.xgrid[i] + simulation.xgrid[i + 1]) / 2;
	}

	for (int i = 0; i < simulation.ynumber; ++i) {
		simulation.middleYgrid[i] = (simulation.ygrid[i] + simulation.ygrid[i + 1]) / 2;
	}

	for (int i = 0; i < simulation.znumber; ++i) {
		simulation.middleZgrid[i] = (simulation.zgrid[i] + simulation.zgrid[i + 1]) / 2;
	}

	readFields(Efile, Bfile, simulation);
	simulation.resetNewTempFields();
	readParticles(particlesFile, simulation);

	return simulation;
}

void readFields(FILE* Efile, FILE* Bfile, Simulation& simulation) {
	for (int i = 0; i < simulation.xnumber; ++i) {
		for (int j = 0; j < simulation.ynumber; ++j) {
			for (int k = 0; k < simulation.znumber; ++k) {
				fscanf(Bfile, "%lf %lf %lf", &(simulation.Bfield[i][j][k].x), &(simulation.Bfield[i][j][k].y), &(simulation.Bfield[i][j][k].z));
			}
		}
	}

	for (int i = 0; i < simulation.xnumber + 1; ++i) {
		for (int j = 0; j < simulation.ynumber + 1; ++j) {
			for (int k = 0; k < simulation.znumber + 1; ++k) {
				fscanf(Efile, "%lf %lf %lf", &(simulation.Efield[i][j][k].x), &(simulation.Efield[i][j][k].y), &(simulation.Efield[i][j][k].z));
			}
		}
	}
}

void readParticles(FILE* particlesFile, Simulation& simulation) {
	for (int i = 0; i < simulation.particlesNumber; ++i) {
		int number = 0;
		double mass = 0;
		double charge = 0;
		double weight = 0;
		int type;
		double x;
		double y;
		double z;
		double px;
		double py;
		double pz;
		double dx;
		double dy;
		double dz;

		fscanf(particlesFile, "%d %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &number, &mass, &charge, &weight, &type, &x, &y, &z, &px, &py, &pz, &dx, &dy, &dz);

		ParticleTypes particleType;
		if (type == 0) {
			particleType = ELECTRON;
		} else if (type == 1) {
			particleType = PROTON;
		} else if (type == 2) {
			particleType = POSITRON;
		} else if (type == 3) {
			particleType = ALPHA;
		} else if (type == 4) {
			particleType = DEUTERIUM;
		} else if (type == 5) {
			particleType = HELIUM3;
		} else {
			printf("particle type must be 0 1 2 3 4 5\n");
			exit(0);
		}
		ParticleTypeContainer typeContainer = simulation.types[type - 1];

		Particle* particle = new Particle(number, mass, typeContainer.chargeCount, typeContainer.charge, weight, particleType, typeContainer, x, y, z, px, py, pz, dx, dy, dz);
		simulation.chargeBalance += particle->chargeCount;
		simulation.particles.push_back(particle);
	}
}

