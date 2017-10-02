#include "stdio.h"
#include "simulation.h"

#include "input.h"

Simulation readInput(FILE* inputFile) {
	int xnumber;
	char ch = ' ';
	fscanf(inputFile, "%d", &xnumber);

	if (xnumber < 0) {
		printf("xnumber must be > 0\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "xnumber must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	int ynumber;
	fscanf(inputFile, "%d", &ynumber);

	if (ynumber < 0) {
		printf("ynumber must be > 0\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "ynumber must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	int znumber;
	fscanf(inputFile, "%d", &znumber);
	if (znumber < 0) {
		printf("znumber must be > 0\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "znumber must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double xsize;
	fscanf(inputFile, "%lf", &xsize);
	if (xsize < 0) {
		printf("xsize must be > 0\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "xsize must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	double ysize;
	fscanf(inputFile, "%lf", &ysize);
	if (ysize < 0) {
		printf("ysize must be > 0\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "ysize must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	double zsize;
	fscanf(inputFile, "%lf", &zsize);
	if (zsize < 0) {
		printf("zsize must be > 0\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
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
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "v/c > 1 in setMomentumByV\n");
		fclose(errorLogFile);
		exit(0);
	}

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double density;
	fscanf(inputFile, "%lf", &density);

	if (density < 0) {
		printf("density must be > 0\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "density must be > 0\n");
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
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
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
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "maxTime must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	int particlesPerBin;
	fscanf(inputFile, "%d", &particlesPerBin);

	if (particlesPerBin < 0) {
		printf("particles per bin must be > 0\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "particles per bin must be > 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	int positronsPerBin;
	fscanf(inputFile, "%d", &positronsPerBin);

	if (positronsPerBin < 0) {
		printf("positrons per bin must be >= 0\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "positrons per bin must be >= 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	int alphaPerBin;
	fscanf(inputFile, "%d", &alphaPerBin);

	if (alphaPerBin < 0) {
		printf("alpha per bin must be >= 0\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "alpha per bin must be >= 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	ch = ' ';
	while (ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	return Simulation(xnumber, xsize, temperature, density, Vx, Vy, Vz, Ex, Ey, Ez, Bx, By, Bz, maxIterations, maxTime, particlesPerBin, positronsPerBin, alphaPerBin);
}

Simulation readBackup(const char* generalFileName, const char* EfileName, const char* BfileName, const char* particlesFileName) {
	Simulation simulation = Simulation();

	int xnumber = 0;
	fscanf(generalFileName, "%d", &xnumber);
	simulation.xnumber = xnumber;
	fscanf(generalFileName, "%d", &simulation.particlesNumber);
	fscanf(generalFileName, "%d", &simulation.electronsPerBin);
	fscanf(generalFileName, "%d", &simulation.protonsPerBin);
	fscanf(generalFileName, "%d", &simulation.positronsPerBin);
	fscanf(generalFileName, "%d", &simulation.alphaPerBin);

	fscanf(generalFileName, "%lf", &simulation.density);
	fscanf(generalFileName, "%lf", &simulation.temperature);
	fscanf(generalFileName, "%lf", &simulation.plasma_period);
	fscanf(generalFileName, "%lf", &simulation.gyroradius);
	fscanf(generalFileName, "%lf", &simulation.fieldScale);

	fscanf(generalFileName, "%lf", &simulation.time);
	fscanf(generalFileName, "%lf", &simulation.maxTime);

	fscanf(generalFileName, "%d", &simulation.currentIteration);
	fscanf(generalFileName, "%d", &simulation.maxIteration);

	fscanf(generalFileName, "%lf", &simulation.xsize);
	fscanf(generalFileName, "%lf", &simulation.theta);
	fscanf(generalFileName, "%lf", &simulation.eta);

	int debugMode = 0;

	fscanf(generalFileName, "%d", &debugMode);

	simulation.debugMode = (debugMode == 1);

	int preserveCharge = 0;

	fscanf(generalFileName, "%d", &preserveCharge);

	simulation.preserveCharge = (preserveCharge == 1);

	int solverType = 0;

	fscanf(generalFileName, "%d", &solverType);

	if (solverType == 1) {
		simulation.solverType = IMPLICIT;
	} else {
		simulation.solverType = EXPLICIT;
	}

	int boundaryConditionType = 0;

	fscanf(generalFileName, "%d", &boundaryConditionType);

	if (boundaryConditionType == 1) {
		simulation.boundaryConditionType = PERIODIC;
	} else {
		simulation.boundaryConditionType = SUPER_CONDUCTOR_LEFT;
	}

	fscanf(generalFileName, "%d", &simulation.maxwellEquationMatrixSize);

	fscanf(generalFileName, "%lf", &simulation.extJ);

	fscanf(generalFileName, "%lf", &simulation.V0.x);
	fscanf(generalFileName, "%lf", &simulation.V0.y);
	fscanf(generalFileName, "%lf", &simulation.V0.z);
	fscanf(generalFileName, "%lf", &simulation.E0.x);
	fscanf(generalFileName, "%lf", &simulation.E0.y);
	fscanf(generalFileName, "%lf", &simulation.E0.z);
	fscanf(generalFileName, "%lf", &simulation.B0.x);
	fscanf(generalFileName, "%lf", &simulation.B0.y);
	fscanf(generalFileName, "%lf", &simulation.B0.z);

	simulation.deltaX = simulation.xsize / (simulation.xnumber);

	simulation.deltaX2 = simulation.deltaX * simulation.deltaX;

	simulation.createArrays();

	simulation.rescaleConstants();

	for (int i = 0; i <= simulation.xnumber; ++i) {
		simulation.xgrid[i] = i * simulation.deltaX;
	}

	for (int i = 0; i < simulation.xnumber; ++i) {
		simulation.middleXgrid[i] = (simulation.xgrid[i] + simulation.xgrid[i + 1]) / 2;
	}

	readFields(EfileName, BfileName, simulation);
	simulation.resetNewTempFields();
	readParticles(particlesFileName, simulation);

	return simulation;
}

void readFields(const char* Efile, const char* Bfile, Simulation& simulation) {
	for (int i = 0; i < simulation.xnumber; ++i) {
		fscanf(Bfile, "%lf %lf %lf", &(simulation.Bfield[i].x), &(simulation.Bfield[i].y), &(simulation.Bfield[i].z));
	}

	for (int i = 0; i < simulation.xnumber + 1; ++i) {
		fscanf(Efile, "%lf %lf %lf", &(simulation.Efield[i].x), &(simulation.Efield[i].y), &(simulation.Efield[i].z));
	}
}

void readParticles(const char* particlesFile, Simulation& simulation) {
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

		fscanf(particlesFile, "%d %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf", &number, &mass, &charge, &weight, &type, &x, &y, &z, &px, &py, &pz, &dx);

		ParticleTypes particleType;
		if (type == 1) {
			particleType = ELECTRON;
		} else if (type == 2) {
			particleType = PROTON;
		} else if (type == 3) {
			particleType = POSITRON;
		} else if (type == 4) {
			particleType = ALPHA;
		} else {
			printf("particle type must be 1 2 3 4\n");
			FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "particle type must be 1 2 3 4\n");
			fclose(errorLogFile);
			exit(0);
		}
		ParticleTypeContainer typeContainer = simulation.types[type - 1];

		Particle* particle = new Particle(number, mass, charge, weight, particleType, typeContainer, x, px, py, pz, dx);
		particle->y = y;
		particle->z = z;

		simulation.particles.push_back(particle);
	}
}

