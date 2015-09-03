#include "stdio.h"
#include "simulation.h"

#include "input.h"

Simulation readInput(FILE* inputFile) {
	int xnumber;
	char ch = ' ';
	fscanf(inputFile, "%d", &xnumber);

	if(xnumber < 0) {
		printf("xnumber must be > 0\n");
		exit(0);
	}

	int ynumber;
	fscanf(inputFile, "%d", &ynumber);

	if(ynumber < 0) {
		printf("ynumber must be > 0\n");
		exit(0);
	}

	int znumber;
	fscanf(inputFile, "%d", &znumber);
	if(znumber < 0) {
		printf("znumber must be > 0\n");
		exit(0);
	}

	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double xsize;
	fscanf(inputFile, "%lf", &xsize);
	if(xsize < 0) {
		printf("xsize must be > 0\n");
		exit(0);
	}

	double ysize;
	fscanf(inputFile, "%lf", &ysize);
	if(ysize < 0) {
		printf("ysize must be > 0\n");
		exit(0);
	}

	double zsize;
	fscanf(inputFile, "%lf", &zsize);
	if(zsize < 0) {
		printf("zsize must be > 0\n");
		exit(0);
	}

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double temperature;
	fscanf(inputFile, "%lf", &temperature);

	if(temperature < 0) {
		printf("temperature must be > 0\n");
		exit(0);
	}

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double density;
	fscanf(inputFile, "%lf", &density);

	if(density < 0) {
		printf("density must be > 0\n");
		exit(0);
	}

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double Ex;
	fscanf(inputFile, "%lf", &Ex);

	double Ey;
	fscanf(inputFile, "%lf", &Ey);

	double Ez;
	fscanf(inputFile, "%lf", &Ez);

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double Bx;
	fscanf(inputFile, "%lf", &Bx);

	double By;
	fscanf(inputFile, "%lf", &By);

	double Bz;
	fscanf(inputFile, "%lf", &Bz);

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	int maxIterations;
	fscanf(inputFile, "%d", &maxIterations);

	if(maxIterations < 0) {
		printf("max iterations must be > 0\n");
		exit(0);
	}

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double maxTime;
	fscanf(inputFile, "%lf", &maxTime);

	if(maxTime < 0) {
		printf("max time must be > 0\n");
		exit(0);
	}

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	int particlesPerBin;
	fscanf(inputFile, "%d", &particlesPerBin);

	if(particlesPerBin < 0) {
		printf("particles per bin must be > 0\n");
		exit(0);
	}

	return Simulation(xnumber, xsize, temperature, density, Ex, Ey, Ez, Bx, By, Bz, maxIterations, maxTime, particlesPerBin);
}