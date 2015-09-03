#include "stdio.h"
#include "math.h"
#include "vector"

#include "output.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "util.h"

void outputGrid(FILE* outFile, double* grid, int number) {
	for(int i = 0; i < number; ++i) {
		fprintf(outFile, "%15.10g %15.10g\n", grid[i], (grid[i+1] + grid[i])/2);
	}
}

void outputFields(FILE* outEfile, FILE* outBfile, Vector3d* Efield, Vector3d* Bfield, int xnumber, double plasma_preiod, double gyroradius) {
	double scale = 1.0/(plasma_preiod*gyroradius);
	for(int i = 0; i < xnumber; ++i) {
				fprintf(outBfile, "%15.10g %15.10g %15.10g\n", scale*Bfield[i].x, scale*Bfield[i].y, scale*Bfield[i].z);
	}

	for(int i = 0; i < xnumber + 1; ++i) {
				fprintf(outEfile, "%15.10g %15.10g %15.10g\n", scale*Efield[i].x, scale*Efield[i].y, scale*Efield[i].z);
	}
}

void outputArrayParameter(FILE* outFile, double*** arrayParameter, int xnumber, int ynumber, int znumber) {
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%g\n", arrayParameter[i][j][k]);
			}
		}
	}
}

void outputGeneral(FILE* outFile, Simulation* simulation) {
	fprintf(outFile, "%d %g %g %g %g %g %g %g %g %g\n", simulation->currentIteration, simulation->time, simulation->time*simulation->plasma_period, simulation->particleEnergy,
		simulation->electricFieldEnergy, simulation->magneticFieldEnergy, simulation->energy, simulation->momentum.x, simulation->momentum.y, simulation->momentum.z);
}

void outputDivergenceError(FILE* outFile, Simulation* simulation) {
	for(int i = 0; i < simulation->xnumber; ++i) {
				double div = simulation->evaluateDivE(i);
				double div2 = simulation->evaluateDivTempE(i);
				fprintf(outFile, "%g %g\n", div, div2);
	}
}