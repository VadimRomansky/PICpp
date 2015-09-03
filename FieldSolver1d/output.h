#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "vector"
#include "stdio.h"

#include "vector3d.h"
#include "simulation.h"

void outputGrid(FILE* outFile, double* grid, int number);
void outputFields(FILE* outEfile, FILE* outBfile, Vector3d* Efield, Vector3d* Bfield, int xnumber, double plasma_priod, double gyroradius);
void outputGeneral(FILE* outFile, Simulation* simulatiom);
void outputDivergenceError(FILE* outFile, Simulation* simulation);
void outputGeneral(FILE* outFile, Simulation* simulation);


#endif