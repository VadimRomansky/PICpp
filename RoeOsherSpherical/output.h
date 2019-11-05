#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include "stdio.h"
#include "simulation.h"

void output(FILE* outFile, Simulation* simulation);
void outputDistribution(FILE* distributionFile, FILE* fullDistributionFile, FILE* coordinateDistributionFile, Simulation* simulation);
void outputDistributionP3(FILE* distributionFile, FILE* fullDistributionFile, FILE* coordinateDistributionFile, Simulation* simulation);
void outputNewGrid(FILE* outFile, Simulation* simulation);
void outMatrix(double* a, double* c, double* b, int N, double* f, double* x);
void outputField(FILE* outFile, FILE* coordinateFile, FILE* outFull, FILE* coefFile, FILE* xfile, FILE* kfile, FILE* pfile, Simulation* simulation);
void serialize(FILE* hydroFile, FILE* distributionFile, FILE* fieldFile, FILE* gridFile, FILE* pgridFile, FILE* kgridFile, FILE* infoFile, Simulation* simulation);

#endif