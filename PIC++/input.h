#ifndef _INPUT_H_
#define _INPUT_H_

#include "stdio.h"
class Simulation;

Simulation readInput(FILE* input_file, MPI_Comm& comm);
Simulation readBackup(const char* generalFileName, const char* EfileName, const char* BfileName, const char* particlesFileName, MPI_Comm& comm);
void readFields(const char* Efile, const char* Bfile, Simulation& simulation);
void readParticles(const char* particlesFile, Simulation& simulation);
int** readTrackedParticlesNumbers(const char* particlesFile, int& number);

#endif 