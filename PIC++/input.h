#ifndef _INPUT_H_
#define _INPUT_H_

#include "stdio.h"
class Simulation;

Simulation readInput(FILE* input_file);
Simulation readBackup(const char* generalFileName, const char* EfileName, const char* BfileName, const char* particlesFileName);
void readFields(const char* Efile, const char* Bfile, Simulation& simulation);
void readParticles(const char* particlesFile, Simulation& simulation);

#endif 