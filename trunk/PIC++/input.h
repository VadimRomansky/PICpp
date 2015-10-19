#ifndef _INPUT_H_
#define _INPUT_H_

#include "stdio.h"
#include "simulation.h"

Simulation readInput(FILE* input_file);
Simulation readBackup(FILE* generalFile, FILE* Efile, FILE* Bfile, FILE* particlesFile);
void readFields(FILE* Efile, FILE* Bfile, Simulation& simulation);
void readParticles(FILE* particlesFile, Simulation& simulation);

#endif 