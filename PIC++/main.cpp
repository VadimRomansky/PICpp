#include "stdio.h"
#include <stdlib.h>
#include <time.h>
//#include <crtdbg.h>
#include <omp.h>

#include "constants.h"
#include "matrix3d.h"
#include "util.h"
#include "input.h"
#include "simulation.h"

int main() {
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_ALWAYS_DF);
	//omp_set_num_threads(numThreads);
	printf("start\n");

	printf("random initialize\n");
	srand(time(NULL));

	bool startNew = true;

	if (startNew) {
		printf("open input\n");
		std::string inputDir = inputDirectory;
		/*FILE* tempFile = fopen((inputDir + "temp.dat").c_str(), "w");
		fclose(tempFile);*/
		FILE* inputFile = fopen((inputDir + "input.dat").c_str(), "r");
		printf("read input\n");
		Simulation simulation = readInput(inputFile);
		printf("close input\n");
		fclose(inputFile);

		printf("simulate\n");
		simulation.simulate();
		printf("end\n");
	} else {
		printf("reading backup\n");

		std::string backupDir = backupDirectory;
		FILE* backupGeneralFile = fopen((backupDir + "general.dat").c_str(), "r");
		FILE* backupEfieldFile = fopen((backupDir + "Efield.dat").c_str(), "r");
		FILE* backupBfieldFile = fopen((backupDir + "Bfield.dat").c_str(), "r");
		FILE* backupParticlesFile = fopen((backupDir + "particles.dat").c_str(), "r");

		Simulation simulation = readBackup(backupGeneralFile, backupEfieldFile, backupBfieldFile, backupParticlesFile);

		fclose(backupGeneralFile);
		fclose(backupEfieldFile);
		fclose(backupBfieldFile);
		fclose(backupParticlesFile);

		printf("simulate\n");
		simulation.simulate();
		printf("end\n");
	}
}

