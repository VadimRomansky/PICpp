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
		FILE* inputFile = fopen("./input/input.dat", "r");
		printf("read input\n");
		Simulation simulation = readInput(inputFile);
		printf("close input\n");
		fclose(inputFile);

		printf("simulate\n");
		simulation.simulate();
		printf("end\n");
	} else {
		printf("reading backup\n");

		FILE* backupGeneralFile = fopen("./backup/general.dat", "r");
		FILE* backupEfieldFile = fopen("./backup/Efield.dat", "r");
		FILE* backupBfieldFile = fopen("./backup/Bfield.dat", "r");
		FILE* backupParticlesFile = fopen("./backup/particles.dat", "r");

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

