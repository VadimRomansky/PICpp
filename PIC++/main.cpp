#include "stdio.h"
#include <stdlib.h>
#include <time.h>
//#include <crtdbg.h>
#include <omp.h>

//#include "memory_debug.h"
#include "constants.h"
#include "particle.h"
#include "matrix3d.h"
#include "util.h"
#include "input.h"
#include "simulation.h"

int main(int argc, char **argv) {
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_ALWAYS_DF);
	int size;
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//printf("hello, world!");
	//exit(0);
	if(rank == 0){
		printf("start version 17.01.2017 23:00\n");
	}

	for(int i = 0; i < size; ++i) {
		MPI_Barrier(MPI_COMM_WORLD);
		if(i == rank) {
			printf("start thread number %d\n", rank);
			fflush(stdout);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0) printf("start\n");
	//fflush(stdout);

	if(rank == 0) printf("random initialize\n");
	//fflush(stdout);
	srand(time(NULL) + rank);

    /*double x = 5;
    double index = 2;
    double a = McDonaldFunction(x, index);

    printf("K(x = %lf, index = %lf) = %15.10g\n", x, index, a);

    return 0;*/

	bool startNew = true;
	bool writelog = false;

	if (startNew) {
		Simulation simulation;
		if(rank == 0) {
			printf("open input\n");
			//fflush(stdout);
			std::string inputDir = inputDirectory;
			FILE *inputFile = fopen((inputDir + "input.dat").c_str(), "r");
			printf("read input\n");
			//fflush(stdout);
			simulation = readInput(inputFile);
			printf("close input\n");
			//fflush(stdout);
			fclose(inputFile);
			simulation.rank = 0;
			sendInput(simulation, size);
		} else {
			//printf("recieve\n");
			simulation = recieveInput();
			simulation.rank = rank;
		}

		if(rank == 0) printf("simulate\n");
		//fflush(stdout);
		simulation.simulate();
		if(rank == 0) printf("end\n");
		fflush(stdout);
	} else {
		//todo parallel
		Simulation simulation;
			printf("reading backup\n");
			fflush(stdout);

			std::string backupDir = backupDirectory;
			const char* backupGeneralFileName = (backupDir + "general.dat").c_str();
			const char* backupEfieldFileName = (backupDir + "Efield.dat").c_str();
			const char* backupBfieldFileName = (backupDir + "Bfield.dat").c_str();
			const char* backupParticlesFileName = (backupDir + "particles.dat").c_str();

			simulation = readBackup(backupGeneralFileName, backupEfieldFileName, backupBfieldFileName, backupParticlesFileName);



		if(rank == 0) printf("simulate\n");
		fflush(stdout);
		simulation.simulate();
		if(rank == 0) printf("end\n");
		fflush(stdout);
	}
	MPI_Finalize();

	//std::string outDir = outputDirectory;
	//const char* memoryFileName = (outDir + "memory.dat").c_str();
	//FILE* memoryFile = fopen(memoryFileName,"w");
	//_CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
	//_CrtSetReportFile( _CRT_WARN, memoryFile );
	//_CrtDumpMemoryLeaks();
	//fclose(memoryFile);
}

