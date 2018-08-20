#include "stdio.h"
#include <stdlib.h>
#include <time.h>
//#include <crtdbg.h>

#include "mpi.h"
#include <omp.h>
//#include "memory_debug.h"
#include "constants.h"
#include "particle.h"
#include "matrix3d.h"
#include "util.h"
#include "input.h"
#include "simulation.h"
#include "mpi_util.h"
#include "paths.h"

int main(int argc, char** argv) {
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_ALWAYS_DF);
	int size;
	int rank;
	MPI_Init(&argc, &argv);

	int dims[MPI_dim];
	int period[MPI_dim];
	period[0] = 0;
	period[1] = 1;
	period[2] = 1;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		std::string inputDir = inputDirectory;
		FILE* MPI_inputFile = fopen((inputDir + "MPI_input.dat").c_str(), "r");
		fscanf(MPI_inputFile, "%d %d %d", &dims[0], &dims[1], &dims[2]);
		fclose(MPI_inputFile);
		if (size != dims[0] * dims[1] * dims[2]) {
			printf("wrong MPI dimensions\n");
			MPI_Finalize();
			exit(0);
		}
	}

	MPI_Bcast(dims, MPI_dim, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Comm comm;
	MPI_Cart_create(MPI_COMM_WORLD, MPI_dim, dims, period, 0, &comm);

	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	//printf("hello, world!");
	//exit(0);
	if (rank == 0) {
		printf("start version 16.11.2017 21:00\n");
	}

	for (int i = 0; i < size; ++i) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (i == rank) {
			printf("start thread number %d\n", rank);
			fflush(stdout);
		}
	}
	//printf("additionalBinNumber = %d\n", additionalBinNumber);
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) printf("start\n");
	//fflush(stdout);

	if (rank == 0) printf("random initialize\n");

	srand(time(NULL) + rank);
	//srand(initialRandom);
	//srand(initialRandom + rank);


	bool startNew = true;
	bool writelog = false;

	if (startNew) {
		Simulation simulation;
		if (rank == 0) {
			printf("open input\n");
			//fflush(stdout);
			std::string inputDir = inputDirectory;
			FILE* inputFile = fopen((inputDir + "input.dat").c_str(), "r");
			printf("read input\n");
			//fflush(stdout);
			simulation = readInput(inputFile, comm);
			printf("close input\n");
			//fflush(stdout);
			fclose(inputFile);
			//simulation.rank = 0;
			sendInput(simulation, size);
		} else {
			//printf("recieve\n");
			simulation = recieveInput(comm);
			//simulation.rank = rank;
		}

		if (rank == 0) printf("simulate\n");
		//fflush(stdout);
		simulation.simulate();
		if (rank == 0) printf("end\n");
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

		simulation = readBackup(backupGeneralFileName, backupEfieldFileName, backupBfieldFileName, backupParticlesFileName,
		                        comm);


		if (rank == 0) printf("simulate\n");
		fflush(stdout);
		simulation.simulate();
		if (rank == 0) printf("end\n");
		fflush(stdout);
	}

	//std::string outDir = outputDirectory;
	//const char* memoryFileName = (outDir + "memory.dat").c_str();
	//const char* memoryFileName = ("C:/users/Vadim/Documents/Visual Studio 2010/Projects/v4/PIC++/output/memory.dat");
	//FILE* memoryFile = fopen((outDir + "memory.dat").c_str(),"w");
	//fclose(memoryFile);
	//_CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
	//_CrtSetReportFile( _CRT_WARN, memoryFile );
	//_CrtDumpMemoryLeaks();
	//fclose(memoryFile);
	MPI_Finalize();
}

