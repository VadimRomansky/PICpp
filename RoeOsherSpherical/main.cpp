//#include <crtdbg.h>
#include <omp.h>
#include <time.h>
#include "constants.h"
#include "input.h"
#include "simulation.h"
#include "complex.h"



int main()
{
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_ALWAYS_DF);
    //omp_set_dynamic(0);      // запретить библиотеке openmp менять число потоков во время исполнения
    omp_set_num_threads(numThreads); // установить число потоков
	Simulation* simulation;
	//bool deserializeParameter = true;
	bool deserializeParameter = false;
	if(deserializeParameter){
		simulation = new Simulation();
		FILE* hydroFile = fopen("./save/hydro.dat", "r");
		FILE* distributionFile = fopen("./save/distribution.dat", "r");
		FILE* fieldFile = fopen("./save/field.dat", "r");
		FILE* gridFile = fopen("./save/grid.dat", "r");
		FILE* pgridFile = fopen("./save/pgrid.dat", "r");
		FILE* kgridFile = fopen("./save/kgrid.dat", "r");
		FILE* infoFile = fopen("./save/info.dat", "r");

		deserialize(hydroFile, distributionFile, fieldFile, gridFile, pgridFile, kgridFile, infoFile, simulation);

		fclose(hydroFile);
		fclose(distributionFile);
		fclose(fieldFile);
		fclose(gridFile);
		fclose(pgridFile);
		fclose(kgridFile);
		fclose(infoFile);
	} else {
		simulation = readInput();
	}
	simulation->simulate();
	return 0;
}
