#include "stdio.h"
#include <stdlib.h>
#include <time.h>
//#include <crtdbg.h>
#include "mpi.h"

#include "constants.h"
#include "complex.h"
#include "simulation.h"

int main(int argc, char **argv){
	int size;
	int rank;
	MPI_Init(&argc, &argv);

	int dims[MPI_dim];
	int period[MPI_dim];
	period[0] = 0;
	period[1] = 1;
	period[2] = 1;
	dims[0] = MPI_nx;
	dims[1] = MPI_ny;
	dims[2] = MPI_nz;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm comm;
	MPI_Cart_create(MPI_COMM_WORLD, MPI_dim, dims, period, 0, &comm);

	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	Simulation* simulation = new Simulation(comm);
	//simulation->randomSimulation();
	simulation->poissonSolving();

	MPI_Finalize();
	return 0;
}