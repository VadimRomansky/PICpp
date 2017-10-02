#ifndef SIMULATION_H
#define SIMULATION_H

#include "mpi.h"

#include "constants.h";

class Complex;

class Simulation{
public:
	int rank;
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Comm cartComm;
	int nprocs;
	int leftRank;
	int rightRank;
	int frontRank;
	int backRank;
	int bottomRank;
	int topRank;

	int xnumber;
	int ynumber;
	int znumber;

	double deltaX;
	double deltaY;
	double deltaZ;

	double deltaKx;
	double deltaKy;
	double deltaKz;

	double* xgrid;
	double* ygrid;
	double* zgrid;

	double* middleXgrid;
	double* middleYgrid;
	double* middleZgrid;

	double* kxgrid;
	double* kygrid;
	double* kzgrid;

	int* xabsoluteIndex;
	int* yabsoluteIndex;
	int* zabsoluteIndex;

	int firstXabsoluteIndex;
	int firstYabsoluteIndex;
	int firstZabsoluteIndex;

	Complex*** function;
	Complex*** fourierImage;
	Complex*** result;


	Simulation(MPI_Comm& comm);
	~Simulation();
	void setSpaceForProc();
	void createArrays();
	void initializeArrays();
	void createFiles();

	void randomSimulation();
	void poissonSolving();

	Complex*** create3dArray();
	void delete3dArray(Complex*** array);
};

#endif