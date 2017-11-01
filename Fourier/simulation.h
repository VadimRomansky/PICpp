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

	int xnumberAdded;
	int ynumberAdded;
	int znumberAdded;

	double deltaX;
	double deltaY;
	double deltaZ;

	double xsize;
	double ysize;
	double zsize;

	double leftX;
	double rightX;
	double leftY;
	double rightY;
	double leftZ;
	double rightZ;

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

	int firstAbsoluteXindex;
	int firstAbsoluteYindex;
	int firstAbsoluteZindex;

	Complex*** function;
	Complex*** fourierImage;
	Complex*** result;

	Complex* function1d;
	Complex* fourierImage1d;
	Complex* result1d;


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