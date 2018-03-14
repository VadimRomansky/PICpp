#ifndef RIGHT_PART_EVALUATOR_H
#define RIGHT_PART_EVALUATOR_H

#include <mpi.h>

#include "vector3d.h"
#include "particle.h"

class RightPartIterativeEvaluator{
public:
	virtual ~RightPartIterativeEvaluator() {}
	int rank;
	int nprocs;
	int xnumberAdded;
	int ynumberAdded;
	int znumberAdded;
	int* cartCoord;
	int* cartDim;
	MPI_Comm cartComm;
	int lnumber;
    RightPartIterativeEvaluator(int xn, int yn, int zn, int ln, int rank, int nprocs, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
    virtual double rightPart(double**** vector, int i, int j, int k, int l) = 0;
    virtual double rightPartInitialNorm() = 0;
    void getError(double**** vector, double**** errorVector);
};

class InitializeShockRightPartEvaluator : public RightPartIterativeEvaluator {
public:
	Vector3d downstreamV;
	Vector3d upstreamV;
	Vector3d downstreamB;
	Vector3d upstreamB;
	Vector3d E;
	double deltaX;
	double c;
	double fieldNorm;
	double velocityNorm;
	double* charge;
	double* mass;
	double* flux;
	int typesNumber;
	InitializeShockRightPartEvaluator(Vector3d v1, Vector3d v2, Vector3d B1, Vector3d B2, Vector3d E1, double dx, double cvalue, double fieldNorm, double velocityNorm, ParticleTypeContainer* types, int t, int xn, int yn, int zn, int ln, int rank, int nprocs, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
	~InitializeShockRightPartEvaluator();
	double evaluateVx(double**** vector, int i, int l);
	double evaluateVy(double**** vector, int i, int l);
	double evaluateVz(double**** vector, int i, int l);
	double evaluateEx(double**** vector, int i);
	double evaluateBy(double**** vector, int i);
	double evaluateBz(double**** vector, int i);
	virtual double rightPart(double**** vector, int i, int j, int k, int l);
	virtual double rightPartInitialNorm();
};

class PoissonRightPartEvaluator : public RightPartIterativeEvaluator {
public:
	PoissonRightPartEvaluator(double*** densityv, double**** tempLargeVector, double dx, double dy, double dz, int xn, int yn, int zn, bool periodicXv, bool periodicYv, bool periodicZv, int rankv, int nprocsv,
	                                         MPI_Comm& cartCommv, int* cartCoordv, int* cartDimv);

	double deltaX;
	double deltaY;
	double deltaZ;

	double**** density;
	bool periodicX;
	bool periodicY;
	bool periodicZ;

	virtual double rightPart(double**** vector, int i, int j, int k, int l);
	virtual double rightPartInitialNorm();
};
#endif