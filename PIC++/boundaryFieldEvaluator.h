#ifndef boundaryFieldEvaluator
#define boundaryFieldEvaluator
#include "vector3d.h"

class Simulation;

class BoundaryFieldEvaluator {
public:
	BoundaryFieldEvaluator();
	virtual ~BoundaryFieldEvaluator();
	virtual Vector3d evaluateEfield(const double& t, int j, int k) = 0;
	virtual Vector3d evaluateBfield(const double& t, int j, int k) = 0;
	virtual void prepareB(const double& t);
	virtual void prepareE(const double& t);
};

/*Vector3d BoundaryFieldEvaluator::evaluateEfield(double& t, int j, int k) {
	return Vector3d(0, 0, 0);
}

Vector3d BoundaryFieldEvaluator::evaluateBfield(double& t, int j, int k) {
	return Vector3d(0, 0, 0);
}*/


class ConstantBoundaryFieldEvaluator : public BoundaryFieldEvaluator{
	Vector3d E;
	Vector3d B;
public:
	ConstantBoundaryFieldEvaluator(Vector3d& E0, Vector3d& B0);

	virtual Vector3d evaluateEfield(const double& t, int j, int k);
	virtual Vector3d evaluateBfield(const double& t, int j, int k);
};

class StripeBoundaryFieldEvaluator : BoundaryFieldEvaluator{
	Vector3d E;
	Vector3d B;
	double lambda;
	double ux;
	double shiftX;

	int halfPeriodNumber(const double& t);
public:
	StripeBoundaryFieldEvaluator(Vector3d& E0, Vector3d& B0, double l, double u, double shift);

	virtual Vector3d evaluateEfield(const double& t, int j, int k);
	virtual Vector3d evaluateBfield(const double& t, int j, int k);
};

class TurbulenceBoundaryFieldEvaluator : public  BoundaryFieldEvaluator {
	Vector3d E0;
	Vector3d B0;
	Vector3d V0;
	double x;
	double speed_of_light_normalized;
	int number;
	double* amplitude;
	double* phase;
	double* omega;
	double* kw;
public:
	TurbulenceBoundaryFieldEvaluator(Vector3d& E0, Vector3d& B0, Vector3d& V0, int number, double* amplitude, double* phase, double* kv, double* omega, double xv, double c);
	~TurbulenceBoundaryFieldEvaluator();

	virtual Vector3d evaluateEfield(const double& t, int j, int k);
	virtual Vector3d evaluateBfield(const double& t, int j, int k);
};

class RandomTurbulenceBoundaryFieldEvaluator : public BoundaryFieldEvaluator {
	int randomSeed;
	int minLengthX, maxLengthX;
	int minLengthY, maxLengthY;
	int minLengthZ, maxLengthZ;
	Simulation* simulation;
	Vector3d V;
	Vector3d E0;
	Vector3d B0;
	double x;
	double deltaX;
	double deltaY;
	double deltaZ;
	int xnumberGeneral;
	int ynumberGeneral;
	int znumberGeneral;
	int xnumberAdded;
	int ynumberAdded;
	int znumberAdded;
	Vector3d** Bfield;
	Vector3d** Efield;
public:
	RandomTurbulenceBoundaryFieldEvaluator(int randomSeedV, int minLengthXV, int maxLengthXV, int minLengthYV, int maxLengthYV, int minLengthZV, int maxLengthZV, Simulation* simuationV, Vector3d V, Vector3d E, Vector3d B, double x, double dx, double dy, double dz, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, int xnumberAdded, int ynumberAdded, int znumberAdded);
	virtual ~RandomTurbulenceBoundaryFieldEvaluator();
	virtual Vector3d evaluateEfield(const double& t, int j, int k);
	virtual Vector3d evaluateBfield(const double& t, int j, int k);
	virtual void prepareB(const double& t);
	virtual void prepareE(const double& t);
};

#endif