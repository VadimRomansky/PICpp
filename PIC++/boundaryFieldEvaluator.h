#ifndef boundaryFieldEvaluator
#define boundaryFieldEvaluator
#include "vector3d.h"

class Simulation;

class BoundaryFieldEvaluator {
public:
	BoundaryFieldEvaluator();
	virtual ~BoundaryFieldEvaluator();
	virtual Vector3d evaluateEfield(double t, int j, int k) = 0;
	virtual Vector3d evaluateBfield(double t, int j, int k) = 0;
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

	virtual Vector3d evaluateEfield(double t, int j, int k);
	virtual Vector3d evaluateBfield(double t, int j, int k);
};

class StripeBoundaryFieldEvaluator : BoundaryFieldEvaluator{
	Vector3d E;
	Vector3d B;
	double lambda;
	double ux;
	double shiftX;

	int halfPeriodNumber(double& t);
public:
	StripeBoundaryFieldEvaluator(Vector3d& E0, Vector3d& B0, double l, double u, double shift);

	virtual Vector3d evaluateEfield(double t, int j, int k);
	virtual Vector3d evaluateBfield(double t, int j, int k);
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

	virtual Vector3d evaluateEfield(double t, int j, int k);
	virtual Vector3d evaluateBfield(double t, int j, int k);
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
public:
	RandomTurbulenceBoundaryFieldEvaluator(int randomSeedV, int minLengthXV, int maxLengthXV, int minLengthYV, int maxLengthYV, int minLengthZV, int maxLengthZV, Simulation* simuationV, Vector3d V, Vector3d E, Vector3d B, double x, double dx, double dy, double dz);
	virtual Vector3d evaluateEfield(double t, int j, int k);
	virtual Vector3d evaluateBfield(double t, int j, int k);
};

#endif