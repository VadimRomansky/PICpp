#ifndef boundaryFieldEvaluator
#define boundaryFieldEvaluator
#include "vector3d.h"

class BaseBoundaryFieldEvaluator {
public:
	BaseBoundaryFieldEvaluator();
	~BaseBoundaryFieldEvaluator();
	virtual Vector3d evaluateEfield(double& t) = 0;
	virtual Vector3d evaluateBfield(double& t) = 0;
};

class ConstantBoundaryFieldEvaluator : BaseBoundaryFieldEvaluator{
	Vector3d E;
	Vector3d B;
public:
	ConstantBoundaryFieldEvaluator(Vector3d& E0, Vector3d& B0);

	virtual Vector3d evaluateEfield(double& t);
	virtual Vector3d evaluateBfield(double& t);
};

class StripeBoundaryFieldEvaluator : BaseBoundaryFieldEvaluator{
	Vector3d E;
	Vector3d B;
	double lambda;
	double ux;
	double shiftX;

	int halfPeriodNumber(double& t);
public:
	StripeBoundaryFieldEvaluator(Vector3d& E0, Vector3d& B0, double l, double u, double shift);

	virtual Vector3d evaluateEfield(double& t);
	virtual Vector3d evaluateBfield(double& t);
};

#endif