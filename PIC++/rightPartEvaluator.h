#ifndef RIGHT_PART_EVALUATOR_H
#define RIGHT_PART_EVALUATOR_H

#include "vector3d.h"
#include "particle.h"

class BaseRightPartEvaluator{
public:
	virtual ~BaseRightPartEvaluator() {
	}
	int rank;
	int nprocs;
	int xnumber;
	int ynumber;
	int znumber;
	int lnumber;
    BaseRightPartEvaluator(int xn, int yn, int zn, int ln);
    virtual double rightPart(double**** vector, int i, int j, int k, int l);
};

class InitializeShockRightPartEvaluator : public BaseRightPartEvaluator {
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
	InitializeShockRightPartEvaluator(Vector3d v1, Vector3d v2, Vector3d B1, Vector3d B2, Vector3d E1, double dx, double cvalue, double fieldNorm, double velocityNorm, ParticleTypeContainer* types, int t, int xn, int yn, int zn, int ln);
	~InitializeShockRightPartEvaluator();
	double evaluateVx(double**** vector, int i, int l);
	double evaluateVy(double**** vector, int i, int l);
	double evaluateVz(double**** vector, int i, int l);
	double evaluateEx(double**** vector, int i);
	double evaluateBy(double**** vector, int i);
	double evaluateBz(double**** vector, int i);
	virtual double rightPart(double**** vector, int i, int j, int k, int l);
};
#endif