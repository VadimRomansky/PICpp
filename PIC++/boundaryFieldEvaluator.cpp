#include "stdlib.h"
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "boundaryFieldEvaluator.h"
#include "paths.h"

BaseBoundaryFieldEvaluator::BaseBoundaryFieldEvaluator() {

}

BaseBoundaryFieldEvaluator::~BaseBoundaryFieldEvaluator() {

}

ConstantBoundaryFieldEvaluator::ConstantBoundaryFieldEvaluator(Vector3d& E0, Vector3d& B0) {
	E = E0;
	B = B0;
}

Vector3d ConstantBoundaryFieldEvaluator::evaluateEfield(double& t) {
	return E;
}

Vector3d ConstantBoundaryFieldEvaluator::evaluateBfield(double& t) {
	return B;
}

int StripeBoundaryFieldEvaluator::halfPeriodNumber(double& t) {
	return ((t * ux * 2 + shiftX) / lambda);
}

StripeBoundaryFieldEvaluator::
StripeBoundaryFieldEvaluator(Vector3d& E0, Vector3d& B0, double l, double u, double shift) {
	E = E0;
	B = B0;
	lambda = 0;
	ux = ux;
	shiftX = shift;
}

Vector3d StripeBoundaryFieldEvaluator::evaluateEfield(double& t) {
	int n = halfPeriodNumber(t);
	return E * n;
}

Vector3d StripeBoundaryFieldEvaluator::evaluateBfield(double& t) {
	int n = halfPeriodNumber(t);
	return B * n;
}








