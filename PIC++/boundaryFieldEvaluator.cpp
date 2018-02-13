#include "stdlib.h"
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "boundaryFieldEvaluator.h"
#include "paths.h"
#include <cmath>
#include "rightPartEvaluator.h"

BoundaryFieldEvaluator::BoundaryFieldEvaluator() {

}

BoundaryFieldEvaluator::~BoundaryFieldEvaluator() {

}

ConstantBoundaryFieldEvaluator::ConstantBoundaryFieldEvaluator(Vector3d& E0, Vector3d& B0) {
	E = E0;
	B = B0;
}

Vector3d ConstantBoundaryFieldEvaluator::evaluateEfield(double t, int j, int k) {
	return E;
}

Vector3d ConstantBoundaryFieldEvaluator::evaluateBfield(double t, int j, int k) {
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

Vector3d StripeBoundaryFieldEvaluator::evaluateEfield(double t, int j, int k) {
	int n = halfPeriodNumber(t);
	return E * n;
}

Vector3d StripeBoundaryFieldEvaluator::evaluateBfield(double t, int j, int k) {
	int n = halfPeriodNumber(t);
	return B * n;
}

TurbulenceBoundaryFieldEvaluator::TurbulenceBoundaryFieldEvaluator(Vector3d& E, Vector3d& B, Vector3d& V, int numberv, double* amplitudev, double* phasev, double* kv, double* omegav, double xv, double c) {
	E0 = E;
	B0 = B;
	V0 = V;
	speed_of_light_normalized = c;
	x = xv;
	number = numberv;
	phase = new double[2*number];
	amplitude = new double[2*number];
	kw = new double[number];
	omega = new double[number];
	for(int i = 0; i < number; ++i) {
		kw[i] = kv[i];
		omega[i] = omegav[i];
	}
	for(int i = 0; i < 2*number; ++i) {
		phase[i] = phasev[i];
		amplitude[i] = amplitudev[i];
	}
}

TurbulenceBoundaryFieldEvaluator::~TurbulenceBoundaryFieldEvaluator() {
	delete[] phase;
	delete[] amplitude;
	delete[] kw;
	delete[] omega;
}

Vector3d TurbulenceBoundaryFieldEvaluator::evaluateBfield(double t, int j, int k) {
	Vector3d result = B0;
	for(int i = 0; i < number; ++i) {
		result.y += amplitude[2*i]*sin(kw[i]*(x - V0.x*t) - omega[i]*t + phase[2*i]);
		result.z += amplitude[2*i + 1]*sin(kw[i]*(x - V0.x*t) - omega[i]*t + phase[2*i + 1]);
	}
	return result;
}

Vector3d TurbulenceBoundaryFieldEvaluator::evaluateEfield(double t, int j, int k) {
	Vector3d result = E0;
	Vector3d B = evaluateBfield(t, j, k);
	result = result - V0.vectorMult(B)/(speed_of_light_normalized*speed_of_light_correction);
	return result;
}










