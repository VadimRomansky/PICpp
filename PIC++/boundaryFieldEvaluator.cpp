#include "stdlib.h"
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "boundaryFieldEvaluator.h"
#include "paths.h"
#include <cmath>
#include "rightPartEvaluator.h"
#include "random.h"
#include "simulation.h"

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

TurbulenceBoundaryFieldEvaluator::TurbulenceBoundaryFieldEvaluator(Vector3d& E, Vector3d& B, Vector3d& V, int numberv, double* amplitudev,
                                                                   double* phasev, double* kv, double* omegav, double xv, double c) {
	E0 = E;
	B0 = B;
	V0 = V;
	speed_of_light_normalized = c;
	x = xv;
	number = numberv;
	phase = new double[2 * number];
	amplitude = new double[2 * number];
	kw = new double[number];
	omega = new double[number];
	for (int i = 0; i < number; ++i) {
		kw[i] = kv[i];
		omega[i] = omegav[i];
	}
	for (int i = 0; i < 2 * number; ++i) {
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
	for (int i = 0; i < number; ++i) {
		result.y += amplitude[2 * i] * sin(kw[i] * (x - V0.x * t) - omega[i] * t + phase[2 * i]);
		result.z += amplitude[2 * i + 1] * sin(kw[i] * (x - V0.x * t) - omega[i] * t + phase[2 * i + 1]);
	}
	return result;
}

Vector3d TurbulenceBoundaryFieldEvaluator::evaluateEfield(double t, int j, int k) {
	Vector3d result = E0;
	Vector3d B = evaluateBfield(t, j, k) - B0;
	result = result - V0.vectorMult(B) / (speed_of_light_normalized * speed_of_light_correction);
	return result;
}

RandomTurbulenceBoundaryFieldEvaluator::RandomTurbulenceBoundaryFieldEvaluator(int randomSeedV, int minLengthXV, int maxLengthXV, int minLengthYV,
                                                                               int maxLengthYV, int minLengthZV, int maxLengthZV,
                                                                               Simulation* simuationV, Vector3d V0, Vector3d E, Vector3d B, double xv,
                                                                               double dx, double dy, double dz) {
	randomSeed = randomSeedV;
	minLengthX = minLengthXV;
	maxLengthX = maxLengthXV;
	minLengthY = minLengthYV;
	maxLengthY = maxLengthYV;
	minLengthZ = minLengthZV;
	maxLengthZ = maxLengthZV;
	simulation = simuationV;
	V = V0;
	E0 = E;
	B0 = B;
	x = xv;
	deltaX = dx;
	deltaY = dy;
	deltaZ = dz;
}

Vector3d RandomTurbulenceBoundaryFieldEvaluator::evaluateBfield(double t, int j, int k) {
	srand(randomSeed);
	int maxKxnumber = maxLengthX / minLengthX;
	int maxKynumber = maxLengthY / minLengthY;
	int maxKznumber = maxLengthZ / minLengthZ;
	Vector3d result = B0;
	for (int ki = 0; ki < maxKxnumber; ki++) {
		for (int kj = 0; kj < maxKynumber; ++kj) {
			for (int kk = 0; kk < maxKznumber; ++kk) {
				if (ki + kj + kk > 0) {
					double kx = 2 * pi * ki / (deltaX * maxLengthX);
					double ky = 2 * pi * kj / (deltaY * maxLengthY);
					double kz = 2 * pi * kk / (deltaZ * maxLengthZ);
					double phase1 = 2 * pi * uniformDistribution();
					double phase2 = 2 * pi * uniformDistribution();

					double kw = sqrt(kx * kx + ky * ky + kz * kz);
					double kyz = sqrt(ky * ky + kz * kz);
					double cosTheta = kx / kw;
					double sinTheta = kyz / kw;
					double cosPhi;
					double sinPhi;
					if (kk + kj > 0) {
						cosPhi = ky / kyz;
						sinPhi = kz / kyz;
					} else {
						cosPhi = 1.0;
						sinPhi = 0.0;
					}
					double Bturbulent = simulation->evaluateTurbulenceFieldAmplitude(kx, ky, kz);
					double kmultr = kx * (x - 0.5 * deltaX - V.x * t) + ky * (simulation->middleYgrid[j] - V.y * t) + kz * (simulation->middleZgrid[k] - V.z * t);
					double localB1 = Bturbulent * sin(kmultr + phase1);
					double localB2 = Bturbulent * sin(kmultr + phase2);
					Vector3d localB = Vector3d(-sinTheta * localB1, cosTheta * cosPhi * localB1 - sinTheta * localB2,
					                           cosTheta * sinPhi * localB1 - cosTheta * localB2);
					result = result + localB;
				}
			}
		}
	}
	return result;
}

Vector3d RandomTurbulenceBoundaryFieldEvaluator::evaluateEfield(double t, int j, int k) {
	srand(randomSeed);
	int maxKxnumber = maxLengthX / minLengthX;
	int maxKynumber = maxLengthY / minLengthY;
	int maxKznumber = maxLengthZ / minLengthZ;
	Vector3d result = E0;
	for (int ki = 0; ki < maxKxnumber; ki++) {
		for (int kj = 0; kj < maxKynumber; ++kj) {
			for (int kk = 0; kk < maxKznumber; ++kk) {
				if (ki + kj + kk > 0) {
					double kx = 2 * pi * ki / (deltaX * maxLengthX);
					double ky = 2 * pi * kj / (deltaY * maxLengthY);
					double kz = 2 * pi * kk / (deltaZ * maxLengthZ);
					double phase1 = 2 * pi * uniformDistribution();
					double phase2 = 2 * pi * uniformDistribution();

					double kw = sqrt(kx * kx + ky * ky + kz * kz);
					double kyz = sqrt(ky * ky + kz * kz);
					double cosTheta = kx / kw;
					double sinTheta = kyz / kw;
					double cosPhi;
					double sinPhi;
					if (kk + kj > 0) {
						cosPhi = ky / kyz;
						sinPhi = kz / kyz;
					} else {
						cosPhi = 1.0;
						sinPhi = 0.0;
					}

					double Bturbulent = simulation->evaluateTurbulenceFieldAmplitude(kx, ky, kz);

					double kmultr = kx * (x - V.x * t) + ky * (simulation->ygrid[j] - V.y * t) + kz * (simulation->zgrid[k] - V.z * t);
					double localB1 = Bturbulent * sin(kmultr + phase1);
					double localB2 = Bturbulent * sin(kmultr + phase2);
					Vector3d localB = Vector3d(-sinTheta * localB1, cosTheta * cosPhi * localB1 - sinTheta * localB2,
					                           cosTheta * sinPhi * localB1 - cosTheta * localB2);
					result = result - V.vectorMult(localB) / (simulation->speed_of_light_normalized);
				}
			}
		}
	}
	return result;
}













