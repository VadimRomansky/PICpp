#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include "constants.h"

#include "distribution.h"

void createParticle(double& px, double& py, double& pz, double temperature, double mass, double* juttnerValue, double* juttnerFunction, int juttnerN){
	double p;
	double thetaParameter = kBoltzman * temperature / (mass*c*c);

	if (thetaParameter < 0.1) {
		//printf("create maxwell particle\n");
		//energy = maxwellDistribution(localTemparature, kBoltzman_normalized);
		px = sqrt(mass * kBoltzman * temperature) * normalDistribution();
		py = sqrt(mass * kBoltzman * temperature) * normalDistribution();
		pz = sqrt(mass * kBoltzman * temperature) * normalDistribution();
		p = sqrt(px * px + py * py + pz * pz);
	} else {
		if (thetaParameter < 2) {
			//printf("create cold juttner particle\n");
			//p = maxwellJuttnerMomentumColdDistribution(temperature, mass, speed_of_light_normalized, kBoltzman_normalized,typeContainer.juttnerFunction, typeContainer.juttnerValue,typeContainer.jutnerNumber);
			p = maxwellJuttnerMomentumColdDistribution(temperature, mass, kBoltzman, juttnerFunction, juttnerValue, juttnerN);
		} else {
			//printf("create hot juttner particle\n");
			//p = maxwellJuttnerMomentumHotDistribution(temperature, mass, speed_of_light_normalized, kBoltzman_normalized);
			p = maxwellJuttnerMomentumHotDistribution(temperature, mass, kBoltzman);
		}


		//p = 0;

		pz = p * (2 * uniformDistribution() - 1);
		double phi = 2 * pi * uniformDistribution();
		double pnormal = sqrt(p * p - pz * pz);
		px = pnormal * cos(phi);
		py = pnormal * sin(phi);
	}
}

double uniformDistribution() {
	return (rand() % randomParameter + 0.5) / randomParameter;
}

double normalDistribution() {
	//Box–Muller transform
	double x = uniformDistribution();
	double y = uniformDistribution();
	return cos(2 * pi * x) * sqrt(-2 * log(y));
}

double maxwellDistribution(const double& temperature, const double& k) {
	double normal = normalDistribution();
	return k * temperature * (0.5 * normal * normal - log(uniformDistribution()));
}

double maxwellJuttnerMomentumColdDistribution(const double& temperature, const double& mass, const double& k, const double* function, const double* xvalue, int number) {
	double theta = k * temperature / (mass * c * c);
	double y = uniformDistribution();
	double x = dichotomySolver(function, 0, number - 1, xvalue, y);
	return mass*c*x;
}

double maxwellJuttnerMomentumHotDistribution(const double& temperature, const double& mass, const double& k) {
	double theta = k * temperature / (mass * c * c);
	double x1 = uniformDistribution();
	double x2 = uniformDistribution();
	double x3 = uniformDistribution();
	double x4 = uniformDistribution();
	double u = - theta*log(x2*x3*x4);
	double exp1 = exp((u - sqrt(1 + u*u))/theta);
	int iterationCount = 0;
	while(exp1 < x1) {
		x2 = uniformDistribution();
		x3 = uniformDistribution();
		x4 = uniformDistribution();
		u = - theta*log(x2*x3*x4);
		exp1 = exp((u - sqrt(1 + u*u))/theta);
		if(iterationCount > 100000){
			printf("infinite cycle in juttner hot distribution\n");
			return mass*c*u;
		}
		iterationCount++;
	}
	return mass*c*u;
}

double dichotomySolver(const double* functionValue, int minIndex, int maxIndex, const double* xValue, double y) {
	if(minIndex + 1 == maxIndex) {
		double y1 = functionValue[minIndex];
		double y2 = functionValue[maxIndex];
		double x1 = xValue[minIndex];
		double x2 = xValue[maxIndex];
		return (x2 - x1)*(y - y1)/(y2 - y1) + x1;
	}
	int middleIndex = (minIndex + maxIndex)/2;
	if(functionValue[middleIndex] > y) {
		return dichotomySolver(functionValue, minIndex, middleIndex, xValue, y);
	} else if(functionValue[middleIndex] < y) {
		return dichotomySolver(functionValue, middleIndex, maxIndex, xValue, y);
	} else {
		return xValue[middleIndex];
	}
}

double juttnerDistribution(const double& u, const double& theta) {
	return exp(-sqrt(1 + u*u)/theta)*u*u;
}

void evaluateJuttnerFunction(double* juttnerValue, double* juttnerFunction, double temperature, double mass, int juttnerNumber) {
	double theta1 = kBoltzman * temperature / (mass*c*c);
	double umax = 20 * theta1;
	double deltaU = umax / (juttnerNumber - 1);
	juttnerFunction[0] = 0;
	juttnerValue[0] = 0;
	int subIter = 10;
	double subDeltaU = deltaU / subIter;
	for (int i = 1; i < juttnerNumber; ++i) {
		juttnerValue[i] = juttnerValue[i - 1] + deltaU;
		juttnerFunction[i] = juttnerFunction[i - 1];
		double localU = juttnerValue[i - 1];
		for (int j = 0; j < subIter; ++j) {
			juttnerFunction[i] += juttnerDistribution(localU, theta1) * subDeltaU;
			localU += subDeltaU;
		}
	}
	for (int i = 0; i < juttnerNumber; ++i) {
		juttnerFunction[i] /= juttnerFunction[juttnerNumber - 1];
	}
	juttnerFunction[juttnerNumber - 1] = 1.0;
}