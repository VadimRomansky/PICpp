#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string>
#include <mpi.h>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "util.h"
#include "constants.h"
#include "random.h"
#include "paths.h"

double uniformDistribution() {
	return (rand() % randomParameter + 0.5) / randomParameter;
}

double normalDistribution() {
	//Box–Muller transform
	double x = uniformDistribution();
	double y = uniformDistribution();
	return cos(2 * pi * x) * sqrt(-2 * log(y));
}

double maxwellDistribution(double temperature, double k) {
	double normal = normalDistribution();
	return k * temperature * (0.5 * normal * normal - log(uniformDistribution()));
}

double maxwellJuttnerDistribution(double temperature, double mass, double c, double k) {
	double theta = k * temperature / (mass * c * c);
	double besselK = McDonaldFunction(1.0 / theta, 2.0);

	double x = uniformDistribution();

	double gamma = solveInverceJuttnerFunction(x, theta, besselK);

	return gamma * mass * c * c;
}

double solveInverceJuttnerFunction(double x, double theta, double besselK) {
	if (x >= 1) {
		printf("distribution function can not be more than 1\n");
		std::string outputDir = outputDirectory;
		//FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "distribution function can not be more than 1\n");
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (x <= 0) {
		return 1;
	}
	double gamma = max2(1, theta);
	double integral = maxwellJuttnerIntegral(gamma, theta, besselK);
	if (integral < x) {
		while (integral < x) {
			gamma *= 2;
			integral = maxwellJuttnerIntegral(gamma, theta, besselK);
		}
	}
	return solveInverceJuttnerFunction(x, theta, besselK, 0.0, gamma);
}

double solveInverceJuttnerFunction(double x, double theta, double besselK, double left, double right) {
	if (right < left) {
		printf("right < left\n");
		std::string outputDir = outputDirectory;
		FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "right < left in solveInverceJuttnerFunction\n");
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	if (right - left < 0.000001 * left) {
		return (right + left) / 2;
	}
	double gamma = (left + right) / 2;

	double integral = maxwellJuttnerIntegral(gamma, theta, besselK);

	if (integral > x) {
		return solveInverceJuttnerFunction(x, theta, besselK, left, gamma);
	}
	return solveInverceJuttnerFunction(x, theta, besselK, gamma, right);
}

double maxwellJuttnerFunction(double gamma, double theta, double besselK) {
	double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
	return gamma * gamma * beta * exp(-gamma / theta) / (theta * besselK);
}

double maxwellJuttnerIntegral(double gamma, double theta, double besselK) {
	if (gamma < 1) {
		return 0;
	}
	int count = 1000;
	double dgamma = (gamma - 1) / count;
	double sum = 0;
	double tempGamma = 1 + dgamma / 2;
	for (int i = 0; i < count; ++i) {
		sum += maxwellJuttnerFunction(tempGamma, theta, besselK) * dgamma;
		tempGamma += dgamma;
	}
	return sum;
}

void anisotropicMaxwellJuttnerDistribution(double& momentumNormal, double& momentumParallel, double alphaNormal,
                                           double alphaParallel, double m_c) {

	double maxMomentumNormal = 10.0 / alphaNormal;
	double maxMomentumParallel = 10.0 / alphaParallel;

	double b = 2 * alphaNormal * alphaNormal * alphaParallel * (2 * alphaNormal + alphaParallel) / sqr(
		alphaNormal - alphaParallel);

	double extremumP2PointSqr = ((4 * sqr(alphaNormal) + 2 * b) / (b * b) - 1);
	double extremumP2Point = 0;
	if (extremumP2PointSqr && (alphaNormal > alphaParallel) > 0) {
		extremumP2Point = sqrt(extremumP2PointSqr);
	}

	double extremumP1Point = sqrt(
		(1 + sqrt(1 + 4 * sqr(alphaNormal) * (sqr(extremumP2Point) + 1))) / (2 * sqr(alphaNormal)));

	double maxFunctionValue = anisotropicMaxwellJuttnerFunction(extremumP1Point, extremumP2Point, alphaNormal,
	                                                            alphaParallel);

	double tempP1 = uniformDistribution() * maxMomentumNormal;
	double tempP2 = uniformDistribution() * maxMomentumParallel;

	double functionValue = anisotropicMaxwellJuttnerFunction(tempP1, tempP2, alphaNormal, alphaParallel);

	if (functionValue > maxFunctionValue) {
		printf("functionValue > maxFunctionValue\n");
		std::string outputDir = outputDirectory;
		FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "functionValue > maxFunctionValue in anisotropicMaxwellJuttnerDistribution\n");
		fprintf(errorLogFile, "functionValue = %g maxValue = %g p1 = %g p2 = %g alpha1 = %g alpha2 = %g\n", functionValue,
		        maxFunctionValue, tempP1, tempP2, alphaNormal, alphaParallel);
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}

	double tempValue = uniformDistribution() * maxFunctionValue;

	while (tempValue > functionValue) {
		tempP1 = uniformDistribution() * maxMomentumNormal;
		tempP2 = uniformDistribution() * maxMomentumParallel;

		functionValue = anisotropicMaxwellJuttnerFunction(tempP1, tempP2, alphaNormal, alphaParallel);

		if (functionValue > maxFunctionValue) {
			printf("functionValue > maxFunctionValue\n");
			std::string outputDir = outputDirectory;
			FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "functionValue > maxFunctionValue in anisotropicMaxwellJuttnerDistribution\n");
			fprintf(errorLogFile, "functionValue = %g maxValue = %g p1 = %g p2 = %g alpha1 = %g alpha2 = %g\n", functionValue,
			        maxFunctionValue, tempP1, tempP2, alphaNormal, alphaParallel);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}

		tempValue = uniformDistribution() * maxFunctionValue;
	}

	momentumNormal = tempP1 * m_c;
	momentumParallel = tempP2 * m_c;
}

double anisotropicMaxwellJuttnerFunction(double p1, double p2, double alpha1, double alpha2) {
	return p1 * exp(-alpha1 * (sqrt(1 + p1 * p1 + p2 * p2) - sqrt(1 + p2 * p2)) - alpha2 * sqrt(1 + p2 * p2));
}
