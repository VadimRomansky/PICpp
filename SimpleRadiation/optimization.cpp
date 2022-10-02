#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include <time.h>
#include <cmath>
#include <string>
#include <omp.h>

#include "startparameters.h"
#include "constants.h"
#include "spectrum.h"
#include "util.h"

#include "optimization.h"

double evaluateOptimizationFunction(double* nu, double* observedInu, double* observedError, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, const double& B, const double& R, const double& fraction, const double& epsilonB, double dopplerBeta){
	double* totalInu = new double[Nnu];
	//evaluateVolumeAndLength(area, length, rmax, fractionSize);
	double err = 0;
	double n = B*B/(4*pi*massProtonReal*speed_of_light2*epsilonB);
	evaluateAllEmissivityAndAbsorption(nu, Inu, Anu, Nnu, Ee, dFe, Np, Nd, Bn, sintheta, thetaIndex, concentrations, n, B, 1.0, dopplerBeta);
		//
		//evaluateSpectrum(nu[i], totalInu, Inu[i], Anu[i], area, length, Nnu, rfactor);
		//
		if(geometry == SPHERICAL){
			evaluateSpectrumSpherical(nu, totalInu, Inu, Anu, R, Nnu, 1.0, fraction);
		} else {
			evaluateSpectrumFlat(nu, totalInu, Inu, Anu, R, Nnu, 1.0, fraction);
		}
		for(int j = 0; j < Nnu; ++j){
			double err1 = 0;
			if(scale == LINEAR){
				err1 = sqr(totalInu[j] - observedInu[j])/sqr(observedError[j]);
			} else {
				//todo
				err1 = sqr(log(totalInu[j]) - log(observedInu[j]))/sqr(observedError[j]/observedInu[j]);
			}
			err = err + err1;
		}
	
	delete[] totalInu;
	
	return err;
}

double evaluateOptimizationFunctionBandR(const double* vector, double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, const double& fraction, const double& epsilonB, double dopplerBeta){
	double r = vector[1]*maxR;
	double B = vector[0]*maxB;

	double* observedError = new double[Nnu];
	for(int i = 0; i < Nnu; ++i){
		observedError[i] = 1.0;
	}
	
	double result = evaluateOptimizationFunction(nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, B, r, fraction, epsilonB, dopplerBeta);
	delete[] observedError;
	return result;
}

double evaluateOptimizationFunctionGeneral(const double* vector, double* nu, double* observedInu, double* observedError, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double dopplerBeta){
	double r = vector[1]*maxR;
	double B = vector[0]*maxB;
	double fraction = vector[2];
	double epsilonB = vector[3];
	if (epsilonB < 0) {
		printf("epsilonB < 0\n");
		exit(0);
	}

	
	double result = evaluateOptimizationFunction(nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, B, r, fraction, epsilonB, dopplerBeta);

	return result;
}

void findMinParametersBandR(double* vector, const double* grad, double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double& currentF, const double& fraction, const double& epsilonB, double dopplerBeta){
	const int Ngrad = 2;
	double minVector[Ngrad];
	double tempVector1[Ngrad];
	double tempVector2[Ngrad];
	minVector[0] = minB/maxB;
	minVector[1] = minR/maxR;

	double maxLambda = sqrt(1.0*Ngrad);
	for(int i = 0; i < Ngrad; ++i) {
		if(grad[i] > 0) {
			maxLambda = min(maxLambda, fabs((vector[i] - minVector[i])/grad[i]));
		}
		if(grad[i] < 0) {
			maxLambda = min(maxLambda, fabs((1.0 - vector[i])/grad[i]));
		}
	}

	double step = maxLambda/2;

	for(int i = 0; i < Ngrad; ++i) {
		tempVector1[i] = vector[i] - grad[i]*step;
		if(tempVector1[i] != tempVector1[i]) {
			printf("tempvector[i] != tempVector[i]\n");
			exit(0);
		}
	}

	double f = evaluateOptimizationFunctionBandR(tempVector1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);
	if(f > currentF){
		int count = 0;
		while (f > currentF && count < 20){
			count++;
			step = step/2;
			
			for(int i = 0; i < Ngrad; ++i) {
				tempVector1[i] = vector[i] - grad[i]*step;
				if(tempVector1[i] != tempVector1[i]) {
					printf("tempvector[i] != tempVector[i]\n");
					exit(0);
				}
			}

			f = evaluateOptimizationFunctionBandR(tempVector1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);
		}

		if(f < currentF){
			//todo!
			for(int i = 0; i < Ngrad; ++i) {
				vector[i] = tempVector1[i];
			}
			return;
		} else {
			return;
		}
	}

	f = evaluateOptimizationFunctionBandR(vector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);

	maxLambda = sqrt(1.0*Ngrad);
	for(int i = 0; i < Ngrad; ++i) {
		if(grad[i] > 0) {
			maxLambda = min(maxLambda, fabs((vector[i] - minVector[i])/grad[i]));
		}
		if(grad[i] < 0) {
			maxLambda = min(maxLambda, fabs((1.0 - vector[i])/grad[i]));
		}
	}

	for(int i = 0; i < Ngrad; ++i) {
		tempVector1[i] = vector[i] - maxLambda*grad[i];
		if(tempVector1[i] > 1.0) {
			tempVector1[i] = 1.0;
		}
		if(tempVector1[i] < minVector[i]) {
			tempVector1[i] = minVector[i];
		}
	}
	double f3 = evaluateOptimizationFunctionBandR(tempVector1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);

	double leftLambda = 0;
	double rightLambda = maxLambda;


	for(int j = 0; j < 10; ++j) {
		double lambda1 = leftLambda + (rightLambda - leftLambda)/3.0;
		double lambda2 = leftLambda + (rightLambda - leftLambda)*2.0/3.0;

		for(int i = 0; i < Ngrad; ++i) {
			tempVector1[i] = vector[i] - lambda1*grad[i];
			tempVector2[i] = vector[i] - lambda2*grad[i];
		}

		double f1 = evaluateOptimizationFunctionBandR(tempVector1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);
		double f2 = evaluateOptimizationFunctionBandR(tempVector2, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);

		if(f1 > f2) {
			leftLambda = lambda1;
		} else {
			rightLambda = lambda2;
		}
	}

	double lambda = (leftLambda + rightLambda)/2.0;
	for(int i = 0; i < Ngrad; ++i) {
		tempVector1[i] = vector[i] - lambda*grad[i];
		if(tempVector1[i] > 1.0) {
			tempVector1[i] = 1.0;
		}
		if(tempVector1[i] < minVector[i]) {
			tempVector1[i] = minVector[i];
		}
	}

	double tempF = evaluateOptimizationFunctionBandR(tempVector1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);

	if(f < tempF) {
		if(f < f3) {
			return;
		} else {
			for(int i = 0; i < Ngrad; ++i) {
				vector[i] = vector[i] - maxLambda*grad[i];
			}
		}
	} else {
		if(tempF < f3) {
			for(int i = 0; i < Ngrad; ++i) {
				vector[i] = vector[i] - lambda*grad[i];
			}
		} else {
			for(int i = 0; i < Ngrad; ++i) {
				vector[i] = vector[i] - maxLambda*grad[i];
			}
		}
	}
	for(int i  = 0; i < Ngrad; ++i){
		if(vector[i] > 1.0) {
			vector[i] = 1.0;
		}
		if(vector[i] < minVector[i]) {
			vector[i] = minVector[i];
		}
	}
}

void optimizeParametersBandR(double* vector,  double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, const double& fraction, const double& epsilonB, double dopplerBeta, FILE* logFile){
	const int Ngrad = 2;
	double minVector[Ngrad];
	double tempVector[Ngrad];
	minVector[0] = minB/maxB;
	minVector[1] = minR/maxR;
	double currentF = evaluateOptimizationFunctionBandR(vector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g rmax = %g\n", vector[0]*maxB, vector[1]*maxR);
	fprintf(logFile, "Bfactor = %grmax = %g\n", vector[0]*maxB, vector[1]*maxR);
	for(int k = 0; k < Niterations; ++k) {
		///randomization;
		/*for(int j = 0; j < 5; ++j){
			for(int i = 0; i < Ngrad; ++i) {
				//tempVector[i] = minVector[i] + (1.0 - minVector[i])*uniformDistribution();
				tempVector[i] = tempVector[i] + 0.2*minVector[i]*(0,5 - uniformDistribution());
				if(tempVector[i] > 1.0) {
					tempVector[i] = 1.0;
				}
				if(tempVector[i] < minVector[i]) {
					tempVector[i] = minVector[i];
				}
			}
			
			double tempF = evaluateOptimizationFunction5(tempVector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);
			if(tempF < currentF){
				currentF = tempF;
				for(int i = 0; i < Ngrad; ++i) {
					vector[i] = tempVector[i];
				}
				printf("random search\n");
				fprintf(logFile, "random search\n");
			}
		}*/
		double prevF = currentF;
		//
		//valley first step
		double valley1[Ngrad];
		for(int i = 0; i < Ngrad; ++i) {
			valley1[i] = vector[i];
		}

		double grad[Ngrad];

		printf("optimization k = %d\n",k);
		fprintf(logFile, "optimization k = %d\n",k);
		fflush(logFile);

		for(int i = 0; i < Ngrad; ++i) {
			for(int j = 0; j < Ngrad; ++j) {
				tempVector[j] = vector[j];
			}
			double dx = fabs(vector[i])/1000;
			tempVector[i] = vector[i] + dx;
			//currentF = evaluateOptimizationFunction5(vector[0]*maxB, vector[1]*maxN, vector[2]*maxFraction, vector[3]*maxR, vector[4]*maxV, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			double f = evaluateOptimizationFunctionBandR(tempVector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);
			tempVector[i] = vector[i] - dx;
			double f1 = evaluateOptimizationFunctionBandR(tempVector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);

			grad[i] = (f - f1)/(2*dx);
			if(grad[i] != grad[i]) {
				printf("grad[i] = NaN\n");
				exit(0);
			}
			if(grad[i] > 0) {
				if(vector[i] - minVector[i] < 0.000001) {
					grad[i] = 0;
				}
			}
			if(grad[i] < 0) {
				if(1.0 - vector[i] < 0.000001) {
					grad[i] = 0;
				}
			}
		}

		double gradNorm = 0;
		for(int i = 0; i < Ngrad; ++i){
			gradNorm = gradNorm + grad[i]*grad[i];
		}
		gradNorm = sqrt(gradNorm);
		if(gradNorm <= 0){
			printf("gradNorm <= 0\n");
			continue;
		}
		if(gradNorm != gradNorm){
			printf("gradNorm = NaN\n");
			exit(0);
		}
		for(int i = 0; i < Ngrad; ++i) {
			grad[i] = grad[i]/gradNorm;
		}


		findMinParametersBandR(vector, grad, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, fraction, epsilonB, dopplerBeta);

		currentF = evaluateOptimizationFunctionBandR(vector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);
		//valley second step

		for(int i = 0; i < Ngrad; ++i) {
			for(int j = 0; j < Ngrad; ++j) {
				tempVector[j] = vector[j];
			}
			double dx = fabs(vector[i])/1000;
			tempVector[i] = vector[i] + dx;
			double f = evaluateOptimizationFunctionBandR(tempVector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);
			tempVector[i] = vector[i] - dx;
			double f1 = evaluateOptimizationFunctionBandR(tempVector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);

			grad[i] = (f - f1)/(2*dx);
			if(grad[i] != grad[i]) {
				printf("grad[i] = NaN\n");
				exit(0);
			}
			if(grad[i] > 0) {
				if(vector[i] - minVector[i] < 0.000001) {
					grad[i] = 0;
				}
			}
			if(grad[i] < 0) {
				if(1.0 - vector[i] < 0.000001) {
					grad[i] = 0;
				}
			}
		}

 		gradNorm = 0;
		for(int i = 0; i < Ngrad; ++i){
			gradNorm = gradNorm + grad[i]*grad[i];
		}
		gradNorm = sqrt(gradNorm);
		if(gradNorm <= 0){
			printf("gradNorm <= 0\n");
			continue;
		}
		if(gradNorm != gradNorm){
			printf("gradNorm = NaN\n");
			exit(0);
		}
		for(int i = 0; i < Ngrad; ++i) {
			grad[i] = grad[i]/gradNorm;
		}

		findMinParametersBandR(vector, grad, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, fraction, epsilonB, dopplerBeta);

		currentF = evaluateOptimizationFunctionBandR(vector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);

		double valley2[Ngrad];
		for(int i = 0; i < Ngrad; ++i) {
			valley2[i] = vector[i];
		}

		//// valley third step

		for(int i = 0; i < Ngrad; ++i) {
			grad[i] = valley1[i] - valley2[i];
		}

		
		gradNorm = 0;
		for(int i = 0; i < Ngrad; ++i){
			gradNorm = gradNorm + grad[i]*grad[i];
		}
		gradNorm = sqrt(gradNorm);
		if(gradNorm <= 0){
			printf("gradNorm <= 0\n");
		}
		if(gradNorm != gradNorm){
			printf("gradNorm = NaN\n");
			exit(0);
		}
		
		if(gradNorm > 0){
			for(int i = 0; i < Ngrad; ++i) {
				grad[i] = grad[i]/gradNorm;
			}

			findMinParametersBandR(vector, grad, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, fraction, epsilonB, dopplerBeta);

		}
		currentF = evaluateOptimizationFunctionBandR(vector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB, dopplerBeta);
		
		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/

		printf("optimization function = %g\n", currentF);
		fprintf(logFile, "optimization function = %g\n", currentF);
		printf("Bfactor = %g rmax = %g\n", vector[0]*maxB, vector[1]*maxR);
		fprintf(logFile, "Bfactor = %g rmax = %g\n", vector[0]*maxB, vector[1]*maxR);
	}
	printf("finish optimization\n");
	fprintf(logFile, "finish optimization\n");
}

void findMinParametersGeneral(double* vector, bool* optPar, const double* grad, double* nu, double* observedInu, double* observedError, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double& currentF, double dopplerBeta){
	const int Ngrad = 4;
	int Npar = 0;
	for (int i = 0; i < Ngrad; ++i) {
		if (optPar[i]) {
			Npar++;
		}
	}
	double minVector[Ngrad];
	double tempVector1[Ngrad];
	double tempVector2[Ngrad];
	minVector[0] = minB/maxB;
	minVector[1] = minR/maxR;
	minVector[2] = minFraction;
	minVector[3] = minSigma;

	double maxLambda = sqrt(1.0*Npar);
	for(int i = 0; i < Ngrad; ++i) {
		if (optPar[i]) {
			if (grad[i] > 0) {
				maxLambda = min(maxLambda, fabs((vector[i] - minVector[i]) / grad[i]));
			}
			if (grad[i] < 0) {
				maxLambda = min(maxLambda, fabs((1.0 - vector[i]) / grad[i]));
			}
		}
	}

	double step = maxLambda/2;

	for(int i = 0; i < Ngrad; ++i){
		tempVector1[i] = vector[i];
	}
	for(int i = 0; i < Ngrad; ++i) {
		if (optPar[i]) {
			tempVector1[i] = vector[i] - grad[i] * step;
			if (tempVector1[i] != tempVector1[i]) {
				printf("tempvector[i] != tempVector[i]\n");
				exit(0);
			}
		}
	}

	double f = evaluateOptimizationFunctionGeneral(tempVector1, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
	if(f > currentF){
		int count = 0;
		while (f > currentF && count < 20){
			count++;
			step = step/2;
			
			for(int i = 0; i < Ngrad; ++i){
				tempVector1[i] = vector[i];
			}
			for(int i = 0; i < Ngrad; ++i) {
				if (optPar[i]) {
					tempVector1[i] = vector[i] - grad[i] * step;
					if (tempVector1[i] != tempVector1[i]) {
						printf("tempvector[i] != tempVector[i]\n");
						exit(0);
					}
				}
			}

			f = evaluateOptimizationFunctionGeneral(tempVector1, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
		}

		if(f < currentF){
			//todo!
			for(int i = 0; i < Ngrad; ++i) {
				vector[i] = tempVector1[i];
			}
			//return;
		} else {
			return;
		}
	}

	f = evaluateOptimizationFunctionGeneral(vector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);

	maxLambda = sqrt(1.0*Npar);
	for(int i = 0; i < Ngrad; ++i) {
		if (optPar[i]) {
			if (grad[i] > 0) {
				maxLambda = min(maxLambda, fabs((vector[i] - minVector[i]) / grad[i]));
			}
			if (grad[i] < 0) {
				maxLambda = min(maxLambda, fabs((1.0 - vector[i]) / grad[i]));
			}
		}
	}

	for(int i = 0; i < Ngrad; ++i){
		tempVector1[i] = vector[i];
	}
	for(int i = 0; i < Ngrad; ++i) {
		if (optPar[i]) {
			tempVector1[i] = vector[i] - maxLambda * grad[i];
			if (tempVector1[i] > 1.0) {
				tempVector1[i] = 1.0;
			}
			if (tempVector1[i] < minVector[i]) {
				tempVector1[i] = minVector[i];
			}
		}
	}
	double f3 = evaluateOptimizationFunctionGeneral(tempVector1, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);

	double leftLambda = 0;
	double rightLambda = maxLambda;


	for(int j = 0; j < 10; ++j) {
		double lambda1 = leftLambda + (rightLambda - leftLambda)/3.0;
		double lambda2 = leftLambda + (rightLambda - leftLambda)*2.0/3.0;

		for(int i = 0; i < Ngrad; ++i){
			tempVector1[i] = vector[i];
			tempVector2[i] = vector[i];
		}
		for(int i = 0; i < Ngrad; ++i) {
			if (optPar[i]) {
				tempVector1[i] = vector[i] - lambda1 * grad[i];
				tempVector2[i] = vector[i] - lambda2 * grad[i];
			}
		}

		double f1 = evaluateOptimizationFunctionGeneral(tempVector1, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
		double f2 = evaluateOptimizationFunctionGeneral(tempVector2, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);

		if(f1 > f2) {
			leftLambda = lambda1;
		} else {
			rightLambda = lambda2;
		}
	}

	double lambda = (leftLambda + rightLambda)/2.0;
	for(int i = 0; i < Ngrad; ++i){
		tempVector1[i] = vector[i];
	}
	for(int i = 0; i < Ngrad; ++i) {
		if (optPar[i]) {
			tempVector1[i] = vector[i] - lambda * grad[i];
			if (tempVector1[i] > 1.0) {
				tempVector1[i] = 1.0;
			}
			if (tempVector1[i] < minVector[i]) {
				tempVector1[i] = minVector[i];
			}
		}
	}

	double tempF = evaluateOptimizationFunctionGeneral(tempVector1, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);

	if(f < tempF) {
		if(f < f3) {
			return;
		} else {
			for(int i = 0; i < Ngrad; ++i) {
				if (optPar[i]) {
					vector[i] = vector[i] - maxLambda * grad[i];
				}
			}
		}
	} else {
		if(tempF < f3) {
			for(int i = 0; i < Ngrad; ++i) {
				if (optPar[i]) {
					vector[i] = vector[i] - lambda * grad[i];
				}
			}
		} else {
			for(int i = 0; i < Ngrad; ++i) {
				if (optPar[i]) {
					vector[i] = vector[i] - maxLambda * grad[i];
				}
			}
		}
	}
	for(int i  = 0; i < Ngrad; ++i){
		if (optPar[i]) {
			if (vector[i] > 1.0) {
				vector[i] = 1.0;
			}
			if (vector[i] < minVector[i]) {
				vector[i] = minVector[i];
			}
		}
	}
}

void optimizeParametersGeneral(double* vector, bool* optPar,  double* nu, double* observedInu, double* observedError, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double dopplerBeta, FILE* logFile){
	const int Ngrad = 4;
	int Npar = 0;
	for (int i = 0; i < Ngrad; ++i) {
		if (optPar[i]) {
			Npar++;
		}
	}
	if(Npar < 2){
		printf("if Npar = 1 use optimizeB\n");
		fprintf(logFile,"if Npar = 1 use optimizeB\n");
		fclose(logFile);
		exit(0);
	}
	if(Npar > 4){
		printf("Npar must be <= 4\n");
		fprintf(logFile,"Npar must be <= 4\n");
		fclose(logFile);
		exit(0);
	}
	double minVector[Ngrad];
	double tempVector[Ngrad];
	double prevVector[Ngrad];
	double currentVector[Ngrad];
	minVector[0] = minB/maxB;
	minVector[1] = minR/maxR;
	minVector[2] = minFraction;
	minVector[3] = minSigma;
	for(int i = 0; i < Ngrad; ++i){
		prevVector[i] = vector[i];
		currentVector[i] = vector[i];
	}
	double currentF = evaluateOptimizationFunctionGeneral(vector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g r = %g fraction = %10.7g epsilonB = %10.7g\n", vector[0]*maxB, vector[1]*maxR, vector[2], vector[3]);
	fprintf(logFile, "Bfactor = %g r = %g  fraction = %10.7g epsilonB = %10.7g\n", vector[0]*maxB, vector[1]*maxR, vector[2], vector[3]);
	for(int k = 0; k < Niterations; ++k) {
		///randomization;
		/*for(int j = 0; j < 5; ++j){
			for(int i = 0; i < Ngrad; ++i) {
				if(optPar[i]){
					//tempVector[i] = minVector[i] + (1.0 - minVector[i])*uniformDistribution();
					tempVector[i] = tempVector[i] + 0.2*minVector[i]*(0,5 - uniformDistribution());
					if(tempVector[i] > 1.0) {
						tempVector[i] = 1.0;
					}
					if(tempVector[i] < minVector[i]) {
						tempVector[i] = minVector[i];
					}
				}
			}
			
			double tempF = evaluateOptimizationFunction5(tempVector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);
			if(tempF < currentF){
				currentF = tempF;
				for(int i = 0; i < Ngrad; ++i) {
					vector[i] = tempVector[i];
				}
				printf("random search\n");
				fprintf(logFile, "random search\n");
			}
		}*/
		double prevF = currentF;
		//
		//valley first step
		double valley1[Ngrad];
		for(int i = 0; i < Ngrad; ++i) {
			valley1[i] = vector[i];
		}

		double grad[Ngrad];

		printf("optimization k = %d\n",k);
		fprintf(logFile, "optimization k = %d\n",k);
		fflush(logFile);

		for(int i = 0; i < Ngrad; ++i) {
			if (optPar[i]) {
				for (int j = 0; j < Ngrad; ++j) {
					tempVector[j] = vector[j];
				}
				double dx = fabs(vector[i]) / 100;
				if (k > 0) {
					dx = min(max(0.000001, fabs(currentVector[i] - prevVector[i]) / 10), fabs(currentVector[i] - minVector[i]) / 2);
				}
				tempVector[i] = vector[i] + dx;
				//currentF = evaluateOptimizationFunction5(vector[0]*maxB, vector[1]*maxN, vector[2]*maxFraction, vector[3]*maxR, vector[4]*maxV, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
				double f = evaluateOptimizationFunctionGeneral(tempVector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
				tempVector[i] = vector[i] - dx;
				double f1 = evaluateOptimizationFunctionGeneral(tempVector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);

				grad[i] = (f - f1) / (2 * dx);
				if (grad[i] != grad[i]) {
					printf("grad[i] = NaN\n");
					exit(0);
				}
				if (grad[i] > 0) {
					if (vector[i] - minVector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
				if (grad[i] < 0) {
					if (1.0 - vector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
			}
		}

		double gradNorm = 0;
		for(int i = 0; i < Ngrad; ++i){
			if (optPar[i]) {
				gradNorm = gradNorm + grad[i] * grad[i];
			}
		}
		gradNorm = sqrt(gradNorm);
		if(gradNorm <= 0){
			printf("gradNorm <= 0\n");
			continue;
		}
		if(gradNorm != gradNorm){
			printf("gradNorm = NaN\n");
			exit(0);
		}
		for(int i = 0; i < Ngrad; ++i) {
			if (optPar[i]) {
				grad[i] = grad[i] / gradNorm;
			}
		}


		findMinParametersGeneral(vector, optPar, grad, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, dopplerBeta);

		currentF = evaluateOptimizationFunctionGeneral(vector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
		//valley second step

		for(int i = 0; i < Ngrad; ++i) {
			if (optPar[i]) {
				for (int j = 0; j < Ngrad; ++j) {
					tempVector[j] = vector[j];
				}
				double dx = fabs(vector[i]) / 100;
				if (k > 0) {
					dx = min(max(0.000001, fabs(currentVector[i] - prevVector[i]) / 10), fabs(currentVector[i] - minVector[i]) / 2);
				}
				tempVector[i] = vector[i] + dx;
				double f = evaluateOptimizationFunctionGeneral(tempVector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
				tempVector[i] = vector[i] - dx;
				double f1 = evaluateOptimizationFunctionGeneral(tempVector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);

				grad[i] = (f - f1) / (2 * dx);
				if (grad[i] != grad[i]) {
					printf("grad[i] = NaN\n");
					exit(0);
				}
				if (grad[i] > 0) {
					if (vector[i] - minVector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
				if (grad[i] < 0) {
					if (1.0 - vector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
			}
		}

 		gradNorm = 0;
		for(int i = 0; i < Ngrad; ++i){
			if (optPar[i]) {
				gradNorm = gradNorm + grad[i] * grad[i];
			}
		}
		gradNorm = sqrt(gradNorm);
		if(gradNorm <= 0){
			printf("gradNorm <= 0\n");
			continue;
		}
		if(gradNorm != gradNorm){
			printf("gradNorm = NaN\n");
			exit(0);
		}
		for(int i = 0; i < Ngrad; ++i) {
			if (optPar[i]) {
				grad[i] = grad[i] / gradNorm;
			}
		}

		findMinParametersGeneral(vector, optPar, grad, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, dopplerBeta);

		currentF = evaluateOptimizationFunctionGeneral(vector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);

		double valley2[Ngrad];
		for(int i = 0; i < Ngrad; ++i) {
			valley2[i] = vector[i];
		}

		//// valley third step

		for(int i = 0; i < Ngrad; ++i) {
			if (optPar[i]) {
				grad[i] = valley1[i] - valley2[i];
			}
		}

		
		gradNorm = 0;
		for(int i = 0; i < Ngrad; ++i){
			if (optPar[i]) {
				gradNorm = gradNorm + grad[i] * grad[i];
			}
		}
		gradNorm = sqrt(gradNorm);
		if(gradNorm <= 0){
			printf("gradNorm <= 0\n");
		}
		if(gradNorm != gradNorm){
			printf("gradNorm = NaN\n");
			exit(0);
		}
		
		if(gradNorm > 0){
			for(int i = 0; i < Ngrad; ++i) {
				if (optPar[i]) {
					grad[i] = grad[i] / gradNorm;
				}
			}

			findMinParametersGeneral(vector, optPar, grad, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, dopplerBeta);

		}
		currentF = evaluateOptimizationFunctionGeneral(vector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
		
		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/
		for(int i = 0; i < Ngrad; ++i){
			prevVector[i] = currentVector[i];
			currentVector[i] = vector[i];
		}
		printf("optimization function = %g\n", currentF);
		fprintf(logFile, "optimization function = %g\n", currentF);
		printf("Bfactor = %g r = %g fraction = %10.7g epsilonB = %10.7g\n", vector[0]*maxB, vector[1]*maxR, vector[2], vector[3]);
		fprintf(logFile, "Bfactor = %g r = %g  fraction = %10.7g epsilonB = %10.7g\n", vector[0]*maxB, vector[1]*maxR, vector[2], vector[3]);
	}
	printf("finish optimization\n");
	fprintf(logFile, "finish optimization\n");
}

void optimizeParametersGeneral(double* vector, bool* optPar,  double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double dopplerBeta, FILE* logFile){
	double* observedError = new double[Nnu];
	for(int i = 0; i < Nnu; ++i){
		observedError[i] = 1.0;
	}

	optimizeParametersGeneral(vector, optPar, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta, logFile);

	delete[] observedError;
}


void optimizeParameterB(double& B, double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, const double& R, const double& fraction, const double& epsilonB, double dopplerBeta, FILE* logFile){
	double* observedError = new double[Nnu];
	for(int i = 0; i < Nnu; ++i){
		observedError[i] = 1.0;
	}
	
	double error = evaluateOptimizationFunction(nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, B, R, fraction, epsilonB, dopplerBeta);
	printf("optimization function = %g\n", error);
	fprintf(logFile, "optimization function = %g\n", error);
	const double Nb = 100;
	double factor = pow(maxB/minB, 1.0/(Nb - 1));
	double curMinB = minB;
	double curB = minB;
	double curMin = evaluateOptimizationFunction(nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, curMinB, R, fraction, epsilonB, dopplerBeta);
	for(int i = 1; i < Nb; ++i){
		curB = curB*factor;
		double curErr = evaluateOptimizationFunction(nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, curB, R, fraction, epsilonB, dopplerBeta);
		if(curErr < curMin){
			curMin = curErr;
			curMinB = curB;
		}
	}
	B = curMinB;
	double leftB = B/factor;
	double rightB = B*factor;
	double B1 = leftB + (rightB - leftB)/3.0;
	double B2 = leftB + 2*(rightB - leftB)/3.0;

	for(int i = 0; i < Niterations; ++i){
		double error1 = evaluateOptimizationFunction(nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, B1, R, fraction, epsilonB, dopplerBeta);
		double error2 = evaluateOptimizationFunction(nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, B2, R, fraction, epsilonB, dopplerBeta);
		if(error1 <= error2){
			rightB = B2;
		} else {
			leftB = B1;
		}
		B1 = leftB + (rightB - leftB)/3.0;
	    B2 = leftB + 2*(rightB - leftB)/3.0;
	}
	B = (B1 + B2)/2;
	error = evaluateOptimizationFunction(nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, B, R, fraction, epsilonB, dopplerBeta);
	printf("optimization function = %g\n", error);
	fprintf(logFile, "optimization function = %g\n", error);
	printf("Bfactor = %g\n", B);
	fprintf(logFile, "Bfactor = %g\n", B);
	delete[] observedError;
}

void stochasticGradientOptimization(double* vector, bool* optPar, double* nu, double* observedInu, double* observedError, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double dopplerBeta, FILE* logFile) {
	const int Ngrad = 4;
	int Npar = 0;
	for (int i = 0; i < Ngrad; ++i) {
		if (optPar[i]) {
			Npar++;
		}
	}
	if (Npar < 2) {
		printf("if Npar = 1 use optimizeB\n");
		fprintf(logFile, "if Npar = 1 use optimizeB\n");
		fclose(logFile);
		exit(0);
	}
	if (Npar > 4) {
		printf("Npar must be <= 4\n");
		fprintf(logFile, "Npar must be <= 4\n");
		fclose(logFile);
		exit(0);
	}

	double minVector[Ngrad];
	double tempVector[Ngrad];
	double prevVector[Ngrad];
	double currentVector[Ngrad];
	minVector[0] = minB / maxB;
	minVector[1] = minR / maxR;
	minVector[2] = minFraction;
	minVector[3] = minSigma;
	for (int i = 0; i < Ngrad; ++i) {
		prevVector[i] = vector[i];
		currentVector[i] = vector[i];
	}
	double currentF = evaluateOptimizationFunctionGeneral(vector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g r = %g fraction = %10.7g epsilonB = %10.7g\n", vector[0] * maxB, vector[1] * maxR, vector[2], vector[3]);
	fprintf(logFile, "Bfactor = %g r = %g  fraction = %10.7g epsilonB = %10.7g\n", vector[0] * maxB, vector[1] * maxR, vector[2], vector[3]);
	for (int k = 0; k < Niterations * Npar; ++k) {
		printf("optimization k = %d\n", k);
		fprintf(logFile, "optimization k = %d\n", k);
		fflush(logFile);
		int i = rand() % Ngrad;
		if (optPar[i]) {
			for (int j = 0; j < Ngrad; ++j) {
				tempVector[j] = vector[j];
			}
			
			double dx = vector[i] * 1E-6;
			tempVector[i] -= dx;
			double f1 = evaluateOptimizationFunctionGeneral(tempVector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
			if (f1 > currentF) {
				findMinAlongAxis(vector, i, vector[i], 1.0, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, dopplerBeta);
			}
			else {
				findMinAlongAxis(vector, i, minVector[i], vector[i], nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, dopplerBeta);
			}
			currentF = evaluateOptimizationFunctionGeneral(vector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);

			printf("optimization function = %g\n", currentF);
			fprintf(logFile, "optimization function = %g\n", currentF);
			printf("Bfactor = %g r = %g fraction = %10.7g epsilonB = %10.7g\n", vector[0] * maxB, vector[1] * maxR, vector[2], vector[3]);
			fprintf(logFile, "Bfactor = %g r = %g  fraction = %10.7g epsilonB = %10.7g\n", vector[0] * maxB, vector[1] * maxR, vector[2], vector[3]);
		}
	}
}

void stochasticGradientOptimization(double* vector, bool* optPar, double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double dopplerBeta, FILE* logFile) {
	double* observedError = new double[Nnu];
	for (int i = 0; i < Nnu; ++i) {
		observedError[i] = 1.0;
	}

	stochasticGradientOptimization(vector, optPar, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta, logFile);

	delete[] observedError;
}

void findMinAlongAxis(double* vector, int axis, double left, double right, double* nu, double* observedInu, double* observedError, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double& currentF, double dopplerBeta) {
	const int Ngrad = 4;
	double tempVector[4][Ngrad];
	for (int i = 0; i < Ngrad; ++i) {
		for (int j = 0; j < 4; ++j) {
			tempVector[j][i] = vector[i];
		}
	}
	tempVector[0][axis] = left;
	tempVector[1][axis] = left + (right - left) / 3;
	tempVector[2][axis] = left + 2 * (right - left) / 3;
	tempVector[3][axis] = right;

	double f = evaluateOptimizationFunctionGeneral(vector, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);

    //double f0 = evaluateOptimizationFunctionGeneral(tempVector[0], nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
	double f1 = evaluateOptimizationFunctionGeneral(tempVector[1], nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
	double f2 = evaluateOptimizationFunctionGeneral(tempVector[2], nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
	//double f3 = evaluateOptimizationFunctionGeneral(tempVector[3], nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);

	if (right - left < 1E-6) {
		if (f1 < f) {
			vector[axis] = tempVector[1][axis];
		}
		else if (f2 < f) {
			vector[axis] = tempVector[2][axis];
		} else {
			printf("x0 is already mininmum along axis %d\n", axis);
		}
		return;
	}

	if (f1 > f2) {
		findMinAlongAxis(vector, axis, tempVector[1][axis], right, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, dopplerBeta);
	} else {
		findMinAlongAxis(vector, axis, left, tempVector[2][axis], nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, dopplerBeta);
	}
}