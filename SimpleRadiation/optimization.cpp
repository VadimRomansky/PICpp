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

double evaluateOptimizationFunction5(const double* vector, double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double fraction, double epsilonB){
	// v[0] = B, v[1] = R
	double* totalInu = new double[Nnu];
	//evaluateVolumeAndLength(area, length, rmax, fractionSize);
	double err = 0;
	double r = vector[1]*maxR;
	double B = vector[0]*maxB;
	double n = B*B/(4*pi*massProtonReal*speed_of_light2*epsilonB);
	evaluateAllEmissivityAndAbsorption(nu, Inu, Anu, Nnu, Ee, dFe, Np, Nd, Bn, sintheta, thetaIndex, concentrations, n, B, 1.0);
		//
		//evaluateSpectrum(nu[i], totalInu, Inu[i], Anu[i], area, length, Nnu, rfactor);
		//
		if(geometry == SPHERICAL){
			evaluateSpectrumSpherical(nu, totalInu, Inu, Anu, r, Nnu, 1.0, fraction);
		} else {
			evaluateSpectrumFlat(nu, totalInu, Inu, Anu, r, Nnu, 1.0, fraction);
		}
		for(int j = 0; j < Nnu; ++j){
			double err1 = 0;
			if(scale == LINEAR){
				err1 = sqr(totalInu[j] - observedInu[j]);
			} else {
				err1 = sqr(log(totalInu[j]) - log(observedInu[j]));
			}
			err = err + err1;
		}
	
	delete[] totalInu;
	
	return err;
}

void findMinParameters5(double* vector, const double* grad, double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double& currentF, double fraction, double epsilonB){
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

	double f = evaluateOptimizationFunction5(tempVector1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);
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

			f = evaluateOptimizationFunction5(tempVector1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);
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

	f = evaluateOptimizationFunction5(vector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);

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
	double f3 = evaluateOptimizationFunction5(tempVector1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);

	double leftLambda = 0;
	double rightLambda = maxLambda;


	for(int j = 0; j < 10; ++j) {
		double lambda1 = leftLambda + (rightLambda - leftLambda)/3.0;
		double lambda2 = leftLambda + (rightLambda - leftLambda)*2.0/3.0;

		for(int i = 0; i < Ngrad; ++i) {
			tempVector1[i] = vector[i] - lambda1*grad[i];
			tempVector2[i] = vector[i] - lambda2*grad[i];
		}

		double f1 = evaluateOptimizationFunction5(tempVector1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);
		double f2 = evaluateOptimizationFunction5(tempVector2, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);

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

	double tempF = evaluateOptimizationFunction5(tempVector1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);

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

void optimizeParameters5(double* vector,  double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double fraction, double epsilonB, FILE* logFile){
	const int Ngrad = 2;
	double minVector[Ngrad];
	double tempVector[Ngrad];
	minVector[0] = minB/maxB;
	minVector[1] = minR/maxR;
	double currentF = evaluateOptimizationFunction5(vector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g rmax = %g\n", vector[0]*maxB, vector[1]*maxR);
	fprintf(logFile, "Bfactor = %grmax = %g\n", vector[0]*maxB, vector[1]*maxR);
	for(int k = 0; k < Niterations; ++k) {
		///randomization;
		/*for(int j = 0; j < 5; ++j){
			for(int i = 0; i < Ngrad; ++i) {
				tempVector[i] = minVector[i] + (1.0 - minVector[i])*uniformDistribution();
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
			double f = evaluateOptimizationFunction5(tempVector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);
			tempVector[i] = vector[i] - dx;
			double f1 = evaluateOptimizationFunction5(tempVector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);

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


		findMinParameters5(vector, grad, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, fraction, epsilonB);

		currentF = evaluateOptimizationFunction5(vector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);
		//valley second step

		for(int i = 0; i < Ngrad; ++i) {
			for(int j = 0; j < Ngrad; ++j) {
				tempVector[j] = vector[j];
			}
			double dx = fabs(vector[i])/1000;
			tempVector[i] = vector[i] + dx;
			double f = evaluateOptimizationFunction5(tempVector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);
			tempVector[i] = vector[i] - dx;
			double f1 = evaluateOptimizationFunction5(tempVector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);

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

		findMinParameters5(vector, grad, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, fraction, epsilonB);

		currentF = evaluateOptimizationFunction5(vector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);

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

			findMinParameters5(vector, grad, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF, fraction, epsilonB);

		}
		currentF = evaluateOptimizationFunction5(vector, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, fraction, epsilonB);
		
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