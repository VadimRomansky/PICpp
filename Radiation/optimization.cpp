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

double evaluateOptimizationFunction5(double* vector, double* time, double** nu, double** observedInu, double** observedError, double* Ee, double**** dFe, int Np, int* Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, double*** psi , int*** thetaIndex, double*** concentrations, double***** nudoppler, double***** Inu, double***** Anu) {
	// v[0] = B, v[1] - N, v[2] - f, v[3] - v v[4] - r0 v[5] B/r^(a-1) v[6] N/r^(b-1) v[7] f/r^(f-1)
	//evaluateVolumeAndLength(area, length, rmax, fractionSize);
	double dopplerBeta =  0.75*vector[3]*maxV/speed_of_light;
	evaluateNuDoppler(nudoppler, Nmonth, Nnu, nu, dopplerBeta);
	double err = 0;
	for(int i = 0; i < Nmonth; ++i){
		double* totalInu = new double[Nnu[i]];
		double r = vector[4]*maxR0 + vector[3]*maxV*times[i];
		//important 0 or Nmonth-1!
		double rfactor = r/(vector[4]*maxR0 + vector[3]*maxV*times[0]);
		updateConcentartions(concentrations, vector[3] * maxV / speed_of_light);
		evaluateAllEmissivityAndAbsorption(nudoppler[i], Inu[i], Anu[i], Nnu[i], Ee, dFe, Np, Nd, Bn, sintheta, psi, thetaIndex, concentrations, vector[1]*maxN, vector[0]*maxB, rfactor, vector[5]*maxBpower, vector[6]*maxNpower, dopplerBeta);
		//
		//evaluateSpectrum(nu[i], totalInu, Inu[i], Anu[i], area, length, Nnu, rfactor);
		//
		double tempFraction = vector[2]*maxFraction/pow(rfactor, vector[7]*maxFpower - 1.0);
		if(geometry == SPHERICAL){
			evaluateSpectrumSpherical(nu[i], nudoppler[i], totalInu, Inu[i], Anu[i], r, Nnu[i], tempFraction, dopplerBeta);
		} else {
			evaluateSpectrumFlat(nu[i], nudoppler[i], totalInu, Inu[i], Anu[i], r, Nnu[i], tempFraction, dopplerBeta);
		}
		for(int j = 0; j < Nnu[i]; ++j){
			double err1 = 0;
			if(scale == LINEAR){
				err1 = sqr(totalInu[j] - observedInu[i][j])/sqr(observedError[i][j]);
			} else {
				err1 = sqr(log(totalInu[j]) - log(observedInu[i][j]));
			}
			//printf("%lf\n",err1);
			err = err + err1;
		}
		//printf("\n");
		delete[] totalInu;
	}
	
	return err;
}


void findMinParametersGeneral(double* vector, bool* optPar, const double* grad, double* time, double** nu, double** observedInu, double** observedError, double* Ee, double**** dFe, int Np, int* Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, double*** psi, int*** thetaIndex, double*** concentrations, double***** nudoppler, double***** Inu, double***** Anu, double& currentF){
	//v[0] = B, v[1] - N, v[2] - f, v[3] - v, v[4] - r0


	const int Ngrad = 8;
	int Npar = 0;
	for(int i = 0; i < Ngrad; ++i){
		if(optPar[i]){
			Npar++;
		}
	}
	if(Npar == 0){
		printf("Npar must be > 0\n");
		//fprintf(logFile, "Npar must be > 0\n");
		//fclose(logFile);
		exit(0);
	}

	double minVector[Ngrad];
	double tempVector1[Ngrad];
	double tempVector2[Ngrad];
	minVector[0] = minB/maxB;
	minVector[1] = minN/maxN;
	minVector[2] = minFraction/maxFraction;
	minVector[3] = minV/maxV;
	minVector[4] = minR0/maxR0;
	minVector[5] = minBpower*1.0/maxBpower;
	minVector[6] = minNpower*1.0/maxNpower;
	minVector[7] = minFpower*1.0/maxFpower;

	double maxLambda = sqrt(1.0*Npar);
	for(int i = 0; i < Ngrad; ++i) {
		if(optPar[i]){
			if(grad[i] > 0) {
				maxLambda = min(maxLambda, fabs((vector[i] - minVector[i])/grad[i]));
			}
			if(grad[i] < 0) {
				maxLambda = min(maxLambda, fabs((1.0 - vector[i])/grad[i]));
			}
		}
	}

	double step = maxLambda/2;

	for(int i = 0; i < Ngrad; ++i){
		tempVector1[i] = vector[i];
	}
	for(int i = 0; i < Ngrad; ++i) {
		if(optPar[i]){
			tempVector1[i] = vector[i] - grad[i]*step;
			if(tempVector1[i] != tempVector1[i]) {
				printf("tempvector[i] != tempVector[i]\n");
				exit(0);
			}
		}
	}

	double f = evaluateOptimizationFunction5(tempVector1, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
	if(f > currentF){
		int count = 0;
		while (f > currentF && count < 20){
			count++;
			step = step/2;
			
			for(int i = 0; i < Ngrad; ++i){
				tempVector1[i] = vector[i];
			}
			for(int i = 0; i < Ngrad; ++i) {
				if(optPar[i]){
					tempVector1[i] = vector[i] - grad[i]*step;
					if(tempVector1[i] != tempVector1[i]) {
						printf("tempvector[i] != tempVector[i]\n");
						exit(0);
					}
				}
			}

			f = evaluateOptimizationFunction5(tempVector1, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
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

	f = evaluateOptimizationFunction5(vector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);

	maxLambda = sqrt(1.0*Npar);
	for(int i = 0; i < Ngrad; ++i) {
		if(optPar[i]){
			if(grad[i] > 0) {
				maxLambda = min(maxLambda, fabs((vector[i] - minVector[i])/grad[i]));
			}
			if(grad[i] < 0) {
				maxLambda = min(maxLambda, fabs((1.0 - vector[i])/grad[i]));
			}
		}
	}

	for(int i = 0; i < Ngrad; ++i){
		tempVector1[i] = vector[i];
	}
	for(int i = 0; i < Ngrad; ++i) {
		if(optPar[i]){
			tempVector1[i] = vector[i] - maxLambda*grad[i];
			if(tempVector1[i] > 1.0) {
				tempVector1[i] = 1.0;
			}
			if(tempVector1[i] < minVector[i]) {
				tempVector1[i] = minVector[i];
			}
		}
	}
	double f3 = evaluateOptimizationFunction5(tempVector1, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);

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
			if(optPar[i]){
				tempVector1[i] = vector[i] - lambda1*grad[i];
				tempVector2[i] = vector[i] - lambda2*grad[i];
			}
		}

		double f1 = evaluateOptimizationFunction5(tempVector1, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
		double f2 = evaluateOptimizationFunction5(tempVector2, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);

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
		if(optPar[i]){
			tempVector1[i] = vector[i] - lambda*grad[i];
			if(tempVector1[i] > 1.0) {
				tempVector1[i] = 1.0;
			}
			if(tempVector1[i] < minVector[i]) {
				tempVector1[i] = minVector[i];
			}
		}
	}

	double tempF = evaluateOptimizationFunction5(tempVector1, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);

	if(f < tempF) {
		if(f < f3) {
			return;
		} else {
			for(int i = 0; i < Ngrad; ++i) {
				if(optPar[i]){
					vector[i] = vector[i] - maxLambda*grad[i];
				}
			}
		}
	} else {
		if(tempF < f3) {
			for(int i = 0; i < Ngrad; ++i) {
				if(optPar[i]){
					vector[i] = vector[i] - lambda*grad[i];
				}
			}
		} else {
			for(int i = 0; i < Ngrad; ++i) {
				if(optPar[i]){
					vector[i] = vector[i] - maxLambda*grad[i];
				}
			}
		}
	}
	for(int i  = 0; i < Ngrad; ++i){
		if(optPar[i]){
			if(vector[i] > 1.0) {
				vector[i] = 1.0;
			}
			if(vector[i] < minVector[i]) {
				vector[i] = minVector[i];
			}
		}
	}
}

void optimizeParametersGeneral(double* vector, bool* optPar, double* time,  double** nu, double** observedInu, double** observedError, double* Ee, double**** dFe, int Np, int* Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, double*** psi, int*** thetaIndex, double*** concentrations, double***** nudoppler, double***** Inu, double***** Anu, FILE* logFile){
	const int Ngrad = 8;
	int Npar = 0;
	for(int i = 0; i < Ngrad; ++i){
		if(optPar[i]){
			Npar++;
		}
	}
	if(Npar < 2){
		printf("if Npar = 1 use optimizeB\n");
		fprintf(logFile,"if Npar = 1 use optimizeB\n");
		fclose(logFile);
		exit(0);
	}
	if(Npar > 8){
		printf("Npar must be <= 8\n");
		fprintf(logFile,"Npar must be <= 8\n");
		fclose(logFile);
		exit(0);
	}
	double minVector[Ngrad];
	double minBoundVector[Ngrad];
	double tempVector[Ngrad];
	double prevVector[Ngrad];
	double currentVector[Ngrad];
	minVector[0] = minB/maxB;
	minVector[1] = minN/maxN;
	minVector[2] = minFraction/maxFraction;
	minVector[3] = minV/maxV;
	minVector[4] = minR0/maxR0;
	minVector[5] = minBpower*1.0/maxBpower;
	minVector[6] = minNpower*1.0/maxNpower;
	minVector[7] = minFpower*1.0/maxFpower;
	minBoundVector[0] = minBbound/maxB;
	minBoundVector[1] = minNbound/maxN;
	minBoundVector[2] = minFbound/maxFraction;
	minBoundVector[3] = minVbound/maxV;
	minBoundVector[4] = minRbound/maxR0;
	minBoundVector[5] = minBpbound/maxBpower;
	minBoundVector[6] = minNpbound/maxNpower;
	minBoundVector[7] = minFpbound/maxFpower;
	for(int i = 0; i < Ngrad; ++i){
		prevVector[i] = vector[i];
		currentVector[i] = vector[i];
	}
	double currentF = evaluateOptimizationFunction5(vector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g n = %g fraction = %10.7g v/c = %10.7g r0 = %10.7g a = %10.7g b = %10.7g fpower = %10.7g\n", vector[0]*maxB, vector[1]*maxN, vector[2]*maxFraction, vector[3]*maxV/speed_of_light, vector[4]*maxR0, vector[5]*maxBpower, vector[6]*maxNpower, vector[7]*maxFpower);
	fprintf(logFile, "Bfactor = %g n = %g  fraction = %10.7g v/c = %10.7g r0 = %10.7g a = %10.7g b = %10.7g fpower = %10.7g\n", vector[0]*maxB, vector[1]*maxN, vector[2]*maxFraction, vector[3]*maxV/speed_of_light, vector[4]*maxR0, vector[5]*maxBpower, vector[6]*maxNpower, vector[7]*maxFpower);
	for(int k = 0; k < Niterations; ++k) {
		///randomization;
		for(int j = 0; j < 5; ++j){
			for(int i = 0; i < Ngrad; ++i) {
				tempVector[i] = vector[i];
				if(optPar[i]){
					//tempVector[i] = minVector[i] + (1.0 - minVector[i])*uniformDistribution();
					tempVector[i] = tempVector[i] + 0.2*minVector[i]*(0.5 - uniformDistribution());
					if(tempVector[i] > 1.0) {
						tempVector[i] = 1.0;
					}
					if(tempVector[i] < minVector[i]) {
						tempVector[i] = minVector[i];
					}
				}
			}
			
			double tempF = evaluateOptimizationFunction5(tempVector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
			if(tempF < currentF){
				currentF = tempF;
				for(int i = 0; i < Ngrad; ++i) {
					vector[i] = tempVector[i];
				}
				printf("random search\n");
				fprintf(logFile, "random search\n");
			}
		}
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
			if(optPar[i]){
				for(int j = 0; j < Ngrad; ++j) {
					tempVector[j] = vector[j];
				}
				double dx = fabs(vector[i])/100;
				if(k > 0){
					dx = max(0.000001,fabs(currentVector[i] - prevVector[i])/10);
				}
				dx = min(0.5*fabs(vector[i] - minBoundVector[i]), dx);
				tempVector[i] = vector[i] + dx;
				//currentF = evaluateOptimizationFunction5(vector[0]*maxB, vector[1]*maxN, vector[2]*maxFraction, vector[3]*maxR, vector[4]*maxV, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
				double f = evaluateOptimizationFunction5(tempVector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
				tempVector[i] = vector[i] - dx;
				double f1 = evaluateOptimizationFunction5(tempVector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);

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
		}

		double gradNorm = 0;
		for(int i = 0; i < Ngrad; ++i){
			if(optPar[i]){
				gradNorm = gradNorm + grad[i]*grad[i];
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
			if(optPar[i]){
				grad[i] = grad[i]/gradNorm;
			}
		}


		findMinParametersGeneral(vector, optPar, grad, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu, currentF);

		currentF = evaluateOptimizationFunction5(vector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
		//valley second step

		for(int i = 0; i < Ngrad; ++i) {
			if(optPar[i]){
				for(int j = 0; j < Ngrad; ++j) {
					tempVector[j] = vector[j];
				}
				double dx = fabs(vector[i])/100;
				if(k > 0){
					dx = max(0.000001,fabs(currentVector[i] - prevVector[i])/10);
				}
				dx = min(0.5*fabs(vector[i] - minBoundVector[i]), dx);
				tempVector[i] = vector[i] + dx;
				double f = evaluateOptimizationFunction5(tempVector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
				tempVector[i] = vector[i] - dx;
				double f1 = evaluateOptimizationFunction5(tempVector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);

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
		}

 		gradNorm = 0;
		for(int i = 0; i < Ngrad; ++i){
			if(optPar[i]){
				gradNorm = gradNorm + grad[i]*grad[i];
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
			if(optPar[i]){
				grad[i] = grad[i]/gradNorm;
			}
		}

		findMinParametersGeneral(vector, optPar, grad,time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu, currentF);

		currentF = evaluateOptimizationFunction5(vector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);

		double valley2[Ngrad];
		for(int i = 0; i < Ngrad; ++i) {
			valley2[i] = vector[i];
		}

		//// valley third step

		for(int i = 0; i < Ngrad; ++i) {
			if(optPar[i]){
				grad[i] = valley1[i] - valley2[i];
			}
		}

		
		gradNorm = 0;
		for(int i = 0; i < Ngrad; ++i){
			if(optPar[i]){
				gradNorm = gradNorm + grad[i]*grad[i];
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
				if(optPar[i]){
					grad[i] = grad[i]/gradNorm;
				}
			}

			findMinParametersGeneral(vector, optPar, grad, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu, currentF);

		}
		currentF = evaluateOptimizationFunction5(vector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
		
		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/
		for(int i = 0; i < Ngrad; ++i){
			prevVector[i] = currentVector[i];
			currentVector[i] = vector[i];
		}
		printf("optimization function = %g\n", currentF);
		fprintf(logFile, "optimization function = %g\n", currentF);
		printf("Bfactor = %g n = %g fraction = %10.7g v/c = %10.7g r0 = %10.7g a = %10.7g b = %10.7g fpower = %10.7g\n", vector[0]*maxB, vector[1]*maxN, vector[2]*maxFraction, vector[3]*maxV/speed_of_light, vector[4]*maxR0, vector[5]*maxBpower, vector[6]*maxNpower, vector[7]*maxFpower);
		fprintf(logFile, "Bfactor = %g n = %g  fraction = %10.7g v/c = %10.7g r0 = %10.7g a = %10.7g b = %10.7g fpower = %10.7g\n", vector[0]*maxB, vector[1]*maxN, vector[2]*maxFraction, vector[3]*maxV/speed_of_light, vector[4]*maxR0, vector[5]*maxBpower, vector[6]*maxNpower, vector[7]*maxFpower);
	}
	printf("finish optimization\n");
	fprintf(logFile, "finish optimization\n");
}

void optimizeParametersGeneral(double* vector, bool* optPar, double* time, double** nu, double** observedInu, double* Ee, double**** dFe, int Np, int* Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, double*** psi, int*** thetaIndex, double*** concentrations, double***** nudoppler, double***** Inu, double***** Anu, FILE* logFile){
	double** observedError = new double*[Nmonth];
	for(int i = 0; i < Nmonth; ++i){
		observedError[i] = new double[Nnu[i]];
	}

	optimizeParametersGeneral(vector, optPar, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu, logFile);

	for(int i = 0; i < Nmonth; ++i){
		delete[] observedError[i];
	}
	delete[] observedError;
}

void stochasticGradientOptimization(double* vector, bool* optPar, double* time, double** nu, double** observedInu, double** observedError, double* Ee, double**** dFe, int Np, int* Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, double*** psi, int*** thetaIndex, double*** concentrations, double***** nudoppler, double***** Inu, double***** Anu, FILE* logFile) {
	const int Ngrad = 8;
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
	if (Npar > 8) {
		printf("Npar must be <= 8\n");
		fprintf(logFile, "Npar must be <= 8\n");
		fclose(logFile);
		exit(0);
	}

	double minVector[Ngrad];
	double tempVector[Ngrad];
	double prevVector[Ngrad];
	double currentVector[Ngrad];
	minVector[0] = minB / maxB;
	minVector[1] = minN / maxN;
	minVector[2] = minFraction / maxFraction;
	minVector[3] = minV / maxV;
	minVector[4] = minR0 / maxR0;
	minVector[5] = minBpower * 1.0 / maxBpower;
	minVector[6] = minNpower * 1.0 / maxNpower;
	minVector[7] = minFpower * 1.0 / maxFpower;
	for (int i = 0; i < Ngrad; ++i) {
		prevVector[i] = vector[i];
		currentVector[i] = vector[i];
	}
	double currentF = evaluateOptimizationFunction5(vector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g n = %g fraction = %10.7g v/c = %10.7g r0 = %10.7g a = %10.7g b = %10.7g fpower = %10.7g\n", vector[0] * maxB, vector[1] * maxN, vector[2] * maxFraction, vector[3] * maxV / speed_of_light, vector[4] * maxR0, vector[5] * maxBpower, vector[6] * maxNpower, vector[7] * maxFpower);
	fprintf(logFile, "Bfactor = %g n = %g  fraction = %10.7g v/c = %10.7g r0 = %10.7g a = %10.7g b = %10.7g fpower = %10.7g\n", vector[0] * maxB, vector[1] * maxN, vector[2] * maxFraction, vector[3] * maxV / speed_of_light, vector[4] * maxR0, vector[5] * maxBpower, vector[6] * maxNpower, vector[7] * maxFpower);
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
			double f1 = evaluateOptimizationFunction5(tempVector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
			if (f1 > currentF) {
				findMinAlongAxis(vector, i, vector[i], 1.0, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu, logFile);
			}
			else {
				findMinAlongAxis(vector, i, minVector[i], vector[i], time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu, logFile);
			}
			currentF = evaluateOptimizationFunction5(vector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);

			printf("optimization function = %g\n", currentF);
			fprintf(logFile, "optimization function = %g\n", currentF);
			printf("Bfactor = %g n = %g fraction = %10.7g v/c = %10.7g r0 = %10.7g a = %10.7g b = %10.7g fpower = %10.7g\n", vector[0] * maxB, vector[1] * maxN, vector[2] * maxFraction, vector[3] * maxV / speed_of_light, vector[4] * maxR0, vector[5] * maxBpower, vector[6] * maxNpower, vector[7] * maxFpower);
			fprintf(logFile, "Bfactor = %g n = %g  fraction = %10.7g v/c = %10.7g r0 = %10.7g a = %10.7g b = %10.7g fpower = %10.7g\n", vector[0] * maxB, vector[1] * maxN, vector[2] * maxFraction, vector[3] * maxV / speed_of_light, vector[4] * maxR0, vector[5] * maxBpower, vector[6] * maxNpower, vector[7] * maxFpower);
		}
	}
}
void findMinAlongAxis(double* vector, int axis, double left, double right, double* time, double** nu, double** observedInu, double** observedError, double* Ee, double**** dFe, int Np, int* Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, double*** psi, int*** thetaIndex, double*** concentrations, double***** nudoppler, double***** Inu, double***** Anu, FILE* logFile) {
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

	double f = evaluateOptimizationFunction5(vector, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);

	//double f0 = evaluateOptimizationFunctionGeneral(tempVector[0], nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);
	double f1 = evaluateOptimizationFunction5(tempVector[1], time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
	double f2 = evaluateOptimizationFunction5(tempVector[2], time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu);
	//double f3 = evaluateOptimizationFunctionGeneral(tempVector[3], nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, dopplerBeta);

	if (right - left < 1E-6) {
		if (f1 < f) {
			vector[axis] = tempVector[1][axis];
		}
		else if (f2 < f) {
			vector[axis] = tempVector[2][axis];
		}
		else {
			printf("x0 is already mininmum along axis %d\n", axis);
		}
		return;
	}

	if (f1 > f2) {
		findMinAlongAxis(vector, axis, tempVector[1][axis], right, time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu, logFile);
	}
	else {
		findMinAlongAxis(vector, axis, left, tempVector[2][axis], time, nu, observedInu, observedError, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, psi, thetaIndex, concentrations, nudoppler, Inu, Anu, logFile);
	}
}