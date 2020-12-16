#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include <time.h>
#include <cmath>
#include <string>
#include <omp.h>

#include "constants.h"
#include "spectrum.h"
#include "util.h"

#include "optimization.h"

double evaluateOptimizationFunction(double Bfactor, double n, double fractionSize, double rmax, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double** Bn, double** sintheta, int** thetaIndex, double** concentrations, double*** Inu, double*** Anu, double** area, double** length, double* Rho, double* Phi, double* totalInu) {
	if(fractionSize > 1.0){
		fractionSize = 1.0;
	}
	evaluateAllEmissivityAndAbsorption1(nu, Inu, Anu, Nnu, Ee, Fe, Np, Nd, Bn, sintheta, thetaIndex, concentrations, n, Bfactor, fractionSize);
	evaluateVolumeAndLength(area, length, rmax, Rho, Phi, fractionSize);
	evaluateSpectrum(nu, totalInu, Inu, Anu, area, length, Nnu, rmax, Rho, Phi);

	double magnetization = sqr(Bfactor*Bn[Ntheta/2][0])/(4*pi*n*concentrations[Ntheta/2][0]*massProtonReal*speed_of_light2);
	double magnErr = 0;
	if(magnetization > maxSigma){
		magnErr = (magnetization)*1E6;
	}
	double consErr = 0;
	if(n*concentrations[Ntheta/2][0] > maxN){
		consErr = (n*concentrations[Ntheta/2][0])*1E3;
	}

	double I0 = totalInu[0] - apry[0];
	double I1 = totalInu[1] - apry[1];
	double I2 = totalInu[2] - apry[2];
	double I3 = totalInu[3] - apry[3];
	//double I4 = totalInu[4] - augy[4];
	//double I5 = totalInu[5] - augy[5];

	double fracErr = 0.001*sqr(1.0/fractionSize) + 0.001*sqr(1.0/(1.0 - fractionSize));

	//double I0 = log(totalInu[0]) - log(augy[0]);
	//double I1 = log(totalInu[1]) - log(augy[1]);
	//double I2 = log(totalInu[2]) - log(augy[2]);
	//double I3 = log(totalInu[3]) - log(augy[3]);
	//double I4 = log(totalInu[4]) - log(augy[4]);

	return I0*I0 + I1*I1 + I2*I2 + I3*I3;// + magnErr;
	//return fabs(I0) + fabs(I1) + fabs(I2) + fabs(I3) + fabs(I4);
	//return fabs(I1);
}

void findMinParameters(double& Bfactor, double& N, double& fractionSize, double rmax, double minLambda, double maxLambda, double gradB, double gradn, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double** Bn, double** sintheta, int** thetaIndex, double** concentrations, double*** Inu, double*** Anu, double** area, double** length, double* Rho, double* Phi, double* totalInu) {
	if(maxLambda - minLambda < 0.01*maxLambda) {
		N = N - maxLambda*gradn;
		Bfactor = Bfactor - (maxLambda + minLambda)*0.5*gradB;
		return;
	}
	double lambda1 = minLambda + (maxLambda - minLambda)/3.0;
	double lambda2 = minLambda + (maxLambda - minLambda)*2.0/3.0;

	double N1 = N - lambda1*gradn;
	double B1 = Bfactor - lambda1*gradB;

	double N2 = N - lambda2*gradn;
	double B2 = Bfactor - lambda2*gradB;

	//double concentration1 = evaluateConcentrationFromB(B*b1, gamma0, sigma);
	//double concentration2 = evaluateConcentrationFromB(B*b2, gamma0, sigma);

	//double f = evaluateOptimizationFunction(Bfactor, N, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, distance, Inu, Anu, X, Y, Z, rmax);
	double f1 = evaluateOptimizationFunction(B1, N1, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	double f2 = evaluateOptimizationFunction(B2, N2, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	if(f1 < f2) {
		findMinParameters(Bfactor, N, fractionSize, rmax, minLambda, lambda2, gradB, gradn, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	} else {
		findMinParameters(Bfactor, N, fractionSize, rmax, lambda1, maxLambda, gradB, gradn, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	}
}

void findMinParameters(double& Bfactor, double& N, double& fractionSize, double rmax, double gradB, double gradN, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double** Bn, double** sintheta, int** thetaIndex, double** concentrations, double*** Inu, double*** Anu, double** area, double** length, double* Rho, double* Phi, double& currentF, double* totalInu) {
	double minLambda = 0;
	double lambdaB = fabs(maxB/gradB);
	double lambdaN = fabs(maxN/gradN);
	double maxLambda = 1.0*min(lambdaN, lambdaB);
	if(N - maxLambda*gradN < 0 ){
		maxLambda = N/gradN;
	}
	if(Bfactor - maxLambda*gradB < 0 ){
		maxLambda = Bfactor/gradB;
	}

	double step = 0.4*min(fabs(Bfactor/gradB), fabs(N/gradN));
	double B1 = Bfactor - gradB*step;
	double N1 = N - gradN*step;
	if(B1 != B1){
		printf("B1 Nan\n");
		exit(0);
	}
	if(N1 != N1){
		printf("N1 Nan\n");
		exit(0);
	}
	double f1 = evaluateOptimizationFunction(B1, N1, fractionSize, rmax,nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	if(f1 > currentF){
		while (f1 > currentF){
			step = step/2;
			B1 = Bfactor - gradB*step;
			N1 = N - gradN*step;
			if(B1 != B1){
				printf("B1 Nan\n");
				exit(0);
			}
			if(N1 != N1){
				printf("N1 Nan\n");
				exit(0);
			}
			f1 = evaluateOptimizationFunction(B1, N1, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		Bfactor = B1;
		N = N1;
		return;
	}
	Bfactor = B1;
	N = N1;
	step = 0.4*min(fabs(Bfactor/gradB), fabs(N/gradN));
	double B2 = B1 - gradB*step;
	double N2 = N1 - gradN*step;
	if(B2 != B2){
		printf("B2 Nan\n");
		exit(0);
	}
	if(N2 != N2){
		printf("N2 Nan\n");
		exit(0);
	}
	double f2 = evaluateOptimizationFunction(B2, N2, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	int iterations = 0;
	while(f2 < f1){
		iterations++;
		Bfactor = B2;
		N = N2;
		double step = 0.4*min(fabs(Bfactor/gradB), fabs(N/gradN));
		B2 = Bfactor - gradB*step;
		N2 = N - gradN*step;
		if(B2 != B2){
			printf("B2 Nan\n");
			exit(0);
		}
		if(N2 != N2){
			printf("N2 Nan\n");
			exit(0);
		}
		f1 = f2;
		f2 = evaluateOptimizationFunction(B2, N2, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	}

	//findMinParameters(Bfactor, N, minLambda, maxLambda, gradB, gradN, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
}

void optimizeParameters(double& Bfactor, double& N, double& fractionSize, double rmax, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double** Bn, double** sintheta, int** thetaIndex, double** concentrations, double*** Inu, double*** Anu, double** area, double** length, double* Rho, double* Phi, FILE* logFile) {
	double* totalInu = new double[Nnu];
	double currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g n = %g\n", Bfactor, N);
	fprintf(logFile, "Bfactor = %g n = %g\n", Bfactor, N);
	for(int i = 0; i < Niterations; ++i) {
		///randomization;
		for(int j = 0; j < 5; ++j){
			double tempN = N + 0.5*N*(uniformDistribution()-0.5);
			double tempB = Bfactor + 0.5*Bfactor*(uniformDistribution() - 0.5);
			double tempF = evaluateOptimizationFunction(tempB, tempN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			if(tempF < currentF){
				currentF = tempF;
				Bfactor = tempB;
				N = tempN;
				printf("random search\n");
				fprintf(logFile, "random search\n");
			}
		}
		double prevF = currentF;
		//
		//valley first step
		double valleyB1 = Bfactor;
		double valleyN1 = N;

		double dxB = fabs(Bfactor)/20;
		double dxN = fabs(N)/20;
		printf("optimization i = %d\n",i);
		fprintf(logFile, "optimization i = %d\n",i);
		fflush(logFile);

		double gradB;
		double gradN;

		double Fb = evaluateOptimizationFunction(Bfactor + dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		double Fb1 = evaluateOptimizationFunction(Bfactor - dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fb > currentF && Fb1 > currentF){
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction(Bfactor + dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fb1 = evaluateOptimizationFunction(Bfactor - dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		double Fn = evaluateOptimizationFunction(Bfactor, N + dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		double Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fn > currentF && Fn1 > currentF){
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction(Bfactor, N + dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		double gradNorm = sqrt(gradB*gradB + gradN*gradN);
		gradB = gradB/gradNorm;
		gradN = gradN/gradNorm;

		findMinParameters(Bfactor, N, fractionSize, rmax, gradB, gradN, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, currentF, totalInu);

		currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		//valley second step

		dxB = fabs(Bfactor)/20;
		dxN = fabs(N)/20;

		Fb = evaluateOptimizationFunction(Bfactor + dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		Fb1 = evaluateOptimizationFunction(Bfactor - dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fb > currentF && Fb1 > currentF){
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction(Bfactor + dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fb1 = evaluateOptimizationFunction(Bfactor - dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		Fn = evaluateOptimizationFunction(Bfactor, N + dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fn > currentF && Fn1 > currentF){
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction(Bfactor, N + dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradN = (Fn - Fn1)/(2*dxN);
		gradNorm = sqrt(gradB*gradB + gradN*gradN);
		gradB = gradB/gradNorm;
		gradN = gradN/gradNorm;

		findMinParameters(Bfactor, N, fractionSize, rmax, gradB, gradN, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, currentF, totalInu);

		currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		double valleyB2 = Bfactor;
		double valleyN2 = N;

		//// valley third step

		
		gradB = (valleyB1 - valleyB2);
		gradN = (valleyN1 - valleyN2);
		gradNorm = sqrt(gradB*gradB + gradN*gradN);
		if(gradNorm > 0){
			gradB = gradB/gradNorm;
			gradN = gradN/gradNorm;

			double step = 0.1*min(Bfactor, N);
			double Fv = evaluateOptimizationFunction(Bfactor - gradB*step, N - gradN*step, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

			if(Fv < currentF){
				findMinParameters(Bfactor, N, fractionSize, rmax, gradB, gradN, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, currentF, totalInu);
			} else {
				findMinParameters(Bfactor, N, fractionSize, rmax, -gradB, -gradN, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, currentF, totalInu);
			}
			currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/

		printf("optimization function = %g\n", currentF);
		fprintf(logFile, "optimization function = %g\n", currentF);
		printf("Bfactor = %g n = %g\n", Bfactor, N);
		fprintf(logFile, "Bfactor = %g n = %g\n", Bfactor, N);
	}
	printf("finish optimization\n");
	fprintf(logFile, "finish optimization\n");
	delete[] totalInu;
}

void findMinParameters3(double& Bfactor, double& N, double& fractionSize, double rmax, double minLambda, double maxLambda, double gradB, double gradn, double gradS, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double** Bn, double** sintheta, int** thetaIndex, double** concentrations, double*** Inu, double*** Anu, double** area, double** length, double* Rho, double* Phi, double* totalInu) {
	if(maxLambda - minLambda < 0.01*maxLambda) {
		N = N - maxLambda*gradn;
		Bfactor = Bfactor - (maxLambda + minLambda)*0.5*gradB;
		return;
	}
	double lambda1 = minLambda + (maxLambda - minLambda)/3.0;
	double lambda2 = minLambda + (maxLambda - minLambda)*2.0/3.0;

	double N1 = N - lambda1*gradn;
	double B1 = Bfactor - lambda1*gradB;
	double S1 = fractionSize - lambda1*gradS;
	if(S1 > 1){
		S1 = 1.0;
		printf("S1 > 1\n");
	}
	if(S1 < 0){
		S1 = 1E-9;
		printf("S1 < 0\n");
	}

	double N2 = N - lambda2*gradn;
	double B2 = Bfactor - lambda2*gradB;
	double S2 = fractionSize - lambda2*gradS;
	if(S2 > 1.0){
		S2 = 1.0;
		printf("S2 > 1\n");
	}
	if(S2 < 0){
		S2 = 1E-9;
		printf("S2 < 0\n");
	}

	//double concentration1 = evaluateConcentrationFromB(B*b1, gamma0, sigma);
	//double concentration2 = evaluateConcentrationFromB(B*b2, gamma0, sigma);

	//double f = evaluateOptimizationFunction(Bfactor, N, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, distance, Inu, Anu, X, Y, Z, rmax);
	double f1 = evaluateOptimizationFunction(B1, N1, S1, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	double f2 = evaluateOptimizationFunction(B2, N2, S2, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	if(f1 < f2) {
		findMinParameters3(Bfactor, N, fractionSize, rmax, minLambda, lambda2, gradB, gradn, gradS, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	} else {
		findMinParameters3(Bfactor, N, fractionSize, rmax, lambda1, maxLambda, gradB, gradn, gradS, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	}
}

void findMinParameters3(double& Bfactor, double& N, double& fractionSize, double rmax, double gradB, double gradN, double gradS, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double** Bn, double** sintheta, int** thetaIndex, double** concentrations, double*** Inu, double*** Anu, double** area, double** length, double* Rho, double* Phi, double& currentF, double* totalInu) {

	double step = 0.4;
	if(fabs(gradB) > 0){
		step = 0.4*fabs(Bfactor/gradB);
	}
	if(fabs(gradN) > 0){
		step = min(step, 0.4*fabs(N/gradN));
	}
	if(gradS > 0){
		step = min(step, 0.4*fabs((maxFraction - fractionSize)/gradS));
	}
	if(gradS < 0){
		step = min(step, 0.4*fabs((fractionSize - minFraction)/gradS));
	}
	double B1 = Bfactor - gradB*step;
	double N1 = N - gradN*step;
	double S1 = fractionSize - gradS*step;
	if(B1 != B1){
		printf("1 Nan\n");
		exit(0);
	}
	if(S1 > maxFraction){
		S1 = maxFraction;
		//printf("S1 > 1\n");
	}
	if(S1 <= minFraction){
		S1 = minFraction;
		//printf("S1 < 0\n");
	}
	double f1 = evaluateOptimizationFunction(B1, N1, S1, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	if(f1 > currentF){
		while (f1 > currentF){
			step = step/2;
			B1 = Bfactor - gradB*step;
			N1 = N - gradN*step;
			S1 = fractionSize - gradS*step;
			if(B1 != B1){
				printf("B1 Nan\n");
				exit(0);
			}
			if(S1 > maxFraction){
				//printf("S1 > 1\n");
				S1 = maxFraction;
			}
			if(S1 <= minFraction){
				//printf("S1 < 0\n");
				S1 = minFraction;
			}
			f1 = evaluateOptimizationFunction(B1, N1, S1, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		Bfactor = B1;
		N = N1;
		fractionSize = S1;
		return;
	}
	Bfactor = B1;
	N = N1;
	fractionSize = S1;
	step = 0.4*min3(fabs(Bfactor/gradB), fabs(N/gradN), fabs(fractionSize/gradS));
	double B2 = B1 - gradB*step;
	double N2 = N1 - gradN*step;
	double S2 = S1 - gradS*step;
	if(B2 != B2){
		printf("B2 Nan\n");
		exit(0);
	}
	if(S2 > maxFraction){
		printf("S2 > 1\n");
		S2 = maxFraction;
	}
	if(S2 < minFraction){
		S2 = minFraction;
	}
	double f2 = evaluateOptimizationFunction(B2, N2, S2, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	int iterations = 0;
	while(f2 < f1){
		iterations++;
		Bfactor = B2;
		N = N2;
		fractionSize = S2;
		double step = 0.4*min4(fabs(Bfactor/gradB), fabs(N/gradN), fabs(fractionSize/gradS), fabs((1.0 - fractionSize)/gradS));
		B2 = Bfactor - gradB*step;
		N2 = N - gradN*step;
		S2 = fractionSize - gradS*step;
		if(B2 != B2){
			printf("B2 Nan\n");
			exit(0);
		}
		if(S2 > maxFraction){
			printf("S2 > 1\n");
			S2 = maxFraction;
			return;
		}
		if(S2 < minFraction){
			S2 = minFraction;
			//printf("S2 < 0\n");
			return;
		}
		f1 = f2;
		f2 = evaluateOptimizationFunction(B2, N2, S2, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	}

	//findMinParameters(Bfactor, N, minLambda, maxLambda, gradB, gradN, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
}

void optimizeParameters3(double& Bfactor, double& N, double& fractionSize, double rmax, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double** Bn, double** sintheta, int** thetaIndex, double** concentrations, double*** Inu, double*** Anu, double** area, double** length, double* Rho, double* Phi, FILE* logFile) {
	double* totalInu = new double[Nnu];
	double currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g n = %g fraction = %g\n", Bfactor, N, fractionSize);
	fprintf(logFile, "Bfactor = %g n = %g fraction = %g\n", Bfactor, N, fractionSize);
	for(int i = 0; i < Niterations; ++i) {
		///randomization;
		for(int j = 0; j < 5; ++j){
			double tempN = N + 0.2*N*(uniformDistribution() - 0.5);
			double tempB = Bfactor + 0.2*Bfactor*(uniformDistribution() - 0.5);
			double tempS  = min(max(minFraction,fractionSize + 0.2*fractionSize*(uniformDistribution() - 0.5)), maxFraction);
			double tempF = evaluateOptimizationFunction(tempB, tempN, tempS, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			if(tempF < currentF){
				currentF = tempF;
				Bfactor = tempB;
				fractionSize = tempS;
				N = tempN;
				printf("random search\n");
				fprintf(logFile, "random search\n");
			}
		}
		double prevF = currentF;
		//
		//valley first step
		double valleyB1 = Bfactor;
		double valleyN1 = N;
		double valleyS1 = fractionSize;

		double dxB = fabs(Bfactor)/20;
		double dxN = fabs(N)/20;
		double dxS = fabs(fractionSize)/20;
		printf("optimization i = %d\n",i);
		fprintf(logFile, "optimization i = %d\n",i);
		fflush(logFile);

		double gradB;
		double gradN;
		double gradS;

		double Fb = evaluateOptimizationFunction(Bfactor + dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		double Fb1 = evaluateOptimizationFunction(Bfactor - dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fb > currentF && Fb1 > currentF){
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction(Bfactor + dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fb1 = evaluateOptimizationFunction(Bfactor - dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		double Fn = evaluateOptimizationFunction(Bfactor, N + dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		double Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fn > currentF && Fn1 > currentF){
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction(Bfactor, N + dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		double Fs = evaluateOptimizationFunction(Bfactor, N, fractionSize + dxS, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		double Fs1 = evaluateOptimizationFunction(Bfactor, N, fractionSize - dxS, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradS = (Fs - Fs1)/(2*dxS);
		gradS = (Fs - Fs1)/(2*dxS);
		if(fractionSize >= maxFraction){
			if(gradS > 0){
				gradS = 0;
			}
		}
		if(fractionSize <= minFraction){
			if(gradS < 0){
				gradS = 0;
			}
		}

		double gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS);
		gradB = gradB/gradNorm;
		gradN = gradN/gradNorm;
		gradS = gradS/gradNorm;

		findMinParameters3(Bfactor, N, fractionSize, rmax, gradB, gradN, gradS, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, currentF, totalInu);

		currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		//valley second step

		dxB = fabs(Bfactor)/20;
		dxN = fabs(N)/20;
		dxS = fabs(fractionSize)/20;

		Fb = evaluateOptimizationFunction(Bfactor + dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		Fb1 = evaluateOptimizationFunction(Bfactor - dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fb > currentF && Fb1 > currentF){
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction(Bfactor + dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fb1 = evaluateOptimizationFunction(Bfactor - dxB, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		Fn = evaluateOptimizationFunction(Bfactor, N + dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fn > currentF && Fn1 > currentF){
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction(Bfactor, N + dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		Fs = evaluateOptimizationFunction(Bfactor, N, fractionSize + dxS, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		Fs1 = evaluateOptimizationFunction(Bfactor, N, fractionSize - dxS, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradS = (Fs - Fs1)/(2*dxS);
		gradS = (Fs - Fs1)/(2*dxS);
		if(fractionSize >= maxFraction){
			if(gradS > 0){
				gradS = 0;
			}
		}
		if(fractionSize <= minFraction){
			if(gradS < 0){
				gradS = 0;
			}
		}

		gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS);
		gradB = gradB/gradNorm;
		gradN = gradN/gradNorm;
		gradS = gradS/gradNorm;

		findMinParameters3(Bfactor, N, fractionSize, rmax, gradB, gradN, gradS, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, currentF, totalInu);

		currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		double valleyB2 = Bfactor;
		double valleyN2 = N;
		double valleyS2 = fractionSize;

		//// valley third step

		
		gradB = (valleyB1 - valleyB2);
		gradN = (valleyN1 - valleyN2);
		gradS = (valleyS1 - valleyS2);
		gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS);
		if(gradNorm > 0){
			gradB = gradB/gradNorm;
			gradN = gradN/gradNorm;
			gradS = gradS/gradNorm;

			double step = 0.1*min3(Bfactor, N, fractionSize);
			double Fv = evaluateOptimizationFunction(Bfactor - gradB*step, N - gradN*step, fractionSize - gradS*step, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		//if(Fv < currentF){
				findMinParameters3(Bfactor, N, fractionSize, rmax, gradB, gradN, gradS, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, currentF, totalInu);
		//} else {
			//findMinParameters3(Bfactor, N, fractionSize, -gradB, -gradN, gradS, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, Z, currentF, totalInu);
		//}
		}
		currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		
		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/

		printf("optimization function = %g\n", currentF);
		fprintf(logFile, "optimization function = %g\n", currentF);
		printf("Bfactor = %g n = %g\n fraction = %g\n", Bfactor, N, fractionSize);
		fprintf(logFile, "Bfactor = %g n = %g fraction = %g\n", Bfactor, N, fractionSize);
	}
	printf("finish optimization\n");
	fprintf(logFile, "finish optimization\n");
	delete[] totalInu;
}

void findMinParameters4(const double& B0, const double& N0, const double& R0, double& Bfactor, double& N, double& fractionSize, double& rmax, double gradB, double gradN, double gradS, double gradR, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double** Bn, double** sintheta, int** thetaIndex, double** concentrations, double*** Inu, double*** Anu, double** area, double** length, double* Rho, double* Phi, double& currentF, double* totalInu) {
	double step = 0.4;
	if(fabs(gradB) > 0){
		step = 0.4*fabs(Bfactor/(B0*gradB));
	}
	if(fabs(gradN) > 0){
		step = min(step, 0.4*fabs(N/(N0*gradN)));
	}
	if(gradS > 0){
		step = min(step, 0.4*fabs((maxFraction - fractionSize)/gradS));
	}
	if(gradS < 0){
		step = min(step, 0.4*fabs((fractionSize - minFraction)/gradS));
	}
	if(fabs(gradR) > 0){
		step = min(step, 0.4*fabs(rmax/(R0*gradR)));
	}
	
	double B1 = Bfactor - gradB*step*B0;
	double N1 = N - gradN*step*N0;
	double S1 = fractionSize - gradS*step;
	double r1 = rmax - gradR*step*R0;
	if(B1 != B1){
		printf("B1 Nan\n");
		exit(0);
	}
	if(r1 != r1){
		printf("eta1 Nan\n");
		exit(0);
	}

	if(S1 > maxFraction){
		S1 = maxFraction;
		//printf("S1 > 1\n");
	}
	if(S1 <= minFraction){
		S1 = minFraction;
		//printf("S1 < 0\n");
	}

	double f1 = evaluateOptimizationFunction(B1, N1, S1, r1, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	if(f1 > currentF){
		while (f1 > currentF){
			step = step/2;
			B1 = Bfactor - gradB*step*B0;
			N1 = N - gradN*step*N0;
			S1 = fractionSize - gradS*step;
			r1 = rmax - gradR*step*R0;
			if(B1 != B1){
				printf("B1 Nan\n");
				exit(0);
			}
			if(S1 > maxFraction){
				//printf("S1 > 1\n");
				S1 = maxFraction;
			}
			if(S1 <= minFraction){
				//printf("S1 < 0\n");
				S1 = minFraction;
			}

			f1 = evaluateOptimizationFunction(B1, N1, S1, r1, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		Bfactor = B1;
		N = N1;
		fractionSize = S1;
		rmax = r1;
		return;
	}
	Bfactor = B1;
	N = N1;
	fractionSize = S1;
	rmax = r1;
	step = 0.4*min3(fabs(Bfactor/gradB), fabs(N/gradN), fabs(fractionSize/gradS));
	step = min(step, 0.4*fabs(rmax/gradR));
	double B2 = B1 - gradB*step*B0;
	double N2 = N1 - gradN*step*N0;
	double S2 = S1 - gradS*step;
	double r2 = r1  - gradR*step*R0;
	if(B2 != B2){
		printf("B2 Nan\n");
		exit(0);
	}
	if(S2 > maxFraction){
		printf("S2 > 1\n");
		S2 = maxFraction;
	}
	if(S2 < minFraction){
		S2 = minFraction;
	}
	double f2 = evaluateOptimizationFunction(B2, N2, S2, r2, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	int iterations = 0;
	while(f2 < f1){
		iterations++;
		Bfactor = B2;
		N = N2;
		fractionSize = S2;
		rmax = r2;
		step = 0.4*min4(fabs(Bfactor/gradB), fabs(N/gradN), fabs(fractionSize/gradS), fabs((1.0 - fractionSize)/gradS));
		step = min(step, 0.4*fabs(rmax/gradR));
		B2 = Bfactor - gradB*step*B0;
		N2 = N - gradN*step*N0;
		S2 = fractionSize - gradS*step;
		r2 = rmax - gradR*step*R0;
		if(B2 != B2){
			printf("B2 Nan\n");
			exit(0);
		}
		if(S2 > maxFraction){
			printf("S2 > 1\n");
			S2 = maxFraction;
			return;
		}
		if(S2 < minFraction){
			S2 = minFraction;
			//printf("S2 < 0\n");
			return;
		}
		f1 = f2;
		f2 = evaluateOptimizationFunction(B2, N2, S2, r2, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	}

	//findMinParameters(Bfactor, N, minLambda, maxLambda, gradB, gradN, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
}

void optimizeParameters4(const double& B0, const double& N0, const double& R0, double& Bfactor, double& N, double& fractionSize, double& rmax, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double** Bn, double** sintheta, int** thetaIndex, double** concentrations, double*** Inu, double*** Anu, double** area, double** length, double* Rho, double* Phi, FILE* logFile) {
	double* totalInu = new double[Nnu];
	double currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g n = %g fraction = %g rmax = %g\n", Bfactor, N, fractionSize, rmax);
	fprintf(logFile, "Bfactor = %g n = %g fraction = %g rmax = %g\n", Bfactor, N, fractionSize, rmax);
	for(int i = 0; i < Niterations; ++i) {
		///randomization;
		for(int j = 0; j < 5; ++j){
			double tempN = N + 0.2*N*(uniformDistribution() - 0.5);
			double tempB = Bfactor + 0.2*Bfactor*(uniformDistribution() - 0.5);
			double tempS  = min(max(minFraction,fractionSize + 0.2*fractionSize*(uniformDistribution() - 0.5)), maxFraction);
			double tempF = evaluateOptimizationFunction(tempB, tempN, tempS, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			if(tempF < currentF){
				currentF = tempF;
				Bfactor = tempB;
				fractionSize = tempS;
				N = tempN;
				printf("random search\n");
				fprintf(logFile, "random search\n");
			}
		}
		double prevF = currentF;
		//
		//valley first step
		double valleyB1 = Bfactor;
		double valleyN1 = N;
		double valleyS1 = fractionSize;
		double valleyR1 = rmax;

		double dxB = fabs(Bfactor/B0)/20;
		double dxN = fabs(N/N0)/20;
		double dxS = fabs(fractionSize)/20;
		double dxR = fabs(rmax/R0)/20;
		printf("optimization i = %d\n",i);
		fprintf(logFile, "optimization i = %d\n",i);
		fflush(logFile);

		double gradB;
		double gradN;
		double gradS;
		double gradR;

		double Fb = evaluateOptimizationFunction(Bfactor + dxB*B0, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		double Fb1 = evaluateOptimizationFunction(Bfactor - dxB*B0, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fb > currentF && Fb1 > currentF){
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction(Bfactor + dxB*B0, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fb1 = evaluateOptimizationFunction(Bfactor - dxB*B0, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		double Fn = evaluateOptimizationFunction(Bfactor, N + dxN*N0, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		double Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN*N0, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fn > currentF && Fn1 > currentF){
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction(Bfactor, N + dxN*N0, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN*N0, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		double Fs = evaluateOptimizationFunction(Bfactor, N, fractionSize + dxS, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		double Fs1 = evaluateOptimizationFunction(Bfactor, N, fractionSize - dxS, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradS = (Fs - Fs1)/(2*dxS);
		if(fractionSize >= maxFraction){
			if(gradS > 0){
				gradS = 0;
			}
		}
		if(fractionSize <= minFraction){
			if(gradS < 0){
				gradS = 0;
			}
		}

		double Fr = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax + dxR*R0, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - currentF)/(dxR);

		double gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS + gradR*gradR);
		gradB = gradB/gradNorm;
		gradN = gradN/gradNorm;
		gradS = gradS/gradNorm;
		gradR = gradR/gradNorm;

		findMinParameters4(B0, N0, R0, Bfactor, N, fractionSize, rmax, gradB, gradN, gradS, gradR, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, currentF, totalInu);

		currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		//valley second step

		dxB = fabs(Bfactor/B0)/20;
		dxN = fabs(N/N0)/20;
		dxS = fabs(fractionSize)/20;
		dxR = fabs(rmax/R0)/20;

		Fb = evaluateOptimizationFunction(Bfactor + dxB*B0, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		Fb1 = evaluateOptimizationFunction(Bfactor - dxB*B0, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fb > currentF && Fb1 > currentF){
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction(Bfactor + dxB*B0, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fb1 = evaluateOptimizationFunction(Bfactor - dxB*B0, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		Fn = evaluateOptimizationFunction(Bfactor, N + dxN*N0, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN*N0, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		while( Fn > currentF && Fn1 > currentF){
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction(Bfactor, N + dxN*N0, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
			Fn1 = evaluateOptimizationFunction(Bfactor, N - dxN*N0, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		Fs = evaluateOptimizationFunction(Bfactor, N, fractionSize + dxS, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		Fs1 = evaluateOptimizationFunction(Bfactor, N, fractionSize - dxS, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradS = (Fs - Fs1)/(2*dxS);
		if(fractionSize >= maxFraction){
			if(gradS > 0){
				gradS = 0;
			}
		}
		if(fractionSize <= minFraction){
			if(gradS < 0){
				gradS = 0;
			}
		}

		Fr = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax + dxR*R0, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		//Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - currentF)/(dxR);

		gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS + gradR*gradR);
		gradB = gradB/gradNorm;
		gradN = gradN/gradNorm;
		gradS = gradS/gradNorm;
		gradR = gradR/gradNorm;

		findMinParameters4(B0, N0, R0, Bfactor, N, fractionSize, rmax, gradB, gradN, gradS, gradR, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, currentF, totalInu);

		currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		double valleyB2 = Bfactor;
		double valleyN2 = N;
		double valleyS2 = fractionSize;
		double valleyR2 = rmax;

		//// valley third step

		
		gradB = (valleyB1 - valleyB2)/B0;
		gradN = (valleyN1 - valleyN2)/N0;
		gradS = (valleyS1 - valleyS2);
		gradR = (valleyR1 - valleyR2)/R0;
		gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS + gradR*gradR);
		if(gradNorm > 0){
			gradB = gradB/gradNorm;
			gradN = gradN/gradNorm;
			gradS = gradS/gradNorm;
			gradR = gradR/gradNorm;

			//double step = 0.1*min3(Bfactor, N, fractionSize);
			//double Fv = evaluateOptimizationFunction(Bfactor - gradB*step, N - gradN*step, fractionSize - gradS*step, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		//if(Fv < currentF){
				findMinParameters4(B0, N0, R0, Bfactor, N, fractionSize, rmax, gradB, gradN, gradS, gradR, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, currentF, totalInu);
		//} else {
			//findMinParameters3(Bfactor, N, fractionSize, -gradB, -gradN, gradS, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, Z, currentF, totalInu);
		//}
		}
		currentF = evaluateOptimizationFunction(Bfactor, N, fractionSize, rmax, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		
		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/

		printf("optimization function = %g\n", currentF);
		fprintf(logFile, "optimization function = %g\n", currentF);
		printf("Bfactor = %g n = %g fraction = %g rmax = %g\n", Bfactor, N, fractionSize, rmax);
		fprintf(logFile, "Bfactor = %g n = %g fraction = %g rmax = %g\n", Bfactor, N, fractionSize, rmax);
	}
	printf("finish optimization\n");
	fprintf(logFile, "finish optimization\n");
	delete[] totalInu;
}

//for fla disc

double evaluateOptimizationFunction5(double Bfactor, double n, double fractionSize, double rmax, double v, double** nu, double** F, double** Ee, double** Fe, int Np, int Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double***** Inu, double***** Anu, double*** area, double*** length){
	if(fractionSize > 1.0){
		fractionSize = 1.0;
	}
	if(v > maxV){
		v = maxV;
	}
	double err = 0;
	double r = rmax;

	double* totalInu = new double[Nnu];
	//evaluateVolumeAndLength(area, length, rmax, fractionSize);
	for(int i = 0; i < Nmonth; ++i){
		r = rmax + v*times[i];
		double rfactor = r/rmax;
		evaluateAllEmissivityAndAbsorption(nu[i], Inu[i], Anu[i], Nnu, Ee, Fe, Np, Nd, Bn, sintheta, thetaIndex, concentrations, n, Bfactor, rfactor);
		//
		//evaluateSpectrum(nu[i], totalInu, Inu[i], Anu[i], area, length, Nnu, rfactor);
		//
		if(geometry == SPHERICAL){
			evaluateSpectrumSpherical(nu[i], totalInu, Inu[i], Anu[i], rmax, Nnu, rfactor, fractionSize);
		} else {
			evaluateSpectrumFlat(nu[i], totalInu, Inu[i], Anu[i], rmax, Nnu, rfactor, fractionSize);
		}
		for(int j = 0; j < Nnu; ++j){
			double err1 = 0;
			if(scale == LINEAR){
				err1 = sqr(totalInu[j] - F[i][j]);
			} else {
				err1 = sqr(log(totalInu[j]) - log(F[i][j]));
			}
			err = err + err1;
		}
	}
	delete[] totalInu;
	
	return err;
}

void findMinParameters5(const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v, double gradB, double gradN, double gradS, double gradR, double gradV, double** nu, double** F, double** Ee, double** Fe, int Np, int Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double***** Inu, double***** Anu, double*** area, double*** length, double& currentF){
	double step = 0.4;
	if(fabs(gradB) > 0){
		step = 0.4*fabs(Bfactor/(B0*gradB));
	}
	if(fabs(gradN) > 0){
		step = min(step, 0.4*fabs(N/(N0*gradN)));
	}
	if(gradS < 0){
		step = min(step, 0.4*fabs((maxFraction - fractionSize)/gradS));
	}
	if(gradS > 0){
		step = min(step, 0.4*fabs((fractionSize - minFraction)/gradS));
	}
	if(fabs(gradR) > 0){
		step = min(step, 0.4*fabs(rmax/(R0*gradR)));
	}
	if(gradV < 0){
		step = min(step, 0.4*fabs((maxV - v)/(V0*gradV)));
	}
	if(gradV > 0){
		step = min(step, 0.4*fabs(v/(V0*gradV)));
	}

	
	double B1 = Bfactor - gradB*step*B0;
	double N1 = N - gradN*step*N0;
	double S1 = fractionSize - gradS*step;
	double r1 = rmax - gradR*step*R0;
	double v1 = v - gradV*step*V0;
	if(B1 != B1){
		printf("B1 Nan\n");
		exit(0);
	}
	if(r1 != r1){
		printf("eta1 Nan\n");
		exit(0);
	}

	if(S1 > maxFraction){
		S1 = maxFraction;
		//printf("S1 > 1\n");
	}
	if(S1 <= minFraction){
		S1 = minFraction;
		//printf("S1 < 0\n");
	}
	if(v1 > maxV){
		v1 = maxV;
	}

	double f1 = evaluateOptimizationFunction5(B1, N1, S1, r1, v1, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	if(f1 > currentF){
		int count = 0;
		while (f1 > currentF && count < 10){
			count++;
			step = step/2;
			B1 = Bfactor - gradB*step*B0;
			N1 = N - gradN*step*N0;
			S1 = fractionSize - gradS*step;
			r1 = rmax - gradR*step*R0;
			v1 = v - gradV*step*V0;
			if(B1 != B1){
				printf("B1 Nan\n");
				exit(0);
			}
			if(S1 > maxFraction){
				//printf("S1 > 1\n");
				S1 = maxFraction;
			}
			if(S1 <= minFraction){
				//printf("S1 < 0\n");
				S1 = minFraction;
			}
			if(v1 > maxV){
				v1 = maxV;
			}

			f1 = evaluateOptimizationFunction5(B1, N1, S1, r1, v1, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		if(f1 < currentF){
			Bfactor = B1;
			N = N1;
			fractionSize = S1;
			rmax = r1;
			v = v1;
		}
		return;
	}
	Bfactor = B1;
	N = N1;
	fractionSize = S1;
	rmax = r1;
	v = v1;
	step = 0.4;
	if(fabs(gradB) > 0){
		step = 0.4*fabs(Bfactor/(B0*gradB));
	}
	if(fabs(gradN) > 0){
		step = min(step, 0.4*fabs(N/(N0*gradN)));
	}
	if(gradS < 0){
		step = min(step, 0.4*fabs((maxFraction - fractionSize)/gradS));
	}
	if(gradS > 0){
		step = min(step, 0.4*fabs((fractionSize - minFraction)/gradS));
	}
	if(fabs(gradR) > 0){
		step = min(step, 0.4*fabs(rmax/(R0*gradR)));
	}
	if(gradV < 0){
		step = min(step, 0.4*fabs((maxV - v)/(V0*gradV)));
	}
	if(gradV > 0){
		step = min(step, 0.4*fabs(v/(V0*gradV)));
	}
	double B2 = B1 - gradB*step*B0;
	double N2 = N1 - gradN*step*N0;
	double S2 = S1 - gradS*step;
	double r2 = r1  - gradR*step*R0;
	double v2 = v1 - gradV*step*V0;
	if(B2 != B2){
		printf("B2 Nan\n");
		exit(0);
	}
	if(S2 > maxFraction){
		printf("S2 > 1\n");
		S2 = maxFraction;
	}
	if(S2 < minFraction){
		S2 = minFraction;
	}
	if(v2 > maxV){
		v2 = maxV;
	}
	double f2 = evaluateOptimizationFunction5(B2, N2, S2, r2, v2, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	int iterations = 0;
	while((f2 < f1) && (iterations < 100)){
		iterations++;
		Bfactor = B2;
		N = N2;
		fractionSize = S2;
		rmax = r2;
		v = v2;
		step = 0.4;
		if(fabs(gradB) > 0){
			step = 0.4*fabs(Bfactor/(B0*gradB));
		}
		if(fabs(gradN) > 0){
			step = min(step, 0.4*fabs(N/(N0*gradN)));
		}
		if(gradS < 0){
			step = min(step, 0.4*fabs((maxFraction - fractionSize)/gradS));
		}
		if(gradS > 0){
			step = min(step, 0.4*fabs((fractionSize - minFraction)/gradS));
		}
		if(fabs(gradR) > 0){
			step = min(step, 0.4*fabs(rmax/(R0*gradR)));
		}
		if(gradV < 0){
			step = min(step, 0.4*fabs((maxV - v)/(V0*gradV)));
		}
		if(gradV > 0){
			step = min(step, 0.4*fabs(v/(V0*gradV)));
		}
		B2 = Bfactor - gradB*step*B0;
		N2 = N - gradN*step*N0;
		S2 = fractionSize - gradS*step;
		r2 = rmax - gradR*step*R0;
		v2 = v - gradV*step*V0;
		if(B2 != B2){
			printf("B2 Nan\n");
			exit(0);
		}
		if(S2 > maxFraction){
			printf("S2 > 1\n");
			S2 = maxFraction;
			return;
		}
		if(S2 < minFraction){
			S2 = minFraction;
			//printf("S2 < 0\n");
			return;
		}
		if(v2 > maxV){
			v2 = maxV;
		}
		f1 = f2;
		f2 = evaluateOptimizationFunction5(B2, N2, S2, r2, v2, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	}
}

void optimizeParameters5(const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v,  double** nu, double** F, double** Ee, double** Fe, int Np, int Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double***** Inu, double***** Anu, double*** area, double*** length, FILE* logFile){
	double currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g n = %g fraction = %g rmax = %g v = %g\n", Bfactor, N, fractionSize, rmax, v);
	fprintf(logFile, "Bfactor = %g n = %g fraction = %g rmax = %g v = %g\n", Bfactor, N, fractionSize, rmax, v);
	for(int i = 0; i < Niterations; ++i) {
		///randomization;
		for(int j = 0; j < 5; ++j){
			double tempN = N + 0.2*N*(uniformDistribution() - 0.5);
			//tempN = N;
			double tempB = Bfactor + 0.2*Bfactor*(uniformDistribution() - 0.5);
			//tempB = Bfactor;
			double tempS  = min(max(minFraction,fractionSize + 0.2*fractionSize*(uniformDistribution() - 0.5)), maxFraction);
			double tempR = rmax + 0.2*rmax*(uniformDistribution() - 0.5);
			double tempV = 0.2*speed_of_light + (maxV - 0.2*speed_of_light)*uniformDistribution();
			double tempF = evaluateOptimizationFunction5(tempB, tempN, tempS, tempR, tempV, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			if(tempF < currentF){
				currentF = tempF;
				Bfactor = tempB;
				fractionSize = tempS;
				N = tempN;
				rmax = tempR;
				v = tempV;
				printf("random search\n");
				fprintf(logFile, "random search\n");
			}
		}
		double prevF = currentF;
		//
		//valley first step
		double valleyB1 = Bfactor;
		double valleyN1 = N;
		double valleyS1 = fractionSize;
		double valleyR1 = rmax;
		double valleyV1 = v;

		double dxB = fabs(Bfactor/B0)/20;
		double dxN = fabs(N/N0)/20;
		double dxS = fabs(fractionSize)/20;
		double dxR = fabs(rmax/R0)/20;
		double dxV = fabs(v/V0)/20;
		printf("optimization i = %d\n",i);
		fprintf(logFile, "optimization i = %d\n",i);
		fflush(logFile);

		double gradB;
		double gradN;
		double gradS;
		double gradR;
		double gradV;

		double Fb = evaluateOptimizationFunction5(Bfactor + dxB*B0, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fb1 = evaluateOptimizationFunction5(Bfactor - dxB*B0, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		int count = 0;
		while( Fb > currentF && Fb1 > currentF && count < 5){
			count ++;
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction5(Bfactor + dxB*B0, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fb1 = evaluateOptimizationFunction5(Bfactor - dxB*B0, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		double Fn = evaluateOptimizationFunction5(Bfactor, N + dxN*N0, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fn1 = evaluateOptimizationFunction5(Bfactor, N - dxN*N0, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fn > currentF && Fn1 > currentF && count < 5){
			count++;
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction5(Bfactor, N + dxN*N0, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fn1 = evaluateOptimizationFunction5(Bfactor, N - dxN*N0, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		double Fs = evaluateOptimizationFunction5(Bfactor, N, fractionSize + dxS, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fs1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize - dxS, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		gradS = (Fs - Fs1)/(2*dxS);
		if(fractionSize >= maxFraction){
			if(gradS > 0){
				gradS = 0;
			}
		}
		if(fractionSize <= minFraction){
			if(gradS < 0){
				gradS = 0;
			}
		}

		double Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while( Fr > currentF && Fr1 > currentF && count < 5){
			count++;
			dxR = dxR/2;
			Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - Fr1)/(2*dxR);

		double Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while(Fv > currentF && Fv1 > currentF && count < 5){
			count++;
			dxV = dxV/2;
			Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradV = (Fv - Fv1)/(2*dxV);
		if(v >= maxV){
			if(gradV > 0){
				gradV = 0;
			}
		}

		//gradB = 0;
		//gradN = 0;

		double gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS + gradR*gradR + gradV*gradV);
		if(gradNorm != gradNorm){
			printf("gradNorm = NaN\n");
			exit(0);
		}
		gradB = gradB/gradNorm;
		gradN = gradN/gradNorm;
		gradS = gradS/gradNorm;
		gradR = gradR/gradNorm;
		gradV = gradV/gradNorm;

		findMinParameters5(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);

		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		//valley second step

		dxB = fabs(Bfactor/B0)/20;
		dxN = fabs(N/N0)/20;
		dxS = fabs(fractionSize)/20;
		dxR = fabs(rmax/R0)/20;
		dxV = fabs(v/V0)/20;

		fflush(logFile);

		Fb = evaluateOptimizationFunction5(Bfactor + dxB*B0, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fb1 = evaluateOptimizationFunction5(Bfactor - dxB*B0, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fb > currentF && Fb1 > currentF && count < 5){
			count++;
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction5(Bfactor + dxB*B0, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fb1 = evaluateOptimizationFunction5(Bfactor - dxB*B0, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		Fn = evaluateOptimizationFunction5(Bfactor, N + dxN*N0, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fn1 = evaluateOptimizationFunction5(Bfactor, N - dxN*N0, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fn > currentF && Fn1 > currentF && count < 5){
			count++;
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction5(Bfactor, N + dxN*N0, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fn1 = evaluateOptimizationFunction5(Bfactor, N - dxN*N0, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		Fs = evaluateOptimizationFunction5(Bfactor, N, fractionSize + dxS, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fs1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize - dxS, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		gradS = (Fs - Fs1)/(2*dxS);
		if(fractionSize >= maxFraction){
			if(gradS > 0){
				gradS = 0;
			}
		}
		if(fractionSize <= minFraction){
			if(gradS < 0){
				gradS = 0;
			}
		}

		Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fr > currentF && Fr1 > currentF && count < 5){
			count++;
			dxR = dxR/2;
			Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - Fr1)/(2*dxR);

		Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while(Fv > currentF && Fv1 > currentF && count < 5){
			count++;
			dxV = dxV/2;
			Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradV = (Fv - Fv1)/(2*dxV);
		if(v >= maxV){
			if(gradV > 0){
				gradV = 0;
			}
		}

		//gradB = 0;
		//gradN = 0;

		gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS + gradR*gradR + gradV*gradV);
		gradB = gradB/gradNorm;
		gradN = gradN/gradNorm;
		gradS = gradS/gradNorm;
		gradR = gradR/gradNorm;
		gradV = gradV/gradNorm;

		findMinParameters5(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);

		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		double valleyB2 = Bfactor;
		double valleyN2 = N;
		double valleyS2 = fractionSize;
		double valleyR2 = rmax;
		double valleyV2 = v;

		//// valley third step

		
		gradB = (valleyB1 - valleyB2)/B0;
		gradN = (valleyN1 - valleyN2)/N0;
		gradS = (valleyS1 - valleyS2);
		gradR = (valleyR1 - valleyR2)/R0;
		gradV = (valleyV1 - valleyV2)/V0;

		//gradB = 0;
		//gradN = 0;
		gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS + gradR*gradR + gradV*gradV);
		if(gradNorm > 0){
			gradB = gradB/gradNorm;
			gradN = gradN/gradNorm;
			gradS = gradS/gradNorm;
			gradR = gradR/gradNorm;
			gradV = gradV/gradNorm;
			//double step = 0.1*min3(Bfactor, N, fractionSize);
			//double Fv = evaluateOptimizationFunction(Bfactor - gradB*step, N - gradN*step, fractionSize - gradS*step, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		//if(Fv < currentF){
				findMinParameters5(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);
		//} else {
			//findMinParameters3(Bfactor, N, fractionSize, -gradB, -gradN, gradS, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, Z, currentF, totalInu);
		//}
		}
		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		
		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/

		printf("optimization function = %g\n", currentF);
		fprintf(logFile, "optimization function = %g\n", currentF);
		printf("Bfactor = %g n = %g fraction = %g rmax = %g v/c = %g\n", Bfactor, N, fractionSize, rmax, v/speed_of_light);
		fprintf(logFile, "Bfactor = %g n = %g fraction = %g rmax = %g v/c =%g\n", Bfactor, N, fractionSize, rmax, v/speed_of_light);
	}
	printf("finish optimization\n");
	fprintf(logFile, "finish optimization\n");
}

void findMinParameters5sigma(const double& sigma, const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v, double gradB, double gradS, double gradR, double gradV, double** nu, double** F, double** Ee, double** Fe, int Np, int Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double***** Inu, double***** Anu, double*** area, double*** length, double& currentF){
	double step = 0.4;
	if(fabs(gradB) > 0){
		step = 0.4*fabs(Bfactor/(B0*gradB));
	}
	if(gradS < 0){
		step = min(step, 0.4*fabs((maxFraction - fractionSize)/gradS));
	}
	if(gradS > 0){
		step = min(step, 0.4*fabs((fractionSize - minFraction)/gradS));
	}
	if(fabs(gradR) > 0){
		step = min(step, 0.4*fabs(rmax/(R0*gradR)));
	}
	if(gradV < 0){
		step = min(step, 0.4*fabs((maxV - v)/(V0*gradV)));
	}
	if(gradV > 0){
		step = min(step, 0.4*fabs(v/(V0*gradV)));
	}

	
	double B1 = Bfactor - gradB*step*B0;
	double N1 = sqr(B1)/(sigma*4*pi*massProtonReal*speed_of_light2);
	double S1 = fractionSize - gradS*step;
	double r1 = rmax - gradR*step*R0;
	double v1 = v - gradV*step*V0;
	if(B1 != B1){
		printf("B1 Nan\n");
		exit(0);
	}
	if(r1 != r1){
		printf("eta1 Nan\n");
		exit(0);
	}

	if(S1 > maxFraction){
		S1 = maxFraction;
		//printf("S1 > 1\n");
	}
	if(S1 <= minFraction){
		S1 = minFraction;
		//printf("S1 < 0\n");
	}
	if(v1 > maxV){
		v1 = maxV;
	}

	double f1 = evaluateOptimizationFunction5(B1, N1, S1, r1, v1, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	if(f1 > currentF){
		int count = 0;
		while (f1 > currentF && count < 10){
			count++;
			step = step/2;
			B1 = Bfactor - gradB*step*B0;
			N1 = sqr(B1)/(sigma*4*pi*massProtonReal*speed_of_light2);
			S1 = fractionSize - gradS*step;
			r1 = rmax - gradR*step*R0;
			v1 = v - gradV*step*V0;
			if(B1 != B1){
				printf("B1 Nan\n");
				exit(0);
			}
			if(S1 > maxFraction){
				//printf("S1 > 1\n");
				S1 = maxFraction;
			}
			if(S1 <= minFraction){
				//printf("S1 < 0\n");
				S1 = minFraction;
			}
			if(v1 > maxV){
				v1 = maxV;
			}

			f1 = evaluateOptimizationFunction5(B1, N1, S1, r1, v1, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		if(f1 < currentF){
			Bfactor = B1;
			N = N1;
			fractionSize = S1;
			rmax = r1;
			v = v1;
		}
		return;
	}
	Bfactor = B1;
	N = N1;
	fractionSize = S1;
	rmax = r1;
	v = v1;
	step = 0.4;
	if(fabs(gradB) > 0){
		step = 0.4*fabs(Bfactor/(B0*gradB));
	}
	if(gradS < 0){
		step = min(step, 0.4*fabs((maxFraction - fractionSize)/gradS));
	}
	if(gradS > 0){
		step = min(step, 0.4*fabs((fractionSize - minFraction)/gradS));
	}
	if(fabs(gradR) > 0){
		step = min(step, 0.4*fabs(rmax/(R0*gradR)));
	}
	if(gradV < 0){
		step = min(step, 0.4*fabs((maxV - v)/(V0*gradV)));
	}
	if(gradV > 0){
		step = min(step, 0.4*fabs(v/(V0*gradV)));
	}
	double B2 = B1 - gradB*step*B0;
	double N2 = sqr(B2)/(sigma*4*pi*massProtonReal*speed_of_light2);
	double S2 = S1 - gradS*step;
	double r2 = r1  - gradR*step*R0;
	double v2 = v1 - gradV*step*V0;
	if(B2 != B2){
		printf("B2 Nan\n");
		exit(0);
	}
	if(S2 > maxFraction){
		printf("S2 > 1\n");
		S2 = maxFraction;
	}
	if(S2 < minFraction){
		S2 = minFraction;
	}
	if(v2 > maxV){
		v2 = maxV;
	}
	double f2 = evaluateOptimizationFunction5(B2, N2, S2, r2, v2, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	int iterations = 0;
	while((f2 < f1) && (iterations < 100)){
		iterations++;
		Bfactor = B2;
		N = N2;
		fractionSize = S2;
		rmax = r2;
		v = v2;
		step = 0.4;
		if(fabs(gradB) > 0){
			step = 0.4*fabs(Bfactor/(B0*gradB));
		}
		if(gradS < 0){
			step = min(step, 0.4*fabs((maxFraction - fractionSize)/gradS));
		}
		if(gradS > 0){
			step = min(step, 0.4*fabs((fractionSize - minFraction)/gradS));
		}
		if(fabs(gradR) > 0){
			step = min(step, 0.4*fabs(rmax/(R0*gradR)));
		}
		if(gradV < 0){
			step = min(step, 0.4*fabs((maxV - v)/(V0*gradV)));
		}
		if(gradV > 0){
			step = min(step, 0.4*fabs(v/(V0*gradV)));
		}
		B2 = Bfactor - gradB*step*B0;
		N2 = sqr(B2)/(sigma*4*pi*massProtonReal*speed_of_light2);
		S2 = fractionSize - gradS*step;
		r2 = rmax - gradR*step*R0;
		v2 = v - gradV*step*V0;
		if(B2 != B2){
			printf("B2 Nan\n");
			exit(0);
		}
		if(S2 > maxFraction){
			printf("S2 > 1\n");
			S2 = maxFraction;
			return;
		}
		if(S2 < minFraction){
			S2 = minFraction;
			//printf("S2 < 0\n");
			return;
		}
		if(v2 > maxV){
			v2 = maxV;
		}
		f1 = f2;
		f2 = evaluateOptimizationFunction5(B2, N2, S2, r2, v2, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	}
}

void optimizeParameters5sigma(const double& sigma, const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v,  double** nu, double** F, double** Ee, double** Fe, int Np, int Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double***** Inu, double***** Anu, double*** area, double*** length, FILE* logFile){
	double currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g n = %g fraction = %g rmax = %g v = %g\n", Bfactor, N, fractionSize, rmax, v);
	fprintf(logFile, "Bfactor = %g n = %g fraction = %g rmax = %g v = %g\n", Bfactor, N, fractionSize, rmax, v);
	for(int i = 0; i < Niterations; ++i) {
		///randomization;
		for(int j = 0; j < 5; ++j){
			//double tempN = N + 0.2*N*(uniformDistribution() - 0.5);
			//tempN = N;
			double tempB = Bfactor + 0.2*Bfactor*(uniformDistribution() - 0.5);
			double tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
			//tempB = Bfactor;
			double tempS  = min(max(minFraction,fractionSize + 0.2*fractionSize*(uniformDistribution() - 0.5)), maxFraction);
			double tempR = rmax + 0.2*rmax*(uniformDistribution() - 0.5);
			double tempV = 0.2*speed_of_light + (maxV - 0.2*speed_of_light)*uniformDistribution();
			double tempF = evaluateOptimizationFunction5(tempB, tempN, tempS, tempR, tempV, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			if(tempF < currentF){
				currentF = tempF;
				Bfactor = tempB;
				fractionSize = tempS;
				N = tempN;
				rmax = tempR;
				v = tempV;
				printf("random search\n");
				fprintf(logFile, "random search\n");
			}
		}
		double prevF = currentF;
		//
		//valley first step
		double valleyB1 = Bfactor;
		double valleyS1 = fractionSize;
		double valleyR1 = rmax;
		double valleyV1 = v;

		double dxB = fabs(Bfactor/B0)/20;
		double dxS = fabs(fractionSize)/20;
		double dxR = fabs(rmax/R0)/20;
		double dxV = fabs(v/V0)/20;
		printf("optimization i = %d\n",i);
		fprintf(logFile, "optimization i = %d\n",i);
		fflush(logFile);

		double gradB;
		double gradS;
		double gradR;
		double gradV;

		double tempB = Bfactor + dxB*B0;
		double tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);

		double Fb = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		tempB = Bfactor - dxB*B0;
		tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
		double Fb1 = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		int count = 0;
		while( Fb > currentF && Fb1 > currentF && count < 5){
			count ++;
			dxB = dxB/2;
			tempB = Bfactor + dxB*B0;
			tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
			Fb = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			tempB = Bfactor - dxB*B0;
			tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
			Fb1 = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		double Fs = evaluateOptimizationFunction5(Bfactor, N, fractionSize + dxS, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fs1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize - dxS, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		gradS = (Fs - Fs1)/(2*dxS);
		if(fractionSize >= maxFraction){
			if(gradS > 0){
				gradS = 0;
			}
		}
		if(fractionSize <= minFraction){
			if(gradS < 0){
				gradS = 0;
			}
		}

		double Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while( Fr > currentF && Fr1 > currentF && count < 5){
			count++;
			dxR = dxR/2;
			Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - Fr1)/(2*dxR);

		double Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while(Fv > currentF && Fv1 > currentF && count < 5){
			count++;
			dxV = dxV/2;
			Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradV = (Fv - Fv1)/(2*dxV);
		if(v >= maxV){
			if(gradV > 0){
				gradV = 0;
			}
		}

		//gradB = 0;
		//gradN = 0;

		double gradNorm = sqrt(gradB*gradB + gradS*gradS + gradR*gradR + gradV*gradV);
		if(gradNorm != gradNorm){
			printf("gradNorm = NaN\n");
			exit(0);
		}
		gradB = gradB/gradNorm;
		gradS = gradS/gradNorm;
		gradR = gradR/gradNorm;
		gradV = gradV/gradNorm;

		findMinParameters5sigma(sigma, B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradS, gradR, gradV, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);

		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		//valley second step

		dxB = fabs(Bfactor/B0)/20;
		dxS = fabs(fractionSize)/20;
		dxR = fabs(rmax/R0)/20;
		dxV = fabs(v/V0)/20;

		fflush(logFile);

		tempB = Bfactor + dxB*B0;
		tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
		Fb = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		tempB = Bfactor - dxB*B0;
		tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
		Fb1 = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fb > currentF && Fb1 > currentF && count < 5){
			count++;
			dxB = dxB/2;
			tempB = Bfactor + dxB*B0;
			tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
			Fb = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			tempB = Bfactor - dxB*B0;
			tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
			Fb1 = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradB = (Fb - Fb1)/(2*dxB);


		Fs = evaluateOptimizationFunction5(Bfactor, N, fractionSize + dxS, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fs1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize - dxS, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		gradS = (Fs - Fs1)/(2*dxS);
		if(fractionSize >= maxFraction){
			if(gradS > 0){
				gradS = 0;
			}
		}
		if(fractionSize <= minFraction){
			if(gradS < 0){
				gradS = 0;
			}
		}

		Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fr > currentF && Fr1 > currentF && count < 5){
			count++;
			dxR = dxR/2;
			Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - Fr1)/(2*dxR);

		Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while(Fv > currentF && Fv1 > currentF && count < 5){
			count++;
			dxV = dxV/2;
			Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradV = (Fv - Fv1)/(2*dxV);
		if(v >= maxV){
			if(gradV > 0){
				gradV = 0;
			}
		}

		//gradB = 0;
		//gradN = 0;

		gradNorm = sqrt(gradB*gradB + gradS*gradS + gradR*gradR + gradV*gradV);
		gradB = gradB/gradNorm;
		gradS = gradS/gradNorm;
		gradR = gradR/gradNorm;
		gradV = gradV/gradNorm;

		findMinParameters5sigma(sigma, B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradS, gradR, gradV, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);

		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		double valleyB2 = Bfactor;
		double valleyS2 = fractionSize;
		double valleyR2 = rmax;
		double valleyV2 = v;

		//// valley third step

		
		gradB = (valleyB1 - valleyB2)/B0;
		gradS = (valleyS1 - valleyS2);
		gradR = (valleyR1 - valleyR2)/R0;
		gradV = (valleyV1 - valleyV2)/V0;

		//gradB = 0;
		//gradN = 0;
		gradNorm = sqrt(gradB*gradB + gradS*gradS + gradR*gradR + gradV*gradV);
		if(gradNorm > 0){
			gradB = gradB/gradNorm;
			gradS = gradS/gradNorm;
			gradR = gradR/gradNorm;
			gradV = gradV/gradNorm;
			//double step = 0.1*min3(Bfactor, N, fractionSize);
			//double Fv = evaluateOptimizationFunction(Bfactor - gradB*step, N - gradN*step, fractionSize - gradS*step, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		//if(Fv < currentF){
				findMinParameters5sigma(sigma, B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradS, gradR, gradV, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);
		//} else {
			//findMinParameters3(Bfactor, N, fractionSize, -gradB, -gradN, gradS, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, Z, currentF, totalInu);
		//}
		}
		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, F, Ee, Fe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		
		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/

		printf("optimization function = %g\n", currentF);
		fprintf(logFile, "optimization function = %g\n", currentF);
		printf("Bfactor = %g n = %g fraction = %g rmax = %g v/c = %g\n", Bfactor, N, fractionSize, rmax, v/speed_of_light);
		fprintf(logFile, "Bfactor = %g n = %g fraction = %g rmax = %g v/c =%g\n", Bfactor, N, fractionSize, rmax, v/speed_of_light);
	}
	printf("finish optimization\n");
	fprintf(logFile, "finish optimization\n");
}

double evaluateOptimizationFunction5simple(double Bfactor, double n, double fractionSize, double rmax, double v, double** Ee, double** Fe, int Np, int Nd, double sintheta, int thetaIndex) {
	if(fractionSize > 1.0){
		fractionSize = 1.0;
	}
	if(v > maxV){
		v = maxV;
	}
	double err = 0;
	double nu = 0;
	double Inu = 0;
	double Anu = 0;
	double r = rmax;
	double B = Bfactor;
	double N = n;
	for(int i = 0; i < 4; ++i){
		double I = 0;
		nu = aprx[i]*1E9;
		evaluateEmissivityAndAbsorptionAtNuSimple(nu, Inu, Anu, Ee[thetaIndex], Fe[thetaIndex], Np, sintheta, B, N);
		evaluateSpectrumAtNuSimple(nu, I, Inu, Anu, r, fractionSize);
		err = err + sqr(I - apry[i]);
	}

	r = rmax + v*times[1];
	B = Bfactor*rmax/r;
	N = n*sqr(rmax/r);
	for(int i = 0; i < 3; ++i){
		double I = 0;
		nu = mayx[i]*1E9;
		evaluateEmissivityAndAbsorptionAtNuSimple(nu, Inu, Anu, Ee[thetaIndex], Fe[thetaIndex], Np, sintheta, B, N);
		evaluateSpectrumAtNuSimple(nu, I, Inu, Anu, r, fractionSize);
		err = err + sqr(I - mayy[i]);
	}

	r = rmax + v*times[2];
	B = Bfactor*rmax/r;
	N = n*sqr(rmax/r);
	for(int i = 0; i < 4; ++i){
		double I = 0;
		nu = junx[i]*1E9;
		evaluateEmissivityAndAbsorptionAtNuSimple(nu, Inu, Anu, Ee[thetaIndex], Fe[thetaIndex], Np, sintheta, B, N);
		evaluateSpectrumAtNuSimple(nu, I, Inu, Anu, r, fractionSize);
		err = err + sqr(I - juny[i]);
	}

	r = rmax + v*times[3];
	B = Bfactor*rmax/r;
	N = n*sqr(rmax/r);
	for(int i = 0; i < 5; ++i){
		double I = 0;
		nu = augx[i]*1E9;
		evaluateEmissivityAndAbsorptionAtNuSimple(nu, Inu, Anu, Ee[thetaIndex], Fe[thetaIndex], Np, sintheta, B, N);
		evaluateSpectrumAtNuSimple(nu, I, Inu, Anu, r, fractionSize);
		err = err + sqr(I - augy[i]);
	}

	return err;
}



void findMinParameters5simple(const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v, double gradB, double gradN, double gradS, double gradR, double gradV, double** Ee, double** Fe, int Np, int Nd, double sintheta, int thetaIndex, double& currentF) {
	double step = 0.4;
	if(fabs(gradB) > 0){
		step = 0.4*fabs(Bfactor/(B0*gradB));
	}
	if(fabs(gradN) > 0){
		step = min(step, 0.4*fabs(N/(N0*gradN)));
	}
	if(gradS > 0){
		step = min(step, 0.4*fabs((maxFraction - fractionSize)/gradS));
	}
	if(gradS < 0){
		step = min(step, 0.4*fabs((fractionSize - minFraction)/gradS));
	}
	if(fabs(gradR) > 0){
		step = min(step, 0.4*fabs(rmax/(R0*gradR)));
	}
	if(gradV > 0){
		step = min(step, 0.4*fabs((maxV - v)/(V0*gradV)));
	}
	if(gradV < 0){
		step = min(step, 0.4*fabs(v/(V0*gradV)));
	}

	
	double B1 = Bfactor - gradB*step*B0;
	double N1 = N - gradN*step*N0;
	double S1 = fractionSize - gradS*step;
	double r1 = rmax - gradR*step*R0;
	double v1 = v - gradV*step*V0;
	if(B1 != B1){
		printf("B1 Nan\n");
		exit(0);
	}
	if(r1 != r1){
		printf("eta1 Nan\n");
		exit(0);
	}

	if(S1 > maxFraction){
		S1 = maxFraction;
		//printf("S1 > 1\n");
	}
	if(S1 <= minFraction){
		S1 = minFraction;
		//printf("S1 < 0\n");
	}
	if(v1 > maxV){
		v1 = maxV;
	}

	double f1 = evaluateOptimizationFunction5simple(B1, N1, S1, r1, v1, Ee, Fe, Np, Nd, sintheta, thetaIndex);
	if(f1 > currentF){
		while (f1 > currentF){
			step = step/2;
			B1 = Bfactor - gradB*step*B0;
			N1 = N - gradN*step*N0;
			S1 = fractionSize - gradS*step;
			r1 = rmax - gradR*step*R0;
			v1 = v - gradV*step*V0;
			if(B1 != B1){
				printf("B1 Nan\n");
				exit(0);
			}
			if(S1 > maxFraction){
				//printf("S1 > 1\n");
				S1 = maxFraction;
			}
			if(S1 <= minFraction){
				//printf("S1 < 0\n");
				S1 = minFraction;
			}
			if(v1 > maxV){
				v1 = maxV;
			}

			f1 = evaluateOptimizationFunction5simple(B1, N1, S1, r1, v1, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		}
		Bfactor = B1;
		N = N1;
		fractionSize = S1;
		rmax = r1;
		v = v1;
		return;
	}
	Bfactor = B1;
	N = N1;
	fractionSize = S1;
	rmax = r1;
	v = v1;
	step = 0.4;
	if(fabs(gradB) > 0){
		step = 0.4*fabs(Bfactor/(B0*gradB));
	}
	if(fabs(gradN) > 0){
		step = min(step, 0.4*fabs(N/(N0*gradN)));
	}
	if(gradS > 0){
		step = min(step, 0.4*fabs((maxFraction - fractionSize)/gradS));
	}
	if(gradS < 0){
		step = min(step, 0.4*fabs((fractionSize - minFraction)/gradS));
	}
	if(fabs(gradR) > 0){
		step = min(step, 0.4*fabs(rmax/(R0*gradR)));
	}
	if(gradV > 0){
		step = min(step, 0.4*fabs((maxV - v)/(V0*gradV)));
	}
	if(gradV < 0){
		step = min(step, 0.4*fabs(v/(V0*gradV)));
	}
	double B2 = B1 - gradB*step*B0;
	double N2 = N1 - gradN*step*N0;
	double S2 = S1 - gradS*step;
	double r2 = r1  - gradR*step*R0;
	double v2 = v1 - gradV*step*V0;
	if(B2 != B2){
		printf("B2 Nan\n");
		exit(0);
	}
	if(S2 > maxFraction){
		printf("S2 > 1\n");
		S2 = maxFraction;
	}
	if(S2 < minFraction){
		S2 = minFraction;
	}
	if(v2 > maxV){
		v2 = maxV;
	}
	double f2 = evaluateOptimizationFunction5simple(B2, N2, S2, r2, v2, Ee, Fe, Np, Nd, sintheta, thetaIndex);
	int iterations = 0;
	while(f2 < f1){
		iterations++;
		Bfactor = B2;
		N = N2;
		fractionSize = S2;
		rmax = r2;
		v = v2;
		step = 0.4;
		if(fabs(gradB) > 0){
			step = 0.4*fabs(Bfactor/(B0*gradB));
		}
		if(fabs(gradN) > 0){
			step = min(step, 0.4*fabs(N/(N0*gradN)));
		}
		if(gradS > 0){
			step = min(step, 0.4*fabs((maxFraction - fractionSize)/gradS));
		}
		if(gradS < 0){
			step = min(step, 0.4*fabs((fractionSize - minFraction)/gradS));
		}
		if(fabs(gradR) > 0){
			step = min(step, 0.4*fabs(rmax/(R0*gradR)));
		}
		if(gradV > 0){
			step = min(step, 0.4*fabs((maxV - v)/(V0*gradV)));
		}
		if(gradV < 0){
			step = min(step, 0.4*fabs(v/(V0*gradV)));
		}
		B2 = Bfactor - gradB*step*B0;
		N2 = N - gradN*step*N0;
		S2 = fractionSize - gradS*step;
		r2 = rmax - gradR*step*R0;
		v2 = v - gradV*step*V0;
		if(B2 != B2){
			printf("B2 Nan\n");
			exit(0);
		}
		if(S2 > maxFraction){
			printf("S2 > 1\n");
			S2 = maxFraction;
			return;
		}
		if(S2 < minFraction){
			S2 = minFraction;
			//printf("S2 < 0\n");
			return;
		}
		if(v2 > maxV){
			v2 = maxV;
		}
		f1 = f2;
		f2 = evaluateOptimizationFunction5simple(B2, N2, S2, r2, v2, Ee, Fe, Np, Nd, sintheta, thetaIndex);
	}

	//findMinParameters(Bfactor, N, minLambda, maxLambda, gradB, gradN, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
}

void optimizeParameters5simple(const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v, double** Ee, double** Fe, int Np, int Nd, double sintheta, int thetaIndex, FILE* logFile) {
	double currentF = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g n = %g fraction = %g rmax = %g v = %g\n", Bfactor, N, fractionSize, rmax, v);
	fprintf(logFile, "Bfactor = %g n = %g fraction = %g rmax = %g v = %g\n", Bfactor, N, fractionSize, rmax, v);
	for(int i = 0; i < Niterations; ++i) {
		///randomization;
		for(int j = 0; j < 5; ++j){
			double tempN = N + 0.2*N*(uniformDistribution() - 0.5);
			double tempB = Bfactor + 0.2*Bfactor*(uniformDistribution() - 0.5);
			double tempS  = min(max(minFraction,fractionSize + 0.2*fractionSize*(uniformDistribution() - 0.5)), maxFraction);
			double tempF = evaluateOptimizationFunction5simple(tempB, tempN, tempS, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
			if(tempF < currentF){
				currentF = tempF;
				Bfactor = tempB;
				fractionSize = tempS;
				N = tempN;
				printf("random search\n");
				fprintf(logFile, "random search\n");
			}
		}
		double prevF = currentF;
		//
		//valley first step
		double valleyB1 = Bfactor;
		double valleyN1 = N;
		double valleyS1 = fractionSize;
		double valleyR1 = rmax;
		double valleyV1 = v;

		double dxB = fabs(Bfactor/B0)/20;
		double dxN = fabs(N/N0)/20;
		double dxS = fabs(fractionSize)/20;
		double dxR = fabs(rmax/R0)/20;
		double dxV = fabs(v/V0)/20;
		printf("optimization i = %d\n",i);
		fprintf(logFile, "optimization i = %d\n",i);
		fflush(logFile);

		double gradB;
		double gradN;
		double gradS;
		double gradR;
		double gradV;

		double Fb = evaluateOptimizationFunction5simple(Bfactor + dxB*B0, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		double Fb1 = evaluateOptimizationFunction5simple(Bfactor - dxB*B0, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);

		while( Fb > currentF && Fb1 > currentF){
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction5simple(Bfactor + dxB*B0, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
			Fb1 = evaluateOptimizationFunction5simple(Bfactor - dxB*B0, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		double Fn = evaluateOptimizationFunction5simple(Bfactor, N + dxN*N0, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		double Fn1 = evaluateOptimizationFunction5simple(Bfactor, N - dxN*N0, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);

		while( Fn > currentF && Fn1 > currentF){
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction5simple(Bfactor, N + dxN*N0, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
			Fn1 = evaluateOptimizationFunction5simple(Bfactor, N - dxN*N0, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		double Fs = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize + dxS, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		double Fs1 = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize - dxS, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		gradS = (Fs - Fs1)/(2*dxS);
		if(fractionSize >= maxFraction){
			if(gradS > 0){
				gradS = 0;
			}
		}
		if(fractionSize <= minFraction){
			if(gradS < 0){
				gradS = 0;
			}
		}

		double Fr = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax + dxR*R0, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - currentF)/(dxR);

		double Fv = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v + dxV*V0, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		double Fv1 = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v - dxV*V0, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradV = (Fv - Fv1)/(2*dxV);
		if(v >= maxV){
			if(gradV > 0){
				gradV = 0;
			}
		}

		double gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS + gradR*gradR + gradV*gradV);
		gradB = gradB/gradNorm;
		gradN = gradN/gradNorm;
		gradS = gradS/gradNorm;
		gradR = gradR/gradNorm;
		gradV = gradV/gradNorm;

		findMinParameters5simple(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, Ee, Fe, Np, Nd, sintheta, thetaIndex, currentF);

		currentF = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		//valley second step

		dxB = fabs(Bfactor/B0)/20;
		dxN = fabs(N/N0)/20;
		dxS = fabs(fractionSize)/20;
		dxR = fabs(rmax/R0)/20;
		dxV = fabs(v/V0)/20;

		fflush(logFile);

		Fb = evaluateOptimizationFunction5simple(Bfactor + dxB*B0, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		Fb1 = evaluateOptimizationFunction5simple(Bfactor - dxB*B0, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);

		while( Fb > currentF && Fb1 > currentF){
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction5simple(Bfactor + dxB*B0, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
			Fb1 = evaluateOptimizationFunction5simple(Bfactor - dxB*B0, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		Fn = evaluateOptimizationFunction5simple(Bfactor, N + dxN*N0, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		Fn1 = evaluateOptimizationFunction5simple(Bfactor, N - dxN*N0, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);

		while( Fn > currentF && Fn1 > currentF){
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction5simple(Bfactor, N + dxN*N0, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
			Fn1 = evaluateOptimizationFunction5simple(Bfactor, N - dxN*N0, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		Fs = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize + dxS, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		Fs1 = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize - dxS, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		gradS = (Fs - Fs1)/(2*dxS);
		if(fractionSize >= maxFraction){
			if(gradS > 0){
				gradS = 0;
			}
		}
		if(fractionSize <= minFraction){
			if(gradS < 0){
				gradS = 0;
			}
		}

		Fr = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax + dxR*R0, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - currentF)/(dxR);

		Fv = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v + dxV*V0, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		Fv1 = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v - dxV*V0, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradV = (Fv - Fv1)/(2*dxV);
		if(v >= maxV){
			if(gradV > 0){
				gradV = 0;
			}
		}

		gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS + gradR*gradR + gradV*gradV);
		gradB = gradB/gradNorm;
		gradN = gradN/gradNorm;
		gradS = gradS/gradNorm;
		gradR = gradR/gradNorm;
		gradV = gradV/gradNorm;

		findMinParameters5simple(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, Ee, Fe, Np, Nd, sintheta, thetaIndex, currentF);

		currentF = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);

		double valleyB2 = Bfactor;
		double valleyN2 = N;
		double valleyS2 = fractionSize;
		double valleyR2 = rmax;
		double valleyV2 = v;

		//// valley third step

		
		gradB = (valleyB1 - valleyB2)/B0;
		gradN = (valleyN1 - valleyN2)/N0;
		gradS = (valleyS1 - valleyS2);
		gradR = (valleyR1 - valleyR2)/R0;
		gradV = (valleyV1 - valleyV2)/V0;
		gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS + gradR*gradR + gradV*gradV);
		if(gradNorm > 0){
			gradB = gradB/gradNorm;
			gradN = gradN/gradNorm;
			gradS = gradS/gradNorm;
			gradR = gradR/gradNorm;
			gradV = gradV/gradNorm;
			//double step = 0.1*min3(Bfactor, N, fractionSize);
			//double Fv = evaluateOptimizationFunction(Bfactor - gradB*step, N - gradN*step, fractionSize - gradS*step, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		//if(Fv < currentF){
				findMinParameters5simple(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, Ee, Fe, Np, Nd, sintheta, thetaIndex, currentF);
		//} else {
			//findMinParameters3(Bfactor, N, fractionSize, -gradB, -gradN, gradS, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, Z, currentF, totalInu);
		//}
		}
		currentF = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v, Ee, Fe, Np, Nd, sintheta, thetaIndex);
		
		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/

		printf("optimization function = %g\n", currentF);
		fprintf(logFile, "optimization function = %g\n", currentF);
		printf("Bfactor = %g n = %g fraction = %g rmax = %g v/c = %g\n", Bfactor, N, fractionSize, rmax, v/speed_of_light);
		fprintf(logFile, "Bfactor = %g n = %g fraction = %g rmax = %g v/c =%g\n", Bfactor, N, fractionSize, rmax, v/speed_of_light);
	}
	printf("finish optimization\n");
	fprintf(logFile, "finish optimization\n");
}