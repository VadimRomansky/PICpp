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

//for fla disc

double evaluateOptimizationFunction5(double Bfactor, double n, double fractionSize, double rmax, double v, double** nu, double** observedInu, double** Ee, double** dFe, int Np, int Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double***** Inu, double***** Anu, double*** area, double*** length){
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
		evaluateAllEmissivityAndAbsorption(nu[i], Inu[i], Anu[i], Nnu, Ee, dFe, Np, Nd, Bn, sintheta, thetaIndex, concentrations, n, Bfactor, rfactor);
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
				err1 = sqr(totalInu[j] - observedInu[i][j]);
			} else {
				err1 = sqr(log(totalInu[j]) - log(observedInu[i][j]));
			}
			err = err + err1;
		}
	}
	delete[] totalInu;
	
	return err;
}

void findMinParameters5(const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v, double gradB, double gradN, double gradS, double gradR, double gradV, double** nu, double** observedInu, double** Ee, double** dFe, int Np, int Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double***** Inu, double***** Anu, double*** area, double*** length, double& currentF){
	double step = 0.4;
	/*if(fabs(gradB) > 0){
		step = 0.4*fabs(Bfactor/(B0*gradB));
	}*/
	if(gradB > 0){
		step = 0.4*fabs((Bfactor-minB)/(B0*gradB));
	}
	if(gradB < 0){
		step = 0.4*fabs((maxB - Bfactor)/(B0*gradB));
	}
	/*if(fabs(gradN) > 0){
		step = min(step, 0.4*fabs(N/(N0*gradN)));
	}*/
	if(gradN > 0){
		step = min(step, 0.4*fabs((N-minN)/(N0*gradN)));
	}
	if(gradN < 0){
		step = min(step, 0.4*fabs((maxN - N)/(N0*gradN)));
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
	if(N1 < 0) {
		printf("N1 < 0\n");
	}
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

	double f1 = evaluateOptimizationFunction5(B1, N1, S1, r1, v1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	if(f1 > currentF){
		int count = 0;
		while (f1 > currentF && count < 10){
			count++;
			step = step/2;
			B1 = Bfactor - gradB*step*B0;
			N1 = N - gradN*step*N0;
			if(N1 < 0) {
				printf("N1 < 0\n");
			}
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

			f1 = evaluateOptimizationFunction5(B1, N1, S1, r1, v1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
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
	if(gradB > 0){
		step = 0.4*fabs((Bfactor-minB)/(B0*gradB));
	}
	if(gradB < 0){
		step = 0.4*fabs((maxB - Bfactor)/(B0*gradB));
	}
	if(gradN > 0){
		step = min(step, 0.4*fabs((N-minN)/(N0*gradN)));
	}
	if(gradN < 0){
		step = min(step, 0.4*fabs((maxN - N)/(N0*gradN)));
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
	if(N2 < 0) {
		printf("N2 < 0\n");
	}
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
	double f2 = evaluateOptimizationFunction5(B2, N2, S2, r2, v2, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	int iterations = 0;
	while((f2 < f1) && (iterations < 100)){
		iterations++;
		Bfactor = B2;
		N = N2;
		fractionSize = S2;
		rmax = r2;
		v = v2;
		step = 0.4;
		if(gradB > 0){
			step = 0.4*fabs((Bfactor-minB)/(B0*gradB));
		}
		if(gradB < 0){
			step = 0.4*fabs((maxB - Bfactor)/(B0*gradB));
		}
		if(gradN > 0){
			step = min(step, 0.4*fabs((N-minN)/(N0*gradN)));
		}
		if(gradN < 0){
			step = min(step, 0.4*fabs((maxN - N)/(N0*gradN)));
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
		if(N2 < 0) {
			printf("N2 < 0\n");
		}
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
		f2 = evaluateOptimizationFunction5(B2, N2, S2, r2, v2, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	}
}

void optimizeParameters5(const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v,  double** nu, double** observedInu, double** Ee, double** dFe, int Np, int Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double***** Inu, double***** Anu, double*** area, double*** length, FILE* logFile){
	double currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	printf("optimization function = %g\n", currentF);
	fprintf(logFile, "optimization function = %g\n", currentF);
	printf("Bfactor = %g n = %g fraction = %g rmax = %g v = %g\n", Bfactor, N, fractionSize, rmax, v);
	fprintf(logFile, "Bfactor = %g n = %g fraction = %g rmax = %g v = %g\n", Bfactor, N, fractionSize, rmax, v);
	for(int i = 0; i < Niterations; ++i) {
		///randomization;
		for(int j = 0; j < 5; ++j){
			double tempN = N + 0.2*N*(uniformDistribution() - 0.5);
			if(tempN < minN){
				tempN = minN;
			}
			if(tempN > maxN){
				tempN = maxN;
			}
			//tempN = N;
			double tempB = Bfactor + 0.2*Bfactor*(uniformDistribution() - 0.5);
			if(tempB < minB){
				tempB = minB;
			}
			if(tempB > maxB){
				tempB = maxB;
			}
			//tempB = Bfactor;
			double tempS  = min(max(minFraction,fractionSize + 0.2*fractionSize*(uniformDistribution() - 0.5)), maxFraction);
			if(tempS < minFraction){
				tempS = minFraction;
			}
			if(tempS > maxFraction){
				tempS = maxFraction;
			}
			double tempR = rmax + 0.2*rmax*(uniformDistribution() - 0.5);
			double tempV = 0.2*speed_of_light + (maxV - 0.2*speed_of_light)*uniformDistribution();
			if(tempV > maxV){
				tempV = maxV;
			}
			double tempF = evaluateOptimizationFunction5(tempB, tempN, tempS, tempR, tempV, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
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

		double Fb = evaluateOptimizationFunction5(Bfactor + dxB*B0, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fb1 = evaluateOptimizationFunction5(Bfactor - dxB*B0, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		int count = 0;
		while( Fb > currentF && Fb1 > currentF && count < 5){
			count ++;
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction5(Bfactor + dxB*B0, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fb1 = evaluateOptimizationFunction5(Bfactor - dxB*B0, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		double Fn = evaluateOptimizationFunction5(Bfactor, N + dxN*N0, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fn1 = evaluateOptimizationFunction5(Bfactor, N - dxN*N0, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fn > currentF && Fn1 > currentF && count < 5){
			count++;
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction5(Bfactor, N + dxN*N0, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fn1 = evaluateOptimizationFunction5(Bfactor, N - dxN*N0, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		double Fs = evaluateOptimizationFunction5(Bfactor, N, fractionSize + dxS, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fs1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize - dxS, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
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

		double Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while( Fr > currentF && Fr1 > currentF && count < 5){
			count++;
			dxR = dxR/2;
			Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - Fr1)/(2*dxR);

		double Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while(Fv > currentF && Fv1 > currentF && count < 5){
			count++;
			dxV = dxV/2;
			Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradV = (Fv - Fv1)/(2*dxV);
		if(v >= maxV){
			if(gradV > 0){
				gradV = 0;
			}
		}

		//gradB = 0;
		//gradN = 0;

		double gradNorm = sqrt(gradB*gradB + gradN*gradN + gradS*gradS + gradR*gradR + gradV*gradV);
		if(gradNorm <= 0){
			printf("gradNorm <= 0\n");
		}
		if(gradNorm != gradNorm){
			printf("gradNorm = NaN\n");
			exit(0);
		}
		gradB = gradB/gradNorm;
		gradN = gradN/gradNorm;
		gradS = gradS/gradNorm;
		gradR = gradR/gradNorm;
		gradV = gradV/gradNorm;

		findMinParameters5(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);

		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		//valley second step

		dxB = fabs(Bfactor/B0)/20;
		dxN = fabs(N/N0)/20;
		dxS = fabs(fractionSize)/20;
		dxR = fabs(rmax/R0)/20;
		dxV = fabs(v/V0)/20;

		fflush(logFile);

		Fb = evaluateOptimizationFunction5(Bfactor + dxB*B0, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fb1 = evaluateOptimizationFunction5(Bfactor - dxB*B0, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fb > currentF && Fb1 > currentF && count < 5){
			count++;
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction5(Bfactor + dxB*B0, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fb1 = evaluateOptimizationFunction5(Bfactor - dxB*B0, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		Fn = evaluateOptimizationFunction5(Bfactor, N + dxN*N0, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fn1 = evaluateOptimizationFunction5(Bfactor, N - dxN*N0, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fn > currentF && Fn1 > currentF && count < 5){
			count++;
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction5(Bfactor, N + dxN*N0, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fn1 = evaluateOptimizationFunction5(Bfactor, N - dxN*N0, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		Fs = evaluateOptimizationFunction5(Bfactor, N, fractionSize + dxS, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fs1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize - dxS, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
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

		Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fr > currentF && Fr1 > currentF && count < 5){
			count++;
			dxR = dxR/2;
			Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - Fr1)/(2*dxR);

		Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while(Fv > currentF && Fv1 > currentF && count < 5){
			count++;
			dxV = dxV/2;
			Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
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

		findMinParameters5(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);

		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

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
			//double Fv = evaluateOptimizationFunction(Bfactor - gradB*step, N - gradN*step, fractionSize - gradS*step, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		//if(Fv < currentF){
				findMinParameters5(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);
		//} else {
			//findMinParameters3(Bfactor, N, fractionSize, -gradB, -gradN, gradS, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, Z, currentF, totalInu);
		//}
		}
		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		
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

void findMinParameters5sigma(const double& sigma, const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v, double gradB, double gradS, double gradR, double gradV, double** nu, double** observedInu, double** Ee, double** dFe, int Np, int Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double***** Inu, double***** Anu, double*** area, double*** length, double& currentF){
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

	double f1 = evaluateOptimizationFunction5(B1, N1, S1, r1, v1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
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

			f1 = evaluateOptimizationFunction5(B1, N1, S1, r1, v1, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
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
	double f2 = evaluateOptimizationFunction5(B2, N2, S2, r2, v2, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
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
		f2 = evaluateOptimizationFunction5(B2, N2, S2, r2, v2, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
	}
}

void optimizeParameters5sigma(const double& sigma, const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v,  double** nu, double** observedInu, double** Ee, double** dFe, int Np, int Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double***** Inu, double***** Anu, double*** area, double*** length, FILE* logFile){
	double currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
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
			double tempF = evaluateOptimizationFunction5(tempB, tempN, tempS, tempR, tempV, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
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

		double Fb = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		tempB = Bfactor - dxB*B0;
		tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
		double Fb1 = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		int count = 0;
		while( Fb > currentF && Fb1 > currentF && count < 5){
			count ++;
			dxB = dxB/2;
			tempB = Bfactor + dxB*B0;
			tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
			Fb = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			tempB = Bfactor - dxB*B0;
			tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
			Fb1 = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		double Fs = evaluateOptimizationFunction5(Bfactor, N, fractionSize + dxS, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fs1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize - dxS, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
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

		double Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while( Fr > currentF && Fr1 > currentF && count < 5){
			count++;
			dxR = dxR/2;
			Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - Fr1)/(2*dxR);

		double Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		double Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while(Fv > currentF && Fv1 > currentF && count < 5){
			count++;
			dxV = dxV/2;
			Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
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

		findMinParameters5sigma(sigma, B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradS, gradR, gradV, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);

		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		//valley second step

		dxB = fabs(Bfactor/B0)/20;
		dxS = fabs(fractionSize)/20;
		dxR = fabs(rmax/R0)/20;
		dxV = fabs(v/V0)/20;

		fflush(logFile);

		tempB = Bfactor + dxB*B0;
		tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
		Fb = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		tempB = Bfactor - dxB*B0;
		tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
		Fb1 = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fb > currentF && Fb1 > currentF && count < 5){
			count++;
			dxB = dxB/2;
			tempB = Bfactor + dxB*B0;
			tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
			Fb = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			tempB = Bfactor - dxB*B0;
			tempN = sqr(tempB)/(sigma*4*pi*massProtonReal*speed_of_light2);
			Fb1 = evaluateOptimizationFunction5(tempB, tempN, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		gradB = (Fb - Fb1)/(2*dxB);


		Fs = evaluateOptimizationFunction5(Bfactor, N, fractionSize + dxS, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fs1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize - dxS, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
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

		Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

		count = 0;
		while( Fr > currentF && Fr1 > currentF && count < 5){
			count++;
			dxR = dxR/2;
			Fr = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax + dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fr1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax - dxR*R0, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - Fr1)/(2*dxR);

		Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		count = 0;
		while(Fv > currentF && Fv1 > currentF && count < 5){
			count++;
			dxV = dxV/2;
			Fv = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v + dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
			Fv1 = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v - dxV*V0, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		}
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
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

		findMinParameters5sigma(sigma, B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradS, gradR, gradV, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);

		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);

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
			//double Fv = evaluateOptimizationFunction(Bfactor - gradB*step, N - gradN*step, fractionSize - gradS*step, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		//if(Fv < currentF){
				findMinParameters5sigma(sigma, B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradS, gradR, gradV, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, currentF);
		//} else {
			//findMinParameters3(Bfactor, N, fractionSize, -gradB, -gradN, gradS, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, Z, currentF, totalInu);
		//}
		}
		currentF = evaluateOptimizationFunction5(Bfactor, N, fractionSize, rmax, v, nu, observedInu, Ee, dFe, Np, Nnu, Nd, Nmonth, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length);
		
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

double evaluateOptimizationFunction5simple(double Bfactor, double n, double fractionSize, double rmax, double v, double** Ee, double** dFe, int Np, int Nd, double sintheta, int thetaIndex) {
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
		evaluateEmissivityAndAbsorptionAtNuSimple(nu, Inu, Anu, Ee[thetaIndex], dFe[thetaIndex], Np, sintheta, B, N);
		evaluateSpectrumAtNuSimple(nu, I, Inu, Anu, r, fractionSize);
		err = err + sqr(I - apry[i]);
	}

	r = rmax + v*times[1];
	B = Bfactor*rmax/r;
	N = n*sqr(rmax/r);
	for(int i = 0; i < 3; ++i){
		double I = 0;
		nu = mayx[i]*1E9;
		evaluateEmissivityAndAbsorptionAtNuSimple(nu, Inu, Anu, Ee[thetaIndex], dFe[thetaIndex], Np, sintheta, B, N);
		evaluateSpectrumAtNuSimple(nu, I, Inu, Anu, r, fractionSize);
		err = err + sqr(I - mayy[i]);
	}

	r = rmax + v*times[2];
	B = Bfactor*rmax/r;
	N = n*sqr(rmax/r);
	for(int i = 0; i < 4; ++i){
		double I = 0;
		nu = junx[i]*1E9;
		evaluateEmissivityAndAbsorptionAtNuSimple(nu, Inu, Anu, Ee[thetaIndex], dFe[thetaIndex], Np, sintheta, B, N);
		evaluateSpectrumAtNuSimple(nu, I, Inu, Anu, r, fractionSize);
		err = err + sqr(I - juny[i]);
	}

	r = rmax + v*times[3];
	B = Bfactor*rmax/r;
	N = n*sqr(rmax/r);
	for(int i = 0; i < 5; ++i){
		double I = 0;
		nu = augx[i]*1E9;
		evaluateEmissivityAndAbsorptionAtNuSimple(nu, Inu, Anu, Ee[thetaIndex], dFe[thetaIndex], Np, sintheta, B, N);
		evaluateSpectrumAtNuSimple(nu, I, Inu, Anu, r, fractionSize);
		err = err + sqr(I - augy[i]);
	}

	return err;
}



void findMinParameters5simple(const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v, double gradB, double gradN, double gradS, double gradR, double gradV, double** Ee, double** dFe, int Np, int Nd, double sintheta, int thetaIndex, double& currentF) {
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

	double f1 = evaluateOptimizationFunction5simple(B1, N1, S1, r1, v1, Ee, dFe, Np, Nd, sintheta, thetaIndex);
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

			f1 = evaluateOptimizationFunction5simple(B1, N1, S1, r1, v1, Ee, dFe, Np, Nd, sintheta, thetaIndex);
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
	double f2 = evaluateOptimizationFunction5simple(B2, N2, S2, r2, v2, Ee, dFe, Np, Nd, sintheta, thetaIndex);
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
		f2 = evaluateOptimizationFunction5simple(B2, N2, S2, r2, v2, Ee, dFe, Np, Nd, sintheta, thetaIndex);
	}

	//findMinParameters(Bfactor, N, minLambda, maxLambda, gradB, gradN, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
}

void optimizeParameters5simple(const double& B0, const double& N0, const double& R0, const double& V0, double& Bfactor, double& N, double& fractionSize, double& rmax, double& v, double** Ee, double** dFe, int Np, int Nd, double sintheta, int thetaIndex, FILE* logFile) {
	double currentF = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
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
			double tempF = evaluateOptimizationFunction5simple(tempB, tempN, tempS, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
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

		double Fb = evaluateOptimizationFunction5simple(Bfactor + dxB*B0, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		double Fb1 = evaluateOptimizationFunction5simple(Bfactor - dxB*B0, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);

		while( Fb > currentF && Fb1 > currentF){
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction5simple(Bfactor + dxB*B0, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
			Fb1 = evaluateOptimizationFunction5simple(Bfactor - dxB*B0, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		double Fn = evaluateOptimizationFunction5simple(Bfactor, N + dxN*N0, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		double Fn1 = evaluateOptimizationFunction5simple(Bfactor, N - dxN*N0, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);

		while( Fn > currentF && Fn1 > currentF){
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction5simple(Bfactor, N + dxN*N0, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
			Fn1 = evaluateOptimizationFunction5simple(Bfactor, N - dxN*N0, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		double Fs = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize + dxS, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		double Fs1 = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize - dxS, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
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

		double Fr = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax + dxR*R0, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - currentF)/(dxR);

		double Fv = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v + dxV*V0, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		double Fv1 = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v - dxV*V0, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
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

		findMinParameters5simple(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, Ee, dFe, Np, Nd, sintheta, thetaIndex, currentF);

		currentF = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		//valley second step

		dxB = fabs(Bfactor/B0)/20;
		dxN = fabs(N/N0)/20;
		dxS = fabs(fractionSize)/20;
		dxR = fabs(rmax/R0)/20;
		dxV = fabs(v/V0)/20;

		fflush(logFile);

		Fb = evaluateOptimizationFunction5simple(Bfactor + dxB*B0, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		Fb1 = evaluateOptimizationFunction5simple(Bfactor - dxB*B0, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);

		while( Fb > currentF && Fb1 > currentF){
			dxB = dxB/2;
			Fb = evaluateOptimizationFunction5simple(Bfactor + dxB*B0, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
			Fb1 = evaluateOptimizationFunction5simple(Bfactor - dxB*B0, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		}
		gradB = (Fb - Fb1)/(2*dxB);

		Fn = evaluateOptimizationFunction5simple(Bfactor, N + dxN*N0, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		Fn1 = evaluateOptimizationFunction5simple(Bfactor, N - dxN*N0, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);

		while( Fn > currentF && Fn1 > currentF){
			dxN = dxN/2;
			Fn = evaluateOptimizationFunction5simple(Bfactor, N + dxN*N0, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
			Fn1 = evaluateOptimizationFunction5simple(Bfactor, N - dxN*N0, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		}
		gradN = (Fn - Fn1)/(2*dxN);

		Fs = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize + dxS, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		Fs1 = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize - dxS, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
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

		Fr = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax + dxR*R0, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		//double Feta1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
		gradR = (Fr - currentF)/(dxR);

		Fv = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v + dxV*V0, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		Fv1 = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v - dxV*V0, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		//double Fv1 = evaluateOptimizationFunction(Bfactor, N, fractionSize, eta - dxEta, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);
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

		findMinParameters5simple(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, Ee, dFe, Np, Nd, sintheta, thetaIndex, currentF);

		currentF = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);

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
			//double Fv = evaluateOptimizationFunction(Bfactor - gradB*step, N - gradN*step, fractionSize - gradS*step, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, totalInu);

		//if(Fv < currentF){
				findMinParameters5simple(B0, N0, R0, V0, Bfactor, N, fractionSize, rmax, v, gradB, gradN, gradS, gradR, gradV, Ee, dFe, Np, Nd, sintheta, thetaIndex, currentF);
		//} else {
			//findMinParameters3(Bfactor, N, fractionSize, -gradB, -gradN, gradS, nu, Ee, dFe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, Inu, Anu, area, length, Rho, Phi, Z, currentF, totalInu);
		//}
		}
		currentF = evaluateOptimizationFunction5simple(Bfactor, N, fractionSize, rmax, v, Ee, dFe, Np, Nd, sintheta, thetaIndex);
		
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