// SphericalTurbulence.cpp : Defines the entry point for the console application.
//
#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include <time.h>
#include <cmath>
#include <string>
#include <omp.h>

const std::string outputfileName = "radiation.dat";
//const std::string fileNameP = "../../tristan-mp-pitp/Pe";
//const std::string fileNameF = "../../tristan-mp-pitp/Fe";
const std::string fileNameP = "Pe";
const std::string fileNameF = "Fe";
const std::string logFileName = "log.dat";
const std::string BFileName = "B.dat";

const int randomParameter = 1024;
const double atomicUnionMass = 1.66053904E-24;
const double massProtonReal = 1.67262177E-24;
const double massAlphaReal = 6.644656E-24;
const double massElectron = 0.910938291E-27;
//const double massElectronFactor = massProtonReal/(10.0*massElectronReal);
//const double massElectronFactor = 100;
const double massDeuteriumReal = 3.34449696893E-24;
const double massHelium3Real = 5.00823792874E-24;
const double massOxygenReal = 26.5601801672E-24;
const double massSiliconReal = 46.4567787264E-24;
const double kBoltzman = 1.3806488E-16;
const double speed_of_light = 2.99792458E10;
const double speed_of_light2 = speed_of_light * speed_of_light;
const double speed_of_light4 = speed_of_light2 * speed_of_light2;
const double electron_charge = 4.803529695E-10;
const double pi = 4*atan2(1.0,1.0);
const double four_pi = 4*pi;

const double emissivityCoef = sqrt(3.0) * electron_charge * electron_charge * electron_charge / (massElectron * speed_of_light2);
const double absorpCoef = 16 * pi * pi * electron_charge / (3 * sqrt(3.0));
const double criticalNuCoef = 3 * electron_charge / (4 * pi * massElectron * massElectron * massElectron * speed_of_light * speed_of_light4);

const int Nx = 20;
const int Ny = 20;
const int Nz = 20;
const int Nk = 10;

const int Napprox = 52;

///      K5/3(x)
const double McDonaldValue[Napprox] = { 1.43E15, 9.8E13, 3.09E13, 2.11E12, 6.65E11, 4.55E10, 1.43E10, 9.8E8,
	3.08E8, 2.11E7, 6.65E6, 4.55E5, 1.43E4, 9802, 3087, 670, 211, 107, 66.3, 33.6, 20.7, 14.1, 11.0, 10.3, 6.26, 4.2, 3.01, 2.25, 1.73, 1.37, 1.10,
	0.737, 0.514, 0.368, 0.269, 0.2, 0.0994, 0.0518, 0.0278, 0.0152, 0.00846, 0.00475, 0.00154, 0.000511, 0.000172, 0.0000589, 0.0000203, 0.00000246, 3.04E-7, 3.81E-8, 4.82E-9, 6.14E-10
};
////     x*int from x to inf K5/3(t)dt
const double UvarovValue[Napprox] = {0.0021495, 0.00367, 0.00463, 0.00792, 0.00997, 0.017, 0.0215, 0.0367,
	0.0461, 0.0791, 0.0995, 0.169, 0.213, 0.358, 0.445, 0.583, 0.702, 0.772, 0.818, 0.874,0.904, 0.917, 0.918, 0.918, 0.901, 0.872, 0.832, 0.788,
	0.742, 0.694, 0.655, 0.566, 0.486, 0.414, 0.354, 0.301, 0.200, 0.130, 0.0845, 0.0541, 0.0339, 0.0214, 0.0085, 0.0033, 0.0013, 0.00050, 0.00019,
	0.0000282, 0.00000409, 5.89E-7, 8.42E-8, 1.19E-8
};
///// x
const double UvarovX[Napprox] = {1.0E-9, 5.0E-9, 1.0E-8, 5.0E-8, 1.0E-7, 5.0E-7, 1.0E-6, 5.0E-6,
	0.00001, 0.00005, 0.0001 ,0.0005, 0.001, 0.005, 0.01, 0.025, 0.050, 0.075, 0.10, 0.15,
	0.20, 0.25, 0.29, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0, 1.2, 1.4, 1.6,
	1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0
};

const double augx[5] = {0.335, 0.625, 1.46, 4.92, 8.57};
const double augy[5] = {3.29, 7.77, 8.53, 2.42, 1.06};

const double augmaxx = 0.886;
const double augmaxy = 11.2;

const double junx[4] = {0.628, 1.45, 4.89, 8.53};
const double juny[4] = {2.98, 12.3, 5.79, 3.15};

const double junmaxx = 1.65;
const double junmaxy = 13.2;

const double mayx[3]= {1.46, 4.94, 8.62};
const double mayy[3] = {4.91, 12.0, 6.67};

const double maymaxx = 2.96;
const double maymaxy = 15.2;

const double aprx[4] = {1.44, 4.91, 8.58, 22.8};
const double apry[4] ={0.993, 13.9, 17.1, 5.11};

const double aprmaxx = 6.50;
const double aprmaxy = 19.3;

const double minB = 0.05;
const double maxB = 3.0;
const double minN = 0.01;
const double maxN = 10;

double uniformDistribution() {
	return (rand() % randomParameter + 0.5) / randomParameter;
}

double sqr(const double& a){
	return a*a;
}

double min(const double& a, const double& b){
	if(a < b) {
		return a;
	} else {
		return b;
	}
}

std::string convertIntToString(int a) {
	if (a == 0) {
		std::string result = "0";
		return result;
	}
	if (a > 0) {
		std::string result = "";
		while (a > 0) {
			int last = a % 10;
			a = a / 10;
			char c = last + '0';
			result = c + result;
		}
		return result;
	}
	a = -a;
	std::string result = "-";
	return result + convertIntToString(a);
}

void findMaxNu(int& nuMaxIndex, double* Inu, int Nnu) {
	double Imax = 0;
	nuMaxIndex = 0;
	for(int i = 0; i < Nnu; ++i) {
		if(Inu[i] > Imax) {
			Imax = Inu[i];
			nuMaxIndex = i;
		}
	}
}

double evaluateTurbB(const double& kx, const double& ky, const double& kz, const double turbNorm){
	double k = sqrt(kx*kx + ky*ky + kz*kz);
	return turbNorm/pow(k, 11.0/3.0);
}

double evaluateTurbNorm(const double& kmax, const int& Nk, const double& B0, const double& energyFraction){
	double dk = kmax/Nk;
	double sum = 0;
	for(int i = 0; i < Nk; ++i){
		for(int j = 0; j < Nk; ++j){
			for(int k = 0; k < Nk; ++k){
				if((i + j + k) > 0){
					double kx = i*dk;
					double ky = j*dk;
					double kz = k*dk;
					double B = evaluateTurbB(kx, ky, kz, 1.0);

					sum = sum + 2*B*B;
				}
			}
		}
	}
	return B0*sqrt(energyFraction/((1.0 - energyFraction)*sum));
}

double evaluateMcDonaldIntegral(const double& nu) {
	int curIndex = 0;
	if (nu < UvarovX[0]) {
		//printf("x < UvarovX[0]\n");
		return 0;
	}
	if (nu > UvarovX[Napprox - 1]) {
		//printf("x > UvarovX[Napprox - 1]\n");
		return 0;
	}
	int leftIndex = 0;
	int rightIndex = Napprox-1;
	while(rightIndex - leftIndex > 1){
		int currentIndex = (rightIndex + leftIndex)/2;
		if(UvarovX[currentIndex] > nu){
			rightIndex = currentIndex;
		} else { 
			leftIndex = currentIndex;
		}
	}

	double result = (UvarovValue[rightIndex]*(nu - UvarovX[leftIndex]) + UvarovValue[leftIndex]*(UvarovX[rightIndex] - nu))/(UvarovX[rightIndex] - UvarovX[leftIndex]);
	//double result = UvarovValue[curIndex - 1] * exp(
	//	log(UvarovValue[curIndex] / UvarovValue[curIndex - 1]) * ((nu - UvarovX[curIndex - 1]) / (UvarovX[curIndex] - UvarovX[curIndex - 1])));
	if (result < 0) {
		printf("result < 0\n");
	}
	return result;
}

double evaluateMcDonaldFunction(const double& nu) {
	int curIndex = 0;
	if (nu < UvarovX[0]) {
		return 0;
	}
	if (nu > UvarovX[Napprox - 1]) {
		return 0;
	}
	int leftIndex = 0;
	int rightIndex = Napprox-1;
	while(rightIndex - leftIndex > 1){
		int currentIndex = (rightIndex + leftIndex)/2;
		if(UvarovX[currentIndex] > nu){
			rightIndex = currentIndex;
		} else { 
			leftIndex = currentIndex;
		}
	}

	double result = (McDonaldValue[rightIndex]*(nu - UvarovX[leftIndex]) + McDonaldValue[leftIndex]*(UvarovX[rightIndex] - nu))/(UvarovX[rightIndex] - UvarovX[leftIndex]);
	//double result = McDonaldValue[curIndex - 1] * exp(
	//	log(McDonaldValue[curIndex] / McDonaldValue[curIndex - 1]) * ((nu - UvarovX[curIndex - 1]) / (UvarovX[curIndex] - UvarovX[curIndex - 1])));
	if (result < 0) {
		printf("result < 0\n");
	}
	return result;
}

double criticalNu(const double& E, const double& sinhi, const double& H) {
	return criticalNuCoef * H * sinhi * E * E;
}

void evaluateLocalEmissivityAndAbsorption(double* nu, double* Inu, double* Anu, int Nnu, double* Ee, double* Fe, int Np, double sinhi, double B, double concentration) {
	Inu[0] = 0;
	Anu[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		Inu[i] = 0;
		Anu[i] = 0;
	}

	if(sinhi == 0.0){
		return;
	}

	double coef = concentration * emissivityCoef;
	double coefAbsorb = concentration * absorpCoef/(B*sinhi);


	for (int i = 0; i < Nnu; ++i) {
		//printf("i = %d\n", i);
		for (int j = 1; j < Np; ++j) {
			//if(Ee[j] < 100*massElectron*speed_of_light2){
				double nuc = criticalNu(Ee[j], sinhi, B);
				double gamma = Ee[j] / (massElectron * speed_of_light2);
				double x = nu[i] / nuc;
				//todo!!! 4pi!!
				//here Fe is (Fe[j] / (4*pi)) * (Ee[j] - Ee[j - 1])
				Inu[i] = Inu[i] + coef * Fe[j] * B * sinhi * evaluateMcDonaldIntegral(nu[i] / nuc);
				Anu[i] = Anu[i] + coefAbsorb * Fe[j] * evaluateMcDonaldFunction(nu[i] / nuc) / (gamma * gamma * gamma * gamma * gamma);
				/*if(Inu[i] != Inu[i]){
					printf("Inu NaN\n");
				}
				if(Anu[i] != Anu[i]){
					printf("Anu Nan\n");
				}*/
			//}
		}
	}

	//for (int i = 0; i < Nnu; ++i) {
	//	Inu[i] = Inu[i] * (1 - exp(-Anu[i] * fractionSize * localSize)) / (Anu[i] * fractionSize * localSize);
	//}
}

void evaluateAllEmissivityAndAbsorption(double* nu, double**** Inu, double**** Anu, int Nnu, double** Ee, double** Fe, int Np, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double concentration, double Bfactor, double* X, double* Y, double* Z, double rmax){
	int i = 0;
	#pragma omp parallel for shared(nu, Inu, Anu, Nnu, Ee, Fe, Np, Nd, Bn, sintheta, thetaIndex, concentration, Bfactor, X, Y, Z, rmax) private(i)
	for(i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			for(int k = 0; k < Nz; ++k){
				double r = sqrt(X[i]*X[i] + Y[j]*Y[j] + Z[k]*Z[k]);
				if(r < rmax){
					evaluateLocalEmissivityAndAbsorption(nu, Inu[i][j][k], Anu[i][j][k], Nnu, Ee[thetaIndex[i][j][k]], Fe[thetaIndex[i][j][k]], Np, sintheta[i][j][k], Bn[i][j][k]*Bfactor, concentration*concentrations[i][j][k]);
				} else {
					for(int l = 0; l < Nnu; ++l){
						Inu[i][j][k][l] = 0;
						Anu[i][j][k][l] = 0;
					}
				}
			}
		}
	}
}

void evaluateOrientationParameters(double*** B, double*** sintheta, int*** thetaIndex, double*** Bx, double*** By, double*** Bz, double* X, double* Y, double* Z, int Nd){
	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			for(int k = 0; k < Nz; ++k){
				double r = sqrt(X[i]*X[i] + Y[j]*Y[j] + Z[k]*Z[k]);
				double Bxy = sqrt(Bx[i][j][k]*Bx[i][j][k] + By[i][j][k]*By[i][j][k]);
				B[i][j][k] = sqrt(Bz[i][j][k]*Bz[i][j][k] + Bx[i][j][k]*Bx[i][j][k] + By[i][j][k]*By[i][j][k]);
				double cosTheta = Bz[i][j][k]/B[i][j][k];
				double sinTheta = Bxy/B[i][j][k];
				sintheta[i][j][k] = sinTheta;

				double cosTheta1 = Z[k]/r;
				double sinTheta1 = sqrt(X[i]*X[i] + Y[j]*Y[j])/r;
				double cosPhi1 = X[i]/sqrt(X[i]*X[i] + Y[j]*Y[j]);
				double sinPhi1 = Y[i]/sqrt(X[i]*X[i] + Y[j]*Y[j]);

				double Br = Bx[i][j][k]*sinTheta1*cosPhi1 + 
						    By[i][j][k]*sinTheta1*sinPhi1 +
							Bz[i][j][k]*cosTheta1;

				double cosTheta2 = Br/B[i][j][k];
				//double sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);


				double theta2 = acos(cosTheta2);
				if(theta2 >= pi/2){
					theta2 = pi - theta2;
				}

				thetaIndex[i][j][k] = floor(theta2 / (pi/(2*Nd))); 
				if(thetaIndex[i][j][k] == Nd){
					printf("thetaIndex == Nd\n");
					printf("%lf\n", By[i][j][k]);
					thetaIndex[i][j][k] = Nd - 1;
				}
				//for debug
				thetaIndex[i][j][k] = 3;
			}
		}
	}
}

void evaluateNu(double* nu, int Nnu, double minEnergy, double maxEnergy, double Bmean){
	double minNu = 0.1 * criticalNu(minEnergy, 1.0, Bmean);
	double maxNu = 100 * criticalNu(maxEnergy, 1.0, Bmean);

	double temp = maxNu / minNu;

	double factor = pow(maxNu / minNu, 1.0 / (Nnu - 1));

	nu[0] = minNu;
	for (int i = 1; i < Nnu; ++i) {
		nu[i] = nu[i - 1] * factor;
	}
}

double findEmissivityAt(double* nu, double* Inu, double currentNu, int Nnu) {
	if(currentNu <= nu[0]) {
		return Inu[0];
	}
	if(currentNu >= nu[Nnu - 1]) {
		return Inu[Nnu - 1];
	}
	int leftIndex = 0;
	int rightIndex = Nnu-1;
	while(rightIndex - leftIndex > 1){
		int currentIndex = (rightIndex + leftIndex)/2;
		if(nu[currentIndex] > currentNu){
			rightIndex = currentIndex;
		} else { 
			leftIndex = currentIndex;
		}
	}
	//return Inu[i - 1] * exp(log(Inu[i] / Inu[i - 1]) * ((currentNu - nu[i-1]) / (nu[i] - nu[i - 1])));
	return (Inu[leftIndex] *(nu[rightIndex] - currentNu) + Inu[rightIndex]*(currentNu - nu[leftIndex]))/ (nu[rightIndex] - nu[leftIndex]);
}

void evaluateSpectrum(double* nu, double* totalInu, double**** Inu, double**** Anu, double distance, int Nnu, double* X, double* Y, double* Z, double rmax){
	double dx = 2*rmax/Nx;
	for(int l = 0; l < Nnu; ++l){
		totalInu[l] = 0;
		for(int i = 0; i < Nx; ++i){
			for(int j = 0; j < Ny; ++j){
				double localInu = 0;
				for(int k = 0; k < Nz; ++k){
					double r = sqrt(X[i]*X[i] + Y[j]*Y[j] + Z[k]*Z[k]);
					/*if(r < rmax){
						printf("aaa\n");
					}*/
					double c = Anu[i][j][k][l]*localInu/(Inu[i][j][k][l]*dx*dx) - 1.0;
					localInu = localInu + Inu[i][j][k][l]*dx*dx*dx;
					if( Anu[i][j][k][l]*dx > 0.01){
						localInu = Inu[i][j][k][l]*dx*dx*(1.0 + c*exp(-Anu[i][j][k][l]*dx))/(Anu[i][j][k][l]);
					} else {
						localInu = localInu*(1.0 - Anu[i][j][k][l]*dx);
					}
				}
				totalInu[l] = totalInu[l] + localInu;
			}
		}
		totalInu[l] = totalInu[l]*1E26/(distance*distance);
	}
}

double evaluateOptimizationFunction(double Bfactor, double n, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double distance, double**** Inu, double**** Anu, double* X, double* Y, double* Z, double rmax) {
	double* totalInu = new double[Nnu];
	evaluateAllEmissivityAndAbsorption(nu, Inu, Anu, Nnu, Ee, Fe, Np, Nd, Bn, sintheta, thetaIndex, concentrations, n, Bfactor, X, Y, Z, rmax);
	evaluateSpectrum(nu, totalInu, Inu, Anu, distance, Nnu, X, Y, Z, rmax);

	//double I0 = findEmissivityAt(nu, totalInu, augx[0]*1E9, Nnu) - augy[0];
	double I1 = findEmissivityAt(nu, totalInu, augx[1]*1E9, Nnu) - augy[1];
	double I2 = findEmissivityAt(nu, totalInu, augx[2]*1E9, Nnu) - augy[2];
	double I3 = findEmissivityAt(nu, totalInu, augx[3]*1E9, Nnu) - augy[3];
	//double I4 = findEmissivityAt(nu, totalInu, augx[4]*1E9, Nnu) - augy[4];

	//return I0*I0 + I1*I1 + I2*I2 + I3*I3 + I4*I4;
	return fabs(I2) + fabs(I3) + fabs(I1);
}

void findMinParameters(double& Bfactor, double& N, double minLambda, double maxLambda, double gradB, double gradn, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double distance, double**** Inu, double**** Anu, double* X, double* Y, double* Z, double rmax) {
	if(maxLambda - minLambda < 0.01*maxLambda) {
		N = N - maxLambda*gradn;
		Bfactor = Bfactor - maxLambda*gradB;
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
	double f1 = evaluateOptimizationFunction(B1, N1, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
	double f2 = evaluateOptimizationFunction(B2, N2, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
	if(f1 < f2) {
		findMinParameters(Bfactor, N, minLambda, lambda2, gradB, gradn, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
	} else {
		findMinParameters(Bfactor, N, lambda1, maxLambda, gradB, gradn, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
	}
}

void findMinParameters(double& Bfactor, double& N, double gradB, double gradN, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double distance, double**** Inu, double**** Anu, double* X, double* Y, double* Z, double rmax) {
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
	findMinParameters(Bfactor, N, minLambda, maxLambda, gradB, gradN, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
}

void optimizeParameters(double& Bfactor, double& N, double* nu, double** Ee, double** Fe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double distance, double**** Inu, double**** Anu, double* X, double* Y, double* Z, double rmax, FILE* logFile) {
	double currentF = evaluateOptimizationFunction(Bfactor, N, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);

	for(int i = 0; i < 10; ++i) {
		double dx = min(Bfactor, N)/100;
		printf("optimiztion i = %d\n",i);
		fprintf(logFile, "optimiztion i = %d\n",i);
		fflush(logFile);
		double Fb = evaluateOptimizationFunction(Bfactor + dx, N, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
		double Fn = evaluateOptimizationFunction(Bfactor, N + dx, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
		double gradB = (Fb - currentF)/dx;
		double gradN = (Fn - currentF)/dx;
		findMinParameters(Bfactor, N, gradB, gradN, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
		currentF = evaluateOptimizationFunction(Bfactor, N, nu, Ee, Fe, Np, Nnu, Nd, Bn, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax);
		//random
		/*for(int j = 0; j < 10; ++j){
			double tempB = 5*b*uniformDistribution();
			double tempN = 5*n*uniformDistribution();
			double tempF = evaluateOptimizationFunction(B*tempB, N*tempN, Ee, Fe, Np, Nnu, minEnergy, maxEnergy, startElectronIndex, sinhi, localSize, normFactor);
			if(tempF < currentF){
				printf("random step\n");
				currentF = tempF;
				b = tempB;
				n = tempN;
			}
		}*/
	}
}

void evaluateFirstBapprox(double& Bfactor, double& concentration, double* nu, double** Ee, double** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double distance, double**** Inu, double**** Anu, double* X, double* Y, double* Z, double rmax, FILE* logFile){
	double* totalInu = new double[Nnu];
	
	double leftB = minB;
	double rightB = maxB;
	double localB = (leftB + rightB)/2;

	double nuMax = 0;
	int nuMaxIndex = 0;

	evaluateAllEmissivityAndAbsorption(nu, Inu, Anu, Nnu, Ee, dFe, Np, Nd, Bn, sintheta, thetaIndex, concentrations, concentration, localB, X, Y, Z, rmax);
	evaluateSpectrum(nu, totalInu, Inu, Anu, distance, Nnu, X, Y, Z, rmax);

	findMaxNu(nuMaxIndex, totalInu, Nnu);
	nuMax = nu[nuMaxIndex];
	int iterations = 0;
	bool converges = false;
	while(!converges) {
		iterations++;
		if(nuMax > augmaxx*1E9) {
			rightB = localB;
			localB = (leftB + rightB)/2;
			evaluateAllEmissivityAndAbsorption(nu, Inu, Anu, Nnu, Ee, dFe, Np, Nd, Bn, sintheta, thetaIndex, concentrations, concentration, localB, X, Y, Z, rmax);
			evaluateSpectrum(nu, totalInu, Inu, Anu, distance, Nnu, X, Y, Z, rmax);
			findMaxNu(nuMaxIndex, totalInu, Nnu);
			//findMaxNu(nuMaxIndex, doplerInu, Nnu);
		} else {
			leftB = localB;
			localB = (leftB + rightB)/2;
			evaluateAllEmissivityAndAbsorption(nu, Inu, Anu, Nnu, Ee, dFe, Np, Nd, Bn, sintheta, thetaIndex, concentrations, concentration, localB, X, Y, Z, rmax);
			evaluateSpectrum(nu, totalInu, Inu, Anu, distance, Nnu, X, Y, Z, rmax);
			findMaxNu(nuMaxIndex, totalInu, Nnu);
			//findMaxNu(nuMaxIndex, doplerInu, Nnu);
		}
		nuMax = nu[nuMaxIndex];
		if((rightB - leftB) < 0.01) {
			converges = true;
		}
	}

	Bfactor = localB;

	delete[] totalInu;
}


int main()
{
	double*** Bx;
	double*** By;
	double*** Bz;

	double*** concentrations;

	double*** B;
	double*** sintheta;
	int*** thetaIndex;

	double thetaObserv = 45;
	double cosThetaObserv = cos(thetaObserv);
	double sinThetaObserv = sin(thetaObserv);

	int Np= 200;
	int Nnu = 100;
	int Nd = 10;
	double** Fe;
	double** dFe;
	double** Ee;

	double**** Inu;
	double**** Anu;
	double* nu;

	double X[Nx];
	double Y[Ny];
	double Z[Nz];
	double distance = 40*3*1.0E24;

	const int Npoints = 4;

	double size[Npoints];

	size[3] = 2.30E17;
	size[2] = 1.10E17;
	size[1] = 6.8E16;
	size[0] = 3.1E16;

	double times[Npoints];

	times[0] = 0;
	times[1] = 2760000;
	times[2] = 5270400;
	times[3] = 10700000;

	double rmax = size[3];

	double dx = 2*rmax/Nx;

	//double rmin = dx*sqrt(3.0)/2;
	double rmin = 1E13;

	//double B0 = 1000;
	double B0 = 1.0;

	double rcorot = rmax/5;

	FILE* logFile = fopen(logFileName.c_str(), "w");

	printf("initialization\n");
	fprintf(logFile, "initialization\n");

	Bx = new double**[Nx];
	By = new double**[Nx];
	Bz = new double**[Nx];

	concentrations = new double**[Nx];

	B = new double**[Nx];
	sintheta = new double**[Nx];
	thetaIndex = new int**[Nx];

	for(int i = 0; i < Nx; ++i){
		double x = dx/2 - rmax + i*dx;
		X[i] = x;
		Bx[i] = new double*[Ny];
		By[i] = new double*[Ny];
		Bz[i] = new double*[Ny];
		concentrations[i] = new double*[Ny];
		B[i] = new double*[Ny];
		sintheta[i] = new double*[Ny];
		thetaIndex[i] = new int*[Ny];
		for(int j = 0; j < Ny; ++j){
			double y = dx/2 - rmax + j*dx;
			Y[j] = y;
			Bx[i][j] = new double[Nz];
			By[i][j] = new double[Nz];
			Bz[i][j] = new double[Nz];
			concentrations[i][j] = new double[Nz];
			B[i][j] = new double[Nz];
			sintheta[i][j] = new double[Nz];
			thetaIndex[i][j] = new int[Nz];
			for(int k = 0; k < Nz; ++k){
				double z = dx/2 - rmax + k*dx;
				Z[k] = z;
				double r = sqrt(x*x + y*y + z*z);

				concentrations[i][j][k] = 1.0;
				//todo
				double adiab = 5.0/3.0;
				concentrations[i][j][k] = pow(r/rmax, 3.0/(adiab - 1.0));

				double y1 = y;
				double z1 = z*cosThetaObserv + x*sinThetaObserv;
				double x1 = x*cosThetaObserv - z*sinThetaObserv;


				double costheta = z1/r;
				double sintheta = sqrt(1.0 - costheta*costheta);

				double sinphi = y1/sqrt(x1*x1 + y1*y1);
				double cosphi = x1/sqrt(x1*x1 + y1*y1);

				double Br = B0*costheta*sqr(rmin/r);
				double Bphi = B0*costheta*((r - rmin)/rcorot)*sqr(rmin/r)*sintheta;
				//double Bphi = 0;

				double Bx1 = Br*sintheta*cosphi - Bphi*sinphi;
				double By1 = Br*sintheta*sinphi + Bphi*cosphi;
				double Bz1 = Br*costheta;

				/*Bx[i][j][k] = Bx1*cosThetaObserv + Bz1*sinThetaObserv;
				By[i][j][k] = By1;
				Bz[i][j][k] = Bz1*cosThetaObserv - Bx1*sinThetaObserv;*/
				Bx[i][j][k] = B0;
				By[i][j][k] = 0;
				Bz[i][j][k] = 0;
			}
		}
	}

	printf("evaluating turbulent field\n");
	fprintf(logFile, "evaluating turbulent field\n");

	srand(time(NULL));

	double kmin = 2*pi*10/rmax;
	double dk = kmin;
	double kmax = Nk*dk;
	double turbNorm = evaluateTurbNorm(kmax, Nk, B0*rmin/rmax, 0.9);

	/*for(int ki = 0; ki < Nk; ++ki){
		printf("%d\n", ki);
		for(int kj = 0; kj < Nk; ++kj){
			for(int kk = 0; kk < Nk; ++kk){
				if ((ki + kj + kk) > 0) {
					double phase1 = 2*pi*uniformDistribution();
					double phase2 = 2*pi*uniformDistribution();


					double kx = ki*dk;
					double ky = kj*dk;
					double kz = kk*dk;

					double kw = sqrt(kx*kx + ky*ky + kz*kz);
					double kyz = sqrt(ky*ky + kz*kz);
					double cosTheta = kx/kw;
					double sinTheta = kyz/kw;
					double cosPhi = 1.0;
					double sinPhi = 0;
					if(kj + kk > 0) {
						cosPhi = ky/kyz;
						sinPhi = kz/kyz;
					} else {
						cosPhi = 1.0;
						sinPhi = 0.0;
					}

					double Bturbulent = evaluateTurbB(kx, ky, kz, turbNorm);

					for(int i = 0; i < Nx; ++i){
						for(int j = 0; j < Ny; ++j){
							for(int k = 0; k < Nz; ++k){

								double kmultr = kx*X[i] + ky*Y[j] + kz*Z[k];
								double localB1 = Bturbulent*sin(kmultr + phase1);
								double localB2 = Bturbulent*sin(kmultr + phase2);

								Bx[i][j][k] = Bx[i][j][k] - localB1*sinTheta;
								By[i][j][k] = By[i][j][k] + (localB1*cosTheta*cosPhi - localB2*sinPhi);
								Bz[i][j][k] = Bz[i][j][k] + (localB1*cosTheta*sinPhi + localB2*cosPhi);
							}
						}
					}
				}

			}
		}
	}*/

	FILE* bFile = fopen(BFileName.c_str(), "w");
	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			for(int k = 0; k < Nz; ++k){
				fprintf(bFile, "%g %g %g\n", Bx[i][j][k], By[i][j][k], Bz[i][j][k]);
			}
		}
	}
	fclose(bFile);

	evaluateOrientationParameters(B, sintheta, thetaIndex, Bx, By, Bz, X, Y, Z, Nd);

	printf("reading input\n");
	fprintf(logFile, "reading input\n");

	Fe = new double*[Nd];
	dFe = new double*[Nd];
	Ee = new double*[Nd];
	for(int j = 0; j < Nd; ++j){
		std::string fileNumber = convertIntToString(j);
		FILE* inputPe = fopen((fileNameP + fileNumber + ".dat").c_str(), "r");
		FILE* inputFe = fopen((fileNameF + fileNumber + ".dat").c_str(), "r");
		Fe[j] = new double[Np];
		dFe[j] = new double[Np];
		Ee[j] = new double[Np];
		double u;
		fscanf(inputPe, "%lf", &u);
		fscanf(inputFe, "%lf", &Fe[j][0]);
		Ee[j][0] = massElectron*speed_of_light2;
		Fe[j][0] = 0;
		dFe[j][0] = 0;
		for (int i = 1; i < Np; ++i) {
			fscanf(inputPe, "%lf", &u);
			fscanf(inputFe, "%lf", &Fe[j][i]);

			//if( u < 3000){
				Ee[j][i] = sqrt(u * u  + 1)*massElectron*speed_of_light2;
				//maxEnergy = Ee[i];
				Fe[j][i] = Fe[j][i] * Ee[j][i]/ (u * u * u * speed_of_light4 *massElectron*massElectron);
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
			//} else {
			//	Fe[i] = 0;
			//}
		
		}

		double norm = 0;
		for (int i = 1; i < Np; ++i) {
			norm = norm + Fe[j][i] * (Ee[j][i] - Ee[j][i - 1]);
		}
		for (int i = 0; i < Np; ++i) {
			Fe[j][i] = Fe[j][i] / norm;
		}
	}

	printf("evaluating local emissivity\n");
	fprintf(logFile, "evaluating local emissivity\n");

	nu = new double[Nnu];
	Inu = new double***[Nx];
	Anu = new double***[Nx];
	for(int i = 0; i < Nx; ++i){
		Inu[i] = new double**[Ny];
		Anu[i] = new double**[Ny];
		for(int j = 0; j < Ny; ++j){
			Inu[i][j]= new double*[Nz];
			Anu[i][j] = new double*[Nz];
			for(int k = 0; k < Nz; ++k){
				Inu[i][j][k] = new double[Nnu];
				Anu[i][j][k] = new double[Nnu];
			}
		}
	}

	//todo chose B
	evaluateNu(nu, Nnu, 1.1*massElectron*speed_of_light2, 1000*massElectron*speed_of_light2, 3*B0*rmin/rmax);

	/////////////////
	//todo concentration!!
	double concentration = 1.0;
	double Bfactor = 0.15;
	////////////////////
	printf("optimizing parameters\n");
	fprintf(logFile, "optimizing parameters\n");
	fflush(logFile);
	evaluateFirstBapprox(Bfactor, concentration, nu, Ee, dFe, Np, Nnu, Nd, B, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax, logFile);
	optimizeParameters(Bfactor, concentration, nu, Ee, dFe, Np, Nnu, Nd, B, sintheta, thetaIndex, concentrations, distance, Inu, Anu, X, Y, Z, rmax, logFile);
	///////////////////
	//concentration = 1.0;
	//Bfactor = 0.15;
	evaluateAllEmissivityAndAbsorption(nu, Inu, Anu, Nnu, Ee, dFe, Np, Nd, B, sintheta, thetaIndex, concentrations, concentration, Bfactor, X, Y, Z, rmax);

	printf("integrating fields\n");
	fprintf(logFile, "integrating Fields\n");
	double* totalInu = new double[Nnu];
	printf("Bfactor = %g, n = %g\n", Bfactor, concentration);
	fprintf(logFile, "Bfactor = %g, n = %g\n", Bfactor, concentration);
	fflush(logFile);

	evaluateSpectrum(nu, totalInu, Inu, Anu, distance, Nnu, X, Y, Z, rmax);

	//////////
	printf("outputing\n");
	fprintf(logFile, "ooutputing\n");
	//char* number = new char[100];
	//itoa(0, number, 10);
	//delete[] number;
	//std::string fileNumber = std::string(number);
	FILE* output = fopen(outputfileName.c_str(), "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(output, "%g %g\n", nu[i]/1E9, totalInu[i]);
	}

	fclose(output);

	////////////////////

	printf("deleting arrays\n");
	fprintf(logFile, "deleting arrays\n");
	delete[] totalInu;

	for(int j = 0; j < Nd; ++j){
		delete[] Fe[j];
		delete[] Ee[j];
	}
	delete[] Fe;
	delete[] Ee;

	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			for(int k = 0; k < Nz; ++k){
				delete[] Inu[i][j][k];
				delete[] Anu[i][j][k];
			}
			delete[] Inu[i][j];
			delete[] Anu[i][j];
			delete[] Bx[i][j];
			delete[] By[i][j];
			delete[] Bz[i][j];
			delete[] concentrations[i][j];
			delete[] B[i][j];
			delete[] sintheta[i][j];
			delete[] thetaIndex[i][j];
		}
		delete[] Inu[i];
		delete[] Anu[i];
		delete[] Bx[i];
		delete[] By[i];
		delete[] Bz[i];
		delete[] concentrations[i];
		delete[] B[i];
		delete[] sintheta[i];
		delete[] thetaIndex[i];
	}
	delete[] Inu;
	delete[] Anu;
	delete[] Bx;
	delete[] By;
	delete[] Bz;
	delete[] concentrations;
	delete[] B;
	delete[] sintheta;
	delete[] thetaIndex;
	
	fclose(logFile);

	return 0;
}

