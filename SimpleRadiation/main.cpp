// SphericalTurbulence.cpp : Defines the entry point for the console application.
//
#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include <time.h>
#include <cmath>
#include <string>
#include <omp.h>

#include "startparameters.h"
#include "constants.h"
#include "util.h"
#include "spectrum.h"
#include "optimization.h"

double phiValue[Nphi];
double sinPhiValue[Nphi];
double cosPhiValue[Nphi];


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
				if((i + j + k) > 1){
					double kx = i*dk;
					double ky = j*dk;
					double kz = k*dk;
					double B = evaluateTurbB(kx, ky, kz, 1.0);

					sum = sum + B*B;
				}
			}
		}
	}
	return B0*sqrt(energyFraction/((1.0 - energyFraction)*sum));
}

void evaluateOrientationParameters(double** B, double** sintheta, int** thetaIndex, double** Bx, double** By, double** Bz, int Nd){
		for(int j = 0; j < Ntheta; ++j){
			for(int k = 0; k < Nphi; ++k){
				//double r = rmax;
				double Bxy = sqrt(Bx[j][k]*Bx[j][k] + By[j][k]*By[j][k]);
				B[j][k] = sqrt(Bz[j][k]*Bz[j][k] + Bx[j][k]*Bx[j][k] + By[j][k]*By[j][k]);
				double cosTheta = Bz[j][k]/B[j][k];
				double sinTheta = Bxy/B[j][k];
				if(sinTheta != sinTheta){
					printf("sintheta NaN\n");
				}
				sintheta[j][k] = sinTheta;

				double cosTheta1 = cosThetaValue[j];
				double sinTheta1 = sinThetaValue[j];
				
				double cosPhi1 = cosPhiValue[k];
				double sinPhi1 = sinPhiValue[k];

				double Br = Bx[j][k]*sinTheta1*cosPhi1 + 
						    By[j][k]*sinTheta1*sinPhi1 +
							Bz[j][k]*cosTheta1;

				double cosTheta2 = Br/B[j][k];
				//double sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);


				double theta2 = acos(cosTheta2);
				if(theta2 >= pi/2){
					theta2 = pi - theta2;
				}

				thetaIndex[j][k] = floor(theta2 / (pi/(2*Nd))); 
				if(thetaIndex[j][k] == Nd){
					printf("thetaIndex == Nd\n");
					printf("%lf\n", By[j][k]);
					thetaIndex[j][k] = Nd - 1;
				}
				//for debug
				//thetaIndex[j][k] = 3;
		}
	}
}

void evaluateOrientationParameters3d(double*** B, double*** sintheta, int*** thetaIndex, double*** Bx, double*** By, double*** Bz, int Nd){
	for(int i = 0; i < Nrho; ++i){
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				//double r = rmax;
				double Bxy = sqrt(Bx[i][j][k]*Bx[i][j][k] + By[i][j][k]*By[i][j][k]);
				B[i][j][k] = sqrt(Bz[i][j][k]*Bz[i][j][k] + Bx[i][j][k]*Bx[i][j][k] + By[i][j][k]*By[i][j][k]);
				double cosTheta = Bz[i][j][k]/B[i][j][k];
				double sinTheta = Bxy/B[i][j][k];
				if(sinTheta != sinTheta){
					printf("sintheta NaN\n");
				}
				sintheta[i][j][k] = sinTheta;

				double r = sqrt((i+0.5)*(i+0.5) + (k + 0.5 - Nz/2.0)*(k + 0.5 - Nz/2.0));
				double z = k + 0.5 - Nz/2.0;

				double cosTheta1 = z/r;
				double sinTheta1 = sqrt(1.0 - cosTheta1*cosTheta1);
				
				double cosPhi1 = cosPhiValue[j];
				double sinPhi1 = sinPhiValue[j];

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
				//thetaIndex[i][j][k] = 9;
			}
		}
	}
}

void evaluateOrientationParameters3dflat(double*** B, double*** sintheta, int*** thetaIndex, double*** Bx, double*** By, double*** Bz, int Nd){
	for(int i = 0; i < Nrho; ++i){
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				//double r = rmax;
				double Bxy = sqrt(Bx[i][j][k]*Bx[i][j][k] + By[i][j][k]*By[i][j][k]);
				B[i][j][k] = sqrt(Bz[i][j][k]*Bz[i][j][k] + Bx[i][j][k]*Bx[i][j][k] + By[i][j][k]*By[i][j][k]);
				double cosTheta = Bz[i][j][k]/B[i][j][k];
				double sinTheta = Bxy/B[i][j][k];
				if(sinTheta != sinTheta){
					printf("sintheta NaN\n");
				}
				sintheta[i][j][k] = sinTheta;

				double theta2 = acos(cosTheta);
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
				//thetaIndex[i][j][k] = 9;
			}
		}
	}
}

void createNu(double* nu, int Nnu, double minNu, double maxNu){
	double temp = maxNu / minNu;

	double factor = pow(maxNu / minNu, 1.0 / (Nnu - 1));

	nu[0] = minNu;
	for (int i = 1; i < Nnu; ++i) {
		nu[i] = nu[i - 1] * factor;
	}
}

void evaluateNu(double* nu, int Nnu, double minEnergy, double maxEnergy, double Bmean){
	double minNu = 0.1 * criticalNu(minEnergy, 1.0, Bmean);
	double maxNu = 100 * criticalNu(maxEnergy, 1.0, Bmean);

	createNu(nu, Nnu, minNu, maxNu);
}

void initializeParker(double*** Bx, double*** By, double*** Bz){
	double thetaObserv = 0;
	double sinThetaObserv = sin(thetaObserv);
	double cosThetaObserv = cos(thetaObserv);
	double rcorot = Nrho/10.0;
	double rmin = Nrho/100.0;
	double Br1 = sqr(rmin/Nrho);
	double Bphi1 = ((Nrho - rmin)/rcorot)*sqr(rmin/Nrho);
	//norm to 1;
	double B0 = 1.0/sqrt(Br1*Br1 + Bphi1*Bphi1);
	for(int i = 0; i < Nrho; ++i){
		double rho = i + 0.5;
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				double z = k + 0.5 - Nz/2.0;
				double r = sqrt(rho*rho + z*z);

				double x = rho*cosPhiValue[j];
				double y = rho*sinPhiValue[j];

				double y1 = y;
				double z1 = z*cosThetaObserv + x*sinThetaObserv;
				double x1 = x*cosThetaObserv - z*sinThetaObserv;

				double costheta = z1/r;
				double sintheta = sqrt(1.0 - costheta*costheta);

				double sinphi = y1/sqrt(x1*x1 + y1*y1);
				double cosphi = x1/sqrt(x1*x1 + y1*y1);

				double Br = B0*costheta*sqr(rmin/r);
				double Bphi = B0*costheta*((r - rmin)/rcorot)*sqr(rmin/r)*sintheta;

				double Bx1 = Br*sintheta*cosphi - Bphi*sinphi;
				double By1 = Br*sintheta*sinphi + Bphi*cosphi;
				double Bz1 = Br*costheta;

				Bx[i][j][k] = Bx1*cosThetaObserv + Bz1*sinThetaObserv;
				By[i][j][k] = By1;
				Bz[i][j][k] = Bz1*cosThetaObserv - Bx1*sinThetaObserv;

				if(Bx[i][j][k] != Bx[i][j][k]){
					printf("Bx = NaN\n");
				}
			}
		}
	}

	double maxB = 0;
	//for(int i = 0; i < Nrho; ++i){
	int i = Nrho - 1;
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				double B = sqrt(Bx[i][j][k]*Bx[i][j][k] + By[i][j][k]*By[i][j][k] + Bz[i][j][k]*Bz[i][j][k]);
				if(B > maxB) {
					maxB = B;
				}
			}
		}
	//}
	printf("maxB in parker = %g\n", maxB);
}

void evaluateAngleWeights(double* weights, int Nd, double B0, double sintheta, double turbulenceFraction){
	const int Nx = 1000;
	const int Nk = 10;
	double x[Nx];
	double y[Nx];
	double z[Nx];
	double Bx[Nx];
	double By[Nx];
	double Bz[Nx];

	for( int i = 0; i < Nx; ++i){
		Bx[i] = B0*sqrt(1.0 - sintheta*sintheta);
		Bz[i] = B0*sintheta;
		By[i] = 0;
		x[i] = uniformDistribution();
		y[i] = uniformDistribution();
		z[i] = uniformDistribution();
	}

	double turbNorm = evaluateTurbNorm(2*pi*10, Nk, B0, turbulenceFraction);

	for(int i = 0; i < Nd; ++i){
		weights[i] = 0;
	}

	double dk = 2*pi;

	for(int ki = 0; ki < Nk; ++ki){
		for(int kj = 0; kj < Nk; ++kj){
			for(int kk = 0; kk < Nk; ++kk){
				if ((ki + kj + kk) > 1) {
					double phase1 = 2*pi*uniformDistribution();
					double phase2 = 2*pi*uniformDistribution();


					double kx = ki*dk;
					double ky = kj*dk;
					double kz = kk*dk;

					double kw = sqrt(kx*kx + ky*ky);
					double kxy = sqrt(ky*ky + kx*kx);
					double cosTheta = 0;
					double sinTheta = 1.0;
					double cosPhi = 1.0;
					double sinPhi = 0;
					if(kj + ki > 0) {
						cosPhi = kx/kxy;
						sinPhi = ky/kxy;
					} else {
						cosPhi = 1.0;
						sinPhi = 0.0;
					}

					double Bturbulent = evaluateTurbB(kx, ky, kz, turbNorm);

					for(int i = 0; i < Nx; ++i){
						double kmultr = (kx*x[i] + ky*y[i] + kz*z[i]);
						double localB1 = Bturbulent*sin(kmultr + phase1);
						double localB2 = Bturbulent*sin(kmultr + phase2);
						//localB1 = 0;
						//localB2 = localB2*sqrt(2.0);

						Bz[i] = Bz[i] - localB1*sinTheta;
						Bx[i] = Bx[i] + (localB1*cosTheta*cosPhi - localB2*sinPhi);
						By[i] = By[i] + (localB1*cosTheta*sinPhi + localB2*cosPhi);
					}
				}
			}
		}
	}

	for(int i = 0; i < Nx; ++i){
		double B = sqrt(Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i]);
		double theta1 = acos(Bx[i]/B);
		if(theta1 >= pi/2){
			theta1 = pi - theta1;
		}

		int thetaIndex = floor(theta1 / (pi/(2*Nd))); 
		weights[thetaIndex] += 1.0/Nx;
		if(weights[thetaIndex] != weights[thetaIndex]){
			printf("weight = NaN\n");
			exit(0);
		}
	}

	for(int i = 0; i < Nd; ++i){
		if(weights[i] != weights[i]){
			printf("weight = NaN\n");
			exit(0);
		}
		if(0*weights[i] != 0*weights[i]){
			printf("weight = NaN\n");
			exit(0);
		}
	}
}

double findFeAt(double* Ee, double* Fe, double currentE, int Np) {
	if(currentE <= Ee[0]) {
		return Fe[0];
	}
	if(currentE >= Ee[Np - 1]) {
		return Fe[Np - 1];
	}
	int leftIndex = 0;
	int rightIndex = Np-1;
	while(rightIndex - leftIndex > 1){
		int currentIndex = (rightIndex + leftIndex)/2;
		if(Ee[currentIndex] > currentE){
			rightIndex = currentIndex;
		} else { 
			leftIndex = currentIndex;
		}
	}
	//return Inu[i - 1] * exp(log(Inu[i] / Inu[i - 1]) * ((currentNu - nu[i-1]) / (nu[i] - nu[i - 1])));
	double result = (Fe[leftIndex] *(Ee[rightIndex] - currentE) + Fe[rightIndex]*(currentE - Ee[leftIndex]))/ (Ee[rightIndex] - Ee[leftIndex]);
	if(result != result){
		printf("result = NaN\n");
		exit(0);
	}
	if(0*result != 0*result){
		printf("result = QNaN\n");
		exit(0);
	}
	return result;
}

double evaluateR(double nupeak, double fpeak, double d, double fraction, double epsilone, double epsilonB, double p){
	double c1 = 6.265E18;
	double c5 = pacholczykC5[5];
	double c6 = pacholczykC6[5];

	double a1 = pow(c6, p+5);
	double a2 = pow(fpeak, p+6);
	double a3 = pow(d, p + 6);

	double d0 = 3.08E24;
	double nu0 = 5E9;
	double f0 = 1E-23;

	double a = 6*(pow(c6, p+5)*pow(d, p + 6))* (pow(fpeak, p+6)*pow(d, p + 6));
	double b = (epsilone/epsilonB)*fraction*(p-2)*pow(pi, p+5)*pow(c5, p+6)*pow(massElectron*speed_of_light2, p-2);

	double result = pow(a/b, 1.0/(2*p + 13))*2*c1/nupeak;

	//double result = 5.1E15*pow(epsilone/epsilonB, -1.0/(2*p +13))*pow(fpeak/f0, (p+6)/(2*p +13))*pow(d/d0, (2*p+12)/(2*p +13))*(nu0/nupeak);

	return result;
}

double evaluateB(double nupeak, double fpeak, double d, double fraction, double epsilone, double epsilonB, double p){
	double c1 = 6.265E18;
	double c5 = pacholczykC5[5];
	double c6 = pacholczykC6[5];

	double d0 = 3.08E24;
	double nu0 = 5E9;
	double f0 = 1E-23;

	double a = 36*pi*pi*pi*c5;
	double b = pow(epsilone/epsilonB, 2)*fraction*fraction*(p-2.0)*(p-2.0)*c6*c6*c6*pow(massElectron*speed_of_light2, 2*(p-2.0))*fpeak*d*d;

	double result = pow(a/b, 2.0/(2*p + 13))*nupeak/(2*c1);


	//double result = 0.3*pow(epsilone/epsilonB, -4.0/(2*p +13))*pow(fpeak/f0, -2.0/(2*p +13))*pow(d/d0, -4.0/(2*p +13))*(nupeak/nu0);

	return result;
}


int main()
{
	double** Bx;
	double** By;
	double** Bz;

	double** area;
	double** length;

	double** concentrations;

	//double** B;
	//double** sintheta;
	//int** thetaIndex;

	double thetaObserv = 0;
	double cosThetaObserv = cos(thetaObserv);
	double sinThetaObserv = sin(thetaObserv);

	//const int Np= 200;
	const int Np= 298;
	const int Nnu = 100;
	const int Nnu1 = 5;
	const int Ndist = 10;

	double** Fe;
	double** dFe;
	double** Ee;


	double* Rho;
	double* Phi;

	double*** B3d;
	double*** Bx3d;
	double*** By3d;
	double*** Bz3d;
	double*** concentrations3d;
	double*** area3d;
	double*** length3d;
	int*** thetaIndex3d;
	double*** sintheta3d;

	double z = 0.2433;
	double d = 816*3.08*1.0E24; //angular diameter distance
	double fraction = 0.5;
	double fpeak = 0.68E-26;
	double nupeak = 22.0E9;
	double Time = 16.5*24*3600;
	double epsilone = 0.33;
	double epsilonB = 0.33;

	//double epsilone = 0.0154;
	//double epsilonB = 0.04;

	//fpeak = 1E-23;
	//nupeak = 5E9;
	//d = 3.085*1.0E24;
	//fraction = 0.5;

	double p = 3.0;

	double B1 = evaluateB(nupeak, fpeak, d, fraction, epsilone, epsilonB, p);
	double R1 = evaluateR(nupeak, fpeak, d, fraction, epsilone, epsilonB, p);
	double t1 = 71*24*3600;
	double v1 = (1+z)*R1/t1;
	double beta1 = v1/speed_of_light;
	double n1 = (B1*B1/(8*pi))/(epsilonB*v1*v1*massProtonReal);

	//double rmin = dx*sqrt(3.0)/2;
	double rmax = size[0];
	double rmin = 1E13;	

	//double B0 = 1000;
	double B0 = 1.0;

	double rcorot = rmax/5;

	FILE* logFile = fopen(logFileName.c_str(), "w");

	printf("initialization\n");
	fprintf(logFile, "initialization\n");
	fflush(logFile);

	Rho = new double[Nrho];
	Phi = new double[Nphi];

	area = new double*[Nrho];
	length = new double*[Nrho];
	double dr = rmax/Nrho;
	for(int i = 0; i < Nrho; ++i){
		Rho[i] = dr/2 + i*dr;
		area[i] = new double[Nphi];
		length[i] = new double[Nphi];
		for(int j = 0; j < Nphi; ++j){
			Phi[j] = dphi/2 + j*dphi;
		}
	}

	for(int i = 0; i < Nphi; ++i){
		phiValue[i] = dphi/2.0 + i*dphi;
		sinPhiValue[i] = sin(phiValue[i]);
		cosPhiValue[i] = cos(phiValue[i]);
	}

	/*Bx = new double*[Ntheta];
	By = new double*[Ntheta];
	Bz = new double*[Ntheta];

	concentrations = new double*[Ntheta];

	B = new double*[Ntheta];
	sintheta = new double*[Ntheta];
	thetaIndex = new int*[Ntheta];

	for(int j = 0; j < Ntheta; ++j){
		Bx[j] = new double[Nphi];
		By[j] = new double[Nphi];
		Bz[j] = new double[Nphi];
		concentrations[j] = new double[Nphi];
		B[j] = new double[Nphi];
		sintheta[j] = new double[Nphi];
		thetaIndex[j] = new int[Nphi];
		for(int k = 0; k < Nphi; ++k){
			double z = rmax*cosThetaValue[j];
			double rho = rmax*sinThetaValue[j];

			double r = rmax;
				
			//change later
			concentrations[j][k] = 1.0;

			double x = rho*cosPhiValue[j];
			double y = rho*sinPhiValue[j];

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

			Bx[j][k] = Bx1*cosThetaObserv + Bz1*sinThetaObserv;
			By[j][k] = By1;
			Bz[j][k] = Bz1*cosThetaObserv - Bx1*sinThetaObserv;
			Bx[j][k] = B0;
			By[j][k] = 0;
			Bz[j][k] = 0;
		}
	}*/

	area3d = new double**[Nrho];
	length3d = new double**[Nrho];
	B3d = new double**[Nrho];
	Bx3d = new double**[Nrho];
	By3d = new double**[Nrho];
	Bz3d = new double**[Nrho];
	concentrations3d = new double**[Nrho];
	sintheta3d = new double**[Nrho];
	thetaIndex3d = new int**[Nrho];
	for(int i = 0; i < Nrho; ++i){
		area3d[i] = new double*[Nphi];
		length3d[i] = new double*[Nphi];
		B3d[i] = new double*[Nphi];
		Bx3d[i] = new double*[Nphi];
		By3d[i] = new double*[Nphi];
		Bz3d[i] = new double*[Nphi];
		concentrations3d[i] = new double*[Nphi];
		sintheta3d[i] = new double*[Nphi];
		thetaIndex3d[i] = new int*[Nphi];
		for(int j = 0; j < Nphi; ++j){
			area3d[i][j] = new double[Nz];
			length3d[i][j] = new double[Nz];
			B3d[i][j] = new double[Nz];
			Bx3d[i][j] = new double[Nz];
			By3d[i][j] = new double[Nz];
			Bz3d[i][j] = new double[Nz];
			concentrations3d[i][j] = new double[Nz];
			sintheta3d[i][j] = new double[Nz];
			thetaIndex3d[i][j] = new int[Nz];
			for(int k = 0; k < Nz; ++k){
				area3d[i][j][k] = 0;
				length3d[i][j][k] = 0;
				double r = sqrt(1.0*(i+0.5)*(i+0.5) + (k + 0.5 - Nz/2.0)*(k + 0.5 - Nz/2.0));
				double rmax = Nrho;
				B3d[i][j][k] = 1.0;
				Bx3d[i][j][k] = 1.0;
				By3d[i][j][k] = 0.0;
				Bz3d[i][j][k] = 0.0;
				if(geometry == SPHERICAL) {
					//concentrations3d[i][j][k] = 1.0*sqr(rmax/r);
					//concentrations3d[i][j][k] = 1.0*sqr(r/rmax);
					concentrations3d[i][j][k] = 1.0;
				} else {
					concentrations3d[i][j][k] = 1.0;
				}
				//sintheta3d[i][j][k] = 1.0;
				//thetaIndex3d[i][j][k] = 9;
			}
		}
	}

	if(parker){
		initializeParker(Bx3d, By3d, Bz3d);
	}


	printf("evaluating turbulent field\n");
	fprintf(logFile, "evaluating turbulent field\n");
	fflush(logFile);

	srand(time(NULL));
	int randomSeed = rand();
	randomSeed = 10;
	printf("random seed = %d\n", randomSeed);
	fprintf(logFile, "random seed = %d\n", randomSeed);
	fflush(logFile);
	srand(randomSeed);

	double kmin = 2*pi*2/(0.5);
	double dk = kmin;
	double kmax = Nk*dk;
	double turbulenceFraction = 0.000000001;
	double turbNorm = evaluateTurbNorm(kmax, Nk, 1.0, turbulenceFraction);

	/*if(turbulence){
		for(int i = 0; i < Nrho; ++i){
			for(int j = 0; j < Nphi; ++j){
				for(int k = 0; k < Nz; ++k){
					Bx3d[i][j][k] = Bx3d[i][j][k]*sqrt(1.0 - turbulenceFraction);
					By3d[i][j][k] = By3d[i][j][k]*sqrt(1.0 - turbulenceFraction);
					Bz3d[i][j][k] = Bz3d[i][j][k]*sqrt(1.0 - turbulenceFraction);
					B3d[i][j][k] = B3d[i][j][k]*sqrt(1.0 - turbulenceFraction);
				}
			}
		}

		const int tempNr = 10000;
		double tempB[tempNr];
		double tempBx[tempNr];
		double tempBy[tempNr];
		double tempBz[tempNr];
		double lambda[tempNr];
		double delta[tempNr];
		double tempr[tempNr];
		for(int i = 0; i < tempNr; ++i){
			tempr[i] = Nrho*(i+0.5)/tempNr;
			tempB[i] = 0;
			tempBx[i] = 0;
			tempBy[i] = 0;
			tempBz[i] = 0;
			lambda[i] = 0;
			delta[i] = 0;
		}

		for(int ki = 0; ki < Nk; ++ki){
			printf("%d\n", ki);
			for(int kj = 0; kj < Nk; ++kj){
				//for(int kk = 0; kk < Nk; ++kk){
					if ((ki + kj) > 4) {
						double phase1 = 2*pi*uniformDistribution();
						double phase2 = 2*pi*uniformDistribution();


						double kx = ki*dk;
						double ky = kj*dk;
						//double kz = kk*dk;

						double kw = sqrt(kx*kx + ky*ky);
						double kxy = sqrt(ky*ky + kx*kx);
						double cosTheta = 0;
						double sinTheta = 1.0;
						double cosPhi = 1.0;
						double sinPhi = 0;
						if(kj + ki > 0) {
							cosPhi = kx/kxy;
							sinPhi = ky/kxy;
						} else {
							cosPhi = 1.0;
							sinPhi = 0.0;
						}

						double Bturbulent = evaluateTurbB(kx, ky, 0, turbNorm);

						for(int i = 0; i < tempNr; ++i){
							double x = tempr[i];
							double kmultr = kx*x;
							double localB1 = Bturbulent*sin(kmultr + phase1);
							double localB2 = Bturbulent*sin(kmultr + phase2);
							localB1 = 0;
							localB2 = localB2*sqrt(2.0);

							tempBz[i] = tempBz[i] - localB1*sinTheta;
							tempBx[i] = tempBx[i] + (localB1*cosTheta*cosPhi - localB2*sinPhi);
							tempBy[i] = tempBy[i] + (localB1*cosTheta*sinPhi + localB2*cosPhi);
						}

						for(int i = 0; i < Nrho; ++i){
							for(int j = 0; j < Nphi; ++j){
								for(int k = 0; k < Nz; ++k){
									double x = (i+0.5)*cos((2*pi/Nphi)*(j+0.5));
									double y = (i+0.5)*sin((2*pi/Nphi)*(j+0.5));
									double z = (k - Nz/2.0 +0.5);

									double kmultr = kx*x + ky*y;
									double localB1 = Bturbulent*sin(kmultr + phase1)*cos(pi*z/(2*Nz));
									double localB2 = Bturbulent*sin(kmultr + phase2)*cos(pi*z/(2*Nz));
									localB1 = 0;
									localB2 = localB2*sqrt(2.0);

									Bz3d[i][j][k] = Bz3d[i][j][k] - localB1*sinTheta;
									Bx3d[i][j][k] = Bx3d[i][j][k] + (localB1*cosTheta*cosPhi - localB2*sinPhi);
									By3d[i][j][k] = By3d[i][j][k] + (localB1*cosTheta*sinPhi + localB2*cosPhi);
								}
							}
						}
					//}
				}
			}
		}

		FILE* tempBfile = fopen("turbulentB.dat","w");
		for(int i = 0; i < tempNr; ++i){
			tempB[i] = sqrt(tempBx[i]*tempBx[i] + tempBy[i]*tempBy[i] + tempBz[i]*tempBz[i]);
			lambda[i] = atan2(tempBy[i],tempBx[i])*180/pi;
			delta[i] = acos(tempBz[i]/tempB[i])*180/pi;
			fprintf(tempBfile, "%g %g %g %g\n", tempr[i], tempB[i], lambda[i], delta[i]);
		}
		fclose(tempBfile);
	}*/


	if(geometry == SPHERICAL){
		evaluateOrientationParameters3d(B3d, sintheta3d, thetaIndex3d, Bx3d, By3d, Bz3d, Ndist);
	} else {
		evaluateOrientationParameters3dflat(B3d, sintheta3d, thetaIndex3d, Bx3d, By3d, Bz3d, Ndist);
	}

	double**** weights = new double***[Nrho];

	for(int i = 0; i < Nrho; ++i){
		weights[i] = new double**[Nphi];
		for(int j = 0; j < Nphi; ++j){
			weights[i][j] = new double*[Nz];
			for(int k = 0; k < Nz; ++k){
				weights[i][j][k] = new double[Ndist];
				if(turbulence){
					evaluateAngleWeights(weights[i][j][k], Ndist, 1.0, sintheta3d[i][j][k], 0.9);
				} else {
					for(int l = 0; l < Ndist; ++l){
						weights[i][j][k][l] = 0;
					}
					weights[i][j][k][thetaIndex3d[i][j][k]] = 1.0;
				}
			}
		}
	}

	/*FILE* bFile = fopen(BFileName.c_str(), "w");
	for(int j = 0; j < Ntheta; ++j){
		for(int k = 0; k < Nphi; ++k){
			fprintf(bFile, "%g %g %g\n", Bx[j][k], By[j][k], Bz[j][k]);
		}
	}
	fclose(bFile);*/


	/*double sigma = 0.04;
	double tempConcentration = sqr(B[Ntheta/2][0])/(sigma*4*pi*massProtonReal*speed_of_light2);
	for(int i = 0; i < Ntheta; ++i){
		for(int j = 0; j < Nphi; ++j){
				concentrations[i][j] = tempConcentration;
		}
	}*/

	printf("reading input\n");
	fprintf(logFile, "reading input\n");
	fflush(logFile);
	int powerlawStart[Ndist] = {140,140,140,140,140,140,140,140,140,140};


	Fe = new double*[Ndist];
	dFe = new double*[Ndist];
	Ee = new double*[Ndist];
	for(int j = 0; j < Ndist; ++j){
		//std::string fileNumber = convertIntToString(j);
		std::string fileNumber = convertIntToString(0);
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

		if(input == TRISTAN){
			for (int i = 1; i < Np; ++i) {
				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				u = u*realMassRelationSqrt/massRelationSqrt;
				//if( u < 3000){
					double gamma = sqrt(u * u  + 1);
					Ee[j][i] = sqrt(u * u  + 1)*massElectron*speed_of_light2;
					//maxEnergy = Ee[i];
					Fe[j][i] = Fe[j][i] * Ee[j][i]/ (u * u * u * speed_of_light4 *massElectron*massElectron);
					/*if(gamma >= 1.0){
						Fe[j][i] = 1.0/pow(Ee[j][i],3);
					} else {
						Fe[j][i] = 0;
					}*/
					dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
				//} else {
				//	Fe[i] = 0;
				//}


				//Fe[j][i] = exp(-gamma/thetae)*gamma*u;
				//dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
			}
		} else if(input == SMILEI){
			for (int i = 1; i < Np; ++i) {
				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				//todo massRelationSqrt?
				double gamma = u*realMassRelationSqrt/massRelationSqrt + 1;
				//if( u < 3000){
				Ee[j][i] = gamma*massElectron*speed_of_light2;
					//maxEnergy = Ee[i];
					Fe[j][i] = Fe[j][i] / (massElectron*speed_of_light2);
					/*if(gamma >= 1.0){
						Fe[j][i] = 1.0/pow(Ee[j][i],3);
					} else {
					Fe[j][i] = 0;
					}*/
					dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
				//} else {
				//	Fe[i] = 0;
				//}


				//Fe[j][i] = exp(-gamma/thetae)*gamma*u;
				//dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
			}
		} else if(input == MAXWELL){
			for (int i = 1; i < Np; ++i) {
				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				double Te = (massRelationSqrt/realMassRelationSqrt)*2.4E11;
				//double Te = 1.6*1E10;
				double thetae = kBoltzman*Te/(massElectron*speed_of_light2);

				//todo massRelationSqrt?
				double gamma = u*realMassRelationSqrt/massRelationSqrt + 1;
				Ee[j][i] = gamma*massElectron*speed_of_light2;


				Fe[j][i] = exp(-gamma/thetae)*gamma*u;
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
			}
		} else if(input == POWERLAW) {
	
			int Npowerlaw = powerlawStart[j];
			Npowerlaw = 1;
			for (int i = 1; i < Np; ++i) {

				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				//todo massRelationSqrt?
				double gamma = u*realMassRelationSqrt/massRelationSqrt + 1;
				Ee[j][i] = gamma*massElectron*speed_of_light2;
				Fe[j][i] = Fe[j][i] / (massElectron*speed_of_light2);
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);

				double minGamma = 6.0;
				double power = 3.5;

				if(gamma >= minGamma){
					Fe[j][i] = 1.0/pow(Ee[j][i],power);
				} else {
					Fe[j][i] = 0;
				}
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
			}
			/*Ee[j][Np-1] = 1E8;
			double factor = pow(Ee[j][Np-1]/Ee[j][Npowerlaw-1], 1.0/(Np-Npowerlaw));
			//double gammae = -log(Fe[j][Npowerlaw-10]/Fe[j][Npowerlaw-1])/log(Ee[j][Npowerlaw-10]/Ee[j][Npowerlaw-1]);
			double gammae = 2.5;
			for(int i = Npowerlaw; i < Np; ++i){
				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				//todo massRelationSqrt?
				double gamma = u*realMassRelationSqrt/massRelationSqrt + 1;
				//Ee[j][i] = gamma*massElectron*speed_of_light2;
				Ee[j][i] = Ee[j][i-1]*factor;
				Fe[j][i] = Fe[j][Npowerlaw-1]*pow(Ee[j][Npowerlaw-1]/Ee[j][i], gammae);
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
			}*/
		} else if (input == COMBINED){
			for (int i = 1; i < Np; ++i) {
				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				//todo massRelationSqrt?
				double gamma = u*realMassRelationSqrt/massRelationSqrt + 1;
				//if( u < 3000){
				Ee[j][i] = gamma*massElectron*speed_of_light2;
				//maxEnergy = Ee[i];
				Fe[j][i] = Fe[j][i] / (massElectron*speed_of_light2);
				//if(i > 137){
				if((i > 137)&&(u < 500)){
					Fe[j][i] = Fe[j][137]*pow(Ee[j][i]/Ee[j][137], -3.5);
				}
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
				if (dFe[j][i] < 0) {
					printf("dFe < 0\n");
				}
			}

		}
		else if (input == MONTECARLO) {
			for (int i = 1; i < Np; ++i) {
				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				double gamma = u;
				//if( u < 3000){
				Ee[j][i] = gamma * massElectron * speed_of_light2;
				//maxEnergy = Ee[i];
				Fe[j][i] = Fe[j][i] / (massElectron * speed_of_light2);
				/*if(gamma >= 1.0){
					Fe[j][i] = 1.0/pow(Ee[j][i],3);
				} else {
				Fe[j][i] = 0;
				}*/
				dFe[j][i] = (Fe[j][i] / (4 * pi)) * (Ee[j][i] - Ee[j][i - 1]);
				//} else {
				//	Fe[i] = 0;
				//}


				//Fe[j][i] = exp(-gamma/thetae)*gamma*u;
				//dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
			}
		}
		else if (input == MONTECARLO_CUT) {
			for (int i = 1; i < Np; ++i) {
				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				double gamma = u;
				//if( u < 3000){
				Ee[j][i] = gamma * massElectron * speed_of_light2;
				//maxEnergy = Ee[i];
				Fe[j][i] = Fe[j][i] / (massElectron * speed_of_light2);
				if (gamma > 5E6) {
					Fe[j][i] = 0;
				}
				/*if(gamma >= 1.0){
					Fe[j][i] = 1.0/pow(Ee[j][i],3);
				} else {
				Fe[j][i] = 0;
				}*/
				dFe[j][i] = (Fe[j][i] / (4 * pi)) * (Ee[j][i] - Ee[j][i - 1]);
				//} else {
				//	Fe[i] = 0;
				//}


				//Fe[j][i] = exp(-gamma/thetae)*gamma*u;
				//dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
			}
		}

		fclose(inputPe);
		fclose(inputFe);

		double norm = 0;
		for (int i = 1; i < Np; ++i) {
			norm = norm + Fe[j][i] * (Ee[j][i] - Ee[j][i - 1]);
		}
		for (int i = 0; i < Np; ++i) {
			Fe[j][i] = Fe[j][i] / norm;
			dFe[j][i] = dFe[j][i]/norm;
		}
	}

	double* weightedEe = new double[Np];
	double**** weightedFe = new double***[Nrho];
	double minE = Ee[0][0];
	double maxE = Ee[0][Np-1];
	for(int i = 0; i < Ndist; ++i){
		if(Ee[i][0] < minE){
			minE = Ee[i][0];
		}
		if(Ee[i][Np-1] > maxE){
			maxE = Ee[i][Np-1];
		}
	}

	/*double factor = pow(maxE/minE, 1.0/(Np - 1));
	weightedEe[0] = minE;
	for(int i = 1; i < Np; ++i){
		weightedEe[i] = weightedEe[i-1]*factor;
	}
	weightedEe[Np-1] = maxE;*/
	for(int i = 0; i < Np; ++i){
		weightedEe[i] = Ee[0][i];
	}

	for(int i = 0; i < Nrho; ++i){
		weightedFe[i] = new double**[Nphi];
		for(int j = 0; j < Nphi; ++j){
			weightedFe[i][j] = new double*[Nz];
			for(int k = 0; k < Nz; ++k){
				weightedFe[i][j][k] = new double[Np];
				for(int l = 0; l < Np; ++l){
					weightedFe[i][j][k][l] = 0;
					for(int m = 0; m < Ndist; ++m){
						weightedFe[i][j][k][l] += weights[i][j][k][m]*findFeAt(Ee[m], Fe[m], weightedEe[l], Np);
						//weightedFe[i][j][k][l] = Fe[0][l];
						if(weightedFe[i][j][k][l] != weightedFe[i][j][k][l]){
							printf("WeightdFe = NaN\n");
							exit(0);
						}
					}
				}
			}
		}
	}


	for(int i = 0; i < Nrho; ++i){
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				double norm = 0;
				weightedFe[i][j][k][0] = 0;
				for(int l = 1; l < Np; ++l){
					/*norm += weightedFe[i][j][k][l]*(weightedEe[l] - weightedEe[l-1]);
					if(norm != norm){
						printf("norm = NaN\n");
						exit(0);
					}*/

					weightedFe[i][j][k][l] = (weightedFe[i][j][k][l]/(4*pi))*(weightedEe[l] - weightedEe[l-1]);
				}
			}
		}
	}

	printf("evaluating local emissivity\n");
	fprintf(logFile, "evaluating local emissivity\n");
	fflush(logFile);

	double* nu = new double[Nnu];



	//todo chose B
	/*double meanB = 0;
	int ncells = 0;
	for(int i = 0; i < Nr; ++i){
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				double r = sqrt(Rho[i]*Rho[i] + Z[k]*Z[k]);
				if((r < rmax) && (r > rmax*(1.0 - fractionSize))){
					ncells++;
					meanB = meanB + sqrt(Bx[i][j][k]*Bx[i][j][k] + By[i][j][k]*By[i][j][k] + Bz[i][j][k]*Bz[i][j][k]);
				}
			}
		}
	}
	meanB = meanB/ncells;*/

	//evaluateNu(nu, Nnu, 1.1*massElectron*speed_of_light2, 1000*massElectron*speed_of_light2, meanB);
	double numax = 10000*1E9;
	numax = 10000*1.6*1E-12/hplank;
	createNu(nu, Nnu, 0.1*1E9, 2E18);
	/*for(int i = 0; i < Nnu1; ++i){
		nu1[i] = aprx[i]*1.0E9;
	}*/

	/////////////////
	//todo concentration!!
	double sigma = 0.004;
	double tempConcentration = sqr(B3d[Ntheta/2][0][0])/(sigma*4*pi*massProtonReal*speed_of_light2);
	double concentration = 1.0*tempConcentration;
	concentration = 11;
	double Bfactor = 0.1;
	//double fractionSize = 1.0 - pow((3.0/(4.0*pi))*(1 - 0.5),1.0/3.0);
	double fractionSize = 0.001;
	double V0 = speed_of_light;
	double v = 0.1*speed_of_light;
	rmax = 3.4E16;
	////////////////////
	printf("initial optimizing parameters\n");
	fprintf(logFile, "initial optimizing parameters\n");
	fflush(logFile);
	printf("optimizing parameters\n");
	fprintf(logFile, "optimizing parameters\n");
	fflush(logFile);
	epsilonB = 0.1;
	v = 0.15*speed_of_light;


	//from Nayan day 138
	rmax = 5.38E16;
	Bfactor = 0.065;
	double Mdot = 4.82E-5;
	concentration = n1;
	fractionSize = 1.0 - pow(0.5, 1.0/3.0);

	rmax = R1;
	Bfactor = B1;
	concentration = n1;
	fractionSize = 1.0 - pow(0.5, 1.0/3.0);
	//fractionSize = 0.5;
	v = v1;

	/*double t = 82.0*24.0*3600.0;
	rmax = 82.0*24.0*3600.0*v;
	double gam = 1.0/sqrt(1.0 - v*v/speed_of_light2);
	double n0 = (9E6)/(82*82/49);
	Bfactor = sqrt(epsilonB*8*pi*(gam - 1.0)*n0*massProtonReal*speed_of_light2);
	concentration = 2*n0;
	fractionSize = 0.1;
	sigma = 0.02;*/
	//concentration = sqr(Bfactor)/(sigma*4*pi*massProtonReal*speed_of_light2);

	//Ho atm2020
	Bfactor = 6;
	rmax = 7E15;
	fractionSize = 0.5;
	concentration = 1.9;
	v = 0.15*speed_of_light;
	z = 0;

	//margutti at2018 16.5 days
	//rmax = 7.7*24*3600*0.1*speed_of_light;
	//css161010 357
	rmax = 357*24*3600*0.3*speed_of_light;
	fractionSize = 0.5;
	v = 0.1*speed_of_light;
	Bfactor = 0.048;
	epsilonB = 0.0012;
	concentration = Bfactor*Bfactor/(4*pi*massProtonReal*speed_of_light2*epsilonB);
	z = 0;


	printf("integrating fields\n");
	fprintf(logFile, "integrating Fields\n");
	fflush(logFile);
	double* totalInu = new double[Nnu];
	double finalSigma = sqr(Bfactor)/(4*pi*concentration*massProtonReal*speed_of_light2);
	

	evaluateVolumeAndLength(area3d, length3d, rmax, fractionSize);

	double**** Inu = new double***[Nnu];
	double**** Anu = new double***[Nnu];
	for(int l = 0; l < Nnu; ++l){
		Inu[l] = new double**[Nrho];
		Anu[l] = new double**[Nrho];
		for(int i = 0; i < Nrho; ++i){
			Inu[l][i] = new double*[Nphi];
			Anu[l][i] = new double*[Nphi];
			for(int j = 0; j < Nphi; ++j){
				Inu[l][i][j] = new double[Nz];
				Anu[l][i][j] = new double[Nz];
			}
		}
	}

	double*** image = new double**[Nrho];
	for(int i = 0; i < Nrho; ++i) {
		image[i] = new double*[Nphi];
		for(int j = 0; j < Nphi; ++j) {
			image[i][j] = new double[Nnu];
			for(int k = 0; k < Nnu; ++k) {
				image[i][j][k] = 0;
			}
		}
	}

	const int Nobs = 4;
	double nu1[Nobs];
	double observedInu[Nobs];
	double observedError[Nobs];


	//css161010 t = 357
	/*nu1[0] = 0.33E9;
	observedInu[0] = 0.357;
	observedError[0] = 0.09;
	nu1[1] = 0.61E9;
	observedInu[1] = 0.79;
	observedError[1] = 0.09;
	nu1[2] = 1.5E9;
	observedInu[2] = 0.27;
	observedError[2] = 0.07;
	nu1[3] = 3.0E9;
	observedInu[3] = 0.17;
	observedError[3] = 0.03;
	nu1[4] = 6.05E9;
	observedInu[4] = 0.07;
	observedError[4] = 0.01;
	nu1[5] = 10.0E9;
	observedInu[5] = 0.032;
	observedError[5] = 0.008;
	Time = 357*24*3600;*/
	//css1610101 t = 98
	nu1[0] = 1.5E9;
	observedInu[0] = 1.5;
	observedError[0] = 0.1;
	nu1[1] = 3.0E9;
	observedInu[1] = 4.3;
	observedError[1] = 0.2;
	nu1[2] = 6.1E9;
	observedInu[2] = 6.1;
	observedError[2] = 0.3;
	nu1[3] = 9.87E9;
	observedInu[3] = 4.2;
	observedError[3] = 0.2;
	Time = 98 * 24 * 3600;
	//at2018 t = 15
	//nu1[0] = 35E9;
	//observedInu[0] = 8;
	//nu1[1] = 225E9;
	//observedInu[1] = 30.8;
	//nu1[2] = 233E9;
	//observedInu[2] = 28.6;
	//nu1[3] = 241E9;
	//observedInu[3] = 27.4;
	//Time = 16.5*24*3600;
	//at2018 t = 7.7
	//nu1[0] = 15.5E9;
	//observedInu[0] = 0.489;
	//nu1[1] = 214E9;
	//observedInu[1] = 36.425;
	//nu1[2] = 326E9;
	//observedInu[2] = 32.705;
	//Time = 7.7*24*3600;


	double**** Inu1 = new double***[Nobs];
	double**** Anu1 = new double***[Nobs];
	for(int l = 0; l < Nobs; ++l){
		Inu1[l] = new double**[Nrho];
		Anu1[l] = new double**[Nrho];
		for(int i = 0; i < Nrho; ++i){
			Inu1[l][i] = new double*[Nphi];
			Anu1[l][i] = new double*[Nphi];
			for(int j = 0; j < Nphi; ++j){
				Inu1[l][i][j] = new double[Nz];
				Anu1[l][i][j] = new double[Nz];
			}
		}
	}

	//Bfactor = 0.069;
	Bfactor = 0.6;
	//rmax = 3E17;
	rmax = 1.3E17;
	epsilonB = 0.0005;
	fractionSize = 0.0001;

	double vector[4];
	vector[0] = Bfactor/maxB;
	vector[1] = rmax/maxR;
	vector[2] = fractionSize;
	vector[3] = epsilonB;
	double dopplerBeta = 0.75 * 0.5 * speed_of_light;

	printf("Bfactor = %g, n = %g\n fraction = %g rmax = %g sigma = %g\n", Bfactor, concentration, fractionSize, rmax, finalSigma);
	fprintf(logFile, "Bfactor = %g, n = %g fraction = %g rmax = %g v/c = %g sigma = %g\n", Bfactor, concentration, fractionSize, rmax, v / speed_of_light, finalSigma);
	fflush(logFile);

	bool optPar[4] = { true, true, true, true };

	//optimizeParametersGeneral(vector, optPar, nu1, observedInu, observedError, weightedEe, weightedFe, Np, Nobs, Ndist, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inu1, Anu1, area3d, length3d, dopplerBeta, logFile);
	stochasticGradientOptimization(vector, optPar, nu1, observedInu, observedError, weightedEe, weightedFe, Np, Nobs, Ndist, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inu1, Anu1, area3d, length3d, dopplerBeta, logFile);
	Bfactor = vector[0]*maxB;
	rmax = vector[1]*maxR;
	fractionSize = vector[2];
	epsilonB = vector[3];
	v = rmax/Time;
	concentration = Bfactor*Bfactor/(4*pi*massProtonReal*speed_of_light2*epsilonB);

	//
	//Bfactor = 1.0;
	//rmax = 3E17;
	//rmax = 1.0E17;
	//epsilonB = 0.005;
	//fractionSize = 0.0001;
	//concentration = Bfactor * Bfactor / (4 * pi * massProtonReal * speed_of_light2 * epsilonB);
	//concentration = 100;
	//double tempnu = pow(Bfactor, 0.3/0.155 - 1.25) * pow(concentration, 0.15/0.155 - 0.65);
	//double tempf = pow(Bfactor, 1.25) * pow(concentration, 0.65) * pow(fractionSize, 0.5);

	//Bfactor = 1.0;
	//concentration = pow(tempnu / pow(Bfactor, 0.3 / 0.155 - 1.25), 1.0 / (0.15 / 0.155 - 0.65));
	//fractionSize = pow(tempf / (pow(Bfactor, 1.25) * pow(concentration, 0.65)), 2.0);

	//concentration = sqrt(temp / Bfactor * Bfactor * Bfactor * Bfactor * Bfactor);
	//

	//optimizeParametersBandR(vector, nu1, observedInu, weightedEe, weightedFe, Np, Nobs, Ndist, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inu1, Anu1, area3d, length3d, fraction, epsilonB, logFile);
	//Bfactor = vector[0]*maxB;
	//rmax = vector[1]*maxR;
	//v = rmax/Time;
	//concentration = Bfactor*Bfactor/(4*pi*massProtonReal*speed_of_light2*epsilonB);


	//optimizeParameterB(Bfactor, nu1, observedInu, weightedEe, weightedFe, Np, Nobs, Ndist, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inu1, Anu1, area3d, length3d, rmax, fraction, epsilonB, logFile);
	//concentration = Bfactor*Bfactor/(4*pi*massProtonReal*speed_of_light2*epsilonB);

	for(int l = 0; l < Nobs; ++l){
		for(int i = 0; i < Nrho; ++i){
			for(int j = 0; j < Nphi; ++j){
				delete[] Inu1[l][i][j];
				delete[] Anu1[l][i][j];
			}
			delete[] Inu1[l][i];
			delete[] Anu1[l][i];
		}
		delete[] Inu1[l];
		delete[] Anu1[l];
	}
	delete[] Inu1;
	delete[] Anu1;

	evaluateAllEmissivityAndAbsorption(nu, Inu, Anu, Nnu, weightedEe, weightedFe, Np, Ndist, B3d, sintheta3d, thetaIndex3d, concentrations3d, concentration, Bfactor, 1.0, 0.75*v/speed_of_light);

	FILE* Ifile = fopen("Inu.dat","w");
	FILE* Afile = fopen("Anu.dat","w");
	FILE* Kfile = fopen("Knu.dat","w");
	double nuc = nu[10];
	for(int i = 0; i < Nnu; ++i){
		fprintf(Ifile,"%g %g\n", nu[i], Inu[i][5][5][0]);
		fprintf(Afile,"%g %g\n", nu[i], Anu[i][5][5][0]);
		double K = evaluateMcDonaldIntegral(nu[i]/nuc);
		fprintf(Kfile,"%g %g\n", nu[i], K);
	}
	fclose(Ifile);
	fclose(Afile);
	fclose(Kfile);



		
	if(geometry == SPHERICAL){
		evaluateImageSpherical(image, nu, Inu, Anu, rmax, Nnu, 1.0, fractionSize);
	} else {
		evaluateImageFlat(image, nu, Inu, Anu, rmax, Nnu, 1.0, fractionSize);
	}
		
	//evaluateSpectrum(nu, tempTotalInu[l], tempInu, tempAnu, area3d, length3d, Nnu, rfactor);
	if(geometry == SPHERICAL){
		evaluateSpectrumSpherical(nu, totalInu, Inu, Anu, rmax, Nnu, 1.0, fractionSize);
	} else {
		evaluateSpectrumFlat(nu, totalInu, Inu, Anu, rmax, Nnu, 1.0, fractionSize);
	}

	//evaluate rontgen
	const int rontgenNnu = 10;
	double**** rontgenInu = new double*** [rontgenNnu];
	double**** rontgenAnu = new double*** [rontgenNnu];
	for (int l = 0; l < Nnu; ++l) {
		rontgenInu[l] = new double** [Nrho];
		rontgenAnu[l] = new double** [Nrho];
		for (int i = 0; i < Nrho; ++i) {
			rontgenInu[l][i] = new double* [Nphi];
			rontgenAnu[l][i] = new double* [Nphi];
			for (int j = 0; j < Nphi; ++j) {
				rontgenInu[l][i][j] = new double[Nz];
				rontgenAnu[l][i][j] = new double[Nz];
			}
		}
	}
	double* rontgenTotalInu = new double[rontgenNnu];
	double* rontgenNu = new double[rontgenNnu];
	rontgenNu[0] = 0.3*1000 * 1.6E-12 / hplank;
	double deltaNu = (10 - 0.3) * 1000 * 1.6E-12 / hplank/(rontgenNnu);
	for (int i = 1; i < rontgenNnu; ++i) {
		rontgenNu[i] = rontgenNu[i - 1] + deltaNu;
	}
	evaluateAllEmissivityAndAbsorption(rontgenNu, rontgenInu, rontgenAnu, rontgenNnu, weightedEe, weightedFe, Np, Ndist, B3d, sintheta3d, thetaIndex3d, concentrations3d, concentration, Bfactor, 1.0, 0.75 * v / speed_of_light);
	if (geometry == SPHERICAL) {
		//evaluateSpectrumSpherical(rontgenNu, rontgenTotalInu, rontgenInu, rontgenAnu, rmax, rontgenNnu, 1.0, fractionSize);
		//evaluateSpectrumSpherical(rontgenNu, rontgenTotalInu, rontgenInu, rontgenAnu, rmax, rontgenNnu, 1.0, 0.001);
	}
	else {
		//evaluateSpectrumFlat(rontgenNu, rontgenTotalInu, rontgenInu, rontgenAnu, rmax, rontgenNnu, 1.0, fractionSize);
		//evaluateSpectrumFlat(rontgenNu, rontgenTotalInu, rontgenInu, rontgenAnu, rmax, rontgenNnu, 1.0, 0.001);
	}
	double rontgenFlux = 0;
	for (int i = 0; i < rontgenNnu; ++i) {
		rontgenFlux += rontgenTotalInu[i] * deltaNu / 1E26;
	}
	printf("rontgen flux = %g\n", rontgenFlux);
	fprintf(logFile, "rontgen flux = %g\n", rontgenFlux);
	fflush(logFile);

	//evaluateSpectrumFlatSimple(nu, tempTotalInu[l], Inuflat, Anuflat, Nnu, r, fractionSize);
	

	FILE* imageFile = fopen("image.dat","w");
	for(int k = 0; k < Nnu; ++k) {
		for(int i = 0; i < Nrho; ++i) {
			for(int j = 0; j < Nphi; ++j) {
				fprintf(imageFile, "%g", image[i][j][k]);
				if(j < Nphi - 1) {
					fprintf(imageFile, " ");
				}
			}
			fprintf(imageFile,"\n");
		}
	}
	fclose(imageFile);
	for(int i = 0; i < Nrho; ++i) {
		for(int j = 0; j < Nphi; ++j) {
			delete[] image[i][j];
		}
		delete[] image[i];
	}
	delete[] image;

	

	//////////
	printf("outputing\n");
	fprintf(logFile, "ooutputing\n");
	fflush(logFile);
	//char* number = new char[100];
	//itoa(0, number, 10);
	//delete[] number;
	//std::string fileNumber = std::string(number);
	FILE* output = fopen(outputfileName.c_str(), "w");
	for (int i = 0; i < Nnu; ++i) {	
		fprintf(output, "%g %g\n", nu[i]/((1 + z)*1E9), totalInu[i]/(1 + z));	
	}

	fclose(output);
	


	FILE* outputParam = fopen("parameters.dat","w");
	fprintf(outputParam, "B %g\n", Bfactor);
	fprintf(outputParam, "n %g\n", concentration);
	fprintf(outputParam, "f %g\n", fractionSize);
	fprintf(outputParam, "R %g\n", rmax);
	fprintf(outputParam, "v/c %g\n", v/speed_of_light);
	fclose(outputParam);

	////////////////////

	printf("deleting arrays\n");
	fprintf(logFile, "deleting arrays\n");
	fflush(logFile);
	delete[] totalInu;

	for(int j = 0; j < Ndist; ++j){
		delete[] Fe[j];
		delete[] Ee[j];
	}
	delete[] Fe;
	delete[] Ee;

	for(int i = 0; i < Nrho; ++i){
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				delete[] weightedFe[i][j][k];
			}
			delete[] area3d[i][j];
			delete[] length3d[i][j];
			delete[] B3d[i][j];
			delete[] Bx3d[i][j];
			delete[] By3d[i][j];
			delete[] Bz3d[i][j];
			delete[] concentrations3d[i][j];
			delete[] sintheta3d[i][j];
			delete[] thetaIndex3d[i][j];
			delete[] weights[i][j];
			delete[] weightedFe[i][j];
		}
		delete[] area3d[i];
		delete[] length3d[i];
		delete[] B3d[i];
		delete[] Bx3d[i];
		delete[] By3d[i];
		delete[] Bz3d[i];
		delete[] concentrations3d[i];
		delete[] sintheta3d[i];
		delete[] thetaIndex3d[i];
		delete[] weights[i];
		delete[] weightedFe[i];
	}
	delete[] area3d;
	delete[] length3d;
	delete[] B3d;
	delete[] Bx3d;
	delete[] By3d;
	delete[] Bz3d;
	delete[] concentrations3d;
	delete[] sintheta3d;
	delete[] thetaIndex3d;
	delete[] weights;
	delete[] weightedFe;
	delete[] weightedEe;

	for(int i = 0; i < Ntheta; ++i){
		for(int j = 0; j < Nphi; ++j){
			delete[] Inu[i][j];
			delete[] Anu[i][j];
		}
		delete[] Inu[i];
		delete[] Anu[i];
		//delete[] concentrations[i];
		//delete[] B[i];
		//delete[] sintheta[i];
		//delete[] thetaIndex[i];
	}
	delete[] Inu;
	delete[] Anu;
	//delete[] concentrations;
	//delete[] B;
	//delete[] sintheta;
	//delete[] thetaIndex;

	for(int i = 0; i < Nrho; ++i){
		delete[] area[i];
		delete[] length[i];
	}
	delete[] area;
	delete[] length;

	delete[] Rho;
	delete[] Phi;
	
	fclose(logFile);

	return 0;
}

