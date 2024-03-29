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
#include "optimization.h"
#include "spectrum.h"

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

void evaluateOrientationParameters3d(double*** B, double*** sintheta, double*** psi, int*** thetaIndex, double*** Bx, double*** By, double*** Bz, int Nd){
	for(int i = 0; i < Nrho; ++i){
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				//double r = rmax;
				psi[i][j][k] = atan2(By[i][j][k], Bx[i][j][k]);
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
				printf("thetaa index[%d][%d][%d] = %d\n", i, j, k, thetaIndex[i][j][k]);
				//for debug
				//thetaIndex[i][j][k] = 9;
			}
		}
	}
}

void evaluateOrientationParameters3dflat(double*** B, double*** sintheta, double*** psi, int*** thetaIndex, double*** Bx, double*** By, double*** Bz, int Nd){
	for(int i = 0; i < Nrho; ++i){
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				//double r = rmax;
				psi[i][j][k] = atan2(By[i][j][k], Bx[i][j][k]);
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
	//todo - also depends on p
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
	//todo - also depends on p
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

	double thetaObserv = 0;
	double cosThetaObserv = cos(thetaObserv);
	double sinThetaObserv = sin(thetaObserv);

	const int Np= 200;
	const int Nnu = 100;
	const int Nnu1 = 5;
	const int Ndist = 10;

	double** Fe;
	double** dFe;
	double** Ee;

	double*** Inu;
	double*** Anu;
	double* nu;

	double*** Inu1;
	double*** Anu1;
	double* nu1;


	double* Rho;
	double* Phi;

	int Nnum = 4;

	double***** Inumonth;
	double***** Anumonth;
	double** Fmonth;
	double** Numonth;
	double***** NuDoppler;

	double*** B3d;
	double*** Bx3d;
	double*** By3d;
	double*** Bz3d;
	double*** concentrations3d;
	int*** thetaIndex3d;
	double*** sintheta3d;
	double*** psi3d;

	

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

	double dr = rmax/Nrho;
	for(int i = 0; i < Nrho; ++i){
		Rho[i] = dr/2 + i*dr;
		//area[i] = new double[Nphi];
		//length[i] = new double[Nphi];
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

	B3d = new double**[Nrho];
	Bx3d = new double**[Nrho];
	By3d = new double**[Nrho];
	Bz3d = new double**[Nrho];
	concentrations3d = new double**[Nrho];
	sintheta3d = new double**[Nrho];
	psi3d = new double** [Nrho];
	thetaIndex3d = new int**[Nrho];
	for(int i = 0; i < Nrho; ++i){
		B3d[i] = new double*[Nphi];
		Bx3d[i] = new double*[Nphi];
		By3d[i] = new double*[Nphi];
		Bz3d[i] = new double*[Nphi];
		concentrations3d[i] = new double*[Nphi];
		sintheta3d[i] = new double*[Nphi];
		psi3d[i] = new double* [Nphi];
		thetaIndex3d[i] = new int*[Nphi];
		for(int j = 0; j < Nphi; ++j){
			B3d[i][j] = new double[Nz];
			Bx3d[i][j] = new double[Nz];
			By3d[i][j] = new double[Nz];
			Bz3d[i][j] = new double[Nz];
			concentrations3d[i][j] = new double[Nz];
			sintheta3d[i][j] = new double[Nz];
			psi3d[i][j] = new double[Nz];
			thetaIndex3d[i][j] = new int[Nz];
			for(int k = 0; k < Nz; ++k){
				double r = sqrt(1.0*(i+0.5)*(i+0.5) + (k + 0.5 - Nz/2.0)*(k + 0.5 - Nz/2.0));
				double rmax = Nrho;
				B3d[i][j][k] = 1.0;
				Bx3d[i][j][k] = 1.0;
				By3d[i][j][k] = 0.0;
				Bz3d[i][j][k] = 0.0;
				if(geometry == SPHERICAL) {
					double tempBeta = 0.75;
					concentrations3d[i][j][k] = 1.0*pow(1 - tempBeta, 7.0/4.0)/pow(1 - tempBeta*r/rmax, 7.0/4.0);
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
		evaluateOrientationParameters3d(B3d, sintheta3d, psi3d, thetaIndex3d, Bx3d, By3d, Bz3d, Ndist);
	} else {
		evaluateOrientationParameters3dflat(B3d, sintheta3d, psi3d, thetaIndex3d, Bx3d, By3d, Bz3d, Ndist);
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
	int powerlawStart[Ndist] = {155,155,155,155,155,155,145,145,145,145};


	Fe = new double*[Ndist];
	dFe = new double*[Ndist];
	Ee = new double*[Ndist];
	for(int j = 0; j < Ndist; ++j){
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
			for (int i = 1; i < Npowerlaw; ++i) {

				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				//todo massRelationSqrt?
				double gamma = u*realMassRelationSqrt/massRelationSqrt + 1;
				Ee[j][i] = gamma*massElectron*speed_of_light2;
				Fe[j][i] = Fe[j][i] / (massElectron*speed_of_light2);
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);

				/*double minGamma = 2;
				double power = 3;

				if(gamma >= minGamma){
					Fe[j][i] = 1.0/pow(Ee[j][i],power);
				} else {
					Fe[j][i] = 0;
				}
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);*/
			}
			Ee[j][Np-1] = 1E8;
			double factor = pow(Ee[j][Np-1]/Ee[j][Npowerlaw-1], 1.0/(Np-Npowerlaw));
			double gammae = -log(Fe[j][Npowerlaw-10]/Fe[j][Npowerlaw-1])/log(Ee[j][Npowerlaw-10]/Ee[j][Npowerlaw-1]);
			for(int i = Npowerlaw; i < Np; ++i){
				Ee[j][i] = Ee[j][i-1]*factor;
				Fe[j][i] = Fe[j][Npowerlaw-1]*pow(Ee[j][Npowerlaw-1]/Ee[j][i], gammae);
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
			}
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
				if(i > 137){
				//if((i > 137)&&(u < 500)){
					Fe[j][i] = Fe[j][137]*pow(Ee[j][i]/Ee[j][137], -3);
				}
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
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

	double factor = pow(maxE/minE, 1.0/(Np - 1));
	weightedEe[0] = minE;
	for(int i = 1; i < Np; ++i){
		weightedEe[i] = weightedEe[i-1]*factor;
	}
	weightedEe[Np-1] = maxE;

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

	int Nnumonth[Nmonth];
	//SN2009bb
	/*Nnumonth[0] = 4;
	Nnumonth[1] = 3;
	Nnumonth[2] = 4;
	Nnumonth[3] = 5;
	*/
	//CSS161010
	Nnumonth[0] = 4;
	Nnumonth[1] = 5;
	Nnumonth[2] = 6;

	nu = new double[Nnu];
	Inu = new double**[Ntheta];
	Anu = new double**[Ntheta];
	nu1 = new double[Nnu1];
	Inu1 = new double**[Ntheta];
	Anu1 = new double**[Ntheta];
	for(int i = 0; i < Ntheta; ++i){
		Inu[i] = new double*[Nphi];
		Anu[i] = new double*[Nphi];
		Inu1[i] = new double*[Nphi];
		Anu1[i] = new double*[Nphi];
		for(int j = 0; j < Nphi; ++j){
			Inu[i][j] = new double[Nnu];
			Anu[i][j] = new double[Nnu];
			Inu1[i][j] = new double[Nnu1];
			Anu1[i][j] = new double[Nnu1];
		}
	}

	double* Inuflat = new double[Nnu];
	double* Anuflat = new double[Nnu];

	Numonth = new double*[Nmonth];
	NuDoppler = new double****[Nmonth];
	Fmonth = new double*[Nmonth];
	double** ErrorMonth = new double*[Nmonth];
	Inumonth = new double****[Nmonth];
	Anumonth = new double****[Nmonth];
	for(int m = 0; m < Nmonth; ++m){
		Numonth[m] = new double[Nnumonth[m]];
		NuDoppler[m] = new double***[Nnumonth[m]];
		Fmonth[m] = new double[Nnumonth[m]];
		ErrorMonth[m] = new double[Nnumonth[m]];
		Inumonth[m] = new double***[Nnumonth[m]];
		Anumonth[m] = new double***[Nnumonth[m]];
		for(int l = 0; l < Nnumonth[m]; ++l){
			Numonth[m][l] = 0;
			NuDoppler[m][l] = new double**[Nrho];
			Fmonth[m][l]= 0;
			ErrorMonth[m][l] = 1.0;
			Inumonth[m][l] = new double**[Nrho];
			Anumonth[m][l] = new double**[Nrho];
			for(int i = 0; i < Nrho; ++i){
				NuDoppler[m][l][i] = new double*[Nphi];
				Inumonth[m][l][i] = new double*[Nphi];
				Anumonth[m][l][i] = new double*[Nphi];
				for(int j = 0; j < Nphi; ++j){
					NuDoppler[m][l][i][j] = new double[Nz];
					Inumonth[m][l][i][j] = new double[Nz];
					Anumonth[m][l][i][j] = new double[Nz];
					for(int k = 0; k < Nz; ++k){
						NuDoppler[m][l][i][j][k] = 0;
						Inumonth[m][l][i][j][k] = 0;
						Anumonth[m][l][i][j][k] = 0;
					}
				}
			}
		}
	}

	//SN2009bb
	/*for (int i = 0; i < Nnumonth[0]; ++i) {
		Numonth[0][i] = aprx[i]*1E9;
		Fmonth[0][i] = apry[i];
		ErrorMonth[0][i] = aprError[i];
	}

	for(int i = 0; i < Nnumonth[1]; ++i){
		Numonth[1][i] = mayx[i]*1E9;
		Fmonth[1][i] = mayy[i];
		ErrorMonth[1][i] = mayError[i];
	}

	for(int i = 0; i < Nnumonth[2]; ++i){
		Numonth[2][i] = junx[i]*1E9;
		Fmonth[2][i] = juny[i];
		ErrorMonth[2][i] = junError[i];
	}

	for(int i = 0; i < Nnumonth[3]; ++i){
		Numonth[3][i] = augx[i]*1E9;
		Fmonth[3][i] = augy[i];
		ErrorMonth[3][i] = augError[i];
	}

	if(Nmonth > 4){
		for(int i = 0; i < Nnumonth[4]; ++i){
			Numonth[4][i] = octx[i]*1E9;
			Fmonth[4][i] = octy[i];
			//ErrorMonth[4][i] = octError[i];
		}

		for(int i = 0; i < Nnumonth[5]; ++i){
			Numonth[5][i] = decx[i]*1E9;
			Fmonth[5][i] = decy[i];
			//ErrorMonth[5][i] = decError[i];
		}
	}*/

	//CSS161010
	for (int i = 0; i < Nnumonth[0]; ++i) {
		Numonth[0][i] = cssx1[i] * 1E9;
		Fmonth[0][i] = cssy1[i];
		ErrorMonth[0][i] = cssError1[i];
	}

	for (int i = 0; i < Nnumonth[1]; ++i) {
		Numonth[1][i] = cssx2[i] * 1E9;
		Fmonth[1][i] = cssy2[i];
		ErrorMonth[1][i] = cssError2[i];
	}

	for (int i = 0; i < Nnumonth[2]; ++i) {
		Numonth[2][i] = cssx3[i] * 1E9;
		Fmonth[2][i] = cssy3[i];
		ErrorMonth[2][i] = cssError3[i];
	}


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
	numax = 1000*1.6*1E-12/hplank;
	createNu(nu, Nnu, 0.001*1E9, numax);
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
	double fractionSize = 0.5;
	double V0 = speed_of_light;
	double v = 0.4*speed_of_light;
	double r0 = 0.1*maxR0;
	double a = 2.0; //B power r ^(a-1)
	double b = 3.0; //N power r ^ (b-1);
	double fpower = 1.0;//f-1 power r
	rmax = 3.4E16;
	////////////////////
	printf("initial optimizing parameters\n");
	fprintf(logFile, "initial optimizing parameters\n");
	fflush(logFile);
	printf("optimizing parameters\n");
	fprintf(logFile, "optimizing parameters\n");
	fflush(logFile);
	Bfactor = 0.5*maxB;
	concentration = 0.5*maxN;
	fractionSize = 0.5*maxFraction;
	rmax = 4.0E16;
	v = 0.67*speed_of_light;
	sigma = 0.02;
	//concentration = sqr(Bfactor)/(sigma*4*pi*massProtonReal*speed_of_light2);
	const int Nopt = 8;
	bool optPar[Nopt] = {true, true, true, true, true, false, false, false};
	double vector[Nopt];
	vector[0] = Bfactor/maxB;
	vector[1] = concentration/maxN;
	vector[2] = fractionSize/maxFraction;
	vector[3] = v/maxV;
	vector[4] = r0/maxR0;
	vector[5] = a/maxBpower;
	vector[6] = b/maxNpower;
	vector[7] = fpower/maxFpower;
	

	double timeMoments[Nmonth];
	for(int i = 0; i < Nmonth; ++i){
		timeMoments[i] = times[i];
	}

	double error = evaluateOptimizationFunction5(vector, timeMoments, Numonth, Fmonth, ErrorMonth, weightedEe, weightedFe, Np, Nnumonth, Ndist, Nmonth, B3d, sintheta3d, psi3d, thetaIndex3d, concentrations3d, NuDoppler, Inumonth, Anumonth);
	printf("error = %lf\n", error);
	fprintf(logFile, "error = %lf\n", error);
	const int Nbp = 8;
	const int Nnp = 11;
	const int Nfp = 6;
	const int Nvp = 8;
	const int Nr0 = 4;
	int Na = 3;
	double* apoints;
	if (parker) {
		Na = 1;
		apoints = new double[Na];
		apoints[0] = 2.0;
	} else {
		Na = 3;
		apoints = new double[Na];
		apoints[0] = 1.0;
		apoints[1] = 2.0;
		apoints[2] = 3.0;
	}
	int Nb1 = 3;
	double* bpoints;
	if (parker) {
		Nb1 = 1;
		bpoints = new double[Nb1];
		bpoints[0] = 3.0;
	} else {
		Nb1 = 2;
		bpoints = new double[Nb1];
		bpoints[0] = 1.0;
		bpoints[1] = 3.0;
	}
	const int Nf1 = 2;
	double Bpoints[Nbp] = {0.1, 0.2, 0.5, 1, 2, 5, 10, 20};
	double npoints[Nnp] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000};
	double fpoints[Nfp] = {0.03, 0.06, 0.1,0.3,0.4,0.5};
	double vpoints[Nvp] = {0.3 * speed_of_light, 0.4*speed_of_light, 0.5*speed_of_light, 0.6*speed_of_light, 0.65*speed_of_light, 0.7*speed_of_light, 0.75*speed_of_light, 0.8*speed_of_light};
	double rpoints[Nr0] = {0.1*maxR0, 0.2*maxR0, 0.4*maxR0, 0.5*maxR0};
	double fppoints[Nf1] = {1,2};
	if(initialGridSearch){
		for(int i = 0; i < Nbp; ++i){
			double tempBfactor = Bpoints[i];
			for(int j = 0; j < Nnp; ++j){
				tempConcentration = npoints[j];
				for(int k = 0; k < Nfp; ++k){
					double tempFractionSize = fpoints[k];
					for(int l = 0; l < Nvp; ++l){
						double tempV = vpoints[l];
						for(int m = 0; m < Nr0; ++m){
							double tempR0 = rpoints[m];
							for(int ii = 0; ii < Na; ++ii){
								double tempA = apoints[ii];
								for(int jj = 0; jj < Nb1; ++jj){
									double tempB = bpoints[jj];
									for(int kk = 0; kk < Nf1; ++kk){
										double tempFpower = fppoints[kk];
										vector[0] = tempBfactor/maxB;
										vector[1] = tempConcentration/maxN;
										vector[2] = tempFractionSize/maxFraction;
										vector[3] = tempV/maxV;
										vector[4] = tempR0/maxR0;
										vector[5] = tempA/maxBpower;
										vector[6] = tempB/maxNpower;
										vector[7] = tempFpower/maxFpower;
										double tempError = evaluateOptimizationFunction5(vector, timeMoments, Numonth, Fmonth, ErrorMonth, weightedEe, weightedFe, Np, Nnumonth, Ndist, Nmonth, B3d, sintheta3d, psi3d, thetaIndex3d, concentrations3d, NuDoppler, Inumonth, Anumonth);
										printf("tempError = %lf\n", tempError);
										//fprintf(logFile, "tempError = %lf\n", tempError);
										if(tempError < error){
											error = tempError;
											Bfactor = vector[0]*maxB;
											concentration = vector[1]*maxN;
											fractionSize = vector[2]*maxFraction;
											v = vector[3]*maxV;
											r0 = vector[4]*maxR0;
											a = vector[5]*maxBpower;
											b = vector[6]*maxNpower;
											fpower = vector[7]*maxFpower;
											fprintf(logFile, "tempError = %lf, Bfactor = %lf, concentration = %lf, fraction = %lf, v/c = %lf r0 = %g a = %lf b = %lf fpower = %lf\n", error, Bfactor, concentration, fractionSize, v/speed_of_light, r0, a, b, fpower);
											printf("tempError = %lf, Bfactor = %lf, concentration = %lf, fraction = %lf, v/c = %lf\n r0 = %g a = %lf b = %lf fpower = %lf", error, Bfactor, concentration, fractionSize, v/speed_of_light, r0, a, b, fpower);
											fflush(logFile);
											FILE* initialFile = fopen("initialParameters.dat","w");
											fprintf(initialFile,"%lf %lf %lf %lf %lf %lf %lf %lf", Bfactor, concentration, fractionSize, v/speed_of_light, r0, a, b, fpower);
											fclose(initialFile);
										}
									}
								}
							}
						}
					}
				}
			}
		} 
	} else {
		FILE* initialFile = fopen("initialParameters.dat","r");
		fscanf(initialFile, "%lf %lf %lf %lf %lf %lf %lf %lf", &Bfactor, &concentration, &fractionSize, &v, &r0, &a, &b, &fpower);
		v = v*speed_of_light;
		fclose(initialFile);
		/*Bfactor = 0.362992;
		concentration = 501.865;
		fractionSize = 0.112738;
		r0 = 5.10995e+15;
		v = 0.364802*speed_of_light;
		a = 2.0;
		b = 2.0;
		fpower = 2.56063;*/
	}
	delete[] bpoints;
	delete[] apoints;

	vector[0] = Bfactor/maxB;
	vector[1] = concentration/maxN;
	vector[2] = fractionSize/maxFraction;
	vector[3] = v/maxV;
	vector[4] = r0/maxR0;
	vector[5] = a/maxBpower;
	vector[6] = b/maxNpower;
	vector[7] = fpower/maxFpower;

	error = evaluateOptimizationFunction5(vector, timeMoments, Numonth, Fmonth, ErrorMonth, weightedEe, weightedFe, Np, Nnumonth, Ndist, Nmonth, B3d, sintheta3d, psi3d, thetaIndex3d, concentrations3d, NuDoppler, Inumonth, Anumonth);
	fprintf(logFile, "tempError = %lf, Bfactor = %lf, concentration = %lf, fraction = %lf, v/c = %lf, r0 = %lf a = %lf b = %lf\n", error, Bfactor, concentration, fractionSize, v/speed_of_light, r0, a, b);
	printf("tempError = %lf, Bfactor = %lf, concentration = %lf, fraction = %lf, rmax = %lf, v/c = %lf, r0 = %lf a = %lf b = %lf\n", error, Bfactor, concentration, fractionSize, v/speed_of_light, r0, a, b);


	if(optimization){
		if(geometry == FLAT_SIMPLE){
			//optimizeParameters5simple(1.0, 2000, 3.4E16, V0, Bfactor, concentration, fractionSize, rmax, v, weightedEe, weightedFe[0][0][0], Np, Ndist, 1.0, 8, logFile);
		} else {
			//optimizeParametersGeneral(vector, optPar, timeMoments, Numonth, Fmonth, ErrorMonth, weightedEe, weightedFe, Np, Nnumonth, Ndist, Nmonth, B3d, sintheta3d, psi3d, thetaIndex3d, concentrations3d, NuDoppler, Inumonth, Anumonth, logFile);
			stochasticGradientOptimization(vector, optPar, timeMoments, Numonth, Fmonth, ErrorMonth, weightedEe, weightedFe, Np, Nnumonth, Ndist, Nmonth, B3d, sintheta3d, psi3d, thetaIndex3d, concentrations3d, NuDoppler, Inumonth, Anumonth, logFile);
			Bfactor = vector[0]*maxB;
			concentration = vector[1]*maxN;
			fractionSize = vector[2]*maxFraction;
			//rmax = vector[3]*maxR;
			v = vector[3]*maxV;
			r0 = vector[4]*maxR0;
			a = vector[5]*maxBpower;
			b = vector[6]*maxNpower;
			fpower = vector[7]*maxFpower;
		}
	}
	
	//error = evaluateOptimizationFunction5(Bfactor, concentration, fractionSize, rmax, v, Numonth, Fmonth, weightedEe, weightedFe, Np, Nnum, Ndist, Nmonth, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inumonth, Anumonth, area3d, length3d);

	///////////////////


	error = evaluateOptimizationFunction5(vector, timeMoments, Numonth, Fmonth, ErrorMonth, weightedEe, weightedFe, Np, Nnumonth, Ndist, Nmonth, B3d, sintheta3d, psi3d, thetaIndex3d, concentrations3d, NuDoppler, Inumonth, Anumonth);

	int optimizationNumber = 0;
	for (int i = 0; i < Nopt; ++i) {
		if (optPar[i]) {
			optimizationNumber++;
		}
	}
	int pointsNumber = 0;
	for (int i = 0; i < Nmonth; ++i) {
		pointsNumber += Nnumonth[i];
	}
	printf("integrating fields\n");
	fprintf(logFile, "integrating Fields\n");
	fflush(logFile);
	double* totalInu = new double[Nnu];
	double Rfiducial = 1E16;
	double gamma = 1.0 / sqrt(1.0 - 0.75 * 0.75 * v * v / speed_of_light2);
	double finalSigma = sqr(Bfactor)/(4*pi*gamma*concentration*massProtonReal*speed_of_light2);
	printf("Bfactor = %g, n = %g\n fraction = %g v/c = %g r0 = %g a = %g b = %g fpower = %g sigma = %g\n", Bfactor, concentration, fractionSize, v/speed_of_light, r0, a, b, fpower, finalSigma);
	printf("error = %g\n", error);
	fprintf(logFile, "Bfactor = %g, n = %g fraction = %g v/c = %g r0 = %g a = %g b = %g fpower = %g\n", Bfactor, concentration, fractionSize, v/speed_of_light, r0, a, b, fpower);
	fprintf(logFile, "sigma = %g\n", finalSigma);
	double r1 = r0 + v * timeMoments[0];
	fprintf(logFile, "R at t0 = %g\n",r1);
	fprintf(logFile, "B ~ 1/r^%lf,  N ~ 1/r^%lf\n", a - 1, b - 1);
	fprintf(logFile, "B at %g cm = %g G, N at %g cm = %g cm^-3\n", Rfiducial, Bfactor / pow(Rfiducial / r1, a - 1),Rfiducial, concentration / pow(Rfiducial / r1, b - 1));
	fprintf(logFile, "error = %g\n", error);
	fprintf(logFile, "number of parameters = %d, number of points = %d\n", optimizationNumber, pointsNumber);
	fprintf(logFile, "chi nu = %g\n", error/(pointsNumber - optimizationNumber));
	fprintf(logFile, "approximate mass loss = %g Msun/year\n", 3.14E7*4*pi*massProtonReal*r1*r1*concentration*1.0E8/(2.0E33));
	fflush(logFile);

	double** tempTotalInu = new double*[Nmonth];

	//evaluateVolumeAndLength(area3d, length3d, rmax, fractionSize);

	double**** tempInu = new double***[Nnu];
	double**** tempAnu = new double***[Nnu];
	double**** tempNuDoppler = new double***[Nnu];
	for(int l = 0; l < Nnu; ++l){
		tempInu[l] = new double**[Nrho];
		tempAnu[l] = new double**[Nrho];
		tempNuDoppler[l] = new double**[Nrho];
		for(int i = 0; i < Nrho; ++i){
			tempInu[l][i] = new double*[Nphi];
			tempAnu[l][i] = new double*[Nphi];
			tempNuDoppler[l][i] = new double*[Nphi];
			for(int j = 0; j < Nphi; ++j){
				tempInu[l][i][j] = new double[Nz];
				tempAnu[l][i][j] = new double[Nz];
				tempNuDoppler[l][i][j] = new double[Nz];
			}
		}
	}

	double*** image = new double**[Nrho];
	for(int i = 0; i < Nrho; ++i) {
		image[i] = new double*[Nphi];
		for(int j = 0; j < Nphi; ++j) {
			image[i][j] = new double[Nnum];
			for(int k = 0; k < Nnum; ++k) {
				image[i][j][k] = 0;
			}
		}
	}

	rmin = r0 + v*timeMoments[0];
	evaluateNuDoppler(tempNuDoppler, Nnu, nu, 0.75*v/speed_of_light); 
	for(int l = 0; l < Nmonth; ++l){
		printf("evaluate final spectrum month = %d\n", l);
		fprintf(logFile, "evaluate final spectrum month = %d\n", l);
		fflush(logFile);
		double r = r0 + v*timeMoments[l];
		double rfactor = r/rmin;
		//double locB = Bfactor/(rfactor*rfactor);
		//double locN = concentration/pow(rfactor, 3.6666);

		tempTotalInu[l] = new double[Nnu];

		//evaluateVolumeAndLength(area, length, r, Rho, Phi, fractionSize);
		//evaluateAllEmissivityAndAbsorption1(nu, Inu, Anu, Nnu, Ee, dFe, Np, Ndist, B, sintheta, thetaIndex, concentrations, locN, locB, fractionSize);
		//evaluateSpectrum(nu, tempTotalInu[l], Inu, Anu, area, length, Nnu, r, Rho, Phi);

		double dopplerBeta = 0.75 * v / speed_of_light;
		updateConcentartions(concentrations3d, v / speed_of_light);
		evaluateAllEmissivityAndAbsorption(tempNuDoppler, tempInu, tempAnu, Nnu, weightedEe, weightedFe, Np, Ndist, B3d, sintheta3d, psi3d, thetaIndex3d, concentrations3d, concentration, Bfactor, rfactor, a, b, dopplerBeta);
		double tempFraction = fractionSize/pow(rfactor, fpower - 1.0);
		if(l == 0) {
			if(geometry == SPHERICAL){
				evaluateImageSpherical(image, nu, tempNuDoppler, tempInu, tempAnu, r, Nnum, rfactor, tempFraction, dopplerBeta);
			} else {
				evaluateImageFlat(image, nu, tempNuDoppler, tempInu, tempAnu, r, Nnum, rfactor, tempFraction, dopplerBeta);
			}
		}
		//evaluateSpectrum(nu, tempTotalInu[l], tempInu, tempAnu, area3d, length3d, Nnu, rfactor);
		if(geometry == SPHERICAL){
			evaluateSpectrumSpherical(nu, tempNuDoppler, tempTotalInu[l], tempInu, tempAnu, r, Nnu, tempFraction, dopplerBeta);
		} else {
			evaluateSpectrumFlat(nu, tempNuDoppler, tempTotalInu[l], tempInu, tempAnu, r, Nnu, tempFraction, dopplerBeta);
		}

		//evaluateSpectrumFlatSimple(nu, tempTotalInu[l], Inuflat, Anuflat, Nnu, r, fractionSize);
	}

	FILE* imageFile = fopen("image.dat","w");
	for(int k = 0; k < Nnum; ++k) {
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

	printf("outputing\n");
	fprintf(logFile, "ooutputing\n");
	fflush(logFile);
	//char* number = new char[100];
	//itoa(0, number, 10);
	//delete[] number;
	//std::string fileNumber = std::string(number);
	FILE* output = fopen(outputfileName.c_str(), "w");
	for (int i = 0; i < Nnu; ++i) {
		for (int j = 0; j < Nmonth; ++j) {
			if (tempTotalInu[j][i] != tempTotalInu[j][i]) {
				tempTotalInu[j][i] = 0;
			}
		}
		fprintf(output, "%g ", nu[i] / 1E9);
		for (int j = 0; j < Nmonth; ++j) {
			fprintf(output, "%g ", tempTotalInu[j][i]);
		}
		fprintf(output, "\n");
	}

	fclose(output);

	for (int l = 0; l < Nmonth; ++l) {
		delete[] tempTotalInu[l];
	}
	delete[] tempTotalInu;

	FILE* errorFile = fopen("error.dat","w");
	FILE* Bp = fopen("Bpoints.dat","w");
	FILE* np = fopen("Npoints.dat","w");
	for(int i = 0; i < Nbp; ++i){
		fprintf(Bp, "%g\n", Bpoints[i]);
	}
	for(int i = 0; i < Nnp; ++i){
		fprintf(np, "%g\n", npoints[i]);
	}
	fclose(Bp);
	fclose(np);
	for(int i = 0; i < Nbp; ++i){
		Bfactor = Bpoints[i];
		for(int j = 0; j < Nnp; ++j){
			concentration = npoints[j];
			vector[0] = Bfactor/maxB;
			vector[1] = concentration/maxN;
			vector[2] = fractionSize/maxFraction;
			vector[3] = v/maxV;
			vector[4] = r0/maxR0;
			vector[5] = a/maxBpower;
			vector[6] = b/maxNpower;
			vector[7] = fpower/maxFpower;

			double tempError = evaluateOptimizationFunction5(vector, timeMoments, Numonth, Fmonth, ErrorMonth, weightedEe, weightedFe, Np, Nnumonth, Ndist, Nmonth, B3d, sintheta3d, psi3d, thetaIndex3d, concentrations3d, NuDoppler, Inumonth, Anumonth);
			fprintf(errorFile, "%g ", tempError);
		}
		fprintf(errorFile, "\n");
	}
	fclose(errorFile);

	//////////

	

	//Chevalier compare to SN1987A
	/*const int Nchevalier = 5;
	const int Ntchev = 100;
	double Nuchev = 0.843E9;
	double* tchev = new double[Ntchev];
	double* chevTotalInu = new double[Ntchev];
	tchev[0] = 24*3600;
	factor = pow(1000.0, 1.0/(Ntchev - 1));
	for(int i = 1; i < Ntchev; ++i){
		tchev[i] = tchev[i-1]*factor;
	}
	double dchev = 0.05*3.08E24;
	double rchev = 1.5E15;
	double Bchev = 0.03;
	double nchev = (6E-19)/(1.6E-24);
	nchev = 1;
	double fchev = 0.2;
	double vchev = 5100000000;

	for(int l = 0; l < Ntchev; ++l){
		double r = rchev + vchev*tchev[l];
		double locB = Bchev*rchev/r;
		//locB = Bchev;
		double locN = nchev*pow(tchev[l]/(365*24*3600),-3)*pow(vchev/1000000000,-9);
		locN - nchev*(rchev/r)*(rchev/r);
		//locN = nchev;

		for(int i = 0; i < Nrho; ++i){
			for(int j = 0; j < Nphi; ++j){
				for(int k = 0; k < Nz; ++k){
					evaluateEmissivityAndAbsorptionAtNuSimple(Nuchev, tempInu[0][i][j][k], tempAnu[0][i][j][k], weightedEe, weightedFe[i][j][k], Np, sintheta3d[i][j][k], locB, locN);
				}
			}
		}

		//evaluateAllEmissivityAndAbsorption(nu, tempInu, tempAnu, Nnu, weightedEe, weightedFe, Np, Ndist, B3d, sintheta3d, thetaIndex3d, concentrations3d, nchev[l], Bchev[l], 1.0);


		//evaluateSpectrum(nu, tempTotalInu[l], tempInu, tempAnu, area3d, length3d, Nnu, rfactor);
		if(geometry == SPHERICAL){
			evaluateSpectrumSphericalAtNu(Nuchev, chevTotalInu[l], tempInu[0], tempAnu[0], r, 1.0, fchev, dchev);
		} else {
			evaluateSpectrumFlatAtNu(Nuchev, chevTotalInu[l], tempInu[0], tempAnu[0], r, 1.0, fchev, dchev);
		}
	}

	FILE* chevOutput = fopen("chevalier.dat","w");

	for(int i = 0; i < Nnu; ++i){
		fprintf(chevOutput, "%g", nu[i]/1E9, chevTotalInu[0][i]);
		for(int l = 0; l < Nchev; ++l){
			fprintf(chevOutput, " %g", chevTotalInu[l][i]);
		}
		fprintf(chevOutput, "\n");
	}

	for(int i = 0; i < Ntchev; ++i){
		fprintf(chevOutput, "%g %g\n", tchev[i]/(24*3600), chevTotalInu[i]);
	}

	fclose(chevOutput);

	delete[] tchev;
	delete[] chevTotalInu;*/
	//

	FILE* outputParam = fopen("parameters.dat","w");
	fprintf(outputParam, "B = %g\n", Bfactor);
	fprintf(outputParam, "n = %g\n", concentration);
	fprintf(outputParam, "f = %g\n", fractionSize);
	fprintf(outputParam, "v/c = %g\n", v/speed_of_light);
	fprintf(outputParam, "r0 = %g\n", r0);
	fprintf(outputParam, "a = %g\n", a);
	fprintf(outputParam, "b = %g\n", b);
	fprintf(outputParam, "fpower = %g\n", fpower);
	fclose(outputParam);

	for(int l = 0; l < Nnu; ++l){
		for(int i = 0; i < Nrho; ++i){
			for(int j = 0; j < Nphi; ++j){
				delete[] tempInu[l][i][j];
				delete[] tempAnu[l][i][j];
			}
			delete[] tempInu[l][i];
			delete[] tempAnu[l][i];
		}
		delete[] tempInu[l];
		delete[] tempAnu[l];
	}
	delete[] tempInu;
	delete[] tempAnu;

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

	for(int m = 0; m < Nmonth; ++m){
		for(int l = 0; l < Nnumonth[m]; ++l){
			for(int i = 0; i < Nrho; ++i){
				for(int j = 0; j < Nphi; ++j){
					delete[] Inumonth[m][l][i][j];
					delete[] Anumonth[m][l][i][j];
				}
				delete[] Inumonth[m][l][i];
				delete[] Anumonth[m][l][i];
			}
			delete[] Inumonth[m][l];
			delete[] Anumonth[m][l];
		}
		delete[] Inumonth[m];
		delete[] Anumonth[m];
		delete[] Fmonth[m];
		delete[] Numonth[m];
	}
	delete[] Inumonth;
	delete[] Anumonth;
	delete[] Fmonth;
	delete[] Numonth;

	for(int i = 0; i < Nrho; ++i){
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				delete[] weightedFe[i][j][k];
			}
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
		delete[] Anu[i];;
	}
	delete[] Inu;
	delete[] Anu;

	delete[] Inuflat;
	delete[] Anuflat;

	delete[] Rho;
	delete[] Phi;
	
	fclose(logFile);

	return 0;
}

