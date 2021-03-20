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
				if((i + j + k) > 4){
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
}



int main()
{
	double** Bx;
	double** By;
	double** Bz;

	double** area;
	double** length;

	double** concentrations;

	double** B;
	double** sintheta;
	int** thetaIndex;

	double thetaObserv = 0;
	double cosThetaObserv = cos(thetaObserv);
	double sinThetaObserv = sin(thetaObserv);

	int Np= 200;
	int Nnu = 100;
	int Nnu1 = 5;
	int Ndist = 10;

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

	double*** B3d;
	double*** Bx3d;
	double*** By3d;
	double*** Bz3d;
	double*** concentrations3d;
	double*** area3d;
	double*** length3d;
	int*** thetaIndex3d;
	double*** sintheta3d;

	

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
					concentrations3d[i][j][k] = 1.0*sqr(rmax/r);
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
	double turbulenceFraction = 0.9;
	double turbNorm = evaluateTurbNorm(kmax, Nk, 1.0, turbulenceFraction);

	if(turbulence){
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
	}

	if(geometry == SPHERICAL){
		evaluateOrientationParameters3d(B3d, sintheta3d, thetaIndex3d, Bx3d, By3d, Bz3d, Ndist);
	} else {
		evaluateOrientationParameters3dflat(B3d, sintheta3d, thetaIndex3d, Bx3d, By3d, Bz3d, Ndist);
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
			for (int i = 1; i < Np; ++i) {
				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				double gamma = u*realMassRelationSqrt/massRelationSqrt + 1;
				Ee[j][i] = gamma*massElectron*speed_of_light2;

				double minGamma = 2;
				double power = 3;

				if(gamma >= minGamma){
					Fe[j][i] = 1.0/pow(Ee[j][i],power);
				} else {
					Fe[j][i] = 0;
				}
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

	printf("evaluating local emissivity\n");
	fprintf(logFile, "evaluating local emissivity\n");
	fflush(logFile);

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
	Fmonth = new double*[Nmonth];
	Inumonth = new double****[Nmonth];
	Anumonth = new double****[Nmonth];
	for(int m = 0; m < Nmonth; ++m){
		Numonth[m] = new double[Nnum];
		Fmonth[m] = new double[Nnum];
		Inumonth[m] = new double***[Nnum];
		Anumonth[m] = new double***[Nnum];
		for(int l = 0; l < Nnum; ++l){
			Numonth[m][l] = 0;
			Fmonth[m][l]= 0;
			Inumonth[m][l] = new double**[Nrho];
			Anumonth[m][l] = new double**[Nrho];
			for(int i = 0; i < Nrho; ++i){
				Inumonth[m][l][i] = new double*[Nphi];
				Anumonth[m][l][i] = new double*[Nphi];
				for(int j = 0; j < Nphi; ++j){
					Inumonth[m][l][i][j] = new double[Nz];
					Anumonth[m][l][i][j] = new double[Nz];
					for(int k = 0; k < Nz; ++k){
						Inumonth[m][l][i][j][k] = 0;
						Anumonth[m][l][i][j][k] = 0;
					}
				}
			}
		}
	}
	Numonth[0][0] = aprx[0]*1E9;
	Fmonth[0][0] = apry[0];
	Numonth[0][1] = aprx[1]*1E9;
	Fmonth[0][1] = apry[1];
	Numonth[0][2] = aprx[2]*1E9;
	Fmonth[0][2] = apry[2];
	Numonth[0][3] = aprx[3]*1E9;
	Fmonth[0][3] = apry[3];

	Numonth[1][0] = mayx[0]*1E9;
	Fmonth[1][0] = mayy[0];
	Numonth[1][1] = mayx[0]*1E9;
	Fmonth[1][1] = mayy[0];
	Numonth[1][2] = mayx[1]*1E9;
	Fmonth[1][2] = mayy[1];
	Numonth[1][3] = mayx[2]*1E9;
	Fmonth[1][3] = mayy[2];

	Numonth[2][0] = junx[0]*1E9;
	Fmonth[2][0] = juny[0];
	Numonth[2][1] = junx[1]*1E9;
	Fmonth[2][1] = juny[1];
	Numonth[2][2] = junx[2]*1E9;
	Fmonth[2][2] = juny[2];
	Numonth[2][3] = junx[3]*1E9;
	Fmonth[2][3] = juny[3];

	Numonth[3][0] = augx[1]*1E9;
	Fmonth[3][0] = augy[1];
	Numonth[3][1] = augx[2]*1E9;
	Fmonth[3][1] = augy[2];
	Numonth[3][2] = augx[3]*1E9;
	Fmonth[3][2] = augy[3];
	Numonth[3][3] = augx[4]*1E9;
	Fmonth[3][3] = augy[4];

	if(Nmonth > 4){
		Numonth[4][0] = octx[0]*1E9;
		Fmonth[4][0] = octy[0];
		Numonth[4][1] = octx[0]*1E9;
		Fmonth[4][1] = octy[0];
		Numonth[4][2] = octx[1]*1E9;
		Fmonth[4][2] = octy[1];
		Numonth[4][3] = octx[2]*1E9;
		Fmonth[4][3] = octy[2];

		Numonth[5][0] = decx[0]*1E9;
		Fmonth[5][0] = decy[0];
		Numonth[5][1] = decx[0]*1E9;
		Fmonth[5][1] = decy[0];
		Numonth[5][2] = decx[1]*1E9;
		Fmonth[5][2] = decy[1];
		Numonth[5][3] = decx[2]*1E9;
		Fmonth[5][3] = decy[2];
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
	createNu(nu, Nnu, 0.001*1E9, 10000*1E9);
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
	double v = 0.75*speed_of_light;
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
	double error = evaluateOptimizationFunction5(Bfactor, concentration, fractionSize, rmax, v, Numonth, Fmonth, Ee, dFe, Np, Nnum, Ndist, Nmonth, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inumonth, Anumonth, area3d, length3d);
	printf("error = %lf\n", error);
	fprintf(logFile, "error = %lf\n", error);
	const int Nbp = 10;
	const int Nnp = 12;
	const int Nfp = 7;
	const int Nvp = 5;
	const int Nrp = 5;
	double Bpoints[Nbp] = {0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2};
	double npoints[Nnp] = {0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200};
	double fpoints[Nfp] = {0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2};
	double vpoints[Nvp] = { 0.6*speed_of_light, 0.65*speed_of_light, 0.7*speed_of_light, 0.75*speed_of_light, 0.8*speed_of_light};
	double rpoints[Nrp] = {3.0E16, 3.4E16, 3.6E16, 3.8E16, 4.0E16};
	for(int i = 0; i < Nbp; ++i){
		double tempBfactor = Bpoints[i];
		for(int j = 0; j < Nnp; ++j){
			tempConcentration = npoints[j];
			for(int k = 0; k < Nfp; ++k){
				double tempFractionSize = fpoints[k];
				for(int l = 0; l < Nvp; ++l){
					double tempV = vpoints[l];
					for(int m = 0; m < Nrp; ++m){
						double tempRmax = rpoints[m];
						double tempError = evaluateOptimizationFunction5(tempBfactor, tempConcentration, tempFractionSize, tempRmax, tempV, Numonth, Fmonth, Ee, dFe, Np, Nnum, Ndist, Nmonth, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inumonth, Anumonth, area3d, length3d);
						//printf("tempError = %lf\n", tempError);
						//fprintf(logFile, "tempError = %lf\n", tempError);
						if(tempError < error){
							error = tempError;
							Bfactor = tempBfactor;
							concentration = tempConcentration;
							fractionSize = tempFractionSize;
							rmax = tempRmax;
							v = tempV;
							fprintf(logFile, "tempError = %lf, Bfactor = %lf, concentration = %lf, fraction = %lf, rmax = %lf, v = %lf\n", error, Bfactor, concentration, fractionSize, rmax, v);
							printf("tempError = %lf, Bfactor = %lf, concentration = %lf, fraction = %lf, rmax = %lf, v = %lf\n", error, Bfactor, concentration, fractionSize, rmax, v);
						}
					}
				}
			}
		}
	}
	/*Bfactor = 0.2;
	concentration = 200;
	fractionSize = 0.2;
	rmax = 3E16;
	v = 0.7*speed_of_light;*/

	error = evaluateOptimizationFunction5(Bfactor, concentration, fractionSize, rmax, v, Numonth, Fmonth, Ee, dFe, Np, Nnum, Ndist, Nmonth, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inumonth, Anumonth, area3d, length3d);
	fprintf(logFile, "tempError = %lf, Bfactor = %lf, concentration = %lf, fraction = %lf, rmax = %lf, v = %lf\n", error, Bfactor, concentration, fractionSize, rmax, v);
	printf("tempError = %lf, Bfactor = %lf, concentration = %lf, fraction = %lf, rmax = %lf, v = %lf\n", error, Bfactor, concentration, fractionSize, rmax, v);

	double vector[5];
	vector[0] = Bfactor/maxB;
	vector[1] = concentration/maxN;
	vector[2] = fractionSize/maxFraction;
	vector[3] = rmax/maxR;
	vector[4] = v/maxV;

	//optimizeParameters(Bfactor, concentration, fractionSize, nu1, rmax, Ee, dFe, Np, Nnu1, Nd, B, sintheta, thetaIndex, concentrations, Inu1, Anu1, area, length, Rho, Phi, logFile);
	//optimizeParameters4(1.0, 2000, 3.4E16, Bfactor, concentration, fractionSize, rmax, nu1, Ee, dFe, Np, Nnu1, Ndist, B, sintheta, thetaIndex, concentrations, Inu1, Anu1, area, length, Rho, Phi, logFile);
	//optimizeParameters4(Bfactor, concentration, fractionSize,rmax, nu1, Ee, dFe, Np, Nnu1, Nd, B, sintheta, thetaIndex, concentrations, Inu1, Anu1, area, length, Rho, Phi, logFile);
	if(geometry == FLAT_SIMPLE){
		optimizeParameters5simple(1.0, 2000, 3.4E16, V0, Bfactor, concentration, fractionSize, rmax, v, Ee, dFe, Np, Ndist, 1.0, 8, logFile);
	} else {
		optimizeParameters5(vector, Numonth, Fmonth, Ee, dFe, Np, Nnum, Ndist, Nmonth, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inumonth, Anumonth, area3d, length3d, logFile);
		Bfactor = vector[0]*maxB;
		concentration = vector[1]*maxN;
		fractionSize = vector[2]*maxFraction;
		rmax = vector[3]*maxR;
		v = vector[4]*maxV;
	}
	//optimizeParameters5sigma(sigma, 1.0, N0, 3.4E16,V0, Bfactor, concentration, fractionSize, rmax, v, Numonth, Fmonth, Ee, dFe, Np, Nnum, Ndist, Nmonth, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inumonth, Anumonth, area3d, length3d, logFile);
	error = evaluateOptimizationFunction5(Bfactor, concentration, fractionSize, rmax, v, Numonth, Fmonth, Ee, dFe, Np, Nnum, Ndist, Nmonth, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inumonth, Anumonth, area3d, length3d);

	///////////////////


	//error = evaluateOptimizationFunction5(Bfactor, concentration, fractionSize, rmax, v, Numonth, Fmonth, Ee, dFe, Np, Nnum, Ndist, Nmonth, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inumonth, Anumonth, area3d, length3d);

	printf("integrating fields\n");
	fprintf(logFile, "integrating Fields\n");
	fflush(logFile);
	double* totalInu = new double[Nnu];
	double finalSigma = sqr(Bfactor)/(4*pi*concentration*massProtonReal*speed_of_light2);
	printf("Bfactor = %g, n = %g\n fraction = %g rmax = %g sigma = %g\n", Bfactor, concentration, fractionSize, rmax, finalSigma);
	printf("error = %g\n", error);
	fprintf(logFile, "Bfactor = %g, n = %g fraction = %g rmax = %g v/c = %g sigma = %g\n", Bfactor, concentration, fractionSize, rmax, v/speed_of_light, finalSigma);
	fprintf(logFile, "error = %g\n", error);
	fflush(logFile);

	double** tempTotalInu = new double*[Nmonth];

	evaluateVolumeAndLength(area3d, length3d, rmax, fractionSize);

	double**** tempInu = new double***[Nnu];
	double**** tempAnu = new double***[Nnu];
	for(int l = 0; l < Nnu; ++l){
		tempInu[l] = new double**[Nrho];
		tempAnu[l] = new double**[Nrho];
		for(int i = 0; i < Nrho; ++i){
			tempInu[l][i] = new double*[Nphi];
			tempAnu[l][i] = new double*[Nphi];
			for(int j = 0; j < Nphi; ++j){
				tempInu[l][i][j] = new double[Nz];
				tempAnu[l][i][j] = new double[Nz];
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

	for(int l = 0; l < Nmonth; ++l){
		double r = rmax + v*times[l];
		double rfactor = r/rmax;
		double locB = Bfactor*rmax/r;
		double locN = concentration*sqr(rmax/r);

		tempTotalInu[l] = new double[Nnu];

		//evaluateVolumeAndLength(area, length, r, Rho, Phi, fractionSize);
		//evaluateAllEmissivityAndAbsorption1(nu, Inu, Anu, Nnu, Ee, dFe, Np, Ndist, B, sintheta, thetaIndex, concentrations, locN, locB, fractionSize);
		//evaluateSpectrum(nu, tempTotalInu[l], Inu, Anu, area, length, Nnu, r, Rho, Phi);

		evaluateAllEmissivityAndAbsorption(nu, tempInu, tempAnu, Nnu, Ee, dFe, Np, Ndist, B3d, sintheta3d, thetaIndex3d, concentrations3d, concentration, Bfactor, rfactor);

		if(l == 0) {
			if(geometry == SPHERICAL){
				evaluateImageSpherical(image, Numonth[l], tempInu, tempAnu, rmax, Nnum, rfactor, fractionSize);
			} else {
				evaluateImageFlat(image, Numonth[l], tempInu, tempAnu, rmax, Nnum, rfactor, fractionSize);
			}
		}
		//evaluateSpectrum(nu, tempTotalInu[l], tempInu, tempAnu, area3d, length3d, Nnu, rfactor);
		if(geometry == SPHERICAL){
			evaluateSpectrumSpherical(nu, tempTotalInu[l], tempInu, tempAnu, rmax, Nnu, rfactor, fractionSize);
		} else {
			evaluateSpectrumFlat(nu, tempTotalInu[l], tempInu, tempAnu, rmax, Nnu, rfactor, fractionSize);
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
			double tempError = evaluateOptimizationFunction5(Bfactor, concentration, fractionSize, rmax, v, Numonth, Fmonth, Ee, dFe, Np, Nnum, Ndist, Nmonth, B3d, sintheta3d, thetaIndex3d, concentrations3d, Inumonth, Anumonth, area3d, length3d);
			fprintf(errorFile, "%g ", tempError);
		}
		fprintf(errorFile, "\n");
	}
	fclose(errorFile);

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
		for(int j = 0; j < Nmonth; ++j){
			if(tempTotalInu[j][i] != tempTotalInu[j][i]){
				tempTotalInu[j][i] = 0;
			}
		}
		//fprintf(output, "%g %g %g %g %g %g %g\n", nu[i]/1E9, tempTotalInu[0][i], tempTotalInu[1][i] , tempTotalInu[2][i], tempTotalInu[3][i], tempTotalInu[4][i], tempTotalInu[5][i]);
		fprintf(output, "%g %g %g %g %g\n", nu[i]/1E9, tempTotalInu[0][i], tempTotalInu[1][i] , tempTotalInu[2][i], tempTotalInu[3][i]);
	}

	fclose(output);

	for(int l = 0; l < Nmonth; ++l){
		delete[] tempTotalInu[l];
	}
	delete[] tempTotalInu;
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


	FILE* outputParam = fopen("parameters.dat","w");
	fprintf(outputParam, "%g\n", Bfactor);
	fprintf(outputParam, "%g\n", concentration);
	fprintf(outputParam, "%g\n", fractionSize);
	fprintf(outputParam, "%g\n", rmax);
	fprintf(outputParam, "%g\n", v/speed_of_light);
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

	for(int m = 0; m < Nmonth; ++m){
		for(int l = 0; l < Nnum; ++l){
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
			delete[] area3d[i][j];
			delete[] length3d[i][j];
			delete[] B3d[i][j];
			delete[] Bx3d[i][j];
			delete[] By3d[i][j];
			delete[] Bz3d[i][j];
			delete[] concentrations3d[i][j];
			delete[] sintheta3d[i][j];
			delete[] thetaIndex3d[i][j];
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

	for(int i = 0; i < Ntheta; ++i){
		for(int j = 0; j < Nphi; ++j){
			delete[] Inu[i][j];
			delete[] Anu[i][j];
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

	delete[] Inuflat;
	delete[] Anuflat;

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

