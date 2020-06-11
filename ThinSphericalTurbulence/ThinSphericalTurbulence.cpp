// SphericalTurbulence.cpp : Defines the entry point for the console application.
//
#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include <time.h>
#include <cmath>
#include <string>
#include <omp.h>

#include "constants.h"
#include "util.h"
#include "optimization.h"
#include "spectrum.h"


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
				thetaIndex[j][k] = 3;
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

	

	//double rmin = dx*sqrt(3.0)/2;
	double rmax = size[0];
	double rmin = 1E13;

	//double B0 = 1000;
	double B0 = 1.0;

	double rcorot = rmax/5;

	FILE* logFile = fopen(logFileName.c_str(), "w");

	printf("initialization\n");
	fprintf(logFile, "initialization\n");

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

	Bx = new double*[Ntheta];
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
	}
	//////////////////////////////////////

	printf("evaluating turbulent field\n");
	fprintf(logFile, "evaluating turbulent field\n");

	srand(time(NULL));
	int randomSeed = rand();
	randomSeed = 10;
	printf("random seed = %d\n", randomSeed);
	fprintf(logFile, "random seed = %d\n", randomSeed);
	srand(randomSeed);

	double kmin = 2*pi*2/rmax;
	double dk = kmin;
	double kmax = Nk*dk;
	double turbNorm = evaluateTurbNorm(kmax, Nk, Bx[Ntheta/2][0], 0.9);

	/*for(int ki = 0; ki < Nk; ++ki){
		printf("%d\n", ki);
		for(int kj = 0; kj < Nk; ++kj){
			for(int kk = 0; kk < Nk; ++kk){
				if ((ki + kj + kk) > 4) {
					double phase1 = 2*pi*uniformDistribution();
					double phase2 = 2*pi*uniformDistribution();


					double kx = ki*dk;
					double ky = kj*dk;
					double kz = kk*dk;

					double kw = sqrt(kx*kx + ky*ky + kz*kz);
					double kxy = sqrt(ky*ky + kx*kx);
					double cosTheta = kz/kw;
					double sinTheta = kxy/kw;
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
						for(int j = 0; j < Ny; ++j){
							for(int k = 0; k < Nz; ++k){

								double kmultr = kx*X[i] + ky*Y[j] + kz*Z[k];
								double localB1 = 0.3*Bturbulent*sin(kmultr + phase1);
								double localB2 = 0.3*Bturbulent*sin(kmultr + phase2);

								Bz[i][j][k] = Bz[i][j][k] - localB1*sinTheta;
								Bx[i][j][k] = Bx[i][j][k] + (localB1*cosTheta*cosPhi - localB2*sinPhi);
								By[i][j][k] = By[i][j][k] + (localB1*cosTheta*sinPhi + localB2*cosPhi);
							}
						}
					}
				}

			}
		}
	}

		for(int i = 0; i < Nx; ++i){
			for(int j = 0; j < Ny; ++j){
				for(int k = 0; k < Nz; ++k){
					double r = sqrt(X[i]*X[i] + Y[j]*Y[j] + Z[k]*Z[k]);
					if(r > 2*rmax){
						Bx[i][j][k] = 0;
						By[i][j][k] = 0;
						Bz[i][j][k] = 0;
					}
					if(r < 1.2*rmax){
						Bx[i][j][k] = 0;
						By[i][j][k] = 0;
						Bz[i][j][k] = 0;
					}
				}
			}
		}*/

	FILE* bFile = fopen(BFileName.c_str(), "w");
	for(int j = 0; j < Ntheta; ++j){
		for(int k = 0; k < Nphi; ++k){
			fprintf(bFile, "%g %g %g\n", Bx[j][k], By[j][k], Bz[j][k]);
		}
	}
	fclose(bFile);

	evaluateOrientationParameters(B, sintheta, thetaIndex, Bx, By, Bz, Ndist);

	/*double sigma = 0.04;
	double tempConcentration = sqr(B[Ntheta/2][0])/(sigma*4*pi*massProtonReal*speed_of_light2);
	for(int i = 0; i < Ntheta; ++i){
		for(int j = 0; j < Nphi; ++j){
				concentrations[i][j] = tempConcentration;
		}
	}*/

	printf("reading input\n");
	fprintf(logFile, "reading input\n");

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
		for (int i = 1; i < Np; ++i) {
			fscanf(inputPe, "%lf", &u);
			fscanf(inputFe, "%lf", &Fe[j][i]);

			u = u*massSpectrumFactor;
			//if( u < 3000){
				double gamma = sqrt(u * u  + 1);
				Ee[j][i] = sqrt(u * u  + 1)*massElectron*speed_of_light2;
				//maxEnergy = Ee[i];
				Fe[j][i] = Fe[j][i] * Ee[j][i]/ (u * u * u * speed_of_light4 *massElectron*massElectron);
				if(gamma >= 1.0){
					Fe[j][i] = 1.0/pow(Ee[j][i],3);
				} else {
					Fe[j][i] = 0;
				}
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
			//} else {
			//	Fe[i] = 0;
			//}
		
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
	for(int i = 0; i < Nnu1; ++i){
		nu1[i] = aprx[i]*1.0E9;
	}

	/////////////////
	//todo concentration!!
	double sigma = 0.004;
	double tempConcentration = sqr(B[Ntheta/2][0])/(sigma*4*pi*massProtonReal*speed_of_light2);
	double concentration = 1.0*tempConcentration;
	concentration = 3300;
	double Bfactor = 0.48;
	//double fractionSize = 1.0 - pow((3.0/(4.0*pi))*(1 - 0.5),1.0/3.0);
	double fractionSize = 0.61;
	double V0 = speed_of_light;
	double v = 0.75*speed_of_light;
	rmax = 4.5E16;
	////////////////////
	printf("initial optimizing parameters\n");
	fprintf(logFile, "initial optimizing parameters\n");
	fflush(logFile);
	printf("optimizing parameters\n");
	fprintf(logFile, "optimizing parameters\n");
	//optimizeParameters(Bfactor, concentration, fractionSize, nu1, rmax, Ee, dFe, Np, Nnu1, Nd, B, sintheta, thetaIndex, concentrations, Inu1, Anu1, area, length, Rho, Phi, logFile);
	optimizeParameters4(1.0, 2000, 3.4E16, Bfactor, concentration, fractionSize, rmax, nu1, Ee, dFe, Np, Nnu1, Ndist, B, sintheta, thetaIndex, concentrations, Inu1, Anu1, area, length, Rho, Phi, logFile);
	//optimizeParameters4(Bfactor, concentration, fractionSize,rmax, nu1, Ee, dFe, Np, Nnu1, Nd, B, sintheta, thetaIndex, concentrations, Inu1, Anu1, area, length, Rho, Phi, logFile);
	//optimizeParameters5(1.0, 2000, 3.4E16, V0, Bfactor, concentration, fractionSize, rmax, v, Ee, Fe, Np, Ndist, 1.0, 3, logFile);
	///////////////////
	//concentration = 1.0;
	//Bfactor = 1.0;
	//fractionSize = 0.1;

	printf("integrating fields\n");
	fprintf(logFile, "integrating Fields\n");
	double* totalInu = new double[Nnu];
	double finalSigma = sqr(B[Ntheta/2][0]*Bfactor)/(4*pi*concentration*concentrations[Ntheta/2][0]*massProtonReal*speed_of_light2);
	printf("Bfactor = %g, n = %g\n fraction = %g rmax = %g sigma = %g\n", Bfactor, concentration, fractionSize, rmax, finalSigma);
	fprintf(logFile, "Bfactor = %g, n = %g fraction = %g rmax = %g sigma = %g\n", Bfactor, concentration, fractionSize, rmax, finalSigma);
	fflush(logFile);

	double** tempTotalInu = new double*[4];


	for(int l = 0; l < 4; ++l){
		double r = rmax + v*times[l];
		double locB = Bfactor*rmax/r;
		double locN = concentration*sqr(rmax/r);


		tempTotalInu[l] = new double[Nnu];

		evaluateVolumeAndLength(area, length, r, Rho, Phi, fractionSize);
		evaluateAllEmissivityAndAbsorption1(nu, Inu, Anu, Nnu, Ee, dFe, Np, Ndist, B, sintheta, thetaIndex, concentrations, locN, locB, fractionSize);
		evaluateSpectrum(nu, tempTotalInu[l], Inu, Anu, area, length, Nnu, r, Rho, Phi);

		//evaluateEmissivityAndAbsorptionFlat(nu, Inuflat, Anuflat, Ee[3], dFe[3], Np, Nnu, 1.0, locB, locN);
		//evaluateSpectrumFlat(nu, tempTotalInu[l], Inuflat, Anuflat, Nnu, r, fractionSize);
	}

	//////////
	printf("outputing\n");
	fprintf(logFile, "ooutputing\n");
	//char* number = new char[100];
	//itoa(0, number, 10);
	//delete[] number;
	//std::string fileNumber = std::string(number);
	FILE* output = fopen(outputfileName.c_str(), "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(output, "%g %g %g %g %g\n", nu[i]/1E9, tempTotalInu[0][i], tempTotalInu[1][i] , tempTotalInu[2][i], tempTotalInu[3][i]);
	}

	fclose(output);

	for(int l = 0; l < 4; ++l){
		delete[] tempTotalInu[l];
	}
	delete[] tempTotalInu;

	FILE* outputParam = fopen("parameters.dat","w");
	fprintf(outputParam, "%g\n", Bfactor);
	fprintf(outputParam, "%g\n", concentration);
	fprintf(outputParam, "%g\n", fractionSize);
	fclose(outputParam);

	////////////////////

	printf("deleting arrays\n");
	fprintf(logFile, "deleting arrays\n");
	delete[] totalInu;

	for(int j = 0; j < Ndist; ++j){
		delete[] Fe[j];
		delete[] Ee[j];
	}
	delete[] Fe;
	delete[] Ee;

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

