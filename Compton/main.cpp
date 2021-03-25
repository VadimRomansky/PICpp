// Compton.cpp : Defines the entry point for the console application.
//

#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include <time.h>
#include <cmath>
#include <string>

#include "constants.h"
#include "util.h"
#include "startparameters.h"

void LorentzTransformationPhotonZ(const double& beta, const double& Einit, const double& cosThetaInit, double& Eprime, double& cosThetaPrime){
	double gamma = 1.0/sqrt(1.0 - beta*beta);
	double p = Einit/speed_of_light;
	double pz = p*cosThetaInit;
	double pnorm = sqrt(p*p - pz*pz);

	Eprime = gamma*Einit - beta*gamma*pz*speed_of_light;
	double pzprime = -beta*gamma*Einit/speed_of_light + gamma*pz;
	double pprime = sqrt(pnorm*pnorm + pzprime*pzprime);

	cosThetaPrime = pzprime/pprime;
}

double evaluatePhotonDistribution(const double& energy, int Np, double* Eph, double* Fph){
	if(energy <= Eph[0]){
		return 0;
	} else if(energy >= Eph[Np-1]){
		return 0;
	} else {
		int currentI = 0;
		int nextI = 1;
		for(int i = 1; i < Np; ++i){
			currentI = i-1;
			nextI = i;
			if(Eph[i] > energy){
				break;
			}
		}
		double result = (Fph[currentI]*(Eph[nextI] - energy) + Fph[nextI]*(energy - Eph[currentI]))/(Eph[nextI] - Eph[currentI]);
		return result;
	}
}

void initializeParker(int Nrho, int Ntheta, int Nphi, double*** Bx, double*** By, double*** Bz,double* rho, double* cosTheta, double* sinPhiValue, double* cosPhiValue){
	double thetaObserv = 0;
	double sinThetaObserv = sin(thetaObserv);
	double cosThetaObserv = cos(thetaObserv);
	double rcorot = rho[Nrho-1]/10.0;
	double rmin = rho[Nrho-1]/100.0;
	double Br1 = sqr(rmin/rho[Nrho-1]);
	double Bphi1 = ((rho[Nrho-1] - rmin)/rcorot)*sqr(rmin/rho[Nrho-1]);
	//norm to 1;
	double B0 = 1.0/sqrt(Br1*Br1 + Bphi1*Bphi1);
	for(int i = 0; i < Nrho; ++i){
		double r = rho[i];
		for(int j = 0; j < Ntheta; ++j){
			for(int k = 0; k < Nphi; ++k){
				double z = r*cosTheta[j];
				double rxy = sqrt(r*r - z*z);

				double x = rxy*cosPhiValue[k];
				double y = rxy*sinPhiValue[k];

				double y1 = y;
				double z1 = z*cosThetaObserv + x*sinThetaObserv;
				double x1 = x*cosThetaObserv - z*sinThetaObserv;

				double costheta = z1/r;
				double sintheta = sqrt(1.0 - costheta*costheta);

				double sinphi = sinPhiValue[k];
				double cosphi = cosPhiValue[k];

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
					printLog("Bx = NaN\n");
					exit(0);
				}
			}
		}
	}
}

void evaluateOrientationParameters3d(int Nrho, int Ntheta, int Nphi, double*** B, double*** sintheta, int*** thetaIndex, double*** Bx, double*** By, double*** Bz, int Nd, double* rho, double* cosThetaValue, double* sinPhiValue, double* cosPhiValue){
	for(int i = 0; i < Nrho; ++i){
		for(int j = 0; j < Ntheta; ++j){
			for(int k = 0; k < Nphi; ++k){
				//double r = rmax;
				double Bxy = sqrt(Bx[i][j][k]*Bx[i][j][k] + By[i][j][k]*By[i][j][k]);
				B[i][j][k] = sqrt(Bz[i][j][k]*Bz[i][j][k] + Bx[i][j][k]*Bx[i][j][k] + By[i][j][k]*By[i][j][k]);
				double cosTheta = Bz[i][j][k]/B[i][j][k];
				double sinTheta = Bxy/B[i][j][k];
				if(sinTheta != sinTheta){
					printf("sintheta NaN\n");
					printLog("sinTheta = NaN\n");
					exit(0);
				}
				sintheta[i][j][k] = sinTheta;

				double r = rho[i];

				double cosTheta1 = cosThetaValue[k];
				double sinTheta1 = sqrt(1.0 - cosTheta1*cosTheta1);
				
				double cosPhi1 = cosPhiValue[k];
				double sinPhi1 = sinPhiValue[k];

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
				if(thetaIndex[i][j][k] < 0 || thetaIndex[i][j][k] > Nd){
					FILE* logFile = fopen("log.dat","a");
					printf("thetaIndex = %d\n", thetaIndex[i][j][k]);
					fprintf(logFile, "thetaIndex = %d\n", thetaIndex[i][j][k]);
					fclose(logFile); 
					exit(0);
				}
				//for debug
				//thetaIndex[i][j][k] = 9;
			}
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


int main()
{
	const int Ndist = 10;
	const int Np = 200;
	const int Nnu = 200;

	const double electronConcentration = 1.0;
	const double photonConcentration = 1.0;


	const int Nphi = 10;
	const int Ntheta = 20;
	double cosThetaLeft[Ntheta];
	double cosTheta[Ntheta];
	double phi[Nphi];
	double sinPhiValue[Nphi];
	double cosPhiValue[Nphi];

	double dcosTheta = 2.0/Ntheta;
	double dphi = 2.0*pi/Nphi;
	for(int i = 0; i < Ntheta; ++i){
		cosTheta[i] = 1.0 - (i + 0.5)*dcosTheta;
		cosThetaLeft[i] = 1.0 - (i)*dcosTheta;
	}
	for(int i = 0; i < Nphi; ++i){
		phi[i] = (i + 0.5)*dphi;
		sinPhiValue[i] = sin(phi[i]);
		cosPhiValue[i] = cos(phi[i]);
	}

	double Tphotons = 30;
	double* Fph = new double[Np];
	double* dFph = new double[Np];
	double* Eph = new double[Np];

	double Ephmin = 0.1*kBoltzman*Tphotons;
	double Ephmax = kBoltzman*Tphotons*100;

	FILE* logFile = fopen("log.dat", "w");
	fclose(logFile);

	printLog("start\n");

	double factor = pow(Ephmax/Ephmin, 1.0/(Np - 1));
	Eph[0] = Ephmin;
	Fph[0] = 0;
	dFph[0] = 0;
	//todo check and normalize
	for(int i = 1; i < Np; ++i){
		Eph[i] = Eph[i-1]*factor;
		double theta = Eph[i]/(kBoltzman*Tphotons);
		Fph[i] =  (2*Eph[i]*Eph[i]/(speed_of_light2*hplank*hplank))/(exp(theta)-1.0);
		dFph[i] = Fph[i]*(Eph[i] - Eph[i-1]);
	}

	printLog("read electrons distribution\n");
	double** Fe = new double*[Ndist];
	double** dFe = new double*[Ndist];
	double** Ee = new double*[Ndist];
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
				double power = 2;

				if(gamma >= minGamma){
					Fe[j][i] = 1.0/pow(Ee[j][i],power);
				} else {
					Fe[j][i] = 0;
				}
			}
		}


		fclose(inputPe);
		fclose(inputFe);
	}

	printLog("initialize fields\n");
	const int Nrho = 10;
	double rmax = 3.6E16;
	double fraction = 0.2;
	double rmin = rmax*(1.0 - fraction);
	double dr = (rmax - rmin)/Nrho;
	double rho[Nrho];

	double*** B3d = new double**[Nrho];
	double*** Bx3d = new double**[Nrho];
	double*** By3d = new double**[Nrho];
	double*** Bz3d = new double**[Nrho];
	double*** concentrations3d = new double**[Nrho];
	double*** sintheta3d = new double**[Nrho];
	int*** thetaIndex3d = new int**[Nrho];
	double*** volume = new double**[Nrho];
	for(int i = 0; i < Nrho; ++i){
		double r = rmin + (i+0.5)*dr;
		rho[i] = r;
		volume[i] = new double*[Ntheta];
		B3d[i] = new double*[Ntheta];
		Bx3d[i] = new double*[Ntheta];
		By3d[i] = new double*[Ntheta];
		Bz3d[i] = new double*[Ntheta];
		concentrations3d[i] = new double*[Ntheta];
		sintheta3d[i] = new double*[Ntheta];
		thetaIndex3d[i] = new int*[Ntheta];
		for(int j = 0; j < Ntheta; ++j){
			volume[i][j] = new double[Nphi];
			B3d[i][j] = new double[Nphi];
			Bx3d[i][j] = new double[Nphi];
			By3d[i][j] = new double[Nphi];
			Bz3d[i][j] = new double[Nphi];
			concentrations3d[i][j] = new double[Nphi];
			sintheta3d[i][j] = new double[Nphi];
			thetaIndex3d[i][j] = new int[Nphi];
			for(int k = 0; k < Nphi; ++k){
				volume[i][j][k] = ((pow(r+dr/2,3) - pow(r-dr/2,3))/3.0)*dphi*dcosTheta;
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

	printLog("initialize Parker\n");
	if(parker){
		initializeParker(Nrho, Ntheta, Nphi, Bx3d, By3d, Bz3d, rho, cosTheta, sinPhiValue, cosPhiValue);
	}

	printLog("initialize turbulence\n");
	srand(time(NULL));
	int randomSeed = rand();
	randomSeed = 10;
	printf("random seed = %d\n", randomSeed);
	srand(randomSeed);

	const int Nk = 10;

	double kmin = 2*pi*2/(0.5);
	double dk = kmin;
	double kmax = Nk*dk;
	double turbulenceFraction = 0.9;
	double turbNorm = evaluateTurbNorm(kmax, Nk, 1.0, turbulenceFraction);

	if(turbulence){
		for(int i = 0; i < Nrho; ++i){
			for(int j = 0; j < Ntheta; ++j){
				for(int k = 0; k < Nphi; ++k){
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
						//double cosTheta1 = 0;
						//double sinTheta1 = 1.0;
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
							for(int j = 0; j < Ntheta; ++j){
								for(int k = 0; k < Nphi; ++k){
									double z = rho[i]*cosTheta[j];
									double rxy = sqrt(rho[i]*rho[i] - z*z);
									double x = rxy*cosPhiValue[k];
									double y = rxy*sinPhiValue[k];

									double kmultr = kx*x + ky*y;
									double localB1 = Bturbulent*sin(kmultr + phase1)*cos(pi*z/(2*rho[Nrho-1]));
									double localB2 = Bturbulent*sin(kmultr + phase2)*cos(pi*z/(2*rho[Nrho-1]));
									localB1 = 0;
									localB2 = localB2*sqrt(2.0);

									Bz3d[i][j][k] = Bz3d[i][j][k] - localB1;
									Bx3d[i][j][k] = Bx3d[i][j][k] - localB2*sinPhi;
									By3d[i][j][k] = By3d[i][j][k] + localB2*cosPhi;
								}
							}
						}
					//}
				}
			}
		}
	}

	printLog("evaluate orientation parameters\n");
	evaluateOrientationParameters3d(Nrho, Ntheta, Nphi, B3d, sintheta3d, thetaIndex3d, Bx3d, By3d, Bz3d, Ndist, rho, cosTheta, sinPhiValue, cosPhiValue);

	int iangle = 3;

	double* E = new double[Nnu];
	double* I = new double[Nnu];
	double Emin = 0.000001*kBoltzman*Tphotons;
	double Emax = Ee[iangle][Np-1];
	factor = pow(Emax/Emin, 1.0/(Nnu - 1));
	E[0] = Emin;
	I[0] = 0;
	for(int i = 1; i < Nnu; ++i){
		E[i] = E[i-1]*factor;
		I[i] = 0;
	}

	//for debug
	/*for(int i = 0; i < Np; ++i){
		Fe[iangle][i] = 0;
		Fph[i] = 0;
	}
	Fe[iangle][Np/2] = 1.0;
	Fph[Np/2] = 1.0;*/
	/////
	printLog("evaluate spectrum\n");
	double re2 = sqr(electron_charge*electron_charge/(massElectron*speed_of_light2));

	#pragma omp parallel for shared(logFile)
	for(int i = 0; i < Nnu; ++i){
		logFile = fopen("log.dat","a");
		printf("i nu = %d\n", i);
		fprintf(logFile, "i nu = %d\n", i);
		fclose(logFile);
		double photonFinalEnergy = E[i];
		for(int j = 0; j < Ntheta; ++j){
			//integration by phi changes to 2 pi?
			double photonFinalCosTheta = cosThetaLeft[j];
			//int ir = 0;
			//int itheta = 0;
			//int iphi = 0;
			for(int ir = 0; ir < Nrho; ++ir){
				for(int itheta = 0; itheta < Ntheta; ++itheta){
					for(int iphi = 0; iphi < Nphi; ++iphi){
						iangle = thetaIndex3d[ir][itheta][iphi];
						for(int k = 0; k < Np; ++k){
							double electronInitialEnergy = Ee[iangle][k];
							/*logFile = fopen("log.dat","a");
							printf("electronInitialEnergy = %g\n", electronInitialEnergy);
							fprintf(logFile, "electronInitialEnergy = %g\n", electronInitialEnergy);
							fclose(logFile);*/
							double electronInitialGamma = electronInitialEnergy/(massElectron*speed_of_light2);
							double electronInitialBeta = sqrt(1.0 - 1.0/(electronInitialGamma*electronInitialGamma));
							/*logFile = fopen("log.dat","a");
							printf("electronInitialBeta = %g\n", electronInitialBeta);
							fprintf(logFile, "electronInitialBeta = %g\n", electronInitialBeta);
							fclose(logFile);*/
							double delectronEnergy;
							if(k == 0){
								delectronEnergy = Ee[iangle][1] - Ee[iangle][0];
							} else {
								delectronEnergy = Ee[iangle][k] - Ee[iangle][k-1];
							}

							double photonFinalEnergyPrimed;
							double photonFinalCosThetaPrimed;
							LorentzTransformationPhotonZ(electronInitialBeta, photonFinalEnergy, photonFinalCosTheta, photonFinalEnergyPrimed, photonFinalCosThetaPrimed);
							double photonFinalSinThetaPrimed = sqrt(1.0 - photonFinalCosThetaPrimed*photonFinalCosThetaPrimed);
							for(int l = 0; l < Ntheta; ++l){
								/*logFile = fopen("log.dat","a");
								printf("l photon initial theta = %d\n", l);
								fprintf(logFile, "l photon initial theta = %d\n", l);
								fclose(logFile);*/
								double photonInitialCosTheta = cosThetaLeft[l];
								double photonInitialCosThetaPrimed = (photonInitialCosTheta - electronInitialBeta)/(1.0 - electronInitialBeta*photonInitialCosTheta);
								double photonInitialSinThetaPrimed = sqrt(1.0 - photonInitialCosThetaPrimed*photonInitialCosThetaPrimed);
								for(int m = 0; m < Nphi; ++m){
									double photonInitialPhi = phi[m];
									double cosXiPrimed = photonInitialCosThetaPrimed*photonFinalCosThetaPrimed + photonInitialSinThetaPrimed*photonFinalSinThetaPrimed*cos(photonInitialPhi);

									double photonInitialEnergyPrimed = photonFinalEnergyPrimed/(1.0 - (photonFinalEnergyPrimed/(massElectron*speed_of_light2))*(1.0 - cosXiPrimed));

									double photonInitialEnergy = electronInitialGamma*photonInitialEnergyPrimed + electronInitialBeta*electronInitialGamma*photonInitialEnergyPrimed*photonInitialCosThetaPrimed;

									I[i] += photonConcentration*electronConcentration*concentrations3d[ir][itheta][iphi]*volume[ir][itheta][iphi]*(re2*speed_of_light*(1.0 - electronInitialBeta*photonInitialCosTheta)/(2*electronInitialGamma*(1.0 - electronInitialBeta*photonFinalCosTheta)))*
										(1.0 + cosXiPrimed*cosXiPrimed + sqr(photonFinalEnergyPrimed/(massElectron*speed_of_light2))*sqr(1.0 - cosXiPrimed)/(1.0 - (photonFinalEnergyPrimed/(massElectron*speed_of_light2))*(1.0 - cosXiPrimed)))*
										2*pi*dcosTheta*dphi*dcosTheta*delectronEnergy*Fe[iangle][k]*evaluatePhotonDistribution(photonInitialEnergy, Np, Eph, Fph);
									if(I[i] != I[i]){
										printf("I[i] = NaN\n");
										printLog("I[i] = NaN\n");
										exit(0);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	printLog("output\n");

	FILE* output = fopen("output.dat","w");
	for(int i = 0; i < Nnu; ++i){
		double nu = E[i]/hplank;
		fprintf(output, "%g %g\n", E[i]/(1.6E-12), nu*E[i]*I[i]/sqr(distance));
	}
	fclose(output);

	printLog("delete arrays\n");

	for(int i = 0; i < Nrho; ++i){
		for(int j = 0; j < Nphi; ++j){
			delete[] volume[i][j];
			delete[] B3d[i][j];
			delete[] Bx3d[i][j];
			delete[] By3d[i][j];
			delete[] Bz3d[i][j];
			delete[] concentrations3d[i][j];
			delete[] sintheta3d[i][j];
			delete[] thetaIndex3d[i][j];
		}
		delete[] volume[i];
		delete[] B3d[i];
		delete[] Bx3d[i];
		delete[] By3d[i];
		delete[] Bz3d[i];
		delete[] concentrations3d[i];
		delete[] sintheta3d[i];
		delete[] thetaIndex3d[i];
	}
	delete[] volume;
	delete[] B3d;
	delete[] Bx3d;
	delete[] By3d;
	delete[] Bz3d;
	delete[] concentrations3d;
	delete[] sintheta3d;
	delete[] thetaIndex3d;

	for(int i = 0; i < Ndist; ++i){
		delete[] Fe[i];
		delete[] dFe[i];
		delete[] Ee[i];
	}
	delete[] Fe;
	delete[] dFe;
	delete[] Ee;
	delete[] Fph;
	delete[] dFph;
	delete[] Eph;

	delete[] E;
	delete[] I;

	return 0;
}

