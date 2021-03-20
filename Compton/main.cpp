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


int main()
{
	const int Ndist = 10;
	const int Np = 200;
	const int Nnu = 400;


	const int Nphi = 20;
	const int Ntheta = 20;
	double cosTheta[Ntheta];
	double phi[Nphi];

	double dcosTheta = 2.0/Ntheta;
	double dphi = 2.0*pi/Nphi;
	for(int i = 0; i < Ntheta; ++i){
		cosTheta[i] = 1.0 - (i + 0.5)*dcosTheta;
	}
	for(int i = 0; i < Nphi; ++i){
		phi[i] = (i + 0.5)*dphi;
	}

	double Tphotons = 30;
	double* Fph = new double[Np];
	double* dFph = new double[Np];
	double* Eph = new double[Np];

	double Ephmin = 0.1*kBoltzman*Tphotons;
	double Ephmax = kBoltzman*Tphotons*100;

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

	double re2 = sqr(electron_charge*electron_charge/(massElectron*speed_of_light2));

	for(int i = 0; i < Nnu; ++i){
		printf("i nu = %d\n", i);
		double photonFinalEnergy = E[i];
		for(int j = 0; j < Ntheta; ++j){
			//integration by phi changes to 2 pi?
			double photonFinalCosTheta = cosTheta[j];
			for(int k = 0; k < Np; ++k){
				double electronInitialEnergy = Ee[iangle][k];
				double electronInitialGamma = electronInitialEnergy/(massElectron*speed_of_light2);
				double electronInitialBeta = sqrt(1.0 - 1.0/(electronInitialGamma*electronInitialGamma));
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
					double photonInitialCosTheta = cosTheta[l];
					double photonInitialCosThetaPrimed = (photonInitialCosTheta - electronInitialBeta)/(1.0 - electronInitialBeta*photonInitialCosTheta);
					double photonInitialSinThetaPrimed = sqrt(1.0 - photonInitialCosThetaPrimed*photonInitialCosThetaPrimed);
					for(int m = 0; m < Nphi; ++m){
						double photonInitialPhi = phi[m];
						double cosXiPrimed = photonInitialCosThetaPrimed*photonFinalCosThetaPrimed + photonInitialSinThetaPrimed*photonFinalSinThetaPrimed*cos(photonInitialPhi);

						double photonInitialEnergyPrimed = photonFinalEnergyPrimed/(1.0 - (photonFinalEnergyPrimed/(massElectron*speed_of_light2))*(1.0 - cosXiPrimed));

						double photonInitialEnergy = electronInitialGamma*photonInitialEnergyPrimed + electronInitialBeta*electronInitialGamma*photonInitialEnergyPrimed*photonInitialCosThetaPrimed;

						I[i] += (re2*speed_of_light*(1.0 - electronInitialBeta*photonInitialCosTheta)/(2*electronInitialGamma*(1.0 - electronInitialBeta*photonFinalCosTheta)))*
							(1.0 + cosXiPrimed*cosXiPrimed + sqr(photonFinalEnergyPrimed/(massElectron*speed_of_light2))*sqr(1.0 - cosXiPrimed)/(1.0 - (photonFinalEnergyPrimed/(massElectron*speed_of_light2))*(1.0 - cosXiPrimed)))*
							2*pi*dcosTheta*dphi*dcosTheta*delectronEnergy*Fe[iangle][k]*evaluatePhotonDistribution(photonInitialEnergy, Np, Eph, Fph);
						if(I[i] != I[i]){
							printf("I[i] = NaN\n");
							exit(0);
						}
					}
				}
			}
		}
	}

	FILE* output = fopen("output.dat","w");
	for(int i = 0; i < Nnu; ++i){
		double nu2 = sqr(E[i]/hplank);
		fprintf(output, "%g %g\n", E[i]/(1.6E-12), nu2*I[i]);
	}
	fclose(output);

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

