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

double evaluateMcDonaldIntegral(const double& nu) {
	int curIndex = 0;
	if (nu < UvarovX[0]) {
		printf("x < UvarovX[0]\n");
		return 0;
	}
	if (nu > UvarovX[Napprox - 1]) {
		//printf("x > UvarovX[Napprox - 1]\n");
		return sqrt(nu*pi/2)*exp(-nu);
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

	//double result = (UvarovValue[rightIndex]*(nu - UvarovX[leftIndex]) + UvarovValue[leftIndex]*(UvarovX[rightIndex] - nu))/(UvarovX[rightIndex] - UvarovX[leftIndex]);
	double result = UvarovValue[leftIndex] * exp(log(UvarovValue[rightIndex] / UvarovValue[leftIndex]) * ((nu - UvarovX[leftIndex]) / (UvarovX[rightIndex] - UvarovX[leftIndex])));
	if (result < 0) {
		printf("result < 0\n");
	}
	return result;
}

double evaluateMcDonaldFunction5_3(const double& nu) {
	int curIndex = 0;
	if (nu < UvarovX[0]) {
		return McDonaldValue5_3[0];
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

	double result = (McDonaldValue5_3[rightIndex]*(nu - UvarovX[leftIndex]) + McDonaldValue5_3[leftIndex]*(UvarovX[rightIndex] - nu))/(UvarovX[rightIndex] - UvarovX[leftIndex]);
	//double result = McDonaldValue[curIndex - 1] * exp(
	//	log(McDonaldValue[curIndex] / McDonaldValue[curIndex - 1]) * ((nu - UvarovX[curIndex - 1]) / (UvarovX[curIndex] - UvarovX[curIndex - 1])));
	if (result < 0) {
		printf("result < 0\n");
	}
	return result;
}

void evaluateMcDonaldFunctions(const double& nu, double& K1_3, double& K2_3, double& K4_3, double& K5_3) {
	int curIndex = 0;
	if (nu < UvarovX[0]) {
		K1_3 = McDonaldValue1_3[0];
		K2_3 = McDonaldValue2_3[0];
		K4_3 = McDonaldValue4_3[0];
		K5_3 = McDonaldValue5_3[0];
		return;
	}
	if (nu > UvarovX[Napprox - 1]) {
		K1_3 = 0;
		K2_3 = 0;
		K4_3 = 0;
		K5_3 = 0;
		return;
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

	K1_3 = (McDonaldValue1_3[rightIndex]*(nu - UvarovX[leftIndex]) + McDonaldValue1_3[leftIndex]*(UvarovX[rightIndex] - nu))/(UvarovX[rightIndex] - UvarovX[leftIndex]);
	K2_3 = (McDonaldValue2_3[rightIndex]*(nu - UvarovX[leftIndex]) + McDonaldValue2_3[leftIndex]*(UvarovX[rightIndex] - nu))/(UvarovX[rightIndex] - UvarovX[leftIndex]);
	K4_3 = (McDonaldValue4_3[rightIndex]*(nu - UvarovX[leftIndex]) + McDonaldValue4_3[leftIndex]*(UvarovX[rightIndex] - nu))/(UvarovX[rightIndex] - UvarovX[leftIndex]);
	K5_3 = (McDonaldValue5_3[rightIndex]*(nu - UvarovX[leftIndex]) + McDonaldValue5_3[leftIndex]*(UvarovX[rightIndex] - nu))/(UvarovX[rightIndex] - UvarovX[leftIndex]);
}

double criticalNu(const double& E, const double& sinhi, const double& H) {
	return criticalNuCoef * H * sinhi * E * E;
}

void evaluateLocalEmissivityAndAbsorption(double* nu, double* Inu, double* Anu, int Nnu, double* Ee, double* Fe, int Np, double sinhi, double B, double concentration, double fractionSize) {
	//Anu from quantum cross-section from??
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
	//double coefAbsorb = concentration * absorpCoef/B;
	//todo what if < 0?
	double coshi = sqrt(1.0 - sinhi*sinhi);


	for (int i = 0; i < Nnu; ++i) {
		//printf("i = %d\n", i);
		double oldA = 0;
		for (int j = 1; j < Np; ++j) {
			//if(Ee[j] < 100*massElectron*speed_of_light2){
			if(Fe[j] > 0){
				double nuc = criticalNu(Ee[j], sinhi, B);
				double gamma = Ee[j] / (massElectron * speed_of_light2);
				double gamma4 = gamma*gamma*gamma*gamma;
				double sigmaCoef = concentration*Fe[j]/(B*gamma4);
				//double x = nu[i] / nuc;
				//todo!!! 4pi!!
				//here dFe is (dFe[j] / (4*pi)) * (Ee[j] - Ee[j - 1])
				Inu[i] = Inu[i] + coef * Fe[j] * B * sinhi * evaluateMcDonaldIntegral(nu[i] / nuc);
				if(Inu[i] < 0){
					printf("Inu[i] < 0\n");
					printf("dFe[j] = %g\n", Fe[j]);
					exit(0);
				}
				// integral for sigma
				double sigmaInt = 0;
				//smallAngles
				double psimax = 2*pi/gamma;
				int Nphi = Npsi;
				double dpsi = psimax/Npsi;
				double dalpha = 2*pi/Nphi;

				for(int l = 0; l < Npsi; ++l){
					double psi = dpsi*(0.5 + l);
					double cospsi = cos(psi);
					double sinpsi = sin(psi);
					for(int m = 0; m < Nphi; ++m){
						double cosalpha = cosAlphaValue[m];
						double sinalpha = sinAlphaValue[m];
						double costheta = coshi*cospsi + sinhi*sinpsi*cosalpha;
						if(fabs(costheta) > 1.0){
							printf("costheta > 1\n");
							costheta = 1.0;
						}
						double sintheta = sqrt(1.0 - costheta*costheta);

						double localNuc = criticalNu(Ee[j], sintheta, B);
						double x = nu[i]/localNuc;

						double t = gamma*gamma*psi*psi;
						double y = 0.5*x*pow(1.0 + t, 1.5); 

						if((x < 0.1) && (t < 0.01)){
							sigmaInt = sigmaInt + (sigmaCoef*absorpCoef/(sintheta*sintheta*pow(x,4.0/3.0)))*sinpsi*dpsi*dalpha;
							if(sigmaInt != sigmaInt){
								printf("sigmaInt Nan\n");
								exit(0);
							}
						} else {
							double K1_3 = 0;
							double K2_3 = 0;
							double K4_3 = 0;
							double K5_3 = 0;
							evaluateMcDonaldFunctions(y, K1_3, K2_3, K4_3, K5_3);

							sigmaInt = sigmaInt + ((sigmaCoef*absorpCoef2/(sintheta*sintheta))*sinpsi*dpsi*dalpha)*(K2_3*K2_3*(13*t - 1)*(t + 1)/3 + K1_3*K1_3*t*(11*t-1)/3 -2*y*(t-2)*((1 + t)*K2_3*K5_3 + t*K1_3*K4_3));
							if(sigmaInt != sigmaInt){
								printf("sigmaInt Nan\n");
								exit(0);
							}
						}
					
					}
				}
				double smallSigmaInt = sigmaInt;

				//large angles
				psimax = 2*pi;
				double psimin = 2*pi/gamma;
				Nphi = Npsi;
				dpsi = (psimax-psimin)/Npsi;
				dalpha = 2*pi/Nphi;

				for(int l = 0; l < Npsi; ++l){
					double psi = psimin + dpsi*(0.5 + l);
					double cospsi = cos(psi);
					double sinpsi = sin(psi);
					for(int m = 0; m < Nphi; ++m){
						double cosalpha = cosAlphaValue[m];
						double sinalpha = sinAlphaValue[m];
						double costheta = coshi*cospsi + sinhi*sinpsi*cosalpha;
						if(fabs(costheta) > 1.0){
							printf("costheta > 1\n");
							costheta = 1.0;
						}
						double sintheta = sqrt(1.0 - costheta*costheta);

						double localNuc = criticalNu(Ee[j], sintheta, B);
						double x = nu[i]/localNuc;

						double t = gamma*gamma*psi*psi;
						double y = 0.5*x*pow(1.0 + t, 1.5); 


						double K1_3 = 0;
						double K2_3 = 0;
						double K4_3 = 0;
						double K5_3 = 0;
						evaluateMcDonaldFunctions(y, K1_3, K2_3, K4_3, K5_3);

						sigmaInt = sigmaInt + ((sigmaCoef*absorpCoef2/(sintheta*sintheta))*sinpsi*dpsi*dalpha)*(K2_3*K2_3*(13*t - 1)*(t + 1)/3 + K1_3*K1_3*t*(11*t-1)/3 -2*y*(t-2)*((1 + t)*K2_3*K5_3 + t*K1_3*K4_3));
						if(sigmaInt != sigmaInt){
							printf("sigmaInt Nan\n");
							exit(0);
						}
					}
				}

				Anu[i] = Anu[i] + sigmaInt;
				oldA = oldA + (16*pi*pi*electron_charge/(3.0*sqrt(3.0)))*Fe[j]*evaluateMcDonaldFunction5_3(nu[i]/nuc)/(gamma*gamma*gamma*gamma*gamma*B*sinhi);
				if(Inu[i] != Inu[i]){
					printf("Inu NaN\n");
					exit(0);
				}
				if(Anu[i] != Anu[i]){
					printf("Anu Nan\n");
					exit(0);
				}
			}
		}
		//Anu[i] = oldA;
		//Anu[i] = 0;
	}

	//for (int i = 0; i < Nnu; ++i) {
	//	Inu[i] = Inu[i] * (1 - exp(-Anu[i] * fractionSize * localSize)) / (Anu[i] * fractionSize * localSize);
	//}
}

void evaluateLocalEmissivityAndAbsorption1(double* nu, double* Inu, double* Anu, int Nnu, double* Ee, double* Fe, int Np, double sinhi, double B, double concentration, double fractionSize) {
	//Anu from ghiselini simple
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
	//double coefAbsorb = concentration * absorpCoef/B;
	//todo what if < 0?
	double coshi = sqrt(1.0 - sinhi*sinhi);


	for (int i = 0; i < Nnu; ++i) {
		//printf("i = %d\n", i);
		double oldA = 0;
		for (int j = 1; j < Np; ++j) {
			//if(Ee[j] < 100*massElectron*speed_of_light2){
			if(Fe[j] > 0){
				double nuc = criticalNu(Ee[j], sinhi, B);
				double gamma = Ee[j] / (massElectron * speed_of_light2);
				double gamma4 = gamma*gamma*gamma*gamma;
				double sigmaCoef = concentration*Fe[j]/(B*gamma4);
				//double x = nu[i] / nuc;
				//todo!!! 4pi!!
				//here dFe is (dFe[j] / (4*pi)) * (Ee[j] - Ee[j - 1])
				Inu[i] = Inu[i] + coef * Fe[j] * B * sinhi * evaluateMcDonaldIntegral(nu[i] / nuc);

				if(Inu[i] < 0){
					printf("Inu[i] < 0\n");
					printf("dFe[j] = %g\n", Fe[j]);
					exit(0);
				}

				double tempP = gamma*gamma*coef * B * sinhi * evaluateMcDonaldIntegral(nu[i] / nuc);
				double dg = 0.1*gamma;
				double tempGamma = gamma + dg;
				double tempnuc = criticalNu(massElectron*speed_of_light2*tempGamma, sinhi, B);
				double tempP2 = tempGamma*tempGamma*coef * B * sinhi * evaluateMcDonaldIntegral(nu[i] / tempnuc);
				double Pder = (tempP2 - tempP)/dg;

				

				Anu[i] = Anu[i] + (1.0/(2*massElectron*nu[i]*nu[i]))*Fe[j]*Pder/(gamma*gamma);
				oldA = oldA + (16*pi*pi*electron_charge/(3.0*sqrt(3.0)))*Fe[j]*evaluateMcDonaldFunction5_3(nu[i]/nuc)/(gamma*gamma*gamma*gamma*gamma*B*sinhi);
				if(Inu[i] != Inu[i]){
					printf("Inu NaN\n");
					exit(0);
				}
				if(Anu[i] != Anu[i]){
					printf("Anu Nan\n");
					exit(0);
				}
			}
		}
		//Anu[i] = oldA;
		//Anu[i] = 0;
	}

	//for (int i = 0; i < Nnu; ++i) {
	//	Inu[i] = Inu[i] * (1 - exp(-Anu[i] * fractionSize * localSize)) / (Anu[i] * fractionSize * localSize);
	//}
}

void evaluateAllEmissivityAndAbsorption(double* nu, double**** Inu, double**** Anu, int Nnu, double* Ee, double**** Fe, int Np, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double concentration, double Bfactor, double rfactor){
#pragma omp parallel for shared(nu, Inu, Anu, Ee, Fe, Np, Nd, Bn, sintheta, thetaIndex, concentrations, concentration, Bfactor, rfactor)	
	for(int i = 0; i < Nrho; ++i){
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				for(int l = 0; l < Nnu; ++l){
					evaluateEmissivityAndAbsorptionAtNuSimple(nu[l], Inu[l][i][j][k], Anu[l][i][j][k], Ee, Fe[i][j][k], Np, sintheta[i][j][k], Bfactor*Bn[i][j][k]/rfactor, concentration*concentrations[i][j][k]/(rfactor*rfactor));
				}
			}
		}
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



void evaluateVolumeAndLength(double** area, double** length, double rmax, double* Rho, double* Phi, const double& fractionSize){
	double rmin = rmax*(1.0 - fractionSize);
	double dr = rmax/Nrho;
	for(int i = 0; i < Nrho; ++i){
		double minRho = dr*i;
		double maxRho = dr*(i + 1);
		
		double s = 0;
		double l = 0;

		s = 0.5*dphi*(maxRho*maxRho - minRho*minRho);
		//todo rho or minrho
		if(minRho < rmin){
			l = sqrt(rmax*rmax - minRho*minRho) - sqrt(rmin*rmin - minRho*minRho);
		} else {
			l = sqrt(rmax*rmax - minRho*minRho);
		}

		for(int j = 0; j < Nphi; ++j){
			/*length[i][j] = l;
			double phi = Phi[j];
			if(phi > pi/2){
				phi = pi - phi;
			}
			double phi1 = phi - dphi/2;
			double phi2 = phi + dphi/2;

			double l1 = maxRho*cos(phi1);
			double l2 = minRho*cos(phi1);
			double l3 = maxRho*cos(phi2);
			double l4 = minRho*cos(phi2);

			if(l1 <= l0){
				area[i][j] = s;
			} else if(l4 >= l0){
				area[i][j] = 0;
			} else if((l2 < l0) && (l3 < l0)){
				double h = sqrt(maxRho*maxRho - l0*l0);
				double phi3 = atan2(h, l0);
				//todo what if infinity?
				area[i][j] = 0.5*l0*h - 0.5*l0*l0*tan(phi1) - 0.5*minRho*minRho*(phi3 - phi1);
			} else if(l2 < l0){
				area[i][j] = 0.5*l0*l0*tan(phi2) - 0.5*l0*l0*tan(phi1) - 0.5*minRho*minRho*dphi;
			} else if(l3 < l0){
				double h3 = sqrt(maxRho*maxRho - l0*l0);
				double h4 = sqrt(minRho*minRho - l0*l0);
				double phi3 = atan2(h3, l0);
				double phi4 = atan2(h4, l0);
				area[i][j] = 0.5*maxRho*maxRho*(phi2 - phi3) + 0.5*rmin*rmax*sin(phi3 - phi4) - 0.5*minRho*minRho*(phi2 - phi4);
			} else {
				double l2 = l0/cos(phi2);
				double h = sqrt(minRho*minRho - l0*l0);
				double phi3 = atan2(h, l0);
				area[i][j] = 0.5*l2*rmin*sin(phi2 - phi3) - 0.5*minRho*minRho*(phi2 - phi3);
			}
			if(area[i][j] < 0){
				printf("area < 0\n");
			}*/
			length[i][j] = (4.0/3.0)*rmax*fractionSize;
			area[i][j] = s;
		}
	}
}

void evaluateVolumeAndLength(double*** area, double*** length, double rmax, double fractionSize){
	double rmin = rmax*(1.0 - fractionSize);
	double dr = rmax/Nrho;
	double dz = 4.0*rmax/(3.0*Nz);
	for(int i = 0; i < Nrho; ++i){
		double minRho = dr*i;
		double maxRho = dr*(i + 1);
		double s = 0.5*dphi*(maxRho*maxRho - minRho*minRho);
		for(int j = 0; j < Nphi; ++j){
			for(int k = 0; k < Nz; ++k){
				double z = dz*(k+1);

				double l = 0;


				//todo rho or minrho
				if(z < 4.0*rmin/3.0){
					l = 0;
				} else if (z - dz < 4.0*rmin/3.0) {
					l = z - 4.0*rmin/3.0;
				} else {
					l = dz;
				}

				length[i][j][k] = l;
				area[i][j][k] = s;
			}
		}
	}
}

void evaluateSpectrum(double* nu, double* I, double**** Inu, double**** Anu, double*** area, double*** length, int Nnu, double rfactor){
	for(int l = 0; l < Nnu; ++l){
		I[l] = 0;
		for(int i = 0; i < Nrho; ++i){
			for(int j = 0; j < Nphi; ++j){
				double localInu = 0;
				for(int k = 0; k < Nz; ++k){
					double I0 = localInu;
					double Q = Inu[l][i][j][k]*area[i][j][k]*rfactor*rfactor;
					double tau = Anu[l][i][j][k]*length[i][j][k]*rfactor;
					double S = Q/Anu[l][i][j][k];

					localInu = S + (I0 - S)*exp(-tau);
				}
				I[l] = I[l] + localInu;
			}
		}
		I[l] = I[l]*1E26/(distance*distance);
	}
}

//spherical
void evaluateSpectrumSpherical(double* nu, double* I, double**** Inu, double**** Anu, double rmax, int Nnu, double rfactor, double fractionLength){
	int tempNr = 100;
	double tempRmax = rmax*rfactor;
	double tempRmin = (1.0 - fractionLength)*tempRmax;
	double tempdr = tempRmax/tempNr;
	double dphi = 2*pi/Nphi;
	double dz = 2*rmax*rfactor/Nz;

	for(int l = 0; l < Nnu; ++l){
		I[l] = 0;
	}

	double drho = tempRmax/Nrho;

	double length[Nz];

	for(int i = 0; i < tempNr; ++i){
		double r = (i + 0.5)*tempdr;
		double s = 0.5*dphi*(2*i + 1)*tempdr*tempdr;
		double z1 = -sqrt(tempRmax*tempRmax -r*r);
		double z2 = 0;
		if(tempRmin > r){
			z2 = -sqrt(tempRmin*tempRmin - r*r);
		}
		double z3 = -z2;
		double z4 = -z1;

		int rhoindex = floor(r/drho);
		if(rhoindex >= Nrho){
			printf("rhoindex > Nrho\n");
			rhoindex = Nrho-1;
		}

		for(int k = 0; k < Nz; ++k){
			double z = - tempRmax + (k +0.5)*dz;
			double minz  = - tempRmax + k*dz;
			double maxz = -tempRmax + (k+1)*dz;
			//length
			length[k] = 0;
			if(z < 0){
				if(z1 > maxz){
					length[k] = 0;
				} else if (z2 < minz){
					length[k] = 0;
				} else {
					double lowz = max(minz, z1);
					double topz = min(maxz, z2);
					length[k] = topz - lowz;
				}
			} else {
				if(z4 < minz){
					length[k] = 0;
				} else if(z3 > maxz){
					length[k] = 0;
				} else {
					double lowz = max(minz, z3);
					double topz = min(maxz, z4);
					length[k] = topz - lowz;
				}
			}
			if(i == 0){
				printf("length = %g\n", length[k]);
			}
		}

		for(int j = 0; j < Nphi; ++j){
			for(int l = 0; l < Nnu; ++l){
				double localI = 0;
				for(int k = 0; k < Nz; ++k){				
					double I0 = localI;
					if(length[k] > 0){
						double Q = Inu[l][rhoindex][j][k]*s;
						double tau = Anu[l][rhoindex][j][k]*length[k];
						double S = 0;
						if(Q > 0){
							S = Q/Anu[l][rhoindex][j][k];
						}
						if(fabs(tau) < 1E-15){
							localI = I0*(1.0 - tau) + S*tau;
						} else {
							localI = S + (I0 - S)*exp(-tau);
						}
					}
				}
				I[l] += localI;
			}
			
		}
	}

	//delete[] length;

	for(int l = 0; l < Nnu; ++l){
		I[l] = I[l]*1E26/(distance*distance);
	}
}

void evaluateSpectrumSphericalAtNu(double nu, double& I, double*** Inu, double*** Anu, double rmax, double rfactor, double fractionLength, double d){
	int tempNr = 100;
	double tempRmax = rmax*rfactor;
	double tempRmin = (1.0 - fractionLength)*tempRmax;
	double tempdr = tempRmax/tempNr;
	double dphi = 2*pi/Nphi;
	double dz = 2*rmax*rfactor/Nz;

	I = 0;

	double drho = tempRmax/Nrho;

	double length[Nz];

	for(int i = 0; i < tempNr; ++i){
		double r = (i + 0.5)*tempdr;
		double s = 0.5*dphi*(2*i + 1)*tempdr*tempdr;
		double z1 = -sqrt(tempRmax*tempRmax -r*r);
		double z2 = 0;
		if(tempRmin > r){
			z2 = -sqrt(tempRmin*tempRmin - r*r);
		}
		double z3 = -z2;
		double z4 = -z1;

		int rhoindex = floor(r/drho);
		if(rhoindex >= Nrho){
			printf("rhoindex > Nrho\n");
			rhoindex = Nrho-1;
		}

		for(int k = 0; k < Nz; ++k){
			double z = - tempRmax + (k +0.5)*dz;
			double minz  = - tempRmax + k*dz;
			double maxz = -tempRmax + (k+1)*dz;
			//length
			length[k] = 0;
			if(z < 0){
				if(z1 > maxz){
					length[k] = 0;
				} else if (z2 < minz){
					length[k] = 0;
				} else {
					double lowz = max(minz, z1);
					double topz = min(maxz, z2);
					length[k] = topz - lowz;
				}
			} else {
				if(z4 < minz){
					length[k] = 0;
				} else if(z3 > maxz){
					length[k] = 0;
				} else {
					double lowz = max(minz, z3);
					double topz = min(maxz, z4);
					length[k] = topz - lowz;
				}
			}
		}

		for(int j = 0; j < Nphi; ++j){
				double localI = 0;
				for(int k = 0; k < Nz; ++k){				
					double I0 = localI;
					if(length[k] > 0){
						double Q = Inu[rhoindex][j][k]*s;
						double tau = Anu[rhoindex][j][k]*length[k];
						double S = 0;
						if(Q > 0){
							S = Q/Anu[rhoindex][j][k];
						}
						if(fabs(tau) < 1E-15){
							localI = I0*(1.0 - tau) + S*tau;
						} else {
							localI = S + (I0 - S)*exp(-tau);
						}
					}
				}
				I += localI;
			}
			
		}

	//delete[] length;

	I = I*1E26/(d*d);
}

//flat
void evaluateSpectrumFlat(double* nu, double* I, double**** Inu, double**** Anu, double rmax, int Nnu, double rfactor, double fractionLength){
	double tempRmax = rmax*rfactor;
	double tempRmin = (1.0 - fractionLength)*tempRmax;
	double tempdr = tempRmax/Nrho;
	double dphi = 2*pi/Nphi;
	double dz = (4.0/3.0)*tempRmax/Nz;
	double zmin = (4.0/3.0)*(1.0 - fractionLength)*tempRmax;

	for(int l = 0; l < Nnu; ++l){
		I[l] = 0;
	}

	double drho = tempRmax/Nrho;

	double length[Nz];

	for(int k = 0; k < Nz; ++k){
		double z1  = k*dz;
		double z2 = (k+1)*dz;
		//length
		length[k] = 0;
		if(z2 < zmin){
			length[k] = 0;
		} else if(z1 > zmin){
			length[k] = dz;
		} else {
			length[k] = z2 - zmin;
		}
	}

	for(int i = 0; i < Nrho; ++i){
		double s = 0.5*dphi*(2*i + 1)*tempdr*tempdr;
		for(int j = 0; j < Nphi; ++j){
			for(int l = 0; l < Nnu; ++l){
				double localI = 0;
				for(int k = 0; k < Nz; ++k){				
					double I0 = localI;
					if(length[k] > 0){
						double Q = Inu[l][i][j][k]*s;
						if(Inu[l][i][j][k] != Inu[l][i][j][k]){
							printf("Inu[l][i][j][k] = NaN\n");
							exit(0);
						}
						double tau = Anu[l][i][j][k]*length[k];
						if(Anu[l][i][j][k] != Anu[l][i][j][k]){
							printf("Anu[l][i][j][k] = NaN\n");
							exit(0);
						}
						double S = 0;
						if(Q > 0){
							S = Q/Anu[l][i][j][k];
						}

						if(fabs(tau) < 1E-15){
							localI = I0*(1.0 - tau) + S*tau;
						} else {
							localI = S + (I0 - S)*exp(-tau);
						}
						if(localI != localI){
							printf("Anu = %g Inu = %g Q = %g s = %g length = %g tau = %g I0 = %g\n", Anu[l][i][j][k], Inu[l][i][j][k], Q, s, length[k],tau, I0);
							printf("localI = NaN\n");
							exit(0);
						}
					}
				}
				I[l] += localI;
			}
			
		}
	}

	//delete[] length;

	for(int l = 0; l < Nnu; ++l){
		I[l] = I[l]*1E26/(distance*distance);
	}
}

void evaluateSpectrumFlatAtNu(double nu, double& I, double*** Inu, double*** Anu, double rmax, double rfactor, double fractionLength, double d){
	double tempRmax = rmax*rfactor;
	double tempRmin = (1.0 - fractionLength)*tempRmax;
	double tempdr = tempRmax/Nrho;
	double dphi = 2*pi/Nphi;
	double dz = (4.0/3.0)*tempRmax/Nz;
	double zmin = (4.0/3.0)*(1.0 - fractionLength)*tempRmax;

	I = 0;

	double drho = tempRmax/Nrho;

	double length[Nz];

	for(int k = 0; k < Nz; ++k){
		double z1  = k*dz;
		double z2 = (k+1)*dz;
		//length
		length[k] = 0;
		if(z2 < zmin){
			length[k] = 0;
		} else if(z1 > zmin){
			length[k] = dz;
		} else {
			length[k] = z2 - zmin;
		}
	}

	for(int i = 0; i < Nrho; ++i){
		double s = 0.5*dphi*(2*i + 1)*tempdr*tempdr;
		for(int j = 0; j < Nphi; ++j){
			double localI = 0;
			for(int k = 0; k < Nz; ++k){				
				double I0 = localI;
				if(length[k] > 0){
					double Q = Inu[i][j][k]*s;
					if(Inu[i][j][k] != Inu[i][j][k]){
						printf("Inu[i][j][k] = NaN\n");
						exit(0);
					}
					double tau = Anu[i][j][k]*length[k];
					if(Anu[i][j][k] != Anu[i][j][k]){
						printf("Anu[i][j][k] = NaN\n");
						exit(0);
					}
					double S = 0;
					if(Q > 0){
						S = Q/Anu[i][j][k];
					}

					if(fabs(tau) < 1E-15){
						localI = I0*(1.0 - tau) + S*tau;
					} else {
						localI = S + (I0 - S)*exp(-tau);
					}
					if(localI != localI){
						printf("Anu = %g Inu = %g Q = %g s = %g length = %g tau = %g I0 = %g\n", Anu[i][j][k], Inu[i][j][k], Q, s, length[k],tau, I0);
						printf("localI = NaN\n");
						exit(0);
					}
				}
			}
			I += localI;
		}
			
	}

	//delete[] length;

	I = I*1E26/(d*d);
}

void evaluateSpectrum(double* nu, double* totalInu, double*** Inu, double*** Anu, double** area, double** length, int Nnu, double rmax, double* Rho, double* Phi){
	for(int l = 0; l < Nnu; ++l){
		totalInu[l] = 0;
		int currentThetaIndex = 0;
		double dr = rmax/Nrho;
		for(int i = 0; i < Nrho; ++i){
			double minRho = dr*i;
			double maxRho = dr*(i+1);
			//todo while, -1 or not?
			while((minRho > rmax*sinThetaValue[currentThetaIndex]) && (currentThetaIndex < Ntheta/2 - 1)){
				currentThetaIndex++;
			}
			for(int j = 0; j < Nphi; ++j){
				double localInu = 0;

					//todo find what Inu
					double Q = Inu[Ntheta - 1 - currentThetaIndex][j][l]*area[i][j];
					double A = Anu[Ntheta - 1 - currentThetaIndex][j][l];

					double tau = A*length[i][j];

					if(tau > 0.01){
						localInu = (Q/A)*(1.0 - exp(-tau));
					} else {
						localInu = Q*length[i][j] - A*Q*length[i][j]*length[i][j]/2.0;
					}

					//Q = Inu[currentThetaIndex][j][l]*area[i][j];
					//A = Anu[currentThetaIndex][j][l];
					//double B = localInu*A/Q - 1.0;

					//tau = A*length[i][j];

					//if(tau > 0.01){
					//	localInu = (Q/A)*(1 + B*exp(-tau));
					//} else {
					//	localInu = localInu - (localInu - Q/A)*tau + (localInu - Q/A)*tau*tau/2.0;
					//}

				totalInu[l] = totalInu[l] + localInu;
			}
		}
		totalInu[l] = totalInu[l]*1E26/(distance*distance);
	}
}

void evaluateSpectrumFlatSimple(double* nu, double* totalInu, double* Inu, double* Anu, int Nnu, double rmax, double fraction){
	double area = pi*rmax*rmax;
	double length = 4.0*fraction*rmax/3.0;
	for(int l = 0; l < Nnu; ++l){
		totalInu[l] = 0;


		double localInu = 0;

					//todo find what Inu
		double Q = Inu[l]*area;
		double A = Anu[l];

		double tau = A*length;

		if(tau > 0.01){
			localInu = (Q/A)*(1.0 - exp(-tau));
		} else {
			localInu = Q*length - A*Q*length*length/2.0;
		}

					//Q = Inu[currentThetaIndex][j][l]*area[i][j];
					//A = Anu[currentThetaIndex][j][l];
					//double B = localInu*A/Q - 1.0;

					//tau = A*length[i][j];

					//if(tau > 0.01){
					//	localInu = (Q/A)*(1 + B*exp(-tau));
					//} else {
					//	localInu = localInu - (localInu - Q/A)*tau + (localInu - Q/A)*tau*tau/2.0;
					//}

		totalInu[l] = totalInu[l] + localInu;

		totalInu[l] = totalInu[l]*1E26/(distance*distance);
	}
}

void evaluateSpectrumAtNuSimple(double nu, double& totalInu, double Inu, double Anu, double rmax, double fraction){
	totalInu = 0;
	int currentThetaIndex = 0;

	double area = pi*rmax*rmax;
	double length = 4.0*fraction*rmax/3.0;

	//todo find what Inu
	double Q = Inu*area;
	double A = Anu;
	double localInu = 0;

	double tau = A*length;

	if(tau > 0.01){
		localInu = (Q/A)*(1.0 - exp(-tau));
	} else {
		localInu = Q*length - A*Q*length*length/2.0;
	}

	//Q = Inu[currentThetaIndex][j][l]*area[i][j];
	//A = Anu[currentThetaIndex][j][l];
	//double B = localInu*A/Q - 1.0;

	//tau = A*length[i][j];

	//if(tau > 0.01){
	//	localInu = (Q/A)*(1 + B*exp(-tau));
	//} else {
	//	localInu = localInu - (localInu - Q/A)*tau + (localInu - Q/A)*tau*tau/2.0;
	//}

	totalInu = totalInu + localInu;

	totalInu = totalInu*1E26/(distance*distance);
}

//spherical
void evaluateImageSpherical(double*** image, double* nu, double**** Inu, double**** Anu, double rmax, int Nnu, double rfactor, double fractionLength){
	double tempRmax = rmax*rfactor;
	double tempRmin = (1.0 - fractionLength)*tempRmax;
	double tempdr = tempRmax/Nrho;
	double dphi = 2*pi/Nphi;
	double dz = 2*rmax*rfactor/Nz;

	for(int i = 0; i < Nrho; ++i) {
		for(int j = 0; j < Nphi; ++j){
			for(int l = 0; l < Nnu; ++l){
				image[i][j][l] = 0;
			}
		}
	}


	double length[Nz];

	for(int i = 0; i < Nrho; ++i){
		double r = (i + 0.5)*tempdr;
		double s = 0.5*dphi*(2*i + 1)*tempdr*tempdr;
		double z1 = -sqrt(tempRmax*tempRmax -r*r);
		double z2 = 0;
		if(tempRmin > r){
			z2 = -sqrt(tempRmin*tempRmin - r*r);
		}
		double z3 = -z2;
		double z4 = -z1;


		for(int k = 0; k < Nz; ++k){
			double z = - tempRmax + (k +0.5)*dz;
			double minz  = - tempRmax + k*dz;
			double maxz = -tempRmax + (k+1)*dz;
			//length
			length[k] = 0;
			if(z < 0){
				if(z1 > maxz){
					length[k] = 0;
				} else if (z2 < minz){
					length[k] = 0;
				} else {
					double lowz = max(minz, z1);
					double topz = min(maxz, z2);
					length[k] = topz - lowz;
				}
			} else {
				if(z4 < minz){
					length[k] = 0;
				} else if(z3 > maxz){
					length[k] = 0;
				} else {
					double lowz = max(minz, z3);
					double topz = min(maxz, z4);
					length[k] = topz - lowz;
				}
			}
		}

		for(int j = 0; j < Nphi; ++j){
			for(int l = 0; l < Nnu; ++l){
				double localI = 0;
				for(int k = 0; k < Nz; ++k){				
					double I0 = localI;
					if(length[k] > 0){
						double Q = Inu[l][i][j][k]*s;
						double tau = Anu[l][i][j][k]*length[k];
						double S = 0;
						if(Q > 0){
							S = Q/Anu[l][i][j][k];
						}

						localI = S + (I0 - S)*exp(-tau);
					}
				}
				image[i][j][l] = localI;
			}
			
		}
	}

	//delete[] length;
	for(int i = 0; i < Nrho; ++i) {
		for(int j = 0; j < Nphi; ++j){
			for(int l = 0; l < Nnu; ++l){
				image[i][j][l] = image[i][j][l]*1E26/(distance*distance);
			}
		}
	}
}

//flat
void evaluateImageFlat(double*** image, double* nu, double**** Inu, double**** Anu, double rmax, int Nnu, double rfactor, double fractionLength){
	double tempRmax = rmax*rfactor;
	double tempRmin = (1.0 - fractionLength)*tempRmax;
	double tempdr = tempRmax/Nrho;
	double dphi = 2*pi/Nphi;
	double dz = (4.0/3.0)*tempRmax/Nz;
	double zmin = (4.0/3.0)*(1.0 - fractionLength)*tempRmax;

	for(int i = 0; i < Nrho; ++i){
		for(int j = 0; j < Nphi; ++j){
			for(int l = 0; l < Nnu; ++l){
				image[i][j][l] = 0;
			}
		}
	}

	double length[Nz];

	for(int k = 0; k < Nz; ++k){
		double z1  = k*dz;
		double z2 = (k+1)*dz;
		//length
		length[k] = 0;
		if(z2 < zmin){
			length[k] = 0;
		} else if(z1 > zmin){
			length[k] = dz;
		} else {
			length[k] = z2 - zmin;
		}
	}

	for(int i = 0; i < Nrho; ++i){
		double s = 0.5*dphi*(2*i + 1)*tempdr*tempdr;
		for(int j = 0; j < Nphi; ++j){
			for(int l = 0; l < Nnu; ++l){
				double localI = 0;
				for(int k = 0; k < Nz; ++k){				
					double I0 = localI;
					if(length[k] > 0){
						double Q = Inu[l][i][j][k]*s;
						if(Inu[l][i][j][k] != Inu[l][i][j][k]){
							printf("Inu[l][i][j][k] = NaN\n");
							exit(0);
						}
						double tau = Anu[l][i][j][k]*length[k];
						if(Anu[l][i][j][k] != Anu[l][i][j][k]){
							printf("Anu[l][i][j][k] = NaN\n");
							exit(0);
						}
						double S = 0;
						if(Q > 0){
							S = Q/Anu[l][i][j][k];
						}

						localI = S + (I0 - S)*exp(-tau);
						if(localI != localI){
							printf("Anu = %g Inu = %g Q = %g s = %g length = %g tau = %g I0 = %g\n", Anu[l][i][j][k], Inu[l][i][j][k], Q, s, length[k],tau, I0);
							printf("localI = NaN\n");
							exit(0);
						}
					}
				}
				image[i][j][l] += localI;
			}
			
		}
	}

	//delete[] length;

	for(int i = 0; i < Nrho; ++i) {	
		for(int j = 0; j < Nphi; ++j){
			for(int l = 0; l < Nnu; ++l){
				image[i][j][l] = image[i][j][l]*1E26/(distance*distance);
			}
		}	
	}
}

double evaluateNextdFe(double* Ee, double* dFe, double dg, int j, int Np) {
	double gamma = Ee[j] / (massElectron * speed_of_light2);
	double nextGamma = gamma + dg;
	double nextdFe = dFe[j];
	if(j == 0) {
		double tempE = (gamma + dg)*massElectron*speed_of_light2;
		double F = 0;
		double nextF = dFe[j+1]/(Ee[j+1]-Ee[j]);
		nextdFe = F*(Ee[j+1] - tempE) + nextF*(tempE - Ee[j]);
	} else if(j == Np-1){
		double tempE = (gamma + dg)*massElectron*speed_of_light2;
		double dE = Ee[j] - Ee[j-1];
		double F = dFe[j]/dE;
		double prevF = dFe[j-1]/(Ee[j-1]-Ee[j-2]);
		nextdFe = prevF*(Ee[j] - tempE) + F*(tempE - Ee[j-1]);
	} else {
		double tempE = (gamma + dg)*massElectron*speed_of_light2;
		//printf("tempe = %g\n", tempE);
		//printf("E[j] = %g\n", Ee[j]);
		//printf("E[j+1] = %g\n", Ee[j+1]);
		double dE = Ee[j] - Ee[j-1];
		double F = dFe[j]/dE;
		//printf("F = %lf\n", F);
		double nextF = dFe[j+1]/(Ee[j+1]-Ee[j]);
		//printf("nextF = %lf\n", nextF);
		nextdFe = (F*(Ee[j+1] - tempE) + nextF*(tempE - Ee[j]))/(Ee[j+1]-Ee[j]);
		nextdFe *= dE;
	}

	return nextdFe;
}

void evaluateEmissivityAndAbsorptionAtNuSimple(double nu, double& Inu, double& Anu, double* Ee, double* dFe, int Np, double sinhi, double B, double concentration){
	//Anu from ghiselini simple
	Inu = 0;
	Anu = 0;

	if(sinhi == 0.0){
		return;
	}

	double coef = concentration * emissivityCoef;
	//double coefAbsorb = concentration * absorpCoef/B;
	//todo what if < 0?
	double coshi = sqrt(1.0 - sinhi*sinhi);


	double oldA = 0;
	for (int j = 1; j < Np; ++j) {
		if(dFe[j] > 0){
			double nuc = criticalNu(Ee[j], sinhi, B);
			double gamma = Ee[j] / (massElectron * speed_of_light2);
			double gamma4 = gamma*gamma*gamma*gamma;
			double sigmaCoef = concentration*dFe[j]/(B*gamma4);

			Inu = Inu + coef * dFe[j] * B * sinhi * evaluateMcDonaldIntegral(nu / nuc);
			if(Inu < 0){
				printf("Inu < 0\n");
				printf("dFe[j] = %g\n", dFe[j]);
				exit(0);
			}
			////todo! F(g +dg)
			double tempP = gamma*gamma*coef * B*dFe[j] * sinhi * evaluateMcDonaldIntegral(nu / nuc);
			double dg = 0.1*(Ee[j] - Ee[j-1])/ (massElectron * speed_of_light2);

			double tempGamma = gamma + dg;
			double tempnuc = criticalNu(massElectron*speed_of_light2*tempGamma, sinhi, B);
			double nextdFe = evaluateNextdFe(Ee, dFe, dg, j, Np);
			if(nextdFe < 0){
				printf("nextdFe < 0\n");
			}
			double tempP2 = tempGamma*tempGamma*coef * B*nextdFe * sinhi * evaluateMcDonaldIntegral(nu / tempnuc);
			double Pder = (tempP2 - tempP)/dg;

				

			Anu = Anu + (1.0/(2*massElectron*nu*nu))*Pder/(gamma*gamma);
			/*if(Anu < 0){
				printf("Anu < 0\n");
				printf("dFe[j] = %g\n", dFe[j]);
				//exit(0);
			}*/
			oldA = oldA + (16*pi*pi*electron_charge/(3.0*sqrt(3.0)))*dFe[j]*evaluateMcDonaldFunction5_3(nu/nuc)/(gamma*gamma*gamma*gamma*gamma*B*sinhi);
			if(Inu != Inu){
				printf("Inu NaN\n");
				exit(0);
			}
			if(Anu != Anu){
				printf("Anu Nan\n");
				exit(0);
			}
		}
	}
}

void evaluateEmissivityAndAbsorptionFlatSimple(double* nu, double* Inu, double* Anu, double* Ee, double* Fe, int Np, int Nnu, double sinhi, double B, double concentration){
	for(int i = 0; i < Nnu; ++i){
		Inu[i] = 0;
		Anu[i] = 0;
	}

	if(sinhi == 0.0){
		return;
	}

	for(int i = 0; i < Nnu; ++i){
		evaluateEmissivityAndAbsorptionAtNuSimple(nu[i], Inu[i], Anu[i], Ee, Fe, Np, sinhi, B, concentration);
	}
}