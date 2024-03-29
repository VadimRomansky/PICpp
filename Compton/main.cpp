// Compton.cpp : Defines the entry point for the console application.
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

//transform from one spherical system to rotated. Rotation on phir around z, and then on mur around x' 
void rotationSphericalCoordinates(const double& mur, const double& phir, const double& mu0, const double& phi0, double& mu1, double& phi1) {
	if (mur > 1.0) {
		printf("mu_r = %g > 1.0\n", mur);
		exit(0);
	}
	if (mur < -1.0) {
		printf("mu_r = %g < -1.0\n", mur);
		exit(0);
	}
	if (mu0 > 1.0) {
		printf("mu0 = %g > 1.0\n", mu0);
		exit(0);
	}
	if (mu0 < -1.0) {
		printf("mu0 = %g < -1.0\n", mu0);
		exit(0);
	}
	double sinThetar = sqrt(1.0 - mur * mur);
	double sinTheta0 = sqrt(1.0 - mu0 * mu0);
	double tempPhi = phi0 - phir;
	if (tempPhi < 0) {
		tempPhi = tempPhi + 2 * pi;
	}
	mu1 = mu0 * mur + sinThetar * sinTheta0 * sin(tempPhi);

	if (mu1 > 1.0) {
		printf("mu1 = %g > 1.0\n", mu1);
		exit(0);
	}
	if (mu1 < -1.0) {
		printf("mu1 = %g < -1.0\n", mu1);
		exit(0);
	}

	double sinTheta1 = sqrt(1.0 - mu1 * mu1);
	if (sinTheta1 == 0) {
		phi1 = 0;
		return;
	}
	double coshi = sinTheta0 * cos(tempPhi);
	double cosPhi1 = coshi / sinTheta1;
	if (cosPhi1 > 1.0) {
		printf("cosPhi1 = %g > 1.0\n", cosPhi1);
		exit(0);
	}
	if (cosPhi1 < -1.0) {
		printf("cosPhi1 = %g < -1.0\n", cosPhi1);
		exit(0);
	}
	double scalarMult = sin(tempPhi) * sinTheta0 * mur - mu0 * sinThetar;
	phi1 = acos(cosPhi1);
	if (scalarMult < 0) {
		phi1 = 2 * pi - phi1;
	}
}

void inverseRotationSphericalCoordinates(const double& mur, const double& phir, const double& mu1, const double& phi1, double& mu0, double& phi0) {
	if (mur > 1.0) {
		printf("mu_r = %g > 1.0\n", mur);
		exit(0);
	}
	if (mur < -1.0) {
		printf("mu_r = %g < -1.0\n", mur);
		exit(0);
	}
	if (mu1 > 1.0) {
		printf("mu1 = %g > 1.0\n", mu1);
		exit(0);
	}
	if (mu1 < -1.0) {
		printf("mu1 = %g < -1.0\n", mu1);
		exit(0);
	}
	double sinThetar = sqrt(1.0 - mur * mur);
	double sinTheta1 = sqrt(1.0 - mu1 * mu1);

	mu0 = mu1 * mur + sinThetar * sinTheta1* cos(phi1+pi/2);
	if (mu0 > 1.0) {
		printf("mu0 = %g > 1.0\n", mu0);
		exit(0);
	}
	if (mu0 < -1.0) {
		printf("mu0 = %g < -1.0\n", mu0);
		exit(0);
	}
	double sinTheta0 = sqrt(1.0 - mu0 * mu0);
	if (sinTheta0 == 0) {
		phi0 = 0;
		return;
	}

	double coshi = sinTheta1 * cos(phi1);

	double cosTempPhi = coshi / sinTheta0;
	if (cosTempPhi > 1.0) {
		printf("cosTempPhi = %g > 1.0\n", cosTempPhi);
		exit(0);
	}
	if (cosTempPhi < -1.0) {
		printf("cosTempPhi = %g < -1.0\n", cosTempPhi);
		exit(0);
	}
	double scalarMult = sin(phi1) * sinTheta1 * mur + mu1 * sinThetar;
	double tempPhi = acos(cosTempPhi);
	if (scalarMult < 0) {
		tempPhi = 2 * pi - tempPhi;
	}
	phi0 = tempPhi + phir;
	if (phi0 > 2 * pi) {
		phi0 = phi0 - 2 * pi;
	}
	if (phi0 < 0) {
		phi0 = phi0 + 2 * pi;
	}
}

double evaluatePhotonDistribution(const double& energy, const double& mu, const double& phi, int Np, double* Eph, double* Fph){
	if(energy <= Eph[1]){

		return 0;
	} else if(energy >= Eph[Np-1]){
		return 0;
	}
	else {
		int currentI = 0;
		int nextI = 1;
		for (int i = 1; i < Np; ++i) {
			currentI = i - 1;
			nextI = i;
			if (Eph[i] > energy) {
				break;
			}
		}
		double result;
		if (Fph[currentI] <= 0 || Fph[nextI] <= 0) {
			result = (Fph[currentI] * (Eph[nextI] - energy) + Fph[nextI] * (energy - Eph[currentI])) / (Eph[nextI] - Eph[currentI]);
		} else {
			result = Fph[currentI] * exp(log(Fph[nextI] / Fph[currentI]) * ((energy - Eph[currentI]) / (Eph[nextI] - Eph[currentI])));
		}
		if(result != result){
			printf("result = NaN\n");
			exit(0);
		}
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

				double cosTheta1 = cosThetaValue[j];
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
				thetaIndex[i][j][k] = 3;
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
	const int Nphot = 2000;
	const int Nnu = 200;

	const double Bfactor = 0.29;
	const double epsilonB = 0.0012;
	//const double electronConcentration = 10*Bfactor*Bfactor/(4*pi*massProtonReal*speed_of_light2*epsilonB);
	const double electronConcentration = 150;
	const double photonConcentration = 1.0;
	//double rmax = 0.3 * speed_of_light * 357 * 24 * 3600;
	double rmax = 1.3E17;


	const int Nphi = 4;
	const int NthetaFinal = 20;
	const int NthetaElectron = 20;
	const int NthetaInitial = 20;
	const int NthetaSpace = 4;
	double cosThetaLeftFinal[NthetaFinal];
	double cosThetaFinal[NthetaFinal];
	double ThetaFinal[NthetaFinal];
	double cosThetaLeftInitial[NthetaInitial];
	double cosThetaInitial[NthetaInitial];
	double cosThetaLeftElectron[NthetaElectron];
	double cosThetaElectron[NthetaElectron];
	double ThetaInitial[NthetaInitial];
	double ThetaElectron[NthetaElectron];
	double cosThetaSpace[NthetaSpace];
	double phi[Nphi];
	double sinPhiValue[Nphi];
	double cosPhiValue[Nphi];

	double dcosThetaFinal[NthetaFinal];
	double dcosThetaInitial[NthetaInitial];
	double dcosThetaElectron[NthetaElectron];
	double dcosThetaSpace = 2.0/NthetaSpace;
	double dphi = 2.0*pi/Nphi;
	double gammamax = 10000;
	double thetamin = 0.02/gammamax;
	double dlogtheta = log(pi/thetamin)/(NthetaFinal-1);
	ThetaFinal[0] = 0;
	cosThetaLeftFinal[0] = 1.0;
	cosThetaFinal[0] = cos(thetamin);
	for(int i = 1; i < NthetaFinal; ++i){
		ThetaFinal[i] = thetamin*exp(dlogtheta*(i));
		cosThetaFinal[i] = cos(thetamin*exp(dlogtheta*(i)));
		cosThetaLeftFinal[i] = (cosThetaFinal[i] + cosThetaFinal[i - 1]) / 2.0;
	}
	for(int i = 0; i < NthetaFinal-1; ++i){
		dcosThetaFinal[i] = -(cosThetaLeftFinal[i+1] - cosThetaLeftFinal[i]);
		//cosThetaFinal[i] = (cosThetaLeftFinal[i] + cosThetaLeftFinal[i+1])/2.0;
	}
	dcosThetaFinal[NthetaFinal - 1] = 1.0 + cosThetaLeftFinal[NthetaFinal - 1];
	cosThetaFinal[NthetaFinal - 1] = -1.0;

	ThetaInitial[0] = 0;
	cosThetaLeftInitial[0] = 1.0;
	cosThetaInitial[0] = cos(thetamin);
	for (int i = 1; i < NthetaInitial; ++i) {
		ThetaInitial[i] = thetamin * exp(dlogtheta * (i));
		cosThetaInitial[i] = cos(thetamin * exp(dlogtheta * (i)));
		cosThetaLeftInitial[i] = (cosThetaInitial[i] + cosThetaInitial[i - 1]) / 2.0;
	}
	for (int i = 0; i < NthetaInitial - 1; ++i) {
		dcosThetaInitial[i] = -(cosThetaLeftInitial[i + 1] - cosThetaLeftInitial[i]);
		//cosThetaInitial[i] = (cosThetaLeftInitial[i] + cosThetaLeftInitial[i+1])/2.0;
	}
	dcosThetaInitial[NthetaInitial - 1] = 1.0 + cosThetaLeftInitial[NthetaInitial - 1];
	cosThetaFinal[NthetaInitial - 1] = -1.0;

	ThetaElectron[0] = 0;
	cosThetaLeftElectron[0] = 1.0;
	cosThetaElectron[0] = cos(thetamin);
	for (int i = 1; i < NthetaElectron; ++i) {
		ThetaElectron[i] = thetamin * exp(dlogtheta * (i));
		cosThetaElectron[i] = cos(thetamin * exp(dlogtheta * (i)));
		cosThetaLeftElectron[i] = (cosThetaElectron[i] + cosThetaElectron[i - 1]) / 2.0;
	}
	for (int i = 0; i < NthetaElectron - 1; ++i) {
		dcosThetaElectron[i] = -(cosThetaLeftElectron[i + 1] - cosThetaLeftElectron[i]);
		//cosThetaElectron[i] = (cosThetaLeftElectron[i] + cosThetaLeftElectron[i+1])/2.0;
	}
	dcosThetaElectron[NthetaElectron - 1] = 1.0 + cosThetaLeftElectron[NthetaElectron - 1];
	cosThetaElectron[NthetaElectron - 1] = -1.0;

	for(int i = 0; i < Nphi; ++i){
		phi[i] = (i + 0.5)*dphi;
		sinPhiValue[i] = sin(phi[i]);
		cosPhiValue[i] = cos(phi[i]);
	}
	for(int i = 0; i < NthetaSpace; ++i){
		cosThetaSpace[i] = 1.0 - (i + 0.5)*dcosThetaSpace;
	}
	/*for (int i = 0; i < NthetaInitial; ++i) {
		cosThetaInitial[i] = 1.0 - (i + 0.5)*dcosThetaInitial;
		cosThetaLeftInitial[i] = 1.0 - i*dcosThetaInitial;
	}*/
	/*
	//for check rotation transform
	double mur = 0.5;
	double phir = 0.6;
	double mu0 = 0.7;
	double phi0 = 0.9;
	double mu1, phi1, mu2, phi2;
	rotationSphericalCoordinates(mur, phir, mu0, phi0, mu1, phi1);
	inverseRotationSphericalCoordinates(mur, phir, mu1, phi1, mu2, phi2);*/

	const double intx2plank = 2.4042;
	const double intx3plank = pi*pi*pi*pi/15;
	//energy density = 8*pi*(kT)^4/(hc)^3 * intx3plank
	//a = real density/energy density
	
	double L = 1.0E44;

	double Tphotons1 = 2.7;
	double Tphotons2 = 20;
	double Tphotons3 = 3000;
	double Tphotons4 = 4000;
	double Tphotons5 = 7500;
	double a1 = 1.0;
	//double a1 = 15*L*cube(hplank*speed_of_light)/(32*pi*pi*pi*pi*pi*pi*speed_of_light*rmax*rmax*pow(kBoltzman*Tphotons1,4));
	//a1 = 0;
	double a2 = 4E-4;
	a2 = 0;
	double a3 = 4.0E-13;
	a3 = 0;
	double a4 = 1.65E-13;
	a4 = 0;
	double a5 = 1E-14;
	a5 = 0;
	double* Fph = new double[Nphot];
	double* dFph = new double[Nphot];
	double* Eph = new double[Nphot];

	double Ephmin = 0.01*kBoltzman*Tphotons1;
	double Ephmax = kBoltzman*Tphotons1*100;
	Ephmax = kBoltzman * Tphotons5 * 100;

	FILE* logFile = fopen("log.dat", "w");
	fclose(logFile);

	printLog("start\n");

	double factor = pow(Ephmax/Ephmin, 1.0/(Nphot - 1));
	Eph[0] = Ephmin;
	Fph[0] = 0;
	dFph[0] = 0;
	//todo check and normalize
	//for intergalactic field
	for (int i = 1; i < Nphot; ++i) {
		//todo pi or not to pi?
		Fph[i] = 0;
		dFph[i] = 0;
		Eph[i] = Eph[i-1]*factor;
		double theta = Eph[i]/(kBoltzman*Tphotons1);
		Fph[i] +=  a1*(2*Eph[i]*Eph[i]/cube(hplank*speed_of_light))/(exp(theta) - 1.0);
		theta = Eph[i]/(kBoltzman*Tphotons2);
		Fph[i] +=  a2*(2*Eph[i]*Eph[i]/cube(hplank*speed_of_light))/(exp(theta) - 1.0);
		theta = Eph[i]/(kBoltzman*Tphotons3);
		Fph[i] +=  a3*(2*Eph[i]*Eph[i]/cube(hplank*speed_of_light))/(exp(theta) - 1.0);
		theta = Eph[i] / (kBoltzman * Tphotons4);
		Fph[i] += a4 * (2 * Eph[i] * Eph[i] / cube(hplank * speed_of_light)) / (exp(theta) - 1.0);
		theta = Eph[i] / (kBoltzman * Tphotons5);
		Fph[i] += a5 * (2 * Eph[i] * Eph[i] / cube(hplank * speed_of_light)) / (exp(theta) - 1.0);
		dFph[i] = Fph[i]*(Eph[i] - Eph[i-1]);
	}

	//for cluster
	/*double R = 1E18;
	double numax = 1E15;
	double westerlundD = 4000 * 3E18;
	double westerlundNuFNu = 1E14;
	double westerlundL = (westerlundNuFNu / numax) * 4 * pi * westerlundD * westerlundD;
	double realenergydensity = westerlundL / (4 * pi * speed_of_light * R * R);
	double temperature = hplank*numax/kBoltzman;
	double planckenergydensity = 8 * pi * pow(kBoltzman * temperature, 4) / pow(hplank * speed_of_light, 3) * intx3plank;
	double a = realenergydensity / planckenergydensity;
	if( a > 1){
	    printf("a > 1 in planck spectrum\n");
	}

	for (int i = 1; i < Nphot; ++i) {
		Eph[i] = Eph[i - 1] * factor;
		double theta = Eph[i] / (kBoltzman * temperature);
		Fph[i] = a * (2 * Eph[i] * Eph[i] / cube(hplank * speed_of_light)) / (exp(theta) - 1.0);
		dFph[i] = Fph[i] * (Eph[i] - Eph[i - 1]);
	}*/

	//for star
	/*double R = rmax;
	double starL = 300000*4E33;
	double realenergydensity = starL / (4 * pi * speed_of_light * R * R);
	double temperature = 50000;
	double planckenergydensity = 8 * pi * pow(kBoltzman * temperature, 4) / pow(hplank * speed_of_light, 3) * intx3plank;
	double a = realenergydensity / planckenergydensity;
	if (a > 1) {
		printf("a > 1 in planck spectrum\n");
	}

	for (int i = 1; i < Nphot; ++i) {
		Eph[i] = Eph[i - 1] * factor;
		double theta = Eph[i] / (kBoltzman * temperature);
		Fph[i] = a * (2 * Eph[i] * Eph[i] / cube(hplank * speed_of_light)) / (exp(theta) - 1.0);
		dFph[i] = Fph[i] * (Eph[i] - Eph[i - 1]);
	}*/

	FILE* photons = fopen("photons.dat","w");
	for(int i = 0; i < Nphot; ++i){
		fprintf(photons, "%g %g\n", Eph[i]/1.6E-12, Fph[i]);
	}
	fclose(photons);

	printLog("read electrons distribution\n");
	double** Fe = new double*[Ndist];
	double** dFe = new double*[Ndist];
	double** Ee = new double*[Ndist];
	int powerlawStart[Ndist] = {140,140,140,140,140,140,140,140,140,140};
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
			double Emax = 4.0E3 * massElectron * speed_of_light2;
			double factor = pow(Emax / Ee[j][0], 1.0 / (Np - 1));
			for (int i = 1; i < Np; ++i) {

				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				//todo massRelationSqrt?
				double gamma = u*realMassRelationSqrt/massRelationSqrt + 1;
				//Ee[j][i] = gamma*massElectron*speed_of_light2;
				Ee[j][i] = Ee[j][i-1]*factor;
				Fe[j][i] = Fe[j][i] / (massElectron*speed_of_light2);
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);

				double minGamma = 1.0;
				double power = 3.5;

				//if(gamma >= minGamma){
				if(Ee[j][i] >= 1.0* massElectron * speed_of_light2){
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
				//gamma = u + 1;
				//if( u < 3000){
				Ee[j][i] = gamma*massElectron*speed_of_light2;
				//maxEnergy = Ee[i];
				Fe[j][i] = Fe[j][i] / (massElectron*speed_of_light2);
				if(i > 137){
					Fe[j][i] = Fe[j][137]*pow(Ee[j][i]/Ee[j][137], -3.5);
				}
				dFe[j][i] = (Fe[j][i] / (4*pi)) * (Ee[j][i] - Ee[j][i - 1]);
			}

		} else if (input == LONG_COMBINED) {
			for (int i = 1; i < Np; ++i) {
				fscanf(inputPe, "%lf", &u);
				fscanf(inputFe, "%lf", &Fe[j][i]);

				//todo massRelationSqrt?
				double gamma = u * realMassRelationSqrt / massRelationSqrt + 1;
				//gamma = u + 1;
				//if( u < 3000){
				Ee[j][i] = gamma * massElectron * speed_of_light2;
				//maxEnergy = Ee[i];
				Fe[j][i] = Fe[j][i] / (massElectron * speed_of_light2);
				if (i > 137) {
					Ee[j][i] = Ee[j][i - 1] * 1.2;
					Fe[j][i] = Fe[j][137] * pow(Ee[j][i] / Ee[j][137], -3.5);
				 }
				dFe[j][i] = (Fe[j][i] / (4 * pi)) * (Ee[j][i] - Ee[j][i - 1]);
			}

		}


		fclose(inputPe);
		fclose(inputFe);
	}

	for(int i = 0; i < Ndist; ++i){
		double norm = 0;
		for(int j = 1; j < Np; ++j){
			norm += Fe[i][j]*(Ee[i][j] - Ee[i][j-1]);
		}
		for(int j = 0; j < Np; ++j){
			Fe[i][j] = Fe[i][j]/norm;
			dFe[i][j] = dFe[i][j]/norm;
		}
	}

	printLog("initialize fields\n");
	const int Nrho = 1;
	double fraction = 0.5;
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
		volume[i] = new double*[NthetaSpace];
		B3d[i] = new double*[NthetaSpace];
		Bx3d[i] = new double*[NthetaSpace];
		By3d[i] = new double*[NthetaSpace];
		Bz3d[i] = new double*[NthetaSpace];
		concentrations3d[i] = new double*[NthetaSpace];
		sintheta3d[i] = new double*[NthetaSpace];
		thetaIndex3d[i] = new int*[NthetaSpace];
		for(int j = 0; j < NthetaSpace; ++j){
			volume[i][j] = new double[Nphi];
			B3d[i][j] = new double[Nphi];
			Bx3d[i][j] = new double[Nphi];
			By3d[i][j] = new double[Nphi];
			Bz3d[i][j] = new double[Nphi];
			concentrations3d[i][j] = new double[Nphi];
			sintheta3d[i][j] = new double[Nphi];
			thetaIndex3d[i][j] = new int[Nphi];
			for(int k = 0; k < Nphi; ++k){
				volume[i][j][k] = ((pow(r+dr/2,3) - pow(r-dr/2,3))/3.0)*dphi*dcosThetaSpace;
				B3d[i][j][k] = 1.0;
				Bx3d[i][j][k] = 1.0;
				By3d[i][j][k] = 0.0;
				Bz3d[i][j][k] = 0.0;
				if(geometry == SPHERICAL) {
					//concentrations3d[i][j][k] = 1.0*sqr(rmax/r);
					concentrations3d[i][j][k] = 1.0;
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
		initializeParker(Nrho, NthetaSpace, Nphi, Bx3d, By3d, Bz3d, rho, cosThetaSpace, sinPhiValue, cosPhiValue);
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
			for(int j = 0; j < NthetaSpace; ++j){
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
							for(int j = 0; j < NthetaSpace; ++j){
								for(int k = 0; k < Nphi; ++k){
									double z = rho[i]*cosThetaSpace[j];
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
	evaluateOrientationParameters3d(Nrho, NthetaSpace, Nphi, B3d, sintheta3d, thetaIndex3d, Bx3d, By3d, Bz3d, Ndist, rho, cosThetaSpace, sinPhiValue, cosPhiValue);

	int iangle = 7;

	double* E = new double[Nnu];
	double* I = new double[Nnu];
	double Emin = 0.0001*kBoltzman*Tphotons1;
	double Emax = 1000*(2*Ee[iangle][Np-1] + Eph[Np-1]);
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

	omp_lock_t write_lock;
	omp_init_lock(&write_lock);

	#pragma omp parallel for shared(write_lock, logFile, re2, E, I, concentrations3d, thetaIndex3d, volume, rho, Ee, Fe, Eph, Fph, cosThetaInitial, cosThetaLeftInitial, cosThetaFinal, cosThetaLeftFinal, sinPhiValue, cosPhiValue) private(iangle)
	for(int i = 0; i < Nnu; ++i){
		omp_set_lock(&write_lock);
		logFile = fopen("log.dat","a");
		printf("i nu = %d\n", i);
		fprintf(logFile, "i nu = %d\n", i);
		fclose(logFile);
		omp_unset_lock(&write_lock);
		double photonFinalEnergy = E[i];

		////Uvarov  Bykov Chevalier Elixon Uvarov 1999
		if(solver == Solver::UVAROV){
			for(int ir = 0; ir < Nrho; ++ir){
				for(int itheta = 0; itheta < NthetaSpace; ++itheta){
					for(int iphi = 0; iphi < Nphi; ++iphi){
						iangle = thetaIndex3d[ir][itheta][iphi];
						for(int k = 0; k < Np; ++k){
							double electronInitialEnergy = Ee[iangle][k];
							double delectronEnergy;
							double electronInitialGamma = electronInitialEnergy/(massElectron*speed_of_light2);
							double electronInitialBeta = sqrt(1.0 - 1.0/(electronInitialGamma*electronInitialGamma));
							if(k == 0){
								delectronEnergy = Ee[iangle][1] - Ee[iangle][0];
							} else {
								delectronEnergy = Ee[iangle][k] - Ee[iangle][k-1];
							}
							for(int l = 0; l < Nphot; ++l){
								double photonInitialEnergy = Eph[l];
								double dphotonInitialEnergy;
								if(l == 0){
									dphotonInitialEnergy = Eph[1] - Eph[0];
								} else {
									dphotonInitialEnergy = Eph[l] - Eph[l-1];
								}

								if(electronInitialGamma > photonFinalEnergy/(massElectron*speed_of_light2)){
									double G = 4*electronInitialGamma*photonInitialEnergy/(massElectron*speed_of_light2);
									double q = (photonFinalEnergy/(massElectron*speed_of_light2))/((electronInitialGamma - photonFinalEnergy/(massElectron*speed_of_light2))*G);
									if( q <= 1.0){
										double sigma = 2*pi*re2*(2*q*log(q) + 1 + q - 2*q*q + 0.5*q*q*(1-q)*G*G/(1 + q*G))/(electronInitialGamma*electronInitialGamma*photonInitialEnergy/(massElectron*speed_of_light2));
										//divide by energy to get number of photons
										I[i] += electronConcentration*concentrations3d[ir][itheta][iphi]*volume[ir][itheta][iphi]*sigma*Fph[l]*Fe[iangle][k]*speed_of_light*electronInitialBeta*delectronEnergy*dphotonInitialEnergy/(massElectron*speed_of_light2);
										//I[i] += electronConcentration*concentrations3d[ir][itheta][iphi]*volume[ir][itheta][iphi]*sigma*Fph[l]*Fe[iangle][k]*speed_of_light*electronInitialBeta*delectronEnergy*dphotonInitialEnergy;
										if(I[i] < 0){
											printf("I[i] < 0\n");
											printLog("I[i] < 0\n");
											exit(0);
										}

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

		} else if(solver == Solver::DUBUS){
			////Dubus 2008
		
			for(int j = 0; j < NthetaFinal; ++j){
				//integral by dmu1
				//integration by phi changes to 2 pi?
				double photonFinalCosTheta = cosThetaLeftFinal[j];
				int ir = 0;
				int itheta = 0;
				int iphi = 0;
				for(int ir = 0; ir < Nrho; ++ir){
					//integral by volume dr
					for(int itheta = 0; itheta < NthetaSpace; ++itheta){
						//integral by volume dtheta
						for(int iphi = 0; iphi < Nphi; ++iphi){
							//intehral by volume dphi
							iangle = thetaIndex3d[ir][itheta][iphi];
							for(int k = 0; k < Np; ++k){
								//integral by dE electron
								//integral by dE photon is eaten by delta function
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

								double electronDist = Fe[iangle][k];
								double photonFinalEnergyPrimed;
								double photonFinalCosThetaPrimed;
								LorentzTransformationPhotonZ(electronInitialGamma, photonFinalEnergy, photonFinalCosTheta, photonFinalEnergyPrimed, photonFinalCosThetaPrimed);
								double photonFinalSinThetaPrimed = sqrt(1.0 - photonFinalCosThetaPrimed*photonFinalCosThetaPrimed);
								for(int l = 0; l < NthetaInitial; ++l){
									//integral by dmu0
									/*logFile = fopen("log.dat","a");
									printf("l photon initial theta = %d\n", l);
									fprintf(logFile, "l photon initial theta = %d\n", l);
									fclose(logFile);*/
									double photonInitialCosTheta = -cosThetaLeftInitial[l];
									double photonInitialCosThetaPrimed = (photonInitialCosTheta - electronInitialBeta)/(1.0 - electronInitialBeta*photonInitialCosTheta);
									double photonInitialSinThetaPrimed = sqrt(1.0 - photonInitialCosThetaPrimed*photonInitialCosThetaPrimed);

									double coef = electronConcentration*concentrations3d[ir][itheta][iphi]*volume[ir][itheta][iphi]*(re2*speed_of_light*sqr(1.0 - electronInitialBeta*photonInitialCosTheta)/(2*sqr(electronInitialGamma)*(1.0 - electronInitialBeta*photonFinalCosTheta)));
									double derivativeEE = electronInitialGamma-electronInitialGamma*electronInitialBeta*photonInitialCosTheta;
									//printf("%g\n", derivativeEE);
									//derivativeE = 1.0;
									for(int m = 0; m < Nphi; ++m){
										//integral by dphi0
										double photonInitialPhi = phi[m];
										double cosXiPrimed = photonInitialCosThetaPrimed*photonFinalCosThetaPrimed + photonInitialSinThetaPrimed*photonFinalSinThetaPrimed*cosPhiValue[m];
										double derivativeXi = -sqr(photonFinalEnergyPrimed) / (sqr(1 - photonFinalEnergyPrimed * (1 - cosXiPrimed) / (massElectron * speed_of_light2)) * massElectron * speed_of_light2);
										double derivativeXiMu0 = photonFinalCosThetaPrimed - cosPhiValue[m] * photonFinalSinThetaPrimed * photonInitialCosThetaPrimed / (photonInitialSinThetaPrimed);
										if (photonInitialSinThetaPrimed <= 0) {
											//because derivativeMu0E contains sin^2
											derivativeXiMu0 = 0;
										}
										double derivativeXiMu1 = photonInitialCosThetaPrimed - cosPhiValue[m] * photonInitialSinThetaPrimed * photonFinalCosThetaPrimed / (photonFinalSinThetaPrimed);
										if (photonFinalSinThetaPrimed <= 0) {
											//because derivativeMu1E contains sin^2
											derivativeXiMu1 = 0;
										}
										double derivativeMu0E = ((sqr(photonInitialCosTheta) - 1) / sqr(electronInitialBeta * photonInitialCosTheta - 1)) / (sqr(electronInitialGamma) * sqrt(electronInitialEnergy - massElectron * speed_of_light2));
										if (electronInitialGamma <= 1.0) {
											//i don't know why
											derivativeMu0E = 0;
										}
										double derivativeMu1E = ((sqr(photonFinalCosTheta) - 1) / sqr(electronInitialBeta * photonFinalCosTheta - 1)) / (sqr(electronInitialGamma) * sqrt(electronInitialEnergy - massElectron * speed_of_light2));
										if (electronInitialGamma <= 1.0) {
											//i don't know why
											derivativeMu1E = 0;
										}
										double photonInitialEnergyPrimed = photonFinalEnergyPrimed/(1.0 - (photonFinalEnergyPrimed/(massElectron*speed_of_light2))*(1.0 - cosXiPrimed));

										double photonInitialEnergy = electronInitialGamma*photonInitialEnergyPrimed + electronInitialBeta*electronInitialGamma*photonInitialEnergyPrimed*photonInitialCosThetaPrimed;
										//double derivativeE = (1.0 + (photonInitialEnergyPrimed / (massElectron * speed_of_light2)) * (1.0 - cosXiPrimed)) / sqr(1.0 - (photonInitialEnergyPrimed / (massElectron * speed_of_light2)) * (1.0 - cosXiPrimed));
										double derivativeE = 1.0;

										//derivative??
										I[i] += coef*
											(1.0 + cosXiPrimed*cosXiPrimed + sqr(photonFinalEnergyPrimed/(massElectron*speed_of_light2))*sqr(1.0 - cosXiPrimed)/(1.0 - (photonFinalEnergyPrimed/(massElectron*speed_of_light2))*(1.0 - cosXiPrimed)))*
											//2*pi*dcosThetaFinal[j]*dphi*dcosThetaInitial*delectronEnergy*Fe[iangle][k]*evaluatePhotonDistribution(photonInitialEnergy, Nphot, Eph, Fph)/(derivativeE*derivativeEE - derivativeXi*(derivativeXiMu0*derivativeMu0E + derivativeXiMu1 * derivativeMu1E));
											2*pi*dcosThetaFinal[j]*dphi*dcosThetaInitial[l] * delectronEnergy * Fe[iangle][k] * evaluatePhotonDistribution(photonInitialEnergy, photonInitialCosTheta, photonInitialPhi, Nphot, Eph, Fph);
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
		else if (solver == Solver::THOMPSON) {
			if (input != Input::POWERLAW) {
				printf("Thompson solver works only with powerlaw distribution\n");
				exit(0);
			}
			//double photonConcentration = 3.536638846E7; //for CMB
			double photonConcentration = 400; //for CMB
			double sigmat = 8 * pi * re2 / 3.0;
			double index = 3.5;
			double E0 = 1.0 * massElectron * speed_of_light2;
			double K = (index - 1) * pow(E0, index - 1);
			for (int ir = 0; ir < Nrho; ++ir) {
				for (int itheta = 0; itheta < NthetaSpace; ++itheta) {
					for (int iphi = 0; iphi < Nphi; ++iphi) {
						I[i] += photonConcentration*electronConcentration * concentrations3d[ir][itheta][iphi] * volume[ir][itheta][iphi] * 0.5 * sigmat * pow(massElectron * speed_of_light2, 1.0 - index) * pow(photonFinalEnergy, -(index + 1.0) / 2.0) * pow((4.0/3.0)*2.7012 * kBoltzman * 2.7, (index - 1.0) / 2.0) *K*speed_of_light / (4 * pi);
					}
				}
			}
		}
		else if (solver == Solver::ANISOTROPIC_KN) {
			double mec2 = massElectron * speed_of_light2;
			double totalVolume = 4 * pi * rmax * rmax * rmax / 3;
			double photonFinalCosTheta = 1.0;
			double photonFinalPhi = 0;
			for (int k = 0; k < Np; ++k) {
				double electronInitialEnergy = Ee[iangle][k];
				double delectronEnergy;
				double electronInitialGamma = electronInitialEnergy / (massElectron * speed_of_light2);
				double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));
				if (k == 0) {
					delectronEnergy = Ee[iangle][1] - Ee[iangle][0];
				}
				else {
					delectronEnergy = Ee[iangle][k] - Ee[iangle][k - 1];
				}
				double electronDist = Fe[iangle][k];
				
				for (int imue = 0; imue < NthetaElectron; ++imue) {
					double mu_e = cosThetaElectron[imue];
					for (int iphie = 0; iphie < Nphi; ++iphie) {
						double phi_e = 2 * pi * (iphie + 0.5) / Nphi;
						double dphi_e = 1.0 / Nphi;
						//rotation
						double photonFinalCosThetaRotated;
						double photonFinalPhiRotated;
						rotationSphericalCoordinates(thets_e, phi_e, photonFinalTheta, photonFinalPhi, photonFinalThetaRotated, photonFinalPhiRotated);
						//double photonFinalEnergyPrimed;
						//double photonFinalCosThetaPrimed;
						//LorentzTransformationPhotonZ(electronInitialBeta, photonFinalEnergy, photonFinalCosThetaRotated, photonFinalEnergyPrimed, photonFinalCosThetaPrimed);
						double photonFinalEnergyPrimed = electronInitialGamma * (1 - photonFinalCosThetaRotated * electronInitialBeta) * photonFinalEnergy;
						double photonFinalCosThetaPrimed = (photonFinalCosThetaRotated - electronInitialBeta) / (1.0 - electronInitialBeta * photonFinalCosThetaRotated);
						double photonFinalSinThetaPrimed = sqrt(1.0 - photonFinalCosThetaPrimed * photonFinalCosThetaPrimed);
						for (int imuph = 0; imuph < NthetaInitial; ++imuph) {
							double mu_ph = -cosThetaInitial[imuph];
							for (int iphiph = 0; iphiph < Nphi; ++iphiph) {
								double phi_ph = 2 * pi * (iphiph + 0.5) / Nphi;
								double dphi_ph = 1.0 / Nphi;
								double photonInitialCosThetaPrimed = mu_ph;
								double photonInitialPhiRotated = phi_ph;
								//inverseRotationSphericalCoordinates(mu_e, phi_e, mu_ph, phi_ph, photonInitialCosThetaRotated, photonInitialPhiRotated);
								//double photonInitialCosThetaPrimed = (photonInitialCosThetaRotated - electronInitialBeta) / (1.0 - electronInitialBeta * photonInitialCosThetaRotated);
								
								double photonInitialSinThetaPrimed = sqrt(1.0 - photonInitialCosThetaPrimed * photonInitialCosThetaPrimed);
								
								double cosXiPrimed = photonInitialCosThetaPrimed * photonFinalCosThetaPrimed + photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);

								double photonInitialEnergyPrimed = photonFinalEnergyPrimed / (1.0 - (photonFinalEnergyPrimed / (mec2)) * (1.0 - cosXiPrimed));
								double photonInitialEnergy = electronInitialGamma * photonInitialEnergyPrimed + electronInitialBeta * electronInitialGamma * photonInitialEnergyPrimed * photonInitialCosThetaPrimed;
								double photonInitialCosThetaRotated = (photonInitialCosThetaPrimed + electronInitialBeta) / (1 + photonInitialCosThetaPrimed * electronInitialBeta);

								double photonInitialCosTheta;
								double photonInitialPhi;
								inverseRotationSphericalCoordinates(theta_e, phi_e, photonInitialThetaRotated, photonInitialPhiRotated, photonInitialTheta, photonInitialPhi);

								double dI = electronConcentration * totalVolume * 0.5 * re2 * speed_of_light * Fe[iangle][k] *
									(sqr(1 - photonInitialCosThetaRotated * electronInitialBeta) / (1.0 - photonFinalCosThetaRotated * electronInitialBeta)) *
									(1 + cosXiPrimed * cosXiPrimed + sqr(photonFinalEnergyPrimed / mec2) * sqr(1 - cosXiPrimed) / (1 - (photonFinalEnergyPrimed / mec2) * (1 - cosXiPrimed))) *
									evaluatePhotonDistribution(photonInitialEnergy, photonInitialCosTheta, photonInitialPhi, Nphot, Eph, Fph) *
									dphi_e * dphi_ph * dcosThetaElectron[imue] * dcosThetaInitial[imuph] * delectronEnergy;

								if (dI < 0) {
									printf("dI[i] <  0\n");
									printLog("dI[i] < 0\n");
									exit(0);
								}

								I[i] += dI;
								if (I[i] != I[i]) {
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

	omp_destroy_lock(&write_lock);

	printLog("output\n");

	FILE* output_ev_EFE = fopen("outputE.dat","w");
	FILE* output_GHz_Jansky = fopen("outputNu.dat","w");
	for(int i = 0; i < Nnu; ++i){
		double nu = E[i]/hplank;
		fprintf(output_ev_EFE, "%g %g\n", E[i]/(1.6E-12), E[i]*E[i]*I[i]/sqr(distance));
		fprintf(output_GHz_Jansky, "%g %g\n", nu/1E9, 1E26*hplank*E[i] * I[i] / sqr(distance));
	}
	fclose(output_ev_EFE);
	fclose(output_GHz_Jansky);

	printLog("evaluate total luminosity\n");
	double totalLuminosity = 0;
	double totalFlux = 0;
	double minEev = 0.3 * 1000;
	double maxEev = 10 * 1000;
	for (int i = 0; i < Nnu; ++i) {
		double Eev = E[i] / (1.6E-12);
		if ((Eev > minEev) && (Eev < maxEev)) {
			totalLuminosity += 4 * pi * E[i] * I[i]*(E[i] - E[i-1]);
			totalFlux += E[i] * I[i] * (E[i] - E[i - 1])/sqr(distance);
		}
	}
	printf("total luminosity = %g erg/s \n", totalLuminosity);
	printf("total flux = %g erg/s cm^2 \n", totalFlux);
	printf("99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39\n");
	//CSS161010
	//99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39
	//130 days F = 1.94+-0.74 10^-15 L = 5.0+-2.5 10^39
	//291 days F << 1.31 10^-15 L << 3.4 10^39
	//H density 4.7*10^20 cm^-2 what does it means?
	double theyRatio = 3.4E39 / 1.33E-15;
	double myRatio = totalLuminosity / totalFlux;

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

