#include "stdio.h"
#include <stdlib.h>
#include <time.h>
#include "math.h"
#include <string>

const double kBoltzman = 1.3806488E-16;
const double speed_of_light = 2.99792458E10;
const double speed_of_light2 = speed_of_light * speed_of_light;
const double speed_of_light4 = speed_of_light2 * speed_of_light2;
const double electron_charge = 4.803529695E-10;
const double pi = 4 * atan2(1.0, 1.0);
const double massElectronReal = 0.910938291E-27;
const double massProtonReal = 1.67262177E-24;
const double massRatio = 25;
//const double massElectron = massProtonReal/massRatio;
const double massElectron = massElectronReal;

const int distributionsN = 10;
const int Nx = 200;
const int Np = 200;

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

double power(const double& v, const double& p) {
	return exp(p * log(v));
}

double findDistribution(double* distribution, double* momentum, const double& p){
	if(p < momentum[0]){
		return 0;
	}
	if(p > momentum[Np-1]){
		return 0;
	}

	for(int i = 1; i < Np; ++i) {
		if(p < momentum[i]) {
			//return distribution[i - 1] * exp(log(distribution[i] / distribution[i - 1]) * ((p - momentum[i-1]) / (momentum[i] - momentum[i - 1])));
			return (distribution[i - 1] *(momentum[i] - p) + distribution[i]*(p - momentum[i-1]))/ (momentum[i] - momentum[i - 1]);
		}
	}
	return 0;
}

int main(){
	double angles[distributionsN+1];
	double dtheta = 90.0/distributionsN;
	for(int i = 0; i <= distributionsN; ++i){
		angles[i] = i*90.0/distributionsN;
	}
	double weights[distributionsN];
	for(int i = 0; i < distributionsN; ++i){
		weights[i] = 0;
	}

	double* theta = new double[Nx];
	FILE* thetaFile = fopen("./theta.dat","r");
	for(int i = 0; i < Nx; ++i){
		fscanf(thetaFile, "%lf", &theta[i]);
	}
	fclose(thetaFile);

	for(int i = 0; i < Nx; ++i){
	if(theta[i] > 90){
		theta[i] = 180 - theta[i];
	}
	int index = floor(theta[i]/dtheta);
		weights[index] += 1.0/Nx;
	}
	delete[] theta;

	double** distributions = new double*[distributionsN];
	double** momentums = new double*[distributionsN];
	for(int i = 0; i < distributionsN; ++i){
		distributions[i] = new double[Np];
		momentums[i] = new double[Np];
		std::string fileNumber = convertIntToString(i);
		std::string fFileName = "./Fe";
		std::string pFileName = "./Pe";
		FILE* distributionFile = fopen((fFileName + fileNumber + ".dat").c_str(),"r");
		FILE* momentumFile = fopen((pFileName + fileNumber + ".dat").c_str(),"r");

		for(int j = 0; j < Np; ++j){
			fscanf(momentumFile, "%lf", &momentums[i][j]);
			fscanf(distributionFile, "%lf", &distributions[i][j]);
		}		

		fclose(distributionFile);
		fclose(momentumFile);
	}

	double* distribution = new double[Np];
	double* momentum = new double[Np];

	for(int i = 0; i < Np; ++i){
		distribution[i] = 0;
		momentum[i] = 0;
	}

	//todo reset momentum
	double minMomentum = momentums[0][0];
	double maxMomentum = momentums[0][Np-1];

	for(int i = 0; i < distributionsN; ++i){
		if(momentums[i][0] < minMomentum){
			minMomentum = momentums[i][0];
		}
		if(momentums[i][Np-1] > maxMomentum){
			maxMomentum = momentums[i][Np-1];
		}
	}

	momentum[0] = 0.9*minMomentum;
	momentum[Np-1] = 1.1*maxMomentum;

	double pfactor = power(momentum[Np-1]/momentum[0], 1.0/(Np - 1));
	for(int i = 1; i < Np; ++i){
		momentum[i] = momentum[i - 1]*pfactor;
	}
	for(int j = 0; j < Np; ++j){
		for(int i = 0; i < distributionsN; ++i){
			distribution[j] += findDistribution(distributions[i], momentums[i], momentum[j])*weights[i];
		}
	}

	/*for(int i = 0; i < distributionsN; ++i){
	for(int j = 0 ; j < Np; ++j){
	distribution[j] += distributions[i][j]*weights[i];
	}
	}*/

	for(int i = 0; i < distributionsN; ++i){
		delete[] momentums[i];
		delete[] distributions[i];
	}
	delete[] momentums;
	delete[] distributions;

	FILE* outFile =fopen("./average_distribution.dat","w");
	for(int i = 0; i < Np; ++i){
		fprintf(outFile, "%g %g\n", momentum[i], distribution[i]);
	}
	fclose(outFile);
	delete[] momentum;
	delete[] distribution;

	return 0;
}