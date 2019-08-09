#include "stdio.h"
#include <stdlib.h>
#include <time.h>
#include "math.h"

const double kBoltzman = 1.3806488E-16;
const double speed_of_light = 2.99792458E10;
const double speed_of_light2 = speed_of_light*speed_of_light;
const double speed_of_light4 = speed_of_light2*speed_of_light2;
const double electron_charge = 4.803529695E-10;
const double pi = 4*atan2(1.0,1.0);
const double massElectronReal = 0.910938291E-27;
const double massProtonReal = 1.67262177E-24;
const double massRatio = 25;
//const double massElectron = massProtonReal/massRatio;
const double massElectron = massElectronReal;

const int Napprox  = 40;

const double McDonaldValue[Napprox] = {3.08E8, 2.11E7, 6.65E6, 4.55E5, 1.43E4, 9802, 3087, 670, 211, 107, 66.3, 33.6, 20.7, 14.1, 11.0, 10.3, 6.26, 4.2, 3.01, 2.25, 1.73, 1.37, 1.10, 0.737, 0.514, 0.368, 0.269, 0.2, 0.0994, 0.0518, 0.0278, 0.0152, 0.00846, 0.00475, 0.00154, 0.000511, 0.000172, 0.0000589, 0.0000203, 0.00000246};
const double UvarovValue[Napprox] = {
		0.0461, 0.0791, 0.0995, 0.169, 0.213, 0.358, 0.445, 0.583, 0.702, 0.772, 0.818, 0.874,0.904, 0.917, 0.918, 0.918, 0.901, 0.872, 0.832, 0.788, 
		0.742, 0.694, 0.655, 0.566, 0.486, 0.414, 0.354, 0.301, 0.200, 0.130, 0.0845, 0.0541, 0.0339, 0.0214, 0.0085, 0.0033, 0.0013, 0.00050, 0.00019, 0.0000282};
const double UvarovX[Napprox] =  {0.00001, 0.00005, 0.0001 ,0.0005, 0.001, 0.005, 0.01, 0.025, 0.050, 0.075, 0.10, 0.15,
		0.20, 0.25, 0.29, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0, 1.2, 1.4, 1.6, 
		1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0};

double power(const double& v, const double& p) {
	return exp(p * log(v));
}

double criticalNu(const double& E, const double& sinhi, const double& H) {
	return 3*electron_charge*H*sinhi*E*E/(4*pi*massElectron*massElectron*massElectron*speed_of_light*speed_of_light4);
}

double evaluateMcDonaldIntegral(const double& nu) {
	int curIndex = 0;
	if(nu < UvarovX[0]) {
		return 0;
	}
	if(nu > UvarovX[Napprox - 1]) {
		return 0;
	}
	for(int i = 1; i < Napprox; ++i) {
		if(nu < UvarovX[i]) {
			curIndex = i;
			break;
		}
	}

	//double result = (UvarovValue[curIndex]*(nu - UvarovX[curIndex - 1]) + UvarovValue[curIndex - 1]*(UvarovX[curIndex] - nu))/(UvarovX[curIndex] - UvarovX[curIndex - 1]);
	double result = UvarovValue[curIndex-1]*exp(log(UvarovValue[curIndex]/UvarovValue[curIndex-1])*((nu - UvarovX[curIndex-1])/(UvarovX[curIndex] - UvarovX[curIndex - 1])));
	if(result < 0) {
		printf("result < 0\n");
	}
	return result;
}

double evaluateMcDonaldFunction(const double& nu) {
	int curIndex = 0;
	if(nu < UvarovX[0]) {
		return 0;
	}
	if(nu > UvarovX[Napprox - 1]) {
		return 0;
	}
	for(int i = 1; i < Napprox; ++i) {
		if(nu < UvarovX[i]) {
			curIndex = i;
			break;
		}
	}

	//double result = (UvarovValue[curIndex]*(nu - UvarovX[curIndex - 1]) + UvarovValue[curIndex - 1]*(UvarovX[curIndex] - nu))/(UvarovX[curIndex] - UvarovX[curIndex - 1]);
	double result = McDonaldValue[curIndex-1]*exp(log(McDonaldValue[curIndex]/McDonaldValue[curIndex-1])*((nu - UvarovX[curIndex-1])/(UvarovX[curIndex] - UvarovX[curIndex - 1])));
	if(result < 0) {
		printf("result < 0\n");
	}
	return result;
}

int main(int argc, char** argv){

	int Np = 200;
	FILE* inputPe = fopen("../../tristan-mp-pitp/Pe.dat","r");
	FILE* inputFe = fopen("../../tristan-mp-pitp/Fe.dat","r");

	double* Pe = new double[Np];
	double* Fe = new double[Np];
	double* Ee = new double[Np];

	int Nnu = 100;
	double* nu = new double[Nnu];
	double* Inu = new double[Nnu];
	double* Anu = new double[Nnu];

	double sigma = 0.36;
	double gamma0 = 1.5;
	double concentration = 1;
	double B0 = sqrt(sigma*4*pi*gamma0*concentration*massProtonReal*speed_of_light2);
	double size = 1E16;

	double Bmean = 4*B0;

	Bmean = 0.1;

	concentration = 1;


	for(int i = 0; i < Np; ++i) {
		fscanf(inputPe, "%lf", &Pe[i]);
		Pe[i] = Pe[i]*massElectron*speed_of_light;
		Ee[i] = sqrt(Pe[i]*Pe[i]*speed_of_light2 + massElectron*massElectron*speed_of_light4);
		fscanf(inputFe,"%lf", &Fe[i]);
		Fe[i] = Fe[i]*Ee[i]/(Pe[i]*Pe[i]*Pe[i]*speed_of_light2);
	}

	/*Fe[0] = 1.0;
	for(int i = 1; i < Np; ++i) {
		Fe[i] = Fe[0]*power(Ee[0]/Ee[i],2);
	}*/

	fclose(inputPe);
	fclose(inputFe);

	double norm = 0;
	for(int i = 1; i < Np; ++i) {
		norm = norm + Fe[i]*(Ee[i] - Ee[i-1]);
	}
	for(int i = 0; i < Np; ++i) {
		Fe[i] = Fe[i]/norm;
	}

	double minEnergy = 10*massElectron*speed_of_light2;
	//double minEnergy = Ee[0];
	double maxEnergy = Ee[Np-1];

	double meanE = 200*massElectron*speed_of_light2;

	int startElectronIndex = 0;
	for(int i = 0; i < Np; ++i) {
		if(Ee[i] > minEnergy) {
			startElectronIndex = i;
			break;
		}
	}

	double hi = pi/4;
	double sinhi = sin(hi);

	double minNu = 0.1*criticalNu(minEnergy, sinhi, Bmean);
	double maxNu = 10*criticalNu(maxEnergy, sinhi, Bmean);
	double meanNu = criticalNu(meanE, sinhi, Bmean);

	double temp = maxNu/minNu;

	double factor = power(maxNu/minNu, 1.0/(Nnu-1));

	nu[0] = minNu;
	Inu[0] = 0;
	Anu[0] = 0;
	for(int i = 1; i < Nnu; ++i) {
		nu[i] = nu[i-1]*factor;
		Inu[i] = 0;
		Anu[i] = 0;
	}

	double coef = concentration*sqrt(3.0)*electron_charge*electron_charge*electron_charge/(massElectron*speed_of_light2);
	double coefAbsorb = concentration*16*pi*pi*electron_charge/(3*sqrt(3.0)*Bmean*sinhi);


	for(int i = 0; i < Nnu; ++i) {
		printf("i = %d\n", i);
		for(int j = startElectronIndex; j < Np; ++j) {
			double nuc = criticalNu(Ee[j], sinhi, Bmean);
			double gamma = Ee[j]/(massElectron*speed_of_light2);
			double x = nu[i]/nuc;
			Inu[i] = Inu[i] + coef*Fe[j]*(Ee[j] - Ee[j-1])*Bmean*sinhi*evaluateMcDonaldIntegral(nu[i]/nuc);
			Anu[i] = Anu[i] + coefAbsorb*Fe[j]*(Ee[j] - Ee[j-1])*evaluateMcDonaldFunction(nu[i]/nuc)/(gamma*gamma*gamma*gamma*gamma);
		}
	}

	int absorbtionIndex = 0;
	for(int i = 0; i < Nnu; ++i) {
		if(Anu[i]*size < 1) {
			absorbtionIndex = i;
			break;
		}
	}

	/*for(int i = 0; i < absorbtionIndex; ++i) {
		Inu[i] = Inu[absorbtionIndex]*power(nu[i]/nu[absorbtionIndex], 5.0/2.0);
	}*/

	for(int i = 0; i < Nnu; ++i){
		Inu[i] = Inu[i]*(1 - exp(-Anu[i]*size))/(Anu[i]*size);
	}

	double totalFlux = 0;
	for(int i = 1; i < Nnu; ++i) {
		totalFlux += Inu[i]*(nu[i] - nu[i-1])*size*size*size;
	}


	double nuc = criticalNu(minEnergy, sinhi, Bmean);
	FILE* output = fopen("../../tristan-mp-pitp/radiation.dat", "w");
	for(int i = 0; i < Nnu; ++i) {
		fprintf(output, "%g %g %g\n", nu[i], Inu[i], Anu[i]*size);
	}
	
	fclose(output);

	delete[] Pe;
	delete[] Ee;
	delete[] Fe;

	delete[] nu;
	delete[] Inu;
	delete[] Anu;
}