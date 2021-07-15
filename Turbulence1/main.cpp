// Turbulence1.cpp : Defines the entry point for the console application.
//
#include <cstdlib>
#include <math.h>
#include "stdio.h"

const double pi = 4*atan2(1.0,1.0);

double uniformDistribution() {
	return (rand() % 2048 + 0.5) / 2048;
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


int main(){
	const int Nd = 10;
	const int Nx = 1000;
	const int Nk = 10;
	double x[Nx];
	double y[Nx];
	double z[Nx];
	double Bx[Nx];
	double By[Nx];
	double Bz[Nx];
	double weights[Nd];

	double sintheta = 1.0;
	double turbulenceFraction = 0.9;
	double B0 = 1.0*sqrt(1.0 - turbulenceFraction);

	for( int i = 0; i < Nx; ++i){
		Bx[i] = B0*sqrt(1.0 - sintheta*sintheta);
		By[i] = -B0*sintheta;
		Bz[i] = 0;
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

	FILE* outFile = fopen("turbulentB.dat", "w");

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
		double delta = 90 - acos(Bz[i]/B)*180/pi;
		double lambda = atan2(By[i], Bx[i])*180/pi;
		if(lambda < 0){
			lambda = lambda + 360;
		}
		fprintf(outFile, "%g %g %g %g %g %g %g\n", B, Bx[i], By[i], Bz[i], theta1*180/pi, lambda, delta);
	}

	fclose(outFile);

	return 0;
}

