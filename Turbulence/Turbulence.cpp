// Turbulence.cpp : Defines the entry point for the console application.
//
#include <time.h>

#include "constants.h"
#include "random.h"
#include "util.h"

int seed = 0;

void evaluateTurbulence(const double& x, const double& y, const double& z, const double& Lbox, const double& Cslab, const double& LmaxSlab, const double& powerSlab, const int Nslab, const double& C2d, const double& Lmax2d, const double& power2d, const int N2d, double& Bx, double& By, double& Bz){
	Bx = 0;
	By = 0;
	Bz = 0;

	srand(seed);

	//slab
	for(int i = 1; i <= Nslab; ++i){
		double phase1 = 2*pi*uniformDistribution();
		double phase2 = 2*pi*uniformDistribution();
		double k = 2*pi*i/Lbox;
		if(k > 2*pi/LmaxSlab){
			Bx = Bx + sqrt(Cslab*power(k, powerSlab))*sin(k*z + phase1);
			By = By + sqrt(Cslab*power(k, powerSlab))*sin(k*z + phase2);
		}
	}

	//2d
	for(int i = 0; i <= N2d; ++i){
		for(int j = 0; j <= N2d; ++j){
			if(i + j > 0){
				double phase1 = 2*pi*uniformDistribution();
				double phase2 = 2*pi*uniformDistribution();
				double kx = 2*pi*i/Lbox;
				double ky = 2*pi*j/Lbox;
				double k = sqrt(kx*kx + ky*ky);

				if(k > 2*pi/Lmax2d){
					double B = sqrt(C2d*power(k, power2d));

					double phi = atan2(kx, ky);
					double sinphi = sin(phi);
					double cosphi = cos(phi);
			
					Bx = Bx + B*sin(kx*x + ky*y + phase1)*sinphi;
					By = By - B*sin(kx*x + ky*y + phase1)*cosphi;
					Bz = Bz + B*sin(kx*x + ky*y + phase2);
				}
			}
		}
	}
}

int main()
{
	double B0 = 1.0;
	double fractionSlab = 0.8;
	double fraction2d = 0.2;

	double powerSlab = -5.0/3.0;
	double power2d = -5.0/3.0;

	double Cslab = 1.0;
	double C2d = 1.0;

	double Lbox = 1E7;
	double LmaxSlab = 1E7;
	double Lmax2d = 1E7;

	int Nslab = 2000;
	int N2d = 2000;

	double integralSlab = (powerSlab + 1.0)*(-power(2*pi/LmaxSlab,powerSlab + 1.0) + power(2*pi*Nslab/LmaxSlab, powerSlab + 1.0));
	double integral2d = (power2d + 1.0)*(-power(2*pi/Lmax2d,power2d + 1.0) + power(2*pi*N2d/Lmax2d, power2d + 1.0));

	Cslab = fractionSlab*B0*B0/(integralSlab*4*pi);
	C2d = fraction2d*B0*B0/(integral2d*4*pi);

	double x = 0;
	double y = Lmax2d*uniformDistribution();
	double z = LmaxSlab*uniformDistribution();
	int Nx = N2d;
	double dx = Lbox/Nx;

	double Bx = 0;
	double By = 0;
	double Bz = 0;

	srand(time(NULL));
	seed = rand();

	FILE* outFile = fopen("B.dat","w");

	for(int i = 0; i < Nx; ++i){
		printf("%d\n",i);
		x = x + dx;
		Bx = 0;
		By = 0;
		Bz = 0;
		evaluateTurbulence(x, y, z, Lbox, Cslab, LmaxSlab, powerSlab, Nslab, C2d, Lmax2d, power2d, N2d, Bx, By, Bz);
		Bz = Bz + B0;
		double theta = acos(Bx/sqrt(Bx*Bx + By*By + Bz*Bz))*180/pi;
		fprintf(outFile, "%lf %lf %lf %lf\n",Bx, By, Bz, theta);
	}
	fclose(outFile);

	/*FILE* outBx = fopen("Bx.dat","w");
	FILE* outBy = fopen("By.dat","w");
	FILE* outBz = fopen("Bz.dat","w");
	FILE* outTheta = fopen("Theta.dat","w");

	x = 0;
	y = 0;
	double dy = dx;
	for(int i = 0; i < Nx; ++i){
		printf("%d\n",i);
		x= x + dx;
		for(int j = 0; j < Nx; ++j){
			y = y + dy;
			Bx = 0;
			By = 0;
			Bz = 0;
			evaluateTurbulence(x, y, z, Lbox, Cslab, LmaxSlab, powerSlab, Nslab, C2d, Lmax2d, power2d, N2d, Bx, By, Bz);
			Bz = Bz + B0;
			double theta = acos(Bx/sqrt(Bx*Bx + By*By + Bz*Bz))*180/pi;
			if(theta > 90.0){
				theta = 180 - theta;
			}
			fprintf(outBx, "%lf ", Bx);
			fprintf(outBy, "%lf ", By);
			fprintf(outBz, "%lf ", Bz);
			fprintf(outTheta, "%lf ", theta);
		}
		fprintf(outBx, "\n");
		fprintf(outBy, "\n");
		fprintf(outBz, "\n");
		fprintf(outTheta, "\n");
	}
	fclose(outBx);
	fclose(outBy);
	fclose(outBz);
	fclose(outTheta);*/

	return 0;
}

