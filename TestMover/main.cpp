#include "stdio.h"
#include <stdlib.h>
#include <time.h>
#include "math.h"

const double c = 2.99792458E10;
const double m = 0.910938291E-27;
const double q = 4.803529695E-10;
const double pi = 4*atan(1.0);
const double gamma = 10;
const double Bmeansqr = 1.0;
const double turbulenceFraction = 0.9;
const int randomParameter = 1024;
double turbulenceLength;

double uniformDistribution() {
	return (rand() % randomParameter + 0.5) / randomParameter;
}

double power(const double& v, const double& p) {
	return exp(p * log(v));
}

double evaluateTurbulentB(const double& kx, const double& ky, const double& kz){
		double kw, v;

		kw = sqrt(kx*kx + ky*ky + kz*kz);

		return 1.0/sqrt(power(kw, 8.0/3.0));

		//return 1.0/sqrt(power(kw, 11.0/3.0));

}

void getBfield(const double& x, const double& y, const double& z, double& Bx, double& By, double& Bz, double& B0x, double& B0y, double& B0z){
	int randomseed = 239;
	srand(randomseed);
	int Nx = 4;
	int Ny = 4;
	double dk = 2*pi/turbulenceLength;
	Bx = B0x;
	By = B0y;
	Bz = B0z;
	double B0sqr = B0x*B0x + B0y*B0y + B0z*B0z;
	double Bsqr = 0;
	double turbulentFieldCorrection = 1.0;
	for(int i = 0; i <= Nx; ++i){
		for(int j = 0; j <= Nx; ++j){
			if((i + j) > 0){
				double kx = i*dk;
				double ky = j*dk;
				double kz = 0;
				double B = evaluateTurbulentB(kx, ky, kz);
				Bsqr = Bsqr + B*B;
			}
		}
	}
	if(Bsqr > 0){
		turbulentFieldCorrection = sqrt(turbulenceFraction*B0sqr/(Bsqr*(1.0 - turbulenceFraction)));
	}
	for(int i = 0; i <= Nx; ++i){
		for(int j = 0; j <= Nx; ++j){
			if((i + j) > 0){
				double kx = i*dk;
				double ky = j*dk;
				double kz = 0;
				double B = evaluateTurbulentB(kx, ky, kz)*turbulentFieldCorrection;
				
				double phase1 = 2*pi*uniformDistribution();
				double phase2 = 2*pi*uniformDistribution();

				double kw = sqrt(kx*kx + ky*ky + kz*kz);
				double kyz = sqrt(ky*ky + kz*kz);
				double cosTheta = kx/kw;
				double sinTheta = kyz/kw;
				double cosPhi = 1.0;
				double sinPhi = 0.0;
				if(ky + kz > 0) {
					cosPhi = ky/kyz;
					sinPhi = kz/kyz;
				} else {
					cosPhi = 1.0;
					sinPhi = 0.0;
				}

				double kmultr = kx*x + ky*y + kz*z;
				double localB1 = B*sin(kmultr + phase1);
				double localB2 = B*sin(kmultr + phase2);

				Bx = Bx - localB1*sinTheta;
				By = By + (localB1*cosTheta*cosPhi - localB2*sinPhi);
				Bz = Bz + (localB1*cosTheta*sinPhi + localB2*cosPhi);

			}
		}
	}
}

void move(double& x, double& y, double& z, double& vx, double& vy, double& vz, const double& Bx, const double& By, const double& Bz, const double& dt){
	double	Bnorm = sqrt(Bx*Bx + By*By + Bz*Bz);
	x = x + vx*dt;
	y = y + vy*dt;
	z = z + vz*dt;
	double tempvx = vx;
	double tempvy = vy;
	double tempvz = vz;
	if (Bnorm > 0) {
		double ex = Bx/Bnorm;
		double ey = By/Bnorm;
		double ez = Bz/Bnorm;

		double dphi = -q*Bnorm*dt/(gamma*m*c);
				
		double sindphi = sin(dphi);
		double cosdphi = cos(dphi);

		double mxx = cosdphi + (1-cosdphi)*ex*ex;
		double mxy = (1-cosdphi)*ex*ey-sindphi*ez;
		double mxz = (1-cosdphi)*ex*ez+sindphi*ey;

		double myx = (1-cosdphi)*ex*ey+sindphi*ez;
		double myy = cosdphi + (1-cosdphi)*ey*ey;
		double myz = (1-cosdphi)*ey*ez-sindphi*ex;

		double mzx = (1-cosdphi)*ex*ez-sindphi*ey;
		double mzy = (1-cosdphi)*ey*ez+sindphi*ex;
		double mzz = cosdphi + (1-cosdphi)*ez*ez;

		tempvx = mxx*vx + mxy*vy + mxz*vz;
		tempvy = myx*vx + myy*vy + myz*vz;
		tempvz = mzx*vx + mzy*vy + mzz*vz;

	}
	//x = x + vx*dt/2;
	//y = y + vy*dt/2;
	//z = z + vz*dt/2;

	vx = tempvx;
	vy = tempvy;
	vz = tempvz;

	double v = sqrt(vx*vx + vy*vy + vz*vz);
	if(v > c){
		printf("v > c\n");
	}
}

int main(int argc, char** argv){
	const int Nt = 1000000;
	double dt;


	double vx, vy, vz;
	double x, y, z;
	double Bx, By, Bz;
	double B0x, B0y, B0z;

	B0x = 0;
	B0y = 0;
	B0z = Bmeansqr*sqrt(1.0 - turbulenceFraction);

	double v = c*sqrt(1.0 - 1.0/(gamma*gamma));
	//double Bmean = sqrt((B0x*B0x + B0y*B0y + B0z*B0z)/(1.0 - turbulenceFraction));
	double omega = Bmeansqr*q/(gamma*m*c);
	double rg = v/omega;
	dt = 0.1/omega;

	turbulenceLength = 10*rg;

	
	srand(time(NULL));
	double theta = pi*uniformDistribution();
	double phi = 2*pi*uniformDistribution();

	vx = v*cos(theta);
	vy = v*sin(theta)*cos(phi);
	vz = v*sin(theta)*sin(phi);
	x = 0;
	y = 0;
	z = 0;

	FILE* information = fopen("information.dat","w");
	fprintf(information, "rg = %g\n", rg);
	fprintf(information, "omega = %g\n", omega);
	fprintf(information, "dt = %g\n", dt);
	fprintf(information, "lambda/rg = %g\n", turbulenceLength/rg);
	fclose(information);

	FILE* out = fopen("trajectory.dat","w");
	int writeParameter = 1;

	for(int i = 0; i < Nt; ++i){
		if(i%1000 == 0){
			printf("interation %d\n", i);
		}
		v = c*sqrt(1.0 - 1.0/(gamma*gamma));
		if(i%writeParameter == 0){
			fprintf(out, "%g %g %g %g %g %g %g\n", x, y, z, vx, vy, vz, v);
		}
		getBfield(x, y, z, Bx, By, Bz, B0x, B0y, B0z);
		move(x, y, z, vx, vy, vz, Bx, By, Bz, dt);
	}

	fclose(out);

	FILE* field = fopen("field.dat","w");

	double dx = rg/10;
	int Nx = 500;
	int Ny =500;
	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j< Ny; ++j){
			x = i*dx;
			y = j*dx;
			z = 0;
			getBfield(x, y, z, Bx, By, Bz, B0x, B0y, B0z);
			fprintf(field, "%g %g %g\n", Bx, By, Bz);
		}
	}

	fclose(field);

	return 0;
}