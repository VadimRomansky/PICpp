#include "stdio.h"
#include <stdlib.h>
#include <time.h>
#include "math.h"
#include <string>

#include "constants.h"
#include "output.h"
#include "distribution.h"

const double q = 4.803529695E-10;
const double gamma = 100;
const double Bmeansqr = 1.0;
const double turbulenceFraction = 0.9;
double turbulenceLength;
int randomSeed;

double power(const double& v, const double& p) {
	return exp(p * log(v));
}

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

double evaluateTurbulentB(const double& kx, const double& ky, const double& kz){
		double kw, v;

		kw = sqrt(kx*kx + ky*ky + kz*kz);

		return 1.0/sqrt(power(kw, 8.0/3.0));

		//return 1.0/sqrt(power(kw, 11.0/3.0));

}

void getBfield(const double& x, const double& y, const double& z, double& Bx, double& By, double& Bz, double& B0x, double& B0y, double& B0z, const int Nxmodes, const int Nymodes, double*** phases){
	//srand(randomSeed);
	double dk = 2*pi/turbulenceLength;
	Bx = B0x;
	By = B0y;
	Bz = B0z;
	double B0sqr = B0x*B0x + B0y*B0y + B0z*B0z;
	double Bsqr = 0;
	double turbulentFieldCorrection = 1.0;
	for(int i = 0; i <= Nxmodes; ++i){
		for(int j = 0; j <= Nymodes; ++j){
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
	for(int i = 0; i <= Nxmodes; ++i){
		for(int j = 0; j <= Nxmodes; ++j){
			if((i + j) > 0){
				double kx = i*dk;
				double ky = j*dk;
				double kz = 0;
				double B = evaluateTurbulentB(kx, ky, kz)*turbulentFieldCorrection;
				
				//double phase1 = 2*pi*uniformDistribution();
				//double phase2 = 2*pi*uniformDistribution();

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
				double localB1 = B*sin(kmultr + phases[i][j][0]);
				double localB2 = B*sin(kmultr + phases[i][j][1]);

				Bx = Bx - localB1*sinTheta;
				By = By + (localB1*cosTheta*cosPhi - localB2*sinPhi);
				Bz = Bz + (localB1*cosTheta*sinPhi + localB2*cosPhi);

			}
		}
	}
}


void getEandBfield(const double& x, const double& y, const double& z, double& Ex, double& Ey, double& Ez, double& Bx, double& By, double& Bz, double*** downstreamE, double*** downstreamB, double*** middleE, double*** middleB, double*** upstreamE, double*** upstreamB, int downstreamNx, int middleNx, int upstreamNx, int Ny, double dx){
	double tempy = y;
	double tempx = x;
	double Ly = Ny*dx;
	if(y > Ly){
		tempy = fmod(y, Ly);
	} else if(y < 0){
		tempy = Ly - fmod(-y,Ly);
	}

	int indexy = int(tempy/dx);

	double upstreamLx = upstreamNx*dx;
	double middleLx = middleNx*dx;
	double downstreamLx = downstreamNx*dx;

	if(x > middleLx){
		tempx = fmod(x - middleLx,upstreamLx);
		int indexx = int(tempx/dx);
		Ex = upstreamE[indexx][indexy][0];
		Ey = upstreamE[indexx][indexy][1];
		Ez = upstreamE[indexx][indexy][2];
		Bx = upstreamB[indexx][indexy][0];
		By = upstreamB[indexx][indexy][1];
		Bz = upstreamB[indexx][indexy][2];
	} else if(x >= 0){
		int indexx = int(x/dx);
		Ex = middleE[indexx][indexy][0];
		Ey = middleE[indexx][indexy][1];
		Ez = middleE[indexx][indexy][2];
		Bx = middleB[indexx][indexy][0];
		By = middleB[indexx][indexy][1];
		Bz = middleB[indexx][indexy][2];
	} else {
		tempx = downstreamLx - fmod(-x,downstreamLx);
		int indexx = int(tempx/dx);
		Ex = downstreamE[indexx][indexy][0];
		Ey = downstreamE[indexx][indexy][1];
		Ez = downstreamE[indexx][indexy][2];
		Bx = downstreamB[indexx][indexy][0];
		By = downstreamB[indexx][indexy][1];
		Bz = downstreamB[indexx][indexy][2];
	}

}

void move(double& x, double& y, double& z, double& vx, double& vy, double& vz, const double& Bx, const double& By, const double& Bz, const double& dt, double m){
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

void moveBoris(double& x, double& y, double& z, double& px, double& py, double& pz, const double& Bx, const double& By, const double& Bz, const double& Ex, const double& Ey, const double& Ez, const double& dt, double m){
	double qdt2 = q * dt * 0.5;

	double dpx = Ex*qdt2;
	double dpy = Ey*qdt2;
	double dpz = Ez*qdt2;


	px = px + dpx;
	py = py + dpy;
	pz = pz + dpz;

	double p2 = px * px + py * py + pz * pz;

	double gamma = sqrt(1 + p2/(m*m*c*c));

	double beta  = q*dt/(2*m*c);//todo???

	double betaShift = beta / gamma;
	double B2 = sqrt(Bx*Bx + By*By + Bz*Bz);
	double f = 2.0 / sqrt(1.0 + betaShift * betaShift * B2);

	double tempMomentumX = px + (py*Bz - pz*By)*betaShift;
	double tempMomentumY = py + (pz*Bx - px*Bz)*betaShift;
	double tempMomentumZ = pz + (px*By - py*Bx)*betaShift;


	double fBetaShift = f*betaShift;
	px = px + (tempMomentumY*Bz - tempMomentumZ*By)*fBetaShift + dpx;
	py = py + (tempMomentumZ*Bx - tempMomentumX*Bz)*fBetaShift + dpy;
	pz = pz + (tempMomentumX*By - tempMomentumY*Bx)*fBetaShift + dpz;

	p2 = px * px + py * py + pz * pz;
	gamma = sqrt(p2/(m*m*c*c) + 1);

	double vx = px/(m*gamma);
	double vy = py/(m*gamma);
	double vz = pz/(m*gamma);

	x = x + vx*dt;
	y = y + vy*dt;
	z = z + vz*dt;
}

int main(int argc, char** argv){
	const int Nt = 10000;
	const int chch = 10000;

	const int Nxmodes = 4;
	const int Nymodes = 4;
	double dt;
	
	double*** B;
	double*** E;
	double*** middleB;
	double*** middleE;
	double*** upstreamB;
	double*** upstreamE;
	double*** downstreamB;
	double*** downstreamE;

	double* juttnerValue;
	double* juttnerFunction;
	int juttnerN = 100;

	int Ny;
	int Nx;
	int middleNx;
	int upstreamNx;
	int downstreamNx;

	double vx, vy, vz;
	double x, y, z;
	double Bx, By, Bz;
	double B0x, B0y, B0z;

	B0x = 0;
	B0y = 0;
	B0z = Bmeansqr*sqrt(1.0 - turbulenceFraction);

	double v = c*sqrt(1.0 - 1.0/(gamma*gamma));
	double gammaFrame = 1.5;
	double vframe = c*sqrt(1.0 -1.0/(gammaFrame*gammaFrame));
	//double Bmean = sqrt((B0x*B0x + B0y*B0y + B0z*B0z)/(1.0 - turbulenceFraction));
	double omega = Bmeansqr*q/(gamma*massElectron*c);
	double rg = v/omega;
	dt = 0.1/omega;
	double dx = 1E6;

	turbulenceLength = 0.1*rg;

	Nx = 5000;
	Ny = 100;
	downstreamNx = 1000;
	middleNx = 1000;
	upstreamNx = 1000;
	int startDownstreamIndex = 100;
	int startMiddleIndex = startDownstreamIndex + downstreamNx;
	int startUpstreamIndex = startMiddleIndex + middleNx;


	B = new double**[Nx];
	E = new double**[Nx];
	for(int i = 0; i < Nx; ++i){
		B[i] = new double*[Ny];
		E[i] = new double*[Ny];
		for(int j = 0; j < Ny; ++j){
			B[i][j] = new double[3];
			E[i][j] = new double[3];
			for(int k = 0; k < 3; ++k){
				B[i][j][k] = 0;
				E[i][j][k] = 0;
			}
		}
	}

	printf("reading input\n");

	double fieldScale = 1.0;

	FILE* Bxfile = fopen("Bx.dat","r");
	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			fscanf(Bxfile,"%lf",&B[i][j][0]);
			B[i][j][0] *= fieldScale;
		}
	}
	fclose(Bxfile);
	FILE* Byfile = fopen("By.dat","r");
	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			fscanf(Byfile,"%lf",&B[i][j][1]);
			B[i][j][1] *= fieldScale;
		}
	}
	fclose(Byfile);
	FILE* Bzfile = fopen("Bz.dat","r");
	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			fscanf(Bzfile,"%lf",&B[i][j][2]);
			B[i][j][2] *= fieldScale;
		}
	}
	fclose(Bzfile);

	FILE* Exfile = fopen("Ex.dat","r");
	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			fscanf(Exfile,"%lf",&E[i][j][0]);
			E[i][j][0] *= fieldScale;
		}
	}
	fclose(Exfile);
	FILE* Eyfile = fopen("Ey.dat","r");
	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			fscanf(Eyfile,"%lf",&E[i][j][1]);
			E[i][j][1] *= fieldScale;
		}
	}
	fclose(Eyfile);
	FILE* Ezfile = fopen("Ez.dat","r");
	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			fscanf(Ezfile,"%lf",&E[i][j][2]);
			E[i][j][2] *= fieldScale;
		}
	}
	fclose(Ezfile);

	printf("initializing fields\n");

	upstreamB = new double**[upstreamNx];
	upstreamE = new double**[upstreamNx];
	for(int i = 0; i < upstreamNx; ++i){
		upstreamB[i] = new double*[Ny];
		upstreamE[i] = new double*[Ny];
		for(int j = 0; j < Ny; ++j){
			upstreamB[i][j] = new double[3];
			upstreamE[i][j] = new double[3];
			for(int k = 0; k < 3; ++k){
				upstreamB[i][j][k] = B[startUpstreamIndex + i][j][k];
				upstreamE[i][j][k] = E[startUpstreamIndex + i][j][k];
			}
		}
	}

	downstreamB = new double**[downstreamNx];
	downstreamE = new double**[downstreamNx];
	for(int i = 0; i < downstreamNx; ++i){
		downstreamB[i] = new double*[Ny];
		downstreamE[i] = new double*[Ny];
		for(int j = 0; j < Ny; ++j){
			downstreamB[i][j] = new double[3];
			downstreamE[i][j] = new double[3];
			for(int k = 0; k < 3; ++k){
				downstreamB[i][j][k] = B[startDownstreamIndex + i][j][k];
				downstreamE[i][j][k] = E[startDownstreamIndex + i][j][k];
			}
		}
	}

	middleB = new double**[middleNx];
	middleE = new double**[middleNx];
	for(int i = 0; i < middleNx; ++i){
		middleB[i] = new double*[Ny];
		middleE[i] = new double*[Ny];
		for(int j = 0; j < Ny; ++j){
			middleB[i][j] = new double[3];
			middleE[i][j] = new double[3];
			for(int k = 0; k < 3; ++k){
				middleB[i][j][k] = B[startMiddleIndex + i][j][k];
				middleE[i][j][k] = E[startMiddleIndex + i][j][k];
			}
		}
	}

	outputField("outBx.dat","outBy.dat","outBz.dat", downstreamB, middleB, upstreamB, downstreamNx, middleNx, upstreamNx, Ny);
	outputField("outEx.dat","outEy.dat","outEz.dat", downstreamE, middleE, upstreamE, downstreamNx, middleNx, upstreamNx, Ny);

	double* meanSqrX = new double[Nt];

	for(int i = 0; i < Nt; ++i) {
		meanSqrX[i] = 0;
	}

	double*** phases = new double**[Nxmodes+1];
	for(int i = 0; i <= Nxmodes; ++i) {
		phases[i] = new double*[Nymodes+1];
		for(int j = 0; j <= Nymodes; ++j) {
			phases[i][j] = new double[2];
			phases[i][j][0] = 0;
			phases[i][j][1] = 0;
		}
	}

	

	FILE* information = fopen("information.dat","w");
	fprintf(information, "rg = %g\n", rg);
	fprintf(information, "omega = %g\n", omega);
	fprintf(information, "dt = %g\n", dt);
	fprintf(information, "lambda/rg = %g\n", turbulenceLength/rg);
	fprintf(information, "randomSeed = %d\n", randomSeed);
	fclose(information);

	//FILE* out = fopen("trajectory.dat","w");
	int writeParameter = 100;

	srand(time(NULL));
	randomSeed = rand();
	//randomSeed = 4935;
	printf("random seed = %d\n",randomSeed);
	srand(randomSeed);

	for(int i = 0; i <= Nxmodes; ++i) {
		for(int j = 0; j <= Nymodes; ++j) {
			phases[i][j][0] = 2*pi*uniformDistribution();
			phases[i][j][1] = 2*pi*uniformDistribution();
		}
	}

	double temperature = 10*massElectron*c*c/kBoltzman;
	juttnerValue = new double[juttnerN];
	juttnerFunction = new double[juttnerN];
	for(int i = 0; i < juttnerN; ++i){
		juttnerValue[i] = 0;
		juttnerFunction[i] = 0;
	}
	evaluateJuttnerFunction(juttnerValue, juttnerFunction, temperature, massElectron, juttnerN);

	double** coordinates = new double*[chch];
	double** momentum = new double*[chch];
	for(int i = 0; i < chch; ++i){
		coordinates[i] = new double[3];
		momentum[i] = new double[3];
		for(int j = 0; j < 3; ++j){
			coordinates[i][j] = 0;
			momentum[i][j] = 0;
		}
	}

	FILE* outTrajectory = fopen("trajectories.dat", "w");
	for(int pcount = 0; pcount < chch; ++pcount){
		printf("particle %d\n", pcount);

		x = uniformDistribution()*middleNx*dx;
		y = uniformDistribution()*Ny*dx;
		z = 0;

		double px = 0;
		double py = 0;
		double pz = 0;

		createParticle(px, py, pz, temperature, massElectron, juttnerValue, juttnerFunction, juttnerN);

		double p2 = px*px + py*py + pz*pz;
		double gamma = sqrt(1 + p2/(massElectron*massElectron*c*c));

		vx = px/(gamma*massElectron);
		vy = py/(gamma*massElectron);
		vz = pz/(gamma*massElectron);

		coordinates[pcount][0] = x;
		coordinates[pcount][1] = y;
		coordinates[pcount][2] = z;
		momentum[pcount][0] = px;
		momentum[pcount][1] = py;
		momentum[pcount][2] = pz;
	}

	int currentWriteNumber = 0;
	
	for(int i = 0; i < Nt; ++i){
		printf("iteration %d\n", i);
		if(i%writeParameter == 0){
			printf("outputing %d\n",currentWriteNumber);
			std::string fileNumber = std::string("_") + convertIntToString(currentWriteNumber);
			outputDistribution(("distribution" + fileNumber + ".dat").c_str(), momentum, chch);
			currentWriteNumber++;
		}
		for(int pcount = 0; pcount < chch; ++pcount){

			x = coordinates[pcount][0];
			y = coordinates[pcount][1];
			z = coordinates[pcount][2];
			double px = momentum[pcount][0];
			double py = momentum[pcount][1];
			double pz = momentum[pcount][2];

			meanSqrX[i] += x*x/chch;
			//if(i%1000 == 0){
			//	printf("interation %d\n", i);
			//}
			v = c*sqrt(1.0 - 1.0/(gamma*gamma));
			if(i%writeParameter == 0){
				double xshift = (x - vframe*i*dt)*gammaFrame;
				fprintf(outTrajectory, "%g ", x);
			}
			double Ex, Ey, Ez;
			//getBfield(x, y, z, Bx, By, Bz, B0x, B0y, B0z, Nxmodes, Nymodes, phases);
			getEandBfield(x, y, z, Ex, Ey, Ez, Bx, By, Bz, downstreamE, downstreamB, middleE, middleB, upstreamE, upstreamB, downstreamNx, middleNx, upstreamNx, Ny, dx);
			//move(x, y, z, vx, vy, vz, Bx, By, Bz, dt, massElectron);

			moveBoris(x, y, z, px, py, pz, Bx, By, Bz, Ex, Ey, Ez, dt,massElectron);

			coordinates[pcount][0] = x;
			coordinates[pcount][1] = y;
			coordinates[pcount][2] = z;
			momentum[pcount][0] = px;
			momentum[pcount][1] = py;
			momentum[pcount][2] = pz;
		}
		fprintf(outTrajectory, "\n");
	}

	fclose(outTrajectory);
	//fclose(out);

	FILE* out = fopen("meanx.dat","w");
	for(int i = 0; i < chch; ++i) {
		fprintf(out, "%g %g\n", dt*i, meanSqrX[i]);
	}
	fclose(out);
	delete[] meanSqrX;

	for(int i = 0;i < chch; ++i){
		delete[] coordinates[i];
		delete[] momentum[i];
	}
	delete[] coordinates;
	delete[] momentum;

	for(int i = 0; i <= Nxmodes; ++i) {
		for(int j = 0; j <= Nymodes; ++j) {
			delete[] phases[i][j];
		}
		delete[] phases[i];
	}
	delete[] phases;

	delete[] juttnerValue;
	delete[] juttnerFunction;

	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			delete[] B[i][j];
			delete[] E[i][j];
		}
		delete[] B[i];
		delete[] E[i];
	}
	delete[] B;
	delete[] E;

	for(int i = 0; i < upstreamNx; ++i){
		for(int j = 0; j < Ny; ++j){
			delete[] upstreamB[i][j];
			delete[] upstreamE[i][j];
		}
		delete[] upstreamB[i];
		delete[] upstreamE[i];
	}
	delete[] upstreamB;
	delete[] upstreamE;

	for(int i = 0; i < downstreamNx; ++i){
		for(int j = 0; j < Ny; ++j){
			delete[] downstreamB[i][j];
			delete[] downstreamE[i][j];
		}
		delete[] downstreamB[i];
		delete[] downstreamE[i];
	}
	delete[] downstreamB;
	delete[] downstreamE;

	for(int i = 0; i < middleNx; ++i){
		for(int j = 0; j < Ny; ++j){
			delete[] middleB[i][j];
			delete[] middleE[i][j];
		}
		delete[] middleB[i];
		delete[] middleE[i];
	}
	delete[] middleB;
	delete[] middleE;

	/*FILE* field = fopen("field.dat","w");

	double dx = rg/10;
	int Nx = 500;
	int Ny = 500;
	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j< Ny; ++j){
			x = i*dx;
			y = j*dx;
			z = 0;
			getBfield(x, y, z, Bx, By, Bz, B0x, B0y, B0z);
			fprintf(field, "%g %g %g\n", Bx, By, Bz);
		}
	}

	fclose(field);*/

	return 0;
}