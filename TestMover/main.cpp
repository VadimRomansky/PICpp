#include "stdio.h"
#include <stdlib.h>
#include <time.h>
#include "math.h"
#include <omp.h>
#include <string>

#include "constants.h"
#include "output.h"
#include "distribution.h"

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
	double turbulenceLength = 1.0;////todo!!
	double turbulenceFraction = 0.9;

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


void getDipoleField(const double& x, const double& y, const double& z, double& Bx, double& By, double& Bz, const double& mux, const double& muy, const double& muz){
	double r = sqrt(x*x + y*y + z*z);
	double nx = x/r;
	double ny = y/r;
	double nz = z/r;

	double mun = nx*mux + ny*muy + nz*muz;

	double r3 = r*r*r;
	
	Bx = (3*mun*nx - mux)/r3;
	By = (3*mun*ny - muy)/r3;
	Bz = (3*mun*nz - muz)/r3;
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

	/*Ex = 0.0;
	Ey = 0;
	Ez = 0;
	Bx = 0;
	By = 0;
	Bz = 1.0;*/
}

void move(double& x, double& y, double& z, double& px, double& py, double& pz, const double& Bx, const double& By, const double& Bz, const double& dt, double m){
	double	Bnorm = sqrt(Bx*Bx + By*By + Bz*Bz);
	double p2 = px*px + py*py + pz*pz;
	double gamma = sqrt(1.0 + p2/(m*m*c*c));
	double vx = px/(m*gamma);
	double vy = py/(m*gamma);
	double vz = pz/(m*gamma);
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

		double dphi = -electron_charge*Bnorm*dt/(gamma*m*c);
				
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

	px = vx*m*gamma;
	py = vy*m*gamma;
	pz = vz*m*gamma;
}

void moveBoris(double& x, double& y, double& z, double& px, double& py, double& pz, const double& Bx, const double& By, const double& Bz, const double& Ex, const double& Ey, const double& Ez, const double& dt, double m){
	
	double qdt2 = -electron_charge * dt * 0.5;

	double dpx = Ex*qdt2;
	double dpy = Ey*qdt2;
	double dpz = Ez*qdt2;


	px = px + dpx;
	py = py + dpy;
	pz = pz + dpz;

	double p2 = px * px + py * py + pz * pz;

	double gamma = sqrt(1.0 + p2/(m*m*c*c));

	double beta  = -electron_charge*dt/(2*m*c);//todo???

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
	gamma = sqrt(p2/(m*m*c*c) + 1.0);

	double vx = px/(m*gamma);
	double vy = py/(m*gamma);
	double vz = pz/(m*gamma);

	x = x + vx*dt;
	y = y + vy*dt;
	z = z + vz*dt;
}

void LorentzTransformationFields(double*** E, double*** B, double u, int Nx, int Ny){
	double gamma = 1.0/sqrt(1 - u*u/(c*c));
	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			double tempBx = B[i][j][0];
			double tempBy = gamma*(B[i][j][1] + u*E[i][j][2]/c);
			double tempBz = gamma*(B[i][j][2] - u*E[i][j][1]/c);

			double tempEx = E[i][j][0];
			double tempEy = gamma*(E[i][j][1] - u*B[i][j][2]/c);
			double tempEz = gamma*(E[i][j][2] + u*B[i][j][1]/c);

			B[i][j][0] = tempBx;
			B[i][j][1] = tempBy;
			B[i][j][2] = tempBz;

			E[i][j][0] = tempEx;
			E[i][j][1] = tempEy;
			E[i][j][2] = tempEz;
		}
	}
}

int main(int argc, char** argv){
	//omp_set_num_threads(28);
	const int Nt = 100000;
	const int chch = 500;
	const int writeParameter = 10000;

	const int Nxmodes = 4;
	const int Nymodes = 4;
	
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

	double sigma = 4.0;
	double gammaFrame = 1.5;
	double n = 1;
	double ntristan = 2;//0.5*ppc0
	double ctristan = 0.45;
	double comp = 5;
	double omp = ctristan/comp;
	double metristan = omp*omp*gammaFrame/(ntristan*(1 + massElectron/massProtonReal));
	double omega_pe = sqrt(4*pi*n*electron_charge*electron_charge/(gammaFrame*massElectron));
	double vframe = c*sqrt(1.0 - 1.0/(gammaFrame*gammaFrame));
	double Bmeansqr = sqrt(gammaFrame*n*(1+massElectron/massProtonReal)*c*c*(massElectron)*sigma);
	double omega_ce = Bmeansqr*electron_charge/(gammaFrame*massElectron*c);
	double fieldScale = sqrt(4*pi*(n/ntristan)*(massElectron/metristan)*c*c/(ctristan*ctristan));

	int sampling = 20;
	double dx = 0.2*c/omega_pe;

	double dt = 0.1*0.2/omega_pe;
	dt = 1E-5;

	printf("start\n");
	Nx = 5000;
	Ny = 100;
	downstreamNx = 300;
	middleNx = 1000;
	upstreamNx = 2000;
	int startDownstreamIndex = 100;
	int startMiddleIndex = startDownstreamIndex + downstreamNx;
	int startUpstreamIndex = startMiddleIndex + middleNx;

	printf("init array Nx = %d Ny = %d\n", Nx, Ny);


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

	LorentzTransformationFields(E, B, 0.18*c, Nx, Ny);

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

	outputField("./output/outBx.dat","./output/outBy.dat","./output/outBz.dat", downstreamB, middleB, upstreamB, downstreamNx, middleNx, upstreamNx, Ny);
	outputField("./output/outEx.dat","./output/outEy.dat","./output/outEz.dat", downstreamE, middleE, upstreamE, downstreamNx, middleNx, upstreamNx, Ny);

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

	



	//FILE* out = fopen("trajectory.dat","w");

	srand(time(NULL));
	randomSeed = rand();
	randomSeed = 4935;
	srand(randomSeed);
	FILE* information = fopen("./output/information.dat","w");
	fprintf(information, "dt = %g\n", dt);
	fprintf(information, "randomSeed = %d\n", randomSeed);
	fclose(information);
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

	FILE* outTrajectory = fopen("./output/trajectories.dat", "w");
	for(int pcount = 0; pcount < chch; ++pcount){
		if(pcount%100 == 0){
			printf("particle %d\n", pcount);
		}

		x = uniformDistribution()*middleNx*dx*sampling;
		y = uniformDistribution()*Ny*dx*sampling;
		x = 1E10;
		y = 0;
		z = 0;

		double px = 0;
		double py = 0;
		double pz = 0;

		//createParticle(px, py, pz, temperature, massElectron, juttnerValue, juttnerFunction, juttnerN);
		createFastParticle(px,py,pz, massElectron, 2.0);
		pz = px;

		double p2 = px*px + py*py + pz*pz;
		double gamma = sqrt(1.0 + p2/(massElectron*massElectron*c*c));

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

	const int partWrite = 10;
	int writePartNumber = 50;
	int numbers[partWrite];
	numbers[0] = 0;
	numbers[1] = 1;
	numbers[2] = 2;
	numbers[3] = 3;
	numbers[4] = 4;
	numbers[5] = 5;
	numbers[6] = 6;
	numbers[7] = 7;
	numbers[8] = 8;
	numbers[9] = 9;

	bool readNumbers = false;
	if(readNumbers){
		FILE* file = fopen("numbers_10.dat","r");
		int n;
		double p;
		for(int i = 0; i < partWrite; ++i){
			fscanf(file, "%d %lf", &n, &p); 
			numbers[i] = n;
		}
		fclose(file);
	}

	for(int k = 0; k < partWrite; ++k){
		std::string fileNumber = std::string("_") + convertIntToString(k);
		FILE* file = fopen(("./output/trajectory" + fileNumber + ".dat").c_str(),"w");
		fclose(file);
	}
	
	for(int i = 0; i < Nt; ++i){
		printf("iteration %d\n", i);
		if(i%writeParameter == 0){
			printf("outputing %d\n",currentWriteNumber);
			std::string fileNumber = std::string("_") + convertIntToString(currentWriteNumber);
			outputDistribution(("./output/distribution" + fileNumber + ".dat").c_str(), momentum, chch);
			currentWriteNumber++;
			writeMaxParticles(("./output/numbers" + fileNumber + ".dat").c_str(), momentum, chch);
		}
		int pcount = 0;
		#pragma omp parallel for private(pcount)
		for(pcount = 0; pcount < chch; ++pcount){

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
			
			double Ex, Ey, Ez;
			//getBfield(x, y, z, Bx, By, Bz, B0x, B0y, B0z, Nxmodes, Nymodes, phases);
			//getEandBfield(x, y, z, Ex, Ey, Ez, Bx, By, Bz, downstreamE, downstreamB, middleE, middleB, upstreamE, upstreamB, downstreamNx, middleNx, upstreamNx, Ny, dx*sampling);
			getDipoleField(x, y, z, Bx, By, Bz, 0, 0, 1E27);

			double p2 = px*px + py*py + pz*pz;
			double gamma = sqrt(1.0 + p2/(massElectron*massElectron*c*c));
			double p = sqrt(p2);
			double B = sqrt(Bx*Bx + By*By + Bz*Bz);

			double rg = p*c/(electron_charge*B);
			double r = sqrt(x*x + y*y + z*z);
			double relation = rg/r;
			double omega = electron_charge*B/(gamma*massElectron*c);
			double timeRelation = dt*omega;

			if(i %writePartNumber == 0){
				for(int k = 0; k < partWrite; ++k){
					if(pcount == numbers[k]){
						std::string fileNumber = std::string("_") + convertIntToString(k);
						FILE* file = fopen(("./output/trajectory" + fileNumber + ".dat").c_str(),"a");
						double pb = px*Bx + py*By + pz*Bz;
						double adiabaticInvariant = (p2 - (pb*pb/(B*B)))/B;
						fprintf(file, "%d %g %g %g %g %g %g %g %g %g %g\n", i, x, y, z, px, py, pz, gamma, p2, B, adiabaticInvariant);
						fclose(file);
					}
				}
			}
			if(i%writeParameter == 0){
				double xshift = (x - vframe*i*dt)*gammaFrame;
				fprintf(outTrajectory, "%g ", x);
			}

			//moveBoris(x, y, z, px, py, pz, Bx, By, Bz, Ex, Ey, Ez, dt,massElectron);
			move(x, y, z, px, py, pz, Bx, By, Bz, dt, massElectron);
			Ex = 0;
			Ey = 0;
			Ez = 0;
			/*if(x < -downstreamNx*dx){
				px = fabs(px);
				x = -2*downstreamNx*dx - x;
			}*/

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

	FILE* out = fopen("./output/meanx.dat","w");
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