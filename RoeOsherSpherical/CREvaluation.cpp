#include <time.h>
#include "math.h"
#include "stdio.h"
#include <stdlib.h>
#include <omp.h>
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"
#include "progon.h"

void Simulation::updateDiffusionCoef(){
	int i;
//#pragma omp parallel for private(i)
		for(i = 0; i < rgridNumber; ++i){
			int prevK = 0;
			for(int j = pgridNumber-1; j >= 0; --j){
				double p = pgrid[j];
				double B = B0;
				/*for(int k = 0; k < kgridNumber; ++k){
					if(kgrid[k] > electron_charge*largeScaleField[i][k]/(speed_of_light*p) ){
						B = largeScaleField[i][k];
						prevK = k;
						break;
					}
				}*/
				B = largeScaleField[i][kgridNumber-1];
				double coef = p*speed_of_light*speed_of_light/(3*electron_charge*B);
				double dx = deltaR[i];
				double lambda = coef/speed_of_light;
				diffusionCoef[i][j] = coef;
			}
		
	}
}

//инжекционный член
double Simulation::injection(int i){
	double pf = pgrid[injectionMomentum];
	double dp = (pgrid[injectionMomentum + 1] - pgrid[injectionMomentum - 1])/2;
	//double xi = 5;
	double xi = pgrid[injectionMomentum]*speed_of_light/(kBoltzman*temperatureIn(i+1));
	double eta = cube(xi)*exp(-xi);
    return (1E-3)*middleDensity[i+1]*abs2(middleVelocity[i-1])*pf/(massProton*dp*deltaR[i]);
    return (1E-1)*middleDensity[i+1]*abs2(middleVelocity[i-1])*pf/(massProton*dp*deltaR[i]);
    //return (1E-3)*pf*pf*pf*middleDensity[i+1]/(massProton*dp*volume(i)*deltaT);
}


//расчет космических лучей

void Simulation::evaluateCR(){
	//printf("solve CR\n");
	/*if(shockWavePoint > 0 && shockWavePoint < rgridNumber){
		distributionFunction[shockWavePoint][injectionMomentum] += injection()*deltaT;
	}*/
	double mc2 = massProton*sqr(speed_of_light);
	double* upper = new double[rgridNumber + 1];
	double* middle = new double[rgridNumber + 1];
	double* lower = new double[rgridNumber + 1];

	double* f = new double[rgridNumber + 1];
	double* x = new double[rgridNumber + 1];

	double* alpha = new double[rgridNumber];
	double* beta = new double[rgridNumber];
	int k;
	//todo?? private pointer?
#pragma omp parallel for private(k, upper, middle, lower, f, x, alpha, beta)
		for(k = 0; k < pgridNumber; ++k){


			double y = logPgrid[k];
			double p = pgrid[k];
			double gkp = distributionFunction[0][k];
			double gkm=0.0;
			if (k == 0) {
				gkm=0.0;
			} else {
				gkm = distributionFunction[0][k-1];
			}
			double gkpp;
			if(k < pgridNumber-1){
				gkpp = distributionFunction[0][k+1];
			} else {
				gkpp = 0;
			}
			//double dx = (grid[1] + upstreamR/2)/2;
			double dxp = grid[1] - grid[0];
			//double dxm = grid[0] + upstreamR/2;
			double dxm = grid[1] - grid[0];
			double xp = (grid[1] + grid[0])/2;
			//double xm = (grid[0] - upstreamR/2)/2;
			double xm = grid[0];
			double gip = (distributionFunction[1][k] + distributionFunction[0][k])/2;
			double gim = distributionFunction[0][k]/2;
			double dV;
			middle[0] = 1;// + (deltaT/(2*dx))*(diffusionCoef(0,p)/dxp + diffusionCoef(0,p)/dxm);
			upper[0] = 0;// -(deltaT/(2*dx))*(diffusionCoef(0,p)/dxp);
			f[0]=distributionFunction[0][k];// + (deltaT/(2*dx))*(diffusionCoef(0,p)*distributionFunction[1][k]/dxp)
						 // - (deltaT/(2*dx))*(diffusionCoef(0,p)/dxp+diffusionCoef(0,p)/dxm)*distributionFunction[0][k] 
						 // - (deltaT/dx)*(middleVelocity[1]*gip - middleVelocity[0]*gim)
						 // + (deltaT/3)*((middleVelocity[1] - middleVelocity[0])/dx)*((gkp-gkm)/deltaLogP);
			for(int i = 1; i < rgridNumber-1; ++i){
				gkp = distributionFunction[i][k];
				if (k == 0) {
					gkm=0;
				} else {
					gkm=tempDistributionFunction[i][k-1];
				}
				if(k < pgridNumber-1){
					gkpp = distributionFunction[i][k+1];
				} else {
					gkpp = 0;
				}
				//dx = (grid[i+1] - grid[i-1])/2;
				dxp=grid[i+1]-grid[i];
				dxm=grid[i]-grid[i-1];
				xp=(grid[i+1]+grid[i])/2;
				xm=(grid[i]+grid[i-1])/2;
				double localx = grid[i];
				dV = (xp*xp*xp - xm*xm*xm)/3;
				gip=(distributionFunction[i+1][k] + distributionFunction[i][k])/2;
				gim=(distributionFunction[i][k] + distributionFunction[i-1][k])/2;
				lower[i-1] = -(deltaT/(2*dV))*(xm*xm*diffusionCoef[i-1][k]/dxm);
				middle[i] = 1 + (deltaT/(2*dV))*(xp*xp*diffusionCoef[i-1][k]/dxp + xm*xm*diffusionCoef[i][k]/dxm);
				upper[i] = -(deltaT/(2*dV))*(xp*xp*diffusionCoef[i][k]/dxp);
				//f[i] = distributionFunction[i][k] + (deltaT/(2*dV))*(xp*xp*diffusionCoef[i][k]*(distributionFunction[i+1][k] - distributionFunction[i][k])/dxp
								//- xm*xm*diffusionCoef[i-1][k]*(distributionFunction[i][k] - distributionFunction[i-1][k])/dxm)
								//- (deltaT/dV)*(xp*xp*middleVelocity[i]*distributionFunction[i][k] - xm*xm*middleVelocity[i-1]*distributionFunction[i-1][k]);
				if(i == shockWavePoint){
					//double v2 = middleVelocity[i+1] + vscattering[i+1];
					double v2 = middleVelocity[i+1];
					//double v1 = middleVelocity[i] + vscattering[i];
					double v1 = middleVelocity[i];
					f[i] = distributionFunction[i][k]  + deltaT*((1/(2*dV))*(xp*xp*diffusionCoef[i][k]*(distributionFunction[i+1][k] - distributionFunction[i][k])/dxp - xm*xm*diffusionCoef[i-1][k]*(distributionFunction[i][k] - distributionFunction[i-1][k])/dxm) - (1/dV)*localx*localx*distributionFunction[i][k]*(v2 - v1));
				} else if(i == shockWavePoint-1){
					//double v2 = middleVelocity[i+1] + vscattering[i+1];
					double v2 = middleVelocity[i+1];
					//double v2 = middleVelocity[i+1] + vscattering[i+1];
					double v1 = middleVelocity[i];
					f[i] = distributionFunction[i][k]  + deltaT*((1/(2*dV))*(xp*xp*diffusionCoef[i][k]*(distributionFunction[i+1][k] - distributionFunction[i][k])/dxp - xm*xm*diffusionCoef[i-1][k]*(distributionFunction[i][k] - distributionFunction[i-1][k])/dxm) - (1/dV)*(xp*xp*v2*distributionFunction[i+1][k] - xm*xm*v1*distributionFunction[i-1][k]));
				} else {
					//double v2 = middleVelocity[i+1] + vscattering[i+1];
					double v2 = middleVelocity[i+1];
					//double v1 = middleVelocity[i-1] + vscattering[i-1];
					double v1 = middleVelocity[i-1];
					f[i] = distributionFunction[i][k]  + deltaT*((1/(2*dV))*(xp*xp*diffusionCoef[i][k]*(distributionFunction[i+1][k] - distributionFunction[i][k])/dxp - xm*xm*diffusionCoef[i-1][k]*(distributionFunction[i][k] - distributionFunction[i-1][k])/dxm) - (1/dV)*0.5*(gridsquare[i+1]*v2*distributionFunction[i+1][k] - gridsquare[i-1]*v1*distributionFunction[i-1][k]));
				}
				//f[i] = distributionFunction[i][k]  + deltaT*((1/(2*dV))*(xp*xp*diffusionCoef[i][k]*(distributionFunction[i+1][k] - distributionFunction[i][k])/dxp - xm*xm*diffusionCoef[i-1][k]*(distributionFunction[i][k] - distributionFunction[i-1][k])/dxm) - (1/dV)*(middleVelocity[i-1]+vscattering[i-1])*0.5*(distributionFunction[i+1][k] - distributionFunction[i-1][k]));
				if(f[i] < 0){
					//printf("f[i] < 0\n");
					f[i] = 0;
				}
				if((xp*xp*middleVelocity[i+1] - xm*xm*middleVelocity[i]) < 0){
					//double v2 = middleVelocity[i+1] + vscattering[i+1];
					double v2 = middleVelocity[i+1];
					//double v1 = middleVelocity[i] + vscattering[i];
					double v1 = middleVelocity[i];
					if(gkp - gkm < 0)
						//f[i] += (deltaT/3)*((xp*xp*middleVelocity[i] - xm*xm*middleVelocity[i-1])/dV)*((gkp - gkm)/deltaLogP);
						f[i] += -(deltaT/3)*((xp*xp*v2 - xm*xm*v1)/dV)*((gkm)/deltaLogP);
						middle[i] += -(deltaT/3)*((xp*xp*v2 - xm*xm*v1)/dV)/deltaLogP;
				} else {
					//double v2 = middleVelocity[i+1] + vscattering[i+1];
					double v2 = middleVelocity[i+1];
					//double v1 = middleVelocity[i] + vscattering[i];
					double v1 = middleVelocity[i];
					if(gkpp - gkp > 0)
						f[i] += (deltaT/3)*((xp*xp*v2 - xm*xm*v1)/dV)*((gkpp - gkp)/deltaLogP);
				}
				//f[i] = distributionFunction[i][k];
                if(abs2(i - shockWavePoint)<2 && abs2(k - injectionMomentum) < 1){
					double inj = injection(i);
					double E = sqrt(sqr(mc2) + sqr(pgrid[injectionMomentum])*speed_of_light) - mc2;
					double dE = deltaT*inj*E*deltaLogP;
					if(dE > tempEnergy[i]){
						printf("dE < tempEnergy[i]\n");
						inj *= 0.5*tempEnergy[i]/dE;
					}
					f[i] += deltaT*inj;
					//todo shift volume to 1/2
					injectedParticles += inj*deltaT*volume(i)*deltaLogP;
					//tempDensity[i] -= deltaT*inj*massProton*deltaLogP;
					//tempMomentum[i] -= deltaT*inj*massProton*deltaLogP*middleVelocity[i];
					//tempEnergy[i] -= deltaT*inj*E*deltaLogP;
					injectedEnergy += deltaT*inj*E*deltaLogP*volume(i);

					/*double inj = injection(i);
					double E = sqrt(sqr(mc2) + sqr(pgrid[injectionMomentum])*speed_of_light) - mc2;
					double dE = deltaT*inj*E*deltaLogP;
					if(dE > tempEnergy[i]){
						printf("dE < tempEnergy[i]\n");
						inj *= 0.5*tempEnergy[i]/dE;
					}
					//f[i] += deltaT*inj;
					f[i] += deltaT*inj/cube(pgrid[injectionMomentum]);
					//todo shift volume to 1/2
					//injectedParticles += inj*deltaT*volume(i)*deltaLogP;
					injectedParticles += inj*deltaT*volume(i)*deltaLogP/cube(pgrid[injectionMomentum]);
					tempDensity[i] -= deltaT*inj*massProton*deltaLogP;
					tempMomentum[i] -= deltaT*inj*massProton*deltaLogP*middleVelocity[i];
					tempEnergy[i] -= deltaT*inj*E*deltaLogP;
					injectedEnergy += deltaT*inj*E*deltaLogP*volume(i);*/
					if(tempDensity[i] < 0){
						printf("tempDensity[i] < 0 by CR\n");
						exit(0);
					}
					if(tempEnergy[i] < 0){
						printf("tempEnergy[i] < 0 by CR\n");
						exit(0);
					}
					alertNaNOrInfinity(tempDensity[i], "density = NaN");
					alertNaNOrInfinity(tempMomentum[i], "momentum = NaN");
					alertNaNOrInfinity(tempEnergy[i], "energy = NaN");
				}
			}
			gkp = distributionFunction[rgridNumber-1][k];
			if (k==0) {
				gkm=0;
			} else {
				gkm = distributionFunction[rgridNumber-1][k-1];
			}
			//dx = (grid[rgridNumber-1] - grid[rgridNumber - 2])/2;
			dxm = (grid[rgridNumber-1] - grid[rgridNumber - 2]);
			xp = grid[rgridNumber-1];
			xm = (grid[rgridNumber-1] + grid[rgridNumber-2])/2;
			dV = (xp*xp*xp - xm*xm*xm)/3;
			gip = distributionFunction[rgridNumber-1][k];
			gim = (distributionFunction[rgridNumber-1][k] + distributionFunction[rgridNumber - 2][k])/2;
			lower[rgridNumber-2] = -(deltaT/(2*dV))*(xp*xp*diffusionCoef[rgridNumber-1][k]/dxp);
			middle[rgridNumber-1] = 1 + (deltaT/(2*dV))*(xm*xm*diffusionCoef[rgridNumber-1][k]/dxm);
			f[rgridNumber-1] = distributionFunction[rgridNumber-1][k] - (deltaT/(2*dV))*(xm*xm*diffusionCoef[rgridNumber-2][k]*(distributionFunction[rgridNumber-1][k] - distributionFunction[rgridNumber - 2][k])/dxm)
								- (deltaT/dV)*grid[rgridNumber-1]*grid[rgridNumber-1]*(middleVelocity[rgridNumber-1]*gip - middleVelocity[rgridNumber-2]*gim)
								 + (deltaT/3)*((xp*xp*middleVelocity[rgridNumber-1] - xm*xm*middleVelocity[rgridNumber-2])/dV)*((gkp - gkm)/deltaLogP);
			progon(lower,middle, upper,rgridNumber-1,f,x, alpha, beta);

			for(int i = 0; i < rgridNumber; ++i){
				//alertNegative(x[i],"tempDistribution < 0");
				alertNaNOrInfinity(x[i],"tempDistribution <= NaN");
				tempDistributionFunction[i][k]= x[i];
				if(x[i] < 0){
					tempDistributionFunction[i][k] = 0;
                    /*if(abs2(x[i]) > 1E-50){
						printf("distribution[i] < 0\n");
					}*/
				}
			}
			//tempDistributionFunction[rgridNumber][k] =x[rgridNumber-1];
			tempDistributionFunction[rgridNumber][k] = 0;
		
	}

	delete[] upper;
	delete[] middle;
	delete[] lower;
	delete[] f;
	delete[] x;
	delete[] alpha;
	delete[] beta;
	
}

//решение трёх диагональной матрицы
void Simulation::solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x, double* alpha, double* beta){

	alpha[1] = -upper[0]/middle[0];
	beta[1] = f[0]/middle[0];
	for(int i = 2; i < rgridNumber; ++i){
		double temp = lower[i-1]*alpha[i-1] + middle[i-1];
		alpha[i] = -upper[i-1]/temp;
		beta[i] = (f[i-1] - lower[i-1]*beta[i-1])/temp;
	}

	x[rgridNumber - 1] = (f[rgridNumber-1] - lower[rgridNumber-1]*beta[rgridNumber-1])/(lower[rgridNumber-1]*alpha[rgridNumber-1] + middle[rgridNumber-1]);
	//alertNaNOrInfinity(x[rgridNumber-1],"x = NaN");
	//alertNegative(x[rgridNumber-1],"x < 0");

	for(int i = rgridNumber - 2; i >= 0; --i){
		x[i] = alpha[i+1]*x[i+1] + beta[i+1];
		//alertNaNOrInfinity(x[i],"x = NaN");
		//alertNegative(x[i],"x < 0");
	}
}


//вычисление давления космических лучей

void Simulation::evaluateCosmicRayPressure(){
	double* partPressure = new double[pgridNumber];
	double deltaLogP = logPgrid[1] - logPgrid[0];
	for(int j = 0; j < pgridNumber; ++j){
		double p = pgrid[j];
		double v = p/sqrt(massProton*massProton + p*p/(speed_of_light*speed_of_light));
		partPressure[j] = p*v*deltaLogP/3;
		//partPressure[j] = p*v*deltaLogP*cube(p);
	}
	for(int i = 0; i < rgridNumber; ++i){
		double pressure = 0;
		double concentration = 0;
		for(int j = 0; j < pgridNumber; ++j){
			if(j >= goodMomentum){
				pressure += distributionFunction[i][j]*partPressure[j];
			}
			concentration += distributionFunction[i][j]*deltaLogP;
			//concentration += distributionFunction[i][j]*deltaLogP*cube(pgrid[i]);
		}
		//4pi?
		cosmicRayPressure[i] = pressure;
		cosmicRayConcentration[i] = concentration;
	}

	delete[] partPressure;
}
