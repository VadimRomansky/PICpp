#include <list>
#include <time.h>
#include "math.h"
#include "stdio.h"
#include <stdlib.h>
#include <omp.h>
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

//конструктор
Simulation::Simulation(){
	serialized = false;
	initialEnergy = 10E50;
	myTime = 0;
	tracPen = true;
	shockWavePoint = -1;
	shockWaveMoved = false;
	injectedParticles = 0;
	totalEnergy = 0;
	totalKineticEnergy = 0;
	totalTermalEnergy = 0;
	totalMagneticEnergy = 0;
	totalParticleEnergy = 0;
	totalMomentum = 0;
	totalParticles = 0;
	currentIteration = 0;
	injectedEnergy = 0;
	uGradPEnergy = 0;
	stopAmplification = false;
}

//деструктор
Simulation::~Simulation(){
	delete[] pgrid;
	delete[] logPgrid;
	delete[] grid;
	delete[] middleGrid;
	delete[] deltaR;
	delete[] middleDeltaR;
	delete[] tempGrid;
	delete[] gridsquare;
	delete[] pointDensity;
	delete[] pointVelocity;
	delete[] pointEnthalpy;
	delete[] pointSoundSpeed;
	delete[] middleDensity;
	delete[] middleVelocity;
	delete[] middlePressure;
	delete[] cosmicRayPressure;
	delete[] cosmicRayConcentration;
	delete[] tempDensity;
	delete[] tempMomentum;
	delete[] tempEnergy;

	for(int i = 0; i < rgridNumber; ++i){
		delete[] dFluxPlus[i];
		delete[] dFluxMinus[i];
		delete[] mFluxPlus[i];
		delete[] mFluxMinus[i];
		delete[] eFluxPlus[i];
		delete[] eFluxMinus[i];
	}

	delete[] dFlux;
	delete[] mFlux;
	delete[] eFlux;
	delete[] dFluxPlus;
	delete[] mFluxPlus;
	delete[] eFluxPlus;
	delete[] dFluxMinus;
	delete[] mFluxMinus;
	delete[] eFluxMinus;

	delete[] kgrid;
	delete[] logKgrid;

	delete[] tempU;

	for(int i = 0; i <= rgridNumber; ++i){
		delete[] distributionFunction[i];
		delete[] tempDistributionFunction[i];
		delete[] largeScaleField[i];
		delete[] magneticField[i];
		delete[] tempMagneticField[i];
		delete[] crflux[i];
		delete[] growth_rate[i];
		delete[] diffusionCoef[i];
	}
	delete[] magneticEnergy;

	delete[] integratedFlux;
	delete[] magneticInductionSum;
	delete[] tempMagneticField;
	delete[] largeScaleField;
	delete[] crflux;
	delete[] growth_rate;
	delete[] magneticField;
	delete[] distributionFunction;
	delete[] tempDistributionFunction;
	delete[] diffusionCoef;
}

//инициализация профиля после считывания данных

void Simulation::initializeProfile(){
	downstreamR = 0;
	//upstreamR = 20000;

	prevShockWavePoint = -1;
	shockWavePoint = -1;

	initializeArrays();

	minP = massProton*speed_of_light;
	maxP = minP*1E8;

	double r = downstreamR;
	double pressure0 = density0*kBoltzman*temperature/massProton;

	//minP = massProton*speed_of_light/10;


	//deltaR0 = (upstreamR - downstreamR)/rgridNumber;
	//grid[0] = -upstreamR/2 + deltaR0;
	//for(int i = 1; i <= rgridNumber; ++i){
		//grid[i] = grid[i-1] + deltaR0;
	//}

	int leftNumber = rgridNumber/5;
	double shockWaveR = upstreamR/20;
	double a = shockWaveR;
	double b = upstreamR - shockWaveR;
	double R1 = a/10;
	double R2 = b/100;
	double h1=leftNumber/log(1.0+a/R1);
	double h2=(rgridNumber + 1 - leftNumber)/log(1.0+b/R2);
	grid[0]=0;
	for(int i = 1; i < leftNumber; ++ i){ 
		grid[i] = R1*(1 - exp(-(1.0*(i+1)-leftNumber)/h1)) + shockWaveR;
	}
	for(int i = leftNumber; i < rgridNumber; ++i){
		grid[i] = R2*(exp((1.0*(i+1)- leftNumber)/h2)-1.0) + shockWaveR;
	}
	grid[rgridNumber] = upstreamR;

	for(int i = 1; i <= rgridNumber; ++i){
		if(grid[i] <= grid[i-1]){
			printf("grid[i] <= grid[i-1]\n");
			exit(0);
		}
	}

	double mindR = upstreamR;
	for(int i = 0; i < rgridNumber; ++i){
		middleGrid[i] = (grid[i] + grid[i+1])/2;
		tempGrid[i] = grid[i];
		deltaR[i] = grid[i+1] - grid[i];
		if(deltaR[i] < mindR){
			mindR = deltaR[i];
		}
		gridsquare[i] = sqr(grid[i]);
	}
	tempGrid[rgridNumber] = grid[rgridNumber];

	minK = (1E-6)*electron_charge*B0/(speed_of_light*maxP);
	minK = max2(minK, 0.0001/mindR);
	maxK = (1E5)*electron_charge*B0/(speed_of_light*minP);
	deltaLogK = (log(maxK) - log(minK))/(kgridNumber - 1);
	for(int i = 0; i < kgridNumber; ++i){
		logKgrid[i] = log(minK) + i*deltaLogK;
		kgrid[i] = exp(logKgrid[i]);
	}

	for(int i = 0; i < rgridNumber; ++i){
		magneticEnergy[i] = 0;
		maxRate[i] = 0;
		for(int k = 0; k < kgridNumber; ++k){
			magneticField[i][k] = (1E-7)*B0*B0*power(1/kgrid[k], 5.0/3.0)*power(kgrid[0],2.0/3.0);
			tempMagneticField[i][k] = magneticField[i][k];
			if(k == 0){
				largeScaleField[i][k] = sqrt(4*pi*magneticField[i][k]*kgrid[k]*deltaLogK + B0*B0);
			} else {
				largeScaleField[i][k] = sqrt(4*pi*magneticField[i][k]*kgrid[k]*deltaLogK + sqr(largeScaleField[i][k-1]));
			}
			growth_rate[i][k] = 0;
			alertNaNOrInfinity(magneticField[i][k], "magnetic field = NaN");
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		double magneticEnergy = 0;
		for(int k = 0; k < kgridNumber; ++k){
			magneticEnergy += magneticField[i][k]*kgrid[k]*deltaLogK;
		}
		magneticInductionSum[i] = sqrt(4*pi*magneticEnergy + B0*B0);
	}

	for(int i = 0; i < rgridNumber; ++i){
		switch(simulationType){
		//различные варианты профиля
		case 1 :
			middleDensity[i] = density0;
			if(i < rgridNumber/10){
				middleVelocity[i] = U0;
			} else {
				middleVelocity[i] = 0;
			}
			middlePressure[i] = pressure0;
			break;
		case 2 :
			middleDensity[i] = density0 + density0*0.01*sin(i*10*2*pi/rgridNumber);
			middleVelocity[i] = sqrt(_gamma*pressure0/density0)*density0*0.01*sin(i*10*2*pi/rgridNumber)/density0;
			middlePressure[i] = pressure0 + (_gamma*pressure0/density0)*density0*0.01*sin(i*10*2*pi/rgridNumber);
			break;
		case 3 :
			middleDensity[i] = density0;
			if(i <= rgridNumber/10){
				middleVelocity[i] = 10*i*U0/rgridNumber;
			} else if(i > rgridNumber/10 && i < 2*rgridNumber/10){
				middleVelocity[i] = U0;
			} else {
				middleVelocity[i] = 0;
			}
			middlePressure[i] = pressure0;
			break;
		case 4 :
			middleDensity[i] = density0;
			middleVelocity[i] = 0;
			{
				int count = 20;
				if(i < count){
					middlePressure[i] = (_gamma- 1)*initialEnergy/(4*pi*cube(grid[count])/3);
				} else {
					middlePressure[i] = pressure0;
				}
				shockWavePoint = count;
				shockWaveMoved = true;
			}
			break;
		case 5 :
			{
				double sigma = 4;
				double pressure = density0*U0*U0/sigma;
				int count = rgridNumber/2 - 1;
				if(i < count){
					middleDensity[i] = density0/sigma;//*sqr(middleGrid[i]/middleGrid[count-1]);
					middleVelocity[i] = U0;
					//middleVelocity[i] = 1;
					middlePressure[i] = pressure*0.0000000000001;
				} else {
					middleDensity[i] = density0;//*sqr(middleGrid[i]/middleGrid[count]);
					middleVelocity[i] = U0/sigma;
					//middleVelocity[i] = 0.25;
					middlePressure[i] = pressure*0.75;
				}
				shockWavePoint = count;
				shockWaveMoved = true;
			}
			break;
		default:
			middleDensity[i] = density0;
			if((i > 3*rgridNumber/10) && (i < 5*rgridNumber/10)){
				middleVelocity[i] = U0;
			} else {
				middleVelocity[i] = 0;
			}
			middlePressure[i] = pressure0;
			break;
		}
	}
	double pstep = exp(log(maxP/minP)/pgridNumber);
	logPgrid[0] = log(minP);
	//logPgrid[pgridNumber - 1] = log(maxP);
	deltaLogP = (log(maxP) - logPgrid[0])/(pgridNumber);
	pgrid[0] = minP;
	for(int i = 1; i < pgridNumber; ++i){
		logPgrid[i] = logPgrid[i-1] + deltaLogP;
		pgrid[i] = exp(logPgrid[i]);
		printf("%lf\n", pgrid[i]);
	}
	pgrid[pgridNumber-1] = maxP;

	/*for(int i = 0; i < rgridNumber/2; ++i){
		middleVelocity[i] = speed_of_light/10;
	}*/

	for(int i = 0; i < rgridNumber; ++i){
		if(i == 0){
			middleDeltaR[i] = middleGrid[i];
		} else {
			middleDeltaR[i] = middleGrid[i] - middleGrid[i-1];
		}
		double p;
		for(int j = 0; j < pgridNumber; ++j){
			p = exp(logPgrid[j]);
			distributionFunction[i][j] = 0;
			//double x = -(sqrt(sqr(massProton*speed_of_light*speed_of_light) + sqr(p*speed_of_light))-massProton*speed_of_light*speed_of_light)/(kBoltzman*temperatureIn(i));
			//distributionFunction[i][j] = exp(x);
			//distributionFunction[i][j] = epsilon;
		}
		pointDensity[i] = middleDensity[i];
		pointVelocity[i] = middleVelocity[i];
		pointEnthalpy[i] = (energy(i) + middlePressure[i])/middleDensity[i];
		tempDensity[i] = middleDensity[i];
		tempMomentum[i] = momentum(i);
		tempEnergy[i] = energy(i);
	}

	pointDensity[rgridNumber] = pointDensity[rgridNumber - 1];
	pointVelocity[rgridNumber] = pointVelocity[rgridNumber - 1];
	pointEnthalpy[rgridNumber] = (energy(rgridNumber-1) + middlePressure[rgridNumber-1])/middleDensity[rgridNumber-1];
	//grid[rgridNumber] = upstreamR;

	shockWaveT = 0;
	shockWaveSpeed = 0;
}

void Simulation::initializeArrays(){
	pgrid = new double[pgridNumber];
	logPgrid = new double[pgridNumber];
	grid = new double[rgridNumber + 1];
	gridsquare = new double[rgridNumber+ 1];
	middleGrid = new double[rgridNumber];
	deltaR = new double[rgridNumber];
	middleDeltaR = new double[rgridNumber];
	tempGrid = new double[rgridNumber];

	pointDensity = new double[rgridNumber + 1];
	pointVelocity = new double[rgridNumber + 1];
	pointEnthalpy = new double[rgridNumber + 1];
	pointSoundSpeed = new double[rgridNumber + 1];
	middleDensity = new double[rgridNumber];
	middleVelocity = new double[rgridNumber];
	middlePressure = new double[rgridNumber];

	tempDensity = new double[rgridNumber];
	tempMomentum = new double[rgridNumber];
	tempEnergy = new double[rgridNumber];

	vscattering = new double[rgridNumber];

	dFlux = new double[rgridNumber];
	dFluxPlus = new double*[rgridNumber];
	dFluxMinus = new double*[rgridNumber];

	mFlux = new double[rgridNumber];
	mFluxPlus = new double*[rgridNumber];
	mFluxMinus = new double*[rgridNumber];

	eFlux = new double[rgridNumber];
	eFluxPlus = new double*[rgridNumber];
	eFluxMinus = new double*[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		vscattering[i] = 0;
		dFluxPlus[i] = new double[3];
		dFluxMinus[i] = new double[3];
		mFluxPlus[i] = new double[3];
		mFluxMinus[i] = new double[3];
		eFluxPlus[i] = new double[3];
		eFluxMinus[i] = new double[3];
	}

	cosmicRayPressure = new double[rgridNumber+1];
	cosmicRayConcentration = new double[rgridNumber+1];
	tempU = new double[rgridNumber];
	integratedFlux = new double[rgridNumber];

	distributionFunction = new double*[rgridNumber+1];
	tempDistributionFunction = new double*[rgridNumber+1];

	for(int i = 0; i <= rgridNumber; i ++){
		cosmicRayPressure[i] = 0;
		distributionFunction[i] = new double[pgridNumber];
		tempDistributionFunction[i] = new double[pgridNumber];
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = 0;
			tempDistributionFunction[i][j] = 0;
		}
	}

	diffusionCoef = new double*[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		diffusionCoef[i] = new double[pgridNumber];
		integratedFlux[i] = 0;
	}

	kgrid = new double[kgridNumber];
	logKgrid = new double[kgridNumber];

	magneticField = new double*[rgridNumber];
	tempMagneticField = new double*[rgridNumber];
	largeScaleField = new double*[rgridNumber];
	growth_rate = new double*[rgridNumber];
	magneticInductionSum = new double[rgridNumber];
	maxRate = new double[rgridNumber];
	magneticEnergy = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		magneticEnergy[i] = 0;
		maxRate[i] = 0;
		magneticField[i] = new double[kgridNumber];
		tempMagneticField[i] = new double[kgridNumber];
		largeScaleField[i] = new double[kgridNumber];
		growth_rate[i] = new double[kgridNumber];
		for(int k = 0; k < kgridNumber; ++k){
			growth_rate[i][k] = 0;
		}
	}

	crflux = new double*[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		crflux[i] = new double[pgridNumber];
		for(int j = 0;j < pgridNumber; ++j){
			crflux[i][j] = 0;
		}
	}
}

void Simulation::updateAfterSerialization(){
	for(int i = 1; i <= rgridNumber; ++i){
		if(grid[i] <= grid[i-1]){
			printf("grid[i] <= grid[i-1]\n");
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		middleGrid[i] = (grid[i] + grid[i+1])/2;
		tempGrid[i] = grid[i];
		deltaR[i] = grid[i+1] - grid[i];
		gridsquare[i] = sqr(grid[i]);
	}
	tempGrid[rgridNumber] = grid[rgridNumber];

	double pstep = exp(log(maxP/minP)/pgridNumber);
	logPgrid[0] = log(minP);
	deltaLogP = (log(maxP) - logPgrid[0])/(pgridNumber);
	pgrid[0] = minP;
	for(int i = 1; i < pgridNumber; ++i){
		logPgrid[i] = logPgrid[i-1] + deltaLogP;
		pgrid[i] = exp(logPgrid[i]);
	}
	pgrid[pgridNumber-1] = maxP;
	/*for(int i = 0; i < pgridNumber; ++i){
		logPgrid[i] = log(pgrid[i]);
	}

	deltaLogP = logPgrid[1] - logPgrid[0];*/

	deltaLogK = (log(maxK) - log(minK))/(kgridNumber - 1);
	for(int i = 0; i < kgridNumber; ++i){
		logKgrid[i] = log(minK) + i*deltaLogK;
		kgrid[i] = exp(logKgrid[i]);
	}
	/*for(int i = 0; i < kgridNumber; ++i){
		logKgrid[i] = log(kgrid[i]);
	}

	deltaLogK = logKgrid[1] - logKgrid[0];*/

	for(int i = 0; i < rgridNumber; ++i){
		if(i == 0){
			middleDeltaR[i] = middleGrid[i];
		} else {
			middleDeltaR[i] = middleGrid[i] - middleGrid[i-1];
		}

		pointDensity[i] = middleDensity[i];
		pointVelocity[i] = middleVelocity[i];
		pointEnthalpy[i] = (energy(i) + middlePressure[i])/middleDensity[i];
		tempDensity[i] = middleDensity[i];
		tempMomentum[i] = momentum(i);
		tempEnergy[i] = energy(i);
	}

	pointDensity[rgridNumber] = pointDensity[rgridNumber - 1];
	pointVelocity[rgridNumber] = pointVelocity[rgridNumber - 1];
	pointEnthalpy[rgridNumber] = (energy(rgridNumber-1) + middlePressure[rgridNumber-1])/middleDensity[rgridNumber-1];

	for(int i = 0; i < rgridNumber; ++i){
		magneticEnergy[i] = 0;
		maxRate[i] = 0;
		for(int k = 0; k < kgridNumber; ++k){
			tempMagneticField[i][k] = magneticField[i][k];
			if(k == 0){
				largeScaleField[i][k] = sqrt(4*pi*magneticField[i][k]*kgrid[k]*deltaLogK + B0*B0);
			} else {
				largeScaleField[i][k] = sqrt(4*pi*magneticField[i][k]*kgrid[k]*deltaLogK + sqr(largeScaleField[i][k-1]));
			}
			growth_rate[i][k] = 0;
			alertNaNOrInfinity(magneticField[i][k], "magnetic field = NaN");
		}
	}
}