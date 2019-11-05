#include <time.h>
#include "math.h"
#include "stdio.h"
#include <stdlib.h>
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

//изменение сетки

void Simulation::updateGrid(){
	if ((shockWavePoint < 1) || (shockWavePoint > rgridNumber - 1)) return;
	if( !shockWaveMoved) {
		return;
	}
	printf("updating grid\n");
	double shockWaveR = grid[shockWavePoint];
	double rightR = grid[rgridNumber] - shockWaveR;

	int indent = 3; 

	int leftPoints = shockWavePoint - indent - 2;
	int rightPonts = rgridNumber - indent - shockWavePoint + 1;

	double R1 = grid[shockWavePoint - indent - 1];
	double R2 = grid[rgridNumber] - grid[shockWavePoint+3];
	double a = 100;
	double b = 1000;
	double h1 = leftPoints/log(1.0+a);
	double h2 = rightPonts/log(1.0+b);

	tempGrid[0] = 0;
	for(int i = 1; i < shockWavePoint - indent - 1; ++i){
		tempGrid[i] = (R1/a)*(1 - exp(-(i-leftPoints)/h1)) + R1;
	}

	for(int i = shockWavePoint - indent-1; i < shockWavePoint + indent-1; ++i){
		tempGrid[i] = grid[i+1];
	}

	for(int i = shockWavePoint + indent-1; i < rgridNumber; ++i){
		tempGrid[i] = (R2/b)*(exp((i+1-shockWavePoint-indent)/h2)-1.0) + grid[shockWavePoint+3];
	}
	tempGrid[rgridNumber] = grid[rgridNumber];

	shockWavePoint -= 1;

	for(int i = 1; i <= rgridNumber; ++i){
		if(tempGrid[i] < tempGrid[i-1]){
			printf("grid[i] < grid[i-1]\n");
		}
	}

	redistributeValues();
}


//перераспределение величин между ячейками новой сетки

void Simulation::redistributeValues(){
	int oldCount = 1;
	double tDensity = 0;
	double tMomentum = 0;
	double tEnergy = 0;

	double** newDistributionFunction = new double*[rgridNumber];
	double* tempFunction = new double[pgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		newDistributionFunction[i] = new double[pgridNumber];
		double newVolume = cube(tempGrid[i+1]) - cube(tempGrid[i]);
		bool oldChanged = false;
		if(tempGrid[i+1] > grid[oldCount]){
			tDensity = middleDensity[oldCount - 1]*(cube(grid[oldCount]) - cube(tempGrid[i]))/(newVolume);
			tMomentum = momentum(oldCount - 1)*(cube(grid[oldCount]) - cube(tempGrid[i]))/(newVolume);
			tEnergy = energy(oldCount - 1)*(cube(grid[oldCount]) - cube(tempGrid[i]))/(newVolume);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j]*(cube(grid[oldCount]) - cube(tempGrid[i]))/(newVolume);
			}
			++oldCount;
			oldChanged = true;
		} else {
			tDensity = middleDensity[oldCount - 1];
			tMomentum = momentum(oldCount - 1);
			tEnergy = energy(oldCount - 1);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j];
			}
		}
		while(grid[oldCount] < tempGrid[i+1]){
			tDensity += middleDensity[oldCount - 1]*volume(oldCount - 1)/(4*pi*newVolume/3);
			tMomentum += momentum(oldCount - 1)*volume(oldCount - 1)/(4*pi*newVolume/3);
			tEnergy += energy(oldCount - 1)*volume(oldCount - 1)/(4*pi*newVolume/3);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j]*volume(oldCount - 1)/(4*pi*newVolume/3);
			}
			++oldCount;
			oldChanged = true;
			if(oldCount > rgridNumber + 1){
				printf("oldCount > rgridNUmber + 1\n");
			}
		}
		if(oldChanged){
			tDensity += middleDensity[oldCount - 1]*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(newVolume);
			tMomentum += momentum(oldCount - 1)*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(newVolume);
			tEnergy += energy(oldCount - 1)*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(newVolume);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j]*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(newVolume);
			}
		}
		tempDensity[i] = tDensity;
		alertNaNOrInfinity(tempDensity[i], "newDensity = NaN");
		alertNegative(tempDensity[i], "newDensity < 0");
		tempMomentum[i] = tMomentum;
		alertNaNOrInfinity(tempMomentum[i], "newMomentum = NaN");
		tempEnergy[i] = tEnergy;
		alertNaNOrInfinity(tempEnergy[i], "newEnergy = NaN");
		alertNegative(tempEnergy[i], "newEnergy < 0");
		for(int j = 0; j < pgridNumber; ++j){
			newDistributionFunction[i][j] = tempFunction[j];
			alertNaNOrInfinity(newDistributionFunction[i][j], "newDistribution = NaN");
			alertNegative(newDistributionFunction[i][j], "newDistribution < 0");
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		grid[i + 1] = tempGrid[i + 1];
		middleGrid[i] = (grid[i] + grid[i+1])/2;
		deltaR[i] = grid[i+1] - grid[i];
		if(i == 0){
			middleDeltaR[i] = middleGrid[i];
		} else {
			middleDeltaR[i] = middleGrid[i] - middleGrid[i-1];
		}

		middleDensity[i] = tempDensity[i];
		if(tempDensity[i] <= epsilon*density0){
			middleVelocity[i] = 0;
		} else {
			middleVelocity[i] = tempMomentum[i]/tempDensity[i];
		}
		double tempPressure = (tempEnergy[i] - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(_gamma - 1);
		
		if(tempPressure < 0){
			middlePressure[i] = 0.01*min2(middlePressure[i+1],middlePressure[i]);
			printf("pressure < 0\n");
		} else {
			middlePressure[i] = tempPressure;
		}
		alertNegative(middlePressure[i], "middlePressure < 0");
		alertNaNOrInfinity(middlePressure[i], "middlePressure = NaN");

		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = newDistributionFunction[i][j];
		}
	}
	if(tempDensity[rgridNumber - 1] < middleDensity[rgridNumber - 1]){
		printf("aaa\n");
	}

	for(int i = 0; i < rgridNumber; ++i){
		delete[] newDistributionFunction[i];
	}
	delete[] newDistributionFunction;
	delete[] tempFunction;
}