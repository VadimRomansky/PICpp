#include "stdio.h"
#include <stdlib.h>
#include "output.h"
#include "constants.h"
#include "util.h"

void output(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		double t = simulation->temperatureIn(i);
		double v = simulation->middleVelocity[i]-simulation->shockWaveSpeed;
		double shockWaveR = 0;
		if(simulation->shockWavePoint > 0){
			shockWaveR = simulation->middleGrid[simulation->shockWavePoint];
		}
		double tau = abs2((simulation->middleGrid[i]-shockWaveR)*simulation->maxRate[i]/v);
		if(abs(v)< 1){
			tau = 0;
		}
		fprintf(outFile, "%g %g %g %g %g %g %g %g %g %g %g", simulation->grid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], simulation->cosmicRayPressure[i], t, simulation->magneticInductionSum[i]-simulation->B0, simulation->magneticEnergy[i], tau, simulation->cosmicRayConcentration[i], simulation->vscattering[i]);
		fprintf(outFile, "\n");
	}
}

void serialize(FILE* hydroFile, FILE* distributionFile, FILE* fieldFile, FILE* gridFile, FILE* pgridFile, FILE* kgridFile, FILE* infoFile, Simulation* simulation){
	fprintf(infoFile,"%d\n", simulation->iterationNumber);
    fprintf(infoFile,"%d\n", simulation->particlesNumber);
	fprintf(infoFile,"%60.50lf\n", simulation->U0);
	fprintf(infoFile,"%60.50lf\n", simulation->density0);
    fprintf(infoFile,"%60.50lf\n", simulation->B0);
    fprintf(infoFile,"%60.50lf\n", simulation->temperature);
	fprintf(infoFile,"%60.50lf\n", simulation->upstreamR);
	fprintf(infoFile,"%d\n", simulation->rgridNumber);
	fprintf(infoFile,"%d\n", simulation->simulationType);
	fprintf(infoFile,"%60.50lf\n", simulation->maxTime);
	fprintf(infoFile, "%60.50lf\n", simulation->myTime);
	fprintf(infoFile, "%d\n", simulation->currentIteration);
	fprintf(infoFile, "%60.50lf\n", simulation->injectedParticles);
	fprintf(infoFile, "%d\n", simulation->shockWavePoint);
	fprintf(infoFile, "%d\n", simulation->shockWaveMoved);
	fprintf(infoFile, "%d\n", simulation->prevShockWavePoint);
	fprintf(infoFile, "%60.50lf\n", simulation->shockWaveSpeed);
	fprintf(infoFile, "%60.50lf\n", simulation->shockWaveT);
	fprintf(infoFile, "%60.50lf\n", simulation->minP);
	fprintf(infoFile, "%60.50lf\n", simulation->maxP);
	fprintf(infoFile, "%60.50lf\n", simulation->minK);
	fprintf(infoFile, "%60.50lf\n", simulation->maxK);

	for(int j = 0; j < pgridNumber; ++j){
		fprintf(pgridFile, "%60.50lf\n",simulation->pgrid[j]);
	}

	for(int j = 0; j < kgridNumber; ++j){
		fprintf(kgridFile, "%60.50lf\n",simulation->kgrid[j]);
	}

	for(int i = 0; i <= simulation->rgridNumber; ++i){
		fprintf(gridFile,"%60.50lf\n", simulation->grid[i]);
		for(int j = 0; j < pgridNumber; ++j){
			fprintf(distributionFile, "%60.50lf ", simulation->distributionFunction[i][j]);
		}
		fprintf(distributionFile, "\n");
	}

	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(hydroFile, "%60.50lf %60.50lf %60.50lf\n", simulation->middleDensity[i], simulation->middleVelocity[i], simulation->middlePressure[i]);
		for(int k = 0; k < kgridNumber; ++k){
			fprintf(fieldFile, "%60.50lf ", simulation->magneticField[i][k]);
		}
		fprintf(fieldFile, "\n");
	}
}

void outputDistribution(FILE* distributionFile, FILE* fullDistributionFile, FILE* coordinateDistributionFile, Simulation* simulation){
	double* fullDistribution = new double[pgridNumber];
	double volume = simulation->upstreamR;
	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] = 0;
	}
	for(int i = 0; i < simulation->rgridNumber; ++i){
		for(int j = 0; j < pgridNumber - 1; ++j){
			double p = simulation->pgrid[j];
			fullDistribution[j] += simulation->volume(i)*simulation->distributionFunction[i][j];
		}
		fprintf(coordinateDistributionFile, "%20.10lf %g %g %g %g %g\n", simulation->grid[i], simulation->distributionFunction[i][injectionMomentum], simulation->distributionFunction[i][pgridNumber-1], simulation->crflux[i][injectionMomentum], simulation->crflux[i][pgridNumber-1], simulation->integratedFlux[i]);
	}

	if(simulation->shockWavePoint > 0 && simulation->shockWavePoint < simulation->rgridNumber){
		for(int j = 0; j < pgridNumber; ++j){
			fprintf(distributionFile, "%g %g %g\n", simulation->pgrid[j], simulation->distributionFunction[simulation->shockWavePoint][j], simulation->crflux[simulation->shockWavePoint][j]);
		}
	}

	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] /= volume;
		double p = simulation->pgrid[j];
		fprintf(fullDistributionFile, "%g %g\n", p, fullDistribution[j]);
	}
	delete[] fullDistribution;
}

void outputDistributionP3(FILE* distributionFile, FILE* fullDistributionFile, FILE* coordinateDistributionFile, Simulation* simulation){
	double* fullDistribution = new double[pgridNumber];
	double volume = simulation->upstreamR;
	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] = 0;
	}
	for(int i = 0; i < simulation->rgridNumber; ++i){
		for(int j = 0; j < pgridNumber - 1; ++j){
			double p = simulation->pgrid[j];
			fullDistribution[j] += simulation->volume(i)*simulation->distributionFunction[i][j];
		}
		fprintf(coordinateDistributionFile, "%20.10lf %g %g %g %g %g\n", simulation->grid[i], simulation->distributionFunction[i][injectionMomentum]/(cube(simulation->pgrid[injectionMomentum])), simulation->distributionFunction[i][pgridNumber-1]/(cube(simulation->pgrid[pgridNumber-1])), simulation->crflux[i][injectionMomentum], simulation->crflux[i][pgridNumber-1], simulation->integratedFlux[i]);
	}

	if(simulation->shockWavePoint > 0 && simulation->shockWavePoint < simulation->rgridNumber){
		for(int j = 0; j < pgridNumber; ++j){
			fprintf(distributionFile, "%g %g %g %g\n", simulation->pgrid[j], simulation->distributionFunction[simulation->shockWavePoint][j]/cube(simulation->pgrid[j]), simulation->crflux[simulation->shockWavePoint][j]/(cube(simulation->pgrid[j])*simulation->deltaLogP), simulation->distributionFunction[simulation->rgridNumber-10][j]/cube(simulation->pgrid[j]));
		}
	}

	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] /= volume;
		double p = simulation->pgrid[j];
		fprintf(fullDistributionFile, "%g %g\n", p, fullDistribution[j]/(cube(p)));
	}
	delete[] fullDistribution;
}


void outputNewGrid(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(outFile, "%d %17.12lf\n", i, simulation->tempGrid[i]);
	}
}

void outMatrix(double* a, double* c, double* b, int N, double* f, double* x){
	FILE* file = fopen("output/matrix.dat","w");
	fprintf(file, "%g %g %g %g %g\n" , 0.0, c[0], 0.0, f[0], x[0]);
	for(int i = 1; i <= N; ++i){
		fprintf(file, "%g %g %g %g %g\n" , a[i-1], c[i], b[i-1], f[i], x[i]);
	}
	fclose(file);
}

void outputField(FILE* outFile, FILE* coordinateFile, FILE* outFull, FILE* coefFile, FILE* xfile, FILE* kfile, FILE* pfile,  Simulation* simulation){
	double* integralField = new double[kgridNumber];
	for(int k = 0; k < kgridNumber; ++k){
		integralField[k] = 0;
	}

	double volume = 0; 
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(xfile, "%g\n", simulation->grid[i]);
		fprintf(coordinateFile, "%g %g\n", simulation->grid[i],simulation->magneticField[i][35]);
        if(abs2(i-simulation->shockWavePoint) < 20){
			volume += simulation->volume(i);
		}
		for(int k = 0; k < kgridNumber; ++k){
			fprintf(outFull, "%g ", simulation->magneticField[i][k]);
            if(abs2(i-simulation->shockWavePoint) < 20){
				integralField[k] += simulation->magneticField[i][k]*simulation->volume(i);
			}
		}
		fprintf(outFull, "\n");
	}

	for(int j = 0; j < pgridNumber; ++j){
		fprintf(pfile, "%g\n", simulation->pgrid[j]);
	}

	for(int k = 0; k < kgridNumber; ++k){
		integralField[k] /= volume;
	}

	if(simulation->shockWavePoint > 0){
		int shockWavePoint = simulation->shockWavePoint;
		for(int k = 0; k < kgridNumber; ++k){
			fprintf(kfile, "%g\n", simulation->kgrid[k]);
			fprintf(outFile, "%g %g %g %g\n", simulation->kgrid[k], simulation->magneticField[shockWavePoint][k], simulation->growth_rate[shockWavePoint][k], integralField[k]);
		}

		for(int j = 0; j < pgridNumber; ++j){
			fprintf(coefFile, "%g %g %g %g\n", simulation->pgrid[j], simulation->diffusionCoef[shockWavePoint-1][j], simulation->diffusionCoef[shockWavePoint][j], simulation->diffusionCoef[shockWavePoint+1][j]);
		}
	}


	delete[] integralField;
}
