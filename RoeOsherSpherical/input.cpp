#include <stdio.h>
#include "input.h"
#include "simulation.h"

Simulation* readInput(){
	FILE* infile = fopen("tamc_in.dat","r");
	if (infile != NULL){
		char* s = NULL;
		Simulation* simulation = new Simulation();
		fscanf(infile,"%d",&simulation->iterationNumber);
		char ch = ' ';
		while (ch != '\n') {
			fscanf(infile, "%c", &ch);
		}
        fscanf(infile,"%d",&simulation->particlesNumber);
		ch = ' ';
		while (ch != '\n') {
			fscanf(infile, "%c", &ch);
		}
		fscanf(infile,"%lf",&simulation->U0);
		ch = ' ';
		while (ch != '\n') {
			fscanf(infile, "%c", &ch);
		}
		fscanf(infile,"%lf",&simulation->density0);
		ch = ' ';
		while (ch != '\n') {
			fscanf(infile, "%c", &ch);
		}
        fscanf(infile,"%lf",&simulation->B0);
		ch = ' ';
		while (ch != '\n') {
			fscanf(infile, "%c", &ch);
		}
        fscanf(infile,"%lf",&simulation->temperature);
		ch = ' ';
		while (ch != '\n') {
			fscanf(infile, "%c", &ch);
		}
		fscanf(infile,"%lf",&simulation->upstreamR);
		ch = ' ';
		while (ch != '\n') {
			fscanf(infile, "%c", &ch);
		}
		fscanf(infile,"%d",&simulation->rgridNumber);
		ch = ' ';
		while (ch != '\n') {
			fscanf(infile, "%c", &ch);
		}
		fscanf(infile,"%d",&simulation->simulationType);
		ch = ' ';
		while (ch != '\n') {
			fscanf(infile, "%c", &ch);
		}
		fscanf(infile,"%lf",&simulation->maxTime);
		fclose(infile);
		return simulation;
	} else {
		return NULL;
	}
}

void deserialize(FILE* hydroFile, FILE* distributionFile, FILE* fieldFile, FILE* gridFile, FILE* pgridFile, FILE* kgridFile, FILE* infoFile, Simulation* simulation){
	fscanf(infoFile,"%d", &simulation->iterationNumber);
    fscanf(infoFile,"%d", &simulation->particlesNumber);
	fscanf(infoFile,"%lf", &simulation->U0);
	fscanf(infoFile,"%lf", &simulation->density0);
    fscanf(infoFile,"%lf", &simulation->B0);
    fscanf(infoFile,"%lf", &simulation->temperature);
	fscanf(infoFile,"%lf", &simulation->upstreamR);
	fscanf(infoFile,"%d", &simulation->rgridNumber);
	fscanf(infoFile,"%d", &simulation->simulationType);
	fscanf(infoFile,"%lf", &simulation->maxTime);
	fscanf(infoFile, "%lf", &simulation->myTime);
	fscanf(infoFile, "%d", &simulation->currentIteration);
	fscanf(infoFile, "%lf", &simulation->injectedParticles);
	fscanf(infoFile, "%d", &simulation->shockWavePoint);

	int a;
	fscanf(infoFile, "%d", &a);
	simulation->shockWaveMoved = (a==1);

	fscanf(infoFile, "%d", &simulation->prevShockWavePoint);
	fscanf(infoFile, "%lf", &simulation->shockWaveSpeed);
	fscanf(infoFile, "%lf", &simulation->shockWaveT);

	fscanf(infoFile, "%lf", &simulation->minP);
	fscanf(infoFile, "%lf", &simulation->maxP);
	fscanf(infoFile, "%lf", &simulation->minK);
	fscanf(infoFile, "%lf", &simulation->maxK);

	simulation->initializeArrays();

	/*for(int j = 0; j < pgridNumber; ++j){
		double p;
		fprintf(pgridFile, "%60.50lf", &p);
		printf("%lf\n",p);
	}

	for(int j = 0; j < kgridNumber; ++j){
		int k;
		fprintf(kgridFile, "%60.50lf",&simulation->kgrid[j]);
	}*/

	for(int i = 0; i <= simulation->rgridNumber; ++i){
		fscanf(gridFile,"%lf", &simulation->grid[i]);
		for(int j = 0; j < pgridNumber; ++j){
			fscanf(distributionFile, "%lf", &simulation->distributionFunction[i][j]);
		}
	}

	for(int i = 0; i < simulation->rgridNumber; ++i){
		fscanf(hydroFile, "%lf %lf %lf", &simulation->middleDensity[i], &simulation->middleVelocity[i], &simulation->middlePressure[i]);
		for(int k = 0; k < kgridNumber; ++k){
			fscanf(fieldFile, "%lf", &simulation->magneticField[i][k]);
		}
	}

	simulation->updateAfterSerialization();

	simulation->serialized= true;
}
	