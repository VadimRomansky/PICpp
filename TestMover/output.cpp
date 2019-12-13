#include "stdio.h"
#include <stdlib.h>
#include "math.h"

#include "constants.h"

#include "output.h"

void outputDistribution(const char* outFileName, double** momentum, int chch){
	double minMomentum = 0; //todo something else
	double maxMomentum = 0;
	for (int i = 0; i < chch; ++i) {
		double p = sqrt(momentum[i][0]*momentum[i][0] + momentum[i][1]*momentum[i][1] + momentum[i][2]*momentum[i][2]);
		if (minMomentum <= 0) {
			minMomentum = p;
		}
		if ((p  < minMomentum) && (p > 0)) {
			minMomentum = p;
		} else {
			if (p > maxMomentum) {
				maxMomentum = p;
			}
		}
	}
	maxMomentum = maxMomentum*1.1;

	if(maxMomentum <= 0) {
		maxMomentum = massElectron*c;
		minMomentum = maxMomentum/2;
	}


	if (maxMomentum - minMomentum < 1E-8 * maxMomentum) {
		maxMomentum = maxMomentum * 2;
		minMomentum = minMomentum / 2;
	}


	double pgrid[pnumber + 1];
	double distribution[pnumber];
	double tempDistribution[pnumber];
	double logMinMomentum = log(minMomentum);
	pgrid[0] = minMomentum;
	distribution[0] = 0;
	double deltaLogP = (log(maxMomentum) - log(minMomentum)) / (pnumber);
	for (int i = 1; i < pnumber; ++i) {
		distribution[i] = 0;
		pgrid[i] = exp(logMinMomentum + i * deltaLogP);
	}
	pgrid[pnumber] = maxMomentum;

	double weight;
	weight = 0;

	for (int i = 0; i < chch; ++i) {
		double p = sqrt(momentum[i][0]*momentum[i][0] + momentum[i][1]*momentum[i][1] + momentum[i][2]*momentum[i][2]);
		int j = (log(p) - logMinMomentum) / deltaLogP;
		if (j >= 0 && j < pnumber) {
			distribution[j] += 1.0;
			weight += 1.0;
		}
	}


	for (int i = 0; i < pnumber; ++i) {
		if(weight > 0){
		distribution[i] /= (weight * (pgrid[i + 1] - pgrid[i]));
		} else {
			distribution[i] = 0;
		}
		//todo check
		if (maxMomentum <= 0 || minMomentum <= 0) {
			distribution[i] = 0;
			//pgrid[i] = 0;
		}
	}

	FILE* outFile = NULL;
	outFile = fopen(outFileName, "w");
	for (int i = 0; i < pnumber; ++i) {
		fprintf(outFile, "%22.15g %22.15g\n", (pgrid[i] + pgrid[i + 1]) / 2, distribution[i]);
	}
	fclose(outFile);
}

void writeMaxParticles(const char* fileName, double** momentum, int N){
	FILE* outFile = fopen(fileName, "w");
	const int Nwrite = 10;
	int numbers[Nwrite];
	for(int i = 0; i < Nwrite; ++i){
		numbers[i] = i;
	}
	for(int i = 0; i < Nwrite - 1; ++i){
		for(int j = 0; j < Nwrite - i - 1; ++j){
			double p1 = momentum[numbers[Nwrite - i - 1]][0]*momentum[numbers[Nwrite - i - 1]][0]+
				momentum[numbers[Nwrite - i - 1]][1]*momentum[numbers[Nwrite - i - 1]][1]+
				momentum[numbers[Nwrite - i - 1]][2]*momentum[numbers[Nwrite - i - 1]][2];
			double p2 = momentum[numbers[j]][0]*momentum[numbers[j]][0]+
				momentum[numbers[j]][1]*momentum[numbers[j]][1]+
				momentum[numbers[j]][2]*momentum[numbers[j]][2];
			if(p2 > p1){
				int temp = numbers[Nwrite - i - 1];
				numbers[Nwrite - i - 1] = numbers[j];
				numbers[j] = temp;
			}
		}
	}

	for(int k = Nwrite; k < N; ++k){
		double p1 = momentum[k][0]*momentum[k][0]+
				momentum[k][1]*momentum[k][1]+
				momentum[k][2]*momentum[k][2];
		double p2 = momentum[numbers[0]][0]*momentum[numbers[0]][0]+
				momentum[numbers[0]][1]*momentum[numbers[0]][1]+
				momentum[numbers[0]][2]*momentum[numbers[0]][2];
		if(p1 > p2){
			numbers[0] = k;
			for(int i = 0; i < Nwrite - 1; ++i){
				p1 = momentum[numbers[i]][0]*momentum[numbers[i]][0]+
					momentum[numbers[i]][1]*momentum[numbers[i]][1]+
					momentum[numbers[i]][2]*momentum[numbers[i]][2];
				p2 = momentum[numbers[i+1]][0]*momentum[numbers[i+1]][0]+
					momentum[numbers[i+1]][1]*momentum[numbers[i+1]][1]+
					momentum[numbers[i+1]][2]*momentum[numbers[i+1]][2];
				if(p1 > p2){
					int temp = numbers[i];
					numbers[i] = numbers[i+1];
					numbers[i+1] = temp;
				}
			}
		}
	}

	for(int i = 0; i < Nwrite; ++i){
		double p1 = momentum[numbers[i]][0]*momentum[numbers[i]][0]+
				momentum[numbers[i]][1]*momentum[numbers[i]][1]+
				momentum[numbers[i]][2]*momentum[numbers[i]][2];
		fprintf(outFile, "%d %g\n",numbers[i], sqrt(p1));
	}
	fclose(outFile);
}

void outputField(const char* fileNameX, const char* fileNameY, const char* fileNameZ, double*** downstreamField, double*** middleField, double*** upstreamField, int downstreamNx, int middleNx, int upstreamNx, int Ny){
	FILE* outFileX = fopen(fileNameX, "w");
	FILE* outFileY = fopen(fileNameY, "w");
	FILE* outFileZ = fopen(fileNameZ, "w");
	int maxK = 3;
	for(int k = 0; k < maxK; ++k){
		for(int i = 0; i < downstreamNx; ++i){
			for(int j = 0; j < Ny; ++j){
				fprintf(outFileX,"%g ", downstreamField[i][j][0]);
				fprintf(outFileY,"%g ", downstreamField[i][j][1]);
				fprintf(outFileZ,"%g ", downstreamField[i][j][2]);
			}
			fprintf(outFileX,"\n");
			fprintf(outFileY,"\n");
			fprintf(outFileZ,"\n");
		}
	}
	for(int i = 0; i < middleNx; ++i){
		for(int j = 0; j < Ny; ++j){
			fprintf(outFileX,"%g ", middleField[i][j][0]);
			fprintf(outFileY,"%g ", middleField[i][j][1]);
			fprintf(outFileZ,"%g ", middleField[i][j][2]);
		}
		fprintf(outFileX,"\n");
		fprintf(outFileY,"\n");
		fprintf(outFileZ,"\n");
	}
	for(int k = 0; k < maxK; ++k){
		for(int i = 0; i < upstreamNx; ++i){
			for(int j = 0; j < Ny; ++j){
				fprintf(outFileX,"%g ", upstreamField[i][j][0]);
				fprintf(outFileY,"%g ", upstreamField[i][j][1]);
				fprintf(outFileZ,"%g ", upstreamField[i][j][2]);
			}
			fprintf(outFileX,"\n");
			fprintf(outFileY,"\n");
			fprintf(outFileZ,"\n");
		}
	}
	fclose(outFileX);
	fclose(outFileY);
	fclose(outFileZ);
}