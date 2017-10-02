// test.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <mpi.h>

#include "specialmath.h"



void testGMRES(){
	const int generalNumber = 8;
	double** generalMatrix = new double*[generalNumber];
	double* generalRightPart = new double[generalNumber];

	///only 3-diagonal generalMatrix!!!

	for(int i = 0; i < generalNumber; ++i){
		generalMatrix[i] = new double[generalNumber];
		for(int j = 0; j < generalNumber; ++j){
			if(i == j){
				generalMatrix[i][j] = 1;
			} else {
				generalMatrix[i][j] = 0;
			}
		}
	}

    generalMatrix[1][2] = 0.2;
    generalMatrix[2][3] = 0.1;
    generalMatrix[3][1] = 0.3;
    generalMatrix[3][4] = 0.1;
    generalMatrix[4][5] = -0.5;
    generalMatrix[5][6] = 0.1;
    generalMatrix[6][5] = 0.1;

	int rank;
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int localNumber = ((generalNumber-2)/nprocs) + 2;


	for(int i = 0; i < generalNumber; ++i){
		generalRightPart[i] = 0;
	}
    generalRightPart[1] = 1;
    generalRightPart[2] = 1;
    generalRightPart[5] = 1;
    generalRightPart[6] = 1;
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0) {
		for (int i = 0; i < generalNumber; ++i) {
			for (int j = 0; j < generalNumber; ++j) {
				printf("%lf     ", generalMatrix[i][j]);
			}
			printf("   %lf\n", generalRightPart[i]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	double** localMatrix = new double*[localNumber];
	double* localRightPart = new double[localNumber];

	for(int i = 0; i < localNumber; ++i){
		localMatrix[i] = new double[localNumber];
		int globalI = i + (localNumber - 2)*rank;
		for(int j = 0; j < localNumber; ++j){
			int globalJ = j + (localNumber - 2)*rank;
			localMatrix[i][j] = generalMatrix[globalI][globalJ];
		}
		localRightPart[i] = generalRightPart[globalI];
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for(int j = 0; j < nprocs; ++j) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (j == rank) {
			printf("rank = %d\n", rank);
			for (int i = 0; i < localNumber; ++i) {
				for(int k = 0; k < localNumber; ++k){
					printf("%g   ", localMatrix[i][k]);
				}
				printf("%g\n", localRightPart[i]);
			}
		}
	}

	double* result = generalizedMinimalResidualMethod(localMatrix, localRightPart, localNumber);
	MPI_Barrier(MPI_COMM_WORLD);

	for(int j = 0; j < nprocs; ++j) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (j == rank) {
			printf("rank = %d\n", rank);
			for (int i = 0; i < localNumber; ++i) {
				printf("%15.10lf\n", result[i]);
			}
		}
	}
	/*if(rank == 0){
		printf("rank = %d\n", rank);
		for (int i = 0; i < localNumber; ++i) {
			printf("%15.10lf\n", result[i]);
		}
	}*/

	MPI_Barrier(MPI_COMM_WORLD);

	for(int j = 0; j < nprocs; ++j) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (j == rank) {
			printf("rank = %d\n", rank);
			printf("error\n");
			for (int i = 0; i < localNumber; ++i) {
				double value = -localRightPart[i];
				for (int j = 0; j < localNumber; ++j) {
					value += localMatrix[i][j] * result[j];
				}
				printf("%15.10lf\n", value);
			}
		}
	}

	printf("end test GMRES\n");
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	testGMRES();

	MPI_Finalize();

	return 0;
}

