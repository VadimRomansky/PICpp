// test.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <mpi.h>

#include "specialmath.h"
#include "complex.h"

void testFastFourier() {
	const int number = 2*2*2*2*2;
	double* a = new double[number];
	for(int i = 0; i < number; ++i) {
		a[i] = cos(2*pi*i/number) + cos(16*pi*i/number);
		//a[i] = rand()%number;
		//a[i] = i;
	}
	//sortInputFastFourier(a, number);
	time_t time0;
	time_t time1;
	time_t time2;

	printf("evaluating fourier\n");
	time(&time0);
	Complex* result =evaluateFourierTranslation(a, number);
	time(&time1);
	time1 = time1 - time0;
	printf("evaluating fast fourier\n");
	time(&time0);
	Complex* resultFast = fastFourierTransition(a, number);
	time(&time2);
	time2 = time2 - time0;

	for(int i = 0; i < number; ++i) {
		printf("%g %g             %g %g\n", result[i].re, result[i].im, resultFast[i].re, resultFast[i].im);
	}

	printf("time fourier = %lld\n", time1);
	printf("time fast fourier = %lld\n", time2);

	printf("evaluating reverce fourier\n");
	time(&time0);
	double* b = evaluateReverceFourierTranslation(result, number);
	time(&time1);
	time1 = time1 - time0;
	printf("evaluating reverce fast fourier\n");
	time(&time0);
	double* c = fastFourierReverceTransition(resultFast, number);
	time(&time2);
	time2 = time2 - time0;

	for(int i = 0; i < number; ++i) {
		printf("%g %g %g\n", a[i], b[i], c[i]);
	}

	printf("time fourier = %lld\n", time1);
	printf("time fast fourier = %lld\n", time2);

	delete[] result;
	delete[] resultFast;

	delete[] a;
	delete[] b;
	delete[] c;
}

void testGMRES(){
	int rank;
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//bool periodic = false;
	bool periodic = true;
	const int numberGeneral = 10;
	int number = ((numberGeneral - 2)/nprocs) + 2;

	double** generalMatrix = new double*[numberGeneral];
	double** matrix = new double*[number];
	double* generalRightPart = new double[numberGeneral];
	double* rightPart = new double[number];

	for(int i = 0; i < numberGeneral; ++i){
		generalMatrix[i] = new double[numberGeneral];
		for(int j = 0; j < numberGeneral; ++j){
			if(i == j){
				generalMatrix[i][j] = -2;
			} else if((i == j-1) || (i == j+1)){
				generalMatrix[i][j] = 1;
			} else {
				generalMatrix[i][j] = 0;
			}
		}
	}
	generalMatrix[0][9] = 1;
	generalMatrix[9][0] = 1;
	/*generalMatrix[1][2] = 0.2;
	generalMatrix[2][3] = 0.1;
    generalMatrix[3][1] = 0.3;
    generalMatrix[3][4] = 0.1;
    generalMatrix[4][5] = -0.5;
    generalMatrix[5][6] = 0.1;
    generalMatrix[6][5] = 0.1;*/

	for(int i = 0; i < numberGeneral; ++i){
		generalRightPart[i] = 0;
	}
	generalRightPart[0] = -1;
	generalRightPart[1] = 1;
	generalRightPart[2] = -1;
	generalRightPart[3] = 1;
	generalRightPart[4] = -1;
	generalRightPart[5] = 1;
	generalRightPart[6] = -1;
	generalRightPart[7] = 1;
	generalRightPart[8] = -1;
	generalRightPart[9] = 1;

	if(rank == 0){
	for(int i = 0; i < numberGeneral; ++i) {
		for(int j = 0; j < numberGeneral; ++j) {
			printf("%g     ", generalMatrix[i][j]);
		}
		printf("   %g\n", generalRightPart[i]);
	}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for(int i = 0; i < number; ++i) {
		matrix[i] = new double[number];
        if((i == 0) && ((rank > 0)||(periodic))){
            for (int j = 0; j < number; ++j) {
                matrix[i][j] = 0;
            }
            matrix[i][i] = 1.0;
            rightPart[i] = 0;
        } else if((i == number - 1) && ((rank <nprocs - 1)||(periodic))){
            for (int j = 0; j < number; ++j) {
                matrix[i][j] = 0;
            }
            matrix[i][i] = 1.0;
            rightPart[i] = 0;
        } else{
            for (int j = 0; j < number; ++j) {
                matrix[i][j] = generalMatrix[i + rank * (number - 2)][j + rank * (number - 2)];
            }
            rightPart[i] = generalRightPart[i + rank * (number - 2)];
        }
	}

	for(int threadNumber = 0; threadNumber < nprocs; ++threadNumber) {
		if(threadNumber == rank) {
			for(int i = 0; i < number; ++i) {
				for(int j = 0; j < number; ++j) {
					printf("%g     ", matrix[i][j]);
				}
				printf("   %g\n", rightPart[i]);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	double* result = generalizedMinimalResidualMethod(matrix, rightPart, number, periodic, rank, nprocs);

	if(rank == 0) printf("result\n");
	MPI_Barrier(MPI_COMM_WORLD);
	for(int threadNumber = 0; threadNumber < nprocs; ++threadNumber) {
		if(threadNumber == rank) {
			for(int i = 0; i < number; ++i) {
				printf("%g\n", result[i]);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if(rank == 0) printf("error\n");
	MPI_Barrier(MPI_COMM_WORLD);
	for(int threadNumber = 0; threadNumber < nprocs; ++threadNumber) {
		if(threadNumber == rank) {
			for(int i = 0; i < number; ++i){
				double value = -rightPart[i];
				for(int j = 0; j < number; ++j){
					value += matrix[i][j]*result[j];
				}
				printf("%g\n", value);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

void testConjugate(){
	int rank;
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	bool periodic = false;
	const int number = 7;
	double** matrix = new double*[number];
	double* rightPart = new double[number];

	for(int i = 0; i < number; ++i){
		matrix[i] = new double[number];
		for(int j = 0; j < number; ++j){
			if(i == j){
				matrix[i][j] = 1;
			} else {
				matrix[i][j] = 0;
			}
		}
	}
	matrix[0][1] = 1;
	matrix[1][2] = 2;
	matrix[2][3] = 1;
	matrix[3][4] = 1;
	matrix[4][5] = -5;
	matrix[5][6] = 1;
	matrix[3][5] = 1;
	matrix[4][1] = 3;

	for(int i = 0; i < number; ++i){
		rightPart[i] = 0;
	}
	rightPart[number-1] = 1;


	for(int i = 0; i < number; ++i) {
		for(int j = 0; j < number; ++j) {
			printf("%lf     ", matrix[i][j]);
		}
		printf("   %lf\n", rightPart[i]);
	}

	double* result = new double[number];
	conjugateGradientMethod(matrix, rightPart, result, number, 1E-7, number, periodic, rank, nprocs);

	printf("result\n");
	for(int i = 0; i < number; ++i) {
		printf("%15.10lf\n", result[i]);
	}

	printf("error\n");
	for(int i = 0; i < number; ++i){
		double value = -rightPart[i];
		for(int j = 0; j < number; ++j){
			value += matrix[i][j]*result[j];
		}
		printf("%15.10lf\n", value);
	}
	delete[] result;
}

void testBiconjugate(){
	int rank;
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	bool periodic = false;
	const int number = 4;
	double** matrix = new double*[number];
	double* rightPart = new double[number];

	for(int i = 0; i < number; ++i){
		matrix[i] = new double[number];
		for(int j = 0; j < number; ++j){
			if(i == j){
				matrix[i][j] = 1;
			} else {
				matrix[i][j] = 0;
			}
		}
	}
	matrix[0][1] = 1;
	matrix[1][2] = 2;
	matrix[2][3] = 1;

	for(int i = 0; i < number; ++i){
		rightPart[i] = 0;
	}
	rightPart[number-1] = 1;


	for(int i = 0; i < number; ++i) {
		for(int j = 0; j < number; ++j) {
			printf("%lf     ", matrix[i][j]);
		}
		printf("   %lf\n", rightPart[i]);
	}

	double* result = new double[number];
	biconjugateGradientMethod(matrix, rightPart, result, number, 1E-7, number, periodic, rank, nprocs);

	printf("result\n");
	for(int i = 0; i < number; ++i) {
		printf("%15.10lf\n", result[i]);
	}

	printf("error\n");
	for(int i = 0; i < number; ++i){
		double value = -rightPart[i];
		for(int j = 0; j < number; ++j){
			value += matrix[i][j]*result[j];
		}
		printf("%15.10lf\n", value);
	}
	delete[] result;
}

void testBiconjugateStabilized(){
	int rank;
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	bool periodic = false;
	const int number = 7;
	double** matrix = new double*[number];
	double* rightPart = new double[number];

	for(int i = 0; i < number; ++i){
		matrix[i] = new double[number];
		for(int j = 0; j < number; ++j){
			if(i == j){
				matrix[i][j] = 1;
			} else {
				matrix[i][j] = 0;
			}
		}
	}
	matrix[0][1] = 1;
	matrix[1][2] = 2;
	matrix[2][3] = 1;
	matrix[3][4] = 1;
	matrix[4][5] = -5;
	matrix[5][6] = 1;
	matrix[3][5] = 1;
	matrix[4][1] = 3;

	for(int i = 0; i < number; ++i){
		rightPart[i] = 0;
	}
	rightPart[number-1] = 1;


	for(int i = 0; i < number; ++i) {
		for(int j = 0; j < number; ++j) {
			printf("%lf     ", matrix[i][j]);
		}
		printf("   %lf\n", rightPart[i]);
	}

	double* result = new double[number];
	biconjugateStabilizedGradientMethod(matrix, rightPart, result, number, 1E-7, number+1, periodic, rank, nprocs);

	printf("result\n");
	for(int i = 0; i < number; ++i) {
		printf("%15.10lf\n", result[i]);
	}

	printf("error\n");
	for(int i = 0; i < number; ++i){
		double value = -rightPart[i];
		for(int j = 0; j < number; ++j){
			value += matrix[i][j]*result[j];
		}
		printf("%15.10lf\n", value);
	}
	delete[] result;
}

void testSolve4order(){
	double a4 = 1;
	double a3 = 2;
	double a2 = 3;
	double a1 = 4;
	double a0 = -1;

	double startX = 1;

	printf("%lf * x^4 + %lf * x^3 + %lf * x^2 + %lf * x + %lf = 0\n", a4, a3, a2, a1, a0);

	double x = solve4orderEquation(a4, a3, a2, a1, a0, 2.0);

	printf("x = %lf\n", x);
}

void testFourier(){
	int xnumber = 32;
	int ynumber = 1;
	int znumber = 1;

	double xsize = xnumber;
	double ysize = ynumber;
	double zsize = znumber;

	double deltaX = xsize/xnumber;
	double deltaY = ysize/ynumber;
	double deltaZ = zsize/znumber;

	double*** a = new double**[xnumber];
	double*** da2 = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		a[i] = new double*[ynumber];
		da2[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			a[i][j] = new double[znumber];
			da2[i][j] = new double[znumber];
			for(int k = 0; k < znumber; ++k){
				a[i][j][k] = 0;
				a[i][j][k] = sin(2*pi*i*1.0/xnumber);
				da2[i][j][k] = 0;
			}
		}
	}
	//a[1][0][0] = 1;
	//a[0][0][0] = 1;
	//a[7][0][0] = 1;

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				int prevI = i - 1;
				if(prevI < 0){
					prevI = xnumber - 1;
				}
				int nextI = i + 1;
				if(nextI >= xnumber){
					nextI = 0;
				}
				int prevJ = j - 1;
				if(prevJ < 0){
					prevJ = ynumber - 1;
				}
				int nextJ = j + 1;
				if(nextJ >= ynumber){
					nextJ = 0;
				}
				int prevK = k - 1;
				if(prevK < 0){
					prevK = znumber - 1;
				}
				int nextK = k + 1;
				if(nextK >= znumber){
					nextK = 0;
				}

				//da2[i][j][k] = (-(2.0*a[i][j][k]/deltaX) - (2.0*a[i][j][k]/deltaY) - (2.0*a[i][j][k]/deltaZ) + ((a[prevI][j][k] + a[nextI][j][k])/deltaX) + ((a[i][prevJ][k] + a[i][nextJ][k])/deltaY) + ((a[i][j][prevK] + a[i][j][nextK])/deltaZ));
				da2[i][j][k] = (-(2.0*a[i][j][k]/deltaX) + ((a[prevI][j][k] + a[nextI][j][k])/deltaX));
			}
		}
	}

	Complex*** fourier = evaluateFourierTranslation(a, xnumber, ynumber, znumber);

	Complex*** fourier_da2 = new Complex**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		fourier_da2[i] = new Complex*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			fourier_da2[i][j] = new Complex[znumber];
			for(int k = 0; k < znumber; ++k){
				if(i < xnumber/2.0){
					fourier_da2[i][j][k] = fourier[i][j][k]*(-4*pi*pi*((i*i*1.0/(xsize*xsize))  + (j*j*1.0/(ysize*ysize)) + (k*k*1.0/(zsize*zsize))));
				} else {
					fourier_da2[i][j][k] = fourier[i][j][k]*(-4*pi*pi*(((xnumber - i)*(xnumber - i)*1.0/(xsize*xsize))  + (j*j*1.0/(ysize*ysize)) + (k*k*1.0/(zsize*zsize))));
				}
			}
		}
	}

	double*** b = evaluateReverceFourierTranslation(fourier, xnumber, ynumber, znumber);
	double*** db2 = evaluateReverceFourierTranslation(fourier_da2, xnumber, ynumber, znumber);
	

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				printf("a = %g Fa = %g b = %g da2 = %g db2 = %g", a[i][j][k], fourier[i][j][k].module(), b[i][j][k], da2[i][j][k], db2[i][j][k]);
			}
		}
		printf("\n");
	}


	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			delete[] a[i][j];
			delete[] da2[i][j];
			delete[] b[i][j];
			delete[] db2[i][j];
			delete[] fourier[i][j];
			delete[] fourier_da2[i][j];
		}
		delete[] a[i];
		delete[] da2[i];
		delete[] b[i];
		delete[] db2[i];
		delete[] fourier[i];
		delete[] fourier_da2[i];
	}
	delete[] a;
	delete[] da2;
	delete[] b;
	delete[] db2;
	delete[] fourier;
	delete[] fourier_da2;
}

void testFourierSolvePoison(){
	const int number = 2*2*2*2*2*2*2;
	double size = 1000;
	double* a = new double[number];
	double* b = new double[number];
	double* c = new double[number];


	for(int i = 0; i < number; ++i){
		//a[i] = cos(2*pi*i/number);
		a[i] = rand()%16 - 7.5;
		b[i] = - (4*pi*pi/(size*size))*cos(2*pi*i/number);
	}
	double delta = size/number;
	c[0] = (a[1] - 2*a[0] + a[number - 1])/(delta*delta);
	for(int i = 1; i < number - 1; ++i){
		c[i] = (a[i+1] - 2*a[i] + a[i-1])/(delta*delta);
	}
	c[number - 1] = (a[0] - 2*a[number - 1] + a[number - 2])/(delta*delta);

	Complex* rightPartFourier = fastFourierTransition(c,number);

	for(int i = 0; i < number; ++i){
		if(i == 0) {
			rightPartFourier[i] = Complex(0, 0);
			//} else if(i < number/2){
		}else {
			rightPartFourier[i] = rightPartFourier[i]*2 / (-4 * pi * pi * i * i / (size * size));
		}
		/*} else {
			rightPartFourier[i] = rightPartFourier[i]/(-4*pi*pi*(number - 1)*(number - 1)/(size*size));
		}*/
	}

	double* fourierResult = fastFourierReverceTransition(rightPartFourier, number);

	for(int i = 0; i < number; ++i){
		printf("%g %g\n", a[i], fourierResult[i]);
	}

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] rightPartFourier;
	delete[] fourierResult;
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	testGMRES();
	//testBiconjugate();
	//testBiconjugateStabilized();

	//test solution of a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0
	//testSolve4order();

	//testFourier();
	//testFastFourier();
	//testFourierSolvePoison();
	MPI_Finalize();
	return 0;
}

