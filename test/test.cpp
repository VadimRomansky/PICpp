// test.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "specialmath.h"
#include "complex.h"

void testGMRES(){
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

	double* result = generalizedMinimalResidualMethod(matrix, rightPart);

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
	int xnumber = 25;
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
				da2[i][j][k] = 0;
			}
		}
	}
	a[1][0][0] = 1;
	a[0][0][0] = 1;
	a[7][0][0] = 1;

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

				da2[i][j][k] = (-(2.0*a[i][j][k]/deltaX) - (2.0*a[i][j][k]/deltaY) - (2.0*a[i][j][k]/deltaZ) + ((a[prevI][j][k] + a[nextI][j][k])/deltaX) + ((a[i][prevJ][k] + a[i][nextJ][k])/deltaY) + ((a[i][j][prevK] + a[i][j][nextK])/deltaZ));
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
				fourier_da2[i][j][k] = fourier[i][j][k]*(-4*pi*pi*((i*i*1.0/(xsize*xsize))  + (j*j*1.0/(ysize*ysize)) + (k*k*1.0/(zsize*zsize))));
			}
		}
	}

	double*** b = evaluateReverceFourierTranslation(fourier, xnumber, ynumber, znumber);
	double*** db2 = evaluateReverceFourierTranslation(fourier_da2, xnumber, ynumber, znumber);
	

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				printf("a = %g Fa = %g b = %g da2 = %g db2 = %g", a[i][j][k], fourier[i][j][k].re, b[i][j][k], da2[i][j][k], db2[i][j][k]);
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

int main()
{
	//test GMRES
	//testGMRES();

	//test solution of a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0
	//testSolve4order();

	testFourier();

	return 0;
}

