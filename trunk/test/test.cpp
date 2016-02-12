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
	double*** a = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		a[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			a[i][j] = new double[znumber];
			for(int k = 0; k < znumber; ++k){
				a[i][j][k] = 0;
			}
		}
	}
	a[1][0][0] = 1;
	a[0][0][0] = 1;
	a[7][0][0] = 1;

	Complex*** fourier = evaluateFourierTranslation(a, xnumber, ynumber, znumber);

	double*** b = evaluateReverceFourierTranslation(fourier, xnumber, ynumber, znumber);
	

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				printf("a = %g Fa = %g b = %g", a[i][j][k], fourier[i][j][k].re, b[i][j][k]);
			}
		}
		printf("\n");
	}


	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			delete[] a[i][j];
			delete[] b[i][j];
			delete[] fourier[i][j];
		}
		delete[] a[i];
		delete[] b[i];
		delete[] fourier[i];
	}
	delete[] a;
	delete[] b;
	delete[] fourier;
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

