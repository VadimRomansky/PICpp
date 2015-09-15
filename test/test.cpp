// test.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "specialmath.h"


int main()
{
	//test GMRES
	double** matrix = new double*[number];
	double* rightPart = new double[number];

	for(int i = 0; i < number; ++i){
		matrix[i] = new double[number];
		for(int j = 0; j < number; ++j){
			if(i == j){
				matrix[i][j] = 0.1;
			} else {
				matrix[i][j] = 0;
			}
		}
	}
	matrix[0][1] = 1;
	matrix[1][2] = 1;
	matrix[2][3] = 1;
	matrix[3][4] = 1;
	matrix[4][5] = 1;
	matrix[5][6] = 1;

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


	//test solution of a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0

	/*double a4 = 1;
	double a3 = 2;
	double a2 = 3;
	double a1 = 4;
	double a0 = -1;

	double startX = 1;

	printf("%lf * x^4 + %lf * x^3 + %lf * x^2 + %lf * x + %lf = 0\n", a4, a3, a2, a1, a0);

	double x = solve4orderEquation(a4, a3, a2, a1, a0, 2.0);

	printf("x = %lf\n", x);

	return 0;*/
}

