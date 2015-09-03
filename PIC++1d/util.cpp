#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

double power(double v, double p){
	return exp(p*log(v));
}

double sqr(double v){
	return v*v;
}

double cube(double v){
	return v*v*v;
}

double max2(double a, double b){
	if(a >= b){
		return a;
	} else {
		return b;
	}
}

double max3(double a, double b, double c){
	if(a > b){
		return max2(a, c);
	} else {
		return max2(b, c);
	}
}

double min2(double a, double b){
	if(a >= b){
		return b;
	} else {
		return a;
	}
}

double min3(double a, double b, double c){
	if(a > b){
		return min2(b, c);
	} else {
		return min2(a, c);
	}
}

void alertNaNOrInfinity(double value, const char* s){
	if(value != value || 0*value != 0*value){
		printf("%s", s);
		printf("\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "%s", s);
		fclose(errorLogFile);
		exit(0);
		//Sleep(1000);
	}
}

void alertNotPositive(double value, const char* s){
	if(value <= 0){
		printf("%s", s);
		printf("\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "%s", s);
		fclose(errorLogFile);
		exit(0);
	}
}

void alertNegative(double value, const char* s){
	if(value < 0){
		printf("%s", s);
		printf("\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "%s", s);
		fclose(errorLogFile);
		exit(0);
	}
}

/*int trunc(double value){
	int round = value;
	if(round > value){
		round--;
	}
	return round;
}*/

void solveSpecialMatrix(double** const leftHalf, double* const rightPart, double* const output){
	for(int j = 0; j < 2; ++j){
		leftHalf[2][j] /= leftHalf[2][2];
	}
	rightPart[2] /= leftHalf[2][2];
	leftHalf[2][2] = 1;

	for(int i = 0; i < 2; ++i){
		for(int j = 0; j < 2; ++j){
			leftHalf[i][j] -= leftHalf[2][j]*leftHalf[i][2];
		}
		rightPart[i] -= rightPart[2]*leftHalf[i][2];
		leftHalf[i][2] = 0;
	}

	leftHalf[1][0] /= leftHalf[1][1];
	rightPart[1] /= leftHalf[1][1];
	leftHalf[1][1] = 1;

	leftHalf[0][0] -= leftHalf[1][0]*leftHalf[0][1];
	rightPart[0] -= rightPart[1]*leftHalf[0][1];
	leftHalf[0][1] = 0;

	rightPart[0] /= leftHalf[0][0];

	leftHalf[0][0] = 1;

	output[0] = rightPart[0];
	for(int i = 1; i < 6; ++i){
		output[i] = rightPart[i];
		for(int j = 0; j < min2(i,3); ++j){
			output[i] -= leftHalf[i][j]*output[j];
		}
	}
}

double coordinateDifference(double* const a, double* const b, double dt, double mass){
	double result = 0;
	for(int i = 0; i < 3; ++i){
		result += fabs(a[i] - b[i]);
	}

	for(int i = 3; i < 6; ++i){
		result += fabs(a[i] - b[i])*dt/mass;
	}
	return result;
}


double McDonaldFunction(double x, double index){
	//todo approximation with small x!!!
	if(x < 0){
		printf("x in McDonald < 0\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "x in McDonald < 0\n");
		fclose(errorLogFile);
		exit(0);
	}
	double dt;
	double t;
	double prevT = 0;
	double result = 0;
	double maxT;
	double maxLevel = 10;

	double discr = index*index - 2*x*x + 2*x*maxLevel;
	if(discr < 0){
		maxT = log(maxLevel*maxLevel + sqrt(maxLevel*maxLevel - 1));
	} else {
		maxT = (index + sqrt(discr))/(x);
	}
	dt = maxT/100;

	while(x*cosh(maxT) - index*maxT < 10){
		maxT = maxT + dt;
	}
	dt = maxT/10000;

	t = dt;
	int i = 0;
	while(i < 10000){
		double middleT = 0.5*(t + prevT);
		double dresult = exp(-x*cosh(middleT))*cosh(index*middleT)*dt;
		result += dresult;
		if(dresult < result/1E10) break;
		prevT = t;
		t += dt;
		++i;
	}

	return result;
}

double Bspline(double xcenter, double dx, double xvalue) {
	if(fabs(xcenter - xvalue) > dx) {
		return 0;
	}

	if(xvalue > xcenter + dx/2) {
		return 2*sqr(xcenter + dx - xvalue)/cube(dx);
	}

	if(xvalue < xcenter - dx/2) {
		return 2*sqr(xcenter - dx - xvalue)/cube(dx);
	}

	return (1/dx) - 2*sqr(xvalue - xcenter)/cube(dx);
}


