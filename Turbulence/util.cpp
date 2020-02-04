#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "constants.h"
#include "util.h"

double power(const double& v, const double& p) {
	return exp(p * log(v));
}

double sqr(const double& v) {
	return v * v;
}

double cube(const double& v) {
	return v * v * v;
}

double max2(const double& a, const double& b) {
	if (a >= b) {
		return a;
	}
	return b;
}

double max3(const double& a, const double& b, const double& c) {
	if (a > b) {
		return max2(a, c);
	}
	return max2(b, c);
}

double min2(const double& a, const double& b) {
	if (a >= b) {
		return b;
	}
	return a;
}

double min3(const double& a, const double& b, const double& c) {
	if (a > b) {
		return min2(b, c);
	}
	return min2(a, c);
}

void alertNaNOrInfinity(const double& value, const char* s) {
	if (value != value || 0 * value != 0 * value) {
		printf("%s", s);
		printf("\n");
		//FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		FILE* errorLogFile = fopen("errorLog.dat", "w");
		fprintf(errorLogFile, "%s", s);
		exit(0);
		//Sleep(1000);
	}
}

void alertNotPositive(double& value, const char* s) {
	if (value <= 0) {
		printf("%s", s);
		printf("\n");
		//FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		FILE* errorLogFile = fopen("errorLog.dat", "w");
		fprintf(errorLogFile, "%s", s);
		fclose(errorLogFile);
		exit(0);
	}
}

void alertNegative(double& value, const char* s) {
	if (value < 0) {
		printf("%s", s);
		printf("\n");
		//FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		FILE* errorLogFile = fopen("errorLog.dat", "w");
		fprintf(errorLogFile, "%s", s);
		fclose(errorLogFile);
		exit(0);
	}
}

void alertLargeModlue(double& value, double compareValue, const char* s) {
	if (fabs(value) > compareValue) {
		printf("%s", s);
		printf("\n");
		//FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		FILE* errorLogFile = fopen("errorLog.dat", "w");
		fprintf(errorLogFile, "%s", s);
		int rank;
		exit(0);
		//Sleep(1000);
	}
}

void printErrorAndExit(const char* s) {
	printf("%s", s);
	printf("\n");
	//FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
	FILE* errorLogFile = fopen("errorLog.dat", "w");
	fprintf(errorLogFile, "%s", s);
	fclose(errorLogFile);
	exit(0);
}

void printLog(const char* s) {
	int size;
	int rank;
	//FILE* logFile = fopen("./output/log.dat", "a");
	FILE* logFile = fopen("log.dat", "a");
	fprintf(logFile, "rank = %d, message %s", rank, s);
	fflush(logFile);
	fclose(logFile);
}

/*int trunc(double value){
	int round = value;
	if(round > value){
		round--;
	}
	return round;
}*/


double coordinateDifference(double* const a, double* const b, double dt, double mass) {
	double result = 0;
	for (int i = 0; i < 3; ++i) {
		result += fabs(a[i] - b[i]);
	}

	for (int i = 3; i < 6; ++i) {
		result += fabs(a[i] - b[i]) * dt / mass;
	}
	return result;
}


double McDonaldFunction(double x, double index) {
	//todo approximation with small x!!!
	if (x < 0) {
		printf("x in McDonald < 0\n");
		//FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		FILE* errorLogFile = fopen("errorLog.dat", "w");
		fprintf(errorLogFile, "x in McDonald < 0\n");
		fclose(errorLogFile);
		exit(0);
	}

	if (x > 2 * index * index && x > 10) {
		double result = sqrt(pi / (2 * x)) * exp(-x) * (1 + ((4 * index * index - 1) / (8 * x)) + ((4 * index * index - 1) * (
			4 * index * index - 9) / (128 * x * x)) + ((4 * index * index - 1) * (4 * index * index - 9) * (4 * index * index -
			25) / (6 * cube(8 * x))));
		return result;
	}
	double dt;
	double t;
	double prevT = 0;
	double result = 0;
	double maxT;
	double maxLevel = 10;

	double discr = index * index - 2 * x * x + 2 * x * maxLevel;
	if (discr < 0) {
		maxT = log(maxLevel * maxLevel + sqrt(maxLevel * maxLevel - 1));
	} else {
		maxT = (index + sqrt(discr)) / (x);
	}
	dt = maxT / 100;

	while (x * cosh(maxT) - index * maxT < 10) {
		maxT = maxT + dt;
	}
	dt = maxT / 1000000;

	t = dt;
	int i = 0;
	while (i < 1000000) {
		double middleT = 0.5 * (t + prevT);
		double dresult = exp(-x * cosh(middleT)) * cosh(index * middleT) * dt;
		result += dresult;
		if (dresult < result / 1E14) break;
		prevT = t;
		t += dt;
		++i;
	}

	return result;
}

double Bspline(double xcenter, double dx, double xvalue) {
	if (fabs(xcenter - xvalue) > dx) {
		return 0;
	}

	if (xvalue > xcenter + dx / 2) {
		return 2 * sqr(xcenter + dx - xvalue) / cube(dx);
	}

	if (xvalue < xcenter - dx / 2) {
		return 2 * sqr(xcenter - dx - xvalue) / cube(dx);
	}

	return (1 / dx) - 2 * sqr(xvalue - xcenter) / cube(dx);
}

std::string convertIntToString(int a) {
	if (a == 0) {
		std::string result = "0";
		return result;
	}
	if (a > 0) {
		std::string result = "";
		while (a > 0) {
			int last = a % 10;
			a = a / 10;
			char c = last + '0';
			result = c + result;
		}
		return result;
	}
	a = -a;
	std::string result = "-";
	return result + convertIntToString(a);
}


