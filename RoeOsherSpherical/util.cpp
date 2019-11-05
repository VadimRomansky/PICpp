#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "constants.h"

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

double min2(double a, double b){
	if(a >= b){
		return b;
	} else {
		return a;
	}
}

double abs2(double a){
    if(a >= 0.0){
        return a;
    } else {
        return -a;
    }
}

void alertNaNOrInfinity(double value, const char* s){
	if(value != value || 0*value != 0*value){
		printf("%s", s);
		printf("\n");
		exit(0);
		//Sleep(1000);
	}
}

void alertNotPositive(double value, const char* s){
	if(value <= 0){
		printf("%s", s);
		printf("\n");
		exit(0);
	}
}

void alertNegative(double value, const char* s){
	if(value < 0){
		printf("%s", s);
		printf("\n");
		exit(0);
	}
}

double findExpLevel(double y, int N, double min, double max){
	double middle = (max + min)/2;
	if(max - min < epsilon) return middle;
	if((exp(N*log(middle)))/(middle - 1) > y){
		return findExpLevel(y, N, min, middle);
	} else {
		return findExpLevel(y, N, middle, max);
	}
}
