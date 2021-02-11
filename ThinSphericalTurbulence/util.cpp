#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include <time.h>
#include <cmath>
#include <string>
#include <omp.h>

#include "startparameters.h"
#include "constants.h"

#include "util.h"

double uniformDistribution() {
	return (rand() % randomParameter + 0.5) / randomParameter;
}

double sqr(const double& a){
	return a*a;
}

double min(const double& a, const double& b){
	if(a < b) {
		return a;
	} else {
		return b;
	}
}

double min3(const double& a, const double& b, const double& c){
	if(a < b) {
		return min(a, c);
	} else {
		return min(b, c);
	}
}

double min4(const double& a, const double& b, const double& c, const double& d){
	if(a < b) {
		if(c < d){
			return min(a, c);
		} else {
			return min(a, d);
		}
	} else {
		if(c < d){
			return min(b, c);
		} else {
			return min(b, d);
		}
	}
}

double max(const double& a, const double& b){
	if(a > b) {
		return a;
	} else {
		return b;
	}
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