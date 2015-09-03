#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include "util.h"
#include "constants.h"
#include "random.h"

double uniformDistribution(){
	return (rand()%randomSeed + 0.5)/randomSeed;
}

double normalDistribution(){
	double x = uniformDistribution();
	double y = uniformDistribution();
	return cos(2*pi*x)*sqrt(-2*log(y));
}

double maxwellDistribution(double temperature, double k){
	double normal = normalDistribution();
	return k*temperature*(0.5*normal*normal - log(uniformDistribution()));
}

double maxwellJuttnerDistribution(double temperature, double mass, double c, double k){
	double theta = k*temperature/(mass*c*c);
	double besselK = McDonaldFunction(1.0/theta, 2.0);

	double x = uniformDistribution();

	double gamma = solveInverceJuttnerFunction(x, theta, besselK);

	return gamma*mass*c*c;
}

double solveInverceJuttnerFunction(double x, double theta, double besselK){
	if(x >= 1){
		printf("distribution function can not be more then 1\n");
		exit(0);
	}
	if(x <= 0) {
		return 1;
	}
	double gamma = max2(1,theta);
	double integral = maxwellJuttnerIntegral(gamma, theta, besselK);
	if(integral < x){
		while(integral < x){
			gamma *= 2;
			integral = maxwellJuttnerIntegral(gamma, theta, besselK);
		}
	} 
	return solveInverceJuttnerFunction(x, theta, besselK, 0.0, gamma);
}

double solveInverceJuttnerFunction(double x, double theta, double besselK, double left, double right){
	if(right < left){
		printf("right < left\n");
		exit(0);
	}
	if( right - left < 0.000001*left){
		return (right + left)/2;
	}
	double gamma = (left + right)/2;

	double integral = maxwellJuttnerIntegral(gamma, theta, besselK);

	if(integral > x){
		return solveInverceJuttnerFunction(x, theta, besselK, left, gamma);
	} else {
		return solveInverceJuttnerFunction(x, theta, besselK, gamma, right);
	}
}

double maxwellJuttnerFunction(double gamma, double theta, double besselK){
	double beta = sqrt(1.0 - 1.0/(gamma*gamma));
	return gamma*gamma*beta*exp(-gamma/theta)/(theta*besselK);
}

double maxwellJuttnerIntegral(double gamma, double theta, double besselK){
	if(gamma < 1){
		return 0;
	}
	int count = 1000;
	double dgamma = (gamma - 1)/count;
	double sum = 0;
	double tempGamma = 1 + dgamma/2;
	for(int i = 0; i < count; ++i){
		sum += maxwellJuttnerFunction(tempGamma, theta, besselK)*dgamma;
		tempGamma += dgamma;
	}
	return sum;
}