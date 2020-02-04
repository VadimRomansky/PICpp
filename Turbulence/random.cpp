#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "util.h"
#include "constants.h"
#include "random.h"

double uniformDistribution() {
	return (rand() % randomParameter + 0.5) / randomParameter;
}

double normalDistribution() {
	//Box–Muller transform
	double x = uniformDistribution();
	double y = uniformDistribution();
	return cos(2 * pi * x) * sqrt(-2 * log(y));
}

double maxwellDistribution(const double& temperature, const double& k) {
	double normal = normalDistribution();
	return k * temperature * (0.5 * normal * normal - log(uniformDistribution()));
}
