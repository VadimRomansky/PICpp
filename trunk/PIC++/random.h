#ifndef _RANDOM_H_
#define _RANDOM_H_

double uniformDistribution();
double normalDistribution();
double maxwellDistribution(double temperature, double k);
double maxwellJuttnerDistribution(double temperature, double mass, double c, double k);
double solveInverceJuttnerFunction(double x, double theta, double besselK);
double solveInverceJuttnerFunction(double x, double theta, double besselK, double left, double right);
double maxwellJuttnerFunction(double gamma, double theta, double besselK);
double maxwellJuttnerIntegral(double gamma, double theta, double besselK);

#endif