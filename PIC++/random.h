#ifndef _RANDOM_H_
#define _RANDOM_H_

double uniformDistribution();
double normalDistribution();
double maxwellDistribution(double temperature, double k);
double maxwellJuttnerDistribution(double temperature, double mass, double c, double k);
double maxwellJuttnerMomentumDistribution(double temperature, double mass, double c, double k);
double maxwellJuttnerMomentumColdDistribution(double temperature, double mass, double c, double k, const double* function, const double* xvalue, int number);
double maxwellJuttnerMomentumHotDistribution(double temperature, double mass, double c, double k);
double solveInverceJuttnerFunction(double x, double theta, double besselK);
double solveInverceJuttnerFunction(double x, double theta, double besselK, double left, double right);
double maxwellJuttnerFunction(double gamma, double theta, double besselK);
double maxwellJuttnerIntegral(double gamma, double theta, double besselK);
void anisotropicMaxwellJuttnerDistribution(double &momentumNormal, double &momentumParallel, double alphaNormal,
                                           double alphaParallel, double m_c);
double anisotropicMaxwellJuttnerFunction(double p1, double p2, double alpha1, double alpha2);

#endif