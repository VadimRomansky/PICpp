#ifndef _RANDOM_H_
#define _RANDOM_H_

double uniformDistribution();
double normalDistribution();
double maxwellDistribution(const double& temperature, const double& k);
double maxwellJuttnerDistribution(const double& temperature, const double& mass, const double& c, const double& k);
double maxwellJuttnerMomentumDistribution(const double& temperature, const double& mass, const double& c, const double& k);
double maxwellJuttnerMomentumColdDistribution(const double& temperature, const double& mass, const double& c, const double& k, const double* function, const double* xvalue, int number);
double maxwellJuttnerMomentumHotDistribution(const double& temperature, const double& mass, const double& c, const double& k);
double solveInverceJuttnerFunction(const double& x, const double& theta, const double& besselK);
double solveInverceJuttnerFunction(const double& x, const double& theta, const double& besselK, const double& left, const double& right);
double maxwellJuttnerFunction(const double& gamma, const double& theta, const double& besselK);
double maxwellJuttnerIntegral(const double& gamma, const double& theta, const double& besselK);
void anisotropicMaxwellJuttnerDistribution(double &momentumNormal, double &momentumParallel, double alphaNormal,
                                           double alphaParallel, double m_c);
double anisotropicMaxwellJuttnerFunction(const double& p1, const double& p2, const double& alpha1, const double& alpha2);

#endif