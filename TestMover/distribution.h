#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

void createParticle(double& px, double& py, double& pz, double temperature, double mass, double* juttnerValue, double* juttnerFunction, int juttnerN);
void createFastParticle(double& px, double& py, double& pz, double mass, double gamma);
double uniformDistribution();
double normalDistribution();
double maxwellDistribution(const double& temperature, const double& k);
double maxwellJuttnerMomentumColdDistribution(const double& temperature, const double& mass, const double& k, const double* function, const double* xvalue, int number);
double maxwellJuttnerMomentumHotDistribution(const double& temperature, const double& mass, const double& k);
double dichotomySolver(const double* functionValue, int minIndex, int maxIndex, const double* xValue, double y);
double juttnerDistribution(const double& u, const double& theta);
void evaluateJuttnerFunction(double* juttnerValue, double* juttnerFunction, double temperature, double mass, int juttnerNumber);

#endif