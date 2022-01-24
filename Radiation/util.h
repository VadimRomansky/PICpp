#ifndef UTIL_H
#define UTIL_H

#include <string>

double uniformDistribution();
double sqr(const double& a);
double cube(const double& a);
double min(const double& a, const double& b);
double min3(const double& a, const double& b, const double& c);
double min4(const double& a, const double& b, const double& c, const double& d);
double max(const double& a, const double& b);
std::string convertIntToString(int a);

#endif