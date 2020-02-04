#ifndef UTIL_H
#define UTIL_H

#include <string>

double power(const double & v, const double & p);
double sqr(const double& v);
double cube(const double & v);
double max2(const double & a, const double & b);
double max3(const double & a, const double & b, const double & c);
double min2(const double & a, const double & b);
double min3(const double & a, const double & b, const double & c);
void alertNaNOrInfinity(const double & value, const char* s);
void alertNotPositive(double& value, const char* s);
void alertNegative(double& value, const char* s);
void alertLargeModlue(double& value, double compareValue, const char* s);
void printErrorAndExit(const char* s);
void printLog(const char* s);
//int trunc(double value);
double Bspline(double xcenter, double dx, double xvalue);

std::string convertIntToString(int a);


#endif
