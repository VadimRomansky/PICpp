#ifndef UTIL_H
#define UTIL_H

double power(double v, double p);
double sqr(double v);
double cube(double v);
double max2(double a, double b);
double max3(double a, double b, double c);
double min2(double a, double b);
double min3(double a, double b, double c);
void alertNaNOrInfinity(double value, const char* s);
void alertNotPositive(double value, const char* s);
void alertNegative(double value, const char* s);
//int trunc(double value);
double Bspline(double xcenter, double dx, double xvalue);

double coordinateDifference(double* const a, double* const b, double dt, double mass);
void solveSpecialMatrix(double** const leftHalf, double* const rightPart, double* const output);

double* linearLeastSquares(double** matrix, int n);

double McDonaldFunction(double x, double index);


#endif
