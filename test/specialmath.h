#ifndef _SPECIAL_MATH_H_
#define _SPECIAL_MATH_H_

#include "complex.h"

const int number = 7;

const double pi = 3.1415926535897932384626433832795028841971693993751;

double** arnoldiIterations(double** matrix, double** outHessenbergMatrix, int n, double** prevBasis, double** prevHessenbergMatrix);
double* generalizedMinimalResidualMethod(double** matrix, double* rightPart);
double scalarMultiplyLargeVectors(double* a, double* b);
double* multiplyMatrixVector(double** matrix, double* vector);

double solve4orderEquation(double a4, double a3, double a2, double a1, double a0, double startX);
double polynomValue(double a4, double a3, double a2, double a1, double a0, double x);
double polynomDerivativeValue(double a4, double a3, double a2, double a1, double x);

Complex*** evaluateFourierTranslation(double*** a, int xnumber, int ynumber, int znumber);
double*** evaluateReverceFourierTranslation(Complex*** a, int xnumber, int ynumber, int znumber);

#endif