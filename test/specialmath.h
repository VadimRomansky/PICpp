#ifndef _SPECIAL_MATH_H_
#define _SPECIAL_MATH_H_

const int number = 7;

double** arnoldiIterations(double** matrix, double** outHessenbergMatrix, int n, double** prevBasis, double** prevHessenbergMatrix);
double* generalizedMinimalResidualMethod(double** matrix, double* rightPart);
double scalarMultiplyLargeVectors(double* a, double* b);
double* multiplyMatrixVector(double** matrix, double* vector);

double solve4orderEquation(double a4, double a3, double a2, double a1, double a0, double startX);
double polynomValue(double a4, double a3, double a2, double a1, double a0, double x);
double polynomDerivativeValue(double a4, double a3, double a2, double a1, double x);

#endif