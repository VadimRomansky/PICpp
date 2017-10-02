#ifndef _SPECIAL_MATH_H_
#define _SPECIAL_MATH_H_

#include "complex.h"


const double pi = 3.1415926535897932384626433832795028841971693993751;

double** arnoldiIterations(double** matrix, double** outHessenbergMatrix, int n, double** prevBasis, double** prevHessenbergMatrix, int number);
double* generalizedMinimalResidualMethod(double** matrix, double* rightPart, int number);
double scalarMultiplyLargeVectors(double* a, double* b, int number);
double* multiplyMatrixVector(double** matrix, double* vector, int number);
void multiplyMatrixVector(double* result, double** matrix, double* vector, int number);

#endif