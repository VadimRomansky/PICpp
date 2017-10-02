#ifndef _SPECIAL_MATH_H_
#define _SPECIAL_MATH_H_

#include "complex.h"


const double pi = 3.1415926535897932384626433832795028841971693993751;

double** arnoldiIterations(double** matrix, double** outHessenbergMatrix, int n, double** prevBasis, double** prevHessenbergMatrix, int number, bool periodic, int rank, int nprocs);
double* generalizedMinimalResidualMethod(double** matrix, double* rightPart, int number, bool periodic, int rank, int nprocs);
double scalarMultiplyLargeVectors(double* a, double* b, int number, bool periodic, int rank, int nprocs);
double* multiplyMatrixVector(double** matrix, double* vector, int number, bool periodic, int rank, int nprocs);
void multiplyMatrixVector(double* result, double** matrix, double* vector, int number, bool periodic, int rank, int nprocs);

void conjugateGradientMethod(double** matrix, double* rightPart, double* outVector, int number, double precision, int maxIterationr, bool periodic, int rank, int nprocs);
void biconjugateGradientMethod(double** matrix, double* rightPart, double* outVector, int number, double precision, int maxIterationr, bool periodic, int rank, int nprocs);
void biconjugateStabilizedGradientMethod(double** matrix, double* rightPart, double* outVector, int number, double precision, int maxIterationr, bool periodic, int rank, int nprocs);
void transposeMatrix(double** result, double** matrix, int number);

double solve4orderEquation(double a4, double a3, double a2, double a1, double a0, double startX);
double polynomValue(double a4, double a3, double a2, double a1, double a0, double x);
double polynomDerivativeValue(double a4, double a3, double a2, double a1, double x);

Complex*** evaluateFourierTranslation(double*** a, int xnumber, int ynumber, int znumber);
double*** evaluateReverceFourierTranslation(Complex*** a, int xnumber, int ynumber, int znumber);

Complex* evaluateFourierTranslation(double* a, int number);
double* evaluateReverceFourierTranslation(Complex* a, int number);

Complex* fastFourierTransition(double* a, int n);
void sortInputFastFourier(double* a, int n);

double* fastFourierReverceTransition(Complex* a, int n);
void sortInputFastFourierReverce(Complex* a, int n);

#endif