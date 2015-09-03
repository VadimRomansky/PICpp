#ifndef _SPECIAL_MATH_H_
#define _SPECIAL_MATH_H_

#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "matrix3d.h"
#include "matrixElement.h"

void generalizedMinimalResidualMethod(std::vector<MatrixElement>** matrix, double** rightPart, double** outvector, int xnumber, int lnumber);
double*** arnoldiIterations(std::vector<MatrixElement>** matrix,double** outHessenbergMatrix, int n, double*** prevBasis, double** prevHessenbergMatrix, int xnumber, int lnumber);
double** multiplySpecialMatrixVector(std::vector<MatrixElement>** matrix, double** vector, int xnumber, int lnumber);
double** multiplySpecialMatrixVector(std::vector<MatrixElement>** matrix, Vector3d* vector, int xnumber, int lnumber);
double evaluateError(double** hessenbergMatrix, double* vector, double beta, int n);
double scalarMultiplyLargeVectors(double** a, double** b, int xnumber, int lnumber);
double scalarMultiplyLargeVectors(Vector3d* a, Vector3d* b, int xnumber, int lnumber);

double solve4orderEquation(double a4, double a3, double a2, double a1, double a0, double startX);
double polynom4Value(double a4, double a3, double a2, double a1, double a0, double x);
double polynom4DerivativeValue(double a4, double a3, double a2, double a1, double x);

double solveBigOrderEquation(const double* coefficients, int n, double startX);
double polynomValue(const double* coefficients, double x, int n);
double polynomDerivativeValue(const double* coefficients, double x, int n);

#endif