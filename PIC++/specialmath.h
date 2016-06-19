#ifndef _SPECIAL_MATH_H_
#define _SPECIAL_MATH_H_

#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "matrix3d.h"
#include "matrixElement.h"

void generalizedMinimalResidualMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outvector, int xnumber, int ynumber, int znumber, int lnumber, double precision, int maxIteration);
double***** arnoldiIterations(std::vector<MatrixElement>**** matrix,double** outHessenbergMatrix, int n, double***** prevBasis, double** prevHessenbergMatrix, int xnumber, int ynumber, int znumber, int lnumber);
double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, double**** vector, int xnumber, int ynumber, int znumber, int lnumber);
double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, Vector3d*** vector, int xnumber, int ynumber, int znumber, int lnumber);
void multiplySpecialMatrixVector(double**** result, std::vector<MatrixElement>**** matrix, Vector3d*** vector, int xnumber, int ynumber, int znumber, int lnumber);
void multiplySpecialMatrixVector(double**** result, std::vector<MatrixElement>**** matrix, double**** vector, int xnumber, int ynumber, int znumber, int lnumber);
double evaluateError(double**** hessenbergMatrix, double*** vector, double beta, int n);
double scalarMultiplyLargeVectors(double**** a, double**** b, int xnumber, int ynumber, int znumber, int lnumber);
double scalarMultiplyLargeVectors(Vector3d*** a, Vector3d*** b, int xnumber, int ynumber, int znumber, int lnumber);

void conjugateGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int xnumber, int ynumber, int znumber, int lnumber, double precision, int maxIteration);

void biconjugateGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int xnumber, int ynumber, int znumber, int lnumber, double precision, int maxIteration);
void biconjugateStabilizedGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int xnumber, int ynumber, int znumber, int lnumber, double precision, int maxIteration);
void transposeSpecialMatrix(std::vector<MatrixElement>**** result, std::vector<MatrixElement>**** matrix, int xnumber, int ynumber, int znumber, int lnumber);

double solve4orderEquation(double a4, double a3, double a2, double a1, double a0, double startX);
double polynom4Value(double a4, double a3, double a2, double a1, double a0, double x);
double polynom4DerivativeValue(double a4, double a3, double a2, double a1, double x);

double solveBigOrderEquation(const double* coefficients, int n, double startX);
double polynomValue(const double* coefficients, double x, int n);
double polynomDerivativeValue(const double* coefficients, double x, int n);

void sortInputFastFourierX(Complex ***a, int xnumber, int ynumber, int znumber);
void sortInputFastFourierY(Complex ***a, int xnumber, int ynumber, int znumber);
void sortInputFastFourierZ(Complex ***a, int xnumber, int ynumber, int znumber);
void sortInputFastFourierX(double ***a, int xnumber, int ynumber, int znumber);
void sortInputFastFourierY(double ***a, int xnumber, int ynumber, int znumber);
void sortInputFastFourierZ(double ***a, int xnumber, int ynumber, int znumber);
double ***fastFourierReverceTransition(Complex ***a, int xnumber, int ynumber, int znumber);
Complex ***fastFourierTransition(double ***a, int xnumber, int ynumber, int znumber);

#endif