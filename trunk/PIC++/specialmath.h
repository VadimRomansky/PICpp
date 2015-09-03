#ifndef _SPECIAL_MATH_H_
#define _SPECIAL_MATH_H_

#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "matrix3d.h"
#include "matrixElement.h"

void generalizedMinimalResidualMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outvector, int xnumber, int ynumber, int znumber, int lnumber);
double***** arnoldiIterations(std::vector<MatrixElement>**** matrix,double** outHessenbergMatrix, int n, double***** prevBasis, double** prevHessenbergMatrix, int xnumber, int ynumber, int znumber, int lnumber);
double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, double**** vector, int xnumber, int ynumber, int znumber, int lnumber);
double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, Vector3d*** vector, int xnumber, int ynumber, int znumber, int lnumber);
double evaluateError(double** hessenbergMatrix, double* vector, double beta, int n);
double scalarMultiplyLargeVectors(double**** a, double**** b, int xnumber, int ynumber, int znumber, int lnumber);
double scalarMultiplyLargeVectors(Vector3d*** a, Vector3d*** b, int xnumber, int ynumber, int znumber, int lnumber);

#endif