#ifndef _SPECIAL_MATH_H_
#define _SPECIAL_MATH_H_

#include "stdlib.h"
#include "stdio.h"
#include "vector"

class MatrixElement;
class LargeVectorBasis;

void generalizedMinimalResidualMethod(std::vector<MatrixElement> **** matrix, double ****rightPart, double ****outvector, int Nx, int Ny, int Nz, int Nmomentum, double precision, int maxIteration, int verbosity, LargeVectorBasis* gmresBasis);
void generalizedMinimalResidualMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outvector, double**** initialVector, int Nx, int Ny, int Nz, int Nmomentum, double precision, int maxIteration, int verbosity, LargeVectorBasis* gmresBasis);
void arnoldiIterations(std::vector<MatrixElement> ****matrix, double **outHessenbergMatrix, int n, LargeVectorBasis* gmresBasis, double **prevHessenbergMatrix, int Nx, int Ny, int Nz, int Nmomentum);
double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, double**** vector, int Nx, int Ny, int Nz, int Nmomentum);

void multiplySpecialMatrixVector(double**** result, std::vector<MatrixElement>**** matrix, double**** vector, int Nx, int Ny, int Nz, int Nmomentum);
double evaluateError(double** hessenbergMatrix, double* vector, double beta, int n);
double scalarMultiplyLargeVectors(double ****a, double ****b, int Nx, int Ny, int Nz, int Nmomentum);

void biconjugateStabilizedGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int Nx, int Ny, int Nz, int Nmomentum, double precision, int maxIteration, int verbosity, bool& converges);
void biconjugateStabilizedGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int Nx, int Ny, int Nz, int Nmomentum, double precision, int maxIteration, int verbosity, double**** residual, double**** firstResidual, double**** v, double**** p, double**** s, double**** t, bool& converges);

//void simpleIterationSolver(double* outVector, double* tempVector, int N, double precision, int maxIteration);
//double normDifferenceLargeVectors(double* a, double* b, int N);

//void conjugateGradientMethod(std::vector<MatrixElement>* matrix, double* rightPart, double* outVector, int N, double precision, int maxIteration, int verbosity);

//void biconjugateGradientMethod(std::vector<MatrixElement>* matrix, double* rightPart, double* outVector, int N, double precision, int maxIteration, int verbosity);
//void biconjugateStabilizedGradientMethod(std::vector<MatrixElement>* matrix, double* rightPart, double* outVector, int N, double precision, int maxIteration, int verbosity, double* residual, double* fiestResidual, double* v, double* p, double* s, double* t, bool& converges);
//void transposeSpecialMatrix(std::vector<MatrixElement>* result, std::vector<MatrixElement>* matrix, int N);

//void gaussSeidelMethod(std::vector<MatrixElement>* matrix, double* rightPart, double* outVector, int N, double precision, int maxIteration, int verbosity);


double solve4orderEquation(double a4, double a3, double a2, double a1, double a0, double startX);
double polynom4Value(double a4, double a3, double a2, double a1, double a0, double x);
double polynom4DerivativeValue(double a4, double a3, double a2, double a1, double x);

double solveBigOrderEquation(const double* coefficients, int n, double startX);
double polynomValue(const double* coefficients, double x, int n);
double polynomDerivativeValue(const double* coefficients, double x, int n);

double dichotomySolver(const double* functionValue, int minIndex, int maxIndex, const double* xValue, double y);
double juttnerFunction(const double& u, const double& theta);

//void sortInputFastFourierX(Complex ***a, int xnumber, int ynumber, int znumber);
//void sortInputFastFourierY(Complex ***a, int xnumber, int ynumber, int znumber);
//void sortInputFastFourierZ(Complex ***a, int xnumber, int ynumber, int znumber);
//void sortInputFastFourierX(double ***a, int xnumber, int ynumber, int znumber);
//void sortInputFastFourierY(double ***a, int xnumber, int ynumber, int znumber);
//void sortInputFastFourierZ(double ***a, int xnumber, int ynumber, int znumber);
//double ***fastFourierReverceTransition(Complex ***a, int xnumber, int ynumber, int znumber);
//Complex ***fastFourierTransition(double ***a, int xnumber, int ynumber, int znumber);

void sequentialThreeDiagonalSolverX(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum);
void sequentialThreeDiagonalSolverY(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum);
void sequentialThreeDiagonalSolverZ(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum);

#endif