#ifndef _SPECIAL_MATH_H_
#define _SPECIAL_MATH_H_

#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "rightPartEvaluator.h"

class MatrixElement;
class Vector3d;
class Complex;
class LargeVectorBasis;

void generalizedMinimalResidualMethod(std::vector<MatrixElement> ****matrix, double ****rightPart, double ****outvector,
                                      int xnumberAdded, int ynumberAdded, int znumberAdded, int lnumber, int xnumberGeneral,
                                      int znumberGeneral, int ynumberGeneral, int additionalBinNumber, double precision,
                                      int maxIteration, bool periodicX, bool periodicY, bool periodicZ, int verbocity, double* leftOutGmresBuffer, double* rightOutGmresBuffer, double* leftInGmresBuffer, double* rightInGmresBuffer, double* frontOutGmresBuffer, double* backOutGmresBuffer, double* frontInGmresBuffer, double* backInGmresBuffer, double* bottomOutGmresBuffer, double* topOutGmresBuffer, double* bottomInGmresBuffer, double* topInGmresBuffer, LargeVectorBasis* gmresBasis, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
void arnoldiIterations(std::vector<MatrixElement> ****matrix, double **outHessenbergMatrix, int n,
                       LargeVectorBasis* gmresBasis, double **prevHessenbergMatrix, int xnumberAdded, int ynumberAdded,
                       int znumberAdded, int additionalBinNumber, int lnumber, bool periodicX, bool periodicY, bool periodicZ, int rank, int nprocs, double* leftOutGmresBuffer, double* rightOutGmresBuffer, double* leftInGmresBuffer, double* rightInGmresBuffer, double* frontOutGmresBuffer, double* backOutGmresBuffer, double* frontInGmresBuffer, double* backInGmresBuffer, double* bottomOutGmresBuffer, double* topOutGmresBuffer, double* bottomInGmresBuffer, double* topInGmresBuffer, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, double**** vector, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber, bool periodicX, bool
                                       periodicY, bool periodicZ, int rank, int nprocs, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, Vector3d*** vector, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber, bool periodicX, bool
                                       periodicY, bool periodicZ, int rank, int nprocs, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
void multiplySpecialMatrixVector(double**** result, std::vector<MatrixElement>**** matrix, Vector3d*** vector, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber, bool periodicX, bool
                                 periodicY, bool periodicZ, int rank, int nprocs, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
void multiplySpecialMatrixVector(double**** result, std::vector<MatrixElement>**** matrix, double**** vector, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber, bool periodicX, bool
                                 periodicY, bool periodicZ, int rank, int nprocs, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
double evaluateError(double**** hessenbergMatrix, double*** vector, double beta, int n);
double scalarMultiplyLargeVectors(double ****a, double ****b, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                                  int lnumber, bool periodicX, bool periodicY, bool periodicZ, int rank, int nprocs, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
double scalarMultiplyLargeVectors(Vector3d ***a, Vector3d ***b, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                                  int lnumber, bool periodicX, bool periodicY, bool periodicZ, int rank, int nprocs, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
void simpleIterationSolver(double**** outVector, double**** tempVector, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber, int rank, int nprocs,
                           int xnumberGeneral, int znumberGeneral, int ynumberGeneral, double precision,
                           int maxIteration, bool periodicX, bool periodicY, bool periodicZ, int verbocity, double* leftOutBuffer, double* rightOutBuffer, double* leftInBuffer, double* rightInBuffer, double* frontOutBuffer, double* backOutBuffer, double* frontInBuffer, double* backInBuffer, double* bottomOutBuffer, double* topOutBuffer, double* bottomInBuffer, double* topInBuffer, RightPartIterativeEvaluator* rightPartEvaluator, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
double normDifferenceLargeVectors(double**** a, double**** b, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber, bool periodicX, bool
                                  periodicY, bool periodicZ, int rank, int nprocs, MPI_Comm& cartComm, int* cartCoord, int* cartDim);

void conjugateGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber, double precision, int maxIteration, bool periodicX, bool
                             periodicY, bool periodicZ, int verbosity, MPI_Comm& cartComm, int* cartCoord, int* cartDim);

void biconjugateGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber, double precision, int maxIteration, bool periodicX, bool
                               periodicY, bool periodicZ, int verbosity, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
void biconjugateStabilizedGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, double precision, int maxIteration, bool periodicX, bool
                                         periodicY, bool periodicZ, int verbosity, MPI_Comm& cartComm, int* cartCoord, int* cartDim, double**** residual, double**** fiestResidual, double**** v, double**** p, double**** s, double**** t, double* leftOutGmresBuffer, double* rightOutGmresBuffer, double* leftInGmresBuffer, double* rightInGmresBuffer, double* frontOutGmresBuffer, double* backOutGmresBuffer, double* frontInGmresBuffer, double* backInGmresBuffer, double* bottomOutGmresBuffer, double* topOutGmresBuffer, double* bottomInGmresBuffer, double* topInGmresBuffer, bool& converges);
void transposeSpecialMatrix(std::vector<MatrixElement>**** result, std::vector<MatrixElement>**** matrix, int xnumber, int ynumber, int znumber, int lnumber);

void gaussSeidelMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber,
                       int xnumberGeneral, int znumberGeneral, int ynumberGeneral, double precision,
                       int maxIteration, bool periodicX, bool periodicY, bool periodicZ, bool startFromZero, int verbocity, MPI_Comm& cartComm, int* cartCoord, int* cartDim);

bool indexLower(const MatrixElement& element, int i, int j, int k, int l);
bool indexEqual(const MatrixElement& element, int i, int j, int k, int l);
bool indexUpper(const MatrixElement& element, int i, int j, int k, int l);

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