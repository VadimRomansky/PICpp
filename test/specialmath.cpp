#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "specialmath.h"
#include "complex.h"

double** arnoldiIterations(double** matrix, double** outHessenbergMatrix, int n, double** prevBasis, double** prevHessenbergMatrix){
	double** resultBasis = new double*[n];
	for (int m = 0; m < n - 1; ++m) {
		resultBasis[m] = prevBasis[m];
	}
	delete[] prevBasis;

	for (int i = 0; i < n; ++i) {
		if (i < n - 1) {
			for (int j = 0; j < n - 2; ++j) {
				outHessenbergMatrix[i][j] = prevHessenbergMatrix[i][j];
			}
			delete[] prevHessenbergMatrix[i];
		} else {
			for (int j = 0; j < n - 2; ++j) {
				outHessenbergMatrix[i][j] = 0;
			}
		}
	}
	delete[] prevHessenbergMatrix;

	double* tempVector = multiplyMatrixVector(matrix, resultBasis[n - 2]);
	double b = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector));


	for (int m = 0; m < n - 1; ++m) {
		double a = scalarMultiplyLargeVectors(resultBasis[m], tempVector);
		outHessenbergMatrix[m][n - 2] = scalarMultiplyLargeVectors(resultBasis[m], tempVector);
		int i = 0;
		#pragma omp parallel for shared(tempVector, outHessenbergMatrix, resultBasis, lnumber, m, n) private(i)
		for (i = 0; i < number; ++i) {
			tempVector[i] -= outHessenbergMatrix[m][n - 2] * resultBasis[m][i];
			//alertNaNOrInfinity(tempVector[i], "tempVector = NaN\n");
		}
	}
	outHessenbergMatrix[n - 1][n - 2] = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector));
	if (outHessenbergMatrix[n - 1][n - 2] > 0) {
		for (int i = 0; i < number ; ++i) {
			tempVector[i] /= outHessenbergMatrix[n - 1][n - 2];
			//alertNaNOrInfinity(tempVector[i], "tempVector = NaN\n");
		}
	} else {
		printf("outHessenbergMatrix[n-1][n-2] == 0\n");
	}

	resultBasis[n - 1] = tempVector;

	return resultBasis;
}

double* generalizedMinimalResidualMethod(double** matrix, double* rightPart){
	printf("start GMRES\n");
	double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart));
	double* outvector = new double[number];

	if(norm == 0) {
		for(int i = 0; i < number; ++i) {
			outvector[i] = 0;
		}
		return outvector;
	}

	//#pragma omp parallel for
	for (int i = 0; i < number; ++i) {
		rightPart[i] /= norm;
		for (int m = 0; m < number; ++m) {
			double value = matrix[i][m];
			//matrix[i][l][m].value /= norm;
			value = matrix[i][m];
		}
				
	}

	int matrixDimension = number;

	double** hessenbergMatrix;
	double** newHessenbergMatrix;
	hessenbergMatrix = new double*[1];
	hessenbergMatrix[0] = new double[1];

	double** Qmatrix = new double*[2];
	double** Rmatrix = new double*[2];
	double** oldQmatrix = new double*[2];
	double** oldRmatrix = new double*[2];

	for (int i = 0; i < 2; ++i) {
		Qmatrix[i] = new double[2];
		oldQmatrix[i] = new double[2];
	}

	Rmatrix[0] = new double[1];
	Rmatrix[1] = new double[1];
	oldRmatrix[0] = new double[1];
	oldRmatrix[1] = new double[1];

	double** basis = new double*[1];
	basis[0] = new double[number];
	for (int i = 0; i < number; ++i) {
		basis[0][i] = rightPart[i];
	}
	double** newBasis;

	int n = 2;
	double beta = 1.0;
	double error = beta;
	double* y = new double[1];

	double rho;
	double sigma;
	double cosn;
	double sinn;
	double module;

	double relativeError = 1;
	double maxRelativeError = 1E-12/(matrixDimension);

	while (relativeError > maxRelativeError  && n < matrixDimension + 2) {
		printf("GMRES iteration %d\n", n);
		newHessenbergMatrix = new double*[n];
		for (int i = 0; i < n; ++i) {
			newHessenbergMatrix[i] = new double[n - 1];
		}
		newBasis = arnoldiIterations(matrix, newHessenbergMatrix, n, basis, hessenbergMatrix);

		hessenbergMatrix = newHessenbergMatrix;
		basis = newBasis;

		if (n == 2) {
			rho = hessenbergMatrix[0][0];
			sigma = hessenbergMatrix[1][0];

			module = sqrt(rho * rho + sigma * sigma);

			cosn = rho / module;
			sinn = sigma / module;

			Qmatrix[0][0] = cosn;
			Qmatrix[0][1] = sinn;
			Qmatrix[1][0] = -sinn;
			Qmatrix[1][1] = cosn;

			oldQmatrix[0][0] = Qmatrix[0][0];
			oldQmatrix[0][1] = Qmatrix[0][1];
			oldQmatrix[1][0] = Qmatrix[1][0];
			oldQmatrix[1][1] = Qmatrix[1][1];

			Rmatrix[0][0] = module;
			Rmatrix[1][0] = 0;

			oldRmatrix[0][0] = Rmatrix[0][0];
			oldRmatrix[1][0] = Rmatrix[1][0];

		} else {
			Rmatrix = new double*[n];
			for (int i = 0; i < n; ++i) {
				Rmatrix[i] = new double[n - 1];
				if (i < n - 1) {
					for (int j = 0; j < n - 2; ++j) {
						Rmatrix[i][j] = oldRmatrix[i][j];
					}
				} else {
					for (int j = 0; j < n - 2; ++j) {
						Rmatrix[i][j] = 0;
					}
				}
			}

			Qmatrix = new double*[n];
			for (int i = 0; i < n;++i) {
				Qmatrix[i] = new double[n];
				if (i < n - 1) {
					for (int j = 0; j < n - 1; ++j) {
						Qmatrix[i][j] = oldQmatrix[i][j];
					}
					Qmatrix[i][n - 1] = 0;
				} else {
					for (int j = 0; j < n - 1; ++j) {
						Qmatrix[i][j] = 0;
					}
					Qmatrix[n - 1][n - 1] = 1;
				}
			}

			for (int i = 0; i < n; ++i) {
				Rmatrix[i][n - 2] = 0;
				for (int j = 0; j < n; ++j) {
					Rmatrix[i][n - 2] += Qmatrix[i][j] * hessenbergMatrix[j][n - 2];
				}
			}
			rho = Rmatrix[n - 2][n - 2];
			sigma = Rmatrix[n - 1][n - 2];

			module = sqrt(rho * rho + sigma * sigma);

			cosn = rho / module;
			sinn = sigma / module;

			Rmatrix[n - 2][n - 2] = module;
			Rmatrix[n - 1][n - 2] = 0;

			for (int j = 0; j < n - 1; ++j) {
				Qmatrix[n - 2][j] = cosn * oldQmatrix[n - 2][j];
				Qmatrix[n - 1][j] = -sinn * oldQmatrix[n - 2][j];
			}
			Qmatrix[n - 2][n - 1] = sinn;
			Qmatrix[n - 1][n - 1] = cosn;
		}

		delete[] y;
		y = new double[n - 1];

		for (int i = n - 2; i >= 0; --i) {
			y[i] = beta * Qmatrix[i][0];
			for (int j = n - 2; j > i; --j) {
				y[i] -= Rmatrix[i][j] * y[j];
			}
			if(Rmatrix[i][i] > 0){
				y[i] /= Rmatrix[i][i];
			} else {
				y[i] = 0;
				printf("Rmatrix[%d][%d] = 0\n", i, i);
			}
			//alertNaNOrInfinity(y[i], "y = NaN\n");
		}

		error = fabs(beta * Qmatrix[n - 1][0]);
		for (int i = 0; i < number; ++i) {
			outvector[i] = 0;
			for (int m = 0; m < n-1; ++m) {
				outvector[i] += basis[m][i] * y[m] *norm;
				//outvector[i][l] += basis[m][i][l] * y[m];
			}				
		}

		double normRightPart = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart));
		relativeError = error/normRightPart;

		for (int i = 0; i < n - 1; ++i) {
			delete[] oldQmatrix[i];
			delete[] oldRmatrix[i];
		}
		delete[] oldQmatrix;
		delete[] oldRmatrix;

		oldQmatrix = Qmatrix;
		oldRmatrix = Rmatrix;

		n++;
	}

	n = n - 1;

	//out result

	for (int i = 0; i < number; ++i) {
		outvector[i] = 0;
		for (int m = 0; m < n-1; ++m) {
			outvector[i] += basis[m][i] * y[m]*norm;
		}
	}

	for (int i = 0; i < n; ++i) {
		delete[] Qmatrix[i];
		delete[] Rmatrix[i];
		delete[] hessenbergMatrix[i];
	}
	delete[] Qmatrix;
	delete[] Rmatrix;
	delete[] hessenbergMatrix;

	for (int m = 0; m < n-1; ++m) {
		delete[] basis[m];
	}
	delete[] basis;

	delete[] y;

	for (int i = 0; i < number; ++i) {
		rightPart[i] *= norm;
		for (int m = 0; m < number; ++m) {
			double value = matrix[i][m];
			//matrix[i][l][m].value *= norm;
			value = matrix[i][m];
		}
	}

	return outvector;
}

double scalarMultiplyLargeVectors(double* a, double* b){
	double result = 0;
	for(int i = 0; i < number; ++i){
		result += a[i]*b[i];
	}
	return result;
}

double* multiplyMatrixVector(double** matrix, double* vector){
	double* result = new double[number];
	for(int i = 0; i < number; ++i){
		result[i] = 0;
		for(int j = 0; j < number; ++j){
			result[i] += matrix[i][j]*vector[j];
		}
	}

	return result;
}

//a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0
double solve4orderEquation(double a4, double a3, double a2, double a1, double a0, double startX){
	double nextX = startX;
	double prevX = startX;

	int iterationCount = 0;
	while((abs(nextX - prevX) > startX*1E-6 || iterationCount == 0) && (iterationCount < 100)){
		prevX = nextX;
		nextX = prevX - polynomValue(a4, a3, a2, a1, a0, prevX)/polynomDerivativeValue(a4, a3, a2, a1, prevX);
		iterationCount++;
	}

	return nextX;
}

double polynomValue(double a4, double a3, double a2, double a1, double a0, double x){
	return (((a4*x + a3)*x + a2)*x + a1)*x + a0;
}

double polynomDerivativeValue(double a4, double a3, double a2, double a1, double x){
	return  ((4*a4*x + 3*a3)*x + 2*a2)*x + a1;
}

Complex*** evaluateFourierTranslation(double*** a, int xnumber, int ynumber, int znumber){
	Complex*** result = new Complex**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		result[i] = new Complex*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			result[i][j] = new Complex[znumber];
			for(int k = 0; k < znumber; ++k){
				result[i][j][k] = Complex(0, 0);

			for(int tempi = 0; tempi < xnumber; ++tempi){
					for(int tempj = 0; tempj < ynumber; ++tempj){
						for(int tempk = 0; tempk < znumber; ++tempk){
							result[i][j][k] += complexExp(-2*pi*((i*tempi*1.0/xnumber) + (j*tempj*1.0/ynumber) + (k*tempk*1.0/znumber)))*a[tempi][tempj][tempk];
						}
					}
				}

				result[i][j][k] = result[i][j][k]/(xnumber*ynumber*znumber);
			}
		}
	}

	return result;
}

double*** evaluateReverceFourierTranslation(Complex*** a, int xnumber, int ynumber, int znumber){
	double*** result = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		result[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			result[i][j] = new double[znumber];
			for(int k = 0; k < znumber; ++k){
				result[i][j][k] = 0;

				for(int tempi = 0; tempi < xnumber; ++tempi){
					for(int tempj = 0; tempj < ynumber; ++tempj){
						for(int tempk = 0; tempk < znumber; ++tempk){
							if(tempi < xnumber/2.0){
								result[i][j][k] += (complexExp(2*pi*((i*tempi*1.0/xnumber) + (j*tempj*1.0/ynumber) + (k*tempk*1.0/znumber)))*a[tempi][tempj][tempk]).re;
							} else {
								result[i][j][k] += (complexExp(-2*pi*((i*(tempi - xnumber + 2)*1.0/xnumber) + (j*tempj*1.0/ynumber) + (k*tempk*1.0/znumber)))*a[tempi][tempj][tempk]).re;
							}
						}
					}
				}
			}
		}
	}

	return result;
}