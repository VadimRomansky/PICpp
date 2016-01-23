#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "specialmath.h"

double evaluateError(double** hessenbergMatrix, double* vector, double beta, int n) {
	double* resVector = new double[n + 1];

	for (int i = 0; i < n + 1; ++i) {
		resVector[i] = 0;
		for (int j = 0; j < n; ++j) {
			resVector[i] += hessenbergMatrix[i][j] * vector[j];
		}
		if (i == 0) {
			resVector[i] -= beta;
		}
	}

	double norm = 0;
	for (int i = 0; i < n + 1; ++i) {
		norm += resVector[i] * resVector[i];
	}

	delete[] resVector;

	return sqrt(norm);
}

double** multiplySpecialMatrixVector(std::vector<MatrixElement>** matrix, Vector3d* vector, int xnumber, int lnumber) {
	double** result = new double*[xnumber];
	int i = 0;
#pragma omp parallel for shared(result, matrix, vector, xnumber, lnumber) private(i)
	for (i = 0; i < xnumber; ++i) {
		result[i] = new double[lnumber];
		for (int l = 0; l < lnumber; ++l) {
			result[i][l] = 0;
			for (int m = 0; m < matrix[i][l].size(); ++m) {
				MatrixElement element = matrix[i][l][m];

				result[i][l] += element.value * vector[element.i][element.l];
			}
		}
	}

	return result;
}

double** multiplySpecialMatrixVector(std::vector<MatrixElement>** matrix, double** vector, int xnumber, int lnumber) {
	double** result = new double*[xnumber];
	int i = 0;
#pragma omp parallel for shared(matrix, vector, xnumber, lnumber) private(i)
	for (i = 0; i < xnumber; ++i) {
		result[i] = new double[lnumber];
		for (int l = 0; l < lnumber; ++l) {
			result[i][l] = 0;
			for (int m = 0; m < matrix[i][l].size(); ++m) {
				MatrixElement element = matrix[i][l][m];

				result[i][l] += element.value * vector[element.i][element.l];
			}
		}
	}

	return result;
}

double*** arnoldiIterations(std::vector<MatrixElement>** matrix, double** outHessenbergMatrix, int n, double*** prevBasis, double** prevHessenbergMatrix, int xnumber, int lnumber) {
	double*** resultBasis = new double**[n];
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

	double** tempVector = multiplySpecialMatrixVector(matrix, resultBasis[n - 2], xnumber, lnumber);
	double b = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector, xnumber, lnumber));


	for (int m = 0; m < n - 1; ++m) {
		double a = scalarMultiplyLargeVectors(resultBasis[m], tempVector, xnumber, lnumber);
		outHessenbergMatrix[m][n - 2] = scalarMultiplyLargeVectors(resultBasis[m], tempVector, xnumber, lnumber);
		int i = 0;
#pragma omp parallel for shared(tempVector, outHessenbergMatrix, resultBasis, lnumber, m, n) private(i)
		for (i = 0; i < xnumber; ++i) {
			for (int l = 0; l < lnumber; ++l) {
				tempVector[i][l] -= outHessenbergMatrix[m][n - 2] * resultBasis[m][i][l];
				alertNaNOrInfinity(tempVector[i][l], "tempVector = NaN\n");
			}
		}
	}
	outHessenbergMatrix[n - 1][n - 2] = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector, xnumber, lnumber));
	if (outHessenbergMatrix[n - 1][n - 2] > 0) {
		for (int i = 0; i < xnumber; ++i) {
			for (int l = 0; l < lnumber; ++l) {
				tempVector[i][l] /= outHessenbergMatrix[n - 1][n - 2];
				alertNaNOrInfinity(tempVector[i][l], "tempVector = NaN\n");
			}
		}
	} else {
		printf("outHessenbergMatrix[n-1][n-2] == 0\n");
	}

	resultBasis[n - 1] = tempVector;

	return resultBasis;
}

void generalizedMinimalResidualMethod(std::vector<MatrixElement>** matrix, double** rightPart, double** outvector, int xnumber, int lnumber) {
	printf("start GMRES\n");
	double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, lnumber));

	if (norm == 0) {
		for (int i = 0; i < xnumber; ++i) {
			for (int l = 0; l < lnumber; ++l) {
				outvector[i][l] = 0;
			}
		}
		return;
	}

	//#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		for (int l = 0; l < lnumber; ++l) {
			rightPart[i][l] /= norm;
			for (int m = 0; m < matrix[i][l].size(); ++m) {
				double value = matrix[i][l][m].value;
				//matrix[i][l][m].value /= norm;
				value = matrix[i][l][m].value;
			}
		}
	}

	int matrixDimension = lnumber * xnumber;

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

	double*** basis = new double**[1];
	basis[0] = new double*[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		basis[0][i] = new double[lnumber];
		for (int l = 0; l < lnumber; ++l) {
			basis[0][i][l] = rightPart[i][l];
		}
	}
	double*** newBasis;

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
	double maxRelativeError = maxErrorLevel / (matrixDimension);

	while ((relativeError > maxRelativeError && n < min2(maxGMRESIterations, matrixDimension + 3)) || (n <= 4)) {
		printf("GMRES iteration %d\n", n);
		newHessenbergMatrix = new double*[n];
		for (int i = 0; i < n; ++i) {
			newHessenbergMatrix[i] = new double[n - 1];
		}
		newBasis = arnoldiIterations(matrix, newHessenbergMatrix, n, basis, hessenbergMatrix, xnumber, lnumber);

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
			for (int i = 0; i < n; ++i) {
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
			if (Rmatrix[i][i] > 0) {
				y[i] /= Rmatrix[i][i];
			} else {
				y[i] = 0;
				printf("Rmatrix[%d][%d] = 0\n", i, i);
			}
			alertNaNOrInfinity(y[i], "y = NaN\n");
		}

		error = fabs(beta * Qmatrix[n - 1][0]);
		for (int i = 0; i < xnumber; ++i) {
			for (int l = 0; l < lnumber; ++l) {
				outvector[i][l] = 0;
				for (int m = 0; m < n - 1; ++m) {
					outvector[i][l] += basis[m][i][l] * y[m] * norm;
					//outvector[i][l] += basis[m][i][l] * y[m];
				}
			}
		}

		double normRightPart = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, lnumber));
		relativeError = error / normRightPart;

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

	for (int i = 0; i < xnumber; ++i) {
		for (int l = 0; l < lnumber; ++l) {
			outvector[i][l] = 0;
			for (int m = 0; m < n - 1; ++m) {
				outvector[i][l] += basis[m][i][l] * y[m] * norm;
				//outvector[i][l] += basis[m][i][l] * y[m];
			}
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

	for (int m = 0; m < n; ++m) {
		for (int i = 0; i < xnumber; ++i) {
			delete[] basis[m][i];
		}
		delete[] basis[m];
	}
	delete[] basis;

	delete[] y;

	for (int i = 0; i < xnumber; ++i) {
		for (int l = 0; l < lnumber; ++l) {
			rightPart[i][l] *= norm;
			for (int m = 0; m < matrix[i][l].size(); ++m) {
				double value = matrix[i][l][m].value;
				//matrix[i][l][m].value *= norm;
				value = matrix[i][l][m].value;
			}
		}
	}
}

double scalarMultiplyLargeVectors(double** a, double** b, int xnumber, int lnumber) {
	double result = 0;
	for (int i = 0; i < xnumber; ++i) {
		for (int l = 0; l < lnumber; ++l) {
			result += a[i][l] * b[i][l];
		}
	}
	return result;
}

double scalarMultiplyLargeVectors(Vector3d* a, Vector3d* b, int xnumber, int lnumber) {
	double result = 0;
	for (int i = 0; i < xnumber; ++i) {
		for (int l = 0; l < lnumber; ++l) {
			result += a[i][l] * b[i][l];
		}
	}
	return result;
}

//a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0
double solve4orderEquation(double a4, double a3, double a2, double a1, double a0, double startX) {
	double nextX = startX;
	double prevX = startX;

	int iterationCount = 0;
	while ((fabs(nextX - prevX) > fabs(prevX * 1.0E-20) || iterationCount == 0) && (iterationCount < 10000)) {
		prevX = nextX;
		nextX = prevX - polynom4Value(a4, a3, a2, a1, a0, prevX) / polynom4DerivativeValue(a4, a3, a2, a1, prevX);
		iterationCount++;
	}

	return nextX;
}

double polynom4Value(double a4, double a3, double a2, double a1, double a0, double x) {
	return (((a4 * x + a3) * x + a2) * x + a1) * x + a0;
}

double polynom4DerivativeValue(double a4, double a3, double a2, double a1, double x) {
	return ((4.0 * a4 * x + 3.0 * a3) * x + 2.0 * a2) * x + a1;
}

double solveBigOrderEquation(const double* coefficients, int n, double startX) {
	double nextX = startX;
	double prevX = startX;

	int iterationCount = 0;
	while ((fabs(nextX - prevX) > fabs(prevX * 1.0E-20) || iterationCount == 0) && (iterationCount < 10000)) {
		prevX = nextX;
		nextX = prevX - polynomValue(coefficients, prevX, n) / polynomDerivativeValue(coefficients, prevX, n);
		iterationCount++;
	}

	return nextX;
}

double polynomValue(const double* coefficients, double x, int n) {
	double result = coefficients[0];
	for (int i = 1; i <= n; ++i) {
		result = result * x + coefficients[i];
	}
	return result;
}

double polynomDerivativeValue(const double* coefficients, double x, int n) {
	double result = n * coefficients[0];
	for (int i = 1; i < n; ++i) {
		result = result * x + (n - i) * coefficients[i];
	}
	return result;
}

