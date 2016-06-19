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

double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, Vector3d*** vector, int xnumber, int ynumber, int znumber, int lnumber) {
	double**** result = new double***[xnumber];
	int i = 0;
#pragma omp parallel for shared(result, matrix, vector, xnumber, lnumber) private(i)
	for (i = 0; i < xnumber; ++i) {
		result[i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			result[i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				result[i][j][k] = new double[lnumber];
				for (int l = 0; l < lnumber; ++l) {
					result[i][j][k][l] = 0;
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];

						result[i][j][k][l] += element.value * vector[element.i][element.j][element.k][element.l];
					}
				}
			}
		}
	}

	return result;
}

double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, double**** vector, int xnumber, int ynumber, int znumber, int lnumber) {
	double**** result = new double***[xnumber];
	int i = 0;
#pragma omp parallel for shared(result, matrix, vector, xnumber, lnumber) private(i)
	for (i = 0; i < xnumber; ++i) {
		result[i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			result[i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				result[i][j][k] = new double[lnumber];
				for (int l = 0; l < lnumber; ++l) {
					result[i][j][k][l] = 0;
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];

						result[i][j][k][l] += element.value * vector[element.i][element.j][element.k][element.l];
					}
				}
			}
		}
	}

	return result;
}

void multiplySpecialMatrixVector(double**** result, std::vector<MatrixElement>**** matrix, Vector3d*** vector, int xnumber, int ynumber, int znumber, int lnumber) {
	int i = 0;
#pragma omp parallel for shared(result, matrix, vector, xnumber, lnumber) private(i)
	for (i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result[i][j][k][l] = 0;
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];

						result[i][j][k][l] += element.value * vector[element.i][element.j][element.k][element.l];
					}
				}
			}
		}
	}
}

void multiplySpecialMatrixVector(double**** result, std::vector<MatrixElement>**** matrix, double**** vector, int xnumber, int ynumber, int znumber, int lnumber) {
	int i = 0;
#pragma omp parallel for shared(result, matrix, vector, xnumber, lnumber) private(i)
	for (i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) { ;
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result[i][j][k][l] = 0;
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];

						result[i][j][k][l] += element.value * vector[element.i][element.j][element.k][element.l];
					}
				}
			}
		}
	}
}

double***** arnoldiIterations(std::vector<MatrixElement>**** matrix, double** outHessenbergMatrix, int n, double***** prevBasis, double** prevHessenbergMatrix, int xnumber, int ynumber, int znumber, int lnumber) {
	double***** resultBasis = new double****[n];
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

	double**** tempVector = multiplySpecialMatrixVector(matrix, resultBasis[n - 2], xnumber, ynumber, znumber, lnumber);
	/*for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				for(int l = 0; l < lnumber; ++l){
					printf("%g\n", tempVector[i][j][k][l]);
				}
			}
		}
	}*/


	for (int m = 0; m < n - 1; ++m) {
		double a = scalarMultiplyLargeVectors(resultBasis[m], tempVector, xnumber, ynumber, znumber, lnumber);
		outHessenbergMatrix[m][n - 2] = scalarMultiplyLargeVectors(resultBasis[m], tempVector, xnumber, ynumber, znumber, lnumber);
		//printf("outHessenbergMatrix[%d][%d] = %g\n", m, n-2, outHessenbergMatrix[m][n - 2]);
		int i = 0;
#pragma omp parallel for shared(tempVector, outHessenbergMatrix, resultBasis, lnumber, ynumber, znumber, m, n) private(i)
		for (i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][j][k][l] -= outHessenbergMatrix[m][n - 2] * resultBasis[m][i][j][k][l];
						alertNaNOrInfinity(tempVector[i][j][k][l], "tempVector = NaN\n");
						//printf("tempvector[%d][%d][%d][%d] = %g\n", i, j, k, l, tempVector[i][j][k][l]);
					}
				}
			}
		}
	}
	outHessenbergMatrix[n - 1][n - 2] = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector, xnumber, ynumber, znumber, lnumber));
	//printf("outHessenbergMatrix[%d][%d] = %g\n", n-1, n-2, outHessenbergMatrix[n - 1][n - 2]);
	if (outHessenbergMatrix[n - 1][n - 2] > 0) {
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][j][k][l] /= outHessenbergMatrix[n - 1][n - 2];
						alertNaNOrInfinity(tempVector[i][j][k][l], "tempVector = NaN\n");
						//printf("tempvector[%d][%d][%d][%d] = %g\n", i, j, k, l, tempVector[i][j][k][l]);
					}
				}
			}
		}
	} else {
		printf("outHessenbergMatrix[n-1][n-2] == 0\n");
	}

	resultBasis[n - 1] = tempVector;

	return resultBasis;
}

void generalizedMinimalResidualMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outvector, int xnumber, int ynumber, int znumber, int lnumber, double precision, int maxIteration) {
	printf("start GMRES\n");
	double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, ynumber, znumber, lnumber));
	//printf("norm = %g\n", norm);

	if (norm == 0) {
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outvector[i][j][k][l] = 0;
					}
				}
			}
		}
		return;
	}

	//#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					rightPart[i][j][k][l] /= norm;
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						double value = matrix[i][j][k][l][m].value;
						//matrix[i][l][m].value /= norm;
						value = matrix[i][j][k][l][m].value;
					}
				}
			}
		}
	}

	int matrixDimension = lnumber * xnumber * ynumber * znumber;

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

	double***** basis = new double****[1];
	basis[0] = new double***[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		basis[0][i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			basis[0][i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				basis[0][i][j][k] = new double[lnumber];
				for (int l = 0; l < lnumber; ++l) {
					basis[0][i][j][k][l] = rightPart[i][j][k][l];
				}
			}
		}
	}
	double***** newBasis;

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
	//double maxRelativeError = maxErrorLevel / (matrixDimension);
	double maxRelativeError = precision / (matrixDimension);

	while (((relativeError > maxRelativeError && n < min2(maxIteration, matrixDimension + 3)) || (n <= 4))) {
		printf("GMRES iteration %d\n", n);
		newHessenbergMatrix = new double*[n];
		for (int i = 0; i < n; ++i) {
			newHessenbergMatrix[i] = new double[n - 1];
		}
		newBasis = arnoldiIterations(matrix, newHessenbergMatrix, n, basis, hessenbergMatrix, xnumber, ynumber, znumber, lnumber);

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
		//printf("n = %d\n", n);

		for (int i = n - 2; i >= 0; --i) {
			y[i] = beta * Qmatrix[i][0];
			//printf("y[%d] = %g\n", i, y[i]);
			for (int j = n - 2; j > i; --j) {
				y[i] -= Rmatrix[i][j] * y[j];
			}
			if (Rmatrix[i][i] > 0) {
				y[i] /= Rmatrix[i][i];
			} else {
				y[i] = 0;
				printf("Rmatrix[%d][%d] = 0\n", i, i);
			}
			//printf("y[%d] = %g\n", i, y[i]);
			alertNaNOrInfinity(y[i], "y = NaN\n");
		}

		error = fabs(beta * Qmatrix[n - 1][0]);
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outvector[i][j][k][l] = 0;
						for (int m = 0; m < n - 1; ++m) {
							outvector[i][j][k][l] += basis[m][i][j][k][l] * y[m] * norm;
							//printf("outvector[%d][%d][%d][%d] %d = %g\n", i, j, k, l, m, outvector[i][j][k][l]);
							//printf("norm = %g\n", norm);
							//printf("y[%d] = %g\n", m, y[m]);
							//outvector[i][l] += basis[m][i][l] * y[m];
						}
						//printf("outvector[%d][%d][%d][%d] = %g\n", i, j, k, l, outvector[i][j][k][l]);
					}
				}
			}
		}

		double normRightPart = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, ynumber, znumber, lnumber));
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
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					outvector[i][j][k][l] = 0;
					for (int m = 0; m < n - 1; ++m) {
						outvector[i][j][k][l] += basis[m][i][j][k][l] * y[m] * norm;
						//outvector[i][l] += basis[m][i][l] * y[m];
					}
				}
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
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					delete[] basis[m][i][j][k];
				}
				delete[] basis[m][i][j];
			}
			delete[] basis[m][i];
		}
		delete[] basis[m];
	}
	delete[] basis;

	delete[] y;

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					rightPart[i][j][k][l] *= norm;
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						double value = matrix[i][j][k][l][m].value;
						//matrix[i][l][m].value *= norm;
						value = matrix[i][j][k][l][m].value;
					}
				}
			}
		}
	}
}

double scalarMultiplyLargeVectors(double**** a, double**** b, int xnumber, int ynumber, int znumber, int lnumber) {
	double result = 0;
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result += a[i][j][k][l] * b[i][j][k][l];
				}
			}
		}
	}
	return result;
}

double scalarMultiplyLargeVectors(Vector3d*** a, Vector3d*** b, int xnumber, int ynumber, int znumber, int lnumber) {
	double result = 0;
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result += a[i][j][k][l] * b[i][j][k][l];
				}
			}
		}
	}
	return result;
}

void transposeSpecialMatrix(std::vector<MatrixElement>**** result, std::vector<MatrixElement>**** matrix, int xnumber, int ynumber, int znumber, int lnumber) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j =  0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				for(int l = 0; l < lnumber; ++l){
					for(int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];

						result[element.i][element.j][element.k][element.l].push_back(MatrixElement(element.value, i, j, k, l));
					}
				}
			}
		}
	}
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

Complex*** fastFourierTransition(double*** a, int xnumber, int ynumber, int znumber){
	int kx = xnumber;
	while(kx > 1) {
		if(kx%2 != 0) {
			printf("xnumber is not 2^N\n");
		}
		kx = kx/2;
	}
	int ky = ynumber;
	while(ky > 1) {
		if(ky%2 != 0) {
			printf("ynumber is not 2^N\n");
		}
		ky = ky/2;
	}
	int kz = znumber;
	while(kx > 1) {
		if(kx%2 != 0) {
			printf("znumber is not 2^N\n");
		}
		kz = kz/2;
	}

	Complex*** result = new Complex**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		result[i] = new Complex*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			result[i][j] = new Complex[znumber];
			for(int k = 0; k < znumber; ++k){
				result[i][j][k] = Complex(a[i][j][k], 0);
			}
		}
	}
	sortInputFastFourierX(result, xnumber, ynumber, znumber);

	Complex* tempResult;
	if(xnumber > 1) {
		tempResult = new Complex[xnumber];
		for (int ycount = 0; ycount < ynumber; ++ycount) {
			for (int zcount = 0; zcount < znumber; ++zcount) {
				for (int i = 0; i < xnumber / 2; ++i) {
					result[i * 2][ycount][zcount] = result[i * 2][ycount][zcount] + result[i * 2 + 1][ycount][zcount];
					result[i * 2 + 1][ycount][zcount] = result[i * 2][ycount][zcount] - result[i * 2 + 1][ycount][zcount];
				}

				int k = 4;
				while (k <= xnumber) {
					int l = xnumber / k;
					for (int i = 0; i < l; ++i) {
						for (int m = 0; m < k / 2; ++m) {
							tempResult[i * k + m] = result[i * k + m][ycount][zcount] + result[i * k + (k / 2) + m][ycount][zcount] * complexExp(
								-2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[i * k + m][ycount][zcount] - result[i * k + (k / 2) + m][ycount][zcount] * complexExp(
								-2 * pi * m / k);
						}
					}

					for (int i = 0; i < xnumber; ++i) {
						result[i][ycount][zcount] = tempResult[i];
					}
					k = k * 2;
				}
			}
		}
		delete[] tempResult;
	}

	if(ynumber > 1) {
		sortInputFastFourierY(result, xnumber, ynumber, znumber);
		tempResult = new Complex[ynumber];
		for (int xcount = 0; xcount < xnumber; ++xcount) {
			for (int zcount = 0; zcount < znumber; ++zcount) {
				for (int i = 0; i < ynumber / 2; ++i) {
					result[xcount][i * 2][zcount] = result[xcount][i * 2][zcount] + result[xcount][i * 2 + 1][zcount];
					result[xcount][i * 2 + 1][zcount] = result[xcount][i * 2][zcount] - result[xcount][i * 2 + 1][zcount];
				}

				int k = 4;
				while (k <= ynumber) {
					int l = ynumber / k;
					for (int i = 0; i < l; ++i) {
						for (int m = 0; m < k / 2; ++m) {
							tempResult[i * k + m] = result[xcount][i * k + m][zcount] + result[xcount][i * k + (k / 2) + m][zcount] * complexExp(
								-2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[xcount][i * k + m][zcount] - result[xcount][i * k + (k / 2) + m][zcount] * complexExp(
								-2 * pi * m / k);
						}
					}

					for (int i = 0; i < ynumber; ++i) {
						result[xcount][i][zcount] = tempResult[i];
					}
					k = k * 2;
				}
			}
		}
		delete[] tempResult;
	}

	if(znumber > 1) {
		sortInputFastFourierY(result, xnumber, ynumber, znumber);
		tempResult = new Complex[znumber];
		for (int xcount = 0; xcount < xnumber; ++xcount) {
			for (int ycount = 0; ycount < ynumber; ++ycount) {
				for (int i = 0; i < znumber / 2; ++i) {
					result[xcount][ycount][i * 2] = result[xcount][ycount][i * 2] + result[xcount][ycount][i * 2 + 1];
					result[xcount][ycount][i * 2 + 1] = result[xcount][ycount][i * 2] - result[xcount][ycount][i * 2 + 1];
				}

				int k = 4;
				while (k <= znumber) {
					int l = znumber / k;
					for (int i = 0; i < l; ++i) {
						for (int m = 0; m < k / 2; ++m) {
							tempResult[i * k + m] = result[xcount][ycount][i * k + m] + result[xcount][ycount][i * k + (k / 2) + m] * complexExp(
								-2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[xcount][ycount][i * k + m] - result[xcount][ycount][i * k + (k / 2) + m] * complexExp(
								-2 * pi * m / k);
						}
					}

					for (int i = 0; i < znumber; ++i) {
						result[xcount][ycount][i] = tempResult[i];
					}
					k = k * 2;
				}
			}
		}
		delete[] tempResult;
	}


	return result;
}

void sortInputFastFourierX(double*** a, int xnumber, int ynumber, int znumber){
	double* tempA = new double[xnumber];
	for(int ycount = 0; ycount < ynumber; ++ycount) {
		for (int zcount = 0; zcount < znumber; ++zcount) {
			int k = xnumber;

			while (k > 2) {
				int m = xnumber / k;

				for (int i = 0; i < m; ++i) {
					for (int j = 0; j < k / 2; ++j) {
						tempA[i * k + j] = a[i * k + 2 * j][ycount][zcount];
						tempA[i * k + (k / 2) + j] = a[i * k + 2 * j + 1][ycount][zcount];
					}
				}

				for (int i = 0; i < xnumber; ++i) {
					a[i][ycount][zcount] = tempA[i];
				}

				k = k / 2;
			}
		}
	}
	delete[] tempA;
}

void sortInputFastFourierY(double*** a, int xnumber, int ynumber, int znumber){
	double* tempA = new double[ynumber];
	for(int xcount = 0; xcount < xnumber; ++xcount) {
		for (int zcount = 0; zcount < znumber; ++zcount) {
			int k = ynumber;

			while (k > 2) {
				int m = ynumber / k;

				for (int i = 0; i < m; ++i) {
					for (int j = 0; j < k / 2; ++j) {
						tempA[i * k + j] = a[xcount][i * k + 2 * j][zcount];
						tempA[i * k + (k / 2) + j] = a[xcount][i * k + 2 * j + 1][zcount];
					}
				}

				for (int i = 0; i < ynumber; ++i) {
					a[xcount][i][zcount] = tempA[i];
				}

				k = k / 2;
			}
		}
	}
	delete[] tempA;
}

void sortInputFastFourierZ(double*** a, int xnumber, int ynumber, int znumber){
	double* tempA = new double[znumber];
	for(int xcount = 0; xcount < xnumber; ++xcount) {
		for (int ycount = 0; ycount < ynumber; ++ycount) {
			int k = znumber;

			while (k > 2) {
				int m = znumber / k;

				for (int i = 0; i < m; ++i) {
					for (int j = 0; j < k / 2; ++j) {
						tempA[i * k + j] = a[xcount][ycount][i * k + 2 * j];
						tempA[i * k + (k / 2) + j] = a[xcount][ycount][i * k + 2 * j + 1];
					}
				}

				for (int i = 0; i < znumber; ++i) {
					a[xcount][ycount][i] = tempA[i];
				}

				k = k / 2;
			}
		}
	}
	delete[] tempA;
}

double*** fastFourierReverceTransition(Complex ***a, int xnumber, int ynumber, int znumber) {
	int kx = xnumber;
	while(kx > 1) {
		if(kx%2 != 0) {
			printf("xnumber is not 2^N\n");
		}
		kx = kx/2;
	}
	int ky = ynumber;
	while(ky > 1) {
		if(ky%2 != 0) {
			printf("ynumber is not 2^N\n");
		}
		ky = ky/2;
	}
	int kz = znumber;
	while(kx > 1) {
		if(kx%2 != 0) {
			printf("znumber is not 2^N\n");
		}
		kz = kz/2;
	}

	Complex*** result = new Complex**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		result[i] = new Complex*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			result[i][j] = new Complex[znumber];
			for(int k = 0; k < znumber; ++k){
				result[i][j][k] = a[i][j][k];
			}
		}
	}
	sortInputFastFourierX(result, xnumber, ynumber, znumber);

	Complex* tempResult;
	if(xnumber > 1) {
		tempResult = new Complex[xnumber];
		for (int ycount = 0; ycount < ynumber; ++ycount) {
			for (int zcount = 0; zcount < znumber; ++zcount) {
				for (int i = 0; i < xnumber / 2; ++i) {
					result[i * 2][ycount][zcount] = result[i * 2][ycount][zcount] + result[i * 2 + 1][ycount][zcount];
					result[i * 2 + 1][ycount][zcount] = result[i * 2][ycount][zcount] - result[i * 2 + 1][ycount][zcount];
				}

				int k = 4;
				while (k <= xnumber) {
					int l = xnumber / k;
					for (int i = 0; i < l; ++i) {
						for (int m = 0; m < k / 2; ++m) {
							tempResult[i * k + m] = result[i * k + m][ycount][zcount] + result[i * k + (k / 2) + m][ycount][zcount] * complexExp(
								2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[i * k + m][ycount][zcount] - result[i * k + (k / 2) + m][ycount][zcount] * complexExp(
								2 * pi * m / k);
						}
					}

					for (int i = 0; i < xnumber; ++i) {
						result[i][ycount][zcount] = tempResult[i];
					}
					k = k * 2;
				}
			}
		}
		delete[] tempResult;
	}

	if(ynumber > 1) {
		sortInputFastFourierY(result, xnumber, ynumber, znumber);
		tempResult = new Complex[ynumber];
		for (int xcount = 0; xcount < xnumber; ++xcount) {
			for (int zcount = 0; zcount < znumber; ++zcount) {
				for (int i = 0; i < ynumber / 2; ++i) {
					result[xcount][i * 2][zcount] = result[xcount][i * 2][zcount] + result[xcount][i * 2 + 1][zcount];
					result[xcount][i * 2 + 1][zcount] = result[xcount][i * 2][zcount] - result[xcount][i * 2 + 1][zcount];
				}

				int k = 4;
				while (k <= ynumber) {
					int l = ynumber / k;
					for (int i = 0; i < l; ++i) {
						for (int m = 0; m < k / 2; ++m) {
							tempResult[i * k + m] = result[xcount][i * k + m][zcount] + result[xcount][i * k + (k / 2) + m][zcount] * complexExp(
								2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[xcount][i * k + m][zcount] - result[xcount][i * k + (k / 2) + m][zcount] * complexExp(
								2 * pi * m / k);
						}
					}

					for (int i = 0; i < ynumber; ++i) {
						result[xcount][i][zcount] = tempResult[i];
					}
					k = k * 2;
				}
			}
		}
		delete[] tempResult;
	}

	if(znumber > 1) {
		sortInputFastFourierY(result, xnumber, ynumber, znumber);
		tempResult = new Complex[znumber];
		for (int xcount = 0; xcount < xnumber; ++xcount) {
			for (int ycount = 0; ycount < ynumber; ++ycount) {
				for (int i = 0; i < znumber / 2; ++i) {
					result[xcount][ycount][i * 2] = result[xcount][ycount][i * 2] + result[xcount][ycount][i * 2 + 1];
					result[xcount][ycount][i * 2 + 1] = result[xcount][ycount][i * 2] - result[xcount][ycount][i * 2 + 1];
				}

				int k = 4;
				while (k <= znumber) {
					int l = znumber / k;
					for (int i = 0; i < l; ++i) {
						for (int m = 0; m < k / 2; ++m) {
							tempResult[i * k + m] = result[xcount][ycount][i * k + m] + result[xcount][ycount][i * k + (k / 2) + m] * complexExp(
								2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[xcount][ycount][i * k + m] - result[xcount][ycount][i * k + (k / 2) + m] * complexExp(
								2 * pi * m / k);
						}
					}

					for (int i = 0; i < znumber; ++i) {
						result[xcount][ycount][i] = tempResult[i];
					}
					k = k * 2;
				}
			}
		}
		delete[] tempResult;
	}

	double*** realPart = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		realPart[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			realPart[i][j] = new double[znumber];
			for(int k = 0; k < znumber; ++k){
				realPart[i][j][k] = result[i][j][k].re/(xnumber*ynumber*znumber);
			}
			delete[] result[i][j];
		}
		delete[] result[i];
	}

	delete[] result;

	return realPart;
}

void sortInputFastFourierX(Complex ***a, int xnumber, int ynumber, int znumber) {
	Complex* tempA = new Complex[xnumber];
	for(int ycount = 0; ycount < ynumber; ++ycount) {
		for (int zcount = 0; zcount < znumber; ++zcount) {
			int k = xnumber;

			while (k > 2) {
				int m = xnumber / k;

				for (int i = 0; i < m; ++i) {
					for (int j = 0; j < k / 2; ++j) {
						tempA[i * k + j] = a[i * k + 2 * j][ycount][zcount];
						tempA[i * k + (k / 2) + j] = a[i * k + 2 * j + 1][ycount][zcount];
					}
				}

				for (int i = 0; i < xnumber; ++i) {
					a[i][ycount][zcount] = tempA[i];
				}

				k = k / 2;
			}
		}
	}
	delete[] tempA;
}

void sortInputFastFourierY(Complex ***a, int xnumber, int ynumber, int znumber) {
	Complex* tempA = new Complex[ynumber];
	for(int xcount = 0; xcount < xnumber; ++xcount) {
		for (int zcount = 0; zcount < znumber; ++zcount) {
			int k = ynumber;

			while (k > 2) {
				int m = ynumber / k;

				for (int i = 0; i < m; ++i) {
					for (int j = 0; j < k / 2; ++j) {
						tempA[i * k + j] = a[xcount][i * k + 2 * j][zcount];
						tempA[i * k + (k / 2) + j] = a[xcount][i * k + 2 * j + 1][zcount];
					}
				}

				for (int i = 0; i < ynumber; ++i) {
					a[xcount][i][zcount] = tempA[i];
				}

				k = k / 2;
			}
		}
	}
	delete[] tempA;
}

void sortInputFastFourierZ(Complex ***a, int xnumber, int ynumber, int znumber) {
	Complex* tempA = new Complex[znumber];
	for(int xcount = 0; xcount < xnumber; ++xcount) {
		for (int ycount = 0; ycount < ynumber; ++ycount) {
			int k = znumber;

			while (k > 2) {
				int m = znumber / k;

				for (int i = 0; i < m; ++i) {
					for (int j = 0; j < k / 2; ++j) {
						tempA[i * k + j] = a[xcount][ycount][i * k + 2 * j];
						tempA[i * k + (k / 2) + j] = a[xcount][ycount][i * k + 2 * j + 1];
					}
				}

				for (int i = 0; i < znumber; ++i) {
					a[xcount][ycount][i] = tempA[i];
				}

				k = k / 2;
			}
		}
	}
	delete[] tempA;
}

void conjugateGradientMethod(std::vector<MatrixElement> ****matrix, double ****rightPart, double ****outVector,
							 int xnumber, int ynumber, int znumber, int lnumber, double precision, int maxIteration) {

	printf("start conjugate gradient\n");
	double**** residual = new double***[xnumber];
	double**** prevResidual = new double***[xnumber];
	double**** z = new double***[xnumber];
	double**** tempVector = new double***[xnumber];

	for(int i = 0; i < xnumber; ++i){
		residual[i] = new double**[ynumber];
		prevResidual[i] = new double**[ynumber];
		z[i] = new double**[ynumber];
		tempVector[i] = new double**[ynumber];
		for(int j = 0; j < ynumber; ++j){
			residual[i][j] = new double*[znumber];
			prevResidual[i][j] = new double*[znumber];
			z[i][j] = new double*[znumber];
			tempVector[i][j] = new double*[znumber];
			for(int k = 0; k < znumber; ++k){
				residual[i][j][k] = new double[lnumber];
				prevResidual[i][j][k] = new double[lnumber];
				z[i][j][k] = new double[lnumber];
				tempVector[i][j][k] = new double[lnumber];
				for(int l = 0; l < lnumber; ++l){
					outVector[i][j][k][l] = 0;
					prevResidual[i][j][k][l] = rightPart[i][j][k][l];
					z[i][j][k][l] = rightPart[i][j][k][l];
					tempVector[i][j][k][l] = 0;
				}
			}
		}
	}


	/*for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];
						prevResidual[i][j][k][l] += element.value * outVector[element.i][element.j][element.k][element.l];
					}
				}
			}
		}
	}*/

	int iteration = 0;


    double prevResidualNorm2 = scalarMultiplyLargeVectors(prevResidual, prevResidual, xnumber, ynumber, znumber, lnumber);
    double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, ynumber, znumber, lnumber);

	double relativeError = sqrt(residualNorm2/rightPartNorm2);

	while((iteration < maxIteration) && (iteration < xnumber*ynumber*znumber*lnumber) && (relativeError > (precision/(xnumber*ynumber*znumber*lnumber)))){
		printf("conjugate gradient iteration %d\n", iteration);

		multiplySpecialMatrixVector(tempVector, matrix, z, xnumber, ynumber, znumber, lnumber);

		double alpha = prevResidualNorm2/scalarMultiplyLargeVectors(tempVector, z, xnumber, ynumber, znumber, lnumber);

		for(int i = 0; i < xnumber; ++i){
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					for(int l = 0; l < lnumber; ++l){
						outVector[i][j][k][l] += alpha*z[i][j][k][l];
						residual[i][j][k][l] = prevResidual[i][j][k][l] - alpha*tempVector[i][j][k][l];
					}
				}
			}
		}

        residualNorm2 = scalarMultiplyLargeVectors(residual, residual, xnumber, ynumber, znumber, lnumber);

		double beta  = residualNorm2/prevResidualNorm2;

		for(int i = 0; i < xnumber; ++i){
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					for(int l = 0; l < lnumber; ++l){
						z[i][j][k][l] = residual[i][j][k][l] + beta*z[i][j][k][l];
					}
				}
			}
		}

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2/rightPartNorm2);
		iteration++;
	}

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				delete[] residual[i][j][k];
				delete[] prevResidual[i][j][k];
				delete[] z[i][j][k];
				delete[] tempVector[i][j][k];
			}
			delete[] residual[i][j];
			delete[] prevResidual[i][j];
			delete[] z[i][j];
			delete[] tempVector[i][j];
		}
		delete[] residual[i];
		delete[] prevResidual[i];
		delete[] z[i];
		delete[] tempVector[i];
	}
	delete[] residual;
	delete[] prevResidual;
	delete[] z;
	delete[] tempVector;
}

void biconjugateGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int xnumber, int ynumber, int znumber, int lnumber, double precision, int maxIteration) {
	printf("start biconjugate gradient\n");
	double**** residual = new double***[xnumber];
	double**** prevResidual = new double***[xnumber];
	double**** z = new double***[xnumber];
	double**** p = new double***[xnumber];
	double**** s = new double***[xnumber];
	double**** tempVector = new double***[xnumber];
	double**** tempVector2 = new double***[xnumber];
	std::vector<MatrixElement>**** transposedMatrix = new std::vector<MatrixElement>***[xnumber];

	for(int i = 0; i < xnumber; ++i){
		residual[i] = new double**[ynumber];
		prevResidual[i] = new double**[ynumber];
		z[i] = new double**[ynumber];
		p[i] = new double**[ynumber];
		s[i] = new double**[ynumber];
		tempVector[i] = new double**[ynumber];
		tempVector2[i] = new double**[ynumber];
		transposedMatrix[i] = new std::vector<MatrixElement>**[ynumber];
		for(int j = 0; j < ynumber; ++j){
			residual[i][j] = new double*[znumber];
			prevResidual[i][j] = new double*[znumber];
			z[i][j] = new double*[znumber];
			p[i][j] = new double*[znumber];
			s[i][j] = new double*[znumber];
			tempVector[i][j] = new double*[znumber];
			tempVector2[i][j] = new double*[znumber];
			transposedMatrix[i][j] = new std::vector<MatrixElement>*[znumber];
			for(int k = 0; k < znumber; ++k){
				residual[i][j][k] = new double[lnumber];
				prevResidual[i][j][k] = new double[lnumber];
				z[i][j][k] = new double[lnumber];
				p[i][j][k] = new double[lnumber];
				s[i][j][k] = new double[lnumber];
				tempVector[i][j][k] = new double[lnumber];
				tempVector2[i][j][k] = new double[lnumber];
				transposedMatrix[i][j][k] = new std::vector<MatrixElement>[lnumber];
				for(int l = 0; l < lnumber; ++l){
					outVector[i][j][k][l] = 0;
					prevResidual[i][j][k][l] = rightPart[i][j][k][l];
					z[i][j][k][l] = rightPart[i][j][k][l];
					p[i][j][k][l] = rightPart[i][j][k][l];
					s[i][j][k][l] = rightPart[i][j][k][l];
					tempVector[i][j][k][l] = 0;
					tempVector2[i][j][k][l] = 0;
				}
			}
		}
	}

	transposeSpecialMatrix(transposedMatrix, matrix, xnumber, ynumber, znumber, lnumber);



	int iteration = 0;


    double prevResidualNorm2 = scalarMultiplyLargeVectors(p, prevResidual, xnumber, ynumber, znumber, lnumber);
    double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, ynumber, znumber, lnumber);

	double relativeError = sqrt(residualNorm2/rightPartNorm2);

	while((iteration < maxIteration) && (iteration < xnumber*ynumber*znumber*lnumber) && (relativeError > (precision/(xnumber*ynumber*znumber*lnumber)))){
		printf("biconjugate gradient iteration %d\n", iteration);

		multiplySpecialMatrixVector(tempVector, matrix, z, xnumber, ynumber, znumber, lnumber);
		multiplySpecialMatrixVector(tempVector2, transposedMatrix, s, xnumber, ynumber, znumber, lnumber);

		double alpha = prevResidualNorm2/scalarMultiplyLargeVectors(tempVector, s, xnumber, ynumber, znumber, lnumber);

		for(int i = 0; i < xnumber; ++i){
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					for(int l = 0; l < lnumber; ++l){
						outVector[i][j][k][l] += alpha*z[i][j][k][l];
						residual[i][j][k][l] = prevResidual[i][j][k][l] - alpha*tempVector[i][j][k][l];
						p[i][j][k][l] = p[i][j][k][l] - alpha*tempVector2[i][j][k][l];
					}
				}
			}
		}

        residualNorm2 = scalarMultiplyLargeVectors(p, residual, xnumber, ynumber, znumber, lnumber);

		double beta  = residualNorm2/prevResidualNorm2;

		for(int i = 0; i < xnumber; ++i){
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					for(int l = 0; l < lnumber; ++l){
						z[i][j][k][l] = residual[i][j][k][l] + beta*z[i][j][k][l];
						s[i][j][k][l] = p[i][j][k][l] + beta*s[i][j][k][l];
					}
				}
			}
		}

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2/rightPartNorm2);
		iteration++;
	}

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				delete[] residual[i][j][k];
				delete[] prevResidual[i][j][k];
				delete[] z[i][j][k];
				delete[] p[i][j][k];
				delete[] s[i][j][k];
				delete[] tempVector[i][j][k];
				delete[] tempVector2[i][j][k];
				delete[] transposedMatrix[i][j][k];
			}
			delete[] residual[i][j];
			delete[] prevResidual[i][j];
			delete[] z[i][j];
			delete[] p[i][j];
			delete[] s[i][j];
			delete[] tempVector[i][j];
			delete[] tempVector2[i][j];
			delete[] transposedMatrix[i][j];
		}
		delete[] residual[i];
		delete[] prevResidual[i];
		delete[] z[i];
		delete[] p[i];
		delete[] s[i];
		delete[] tempVector[i];
		delete[] tempVector2[i];
		delete[] transposedMatrix[i];
	}
	delete[] residual;
	delete[] prevResidual;
	delete[] z;
	delete[] p;
	delete[] s;
	delete[] tempVector;
	delete[] tempVector2;
	delete[] transposedMatrix;
}

void biconjugateStabilizedGradientMethod(std::vector<MatrixElement> ****matrix, double ****rightPart,
										 double ****outVector, int xnumber, int ynumber, int znumber, int lnumber,
										 double precision, int maxIteration) {
	printf("start biconjugate gradient\n");
	double**** residual = new double***[xnumber];
	double**** firstResidual = new double***[xnumber];
	double**** p = new double***[xnumber];
	double**** v = new double***[xnumber];
	double**** s = new double***[xnumber];
	double**** t = new double***[xnumber];

	double alpha = 1;
	double rho = 1;
	double omega = 1;

	for(int i = 0; i < xnumber; ++i){
		residual[i] = new double**[ynumber];
		firstResidual[i] = new double**[ynumber];
		v[i] = new double**[ynumber];
		p[i] = new double**[ynumber];
		s[i] = new double**[ynumber];
		t[i] = new double**[ynumber];
		for(int j = 0; j < ynumber; ++j){
			residual[i][j] = new double*[znumber];
			firstResidual[i][j] = new double*[znumber];
			v[i][j] = new double*[znumber];
			p[i][j] = new double*[znumber];
			s[i][j] = new double*[znumber];
			t[i][j] = new double*[znumber];
			for(int k = 0; k < znumber; ++k){
				residual[i][j][k] = new double[lnumber];
				firstResidual[i][j][k] = new double[lnumber];
				v[i][j][k] = new double[lnumber];
				p[i][j][k] = new double[lnumber];
				s[i][j][k] = new double[lnumber];
				t[i][j][k] = new double[lnumber];
				for(int l = 0; l < lnumber; ++l){
					outVector[i][j][k][l] = 0;
					firstResidual[i][j][k][l] = rightPart[i][j][k][l];
					residual[i][j][k][l] = rightPart[i][j][k][l];
					v[i][j][k][l] = 0;
					p[i][j][k][l] = 0;
					s[i][j][k][l] = 0;
					t[i][j][k][l] = 0;
				}
			}
		}
	}



	int iteration = 0;


	double prevResidualNorm2 = scalarMultiplyLargeVectors(residual, residual, xnumber, ynumber, znumber, lnumber);
	double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, ynumber, znumber, lnumber);

	double relativeError = sqrt(residualNorm2/rightPartNorm2);

	while((iteration < maxIteration) && (iteration < xnumber*ynumber*znumber*lnumber) && (relativeError > (precision/(xnumber*ynumber*znumber*lnumber)))){
		printf("biconjugate gradient iteration %d\n", iteration);

		double newRho = scalarMultiplyLargeVectors(firstResidual, residual, xnumber, ynumber, znumber, lnumber);

		double beta = (newRho/rho)*(alpha/omega);
		rho = newRho;

		for(int i = 0; i < xnumber; ++i){
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					for(int l = 0; l < lnumber; ++l){
						p[i][j][k][l] = residual[i][j][k][l] + beta*(p[i][j][k][l] - omega*v[i][j][k][l]);
					}
				}
			}
		}
		multiplySpecialMatrixVector(v, matrix, p, xnumber, ynumber, znumber, lnumber);

		alpha = rho/scalarMultiplyLargeVectors(firstResidual, v, xnumber, ynumber, znumber, lnumber);

		for(int i = 0; i < xnumber; ++i){
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					for(int l = 0; l < lnumber; ++l){
						s[i][j][k][l] = residual[i][j][k][l] - alpha*v[i][j][k][l];
					}
				}
			}
		}

		multiplySpecialMatrixVector(t, matrix, s, xnumber, ynumber, znumber, lnumber);

		omega = scalarMultiplyLargeVectors(t, s, xnumber, ynumber, znumber, lnumber)/scalarMultiplyLargeVectors(t, t, xnumber, ynumber, znumber, lnumber);

		for(int i = 0; i < xnumber; ++i){
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					for(int l = 0; l < lnumber; ++l){
						outVector[i][j][k][l] = outVector[i][j][k][l] + omega*s[i][j][k][l] + alpha*p[i][j][k][l];
						residual[i][j][k][l] = s[i][j][k][l] - omega*t[i][j][k][l];
					}
				}
			}
		}

		residualNorm2 = scalarMultiplyLargeVectors(residual, residual, xnumber, ynumber, znumber, lnumber);

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2/rightPartNorm2);
		iteration++;
	}

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				delete[] residual[i][j][k];
				delete[] firstResidual[i][j][k];
				delete[] v[i][j][k];
				delete[] p[i][j][k];
				delete[] s[i][j][k];
				delete[] t[i][j][k];
			}
			delete[] residual[i][j];
			delete[] firstResidual[i][j];
			delete[] v[i][j];
			delete[] p[i][j];
			delete[] s[i][j];
			delete[] t[i][j];
		}
		delete[] residual[i];
		delete[] firstResidual[i];
		delete[] v[i];
		delete[] p[i];
		delete[] s[i];
		delete[] t[i];
	}
	delete[] residual;
	delete[] firstResidual;
	delete[] v;
	delete[] p;
	delete[] s;
	delete[] t;
}




