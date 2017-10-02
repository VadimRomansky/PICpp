#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "specialmath.h"
#include "complex.h"
#include "mpi_util.h"


double** arnoldiIterations(double** matrix, double** outHessenbergMatrix, int n, double** prevBasis, double** prevHessenbergMatrix, int number, bool periodic, int rank, int nprocs){
	MPI_Barrier(MPI_COMM_WORLD);
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
	//printf("update hessenberg\n");
	double* tempVector = multiplyMatrixVector(matrix, resultBasis[n - 2], number, periodic, rank, nprocs);
	//printf("mult special matrix");
	double* bufferRightSend = new double[1];
	double* bufferRightRecv = new double[1];
	double* bufferLeftSend = new double[1];
	double* bufferLeftRecv = new double[1];

	//printf("start exchange\n");

	MPI_Barrier(MPI_COMM_WORLD);

	if( (periodic) || (rank > 0)) sendLargeVectorToLeft(tempVector, bufferLeftSend, number);
	if( (periodic) || (rank < nprocs-1)) sendLargeVectorToRight(tempVector, bufferRightSend, number);

	MPI_Barrier(MPI_COMM_WORLD);

	//printf("finish sending\n");

	if( (periodic) || (rank < nprocs-1)) receiveLargeVectorFromRight(tempVector, bufferRightRecv, number);
	if( (periodic) || (rank > 0)) receiveLargeVectorFromLeft(tempVector, bufferLeftRecv, number);


	MPI_Barrier(MPI_COMM_WORLD);

	/*if(rank == 0){
		printf("tempVector = %15.10g\n", tempVector[0][0][0][1]);
		printf("tempVector = %15.10g\n", tempVector[1][0][0][1]);
		printf("tempVector left = %15.10g\n", tempVector[xnumber - 1][0][0][1]);
		printf("tempVector left = %15.10g\n", tempVector[xnumber][0][0][1]);
	} else {

	}*/
	//printf("finish exchange\n");

	delete[] bufferRightSend;
	delete[] bufferRightRecv;
	delete[] bufferLeftSend;
	delete[] bufferLeftRecv;


	for (int m = 0; m < n - 1; ++m) {
		//double a = scalarMultiplyLargeVectors(resultBasis[m], tempVector, xnumber, ynumber, znumber, lnumber);
		outHessenbergMatrix[m][n - 2] = scalarMultiplyLargeVectors(resultBasis[m], tempVector, number, periodic, rank, nprocs);
		//printf("outHessenbergMatrix[%d][%d] = %g\n", m, n-2, outHessenbergMatrix[m][n - 2]);

		//for (int i = 0; i < xnumber+1; ++i) {
        if(rank == 0 && !periodic){
            tempVector[0] -= outHessenbergMatrix[m][n - 2] * resultBasis[m][0];
        }
		for (int i = 1; i < number - 1; ++i) {
			tempVector[i] -= outHessenbergMatrix[m][n - 2] * resultBasis[m][i];
		}
        if(rank == nprocs - 1 && !periodic){
            tempVector[number - 1] -= outHessenbergMatrix[m][n - 2] * resultBasis[m][number - 1];
        }
	}
	//printf("finish orthogonalisation\n");
	outHessenbergMatrix[n - 1][n - 2] = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector, number, periodic, rank, nprocs));
	//printf("outHessenbergMatrix[%d][%d] = %g\n", n-1, n-2, outHessenbergMatrix[n - 1][n - 2]);
	if (outHessenbergMatrix[n - 1][n - 2] > 0) {
		//for (int i = 0; i < xnumber+1; ++i) {
		for (int i = 0; i < number; ++i) {
						tempVector[i] /= outHessenbergMatrix[n - 1][n - 2];
		}
	} else {
		printf("outHessenbergMatrix[n-1][n-2] == 0\n");
	}
	//printf("finish normalization\n");

	resultBasis[n - 1] = tempVector;
	//printf("arnolsi end\n");

	return resultBasis;
}

double* generalizedMinimalResidualMethod(double** matrix, double* rightPart, int number, bool periodic, int rank, int nprocs){
	printf("start GMRES\n");
	double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, number, periodic, rank, nprocs));
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

	int matrixDimension = (number-2)*nprocs+2;

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
	double maxRelativeError = 1E-16/(matrixDimension);

	while (relativeError > maxRelativeError  && n < matrixDimension + 2) {
		printf("GMRES iteration %d\n", n);
		newHessenbergMatrix = new double*[n];
		for (int i = 0; i < n; ++i) {
			newHessenbergMatrix[i] = new double[n - 1];
		}
		newBasis = arnoldiIterations(matrix, newHessenbergMatrix, n, basis, hessenbergMatrix, number, periodic, rank, nprocs);

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

		double normRightPart = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, number, periodic, rank, nprocs));
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

double scalarMultiplyLargeVectors(double* a, double* b, int number,bool periodic, int rank, int nprocs){
	double* result = new double[1];
	result[0] = 0;
	if((rank == 0) && (!periodic)){
		result[0] += a[0] * b[0];
	}
	for (int i = 1; i < number; ++i) {
	    result[0] += a[i] * b[i];
	}
	if((rank == nprocs - 1) && (!periodic)){
		result[0] += a[number] * b[number];
	}
	if(nprocs > 1) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank != 0) {
			MPI_Send(result, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD);
			//printf("send double rank = %d\n", rank);
		} else {
			for (int i = 1; i < nprocs; ++i) {
				double *temp = new double[1];
				MPI_Status status;
				MPI_Recv(temp, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
				result[0] += temp[0];
				delete[] temp;
				//printf("recv double rnk = 0 from %d\n", i);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0) {
			for (int i = 1; i < nprocs; ++i) {
				MPI_Send(result, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD);
				//printf("send double rank = 0 to %d\n", i);
			}
		} else {
			double *temp = new double[1];
			MPI_Status status;
			MPI_Recv(temp, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
			//printf("recv double rank = %d\n", rank);
			result[0] = temp[0];
			delete[] temp;
		}
	}
	double res = result[0];
	delete[] result;
	MPI_Barrier(MPI_COMM_WORLD);
	return res;
}

double* multiplyMatrixVector(double** matrix, double* vector, int number, bool periodic, int rank, int nprocs){
	double* result = new double[number];
	int i = 0;

	for (i = 0; i < number; ++i) {
		result[i] = 0;
		if((i != 0 && i != number) || ((i == 0) && (rank == 0) &&(!periodic)) || ((i == number) && (rank == nprocs - 1)&&(!periodic))) {
			for (int m = 0; m < number; ++m) {
				result[i] += matrix[i][m] * vector[m];
			}
		}
	}

	return result;
}

void multiplyMatrixVector(double* result, double** matrix, double* vector, int number, bool periodic, int rank, int nprocs){
	int i = 0;

	for (i = 0; i < number; ++i) {
		result[i] = 0;
		if((i != 0 && i != number) || ((i == 0) && (rank == 0) &&(!periodic)) || ((i == number) && (rank == nprocs - 1)&&(!periodic))) {
			for (int m = 0; m < number; ++m) {
				result[i] += matrix[i][m] * vector[m];
			}
		}
	}
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

Complex* evaluateFourierTranslation(double* a, int number){
	Complex* result = new Complex[number];
	for(int i = 0; i < number; ++i){
		result[i] = Complex(0,0);
		for(int tempi = 0; tempi < number; ++tempi){
			result[i] += complexExp(-2*pi*(i*tempi*1.0/number))*a[tempi];
		}
		//result[i] = result[i]/number;
	}

	return result;
}

double* evaluateReverceFourierTranslation(Complex* a, int number){
	double* result = new double[number];
	for(int i = 0; i < number; ++i){
		result[i] = 0;
		for(int tempi = 0; tempi < number; ++tempi){
			if(tempi < number/2.0){
				result[i] += (complexExp(2*pi*((i*tempi*1.0/number)))*a[tempi]).re;
			} else {
				result[i] += (complexExp(-2*pi*((i*(tempi - number + 2)*1.0/number) ))*a[tempi]).re;
			}
		}
		result[i] = result[i]/number;
	}

	return result;
}

Complex* fastFourierTransition(double* a, int n){
	int k = n;
	while(k > 1) {
		if(k%2 != 0) {
			printf("n is not 2^N\n");
		}
		k = k/2;
	}
	if(n == 1) {
		Complex* result = new Complex[1];
		result[0] = Complex(a[0], 0);
		return result;
	}
	if(n == 2) {
		Complex* result = new Complex[2];
		result[0] = Complex(a[0] + a[1], 0);
		result[1] = Complex(a[0] - a[1], 0);
		return result;
	}
	double* tempA = new double[n];
	for(int i = 0; i < n; ++i){
		tempA[i] = a[i];
	}
	sortInputFastFourier(tempA, n);
	Complex* result = new Complex[n];
	Complex* tempResult = new Complex[n];

	for(int i = 0; i < n/2; ++i){
		result[i*2] = Complex(tempA[i*2] + tempA[i*2 + 1], 0);
		result[i*2 + 1] = Complex(tempA[i*2] - tempA[i*2 + 1], 0);
	}

	k = 4;
	while(k <= n){
		int l = n/k;
		for(int i = 0; i < l; ++i){
			for(int m = 0; m < k/2; ++m){
				tempResult[i*k + m] = result[i*k + m] + result[i*k + (k/2) + m]*complexExp(-2*pi*m/k);
				tempResult[i*k + (k/2) + m] = result[i*k + m] - result[i*k + (k/2) + m]*complexExp(-2*pi*m/k);
			}
		}

		for(int i = 0; i < n; ++i){
			result[i] = tempResult[i];
		}
		k = k*2;
	}

	delete[] tempA;
	delete[] tempResult;

	return result;
}

void sortInputFastFourier(double* a, int n){
	double* tempA = new double[n];

	int k = n;

	while (k > 2){
		int m = n/k;

		for(int i = 0; i < m; ++i){
			for(int j = 0; j < k/2; ++j){
				tempA[i*k + j] = a[i*k + 2*j];
				tempA[i*k + (k/2) + j] = a[i*k + 2*j + 1];
			}
		}

		for(int i = 0; i < n; ++i){
			a[i] = tempA[i];
		}

		k = k/2;
	}

	delete[] tempA;
}

double *fastFourierReverceTransition(Complex *a, int n) {
	int k = n;
	while(k > 1) {
		if(k%2 != 0) {
			printf("n is not 2^N\n");
		}
		k = k/2;
	}
	if(n == 1) {
		double* result = new double[1];
		result[0] = a[0].re;
		return result;
	}
	if(n == 2) {
		double* result = new double[2];
		result[0] = (a[0] + a[1]).re/n;
		result[1] = (a[0] - a[1]).re/n;
		return result;
	}
	/*for(int i = n/2; i < n; ++i){
		a[i] = Complex(0,0);
	}*/
	sortInputFastFourierReverce(a, n);
	Complex* result = new Complex[n];
	Complex* tempResult = new Complex[n];

	for(int i = 0; i < n/2; ++i){
		result[i*2] = a[i*2] + a[i*2 + 1];
		result[i*2 + 1] = a[i*2] - a[i*2 + 1];
	}

	k = 4;
	while(k <= n){
		int l = n/k;
		for(int i = 0; i < l; ++i){
			for(int m = 0; m < k/2; ++m){
				tempResult[i*k + m] = result[i*k + m] + result[i*k + (k/2) + m]*complexExp(2*pi*m/k);
				tempResult[i*k + (k/2) + m] = result[i*k + m] - result[i*k + (k/2) + m]*complexExp(2*pi*m/k);
			}
		}

		for(int i = 0; i < n; ++i){
			result[i] = tempResult[i];
		}
		k = k*2;
	}

	double* realPart = new double[n];
	for(int i = 0; i < n; ++i){
		realPart[i] = result[i].re/n;
	}

	delete[] tempResult;
	delete[] result;

	return realPart;
}

void sortInputFastFourierReverce(Complex *a, int n) {
	Complex* tempA = new Complex[n];

	int k = n;

	while (k > 2){
		int m = n/k;

		for(int i = 0; i < m; ++i){
			for(int j = 0; j < k/2; ++j){
				tempA[i*k + j] = a[i*k + 2*j];
				tempA[i*k + (k/2) + j] = a[i*k + 2*j + 1];
			}
		}

		for(int i = 0; i < n; ++i){
			a[i] = tempA[i];
		}

		k = k/2;
	}

	delete[] tempA;
}



void conjugateGradientMethod(double**matrix, double *rightPart, double *outVector, int number, double precision, int maxIteration, bool periodic, int rank, int nprocs) {

	printf("start conjugate gradient\n");
	double* residual = new double[number];
	double* prevResidual = new double[number];
	double* z = new double[number];
	double* tempVector = new double[number];

	for(int i = 0; i < number; ++i) {
		outVector[i] = 0;
		prevResidual[i] = rightPart[i];
		z[i] = rightPart[i];
		tempVector[i] = 0;
	}

	int iteration = 0;



    double prevResidualNorm2 = scalarMultiplyLargeVectors(prevResidual, prevResidual, number, periodic, rank, nprocs);
    double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, number, periodic, rank, nprocs);

	double relativeError = sqrt(residualNorm2/rightPartNorm2);

	while((iteration < maxIteration) && (iteration < number) && (relativeError > (precision/number))){
		printf("conjugate gradient iteration %d\n", iteration);

		multiplyMatrixVector(tempVector, matrix, z, number, periodic, rank, nprocs);

		double alpha = prevResidualNorm2/scalarMultiplyLargeVectors(tempVector, z, number, periodic, rank, nprocs);

		for(int i = 0; i < number; ++i){
						outVector[i] += alpha*z[i];
						residual[i] = prevResidual[i] - alpha*tempVector[i];
		}

        residualNorm2 = scalarMultiplyLargeVectors(residual, residual, number, periodic, rank, nprocs);

		double beta  = residualNorm2/prevResidualNorm2;

		for(int i = 0; i < number; ++i){
						z[i] = residual[i] + beta*z[i];
		}

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2/rightPartNorm2);
		iteration++;
	}

	delete[] residual;
	delete[] prevResidual;
	delete[] z;
	delete[] tempVector;
}

void biconjugateGradientMethod(double** matrix, double* rightPart, double* outVector, int number, double precision, int maxIteration, bool periodic, int rank, int nprocs) {
	printf("start biconjugate gradient\n");
	double* residual = new double[number];
	double* prevResidual = new double[number];
	double* z = new double[number];
	double* p = new double[number];
	double* s = new double[number];
	double* tempVector = new double[number];
	double* tempVector2 = new double[number];
	double** transposedMatrix = new double*[number];

	for(int i = 0; i < number; ++i){
		transposedMatrix[i] = new double[number];
		outVector[i] = 0;
		prevResidual[i] = rightPart[i];
		z[i]= rightPart[i];
		p[i] = rightPart[i];
		s[i] = rightPart[i];
		tempVector[i] = 0;
		tempVector2[i] = 0;
	}

	transposeMatrix(transposedMatrix, matrix, number);



	int iteration = 0;


    double prevResidualNorm2 = scalarMultiplyLargeVectors(p, prevResidual, number, periodic, rank, nprocs);
    double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, number, periodic, rank, nprocs);

	double relativeError = sqrt(residualNorm2/rightPartNorm2);

	while((iteration < maxIteration) && (iteration < number) && (relativeError > (precision/number))){
		printf("biconjugate gradient iteration %d\n", iteration);

		multiplyMatrixVector(tempVector, matrix, z, number, periodic, rank, nprocs);
		multiplyMatrixVector(tempVector2, transposedMatrix, s, number, periodic, rank, nprocs);

		double alpha = prevResidualNorm2/scalarMultiplyLargeVectors(tempVector, s, number, periodic, rank, nprocs);

		for(int i = 0; i < number; ++i){
						outVector[i] += alpha*z[i];
						residual[i] = prevResidual[i] - alpha*tempVector[i];
						p[i] = p[i] - alpha*tempVector2[i];
		}

        residualNorm2 = scalarMultiplyLargeVectors(p, residual, number, periodic, rank, nprocs);

		double beta  = residualNorm2/prevResidualNorm2;

		for(int i = 0; i < number; ++i){
						z[i] = residual[i] + beta*z[i];
						s[i] = p[i] + beta*s[i];
		}

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(scalarMultiplyLargeVectors(residual, residual, number, periodic, rank, nprocs)/rightPartNorm2);
		iteration++;
	}

	for(int i = 0; i < number; ++i){
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

void transposeMatrix(double** result, double** matrix, int number) {
	for(int i = 0; i < number; ++i) {
		for(int j = 0; j < number; ++j) {
			result[i][j] = matrix[j][i];
		}
	}
}

void biconjugateStabilizedGradientMethod(double **matrix, double *rightPart, double *outVector, int number,
										 double precision, int maxIteration, bool periodic, int rank, int nprocs) {
	printf("start biconjugate gradient\n");
	double* residual = new double[number];
	double* firstResidual = new double[number];
	double* p = new double[number];
	double* v = new double[number];
	double* s = new double[number];
	double* t = new double[number];

	double alpha = 1;
	double rho = 1;
	double omega = 1;

	for(int i = 0; i < number; ++i){
					outVector[i] = 0;
					firstResidual[i] = rightPart[i];
					residual[i] = rightPart[i];
					v[i] = 0;
					p[i] = 0;
					s[i] = 0;
					t[i] = 0;
	}



	int iteration = 0;


	double prevResidualNorm2 = scalarMultiplyLargeVectors(residual, residual, number, periodic, rank, nprocs);
	double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, number, periodic, rank, nprocs);

	double relativeError = sqrt(residualNorm2/rightPartNorm2);

	while((iteration < maxIteration) && (iteration < number+1) && (relativeError > (precision/(number)))){
		printf("biconjugate gradient iteration %d\n", iteration);

		double newRho = scalarMultiplyLargeVectors(firstResidual, residual, number, periodic, rank, nprocs);

		double beta = (newRho/rho)*(alpha/omega);
		if(fabs(rho)<1E-200){
			if(fabs(newRho)<1E-200){
				beta = alpha/omega;
				//beta = 0;
			} else {
				printf("denominator = 0\n");
			}
		}
		rho = newRho;

		for(int i = 0; i < number; ++i){
						p[i] = residual[i] + beta*(p[i]- omega*v[i]);
		}
		multiplyMatrixVector(v, matrix, p, number, periodic, rank, nprocs);

		double denominator = scalarMultiplyLargeVectors(firstResidual, v, number, periodic, rank, nprocs);
		if(fabs(denominator)  < 1E-200){
			if(fabs(rho)<1E-200){
				alpha = 1;
			} else {
				alpha = 1;
				printf("denominator = 0\n");
			}
		} else {
			alpha = rho/denominator;
		}


		for(int i = 0; i < number; ++i){
						s[i] = residual[i] - alpha*v[i];
		}

		multiplyMatrixVector(t, matrix, s, number, periodic, rank, nprocs);

		omega = scalarMultiplyLargeVectors(t, s, number, periodic, rank, nprocs)/scalarMultiplyLargeVectors(t, t, number, periodic, rank, nprocs);

		for(int i = 0; i < number; ++i){
						outVector[i] = outVector[i] + omega*s[i] + alpha*p[i];
						residual[i] = s[i] - omega*t[i];

		}

		residualNorm2 = scalarMultiplyLargeVectors(residual, residual, number, periodic, rank, nprocs);

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2/rightPartNorm2);
		iteration++;
	}

	delete[] residual;
	delete[] firstResidual;
	delete[] v;
	delete[] p;
	delete[] s;
	delete[] t;
}


