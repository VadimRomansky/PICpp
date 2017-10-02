#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "specialmath.h"
#include "complex.h"
#include "mpi_util.h"

double** arnoldiIterations(double** matrix, double** outHessenbergMatrix, int n, double** prevBasis, double** prevHessenbergMatrix, int number){
	int rank;
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

	double* tempVector = multiplyMatrixVector(matrix, resultBasis[n - 2], number);



	double bufferRightSend[1];
	double bufferRightRecv[1];
	double bufferLeftSend[1];
	double bufferLeftRecv[1];

	//printf("start exchange\n");


	if(rank > 0) sendGMRESTempVectorToLeft(tempVector, bufferLeftSend, number);
	if(rank < size-1) sendGMRESTempVectorToRight(tempVector, bufferRightSend, number);


	//printf("finish sending\n");
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank < size -1 )receiveGMRESTempVectorFromRight(tempVector, bufferRightRecv, number);
	if(rank > 0) receiveGMRESTempVectorFromLeft(tempVector, bufferLeftRecv, number);



	MPI_Barrier(MPI_COMM_WORLD);


	for (int m = 0; m < n - 1; ++m) {
		outHessenbergMatrix[m][n - 2] = scalarMultiplyLargeVectors(resultBasis[m], tempVector, number);
		if(rank == 0) printf("norm = %g\n", outHessenbergMatrix[m][n - 2]);
		int i = 0;
		for (i = 0; i < number; ++i) {
			tempVector[i] -= outHessenbergMatrix[m][n - 2] * resultBasis[m][i];
			//alertNaNOrInfinity(tempVector[i], "tempVector = NaN\n");
		}
	}
	outHessenbergMatrix[n - 1][n - 2] = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector, number));
	if (outHessenbergMatrix[n - 1][n - 2] > 0) {
		for (int i = 0; i < number ; ++i) {
			tempVector[i] /= outHessenbergMatrix[n - 1][n - 2];
			//alertNaNOrInfinity(tempVector[i], "tempVector = NaN\n");
		}
	} else {
		printf("outHessenbergMatrix[n-1][n-2] == 0\n");
	}


	resultBasis[n - 1] = tempVector;


	printf("end Arnoldi\n");
	return resultBasis;
}

double* generalizedMinimalResidualMethod(double** matrix, double* rightPart, int number){
	printf("start GMRES\n");
	int rank;
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, number));
	if(rank == 0) printf("norm = %g\n", norm);
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
		newBasis = arnoldiIterations(matrix, newHessenbergMatrix, n, basis, hessenbergMatrix, number);
		//printf("end arnoldi\n");
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

		for(int j = 0; j < size; ++j) {
			MPI_Barrier(MPI_COMM_WORLD);
			if (j == rank) {
				printf("rank = %d\n", rank);
				for (int i = 0; i < number; ++i) {
					printf("%15.10lf\n", outvector[i]);
				}
			}
		}

		double normRightPart = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, number));
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

	double bufferRightSend[1];
	double bufferRightRecv[1];
	double bufferLeftSend[1];
	double bufferLeftRecv[1];

	printf("start exchange\n");

	MPI_Barrier(MPI_COMM_WORLD);


	if(rank > 0) sendGMRESTempVectorToLeft(outvector, bufferLeftSend, number);
	if(rank < size-1) sendGMRESTempVectorToRight(outvector, bufferRightSend, number);


	//printf("finish sending\n");
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank < size -1 )receiveGMRESTempVectorFromRight(outvector, bufferRightRecv, number);
	if(rank > 0) receiveGMRESTempVectorFromLeft(outvector, bufferLeftRecv, number);



	MPI_Barrier(MPI_COMM_WORLD);

	//MPI_Barrier(MPI_COMM_WORLD);

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
	printf("enf GMRES\n");
	return outvector;
}

double scalarMultiplyLargeVectors(double* a, double* b, int number){
	double result[1];
	int rank;
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	result[0] = 0;
	if(rank == 0){
		result[0] += a[0]*b[0];
	}
	for(int i = 1; i < number - 1; ++i){
		result[0] += a[i]*b[i];
	}
	if(rank == nprocs - 1){
		result[0] += a[number-1]*b[number-1];
	}

	for(int j = 0; j < nprocs; ++j) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (j == rank) {
			printf("rank = %d\n", rank);
			printf("localScalarMult = %g\n", result[0]);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank != 0){
		MPI_Send(result, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD);
		//printf("send double rank = %d\n", rank);
	} else {
		for(int i = 1; i < nprocs; ++i){
			double temp[1];
			MPI_Status status;
			MPI_Recv(temp, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
			result[0] += temp[0];
			//printf("recv double rnk = 0 from %d\n", i);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0) {
		for (int i = 1; i < nprocs; ++i){
			MPI_Send(result, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD);
			//printf("send double rank = 0 to %d\n", i);
		}
	} else {
		double temp[1];
		MPI_Status status;
		MPI_Recv(temp, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
		//printf("recv double rank = %d\n", rank);
		result[0] = temp[0];
	}
	double res = result[0];
	//MPI_Barrier(MPI_COMM_WORLD);
	return res;

}

double* multiplyMatrixVector(double** matrix, double* vector, int number){
	double* result = new double[number];
	for(int i = 0; i < number; ++i){
		result[i] = 0;
		for(int j = 0; j < number; ++j){
			result[i] += matrix[i][j]*vector[j];
		}
	}

	return result;
}

void multiplyMatrixVector(double* result, double** matrix, double* vector, int number){
	for(int i = 0; i < number; ++i){
		result[i] = 0;
		for(int j = 0; j < number; ++j){
			result[i] += matrix[i][j]*vector[j];
		}
	}
}
