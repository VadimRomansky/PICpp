#include <stdio.h>
#include <math.h>
//#include <mpi.h>
//#include <crtdbg.h>

//#include "memory_debug.h"
//#include "complex.h"
//#include "../util.h"
#include "matrixElement.h"
#include "specialmath.h"
#include "largeVectorBasis.h"
#include "util.h"

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

double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, double**** vector, int Nx, int Ny, int Nz, int Nmomentum) {
	double**** result = create4dArray(Nx, Ny, Nz, Nmomentum, 0.0);
	int i;
#pragma omp parallel for private(i) shared(matrix, vector, result, Nx, Ny, Nz, Nmomentum)
	for(i = 0; i < Nx; ++i){
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
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

void multiplySpecialMatrixVector(double**** result, std::vector < MatrixElement >**** matrix, double**** vector, int Nx, int Ny, int Nz, int Nmomentum) {
	int i;
#pragma omp parallel for private(i) shared(matrix, vector, result, Nx, Ny, Nz, Nmomentum)
	for (i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
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

void arnoldiIterations(std::vector<MatrixElement>**** matrix, double** outHessenbergMatrix, int n, LargeVectorBasis* gmresBasis, double** prevHessenbergMatrix, int Nx, int Ny, int Nz, int Nmomentum) {

	//MPI_Barrier(cartComm);

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
	if (n >= gmresBasis->capacity) {
		gmresBasis->resize(2 * n);
	}
	multiplySpecialMatrixVector(gmresBasis->array[n - 1], matrix, gmresBasis->array[n - 2], Nx, Ny, Nz, Nmomentum);
	gmresBasis->size += 1;
	//printf("mult special matrix");

	//printf("start exchange\n");

	//MPI_Barrier(cartComm);

	//exchangeLargeVector(gmresBasis->array[n-1], xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX, periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutGmresBuffer, rightOutGmresBuffer, leftInGmresBuffer, rightInGmresBuffer, frontOutGmresBuffer, backOutGmresBuffer, frontInGmresBuffer, backInGmresBuffer, bottomOutGmresBuffer, topOutGmresBuffer, bottomInGmresBuffer, topInGmresBuffer);


	for (int m = 0; m < n - 1; ++m) {
		//double a = scalarMultiplyLargeVectors(resultBasis[m], tempVector, xnumber, ynumber, znumber, lnumber);
		outHessenbergMatrix[m][n - 2] = scalarMultiplyLargeVectors(gmresBasis->array[m], gmresBasis->array[n - 1], Nx, Ny, Nz, Nmomentum);
		//printf("outHessenbergMatrix[%d][%d] = %g\n", m, n-2, outHessenbergMatrix[m][n - 2]);

		//for (int i = 0; i < xnumber+1; ++i) {
		int i;
#pragma omp parallel for private(i) shared(gmresBasis, outHessenbergMatrix, m, n, Nx, Ny, Nz, Nmomentum)
		for (i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				for (int k = 0; k < Nz; ++k) {
					for (int l = 0; l < Nmomentum; ++l) {
						gmresBasis->array[n - 1][i][j][k][l] -= outHessenbergMatrix[m][n - 2] * gmresBasis->array[m][i][j][k][l];
					}
				}
			}
		}
	}
	//printf("finish orthogonalisation\n");
	outHessenbergMatrix[n - 1][n - 2] = sqrt(scalarMultiplyLargeVectors(gmresBasis->array[n - 1], gmresBasis->array[n - 1], Nx, Ny, Nz, Nmomentum));
	//printf("outHessenbergMatrix[%d][%d] = %g\n", n-1, n-2, outHessenbergMatrix[n - 1][n - 2]);
	if (outHessenbergMatrix[n - 1][n - 2] > 0) {
		int i;
#pragma omp parallel for private(i) shared(gmresBasis, outHessenbergMatrix, n, Nx, Ny, Nz, Nmomentum)
		for (i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				for (int k = 0; k < Nz; ++k) {
					for (int l = 0; l < Nmomentum; ++l) {
						gmresBasis->array[n - 1][i][j][k][l] /= outHessenbergMatrix[n - 1][n - 2];
					}
				}
			}
		}

		//MPI_Barrier(cartComm);
	} else {
		printf("outHessenbergMatrix[n-1][n-2] == 0\n");
		int i;
#pragma omp parallel for private(i) shared(gmresBasis, n, Nx, Ny, Nz, Nmomentum)
		for (i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				for (int k = 0; k < Nz; ++k) {
					for (int l = 0; l < Nmomentum; ++l) {
						gmresBasis->array[n - 1][i][j][k][l] = 0;
					}
				}
			}
		}
	}

	//exchangeLargeVector(gmresBasis->array[n-1], xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX, periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutGmresBuffer, rightOutGmresBuffer, leftInGmresBuffer, rightInGmresBuffer, frontOutGmresBuffer, backOutGmresBuffer, frontInGmresBuffer, backInGmresBuffer, bottomOutGmresBuffer, topOutGmresBuffer, bottomInGmresBuffer, topInGmresBuffer);
}

void generalizedMinimalResidualMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outvector, int Nx, int Ny, int Nz, int Nmomentum, double precision, int maxIteration, int verbosity, LargeVectorBasis* gmresBasis) {
	double**** initialVector = new double*** [Nx];
	for (int i = 0; i < Nx; ++i) {
		initialVector[i] = new double** [Ny];
		for (int j = 0; j < Ny; ++j) {
			initialVector[i][j] = new double* [Nz];
			for (int k = 0; k < Nz; ++k) {
				initialVector[i][j][k] = new double[Nmomentum];
				for (int l = 0; l < Nmomentum; ++l) {
					initialVector[i][j][k][l] = rightPart[i][j][k][l];
				}
			}
		}
	}

	generalizedMinimalResidualMethod(matrix, rightPart, outvector, initialVector, Nx, Ny, Nz, Nmomentum, precision, maxIteration, verbosity, gmresBasis);

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				delete[] initialVector[i][j][k];
			}
			delete[] initialVector[i][j];
		}
		delete[] initialVector[i];
	}
	delete[] initialVector;
}

void generalizedMinimalResidualMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outvector, double**** initialVector, int Nx, int Ny, int Nz, int Nmomentum, double precision, int maxIteration, int verbosity, LargeVectorBasis* gmresBasis) {

	if (verbosity > 0) printf("start GMRES\n");


	//MPI_Barrier(cartComm);
	/*exchangeLargeVector(rightPart, xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX,
	                    periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutGmresBuffer, rightOutGmresBuffer,
	                    leftInGmresBuffer, rightInGmresBuffer, frontOutGmresBuffer, backOutGmresBuffer, frontInGmresBuffer,
	                    backInGmresBuffer, bottomOutGmresBuffer, topOutGmresBuffer, bottomInGmresBuffer, topInGmresBuffer);*/

	double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, Nx, Ny, Nz, Nmomentum));
	//printf("norm = %g\n", norm);
	//alertNaNOrInfinity(norm, "right partnorm = NaN in gmres\n");

	if (norm == 0) {
		int i;
#pragma omp parallel for private(i) shared(outvector, Nx, Ny, Nz, Nmomentum)
		for (i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				for (int k = 0; k < Nz; ++k) {
					for (int l = 0; l < Nmomentum; ++l) {
						outvector[i][j][k][l] = 0.0;
					}
				}
			}
		}
		return;
	}

	//#pragma omp parallel for
	int i;
#pragma omp parallel for private(i) shared(rightPart, norm, Nx, Ny, Nz, Nmomentum)
	for (i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					rightPart[i][j][k][l] /= norm;
				}
			}
		}
	}

	double normRightPart = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, Nx, Ny, Nz, Nmomentum));

/*#pragma omp parallel for private(i) shared(rightPart, norm, Nx, Ny, Nz, Nmomentum)
	for (i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					rightPart[i][j][k][l] /= normRightPart;
				}
			}
		}
	}
	normRightPart = 1.0;*/
	int matrixDimension = Nx*Ny*Nz*Nmomentum;

	double** hessenbergMatrix;
	double** newHessenbergMatrix;
	hessenbergMatrix = new double*[1];
	hessenbergMatrix[0] = new double[1];
	hessenbergMatrix[0][0] = 0;

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
	//oldQmatrix[0] = new double[1];
	oldRmatrix[1] = new double[1];

#pragma omp parallel for private(i) shared(gmresBasis, rightPart, Nx, Ny, Nz, Nmomentum)
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					initialVector[i][j][k][l] = initialVector[i][j][k][l] / norm;
				}
			}
		}
	}

	if (gmresBasis->capacity <= 0) {
		gmresBasis->resize(10);
	}

	//double normInitialVector = sqrt(scalarMultiplyLargeVectors(initialVector, initialVector, Nx, Ny, Nz, Nmomentum));
	multiplySpecialMatrixVector(gmresBasis->array[0], matrix, initialVector, Nx, Ny, Nz, Nmomentum);
#pragma omp parallel for private(i) shared(gmresBasis, rightPart, Nx, Ny, Nz, Nmomentum)
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					gmresBasis->array[0][i][j][k][l] = rightPart[i][j][k][l] - gmresBasis->array[0][i][j][k][l];
				}
			}
		}
	}
	gmresBasis->size = 1;

	double norm1 = sqrt(scalarMultiplyLargeVectors(gmresBasis->array[0], gmresBasis->array[0], Nx, Ny, Nz, Nmomentum));
#pragma omp parallel for private(i) shared(gmresBasis, rightPart, Nx, Ny, Nz, Nmomentum)
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					gmresBasis->array[0][i][j][k][l] = gmresBasis->array[0][i][j][k][l]/norm1;
				}
			}
		}
	}

	int n = 2;
	//double beta = 1.0;
	double beta = norm1;
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
	//double maxRelativeError = precision;

	int count = 0;
	int Nrestart = 100;
	int niteration = n;

	while ((relativeError > max(maxRelativeError, 1E-15) && (niteration < min(maxIteration, matrixDimension + 3)))) {
		niteration++;
		if (verbosity > 1) printf("GMRES iteration %d\n", niteration);
		if (verbosity > 1) printf("relative error %20.15g\n", relativeError);
		newHessenbergMatrix = new double*[n];
		for (int i = 0; i < n; ++i) {
			newHessenbergMatrix[i] = new double[n - 1];
		}
		arnoldiIterations(matrix, newHessenbergMatrix, n, gmresBasis, hessenbergMatrix, Nx, Ny, Nz, Nmomentum);

		hessenbergMatrix = newHessenbergMatrix;

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
			//alertNaNOrInfinity(y[i], "y = NaN\n");
		}

		error = fabs(beta * Qmatrix[n - 1][0]);
		//double error1 = evaluateError(hessenbergMatrix, y, beta, n);

		relativeError = error / normRightPart;

		for (int i = 0; i < n - 1; ++i) {
			delete[] oldQmatrix[i];
			delete[] oldRmatrix[i];
		}
		if (n == 2) {
			delete[] oldQmatrix[1];
			delete[] oldRmatrix[1];
		}
		delete[] oldQmatrix;
		delete[] oldRmatrix;

		oldQmatrix = Qmatrix;
		oldRmatrix = Rmatrix;

		n++;
		count++;
		if (count >= Nrestart) {
			/*count = 0;
#pragma omp parallel for private(i) shared(outvector, gmresBasis, y, norm, Nx, Ny, Nz, Nmomentum)
			for (i = 0; i < Nx; ++i) {
				for (int j = 0; j < Ny; ++j) {
					for (int k = 0; k < Nz; ++k) {
						for (int l = 0; l < Nmomentum; ++l) {
							//outvector[i][j][k][l] = initialVector[i][j][k][l];
							for (int m = 0; m < n - 1; ++m) {
								initialVector[i][j][k][l] += gmresBasis->array[m][i][j][k][l] * y[m];
								//outvector[i][l] += basis[m][i][l] * y[m];
							}
						}
					}
				}
			}

			double outNorm = sqrt(scalarMultiplyLargeVectors(initialVector, initialVector, Nx, Ny, Nz, Nmomentum));

			gmresBasis->clear();
			multiplySpecialMatrixVector(gmresBasis->array[0], matrix, initialVector, Nx, Ny, Nz, Nmomentum);
#pragma omp parallel for private(i) shared(gmresBasis, rightPart, Nx, Ny, Nz, Nmomentum)
			for (int i = 0; i < Nx; ++i) {
				for (int j = 0; j < Ny; ++j) {
					for (int k = 0; k < Nz; ++k) {
						for (int l = 0; l < Nmomentum; ++l) {
							gmresBasis->array[0][i][j][k][l] = rightPart[i][j][k][l] - gmresBasis->array[0][i][j][k][l];
						}
					}
				}
			}
			gmresBasis->size = 1;

			double norm1 = sqrt(scalarMultiplyLargeVectors(gmresBasis->array[0], gmresBasis->array[0], Nx, Ny, Nz, Nmomentum));
#pragma omp parallel for private(i) shared(gmresBasis, rightPart, Nx, Ny, Nz, Nmomentum)
			for (int i = 0; i < Nx; ++i) {
				for (int j = 0; j < Ny; ++j) {
					for (int k = 0; k < Nz; ++k) {
						for (int l = 0; l < Nmomentum; ++l) {
							gmresBasis->array[0][i][j][k][l] = gmresBasis->array[0][i][j][k][l] / norm1;
						}
					}
				}
			}

			hessenbergMatrix = new double* [1];
			hessenbergMatrix[0] = new double[1];
			hessenbergMatrix[0][0] = 0;

			Qmatrix = new double* [2];
			Rmatrix = new double* [2];
			oldQmatrix = new double* [2];
			oldRmatrix = new double* [2];

			for (int i = 0; i < 2; ++i) {
				Qmatrix[i] = new double[2];
				oldQmatrix[i] = new double[2];
			}

			Rmatrix[0] = new double[1];
			Rmatrix[1] = new double[1];
			oldRmatrix[0] = new double[1];
			//oldQmatrix[0] = new double[1];
			oldRmatrix[1] = new double[1];

			n = 2;
			beta = 1.0;
			error = beta;
			y = new double[1];


			relativeError = 1;*/
		}
	}

	n = n - 1;
	if(verbosity > 0) printf("total GMRES iteration = %d\n", n);
	if (verbosity > 0) printf("GMRES relative error = %g\n", relativeError);
	if (verbosity > 0) printf("GMRES relative error * Nx * Ny * Nz * Nmomentum = %g\n", relativeError * matrixDimension);
	//out result

#pragma omp parallel for private(i) shared(outvector, gmresBasis, y, norm, Nx, Ny, Nz, Nmomentum)
	for (i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					outvector[i][j][k][l] = initialVector[i][j][k][l]*norm;
					for (int m = 0; m < n - 1; ++m) {
						outvector[i][j][k][l] += gmresBasis->array[m][i][j][k][l] * y[m] * norm;
						//outvector[i][l] += basis[m][i][l] * y[m];
					}
				}
			}
		}
	}

	multiplySpecialMatrixVector(gmresBasis->array[0], matrix, outvector, Nx, Ny, Nz, Nmomentum);
#pragma omp parallel for private(i) shared(gmresBasis, rightPart, Nx, Ny, Nz, Nmomentum)
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					gmresBasis->array[0][i][j][k][l] = rightPart[i][j][k][l]*norm - gmresBasis->array[0][i][j][k][l];
				}
			}
		}
	}

	double realError = sqrt(scalarMultiplyLargeVectors(gmresBasis->array[0], gmresBasis->array[0], Nx, Ny, Nz, Nmomentum));

	printf("real error = %g\n", realError / norm);

	gmresBasis->clear();

	//MPI_Barrier(cartComm);

	//exchangeLargeVector(outvector, xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX, periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutGmresBuffer, rightOutGmresBuffer, leftInGmresBuffer, rightInGmresBuffer, frontOutGmresBuffer, backOutGmresBuffer, frontInGmresBuffer, backInGmresBuffer, bottomOutGmresBuffer, topOutGmresBuffer, bottomInGmresBuffer, topInGmresBuffer);
	//MPI_Barrier(cartComm);
	

	for (int i = 0; i < n; ++i) {
		delete[] Qmatrix[i];
		delete[] Rmatrix[i];
		delete[] hessenbergMatrix[i];
	}
	delete[] Qmatrix;
	delete[] Rmatrix;
	delete[] hessenbergMatrix;
	delete[] y;

#pragma omp parallel for private(i) shared(rightPart, norm, Nx, Ny, Nz, Nmomentum)
	for (i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					rightPart[i][j][k][l] *= norm;
					/*for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						double value = matrix[i][j][k][l][m].value;
						//matrix[i][l][m].value *= norm;
						value = matrix[i][j][k][l][m].value;
					}*/
				}
			}
		}
	}
}

double scalarMultiplyLargeVectors(double**** a, double**** b, int Nx, int Ny, int Nz, int Nmomentum) {
	
	double result = 0;
	int i;
#pragma omp parallel for private(i) shared(a, b, Nx, Ny, Nz, Nmomentum) reduction(+:result)
	for (i = 0; i < Nx; ++i) {
		double tempResult = 0;
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					tempResult += a[i][j][k][l] * b[i][j][k][l];
					//alertNaNOrInfinity(a[i][j][k][l], "a[i][j][k][l] = NaN in scalarMult\n");
					//alertNaNOrInfinity(result[0], "result[0] = NaN in scalarMult\n");
				}
			}
		}
		result += tempResult;
	}

	//MPI_Allreduce(result, globalResult, 1, MPI_DOUBLE, MPI_SUM, cartComm);
	return result;
}

/*void simpleIterationSolver(double* outVector, double* tempVector, int N, double precision, int maxIteration) {
	
	double maxErrorLevel = precision / (xnumberGeneral * ynumberGeneral * znumberGeneral* lnumber);
	double relativeError = 1;
	int iterationCount = 0;
	double norm = rightPartEvaluator->rightPartInitialNorm();
	rightPartEvaluator->getError(outVector, tempVector);
	exchangeLargeVector(tempVector, xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX, periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutBuffer, rightOutBuffer, leftInBuffer, rightInBuffer, frontOutBuffer, backOutBuffer, frontInBuffer, backInBuffer, bottomOutBuffer, topOutBuffer, bottomInBuffer, topInBuffer);
	double errorNorm = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector, xnumberAdded, ynumberAdded, znumberAdded,
	                                         additionalBinNumber, lnumber, periodicX, periodicY, periodicZ, rank, nprocs,
	                                         cartComm, cartCoord, cartDim));
	relativeError = errorNorm/norm;
	if(norm <= 0) {
		relativeError = 0;
	}
	while (iterationCount < maxIteration && relativeError > maxErrorLevel) {
		iterationCount++;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][j][k][l] = rightPartEvaluator->rightPart(outVector, i, j, k, l);
						//alertNaNOrInfinity(tempVector[i][j][k][l], "tempvector = NaN in simpleIterationSolver\n");
					}
				}
			}
		}
		exchangeLargeVector(tempVector, xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX, periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutBuffer, rightOutBuffer, leftInBuffer, rightInBuffer, frontOutBuffer, backOutBuffer, frontInBuffer, backInBuffer, bottomOutBuffer, topOutBuffer, bottomInBuffer, topInBuffer);
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outVector[i][j][k][l] = tempVector[i][j][k][l];
					}
				}
			}
		}
		rightPartEvaluator->getError(outVector, tempVector);
		exchangeLargeVector(tempVector, xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX, periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutBuffer, rightOutBuffer, leftInBuffer, rightInBuffer, frontOutBuffer, backOutBuffer, frontInBuffer, backInBuffer, bottomOutBuffer, topOutBuffer, bottomInBuffer, topInBuffer);
		errorNorm = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector, xnumberAdded, ynumberAdded, znumberAdded,
	                                         additionalBinNumber, lnumber, periodicX, periodicY, periodicZ, rank, nprocs,
	                                         cartComm, cartCoord, cartDim));
		relativeError = errorNorm/norm;
	}
}*/

/*double normDifferenceLargeVectors(double* a, double* b, int N) {
	int minI = 1 + additionalBinNumber;
	if (cartCoord[0] == 0 && !periodicX) {
		minI = 0;
	}
	int maxI = xnumberAdded - 1 - additionalBinNumber;
	if (cartCoord[0] == cartDim[0] - 1 && !periodicX) {
		maxI = xnumberAdded - 1;
	}
	int minJ = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - 1 - additionalBinNumber;
	int minK = 1 + additionalBinNumber;
	int maxK = znumberAdded - 1 - additionalBinNumber;
	double diff[1];
	double tempDiff[1];
	diff[0] = 0;
	for (int i = minI; i <= maxI; ++i) {
		for (int j = minJ; j <= maxJ; ++j) {
			for (int k = minK; k <= maxK; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					diff[0] += sqr(a[i][j][k][l] - b[i][j][k][l]);
				}
			}
		}
	}

	MPI_Reduce(diff, tempDiff, 1, MPI_DOUBLE, MPI_SUM, 0, cartComm);
	if (rank == 0) {
		tempDiff[0] = sqrt(tempDiff[0]);
	}
	MPI_Bcast(tempDiff, 1, MPI_DOUBLE, 0, cartComm);
	return tempDiff[0];
}*/

/*void transposeSpecialMatrix(std::vector<MatrixElement>* result, std::vector<MatrixElement>* matrix, int N) {
	for (int i = 0; i < N; ++i) {
					for (int m = 0; m < matrix[i].size(); ++m) {
						MatrixElement element = matrix[i][m];

						result[element.i].push_back(MatrixElement(element.value, i));
					}
	}
}*/

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

double dichotomySolver(const double* functionValue, int minIndex, int maxIndex, const double* xValue, double y) {
	if(minIndex + 1 == maxIndex) {
		double y1 = functionValue[minIndex];
		double y2 = functionValue[maxIndex];
		double x1 = xValue[minIndex];
		double x2 = xValue[maxIndex];
		return (x2 - x1)*(y - y1)/(y2 - y1) + x1;
	}
	int middleIndex = (minIndex + maxIndex)/2;
	if(functionValue[middleIndex] > y) {
		return dichotomySolver(functionValue, minIndex, middleIndex, xValue, y);
	} else if(functionValue[middleIndex] < y) {
		return dichotomySolver(functionValue, middleIndex, maxIndex, xValue, y);
	} else {
		return xValue[middleIndex];
	}
}

double juttnerFunction(const double& u, const double& theta) {
	return exp(-sqrt(1 + u*u)/theta)*u*u;
}

/*Complex*** fastFourierTransition(double*** a, int xnumber, int ynumber, int znumber) {
	int kx = xnumber;
	while (kx > 1) {
		if (kx % 2 != 0) {
			printf("xnumber is not 2^N\n");
		}
		kx = kx / 2;
	}
	int ky = ynumber;
	while (ky > 1) {
		if (ky % 2 != 0) {
			printf("ynumber is not 2^N\n");
		}
		ky = ky / 2;
	}
	int kz = znumber;
	while (kx > 1) {
		if (kx % 2 != 0) {
			printf("znumber is not 2^N\n");
		}
		kz = kz / 2;
	}

	Complex*** result = new Complex**[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		result[i] = new Complex*[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			result[i][j] = new Complex[znumber];
			for (int k = 0; k < znumber; ++k) {
				result[i][j][k] = Complex(a[i][j][k], 0);
			}
		}
	}
	sortInputFastFourierX(result, xnumber, ynumber, znumber);

	Complex* tempResult;
	if (xnumber > 1) {
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
							tempResult[i * k + m] = result[i * k + m][ycount][zcount] + result[i * k + (k / 2) + m][ycount][zcount] *
								complexExp(
									-2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[i * k + m][ycount][zcount] - result[i * k + (k / 2) + m][ycount][zcount]
								* complexExp(
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

	if (ynumber > 1) {
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
							tempResult[i * k + m] = result[xcount][i * k + m][zcount] + result[xcount][i * k + (k / 2) + m][zcount] *
								complexExp(
									-2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[xcount][i * k + m][zcount] - result[xcount][i * k + (k / 2) + m][zcount]
								* complexExp(
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

	if (znumber > 1) {
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
							tempResult[i * k + m] = result[xcount][ycount][i * k + m] + result[xcount][ycount][i * k + (k / 2) + m] *
								complexExp(
									-2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[xcount][ycount][i * k + m] - result[xcount][ycount][i * k + (k / 2) + m]
								* complexExp(
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

void sortInputFastFourierX(double*** a, int xnumber, int ynumber, int znumber) {
	double* tempA = new double[xnumber];
	for (int ycount = 0; ycount < ynumber; ++ycount) {
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

void sortInputFastFourierY(double*** a, int xnumber, int ynumber, int znumber) {
	double* tempA = new double[ynumber];
	for (int xcount = 0; xcount < xnumber; ++xcount) {
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

void sortInputFastFourierZ(double*** a, int xnumber, int ynumber, int znumber) {
	double* tempA = new double[znumber];
	for (int xcount = 0; xcount < xnumber; ++xcount) {
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

double*** fastFourierReverceTransition(Complex*** a, int xnumber, int ynumber, int znumber) {
	int kx = xnumber;
	while (kx > 1) {
		if (kx % 2 != 0) {
			printf("xnumber is not 2^N\n");
		}
		kx = kx / 2;
	}
	int ky = ynumber;
	while (ky > 1) {
		if (ky % 2 != 0) {
			printf("ynumber is not 2^N\n");
		}
		ky = ky / 2;
	}
	int kz = znumber;
	while (kx > 1) {
		if (kx % 2 != 0) {
			printf("znumber is not 2^N\n");
		}
		kz = kz / 2;
	}

	Complex*** result = new Complex**[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		result[i] = new Complex*[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			result[i][j] = new Complex[znumber];
			for (int k = 0; k < znumber; ++k) {
				result[i][j][k] = a[i][j][k];
			}
		}
	}
	sortInputFastFourierX(result, xnumber, ynumber, znumber);

	Complex* tempResult;
	if (xnumber > 1) {
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
							tempResult[i * k + m] = result[i * k + m][ycount][zcount] + result[i * k + (k / 2) + m][ycount][zcount] *
								complexExp(
									2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[i * k + m][ycount][zcount] - result[i * k + (k / 2) + m][ycount][zcount]
								* complexExp(
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

	if (ynumber > 1) {
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
							tempResult[i * k + m] = result[xcount][i * k + m][zcount] + result[xcount][i * k + (k / 2) + m][zcount] *
								complexExp(
									2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[xcount][i * k + m][zcount] - result[xcount][i * k + (k / 2) + m][zcount]
								* complexExp(
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

	if (znumber > 1) {
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
							tempResult[i * k + m] = result[xcount][ycount][i * k + m] + result[xcount][ycount][i * k + (k / 2) + m] *
								complexExp(
									2 * pi * m / k);
							tempResult[i * k + (k / 2) + m] = result[xcount][ycount][i * k + m] - result[xcount][ycount][i * k + (k / 2) + m]
								* complexExp(
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
	for (int i = 0; i < xnumber; ++i) {
		realPart[i] = new double*[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			realPart[i][j] = new double[znumber];
			for (int k = 0; k < znumber; ++k) {
				realPart[i][j][k] = result[i][j][k].re / (xnumber * ynumber * znumber);
			}
			delete[] result[i][j];
		}
		delete[] result[i];
	}

	delete[] result;

	return realPart;
}

void sortInputFastFourierX(Complex*** a, int xnumber, int ynumber, int znumber) {
	Complex* tempA = new Complex[xnumber];
	for (int ycount = 0; ycount < ynumber; ++ycount) {
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

void sortInputFastFourierY(Complex*** a, int xnumber, int ynumber, int znumber) {
	Complex* tempA = new Complex[ynumber];
	for (int xcount = 0; xcount < xnumber; ++xcount) {
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

void sortInputFastFourierZ(Complex*** a, int xnumber, int ynumber, int znumber) {
	Complex* tempA = new Complex[znumber];
	for (int xcount = 0; xcount < xnumber; ++xcount) {
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
}*/

/*void conjugateGradientMethod(std::vector<MatrixElement>* matrix, double* rightPart, double* outVector, int N, double precision, int maxIteration, int verbosity) {
	if (verbosity > 0) printf("start conjugate gradient\n");

	double* residual = new double[N];
	double* prevResidual = new double[N];
	double* z = new double[N];
	double* tempVector = new double[N];

	double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, N));
	//printf("norm = %g\n", norm);
	//alertNaNOrInfinity(norm, "right partnorm = NaN in gmres\n");

	if (norm == 0) {
		for (int i = 0; i < N; ++i) {
			outVector[i] = 0;
		}
		return;
	}

	//#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		rightPart[i] /= norm;
	}

	for (int i = 0; i < N; ++i) {
					outVector[i] = 0;
					prevResidual[i] = rightPart[i];
					z[i] = rightPart[i];
					tempVector[i] = 0;
	}

	int iteration = 0;


	double prevResidualNorm2 = scalarMultiplyLargeVectors(prevResidual, prevResidual, N);
	double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, N);

	double relativeError = sqrt(residualNorm2 / rightPartNorm2);

	while ((iteration < maxIteration) && (iteration < N) && (
		relativeError > (precision / (N)))) {
		if (verbosity > 1) printf("conjugate gradient iteration %d\n", iteration);

		multiplySpecialMatrixVector(tempVector, matrix, z, N);

		double alpha = prevResidualNorm2 / scalarMultiplyLargeVectors(tempVector, z, N);

		for (int i = 0; i < N; ++i) {
						outVector[i] += alpha * z[i];
						residual[i] = prevResidual[i] - alpha * tempVector[i];
		}

		residualNorm2 = scalarMultiplyLargeVectors(residual, residual, N);

		double beta = residualNorm2 / prevResidualNorm2;

		for (int i = 0; i < N; ++i) {
						z[i] = residual[i] + beta * z[i];
		}

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2 / rightPartNorm2);
		iteration++;
	}

	for (int i = 0; i < N; ++i) {
		rightPart[i] *= norm;
		outVector[i] *= norm;
	}

	delete[] residual;
	delete[] prevResidual;
	delete[] z;
	delete[] tempVector;
}*/

/*void biconjugateGradientMethod(std::vector<MatrixElement>* matrix, double* rightPart, double* outVector, int N, double precision, int maxIteration, int verbosity) {
	int rank;
	int nprocs;
	MPI_Comm_size(cartComm, &nprocs);
	MPI_Comm_rank(cartComm, &rank);
	if (rank == 0 && verbosity > 0) printf("start biconjugate gradient\n");
	double**** residual = new double***[xnumberAdded + 1];
	double**** prevResidual = new double***[xnumberAdded + 1];
	double**** z = new double***[xnumberAdded + 1];
	double**** p = new double***[xnumberAdded + 1];
	double**** s = new double***[xnumberAdded + 1];
	double**** tempVector = new double***[xnumberAdded + 1];
	double**** tempVector2 = new double***[xnumberAdded + 1];
	std::vector < MatrixElement >**** transposedMatrix = new std::vector < MatrixElement >***[xnumberAdded + 1];

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		residual[i] = new double**[ynumberAdded];
		prevResidual[i] = new double**[ynumberAdded];
		z[i] = new double**[ynumberAdded];
		p[i] = new double**[ynumberAdded];
		s[i] = new double**[ynumberAdded];
		tempVector[i] = new double**[ynumberAdded];
		tempVector2[i] = new double**[ynumberAdded];
		transposedMatrix[i] = new std::vector < MatrixElement >**[ynumberAdded];
		for (int j = 0; j < ynumberAdded; ++j) {
			residual[i][j] = new double*[znumberAdded];
			prevResidual[i][j] = new double*[znumberAdded];
			z[i][j] = new double*[znumberAdded];
			p[i][j] = new double*[znumberAdded];
			s[i][j] = new double*[znumberAdded];
			tempVector[i][j] = new double*[znumberAdded];
			tempVector2[i][j] = new double*[znumberAdded];
			transposedMatrix[i][j] = new std::vector < MatrixElement >*[znumberAdded];
			for (int k = 0; k < znumberAdded; ++k) {
				residual[i][j][k] = new double[lnumber];
				prevResidual[i][j][k] = new double[lnumber];
				z[i][j][k] = new double[lnumber];
				p[i][j][k] = new double[lnumber];
				s[i][j][k] = new double[lnumber];
				tempVector[i][j][k] = new double[lnumber];
				tempVector2[i][j][k] = new double[lnumber];
				transposedMatrix[i][j][k] = new std::vector < MatrixElement >[lnumber];
				for (int l = 0; l < lnumber; ++l) {
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

	transposeSpecialMatrix(transposedMatrix, matrix, xnumberAdded, ynumberAdded, znumberAdded, lnumber);


	int iteration = 0;


	double prevResidualNorm2 = scalarMultiplyLargeVectors(p, prevResidual, xnumberAdded, ynumberAdded, znumberAdded,
	                                                      additionalBinNumber, lnumber, periodicX, periodicY, periodicZ,
	                                                      rank, nprocs, cartComm, cartCoord, cartDim);
	double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, xnumberAdded, ynumberAdded, znumberAdded,
	                                                   additionalBinNumber, lnumber, periodicX, periodicY, periodicZ, rank,
	                                                   nprocs, cartComm, cartCoord, cartDim);

	double relativeError = sqrt(residualNorm2 / rightPartNorm2);

	while ((iteration < maxIteration) && (iteration < xnumberAdded * ynumberAdded * znumberAdded * lnumber) && (
		relativeError > (precision / (xnumberAdded * ynumberAdded * znumberAdded * lnumber)))) {
		if (rank == 0 && verbosity > 1) printf("biconjugate gradient iteration %d\n", iteration);

		multiplySpecialMatrixVector(tempVector, matrix, z, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber,
		                            lnumber, periodicX, periodicY, periodicZ, rank, nprocs, cartComm, cartCoord, cartDim);
		multiplySpecialMatrixVector(tempVector2, transposedMatrix, s, xnumberAdded, ynumberAdded, znumberAdded,
		                            additionalBinNumber, lnumber, periodicX, periodicY, periodicZ, rank, nprocs, cartComm,
		                            cartCoord, cartDim);

		double alpha = prevResidualNorm2 / scalarMultiplyLargeVectors(tempVector, s, xnumberAdded, ynumberAdded, znumberAdded,
		                                                              additionalBinNumber,
		                                                              lnumber, periodicX, periodicY, periodicZ, rank, nprocs,
		                                                              cartComm, cartCoord, cartDim);

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outVector[i][j][k][l] += alpha * z[i][j][k][l];
						residual[i][j][k][l] = prevResidual[i][j][k][l] - alpha * tempVector[i][j][k][l];
						p[i][j][k][l] = p[i][j][k][l] - alpha * tempVector2[i][j][k][l];
					}
				}
			}
		}

		residualNorm2 = scalarMultiplyLargeVectors(p, residual, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber,
		                                           lnumber, periodicX, periodicY, periodicZ, rank, nprocs, cartComm,
		                                           cartCoord, cartDim);

		double beta = residualNorm2 / prevResidualNorm2;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						z[i][j][k][l] = residual[i][j][k][l] + beta * z[i][j][k][l];
						s[i][j][k][l] = p[i][j][k][l] + beta * s[i][j][k][l];
					}
				}
			}
		}

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2 / rightPartNorm2);
		iteration++;
	}

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
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
}*/

void biconjugateStabilizedGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int Nx, int Ny, int Nz, int Nmomentum, double precision, int maxIteration, int verbosity, double**** residual, double**** firstResidual, double**** v, double**** p, double**** s, double**** t, bool& converges) {
	if (verbosity > 0)printf("start biconjugate gradient\n");

	converges = false;

	double alpha = 1;
	double rho = 1;
	double omega = 1;
	int totalNumber = Nx*Ny*Nz*Nmomentum;

	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, Nx, Ny, Nz, Nmomentum);

	double tempNorm = sqrt(rightPartNorm2) / totalNumber;

	if (tempNorm <= 0) {
		for (int i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				for (int k = 0; k < Nz; ++k) {
					for (int l = 0; l < Nmomentum; ++l) {
						outVector[i][j][k][l] = 0;
					}
				}
			}
		}
		return;
	}

	multiplySpecialMatrixVector(t, matrix, outVector, Nx, Ny, Nz, Nmomentum);


	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					//outVector[i][j][k][l] = 0;
					firstResidual[i][j][k][l] = (rightPart[i][j][k][l]) - t[i][j][k][l];
					residual[i][j][k][l] = (rightPart[i][j][k][l]) - t[i][j][k][l];
					v[i][j][k][l] = 0;
					p[i][j][k][l] = 0;
					s[i][j][k][l] = 0;
					t[i][j][k][l] = 0;
				}
			}
		}
	}


	int iteration = 0;


	double prevResidualNorm2 = scalarMultiplyLargeVectors(residual, residual, Nx, Ny, Nz, Nmomentum);
	double residualNorm2 = prevResidualNorm2;

	double relativeError = sqrt(residualNorm2 / prevResidualNorm2);

	if (fabs(rightPartNorm2) < 1E-300) {
		return;
	}

	while ((iteration < maxIteration) && (iteration < totalNumber) && (relativeError > (precision / (totalNumber)))) {
		if (verbosity > 1) printf("biconjugate gradient iteration %d\n", iteration);
		if (verbosity > 1) printf("relative error %g\n", relativeError);

		double newRho = scalarMultiplyLargeVectors(firstResidual, residual, Nx, Ny, Nz, Nmomentum);
		if (fabs(rho) < 1E-100) {
			printf("rho = 0 in biconjugate\n");
			return;
		}
		if (fabs(omega) < 1E-100) {
			printf("omega = 0 in biconjugate\n");
			return;
		}
		double beta = (newRho / rho) * (alpha / omega);
		rho = newRho;

		for (int i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				for (int k = 0; k < Nz; ++k) {
					for (int l = 0; l < Nmomentum; ++l) {
						s[i][j][k][l] = 0;
					}
				}
			}
		}

		//for (int i = 0; i < xnumberAdded; ++i) {
		//for (int j = 0; j < ynumberAdded; ++j) {
		//for (int k = 0; k < znumberAdded; ++k) {
		for (int i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				for (int k = 0; k < Nz; ++k) {
					for (int l = 0; l < Nmomentum; ++l) {
						p[i][j][k][l] = residual[i][j][k][l] + beta * (p[i][j][k][l] - omega * v[i][j][k][l]);
					}
				}
			}
		}

		multiplySpecialMatrixVector(v, matrix, p, Nx, Ny, Nz, Nmomentum);

		double firstRscalarV = scalarMultiplyLargeVectors(firstResidual, v, Nx, Ny, Nz, Nmomentum);
		if (fabs(firstRscalarV) < 1E-100) {
			printf("firstRscalarV = 0 in biconjugate\n");
			return;
		}
		alpha = rho / firstRscalarV;

		//for (int i = 0; i < xnumberAdded; ++i) {
		//for (int j = 0; j < ynumberAdded; ++j) {
		//for (int k = 0; k < znumberAdded; ++k) {
		for (int i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				for (int k = 0; k < Nz; ++k) {
					for (int l = 0; l < Nmomentum; ++l) {
						s[i][j][k][l] = residual[i][j][k][l] - alpha * v[i][j][k][l];
					}
				}
			}
		}

		multiplySpecialMatrixVector(t, matrix, s, Nx, Ny, Nz, Nmomentum);

		double tnorm2 = scalarMultiplyLargeVectors(t, t, Nx, Ny, Nz, Nmomentum);
		if (tnorm2 <= 0) {
			printf("tnorm2 = 0 in biconjugate\n");
			return;
		}
		omega = scalarMultiplyLargeVectors(t, s, Nx, Ny, Nz, Nmomentum) / tnorm2;

		//for (int i = 0; i < xnumberAdded; ++i) {
		//for (int j = 0; j < ynumberAdded; ++j) {
		//for (int k = 0; k < znumberAdded; ++k) {
		for (int i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				for (int k = 0; k < Nz; ++k) {
					for (int l = 0; l < Nmomentum; ++l) {
						outVector[i][j][k][l] = outVector[i][j][k][l] + omega * s[i][j][k][l] + alpha * p[i][j][k][l];
						residual[i][j][k][l] = s[i][j][k][l] - omega * t[i][j][k][l];
					}
				}
			}
		}

		residualNorm2 = scalarMultiplyLargeVectors(residual, residual, Nx, Ny, Nz, Nmomentum);

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2 / rightPartNorm2);
		iteration++;
	}

	converges = true;

	if (relativeError > 0.1) {
		converges = false;
	}

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					//outVector[i][j][k][l] *= tempNorm;
				}
			}
		}
	}
}

void biconjugateStabilizedGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int Nx, int Ny, int Nz, int Nmomentum, double precision, int maxIteration, int verbosity, bool& converges) {
	double**** residual = new double*** [Nx];
	double**** firstResidual = new double*** [Nx];
	double**** v = new double*** [Nx];
	double**** p = new double*** [Nx];
	double**** s = new double*** [Nx];
	double**** t = new double*** [Nx];
	for (int i = 0; i < Nx; ++i) {
		residual[i] = new double** [Ny];
		firstResidual[i] = new double** [Ny];
		v[i] = new double** [Ny];
		p[i] = new double** [Ny];
		s[i] = new double** [Ny];
		t[i] = new double** [Ny];
		for (int j = 0; j < Ny; ++j) {
			residual[i][j] = new double* [Nz];
			firstResidual[i][j] = new double* [Nz];
			v[i][j] = new double* [Nz];
			p[i][j] = new double* [Nz];
			s[i][j] = new double* [Nz];
			t[i][j] = new double* [Nz];
			for (int k = 0; k < Nz; ++k) {
				residual[i][j][k] = new double[Nmomentum];
				firstResidual[i][j][k] = new double[Nmomentum];
				v[i][j][k] = new double[Nmomentum];
				p[i][j][k] = new double[Nmomentum];
				s[i][j][k] = new double[Nmomentum];
				t[i][j][k] = new double[Nmomentum];
				for (int l = 0; l < Nmomentum; ++l) {
					residual[i][j][k][l] = 0;
					firstResidual[i][j][k][l] = 0;
					v[i][j][k][l] = 0;
					p[i][j][k][l] = 0;
					s[i][j][k][l] = 0;
					t[i][j][k][l] = 0;
				}
			}
		}
	}

	biconjugateStabilizedGradientMethod(matrix, rightPart, outVector, Nx, Ny, Nz, Nmomentum, precision, maxIteration, verbosity, residual, firstResidual, v, p, s, t, converges);

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
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

/*void gaussSeidelMethod(std::vector<MatrixElement>* matrix, double* rightPart, double* outVector, int N, double precision, int maxIteration, int verbosity) {



	if (verbosity > 0) printf("start gauss-seidel\n");
	double normRightPart = scalarMultiplyLargeVectors(rightPart, rightPart, N) / (N);
	for (int i = 0; i < N; ++i) {
						outVector[i] = 0;
	}


	int curIteration = 0;
	while (curIteration < maxIteration) {
		if (verbosity > 1) printf("Gauss-Seidel iteration %d\n", curIteration);

		for (int i = 0; i < N; ++i) {
						double sum = rightPart[i];
						double a = 1;
						for (int m = 0; m < matrix[i].size(); ++m) {
							MatrixElement element = matrix[i][m];
							if (element.i != i) {
								sum -= element.value * outVector[element.i];
							} else {
								a = element.value;
							}
						}
						outVector[i] = sum / a;
		}
		//exchangeLargeVector(outVector, xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX, periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutBuffer, rightOutBuffer, leftInBuffer, rightInBuffer, frontOutBuffer, backOutBuffer, frontInBuffer, backInBuffer, bottomOutBuffer, topOutBuffer, bottomInBuffer, topInBuffer);

		curIteration++;
	}
}*/

bool indexLower(const MatrixElement& element, int i) {
	if (element.i < i) {
		return true;
	}
	if (element.i > i) {
		return false;
	}

	return false;
}

bool indexEqual(const MatrixElement& element, int i) {
	return (element.i == i);
}

bool indexUpper(const MatrixElement& element, int i) {
	return ((!indexEqual(element, i)) && (!indexLower(element, i)));
}




