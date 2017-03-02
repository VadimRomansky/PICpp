#include <stdio.h>
#include <math.h>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "complex.h"
#include "particle.h"
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "matrixElement.h"
#include "random.h"
#include "specialmath.h"

void updatePeriodicBoundaries(double**** tempVector, int xnumber, int ynumber, int znumber,
                              int lnumber) {
	for (int j = 0; j < ynumber; ++j) {
		for (int k = 0; k < znumber; ++k) {
			for (int l = 0; l < lnumber; ++l) {
				tempVector[xnumber][j][k][l] = tempVector[1][j][k][l];
			}
		}
	}
	for (int j = 0; j < ynumber; ++j) {
		for (int k = 0; k < znumber; ++k) {
			for (int l = 0; l < lnumber; ++l) {
				tempVector[0][j][k][l] = tempVector[xnumber - 1][j][k][l];
			}
		}
	}
}

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

double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, Vector3d*** vector, int xnumber, int ynumber, int znumber, int lnumber, bool periodic) {
	double**** result = new double***[xnumber + 1];
	int i = 0;

	for (i = 0; i < xnumber + 1; ++i) {
		//for (i = 1; i < xnumber; ++i) {
		result[i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			result[i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				result[i][j][k] = new double[lnumber];
				for (int l = 0; l < lnumber; ++l) {
					result[i][j][k][l] = 0;
					if ((i != 0 && i != xnumber) || ((i == 0) && (!periodic)) || ((i == xnumber) && (!periodic))) {
						for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
							MatrixElement element = matrix[i][j][k][l][m];

							result[i][j][k][l] += element.value * vector[element.i][element.j][element.k][element.l];
						}
					}
				}
			}
		}
	}

	return result;
}

double**** multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, double**** vector, int xnumber, int ynumber, int znumber, int lnumber, bool periodic) {
	double**** result = new double***[xnumber + 1];
	//printf("multiply\n");
	//for (int i = 1; i < xnumber; ++i) {
	for (int i = 0; i < xnumber + 1; ++i) {
		result[i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			result[i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				result[i][j][k] = new double[lnumber];
				for (int l = 0; l < lnumber; ++l) {
					result[i][j][k][l] = 0;
					if ((i != 0 && i != xnumber) || ((i == 0) && (!periodic)) || ((i == xnumber) && (!periodic))) {
						for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
							MatrixElement element = matrix[i][j][k][l][m];

							result[i][j][k][l] += element.value * vector[element.i][element.j][element.k][element.l];
						}
					}
				}
			}
		}
	}

	return result;
}

void multiplySpecialMatrixVector(double**** result, std::vector<MatrixElement>**** matrix, Vector3d*** vector, int xnumber, int ynumber, int znumber, int lnumber, bool periodic) {

	//for (int i = 1; i < xnumber; ++i) {
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result[i][j][k][l] = 0;
					if ((i != 0 && i != xnumber) || ((i == 0) && (!periodic)) || ((i == xnumber) && (!periodic))) {
						for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
							MatrixElement element = matrix[i][j][k][l][m];

							result[i][j][k][l] += element.value * vector[element.i][element.j][element.k][element.l];
						}
					}
				}
			}
		}
	}
}

void multiplySpecialMatrixVector(double**** result, std::vector<MatrixElement>**** matrix, double**** vector, int xnumber, int ynumber, int znumber, int lnumber, bool periodic) {


	for (int i = 0; i < xnumber + 1; ++i) {
		//for (int i = 1; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			;
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result[i][j][k][l] = 0;
					if ((i != 0 && i != xnumber) || ((i == 0) && (!periodic)) || ((i == xnumber) && (!periodic))) {
						for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
							MatrixElement element = matrix[i][j][k][l][m];

							result[i][j][k][l] += element.value * vector[element.i][element.j][element.k][element.l];
						}
					}
				}
			}
		}
	}
}

void arnoldiIterations(std::vector<MatrixElement>**** matrix, double** outHessenbergMatrix, int n,
                       LargeVectorBasis* gmresBasis, double** prevHessenbergMatrix, int xnumber, int ynumber,
                       int znumber, int lnumber, bool periodic, double* leftOutGmresBuffer, double* rightOutGmresBuffer, double* leftInGmresBuffer, double* rightInGmresBuffer) {

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
	multiplySpecialMatrixVector(gmresBasis->array[n - 1], matrix, gmresBasis->array[n - 2], xnumber, ynumber, znumber, lnumber, periodic);
	gmresBasis->size += 1;
	//printf("mult special matrix");

	//printf("start exchange\n");

	if(periodic) updatePeriodicBoundaries(gmresBasis->array[n-1], xnumber, ynumber, znumber, lnumber);

	/*if(rank == 0){
		printf("tempVector = %15.10g\n", tempVector[0][0][0][1]);
		printf("tempVector = %15.10g\n", tempVector[1][0][0][1]);
		printf("tempVector left = %15.10g\n", tempVector[xnumber - 1][0][0][1]);
		printf("tempVector left = %15.10g\n", tempVector[xnumber][0][0][1]);
	} else {

	}*/
	//printf("finish exchange\n");


	for (int m = 0; m < n - 1; ++m) {
		//double a = scalarMultiplyLargeVectors(resultBasis[m], tempVector, xnumber, ynumber, znumber, lnumber);
		outHessenbergMatrix[m][n - 2] = scalarMultiplyLargeVectors(gmresBasis->array[m], gmresBasis->array[n - 1], xnumber, ynumber,
		                                                           znumber, lnumber, periodic);
		//printf("outHessenbergMatrix[%d][%d] = %g\n", m, n-2, outHessenbergMatrix[m][n - 2]);

		//for (int i = 0; i < xnumber+1; ++i) {
		if (!periodic) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						gmresBasis->array[n - 1][0][j][k][l] -= outHessenbergMatrix[m][n - 2] * gmresBasis->array[m][0][j][k][l];
						//alertNaNOrInfinity(gmresBasis->array[n - 1][0][j][k][l], "tempVector = NaN\n");
						//printf("tempvector[%d][%d][%d][%d] = %g\n", i, j, k, l, tempVector[i][j][k][l]);
					}
				}
			}
		} else {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						gmresBasis->array[n - 1][0][j][k][l] = 0;
					}
				}
			}
		}
		for (int i = 1; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						gmresBasis->array[n - 1][i][j][k][l] -= outHessenbergMatrix[m][n - 2] * gmresBasis->array[m][i][j][k][l];
						//alertNaNOrInfinity(gmresBasis->array[n - 1][i][j][k][l], "tempVector = NaN\n");
						//printf("tempvector[%d][%d][%d][%d] = %g\n", i, j, k, l, tempVector[i][j][k][l]);
					}
				}
			}
		}
		if (!periodic) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						gmresBasis->array[n - 1][xnumber][j][k][l] -= outHessenbergMatrix[m][n - 2] * gmresBasis->array[m][xnumber][j][k][l];
						//alertNaNOrInfinity(gmresBasis->array[n - 1][xnumber][j][k][l], "tempVector = NaN\n");
						//printf("tempvector[%d][%d][%d][%d] = %g\n", i, j, k, l, tempVector[i][j][k][l]);
					}
				}
			}
		} else {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						gmresBasis->array[n - 1][xnumber][j][k][l] = 0;
					}
				}
			}
		}
	}
	//printf("finish orthogonalisation\n");
	outHessenbergMatrix[n - 1][n - 2] = sqrt(
		scalarMultiplyLargeVectors(gmresBasis->array[n - 1], gmresBasis->array[n - 1], xnumber, ynumber, znumber, lnumber, periodic));
	//printf("outHessenbergMatrix[%d][%d] = %g\n", n-1, n-2, outHessenbergMatrix[n - 1][n - 2]);
	if (outHessenbergMatrix[n - 1][n - 2] > 0) {
		if (!periodic) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						gmresBasis->array[n - 1][0][j][k][l] /= outHessenbergMatrix[n - 1][n - 2];
						//alertNaNOrInfinity(gmresBasis->array[n - 1][0][j][k][l], "tempVector = NaN\n");
						//printf("tempvector[%d][%d][%d][%d] = %g\n", i, j, k, l, tempVector[0][j][k][l]);
					}
				}
			}
		}
		//for (int i = 0; i < xnumber+1; ++i) {
		for (int i = 1; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						gmresBasis->array[n - 1][i][j][k][l] /= outHessenbergMatrix[n - 1][n - 2];
						//alertNaNOrInfinity(gmresBasis->array[n - 1][i][j][k][l], "tempVector = NaN\n");
						//printf("tempvector[%d][%d][%d][%d] = %g\n", i, j, k, l, tempVector[i][j][k][l]);
					}
				}
			}
		}
		if (!periodic) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						gmresBasis->array[n - 1][xnumber][j][k][l] /= outHessenbergMatrix[n - 1][n - 2];
						//alertNaNOrInfinity(gmresBasis->array[n - 1][xnumber][j][k][l], "tempVector = NaN\n");
						//printf("tempvector[%d][%d][%d][%d] = %g\n", i, j, k, l, tempVector[0][j][k][l]);
					}
				}
			}
		}
	} else {
		printf("outHessenbergMatrix[n-1][n-2] == 0\n");
	}
	//printf("finish normalization\n");
	//printf("arnolsi end\n");
}

void generalizedMinimalResidualMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outvector,
                                      int xnumber, int ynumber, int znumber, int lnumber, int xnumberGeneral,
                                      int znumberGeneral, int ynumberGeneral, double precision, int maxIteration,
                                      bool periodic, int verbocity, double* leftOutGmresBuffer, double* rightOutGmresBuffer, double* leftInGmresBuffer, double* rightInGmresBuffer, LargeVectorBasis* gmresBasis) {

	if ((verbocity > 0)) printf("start GMRES\n");
	if(periodic) updatePeriodicBoundaries(rightPart, xnumber, ynumber, znumber, lnumber);

	double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, ynumber, znumber, lnumber, periodic));
	//printf("norm = %g\n", norm);

	if (norm == 0) {
		for (int i = 0; i < xnumber + 1; ++i) {
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
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					rightPart[i][j][k][l] /= norm;
				}
			}
		}
	}


	int matrixDimension = lnumber * xnumberGeneral * ynumberGeneral * znumberGeneral;

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

	/*double***** basis = new double****[1];
	basis[0] = new double***[xnumber + 1];
	for (int i = 0; i < xnumber + 1; ++i) {
		//for (int i = 1; i < xnumber; ++i) {
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
	double***** newBasis;*/
	if (gmresBasis->capacity <= 0) {
		gmresBasis->resize(10);
	}
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					gmresBasis->array[0][i][j][k][l] = rightPart[i][j][k][l];
				}
			}
		}
	}

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
		if ((verbocity > 1)) printf("GMRES iteration %d\n", n);
		newHessenbergMatrix = new double*[n];
		for (int i = 0; i < n; ++i) {
			newHessenbergMatrix[i] = new double[n - 1];
		}
		arnoldiIterations(matrix, newHessenbergMatrix, n, gmresBasis, hessenbergMatrix, xnumber, ynumber, znumber,
		                  lnumber, periodic, leftOutGmresBuffer, rightOutGmresBuffer, leftInGmresBuffer, rightInGmresBuffer);

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
			alertNaNOrInfinity(y[i], "y = NaN\n");
		}

		error = fabs(beta * Qmatrix[n - 1][0]);
		//for (int i = 0; i < xnumber+1; ++i) {
		/*for (int j = 0; j < ynumber; ++j) {
		    for (int k = 0; k < znumber; ++k) {
		        for (int l = 0; l < lnumber; ++l) {
		            outvector[0][j][k][l] = 0;
		            if(rank == 0 && !periodic){
						for (int m = 0; m < n - 1; ++m) {
							outvector[0][j][k][l] += basis[m][0][j][k][l] * y[m] * norm;
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
		for (int i = 1; i < xnumber; ++i) {
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

		for (int j = 0; j < ynumber; ++j) {
		    for (int k = 0; k < znumber; ++k) {
		        for (int l = 0; l < lnumber; ++l) {
		            outvector[xnumber][j][k][l] = 0;
		            if(rank == nprocs - 1 && !periodic){
						for (int m = 0; m < n - 1; ++m) {
							outvector[xnumber][j][k][l] += basis[m][xnumber][j][k][l] * y[m] * norm;
							//printf("outvector[%d][%d][%d][%d] %d = %g\n", i, j, k, l, m, outvector[i][j][k][l]);
							//printf("norm = %g\n", norm);
							//printf("y[%d] = %g\n", m, y[m]);
							//outvector[i][l] += basis[m][i][l] * y[m];
						}
						//printf("outvector[%d][%d][%d][%d] = %g\n", i, j, k, l, outvector[i][j][k][l]);
					}
				}
			}
		}*/

		double normRightPart = sqrt(
			scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, ynumber, znumber, lnumber, periodic));
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
	}

	n = n - 1;
	//if(rank == 0) printf("total GMRES iteration = %d\n", n);

	//out result

	for (int j = 0; j < ynumber; ++j) {
		for (int k = 0; k < znumber; ++k) {
			for (int l = 0; l < lnumber; ++l) {
				outvector[0][j][k][l] = 0;
				if (!periodic) {
					for (int m = 0; m < n - 1; ++m) {
						outvector[0][j][k][l] += gmresBasis->array[m][0][j][k][l] * y[m] * norm;
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
	for (int i = 1; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					outvector[i][j][k][l] = 0;
					for (int m = 0; m < n - 1; ++m) {
						outvector[i][j][k][l] += gmresBasis->array[m][i][j][k][l] * y[m] * norm;
						//outvector[i][l] += basis[m][i][l] * y[m];
					}
				}
			}
		}
	}
	if (!periodic) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					outvector[xnumber][j][k][l] = 0;
					for (int m = 0; m < n - 1; ++m) {
						outvector[xnumber][j][k][l] += gmresBasis->array[m][xnumber][j][k][l] * y[m] * norm;
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

	gmresBasis->clear();

	//printf("rank = %d outvector[0][0][0][1] = %g\n", rank, outvector[0][0][0][1]);
	//printf("rank = %d outvector[1][0][0][1] = %g\n", rank, outvector[1][0][0][1]);
	//printf("rank = %d outvector[2][0][0][1] = %g\n", rank, outvector[2][0][0][1]);
	//printf("rank = %d outvector[xnumber - 2][0][0][1] = %g\n", rank, outvector[xnumber - 2][0][0][1]);
	//printf("rank = %d outvector[xnumber - 1][0][0][1] = %g\n", rank, outvector[xnumber - 1][0][0][1]);
	//printf("rank = %d outvector[xnumber][0][0][1] = %g\n", rank, outvector[xnumber][0][0][1]);

	if(periodic) updatePeriodicBoundaries(outvector, xnumber, ynumber, znumber, lnumber);


	//printf("rank = %d outvector[0][0][0][1] = %g\n", rank, outvector[0][0][0][1]);
	//printf("rank = %d outvector[1][0][0][1] = %g\n", rank, outvector[1][0][0][1]);
	//printf("rank = %d outvector[2][0][0][1] = %g\n", rank, outvector[2][0][0][1]);
	//printf("rank = %d outvector[xnumber - 2][0][0][1] = %g\n", rank, outvector[xnumber - 2][0][0][1]);
	//printf("rank = %d outvector[xnumber - 1][0][0][1] = %g\n", rank, outvector[xnumber - 1][0][0][1]);
	//printf("rank = %d outvector[xnumber][0][0][1] = %g\n", rank, outvector[xnumber][0][0][1]);

	for (int i = 0; i < n; ++i) {
		delete[] Qmatrix[i];
		delete[] Rmatrix[i];
		delete[] hessenbergMatrix[i];
	}
	delete[] Qmatrix;
	delete[] Rmatrix;
	delete[] hessenbergMatrix;
	delete[] y;

	for (int i = 0; i < xnumber + 1; ++i) {
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

double scalarMultiplyLargeVectors(double**** a, double**** b, int xnumber, int ynumber, int znumber, int lnumber,
                                  bool periodic) {
	double result[1];
	double globalResult[1];
	globalResult[0] = 0;
	result[0] = 0;
	if ((!periodic)) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result[0] += a[0][j][k][l] * b[0][j][k][l];
				}
			}
		}
	}
	for (int i = 1; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result[0] += a[i][j][k][l] * b[i][j][k][l];
				}
			}
		}
	}
	if ((!periodic)) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result[0] += a[xnumber][j][k][l] * b[xnumber][j][k][l];
				}
			}
		}
	}

	return result[0];
	/*if(nprocs > 1) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank != 0) {
			MPI_Send(result, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD);
			//printf("send double rank = %d\n", rank);
		} else {
			for (int i = 1; i < nprocs; ++i) {
				double temp[1];
				MPI_Status status;
				MPI_Recv(temp, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
				result[0] += temp[0];
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
			double temp[1];
			MPI_Status status;
			MPI_Recv(temp, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
			//printf("recv double rank = %d\n", rank);
			result[0] = temp[0];
		}
	}
	double res = result[0];
	MPI_Barrier(MPI_COMM_WORLD);
	return res;*/
}

double scalarMultiplyLargeVectors(Vector3d*** a, Vector3d*** b, int xnumber, int ynumber, int znumber, int lnumber,
                                  bool periodic, int rank, int nprocs) {
	double result[1];
	double globalResult[1];
	result[0] = 0;
	globalResult[0] = 0;
	if (rank == 0 && !periodic) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result[0] += a[0][j][k][l] * b[0][j][k][l];
				}
			}
		}
	}
	for (int i = 1; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result[0] += a[i][j][k][l] * b[i][j][k][l];
				}
			}
		}
	}
	if (rank == nprocs - 1 && !periodic) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					result[0] += a[xnumber][j][k][l] * b[xnumber][j][k][l];
				}
			}
		}
	}

	return result[0];

	/*if (nprocs > 1) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank != 0) {
			MPI_Send(result, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD);
		} else {
			for (int i = 1; i < nprocs; ++i) {
				double temp[1];
				MPI_Status status;
				MPI_Recv(temp, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
				result[0] += temp[0];

			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0) {
			for (int i = 1; i < nprocs; ++i) {
				MPI_Send(result, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD);
			}
		} else {
			double temp[1];
			MPI_Status status;
			MPI_Recv(temp, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_NUMBER_ALL, MPI_COMM_WORLD, &status);
			result[0] = temp[0];
		}
	}
	double res = result[0];
	return res;*/
}

void transposeSpecialMatrix(std::vector<MatrixElement>**** result, std::vector<MatrixElement>**** matrix, int xnumber, int ynumber, int znumber, int lnumber) {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
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

Complex*** fastFourierTransition(double*** a, int xnumber, int ynumber, int znumber) {
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
}

void conjugateGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector,
                             int xnumber, int ynumber, int znumber, int lnumber, double precision, int maxIteration, bool periodic, int verbosity) {
	if (verbosity > 0) printf("start conjugate gradient\n");

	double**** residual = new double***[xnumber + 1];
	double**** prevResidual = new double***[xnumber + 1];
	double**** z = new double***[xnumber + 1];
	double**** tempVector = new double***[xnumber + 1];

	for (int i = 0; i < xnumber + 1; ++i) {
		residual[i] = new double**[ynumber];
		prevResidual[i] = new double**[ynumber];
		z[i] = new double**[ynumber];
		tempVector[i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			residual[i][j] = new double*[znumber];
			prevResidual[i][j] = new double*[znumber];
			z[i][j] = new double*[znumber];
			tempVector[i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				residual[i][j][k] = new double[lnumber];
				prevResidual[i][j][k] = new double[lnumber];
				z[i][j][k] = new double[lnumber];
				tempVector[i][j][k] = new double[lnumber];
				for (int l = 0; l < lnumber; ++l) {
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


	double prevResidualNorm2 = scalarMultiplyLargeVectors(prevResidual, prevResidual, xnumber, ynumber, znumber,
	                                                      lnumber, false);
	double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, ynumber, znumber, lnumber, false);

	double relativeError = sqrt(residualNorm2 / rightPartNorm2);

	while ((iteration < maxIteration) && (iteration < xnumber * ynumber * znumber * lnumber) && (relativeError > (precision / (xnumber * ynumber * znumber * lnumber)))) {
		if (verbosity > 1) printf("conjugate gradient iteration %d\n", iteration);

		multiplySpecialMatrixVector(tempVector, matrix, z, xnumber, ynumber, znumber, lnumber, periodic);

		double alpha = prevResidualNorm2 / scalarMultiplyLargeVectors(tempVector, z, xnumber, ynumber, znumber, lnumber,
		                                                              false);

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outVector[i][j][k][l] += alpha * z[i][j][k][l];
						residual[i][j][k][l] = prevResidual[i][j][k][l] - alpha * tempVector[i][j][k][l];
					}
				}
			}
		}

		residualNorm2 = scalarMultiplyLargeVectors(residual, residual, xnumber, ynumber, znumber, lnumber, false);

		double beta = residualNorm2 / prevResidualNorm2;

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						z[i][j][k][l] = residual[i][j][k][l] + beta * z[i][j][k][l];
					}
				}
			}
		}

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2 / rightPartNorm2);
		iteration++;
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
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

void biconjugateGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int xnumber, int ynumber, int znumber, int lnumber, double precision, int maxIteration, bool periodic, int verbosity) {
	if (verbosity > 0) printf("start biconjugate gradient\n");
	double**** residual = new double***[xnumber + 1];
	double**** prevResidual = new double***[xnumber + 1];
	double**** z = new double***[xnumber + 1];
	double**** p = new double***[xnumber + 1];
	double**** s = new double***[xnumber + 1];
	double**** tempVector = new double***[xnumber + 1];
	double**** tempVector2 = new double***[xnumber + 1];
	std::vector<MatrixElement>**** transposedMatrix = new std::vector<MatrixElement>***[xnumber + 1];

	for (int i = 0; i < xnumber + 1; ++i) {
		residual[i] = new double**[ynumber];
		prevResidual[i] = new double**[ynumber];
		z[i] = new double**[ynumber];
		p[i] = new double**[ynumber];
		s[i] = new double**[ynumber];
		tempVector[i] = new double**[ynumber];
		tempVector2[i] = new double**[ynumber];
		transposedMatrix[i] = new std::vector<MatrixElement>**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			residual[i][j] = new double*[znumber];
			prevResidual[i][j] = new double*[znumber];
			z[i][j] = new double*[znumber];
			p[i][j] = new double*[znumber];
			s[i][j] = new double*[znumber];
			tempVector[i][j] = new double*[znumber];
			tempVector2[i][j] = new double*[znumber];
			transposedMatrix[i][j] = new std::vector<MatrixElement>*[znumber];
			for (int k = 0; k < znumber; ++k) {
				residual[i][j][k] = new double[lnumber];
				prevResidual[i][j][k] = new double[lnumber];
				z[i][j][k] = new double[lnumber];
				p[i][j][k] = new double[lnumber];
				s[i][j][k] = new double[lnumber];
				tempVector[i][j][k] = new double[lnumber];
				tempVector2[i][j][k] = new double[lnumber];
				transposedMatrix[i][j][k] = new std::vector<MatrixElement>[lnumber];
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

	transposeSpecialMatrix(transposedMatrix, matrix, xnumber, ynumber, znumber, lnumber);


	int iteration = 0;


	double prevResidualNorm2 = scalarMultiplyLargeVectors(p, prevResidual, xnumber, ynumber, znumber, lnumber, false);
	double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, ynumber, znumber, lnumber, false);

	double relativeError = sqrt(residualNorm2 / rightPartNorm2);

	while ((iteration < maxIteration) && (iteration < xnumber * ynumber * znumber * lnumber) && (relativeError > (precision / (xnumber * ynumber * znumber * lnumber)))) {
		if (verbosity > 1) printf("biconjugate gradient iteration %d\n", iteration);

		multiplySpecialMatrixVector(tempVector, matrix, z, xnumber, ynumber, znumber, lnumber, periodic);
		multiplySpecialMatrixVector(tempVector2, transposedMatrix, s, xnumber, ynumber, znumber, lnumber, periodic);

		double alpha = prevResidualNorm2 / scalarMultiplyLargeVectors(tempVector, s, xnumber, ynumber, znumber, lnumber,
		                                                              false);

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outVector[i][j][k][l] += alpha * z[i][j][k][l];
						residual[i][j][k][l] = prevResidual[i][j][k][l] - alpha * tempVector[i][j][k][l];
						p[i][j][k][l] = p[i][j][k][l] - alpha * tempVector2[i][j][k][l];
					}
				}
			}
		}

		residualNorm2 = scalarMultiplyLargeVectors(p, residual, xnumber, ynumber, znumber, lnumber, false);

		double beta = residualNorm2 / prevResidualNorm2;

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
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

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
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

void biconjugateStabilizedGradientMethod(std::vector<MatrixElement>**** matrix, double**** rightPart,
                                         double**** outVector, int xnumber, int ynumber, int znumber, int lnumber,
                                         double precision, int maxIteration, bool periodic, int verbosity) {
	if (verbosity > 0)printf("start biconjugate gradient\n");
	double**** residual = new double***[xnumber + 1];
	double**** firstResidual = new double***[xnumber + 1];
	double**** p = new double***[xnumber + 1];
	double**** v = new double***[xnumber + 1];
	double**** s = new double***[xnumber + 1];
	double**** t = new double***[xnumber + 1];

	double alpha = 1;
	double rho = 1;
	double omega = 1;

	for (int i = 0; i < xnumber + 1; ++i) {
		residual[i] = new double**[ynumber];
		firstResidual[i] = new double**[ynumber];
		v[i] = new double**[ynumber];
		p[i] = new double**[ynumber];
		s[i] = new double**[ynumber];
		t[i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			residual[i][j] = new double*[znumber];
			firstResidual[i][j] = new double*[znumber];
			v[i][j] = new double*[znumber];
			p[i][j] = new double*[znumber];
			s[i][j] = new double*[znumber];
			t[i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				residual[i][j][k] = new double[lnumber];
				firstResidual[i][j][k] = new double[lnumber];
				v[i][j][k] = new double[lnumber];
				p[i][j][k] = new double[lnumber];
				s[i][j][k] = new double[lnumber];
				t[i][j][k] = new double[lnumber];
				for (int l = 0; l < lnumber; ++l) {
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


	double prevResidualNorm2 = scalarMultiplyLargeVectors(residual, residual, xnumber, ynumber, znumber, lnumber, false);
	double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, ynumber, znumber, lnumber, false);

	double relativeError = sqrt(residualNorm2 / rightPartNorm2);

	while ((iteration < maxIteration) && (iteration < xnumber * ynumber * znumber * lnumber) && (relativeError > (precision / (xnumber * ynumber * znumber * lnumber)))) {
		if (verbosity > 1) printf("biconjugate gradient iteration %d\n", iteration);

		double newRho = scalarMultiplyLargeVectors(firstResidual, residual, xnumber, ynumber, znumber, lnumber, false);

		double beta = (newRho / rho) * (alpha / omega);
		rho = newRho;

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						p[i][j][k][l] = residual[i][j][k][l] + beta * (p[i][j][k][l] - omega * v[i][j][k][l]);
					}
				}
			}
		}
		multiplySpecialMatrixVector(v, matrix, p, xnumber, ynumber, znumber, lnumber, periodic);

		alpha = rho / scalarMultiplyLargeVectors(firstResidual, v, xnumber, ynumber, znumber, lnumber, false);

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						s[i][j][k][l] = residual[i][j][k][l] - alpha * v[i][j][k][l];
					}
				}
			}
		}

		multiplySpecialMatrixVector(t, matrix, s, xnumber, ynumber, znumber, lnumber, periodic);

		omega = scalarMultiplyLargeVectors(t, s, xnumber, ynumber, znumber, lnumber, false) / scalarMultiplyLargeVectors(
			t, t,
			xnumber,
			ynumber,
			znumber,
			lnumber,
			false);

		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outVector[i][j][k][l] = outVector[i][j][k][l] + omega * s[i][j][k][l] + alpha * p[i][j][k][l];
						residual[i][j][k][l] = s[i][j][k][l] - omega * t[i][j][k][l];
					}
				}
			}
		}

		residualNorm2 = scalarMultiplyLargeVectors(residual, residual, xnumber, ynumber, znumber, lnumber, false);

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2 / rightPartNorm2);
		iteration++;
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
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

void gaussSeidelMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, double**** outVector, int xnumber,
                       int ynumber, int znumber, int lnumber, int xnumberGeneral, int znumberGeneral,
                       int ynumberGeneral, double precision, int maxIteration, bool periodic, int verbocity) {
	double* bufferRightSend = new double[ynumber * znumber * lnumber];
	double* bufferRightRecv = new double[ynumber * znumber * lnumber];
	double* bufferLeftSend = new double[ynumber * znumber * lnumber];
	double* bufferLeftRecv = new double[ynumber * znumber * lnumber];
	if ((verbocity > 0)) printf("start gauss-seidel\n");
	double normRightPart = scalarMultiplyLargeVectors(rightPart, rightPart, xnumber, ynumber, znumber, lnumber, periodic) / (xnumberGeneral * ynumberGeneral * znumberGeneral * lnumber);
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					//outVector[i][j][k][l] = uniformDistribution()*normRightPart/matrix[0][0][0][0][0].vfalue;
					outVector[i][j][k][l] = 0;
				}
			}
		}
	}


	int curIteration = 0;
	while (curIteration < maxIteration) {
		if ((verbocity > 1)) printf("Gauss-Seidel iteration %d\n", curIteration);
		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						double sum = rightPart[i][j][k][l];
						double a = 1;
						for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
							MatrixElement element = matrix[i][j][k][l][m];
							if (!indexEqual(element, i, j, k, l)) {
								sum -= element.value * outVector[element.i][element.j][element.k][element.l];
							} else {
								a = element.value;
							}
						}
						outVector[i][j][k][l] = sum / a;
					}
				}
			}
		}

		curIteration++;
	}

	delete[] bufferLeftSend;
	delete[] bufferRightRecv;
	delete[] bufferRightSend;
	delete[] bufferLeftRecv;
}

bool indexLower(const MatrixElement& element, int i, int j, int k, int l) {
	if (element.i < i) {
		return true;
	}
	if (element.i > i) {
		return false;
	}
	if (element.j < j) {
		return true;
	}
	if (element.j > j) {
		return false;
	}
	if (element.k < k) {
		return true;
	}
	if (element.k > k) {
		return false;
	}
	if (element.l < l) {
		return true;
	}
	if (element.l > l) {
		return false;
	}

	return false;
}

bool indexEqual(const MatrixElement& element, int i, int j, int k, int l) {
	return (element.i == i) && (element.j == j) && (element.k == k) && (element.l == l);
}

bool indexUpper(const MatrixElement& element, int i, int j, int k, int l) {
	return ((!indexEqual(element, i, j, k, l)) && (!indexLower(element, i, j, k, l)));
}




