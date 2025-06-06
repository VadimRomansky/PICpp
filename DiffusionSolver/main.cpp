// DiffusionSolver.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include "stdio.h"
#include "math.h"
#include "mpi.h"

#include "./Math/largeVectorBasis.h"
#include "./Math/matrixElement.h"
#include "./Math/specialmath.h"
#include "./Math/util.h"

/* test matrix
 * 2 1 0 0 0 0 0 0 0
 * 2 2 -1 0 0 0 0 0 0
 * 0 1 1 -2 0 0 0 0 0
 * 0 0 1 2 3 0 0 0 0
 * 0 0 0 -1 1 2 0 0 0
 * 0 0 0 0 0 3 1 0 0
 * 0 0 0 0 0 2 1 -3 0
 * 0 0 0 0 0 0 1 2 2
 * 0 0 0 0 0 0 0 2 1
 *
 * right part
 * 1
 * 2
 * 3
 * 4
 * 5
 * 6
 * 7
 * 8
 * 9
 *
 * x
 *
 * -3.07143
 * 7.14286
 * 6.14286
 * 5.14286
 * -4.14286
 * 7.14286
 * -15.4286
 * -2.71429
 * 14.4286
 */

double evaluateDiffusionCoefficient(double p) {
	return 0.001 * p;
}

void sequentialThreeDiagonalSolverX(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum) {
	//double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);

	for (int j = 0; j < Ny; ++j) {
		for (int k = 0; k < Nz; ++k) {
			for (int l = 0; l < Nmomentum; ++l) {

				double normRightPart = 0;

				for (int i = 0; i < Nx; ++i) {
					normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
				}

				if (normRightPart <= 0) {
					for (int i = 0; i < Nx; ++i) {
						x[k][j][i][l] = 0;
					}

					continue;
				}
				
				double u = a[k][j][0][l]/b[k][j][0][l];
				double v = c[k][j][Nx - 1][l]/b[k][j][Nx - 1][l];
				for (int i = 0; i < Nx; ++i) {
					a[k][j][i][l] /= b[k][j][i][l];
					c[k][j][i][l] /= b[k][j][i][l];
					rightPart[k][j][i][l] /= b[k][j][i][l];
					b[k][j][i][l] = 1.0;
				}
				//double* d = new double[Nx];
				//d[0] = u;
				//d[1] = 0;

				for (int i = 2; i < Nx; ++i) {
					double r = 1.0 / (1 - a[k][j][i][l] * c[k][j][i - 1][l]);
					//d[i] = -r * a[k][j][i][l] * d[i - 1];
					rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j][i - 1][l]);
					a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j][i - 1][l];
					if (i == Nx - 1) {
						a[k][j][i][l] += v * r;
					}
					c[k][j][i][l] = r * c[k][j][i][l];

				}

				for (int i = Nx - 3; i >= 1; i = i - 1) {
					rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j][i + 1][l] * c[k][j][i][l];
					a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j][i + 1][l];
					c[k][j][i][l] = -c[k][j][i][l] * c[k][j][i + 1][l];
				}

				double r = 1.0 / (1.0 - a[k][j][1][l] * c[k][j][0][l]);
				rightPart[k][j][0][l] = r * (rightPart[k][j][0][l] - rightPart[k][j][1][l] * c[k][j][0][l]);
				c[k][j][0][l] = r * (u - c[k][j][0][l] * c[k][j][1][l]);

				double a1 = 1.0;
				double c1 = c[k][j][0][l];
				double d1 = rightPart[k][j][0][l];

				double a2 = a[k][j][Nx - 1][l];
				double c2 = 1.0;
				double d2 = rightPart[k][j][Nx - 1][l];

				double y2 = (d2 - d1 * a2) / (c2 - c1 * a2);
				double y1 = d1 - c1 * y2;

				x[k][j][0][l] = y1;
				x[k][j][Nx - 1][l] = y2;

				for (int i = 1; i < Nx - 1; ++i) {
					x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * y1 - c[k][j][i][l] * y2;
					if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
						printf("x = NaN in solver X, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
						exit(0);
					}
				}

				//delete[] d;
			}
		}
	}
}

void sequentialThreeDiagonalSolverY(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum) {
	//double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);

	for (int i = 0; i < Nx; ++i) {
		for (int k = 0; k < Nz; ++k) {
			for (int l = 0; l < Nmomentum; ++l) {

				double normRightPart = 0;

				for (int j = 0; j < Ny; ++j) {
					normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
				}

				if (normRightPart <= 0) {
					for (int j = 0; j < Ny; ++j) {
						x[k][j][i][l] = 0;
					}

					continue;
				}

				double u = a[k][0][i][l]/b[k][0][i][l];
				double v = c[k][Ny - 1][i][l]/b[k][Ny - 1][i][l];
				for (int j = 0; j < Ny; ++j) {
					a[k][j][i][l] /= b[k][j][i][l];
					c[k][j][i][l] /= b[k][j][i][l];
					rightPart[k][j][i][l] /= b[k][j][i][l];
					b[k][j][i][l] = 1.0;
				}
				//double* d = new double[Ny];
				//d[0] = u;
				//d[1] = 0;

				for (int j = 2; j < Ny; ++j) {
					double r = 1.0 / (1 - a[k][j][i][l] * c[k][j - 1][i][l]);
					//d[j] = -r * a[k][j][i][l] * d[j - 1];
					rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j - 1][i][l]);
					a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j- 1][i][l];
					if (j == Ny - 1) {
						a[k][j][i][l] += v * r;
					}
					c[k][j][i][l] = r * c[k][j][i][l];

				}

				for (int j = Ny - 3; j >= 1; j = j - 1) {
					rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j + 1][i][l] * c[k][j][i][l];
					a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j + 1][i][l];
					c[k][j][i][l] = - c[k][j][i][l] * c[k][j + 1][i][l];
				}

				double r = 1.0 / (1.0 - a[k][1][i][l] * c[k][0][i][l]);
				rightPart[k][0][i][l] = r * (rightPart[k][0][i][l] - rightPart[k][1][i][l] * c[k][0][i][l]);
				c[k][0][i][l] = r * (u - c[k][0][i][l] * c[k][1][i][l]);

				double a1 = 1.0;
				double c1 = c[k][0][i][l];
				double d1 = rightPart[k][0][i][l];

				double a2 = a[k][Ny - 1][i][l];
				double c2 = 1.0;
				double d2 = rightPart[k][Ny - 1][i][l];

				double y2 = (d2 - d1 * a2) / (c2 - c1 * a2);
				double y1 = d1 - c1 * y2;

				x[k][0][i][l] = y1;
				x[k][Ny - 1][i][l] = y2;

				for (int j = 1; j < Ny - 1; ++j) {
					x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * y1 - c[k][j][i][l] * y2;
					if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
						printf("x = NaN in solver Y, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
						exit(0);
					}
				}

				//delete[] d;
			}
		}
	}
}

void sequentialThreeDiagonalSolverZ(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum) {
	//double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int l = 0; l < Nmomentum; ++l) {

				double normRightPart = 0;

				for (int k = 0; k < Nz; ++k) {
					normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
				}

				if (normRightPart <= 0) {
					for (int k = 0; k < Nz; ++k) {
						x[k][j][i][l] = 0;
					}

					continue;
				}
				
				double u = a[0][j][i][l]/b[0][j][i][l];
				double v = c[Nz - 1][j][i][l]/b[Nz - 1][j][i][l];
				for (int k = 0; k < Nz; ++k) {
					a[k][j][i][l] /= b[k][j][i][l];
					c[k][j][i][l] /= b[k][j][i][l];
					rightPart[k][j][i][l] /= b[k][j][i][l];
					b[k][j][i][l] = 1.0;
				}
				//double* d = new double[Nz];
				//d[0] = u;
				//d[1] = 0;

				for (int k = 2; k < Nz; ++k) {
					double r = 1.0 / (1 - a[k][j][i][l] * c[k - 1][j][i][l]);
					//d[k] = -r * a[k][j][i][l] * d[k - 1];
					rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k - 1][j][i][l]);
					a[k][j][i][l] = -r * a[k][j][i][l] * a[k - 1][j][i][l];
					if (k == Nz - 1) {
						a[k][j][i][l] += v * r;
					}
					c[k][j][i][l] = r * c[k][j][i][l];

				}

				for (int k = Nz - 3; k >= 1; k = k - 1) {
					rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k + 1][j][i][l] * c[k][j][i][l];
					a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k + 1][j][i][l];
					c[k][j][i][l] =  - c[k][j][i][l] * c[k + 1][j][i][l];
				}

				double r = 1.0 / (1.0 - a[1][j][i][l] * c[0][j][i][l]);
				rightPart[0][j][i][l] = r * (rightPart[0][j][i][l] - rightPart[1][j][i][l] * c[0][j][i][l]);
				c[0][j][i][l] = r * (u - c[0][j][i][l] * c[1][j][i][l]);

				double a1 = 1.0;
				double c1 = c[0][j][i][l];
				double d1 = rightPart[0][j][i][l];

				double a2 = a[Nz - 1][j][i][l];
				double c2 = 1.0;
				double d2 = rightPart[Nz - 1][j][i][l];

				double y2 = (d2 - d1 * a2) / (c2 - c1 * a2);
				double y1 = d1 - c1 * y2;

				x[0][j][i][l] = y1;
				x[Nz - 1][j][i][l] = y2;

				for (int k = 1; k < Nz - 1; ++k) {
					x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * y1 - c[k][j][i][l] * y2;
					if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
						printf("x = NaN in solver Z, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
						exit(0);
					}
				}

				//delete[] d;
			}
		}
	}
}

void sequentialThreeDiagonalSolverP(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum) {
	//double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);

	for (int j = 0; j < Ny; ++j) {
		for (int k = 0; k < Nz; ++k) {
			for (int i = 0; i < Nx; ++i) {

				double normRightPart = 0;

				for (int l = 0; l < Nmomentum; ++l) {
					normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
				}

				if (normRightPart <= 0) {
					for (int l = 0; l < Nmomentum; ++l) {
						x[k][j][i][l] = 0;
					}

					continue;
				}

				double u = a[k][j][i][0] / b[k][j][i][0];
				double v = c[k][j][i][Nmomentum - 1] / b[k][j][i][Nmomentum - 1];
				for (int l = 0; l < Nmomentum; ++l) {
					a[k][j][i][l] /= b[k][j][i][l];
					c[k][j][i][l] /= b[k][j][i][l];
					rightPart[k][j][i][l] /= b[k][j][i][l];
					b[k][j][i][l] = 1.0;

					if ((a[k][j][i][l] != a[k][j][i][l]) || (0 * a[k][j][i][l] != 0 * a[k][j][i][l])) {
						printf("a = NaN k = %d, j = %d, i = %d, l = %d\n", k, j, i, l);
						exit(0);
					}

					if ((c[k][j][i][l] != c[k][j][i][l]) || (0 * c[k][j][i][l] != 0 * c[k][j][i][l])) {
						printf("c = NaN k = %d, j = %d, i = %d, l = %d\n", k, j, i, l);
						exit(0);
					}

					if ((rightPart[k][j][i][l] != rightPart[k][j][i][l]) || (0 * rightPart[k][j][i][l] != 0 * rightPart[k][j][i][l])){
						printf("rightPart = NaN k = %d, j = %d, i = %d, l = %d\n", k, j, i, l);
						exit(0);
					}
				}
				//double* d = new double[Nx];
				//d[0] = u;
				//d[1] = 0;

				for (int l = 2; l < Nmomentum; ++l) {
					double r = 1.0 / (1 - a[k][j][i][l] * c[k][j][i][l-1]);
					//d[i] = -r * a[k][j][i][l] * d[i - 1];
					rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j][i][l-1]);
					a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j][i][l-1];
					if (l == Nmomentum - 1) {
						a[k][j][i][l] += v * r;
					}
					c[k][j][i][l] = r * c[k][j][i][l];

				}

				for (int l = Nmomentum - 3; l >= 1; l = l - 1) {
					rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j][i][l+1] * c[k][j][i][l];
					a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j][i][l+1];
					c[k][j][i][l] = -c[k][j][i][l] * c[k][j][i][l+1];
				}

				double y1;
				double y2;

				if (a[k][j][i][1] * c[k][j][i][0] == 1.0) {
					rightPart[k][j][i][0] = rightPart[k][j][i][0] - rightPart[k][j][i][1] * c[k][j][i][0];
					c[k][j][i][0] = u - c[k][j][i][0] * c[k][j][i][1];
					double a1 = 0.0;
					double c1 = c[k][j][i][0];
					double d1 = rightPart[k][j][i][0];

					double a2 = a[k][j][i][Nmomentum - 1];
					double c2 = 1.0;
					double d2 = rightPart[k][j][i][Nmomentum - 1];

					y2 = (d2 - d1 * a2) / (c2 - c1 * a2);
					y1 = d1 - c1 * y2;

					x[k][j][i][0] = y1;
					x[k][j][i][Nmomentum - 1] = y2;
				}
				else {

					double r = 1.0 / (1.0 - a[k][j][i][1] * c[k][j][i][0]);
					rightPart[k][j][i][0] = r * (rightPart[k][j][i][0] - rightPart[k][j][i][1] * c[k][j][i][0]);
					c[k][j][i][0] = r * (u - c[k][j][i][0] * c[k][j][i][1]);

					double a1 = 1.0;
					double c1 = c[k][j][i][0];
					double d1 = rightPart[k][j][i][0];

					double a2 = a[k][j][i][Nmomentum - 1];
					double c2 = 1.0;
					double d2 = rightPart[k][j][i][Nmomentum - 1];

					y2 = (d2 - d1 * a2) / (c2 - c1 * a2);
					y1 = d1 - c1 * y2;

					x[k][j][i][0] = y1;
					x[k][j][i][Nmomentum - 1] = y2;
				}


				for (int l = 0; l < Nmomentum; ++l) {
					x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * y1 - c[k][j][i][l] * y2;
					if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
						printf("x = NaN in solver P, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
						exit(0);
					}
				}

				//delete[] d;
			}
		}
	}
}

void parallelThreeDiagonalSolverX(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum, int Nprocs, int rank, MPI_Comm& comm) {
    //double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);
    double**** parallelRightPart = create4dArray(2*Nprocs, Ny, Nz, Nmomentum);
    double**** parallelA = create4dArray(2*Nprocs, Ny, Nz, Nmomentum);
    double**** parallelB = create4dArray(2*Nprocs, Ny, Nz, Nmomentum);
    double**** parallelC = create4dArray(2*Nprocs, Ny, Nz, Nmomentum);
    double**** parallelX = create4dArray(2*Nprocs, Ny, Nz, Nmomentum);

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < Nmomentum; ++l) {

                double normRightPart = 0;

                for (int i = 0; i < Nx; ++i) {
                    normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
                }
                double norm[1];
                double temp[1];
                norm[0] = normRightPart;
                MPI_Allreduce(norm, temp, 1, MPI_DOUBLE, MPI_SUM, comm);
                normRightPart = temp[0];

                if (normRightPart <= 0) {
                    for (int i = 0; i < Nx; ++i) {
                        x[k][j][i][l] = 0;
                    }

                    continue;
                }

                //double u = a[k][j][0][l]/b[k][j][0][l];
                //double v = c[k][j][Nx - 1][l]/b[k][j][Nx - 1][l];
                for (int i = 0; i < Nx; ++i) {
                    a[k][j][i][l] /= b[k][j][i][l];
                    c[k][j][i][l] /= b[k][j][i][l];
                    rightPart[k][j][i][l] /= b[k][j][i][l];
                    b[k][j][i][l] = 1.0;
                }
                //double* d = new double[Nx];
                //d[0] = u;
                //d[1] = 0;

                for (int i = 2; i < Nx; ++i) {
                    double r = 1.0 / (1 - a[k][j][i][l] * c[k][j][i - 1][l]);
                    //d[i] = -r * a[k][j][i][l] * d[i - 1];
                    rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j][i - 1][l]);
                    a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j][i - 1][l];
                    if (i == Nx - 1) {
                        //a[k][j][i][l] += v * r;
                    }
                    c[k][j][i][l] = r * c[k][j][i][l];

                }

                for (int i = Nx - 3; i >= 1; i = i - 1) {
                    rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j][i + 1][l] * c[k][j][i][l];
                    a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j][i + 1][l];
                    c[k][j][i][l] = -c[k][j][i][l] * c[k][j][i + 1][l];
                }

                double r = 1.0 / (1.0 - a[k][j][1][l] * c[k][j][0][l]);
                rightPart[k][j][0][l] = r * (rightPart[k][j][0][l] - rightPart[k][j][1][l] * c[k][j][0][l]);
                //c[k][j][0][l] = r * (u - c[k][j][0][l] * c[k][j][1][l]);
                c[k][j][0][l] = - r * c[k][j][0][l] * c[k][j][1][l];


                double* outcoef = new double[6];
                outcoef[0] = a[k][j][0][l];
                outcoef[1] = c[k][j][0][l];
                outcoef[2] = rightPart[k][j][0][l];
                outcoef[3] = a[k][j][Nx-1][l];
                outcoef[4] = c[k][j][Nx-1][l];
                outcoef[5] = rightPart[k][j][Nx-1][l];

                double* incoef = new double[6*Nprocs];

                //MPI_Gather(outcoef, 6, MPI_DOUBLE, incoef, 6*Nprocs, MPI_DOUBLE, 0, comm);
                MPI_Gather(outcoef, 6, MPI_DOUBLE, incoef, 6, MPI_DOUBLE, 0, comm);
                if(rank == 0){
                    for(int m = 0; m < Nprocs; ++m){
                        parallelA[k][j][2*m][l] = incoef[6*m];
                        parallelC[k][j][2*m][l] = incoef[6*m+1];
                        parallelRightPart[k][j][2*m][l] = incoef[6*m+2];
                        parallelB[k][j][2*m][l] = 1.0;
                        parallelA[k][j][2*m+1][l] = incoef[6*m + 3];
                        parallelC[k][j][2*m+1][l] = incoef[6*m + 4];
                        parallelRightPart[k][j][2*m+1][l] = incoef[6*m + 5];
                        parallelB[k][j][2*m+1][l] = 1.0;
                    }
                }

                delete[] incoef;
                delete[] outcoef;
            }
        }
    }

    MPI_Barrier(comm);

    /*for(int m = 0; m < Nprocs; ++m){
        if(rank == m){
            for (int j = 0; j < Ny; ++j) {
                for (int k = 0; k < Nz; ++k) {
                    for (int l = 0; l < Nmomentum; ++l) {
                        for(int i = 0; i < Nx; ++i){
                            printf("%g %g %g      %g\n", a[k][j][i][l], b[k][j][i][l], c[k][j][i][l], rightPart[k][j][i][l]);
                        }
                    }
                }
            }
        }
        MPI_Barrier(comm);
    }*/
    MPI_Barrier(comm);

    if(rank == 0){
        sequentialThreeDiagonalSolverX(parallelX, parallelRightPart, parallelA, parallelB, parallelC, 2*Nprocs, Ny, Nz, Nmomentum);
        double* send = new double[2*Ny*Nz*Nmomentum];
        for(int m = 1; m < Nprocs; ++m){
            int count = 0;
            for(int k = 0; k < Nz; ++k){
                for(int j = 0; j < Ny; ++j){
                    for(int l = 0; l < Nmomentum; ++l){
                        send[count] = parallelX[k][j][2*m][l];
                        count = count + 1;
                        send[count] = parallelX[k][j][2*m+1][l];
                        count = count + 1;
                    }
                }
            }
            MPI_Send(send, 2*Ny*Nz*Nmomentum, MPI_DOUBLE, m, 0, comm);
        }
        delete[] send;
        for(int k = 0; k < Nz; ++k){
            for(int j = 0; j < Ny; ++j){
                for(int l = 0; l < Nmomentum; ++l){
                    x[k][j][0][l] = parallelX[k][j][0][l];
                    x[k][j][Nx-1][l] = parallelX[k][j][1][l];
                }
            }
        }
    } else {
        double* recv = new double[2*Ny*Nz*Nmomentum];
        MPI_Status status;
        MPI_Recv(recv, 2*Ny*Nz*Nmomentum, MPI_DOUBLE, 0, 0, comm, &status);
        int count = 0;
        for(int k = 0; k < Nz; ++k){
            for(int j = 0; j < Ny; ++j){
                for(int l = 0; l < Nmomentum; ++l){
                    x[k][j][0][l] = recv[count];
                    count = count + 1;
                    x[k][j][Nx - 1][l] = recv[count];
                    count = count + 1;
                }
            }
        }
        delete[] recv;
    }

    for(int k = 0; k < Nz; ++k){
        for(int j = 0; j < Ny; ++j){
            for(int l = 0; l < Nmomentum; ++l){
                for (int i = 1; i < Nx - 1; ++i) {
                    x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * x[k][j][0][l] - c[k][j][i][l] * x[k][j][Nx-1][l];
                    if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
                        printf("x = NaN in solver X, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }
                }

                //delete[] d;
            }
        }
    }

    delete4dArray(parallelA, 2*Nprocs, Ny, Nz, Nmomentum);
    delete4dArray(parallelB, 2*Nprocs, Ny, Nz, Nmomentum);
    delete4dArray(parallelC, 2*Nprocs, Ny, Nz, Nmomentum);
    delete4dArray(parallelRightPart, 2*Nprocs, Ny, Nz, Nmomentum);
    delete4dArray(parallelX, 2*Nprocs, Ny, Nz, Nmomentum);
}

void parallelThreeDiagonalSolverY(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum, int Nprocs, int rank, MPI_Comm& comm) {
    //double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);
    double**** parallelRightPart = create4dArray(Nx, 2*Nprocs, Nz, Nmomentum);
    double**** parallelA = create4dArray(Nx, 2*Nprocs, Nz, Nmomentum);
    double**** parallelB = create4dArray(Nx, 2*Nprocs, Nz, Nmomentum);
    double**** parallelC = create4dArray(Nx, 2*Nprocs, Nz, Nmomentum);
    double**** parallelX = create4dArray(Nx, 2*Nprocs, Nz, Nmomentum);

    for (int i = 0; i < Nx; ++i) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < Nmomentum; ++l) {

                double normRightPart = 0;

                for (int j = 0; j < Ny; ++j) {
                    normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
                }
                double norm[1];
                double temp[1];
                norm[0] = normRightPart;
                MPI_Allreduce(norm, temp, 1, MPI_DOUBLE, MPI_SUM, comm);
                normRightPart = temp[0];

                if (normRightPart <= 0) {
                    for (int j = 0; j < Ny; ++j) {
                        x[k][j][i][l] = 0;
                    }

                    continue;
                }

                //double u = a[k][j][0][l]/b[k][j][0][l];
                //double v = c[k][j][Nx - 1][l]/b[k][j][Nx - 1][l];
                for (int j = 0; j < Ny; ++j) {
                    a[k][j][i][l] /= b[k][j][i][l];
                    c[k][j][i][l] /= b[k][j][i][l];
                    rightPart[k][j][i][l] /= b[k][j][i][l];
                    b[k][j][i][l] = 1.0;
                }
                //double* d = new double[Nx];
                //d[0] = u;
                //d[1] = 0;

                for (int j = 2; j < Ny; ++j) {
                    double r = 1.0 / (1 - a[k][j][i][l] * c[k][j-1][i][l]);
                    //d[i] = -r * a[k][j][i][l] * d[i - 1];
                    rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j-1][i][l]);
                    a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j-1][i][l];
                    if (i == Nx - 1) {
                        //a[k][j][i][l] += v * r;
                    }
                    c[k][j][i][l] = r * c[k][j][i][l];

                }

                for (int j = Ny - 3; j >= 1; j = j - 1) {
                    rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j+1][i][l] * c[k][j][i][l];
                    a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j+1][i][l];
                    c[k][j][i][l] = -c[k][j][i][l] * c[k][j+1][i][l];
                }

                double r = 1.0 / (1.0 - a[k][1][i][l] * c[k][0][i][l]);
                rightPart[k][0][i][l] = r * (rightPart[k][0][i][l] - rightPart[k][1][i][l] * c[k][0][i][l]);
                //c[k][j][0][l] = r * (u - c[k][j][0][l] * c[k][j][1][l]);
                c[k][0][i][l] = - r * c[k][0][i][l] * c[k][1][i][l];


                double* outcoef = new double[6];
                outcoef[0] = a[k][0][i][l];
                outcoef[1] = c[k][0][i][l];
                outcoef[2] = rightPart[k][0][i][l];
                outcoef[3] = a[k][Ny-1][i][l];
                outcoef[4] = c[k][Ny-1][i][l];
                outcoef[5] = rightPart[k][Ny-1][i][l];

                double* incoef = new double[6*Nprocs];

                //MPI_Gather(outcoef, 6, MPI_DOUBLE, incoef, 6*Nprocs, MPI_DOUBLE, 0, comm);
                MPI_Gather(outcoef, 6, MPI_DOUBLE, incoef, 6, MPI_DOUBLE, 0, comm);
                if(rank == 0){
                    for(int m = 0; m < Nprocs; ++m){
                        parallelA[k][2*m][i][l] = incoef[6*m];
                        parallelC[k][2*m][i][l] = incoef[6*m+1];
                        parallelRightPart[k][2*m][i][l] = incoef[6*m+2];
                        parallelB[k][2*m][i][l] = 1.0;
                        parallelA[k][2*m+1][i][l] = incoef[6*m + 3];
                        parallelC[k][2*m+1][i][l] = incoef[6*m + 4];
                        parallelRightPart[k][2*m+1][i][l] = incoef[6*m + 5];
                        parallelB[k][2*m+1][i][l] = 1.0;
                    }
                }

                delete[] incoef;
                delete[] outcoef;
            }
        }
    }

    MPI_Barrier(comm);

    /*for(int m = 0; m < Nprocs; ++m){
        if(rank == m){
            for (int i = 0; i < Nx; ++i) {
                for (int k = 0; k < Nz; ++k) {
                    for (int l = 0; l < Nmomentum; ++l) {
                        for(int j = 0; j < Ny; ++j){
                            printf("%g %g %g      %g\n", a[k][j][i][l], b[k][j][i][l], c[k][j][i][l], rightPart[k][j][i][l]);
                        }
                    }
                }
            }
        }
        MPI_Barrier(comm);
    }*/
    MPI_Barrier(comm);

    if(rank == 0){
        sequentialThreeDiagonalSolverY(parallelX, parallelRightPart, parallelA, parallelB, parallelC, Nx, 2*Nprocs, Nz, Nmomentum);
        double* send = new double[2*Nx*Nz*Nmomentum];
        for(int m = 1; m < Nprocs; ++m){
            int count = 0;
            for(int k = 0; k < Nz; ++k){
                for(int i = 0; i < Nx; ++i){
                    for(int l = 0; l < Nmomentum; ++l){
                        send[count] = parallelX[k][2*m][i][l];
                        count = count + 1;
                        send[count] = parallelX[k][2*m+1][i][l];
                        count = count + 1;
                    }
                }
            }
            MPI_Send(send, 2*Nx*Nz*Nmomentum, MPI_DOUBLE, m, 0, comm);
        }
        delete[] send;
        for(int k = 0; k < Nz; ++k){
            for(int i = 0; i < Nx; ++i){
                for(int l = 0; l < Nmomentum; ++l){
                    x[k][0][i][l] = parallelX[k][0][i][l];
                    x[k][Ny-1][i][l] = parallelX[k][1][i][l];
                }
            }
        }
    } else {
        double* recv = new double[2*Nx*Nz*Nmomentum];
        MPI_Status status;
        MPI_Recv(recv, 2*Nx*Nz*Nmomentum, MPI_DOUBLE, 0, 0, comm, &status);
        int count = 0;
        for(int k = 0; k < Nz; ++k){
            for(int i = 0; i < Nx; ++i){
                for(int l = 0; l < Nmomentum; ++l){
                    x[k][0][i][l] = recv[count];
                    count = count + 1;
                    x[k][Ny-1][i][l] = recv[count];
                    count = count + 1;
                }
            }
        }
        delete[] recv;
    }

    for(int k = 0; k < Nz; ++k){
        for(int i = 0; i < Nx; ++i){
            for(int l = 0; l < Nmomentum; ++l){
                for (int j = 1; j < Ny - 1; ++j) {
                    x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * x[k][0][i][l] - c[k][j][i][l] * x[k][Ny-1][i][l];
                    if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
                        printf("x = NaN in solver X, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }
                }

                //delete[] d;
            }
        }
    }

    delete4dArray(parallelA, Nx, 2*Nprocs, Nz, Nmomentum);
    delete4dArray(parallelB, Nx, 2*Nprocs, Nz, Nmomentum);
    delete4dArray(parallelC, Nx, 2*Nprocs, Nz, Nmomentum);
    delete4dArray(parallelRightPart, Nx, 2*Nprocs, Nz, Nmomentum);
    delete4dArray(parallelX, Nx, 2*Nprocs, Nz, Nmomentum);
}

void parallelThreeDiagonalSolverZ(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum, int Nprocs, int rank, MPI_Comm& comm) {
    //double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);
    double**** parallelRightPart = create4dArray(Nx, Ny, 2*Nprocs, Nmomentum);
    double**** parallelA = create4dArray(Nx, Ny, 2*Nprocs, Nmomentum);
    double**** parallelB = create4dArray(Nx, Ny, 2*Nprocs, Nmomentum);
    double**** parallelC = create4dArray(Nx, Ny, 2*Nprocs, Nmomentum);
    double**** parallelX = create4dArray(Nx, Ny, 2*Nprocs, Nmomentum);

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            for (int l = 0; l < Nmomentum; ++l) {

                double normRightPart = 0;

                for (int k = 0; k < Nz; ++k) {
                    normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
                }
                double norm[1];
                double temp[1];
                norm[0] = normRightPart;
                MPI_Allreduce(norm, temp, 1, MPI_DOUBLE, MPI_SUM, comm);
                normRightPart = temp[0];

                if (normRightPart <= 0) {
                    for (int k = 0; k < Nz; ++k) {
                        x[k][j][i][l] = 0;
                    }

                    continue;
                }

                //double u = a[k][j][0][l]/b[k][j][0][l];
                //double v = c[k][j][Nx - 1][l]/b[k][j][Nx - 1][l];
                for (int k = 0; k < Nz; ++k) {
                    a[k][j][i][l] /= b[k][j][i][l];
                    c[k][j][i][l] /= b[k][j][i][l];
                    rightPart[k][j][i][l] /= b[k][j][i][l];
                    b[k][j][i][l] = 1.0;
                }
                //double* d = new double[Nx];
                //d[0] = u;
                //d[1] = 0;

                for (int k = 2; k < Nz; ++k) {
                    double r = 1.0 / (1 - a[k][j][i][l] * c[k-1][j][i][l]);
                    //d[i] = -r * a[k][j][i][l] * d[i - 1];
                    rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k-1][j][i][l]);
                    a[k][j][i][l] = -r * a[k][j][i][l] * a[k-1][j][i][l];
                    if (k == Nz - 1) {
                        //a[k][j][i][l] += v * r;
                    }
                    c[k][j][i][l] = r * c[k][j][i][l];

                }

                for (int k = Nz - 3; k >= 1; k = k - 1) {
                    rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k+1][j][i][l] * c[k][j][i][l];
                    a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k+1][j][i][l];
                    c[k][j][i][l] = -c[k][j][i][l] * c[k+1][j][i][l];
                }

                double r = 1.0 / (1.0 - a[1][j][i][l] * c[0][j][i][l]);
                rightPart[0][j][i][l] = r * (rightPart[0][j][i][l] - rightPart[1][j][i][l] * c[0][j][i][l]);
                //c[k][j][0][l] = r * (u - c[k][j][0][l] * c[k][j][1][l]);
                c[0][j][i][l] = - r * c[0][j][i][l] * c[1][j][i][l];


                double* outcoef = new double[6];
                outcoef[0] = a[0][j][i][l];
                outcoef[1] = c[0][j][i][l];
                outcoef[2] = rightPart[0][j][i][l];
                outcoef[3] = a[Nz - 1][j][i][l];
                outcoef[4] = c[Nz - 1][j][i][l];
                outcoef[5] = rightPart[Nz - 1][j][i][l];

                double* incoef = new double[6*Nprocs];

                //MPI_Gather(outcoef, 6, MPI_DOUBLE, incoef, 6*Nprocs, MPI_DOUBLE, 0, comm);
                MPI_Gather(outcoef, 6, MPI_DOUBLE, incoef, 6, MPI_DOUBLE, 0, comm);
                if(rank == 0){
                    for(int m = 0; m < Nprocs; ++m){
                        parallelA[2*m][j][i][l] = incoef[6*m];
                        parallelC[2*m][j][i][l] = incoef[6*m+1];
                        parallelRightPart[2*m][j][i][l] = incoef[6*m+2];
                        parallelB[2*m][j][i][l] = 1.0;
                        parallelA[2*m+1][j][i][l] = incoef[6*m + 3];
                        parallelC[2*m+1][j][i][l] = incoef[6*m + 4];
                        parallelRightPart[2*m+1][j][i][l] = incoef[6*m + 5];
                        parallelB[2*m+1][j][i][l] = 1.0;
                    }
                }

                delete[] incoef;
                delete[] outcoef;
            }
        }
    }

    MPI_Barrier(comm);

    /*for(int m = 0; m < Nprocs; ++m){
        if(rank == m){
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    for (int l = 0; l < Nmomentum; ++l) {
                        for(int k = 0; k < Nz; ++k){
                            printf("%g %g %g      %g\n", a[k][j][i][l], b[k][j][i][l], c[k][j][i][l], rightPart[k][j][i][l]);
                        }
                    }
                }
            }
        }
        MPI_Barrier(comm);
    }*/
    MPI_Barrier(comm);

    if(rank == 0){
        sequentialThreeDiagonalSolverZ(parallelX, parallelRightPart, parallelA, parallelB, parallelC, Nx, Ny, 2*Nprocs, Nmomentum);
        double* send = new double[2*Ny*Nx*Nmomentum];
        for(int m = 1; m < Nprocs; ++m){
            int count = 0;
            for(int i = 0; i < Nx; ++i){
                for(int j = 0; j < Ny; ++j){
                    for(int l = 0; l < Nmomentum; ++l){
                        send[count] = parallelX[2*m][j][i][l];
                        count = count + 1;
                        send[count] = parallelX[2*m + 1][j][i][l];
                        count = count + 1;
                    }
                }
            }
            MPI_Send(send, 2*Ny*Nx*Nmomentum, MPI_DOUBLE, m, 0, comm);
        }
        delete[] send;
        for(int i = 0; i < Nx; ++i){
            for(int j = 0; j < Ny; ++j){
                for(int l = 0; l < Nmomentum; ++l){
                    x[0][j][i][l] = parallelX[0][j][i][l];
                    x[Nz-1][j][i][l] = parallelX[1][j][i][l];
                }
            }
        }
    } else {
        double* recv = new double[2*Ny*Nx*Nmomentum];
        MPI_Status status;
        MPI_Recv(recv, 2*Ny*Nx*Nmomentum, MPI_DOUBLE, 0, 0, comm, &status);
        int count = 0;
        for(int i = 0; i < Nx; ++i){
            for(int j = 0; j < Ny; ++j){
                for(int l = 0; l < Nmomentum; ++l){
                    x[0][j][i][l] = recv[count];
                    count = count + 1;
                    x[Nz-1][j][i][l] = recv[count];
                    count = count + 1;
                }
            }
        }
        delete[] recv;
    }

    for(int i = 0; i < Nx; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int l = 0; l < Nmomentum; ++l){
                for (int k = 1; k < Nz - 1; ++k) {
                    x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * x[0][j][i][l] - c[k][j][i][l] * x[Nz-1][j][i][l];
                    if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
                        printf("x = NaN in solver X, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }
                }

                //delete[] d;
            }
        }
    }

    delete4dArray(parallelA, Nx, Ny, 2*Nprocs, Nmomentum);
    delete4dArray(parallelB, Nx, Ny, 2*Nprocs, Nmomentum);
    delete4dArray(parallelC, Nx, Ny, 2*Nprocs, Nmomentum);
    delete4dArray(parallelRightPart, Nx, Ny, 2*Nprocs, Nmomentum);
    delete4dArray(parallelX, Nx, Ny, 2*Nprocs, Nmomentum);
}

void testSequentialThreeDiagonalSolverX() {
	double**** F = create4dArray(9, 1, 1, 1);
	double**** rightPart = create4dArray(9, 1, 1, 1);
	double**** ax = create4dArray(9, 1, 1, 1);
	double**** bx = create4dArray(9, 1, 1, 1);
	double**** cx = create4dArray(9, 1, 1, 1);

    ax[0][0][0][0] = 0.0;
	bx[0][0][0][0] = 2.0;
	cx[0][0][0][0] = 1.0;

	ax[0][0][1][0] = 2.0;
	bx[0][0][1][0] = 2.0;
	cx[0][0][1][0] = -1.0;

	ax[0][0][2][0] = 1.0;
	bx[0][0][2][0] = 1.0;
	cx[0][0][2][0] = -2.0;

	ax[0][0][3][0] = 1.0;
	bx[0][0][3][0] = 2.0;
	cx[0][0][3][0] = 3.0;

	ax[0][0][4][0] = -1.0;
	bx[0][0][4][0] = 1.0;
	cx[0][0][4][0] = 2.0;

	ax[0][0][5][0] = 0.0;
	bx[0][0][5][0] = 3.0;
	cx[0][0][5][0] = 1.0;

	ax[0][0][6][0] = 2.0;
	bx[0][0][6][0] = 1.0;
	cx[0][0][6][0] = -3.0;

	ax[0][0][7][0] = 1.0;
	bx[0][0][7][0] = 2.0;
	cx[0][0][7][0] = 2.0;

	ax[0][0][8][0] = 2.0;
	bx[0][0][8][0] = 1.0;
    cx[0][0][8][0] = 0.0;

	double v1 = 0;
	double v2 = 0;

	for (int i = 0; i < 9; ++i) {
		rightPart[0][0][i][0] = i + 1;
	}

	for (int i = 0; i < 9; ++i) {
		for (int j = 0; j < 9; ++j) {
			if (i == j) {
				printf("%g ", bx[0][0][i][0]);
			}
			else if (i == j - 1) {
				printf("%g ", cx[0][0][i][0]);
			}
			else if (i == j + 1) {
				printf("%g ", ax[0][0][i][0]);
			}
			else {
				printf("0 ");
			}
		}
		printf("       %g\n", rightPart[0][0][i][0]);
	}

	sequentialThreeDiagonalSolverX(F, rightPart, ax, bx, cx, 9, 1, 1, 1);

	for (int i = 0; i < 9; ++i) {
		printf("%g\n", F[0][0][i][0]);
	}
}

void testSequentialThreeDiagonalSolverZ() {
	double**** F = create4dArray(1, 1, 9, 1);
	double**** rightPart = create4dArray(1, 1, 9, 1);
	double**** ax = create4dArray(1, 1, 9, 1);
	double**** bx = create4dArray(1, 1, 9, 1);
	double**** cx = create4dArray(1, 1, 9, 1);

	ax[0][0][0][0] = 0.0;
	bx[0][0][0][0] = 2.0;
	cx[0][0][0][0] = 1.0;

	ax[1][0][0][0] = 2.0;
	bx[1][0][0][0] = 2.0;
	cx[1][0][0][0] = -1.0;

	ax[2][0][0][0] = 1.0;
	bx[2][0][0][0] = 1.0;
	cx[2][0][0][0] = -2.0;

	ax[3][0][0][0] = 1.0;
	bx[3][0][0][0] = 2.0;
	cx[3][0][0][0] = 3.0;

	ax[4][0][0][0] = -1.0;
	bx[4][0][0][0] = 1.0;
	cx[4][0][0][0] = 2.0;

	ax[5][0][0][0] = 0.0;
	bx[5][0][0][0] = 3.0;
	cx[5][0][0][0] = 1.0;

	ax[6][0][0][0] = 2.0;
	bx[6][0][0][0] = 1.0;
	cx[6][0][0][0] = -3.0;

	ax[7][0][0][0] = 1.0;
	bx[7][0][0][0] = 2.0;
	cx[7][0][0][0] = 2.0;

	ax[8][0][0][0] = 2.0;
	bx[8][0][0][0] = 1.0;
	cx[8][0][0][0] = 0.0;

	double v1 = 0;
	double v2 = 0;

	for (int i = 0; i < 9; ++i) {
		rightPart[i][0][0][0] = i + 1;
	}

	for (int i = 0; i < 9; ++i) {
		for (int j = 0; j < 9; ++j) {
			if (i == j) {
				printf("%g ", bx[i][0][0][0]);
			}
			else if (i == j - 1) {
				printf("%g ", cx[i][0][0][0]);
			}
			else if (i == j + 1) {
				printf("%g ", ax[i][0][0][0]);
			}
			else {
				printf("0 ");
			}
		}
		printf("       %g\n", rightPart[i][0][0][0]);
	}

	sequentialThreeDiagonalSolverZ(F, rightPart, ax, bx, cx, 1, 1, 9, 1);

	for (int i = 0; i < 9; ++i) {
		printf("%g\n", F[i][0][0][0]);
	}
}

void testSequentialThreeDiagonalSolverP() {
	double**** F = create4dArray(1, 1, 1, 9);
	double**** rightPart = create4dArray(1, 1, 1, 9);
	double**** ax = create4dArray(1, 1, 1, 9);
	double**** bx = create4dArray(1, 1, 1, 9);
	double**** cx = create4dArray(1, 1, 1, 9);

	ax[0][0][0][0] = 0.0;
	bx[0][0][0][0] = 2.0;
	cx[0][0][0][0] = 1.0;

	ax[0][0][0][1] = 2.0;
	bx[0][0][0][1] = 2.0;
	cx[0][0][0][1] = -1.0;

	ax[0][0][0][2] = 1.0;
	bx[0][0][0][2] = 1.0;
	cx[0][0][0][2] = -2.0;

	ax[0][0][0][3] = 1.0;
	bx[0][0][0][3] = 2.0;
	cx[0][0][0][3] = 3.0;

	ax[0][0][0][4] = -1.0;
	bx[0][0][0][4] = 1.0;
	cx[0][0][0][4] = 2.0;

	ax[0][0][0][5] = 0.0;
	bx[0][0][0][5] = 3.0;
	cx[0][0][0][5] = 1.0;

	ax[0][0][0][6] = 2.0;
	bx[0][0][0][6] = 1.0;
	cx[0][0][0][6] = -3.0;

	ax[0][0][0][7] = 1.0;
	bx[0][0][0][7] = 2.0;
	cx[0][0][0][7] = 2.0;

	ax[0][0][0][8] = 2.0;
	bx[0][0][0][8] = 1.0;
	cx[0][0][0][8] = 0.0;

	double v1 = 0;
	double v2 = 0;

	for (int i = 0; i < 9; ++i) {
		rightPart[0][0][0][i] = i + 1;
	}

	for (int i = 0; i < 9; ++i) {
		for (int j = 0; j < 9; ++j) {
			if (i == j) {
				printf("%g ", bx[0][0][0][i]);
			}
			else if (i == j - 1) {
				printf("%g ", cx[0][0][0][i]);
			}
			else if (i == j + 1) {
				printf("%g ", ax[0][0][0][i]);
			}
			else {
				printf("0 ");
			}
		}
		printf("       %g\n", rightPart[0][0][0][i]);
	}

	sequentialThreeDiagonalSolverP(F, rightPart, ax, bx, cx, 1, 1, 1, 9);

	for (int i = 0; i < 9; ++i) {
		printf("%g\n", F[0][0][0][i]);
	}

	/*for (int i = 0; i < Nmomentum; ++i) {
		rightPart[0][0][0][i] = 0;
		for (int j = 0; j < Nmomentum; ++j) {
			if (i == j) {
				rightPart[0][0][0][i] += bx[0][0][0][i] * F[0][0][0][i];
			}
			else if (i == j - 1) {
				rightPart[0][0][0][i] += cx[0][0][0][i] * F[0][0][0][j];
			}
			else if (i == j + 1) {
				rightPart[0][0][0][i] += ax[0][0][0][i] * F[0][0][0][j];
			}
		}
		printf("       %g\n", rightPart[0][0][0][i]);
	}*/
}

const double D = 50;

void advanceDiffusionStepExplicit(double**** F, double**** rightPart, double* xgrid, double* ygrid, double* zgrid, double* pgrid, double*** u, double dt, int Nx, int Ny, int Nz, int Nmomentum) {
	double dx = xgrid[1] - xgrid[0];

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				double divu = 0;
				if (i == 0) {
					divu = (u[k][j][1] - u[k][j][0]) / dx;
				}
				else {
					divu = (u[k][j][i] - u[k][j][i-1]) / dx;
				}
				for (int l = 0; l < Nmomentum; ++l) {
					rightPart[k][j][i][l] = F[k][j][i][l];
					if (l == 0) {
						//rightPart[k][j][i][l] = rightPart[k][j][i][l] + dt * (pgrid[l] / 3.0) * divu * (F[k][j][i][l] - 0) / (pgrid[1] - pgrid[0]);
					}
					else {
						rightPart[k][j][i][l] = rightPart[k][j][i][l] + dt * (pgrid[l] / 3.0) * divu * (F[k][j][i][l] - F[k][j][i][l - 1]) / (pgrid[l] - pgrid[l - 1]);
					}
					if ((rightPart[k][j][i][l] != rightPart[k][j][i][l]) || (0 * rightPart[k][j][i][l] != 0 * rightPart[k][j][i][l])) {
						printf("rightPart = NaN k = %d, j = %d, i = %d, l = %d\n", k, j, i, l);
						exit(0);
					}
				}
			}
		}
	}

	double**** a = create4dArray(Nx, Ny, Nz, Nmomentum);
	double**** b = create4dArray(Nx, Ny, Nz, Nmomentum);
	double**** c = create4dArray(Nx, Ny, Nz, Nmomentum);

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				for (int l = 0; l < Nmomentum; ++l) {
					if (i == 0) {
						a[k][j][i][l] = 0.0;
						b[k][j][i][l] = 1.0 + u[k][j][i]*dx/D;
						c[k][j][i][l] = -1.0;
						rightPart[k][j][i][l] = 0;
					}
					else if (i == Nx - 1) {
						a[k][j][i][l] = -1.0;
						b[k][j][i][l] = 1.0;
						c[k][j][i][l] = 0;
						rightPart[k][j][i][l] = 0.0;
					}
					else {
						a[k][j][i][l] = -dt*D / (dx * dx) - dt * u[k][j][i - 1] / dx;
						b[k][j][i][l] = 1 + 2 * dt*D / (dx * dx) + dt * u[k][j][i] / dx;
						c[k][j][i][l] = -dt*D / (dx * dx);
					}
				}
			}
		}
	}

	sequentialThreeDiagonalSolverX(F, rightPart, a, b, c, Nx, Ny, Nz, Nmomentum);

	delete4dArray(a, Nx, Ny, Nz, Nmomentum);
	delete4dArray(b, Nx, Ny, Nz, Nmomentum);
	delete4dArray(c, Nx, Ny, Nz, Nmomentum);
}

void advanceDiffusionStepImplicit(double**** F, double**** rightPart, double* xgrid, double* ygrid, double* zgrid, double* pgrid, double*** u, double dt, int Nx, int Ny, int Nz, int Nmomentum) {
	double dx = xgrid[1] - xgrid[0];

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				double divu = 0;
				if (i == 0) {
					divu = (u[k][j][1] - u[k][j][0]) / dx;
				}
				else {
					divu = (u[k][j][i] - u[k][j][i - 1]) / dx;
				}
				for (int l = 0; l < Nmomentum; ++l) {
					rightPart[k][j][i][l] = F[k][j][i][l];
					if (l == 0) {
						//rightPart[k][j][i][l] = rightPart[k][j][i][l] + dt * (pgrid[l] / 3.0) * divu * (F[k][j][i][l] - 0) / (pgrid[1] - pgrid[0]);
					}
					else {
						//rightPart[k][j][i][l] = rightPart[k][j][i][l] + dt * (pgrid[l] / 3.0) * divu * (F[k][j][i][l] - F[k][j][i][l - 1]) / (pgrid[l] - pgrid[l - 1]);
					}
					if ((rightPart[k][j][i][l] != rightPart[k][j][i][l]) || (0 * rightPart[k][j][i][l] != 0 * rightPart[k][j][i][l])) {
						printf("rightPart = NaN k = %d, j = %d, i = %d, l = %d\n", k, j, i, l);
						exit(0);
					}
				}
			}
		}
	}

	double**** a = create4dArray(Nx, Ny, Nz, Nmomentum);
	double**** b = create4dArray(Nx, Ny, Nz, Nmomentum);
	double**** c = create4dArray(Nx, Ny, Nz, Nmomentum);

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				for (int l = 0; l < Nmomentum; ++l) {
					if (i == 0) {
						a[k][j][i][l] = 0.0;
						b[k][j][i][l] = 1.0 + u[k][j][i]*dx/D;
						c[k][j][i][l] = -1.0;
						rightPart[k][j][i][l] = 0;
					}
					else if (i == Nx - 1) {
						a[k][j][i][l] = -1.0;
						b[k][j][i][l] = 1.0;
						c[k][j][i][l] = 0;
						rightPart[k][j][i][l] = 0.0;
					}
					else {
						a[k][j][i][l] = -dt * D / (dx * dx) - dt * u[k][j][i - 1] / dx;
						b[k][j][i][l] = 1 + 2 * dt * D / (dx * dx) + dt * u[k][j][i] / dx;
						c[k][j][i][l] = -dt * D / (dx * dx);
					}
				}
			}
		}
	}

	sequentialThreeDiagonalSolverX(F, rightPart, a, b, c, Nx, Ny, Nz, Nmomentum);

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				double divu = 0;
				if (i == 0) {
					divu = (u[k][j][1] - u[k][j][0]) / dx;
				}
				else {
					divu = (u[k][j][i] - u[k][j][i-1]) / dx;
				}
				for (int l = 0; l < Nmomentum; ++l) {
					rightPart[k][j][i][l] = F[k][j][i][l];
					if (l == 0) {
						double dp = pgrid[1] - pgrid[0];
						a[k][j][i][l] = 0;
						b[k][j][i][l] = 1.0 - dt * pgrid[0] * divu / (3*dp);
						c[k][j][i][l] = 0;
					}
					/*else if (l == Nmomentum - 1) {
						a[k][j][i][l] = 0.0;
						b[k][j][i][l] = 1.0;
						c[k][j][i][l] = 0;
						rightPart[k][j][i][l] = 0.0;
					}*/
					else {
						double dp = pgrid[l] - pgrid[l - 1];
						a[k][j][i][l] = dt * pgrid[l] * divu / (3*dp);
						b[k][j][i][l] = 1 - dt * pgrid[l] * divu / (3*dp);
						c[k][j][i][l] = 0;
					}
				}
			}
		}
	}

	sequentialThreeDiagonalSolverP(F, rightPart, a, b, c, Nx, Ny, Nz, Nmomentum);

	delete4dArray(a, Nx, Ny, Nz, Nmomentum);
	delete4dArray(b, Nx, Ny, Nz, Nmomentum);
	delete4dArray(c, Nx, Ny, Nz, Nmomentum);
}

void advanceDiffusionStepGMRES(double**** F, double**** rightPart, std::vector<MatrixElement>**** matrix, LargeVectorBasis* GMRESbasis, double* xgrid, double* ygrid, double* zgrid, double* pgrid, double*** u, double dt, int Nx, int Ny, int Nz, int Nmomentum) {
	GMRESbasis->clear();
	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				for (int l = 0; l < Nmomentum; ++l) {
					matrix[k][j][i][l].clear();
				}
			}
		}
	}

	double dx = xgrid[1] - xgrid[0];

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				double divu = 0;
				if (i == 0) {
					divu = (u[k][j][1] - u[k][j][0]) / dx;
				}
				else {
					divu = (u[k][j][i] - u[k][j][i - 1]) / dx;
				}
				for (int l = 0; l < Nmomentum; ++l) {
					rightPart[k][j][i][l] = F[k][j][i][l];
					if (l == 0) {
						//rightPart[k][j][i][l] = rightPart[k][j][i][l] + dt * (pgrid[l] / 3.0) * divu * (F[k][j][i][l] - 0) / (pgrid[1] - pgrid[0]);
					}
					else {
						//rightPart[k][j][i][l] = rightPart[k][j][i][l] + dt * (pgrid[l] / 3.0) * divu * (F[k][j][i][l] - F[k][j][i][l - 1]) / (pgrid[l] - pgrid[l - 1]);
					}
					if ((rightPart[k][j][i][l] != rightPart[k][j][i][l]) || (0 * rightPart[k][j][i][l] != 0 * rightPart[k][j][i][l])) {
						printf("rightPart = NaN k = %d, j = %d, i = %d, l = %d\n", k, j, i, l);
						exit(0);
					}
				}
			}
		}
	}

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				double divu = 0;
				if (i == 0) {
					divu = (u[k][j][1] - u[k][j][0]) / dx;
				}
				else {
					divu = (u[k][j][i] - u[k][j][i - 1]) / dx;
				}
				for (int l = 0; l < Nmomentum; ++l) {
					if (i == 0) {
						//matrix[k][j][i][l].push_back(MatrixElement(1.0, k, j, i, l));
						matrix[k][j][i][l].push_back(MatrixElement(1.0 + u[k][j][i]*dx/D, k, j, i, l));
						matrix[k][j][i][l].push_back(MatrixElement(-1.0, k, j, i + 1, l));
						//a[k][j][i][l] = 0;
						//b[k][j][i][l] = 1.0;
						//c[k][j][i][l] = 0.0;
						rightPart[k][j][i][l] = 0;
					}
					else if (i == Nx - 1) {
						matrix[k][j][i][l].push_back(MatrixElement(1.0, k, j, i, l));
						matrix[k][j][i][l].push_back(MatrixElement(-1.0, k, j, i - 1, l));
						//a[k][j][i][l] = 0.0;
						//b[k][j][i][l] = 1.0;
						//c[k][j][i][l] = 0;
						rightPart[k][j][i][l] = 0.0;
					}
					else {
						matrix[k][j][i][l].push_back(MatrixElement(1 + 2 * dt * D / (dx * dx) + dt * u[k][j][i] / dx, k, j, i, l));
						matrix[k][j][i][l].push_back(MatrixElement(-dt * D / (dx * dx) - dt * u[k][j][i-1] / dx, k, j, i - 1, l));
						matrix[k][j][i][l].push_back(MatrixElement(-dt * D / (dx * dx) , k, j, i + 1, l));
						//a[k][j][i][l] = -dt * D / (dx * dx) - dt * u[k][j][i - 1] / dx;
						//b[k][j][i][l] = 1 + 2 * dt * D / (dx * dx) + dt * u[k][j][i] / dx;
						//c[k][j][i][l] = -dt * D / (dx * dx);

						if (l == 0) {
							double dp = pgrid[1] - pgrid[0];
							matrix[k][j][i][l].push_back(MatrixElement(-dt * pgrid[0] * divu / (3*dp), k, j, i, l));
							//a[k][j][i][l] = 0;
							//b[k][j][i][l] = 1.0 - dt * pgrid[0] * divu / dp;
							//c[k][j][i][l] = 0;
						}
						/*else if (l == Nmomentum - 1) {
							matrix[k][j][i][l].clear();
							matrix[k][j][i][l].push_back(MatrixElement(1.0, k, j, i, l));
							//a[k][j][i][l] = 0;
							//b[k][j][i][l] = 1.0;
							//c[k][j][i][l] = 0;
							rightPart[k][j][i][l] = 0.0;
						}*/
						else {
							double dp = pgrid[l] - pgrid[l - 1];
							matrix[k][j][i][l].push_back(MatrixElement(-dt * pgrid[l] * divu / (3.0*dp), k, j, i, l));
							matrix[k][j][i][l].push_back(MatrixElement(dt * pgrid[l] * divu / (3.0*dp), k, j, i, l - 1));
							//a[k][j][i][l] = dt * pgrid[l] * divu / dp;
							//b[k][j][i][l] = 1 - dt * pgrid[l] * divu / dp;
							//c[k][j][i][l] = 0;
						}
					}
				}
			}
		}
	}

	generalizedMinimalResidualMethod(matrix, rightPart, F, Nz, Ny, Nx, Nmomentum, 1E-5, 100, 1, GMRESbasis);

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				for (int l = 0; l < Nmomentum; ++l) {
					if (F[k][j][i][l] < 0) {
						F[k][j][i][l] = 0;
					}
				}
			}
		}
	}
}

void testDiffusion(){
    int Nx = 100;
    int Ny = 1;
    int Nz = 1;
    int Nmomentum = 50;
    double**** F = create4dArray(Nx, Ny, Nz, Nmomentum);
    double**** F1 = create4dArray(Nx, Ny, Nz, Nmomentum);
    double**** F2 = create4dArray(Nx, Ny, Nz, Nmomentum);
    double**** tempF = create4dArray(Nx, Ny, Nz, Nmomentum);
    double**** rightPart = create4dArray(Nx, Ny, Nz, Nmomentum);

    double*** u = create3dArray(Nx, Ny, Nz);

    LargeVectorBasis* GMRESbasis = new LargeVectorBasis(10, Nz, Ny, Nx, Nmomentum);

    std::vector<MatrixElement>**** matrix = new std::vector<MatrixElement>***[Nz];
    for (int k = 0; k < Nz; ++k) {
        matrix[k] = new std::vector<MatrixElement>**[Ny];
        for (int j = 0; j < Ny; ++j) {
            matrix[k][j] = new std::vector<MatrixElement>*[Nx];
            for (int i = 0; i < Nx; ++i) {
                matrix[k][j][i] = new std::vector<MatrixElement>[Nx];
            }
        }
    }

    double* xgrid = new double[Nx];
    for (int i = 0; i < Nx; ++i) {
        xgrid[i] = i;
    }
    double* ygrid = new double[Ny];
    for (int j = 0; j < Ny; ++j) {
        ygrid[j] = j;
    }
    double* zgrid = new double[Nz];
    for (int k = 0; k < Nz; ++k) {
        zgrid[k] = k;
    }
    double pmin = 1.0;
    double pmax = 1E3;
    double factor = 1.0;
    if(Nmomentum > 1){
        factor = pow(pmax / pmin, 1.0 / (Nmomentum - 1.0));
    }
    double* pgrid = new double[Nmomentum];
    pgrid[0] = pmin;
    for (int l = 1; l < Nmomentum; ++l) {
        pgrid[l] = pgrid[l - 1] * factor;
    }

    FILE* xfile = fopen("./output/xgrid.dat", "w");
    for (int i = 0; i < Nx; ++i) {
        fprintf(xfile, "%g\n", xgrid[i]);
    }
    fclose(xfile);
    FILE* yfile = fopen("./output/ygrid.dat", "w");
    for (int i = 0; i < Ny; ++i) {
        fprintf(yfile, "%g\n", ygrid[i]);
    }
    fclose(yfile);
    FILE* zfile = fopen("./output/zgrid.dat", "w");
    for (int i = 0; i < Nz; ++i) {
        fprintf(zfile, "%g\n", zgrid[i]);
    }
    fclose(zfile);
    FILE* pfile = fopen("./output/pgrid.dat", "w");
    for (int i = 0; i < Nmomentum; ++i) {
        fprintf(pfile, "%g\n", pgrid[i]);
    }
    fclose(pfile);


    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                if (i < Nx / 2 - 5) {
                    u[k][j][i] = 1.0;
                }
                else if (i >= Nx / 2 + 5) {
                    u[k][j][i] = 0.25;
                }
                else {
                    u[k][j][i] = 1.0 - 0.75*(i - Nx/2 + 5)/10.0;
                }
            }
        }
    }

    double dx = xgrid[1] - xgrid[0];
    double maxDivU = (1.0 - 0.25) / (xgrid[Nx / 2] - xgrid[Nx / 2 - 1]);
    double maxU = 1.0;

    double advectiveDt = 0.5 * dx / maxU;
    double accelerationDt = 0.5 * (pgrid[1] - pgrid[0]) / (pgrid[1] * maxDivU);


    int Nt = 1000000;
    double dt = 0.1;
    int writeN = 200;
    int writeCount = 0;
    for (int m = 0; m < Nt; ++m) {
        printf("timestep %d\n", m);
        //injection
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                F[k][j][Nx / 2][1] += 1.0*dt;
                F1[k][j][Nx / 2][1] += 1.0 * dt;
                F2[k][j][Nx / 2][1] += 1.0 * dt;
            }
        }

        advanceDiffusionStepExplicit(F, rightPart, xgrid, ygrid, zgrid, pgrid, u, dt, Nx, Ny, Nz, Nmomentum);
        advanceDiffusionStepImplicit(F1, rightPart, xgrid, ygrid, zgrid, pgrid, u, dt, Nx, Ny, Nz, Nmomentum);
        advanceDiffusionStepGMRES(F2, rightPart, matrix, GMRESbasis, xgrid, ygrid, zgrid, pgrid, u, dt, Nx, Ny, Nz, Nmomentum);

        if (m % writeN == 0) {
            std::string fileName = std::string("./output/F") + std::to_string(writeCount) + std::string(".dat");
            FILE* outFile = fopen(fileName.c_str(), "w");
            for (int k = 0; k < Nz; ++k) {
                for (int j = 0; j < Ny; ++j) {
                    for (int i = 0; i < Nx; ++i) {
                        for (int l = 0; l < Nmomentum; ++l) {
                            fprintf(outFile, "%g %g %g\n", F[k][j][i][l], F1[k][j][i][l], F2[k][j][i][l]);
                        }
                    }
                }
            }
            fclose(outFile);
            writeCount = writeCount + 1;
        }
    }
}

void testParallelThreeDiagonalSolverX( int argc, char *argv[] ){
    int size;
    int rank;
    MPI_Init(&argc, &argv);

    int MPI_Nx = 2;
    int MPI_Ny = 1;
    int MPI_Nz = 1;
    const int MPI_dim = 3;

    int dims[MPI_dim];
    int period[MPI_dim];
    period[0] = 1;
    period[1] = 1;
    period[2] = 1;
    dims[0] = MPI_Nx;
    dims[1] = MPI_Ny;
    dims[2] = MPI_Nz;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, MPI_dim, dims, period, 0, &comm);

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    printf("size = %d\n", size);
    printf("rank = %d\n", rank);

    MPI_Barrier(comm);

    int Nxg = 9;
    int Nyg = 1;
    int Nzg = 1;
    int Np = 1;

    double**** F = create4dArray(Nxg, 1, 1, 1);
    double**** rightPart = create4dArray(Nxg, 1, 1, 1);
    double**** ax = create4dArray(Nxg, 1, 1, 1);
    double**** bx = create4dArray(Nxg, 1, 1, 1);
    double**** cx = create4dArray(Nxg, 1, 1, 1);

    ax[0][0][0][0] = 0.0;
    bx[0][0][0][0] = 2.0;
    cx[0][0][0][0] = 1.0;

    ax[0][0][1][0] = 2.0;
    bx[0][0][1][0] = 2.0;
    cx[0][0][1][0] = -1.0;

    ax[0][0][2][0] = 1.0;
    bx[0][0][2][0] = 1.0;
    cx[0][0][2][0] = -2.0;

    ax[0][0][3][0] = 1.0;
    bx[0][0][3][0] = 2.0;
    cx[0][0][3][0] = 3.0;

    ax[0][0][4][0] = -1.0;
    bx[0][0][4][0] = 1.0;
    cx[0][0][4][0] = 2.0;

    ax[0][0][5][0] = 0.0;
    bx[0][0][5][0] = 3.0;
    cx[0][0][5][0] = 1.0;

    ax[0][0][6][0] = 2.0;
    bx[0][0][6][0] = 1.0;
    cx[0][0][6][0] = -3.0;

    ax[0][0][7][0] = 1.0;
    bx[0][0][7][0] = 2.0;
    cx[0][0][7][0] = 2.0;

    ax[0][0][8][0] = 2.0;
    bx[0][0][8][0] = 1.0;
    cx[0][0][8][0] = 0.0;

    double v1 = 0;
    double v2 = 0;

    for (int i = 0; i < Nxg; ++i) {
        rightPart[0][0][i][0] = i + 1;
    }

    if(rank == 0){
    for (int i = 0; i < Nxg; ++i) {
        for (int j = 0; j < Nxg; ++j) {
            if (i == j) {
                printf("%g ", bx[0][0][i][0]);
            }
            else if (i == j - 1) {
                printf("%g ", cx[0][0][i][0]);
            }
            else if (i == j + 1) {
                printf("%g ", ax[0][0][i][0]);
            }
            else {
                printf("0 ");
            }
        }
        printf("       %g\n", rightPart[0][0][i][0]);
    }
    }

    int Nxl0 = Nxg/MPI_Nx;
    int Nxl = Nxl0;
    int Nxmod = Nxg%MPI_Nx;
    if(rank == MPI_Nx - 1){
        Nxl = Nxl + Nxmod;
    }

    double**** axl = create4dArray(Nxl, Nyg, Nzg, Np);
    double**** bxl = create4dArray(Nxl, Nyg, Nzg, Np);
    double**** cxl = create4dArray(Nxl, Nyg, Nzg, Np);
    double**** rightPartl = create4dArray(Nxl, Nyg, Nzg, Np);
    double**** Fl = create4dArray(Nxl, Nyg, Nzg, Np);

    for(int k = 0; k < Nzg; ++k){
        for(int j = 0; j < Nyg; ++j){
            for(int i = 0; i < Nxl; ++i){
                for(int l = 0; l < Np; ++l){
                    axl[k][j][i][l] = ax[k][j][rank*Nxl0 + i][l];
                    bxl[k][j][i][l] = bx[k][j][rank*Nxl0 + i][l];
                    cxl[k][j][i][l] = cx[k][j][rank*Nxl0 + i][l];
                    rightPartl[k][j][i][l] = rightPart[k][j][rank*Nxl0 + i][l];
                }
            }
        }
    }

    /*for(int m = 0; m < MPI_Nx; ++m){
        if(rank == m){
            for (int i = 0; i < Nxl; ++i) {
                if(i == 0){
                    printf("%g ", axl[0][0][i][0]);
                }
                for (int j = 0; j < Nxl; ++j) {
                    if (i == j) {
                        printf("%g ", bxl[0][0][i][0]);
                    }
                    else if (i == j - 1) {
                        printf("%g ", cxl[0][0][i][0]);
                    }
                    else if (i == j + 1) {
                        printf("%g ", axl[0][0][i][0]);
                    }
                    else {
                        printf("0 ");
                    }
                }
                if(i == Nxl - 1){
                    printf("%g ", cxl[0][0][i][0]);
                }
                printf("       %g\n", rightPartl[0][0][i][0]);
            }
        }
        MPI_Barrier(comm);
    }*/

    MPI_Barrier(comm);

    parallelThreeDiagonalSolverX(Fl, rightPartl, axl, bxl, cxl, Nxl, Nyg, Nzg, Np, MPI_Nx, rank, comm);

    MPI_Barrier(comm);

    for(int m = 0; m < MPI_Nx; ++m){
        if(rank == m){
            for (int i = 0; i < Nxl; ++i) {
                printf("%g\n", Fl[0][0][i][0]);
            }
        }
        MPI_Barrier(comm);
    }

    delete4dArray(axl, Nxl, Nyg, Nzg, Np);
    delete4dArray(bxl, Nxl, Nyg, Nzg, Np);
    delete4dArray(cxl, Nxl, Nyg, Nzg, Np);
    delete4dArray(rightPartl, Nxl, Nyg, Nzg, Np);
    delete4dArray(Fl, Nxl, Nyg, Nzg, Np);

    delete4dArray(ax, Nxg, Nyg, Nzg, Np);
    delete4dArray(bx, Nxg, Nyg, Nzg, Np);
    delete4dArray(cx, Nxg, Nyg, Nzg, Np);
    delete4dArray(rightPart, Nxg, Nyg, Nzg, Np);
    delete4dArray(F, Nxg, Nyg, Nzg, Np);

    printf("before finalize\n");

    MPI_Finalize();
}

void testParallelThreeDiagonalSolverY( int argc, char *argv[] ){
    int size;
    int rank;
    MPI_Init(&argc, &argv);

    int MPI_Nx = 1;
    int MPI_Ny = 2;
    int MPI_Nz = 1;
    const int MPI_dim = 3;

    int dims[MPI_dim];
    int period[MPI_dim];
    period[0] = 1;
    period[1] = 1;
    period[2] = 1;
    dims[0] = MPI_Nx;
    dims[1] = MPI_Ny;
    dims[2] = MPI_Nz;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, MPI_dim, dims, period, 0, &comm);

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    printf("size = %d\n", size);
    printf("rank = %d\n", rank);

    MPI_Barrier(comm);

    int Nxg = 1;
    int Nyg = 9;
    int Nzg = 1;
    int Np = 1;

    double**** F = create4dArray(Nxg, Nyg, Nzg, 1);
    double**** rightPart = create4dArray(Nxg, Nyg, Nzg, 1);
    double**** ax = create4dArray(Nxg, Nyg, Nzg, 1);
    double**** bx = create4dArray(Nxg, Nyg, Nzg, 1);
    double**** cx = create4dArray(Nxg, Nyg, Nzg, 1);

    ax[0][0][0][0] = 0.0;
    bx[0][0][0][0] = 2.0;
    cx[0][0][0][0] = 1.0;

    ax[0][1][0][0] = 2.0;
    bx[0][1][0][0] = 2.0;
    cx[0][1][0][0] = -1.0;

    ax[0][2][0][0] = 1.0;
    bx[0][2][0][0] = 1.0;
    cx[0][2][0][0] = -2.0;

    ax[0][3][0][0] = 1.0;
    bx[0][3][0][0] = 2.0;
    cx[0][3][0][0] = 3.0;

    ax[0][4][0][0] = -1.0;
    bx[0][4][0][0] = 1.0;
    cx[0][4][0][0] = 2.0;

    ax[0][5][0][0] = 0.0;
    bx[0][5][0][0] = 3.0;
    cx[0][5][0][0] = 1.0;

    ax[0][6][0][0] = 2.0;
    bx[0][6][0][0] = 1.0;
    cx[0][6][0][0] = -3.0;

    ax[0][7][0][0] = 1.0;
    bx[0][7][0][0] = 2.0;
    cx[0][7][0][0] = 2.0;

    ax[0][8][0][0] = 2.0;
    bx[0][8][0][0] = 1.0;
    cx[0][8][0][0] = 0.0;

    double v1 = 0;
    double v2 = 0;

    for (int i = 0; i < Nyg; ++i) {
        rightPart[0][i][0][0] = i + 1;
    }

    if(rank == 0){
    for (int i = 0; i < Nyg; ++i) {
        for (int j = 0; j < Nyg; ++j) {
            if (i == j) {
                printf("%g ", bx[0][i][0][0]);
            }
            else if (i == j - 1) {
                printf("%g ", cx[0][i][0][0]);
            }
            else if (i == j + 1) {
                printf("%g ", ax[0][i][0][0]);
            }
            else {
                printf("0 ");
            }
        }
        printf("       %g\n", rightPart[0][i][0][0]);
    }
    }

    int Nyl0 = Nyg/MPI_Ny;
    int Nyl = Nyl0;
    int Nymod = Nyg%MPI_Ny;
    if(rank == MPI_Ny - 1){
        Nyl = Nyl + Nymod;
    }

    double**** axl = create4dArray(Nxg, Nyl, Nzg, Np);
    double**** bxl = create4dArray(Nxg, Nyl, Nzg, Np);
    double**** cxl = create4dArray(Nxg, Nyl, Nzg, Np);
    double**** rightPartl = create4dArray(Nxg, Nyl, Nzg, Np);
    double**** Fl = create4dArray(Nxg, Nyl, Nzg, Np);

    for(int k = 0; k < Nzg; ++k){
        for(int j = 0; j < Nyl; ++j){
            for(int i = 0; i < Nxg; ++i){
                for(int l = 0; l < Np; ++l){
                    axl[k][j][i][l] = ax[k][rank*Nyl0 + j][i][l];
                    bxl[k][j][i][l] = bx[k][rank*Nyl0 + j][i][l];
                    cxl[k][j][i][l] = cx[k][rank*Nyl0 + j][i][l];
                    rightPartl[k][j][i][l] = rightPart[k][rank*Nyl0 + j][i][l];
                }
            }
        }
    }

    for(int m = 0; m < MPI_Ny; ++m){
        if(rank == m){
            for (int i = 0; i < Nyl; ++i) {
                if(i == 0){
                    printf("%g ", axl[0][i][0][0]);
                }
                for (int j = 0; j < Nyl; ++j) {
                    if (i == j) {
                        printf("%g ", bxl[0][i][0][0]);
                    }
                    else if (i == j - 1) {
                        printf("%g ", cxl[0][i][0][0]);
                    }
                    else if (i == j + 1) {
                        printf("%g ", axl[0][i][0][0]);
                    }
                    else {
                        printf("0 ");
                    }
                }
                if(i == Nyl - 1){
                    printf("%g ", cxl[0][i][0][0]);
                }
                printf("       %g\n", rightPartl[0][i][0][0]);
            }
        }
        MPI_Barrier(comm);
    }

    MPI_Barrier(comm);

    parallelThreeDiagonalSolverY(Fl, rightPartl, axl, bxl, cxl, Nxg, Nyl, Nzg, Np, MPI_Ny, rank, comm);

    MPI_Barrier(comm);

    for(int m = 0; m < MPI_Ny; ++m){
        if(rank == m){
            for (int i = 0; i < Nyl; ++i) {
                printf("%g\n", Fl[0][i][0][0]);
            }
        }
        MPI_Barrier(comm);
    }

    delete4dArray(axl, Nxg, Nyl, Nzg, Np);
    delete4dArray(bxl, Nxg, Nyl, Nzg, Np);
    delete4dArray(cxl, Nxg, Nyl, Nzg, Np);
    delete4dArray(rightPartl, Nxg, Nyl, Nzg, Np);
    delete4dArray(Fl, Nxg, Nyl, Nzg, Np);

    delete4dArray(ax, Nxg, Nyg, Nzg, Np);
    delete4dArray(bx, Nxg, Nyg, Nzg, Np);
    delete4dArray(cx, Nxg, Nyg, Nzg, Np);
    delete4dArray(rightPart, Nxg, Nyg, Nzg, Np);
    delete4dArray(F, Nxg, Nyg, Nzg, Np);

    printf("before finalize\n");

    MPI_Finalize();
}

void testParallelThreeDiagonalSolverZ( int argc, char *argv[] ){
    int size;
    int rank;
    MPI_Init(&argc, &argv);

    int MPI_Nx = 1;
    int MPI_Ny = 1;
    int MPI_Nz = 2;
    const int MPI_dim = 3;

    int dims[MPI_dim];
    int period[MPI_dim];
    period[0] = 1;
    period[1] = 1;
    period[2] = 1;
    dims[0] = MPI_Nx;
    dims[1] = MPI_Ny;
    dims[2] = MPI_Nz;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, MPI_dim, dims, period, 0, &comm);

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    printf("size = %d\n", size);
    printf("rank = %d\n", rank);

    MPI_Barrier(comm);

    int Nxg = 1;
    int Nyg = 1;
    int Nzg = 9;
    int Np = 1;

    double**** F = create4dArray(Nxg, Nyg, Nzg, 1);
    double**** rightPart = create4dArray(Nxg, Nyg, Nzg, 1);
    double**** ax = create4dArray(Nxg, Nyg, Nzg, 1);
    double**** bx = create4dArray(Nxg, Nyg, Nzg, 1);
    double**** cx = create4dArray(Nxg, Nyg, Nzg, 1);

    ax[0][0][0][0] = 0.0;
    bx[0][0][0][0] = 2.0;
    cx[0][0][0][0] = 1.0;

    ax[1][0][0][0] = 2.0;
    bx[1][0][0][0] = 2.0;
    cx[1][0][0][0] = -1.0;

    ax[2][0][0][0] = 1.0;
    bx[2][0][0][0] = 1.0;
    cx[2][0][0][0] = -2.0;

    ax[3][0][0][0] = 1.0;
    bx[3][0][0][0] = 2.0;
    cx[3][0][0][0] = 3.0;

    ax[4][0][0][0] = -1.0;
    bx[4][0][0][0] = 1.0;
    cx[4][0][0][0] = 2.0;

    ax[5][0][0][0] = 0.0;
    bx[5][0][0][0] = 3.0;
    cx[5][0][0][0] = 1.0;

    ax[6][0][0][0] = 2.0;
    bx[6][0][0][0] = 1.0;
    cx[6][0][0][0] = -3.0;

    ax[7][0][0][0] = 1.0;
    bx[7][0][0][0] = 2.0;
    cx[7][0][0][0] = 2.0;

    ax[8][0][0][0] = 2.0;
    bx[8][0][0][0] = 1.0;
    cx[8][0][0][0] = 0.0;

    double v1 = 0;
    double v2 = 0;

    for (int i = 0; i < Nzg; ++i) {
        rightPart[i][0][0][0] = i + 1;
    }

    if(rank == 0){
    for (int i = 0; i < Nzg; ++i) {
        for (int j = 0; j < Nzg; ++j) {
            if (i == j) {
                printf("%g ", bx[i][0][0][0]);
            }
            else if (i == j - 1) {
                printf("%g ", cx[i][0][0][0]);
            }
            else if (i == j + 1) {
                printf("%g ", ax[i][0][0][0]);
            }
            else {
                printf("0 ");
            }
        }
        printf("       %g\n", rightPart[i][0][0][0]);
    }
    }

    int Nzl0 = Nzg/MPI_Nz;
    int Nzl = Nzl0;
    int Nzmod = Nzg%MPI_Nz;
    if(rank == MPI_Nz - 1){
        Nzl = Nzl + Nzmod;
    }

    double**** axl = create4dArray(Nxg, Nyg, Nzl, Np);
    double**** bxl = create4dArray(Nxg, Nyg, Nzl, Np);
    double**** cxl = create4dArray(Nxg, Nyg, Nzl, Np);
    double**** rightPartl = create4dArray(Nxg, Nyg, Nzl, Np);
    double**** Fl = create4dArray(Nxg, Nyg, Nzl, Np);

    for(int k = 0; k < Nzl; ++k){
        for(int j = 0; j < Nyg; ++j){
            for(int i = 0; i < Nxg; ++i){
                for(int l = 0; l < Np; ++l){
                    axl[k][j][i][l] = ax[rank*Nzl0 + k][j][i][l];
                    bxl[k][j][i][l] = bx[rank*Nzl0 + k][j][i][l];
                    cxl[k][j][i][l] = cx[rank*Nzl0 + k][j][i][l];
                    rightPartl[k][j][i][l] = rightPart[rank*Nzl0 + k][j][i][l];
                }
            }
        }
    }

    /*for(int m = 0; m < MPI_Nx; ++m){
        if(rank == m){
            for (int i = 0; i < Nxl; ++i) {
                if(i == 0){
                    printf("%g ", axl[0][0][i][0]);
                }
                for (int j = 0; j < Nxl; ++j) {
                    if (i == j) {
                        printf("%g ", bxl[0][0][i][0]);
                    }
                    else if (i == j - 1) {
                        printf("%g ", cxl[0][0][i][0]);
                    }
                    else if (i == j + 1) {
                        printf("%g ", axl[0][0][i][0]);
                    }
                    else {
                        printf("0 ");
                    }
                }
                if(i == Nxl - 1){
                    printf("%g ", cxl[0][0][i][0]);
                }
                printf("       %g\n", rightPartl[0][0][i][0]);
            }
        }
        MPI_Barrier(comm);
    }*/

    MPI_Barrier(comm);

    parallelThreeDiagonalSolverZ(Fl, rightPartl, axl, bxl, cxl, Nxg, Nyg, Nzl, Np, MPI_Nz, rank, comm);

    MPI_Barrier(comm);

    for(int m = 0; m < MPI_Nz; ++m){
        if(rank == m){
            for (int i = 0; i < Nzl; ++i) {
                printf("%g\n", Fl[i][0][0][0]);
            }
        }
        MPI_Barrier(comm);
    }

    delete4dArray(axl, Nxg, Nyg, Nzl, Np);
    delete4dArray(bxl, Nxg, Nyg, Nzl, Np);
    delete4dArray(cxl, Nxg, Nyg, Nzl, Np);
    delete4dArray(rightPartl, Nxg, Nyg, Nzl, Np);
    delete4dArray(Fl, Nxg, Nyg, Nzl, Np);

    delete4dArray(ax, Nxg, Nyg, Nzg, Np);
    delete4dArray(bx, Nxg, Nyg, Nzg, Np);
    delete4dArray(cx, Nxg, Nyg, Nzg, Np);
    delete4dArray(rightPart, Nxg, Nyg, Nzg, Np);
    delete4dArray(F, Nxg, Nyg, Nzg, Np);

    printf("before finalize\n");

    MPI_Finalize();
}

int main( int argc, char *argv[] )
{

    //testSequentialThreeDiagonalSolverX();
	//testSequentialThreeDiagonalSolverZ();
	//testSequentialThreeDiagonalSolverP();

    //testParallelThreeDiagonalSolverX(argc, argv);
    //testParallelThreeDiagonalSolverY(argc, argv);
    testParallelThreeDiagonalSolverZ(argc, argv);

    //testDiffusion();

}
