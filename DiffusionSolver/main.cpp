// DiffusionSolver.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include "stdio.h"
#include "math.h"

#include "./Math/largeVectorBasis.h"
#include "./Math/matrixElement.h"
#include "./Math/specialmath.h"
#include "./Math/util.h"

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

void testSequentialThreeDiagonalSolverX() {
	double**** F = create4dArray(9, 1, 1, 1);
	double**** rightPart = create4dArray(9, 1, 1, 1);
	double**** ax = create4dArray(9, 1, 1, 1);
	double**** bx = create4dArray(9, 1, 1, 1);
	double**** cx = create4dArray(9, 1, 1, 1);

	ax[0][0][0][0] = 1.0;
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
	cx[0][0][8][0] = 1.0;

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

const double D = 5;

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

int main()
{
	int Nx = 100;
	int Ny = 1;
	int Nz = 1;
	int Nmomentum = 100;
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
	double relp = (pgrid[1] - pgrid[0]) / pgrid[1];

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
				if (i < Nx / 2) {
					u[k][j][i] = 1.0;
				}
				else if (i == Nx / 2) {
					u[k][j][i] = 0.7;
				}
				else {
					u[k][j][i] = 0.25;
				}
			}
		}
	}

	//testSequentialThreeDiagonalSolverX();
	//testSequentialThreeDiagonalSolverZ();
	//testSequentialThreeDiagonalSolverP();

	double dx = xgrid[1] - xgrid[0];
	double maxDivU = (1.0 - 0.25) / (xgrid[Nx / 2] - xgrid[Nx / 2 - 1]);
	double maxU = 1.0;
	
	double advectiveDt = 0.5 * dx / maxU;
	double accelerationDt = 0.5 * (pgrid[1] - pgrid[0]) / (pgrid[1] * maxDivU);


	int Nt = 1000000;
	double dt = 0.02;
	int writeN = 2500;
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