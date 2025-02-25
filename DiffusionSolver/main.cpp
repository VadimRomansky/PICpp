// DiffusionSolver.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "stdio.h"
#include "math.h"

const int Nx = 9;
const int Ny = 1;
const int Nz = 1;
const int Nmomentum = 1;

double*** create3Darray(int x, int y, int z) {
	double*** u = new double** [z];
	for (int i = 0; i < z; ++i) {
		u[i] = new double* [y];
		for (int j = 0; j < y; ++j) {
			u[i][j] = new double [x];
			for (int k = 0; k < x; ++k) {
				u[i][j][k] = 0;
			}
		}
	}
	return u;
}

double**** create4Darray(int x, int y, int z, int l) {
	double**** u = new double*** [z];
	for (int i = 0; i < z; ++i) {
		u[i] = new double** [y];
		for (int j = 0; j < y; ++j) {
			u[i][j] = new double* [x];
			for (int k = 0; k < x; ++k) {
				u[i][j][k] = new double[l];
				for (int m = 0; m < l; ++m) {
					u[i][j][k][m] = 0;
				}
			}
		}
	}
	return u;
}

void delete3Darray(double*** u, int x, int y, int z) {
	for (int i = 0; i < z; ++i) {
		for (int j = 0; j < y; ++j) {
			delete[] u[i][j];
		}
		delete[] u[i];
	}
	delete[] u;
}

void delete4Darray(double**** u, int x, int y, int z, int l) {
	for (int i = 0; i < z; ++i) {
		for (int j = 0; j < y; ++j) {
			for (int k = 0; k < x; ++k) {
				delete[] u[i][j][k];
			}
			delete[] u[i][j];
		}
		delete[] u[i];
	}
	delete[] u;
}

double evaluateDiffusionCoefficient(double p) {
	return 0.001 * p;
}

void advanceDiffusionStep(double**** F, double**** rightPart, double* xgrid, double* ygrid, double* zgrid, double* pgrid, double dt) {
	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				for (int l = 0; l < Nmomentum; ++l) {
					rightPart[k][j][i][l] = F[k][j][i][l];
				}
			}
		}
	}

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				for (int l = 0; l < Nmomentum; ++l) {

				}
			}
		}
	}
}

void sequentialThreeDiagonalSolverX(double**** x, double**** rightPart, double**** a, double**** b, double**** c, double u, double v) {
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

				u /= b[k][j][0][l];
				v /= b[k][j][Nx - 1][l];
				for (int i = 0; i < Nx; ++i) {
					a[k][j][i][l] /= b[k][j][i][l];
					c[k][j][i][l] /= b[k][j][i][l];
					rightPart[k][j][i][l] /= b[k][j][i][l];
					b[k][j][i][l] = 1.0;
				}
				double* d = new double[Nx];
				d[0] = u;
				d[1] = 0;

				for (int i = 2; i < Nx; ++i) {
					double r = 1.0 / (1 - a[k][j][i][l] * c[k][j][i - 1][l]);
					d[i] = -r * a[k][j][i][l] * d[i - 1];
					rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j][i - 1][l]);
					a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j][i - 1][l];
					c[k][j][i][l] = r * c[k][j][i][l];

				}

				a[k][j][Nx - 1][l] += v;

				for (int i = Nx - 3; i >= 1; i = i - 1) {
					rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j][i + 1][l] * c[k][j][i][l];
					c[k][j][i][l] = d[i] - c[k][j][i][l] * c[k][j][i + 1][l];
					a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j][i + 1][l];
				}

				double r = 1.0 / (1.0 - a[k][j][1][l] * c[k][j][0][l]);
				rightPart[k][j][0][l] = r * (rightPart[k][j][0][l] - rightPart[k][j][1][l] * c[k][j][0][l]);
				c[k][j][0][l] = r * (d[0] - c[k][j][0][l] * c[k][j][1][l]);

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
				}

				delete[] d;
			}
		}
	}
}

int main()
{
	double**** F = create4Darray(Nx, Ny, Nz, Nmomentum);
	double**** tempF = create4Darray(Nx, Ny, Nz, Nmomentum);
	double**** rightPart = create4Darray(Nx, Ny, Nz, Nmomentum);

	double*** u = create3Darray(Nx, Ny, Nz);
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
	double pmax = 1E7;
	double factor = 1.0;
	if(Nmomentum > 1){
		//factor = pow(pmax / pmin, 1.0 / (Nmomentum - 1.0));
	}
	double* pgrid = new double[Nmomentum];
	pgrid[0] = pmin;
	for (int l = 1; l < Nmomentum; ++l) {
		pgrid[l] = pgrid[l - 1] * factor;
	}

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				if (i < Nx / 2) {
					u[k][j][i] = 1.0;
				}
				else {
					u[k][j][i] = 0.25;
				}
			}
		}
	}

	double**** ax = create4Darray(Nx, Ny, Nz, Nmomentum);
	double**** bx = create4Darray(Nx, Ny, Nz, Nmomentum);
	double**** cx = create4Darray(Nx, Ny, Nz, Nmomentum);

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

	for (int i = 0; i < Nx; ++i) {
		rightPart[0][0][i][0] = i + 1;
	}

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nx; ++j) {
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

	sequentialThreeDiagonalSolverX(F, rightPart, ax, bx, cx, v1, v2);

	for (int i = 0; i < Nx; ++i) {
		printf("%g\n", F[0][0][i][0]);
	}


	/*int Nt = 1000000;
	double dt = 0.1;
	for (int m = 0; m < Nt; ++m) {
		printf("timestep %d\n", m);
		//injection
		for (int k = 0; k < Nz; ++k) {
			for (int j = 0; j < Ny; ++j) {
				F[k][j][Nx / 2][0] += 1.0*dt;
			}
		}

		advanceDiffusionStep(F, rightPart, xgrid, ygrid, zgrid, pgrid, dt);
	}*/
}