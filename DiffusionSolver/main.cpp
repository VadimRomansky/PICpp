// DiffusionSolver.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "stdio.h"
#include "math.h"

const int Nx = 100;
const int Ny = 1;
const int Nz = 1;
const int Nmomentum = 200;

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
	double factor = pow(pmax / pmin, 1.0 / (Nmomentum - 1.0));
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

	int Nt = 1000000;
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
	}
}