// DiffusionSolver.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <string>
#include "stdio.h"
#include "math.h"

const int Nx = 100;
const int Ny = 1;
const int Nz = 1;
const int Nmomentum = 100;

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

void sequentialThreeDiagonalSolverX(double**** x, double**** rightPart, double**** a, double**** b, double**** c) {
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
				}

				//delete[] d;
			}
		}
	}
}

void sequentialThreeDiagonalSolverY(double**** x, double**** rightPart, double**** a, double**** b, double**** c) {
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
				}

				//delete[] d;
			}
		}
	}
}

void sequentialThreeDiagonalSolverZ(double**** x, double**** rightPart, double**** a, double**** b, double**** c) {
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
				}

				//delete[] d;
			}
		}
	}
}

void sequentialThreeDiagonalSolverP(double**** x, double**** rightPart, double**** a, double**** b, double**** c) {
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
				}
				//double* d = new double[Nx];
				//d[0] = u;
				//d[1] = 0;

				for (int l = 2; l < Nmomentum; ++l) {
					double r = 1.0 / (1 - a[k][j][i][l] * c[k][j][i - 1][l]);
					//d[i] = -r * a[k][j][i][l] * d[i - 1];
					rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j][i - 1][l]);
					a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j][i - 1][l];
					if (i == Nx - 1) {
						a[k][j][i][l] += v * r;
					}
					c[k][j][i][l] = r * c[k][j][i][l];

				}

				for (int l = Nmomentum - 3; l >= 1; l = l - 1) {
					rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j][i + 1][l] * c[k][j][i][l];
					a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j][i + 1][l];
					c[k][j][i][l] = -c[k][j][i][l] * c[k][j][i + 1][l];
				}

				double r = 1.0 / (1.0 - a[k][j][i][1] * c[k][j][i][0]);
				rightPart[k][j][i][0] = r * (rightPart[k][j][i][0] - rightPart[k][j][i][1] * c[k][j][i][0]);
				c[k][j][i][0] = r * (u - c[k][j][i][0] * c[k][j][i][1]);

				double a1 = 1.0;
				double c1 = c[k][j][i][0];
				double d1 = rightPart[k][j][i][0];

				double a2 = a[k][j][i][Nmomentum - 1];
				double c2 = 1.0;
				double d2 = rightPart[k][j][i][Nmomentum - 1];

				double y2 = (d2 - d1 * a2) / (c2 - c1 * a2);
				double y1 = d1 - c1 * y2;

				x[k][j][i][0] = y1;
				x[k][j][i][Nmomentum - 1] = y2;

				for (int l = 0; l < Nmomentum; ++l) {
					x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * y1 - c[k][j][i][l] * y2;
				}

				//delete[] d;
			}
		}
	}
}

void testSequentialThreeDiagonalSolver() {
	double**** F = create4Darray(Nx, Ny, Nz, Nmomentum);
	double**** rightPart = create4Darray(Nx, Ny, Nz, Nmomentum);
	double**** ax = create4Darray(Nx, Ny, Nz, Nmomentum);
	double**** bx = create4Darray(Nx, Ny, Nz, Nmomentum);
	double**** cx = create4Darray(Nx, Ny, Nz, Nmomentum);

	ax[0][0][0][0] = 1.0;
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
	cx[8][0][0][0] = 1.0;

	double v1 = 0;
	double v2 = 0;

	for (int i = 0; i < Nz; ++i) {
		rightPart[i][0][0][0] = i + 1;
	}

	for (int i = 0; i < Nz; ++i) {
		for (int j = 0; j < Nz; ++j) {
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

	sequentialThreeDiagonalSolverZ(F, rightPart, ax, bx, cx);

	for (int i = 0; i < Nz; ++i) {
		printf("%g\n", F[i][0][0][0]);
	}
}

void advanceDiffusionStep(double**** F, double**** rightPart, double* xgrid, double* ygrid, double* zgrid, double* pgrid, double*** u, double dt) {
	double dx = xgrid[1] - xgrid[0];

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				for (int l = 0; l < Nmomentum; ++l) {
					rightPart[k][j][i][l] = F[k][j][i][l];
				}
			}
		}
	}

	double**** a = create4Darray(Nx, Ny, Nz, Nmomentum);
	double**** b = create4Darray(Nx, Ny, Nz, Nmomentum);
	double**** c = create4Darray(Nx, Ny, Nz, Nmomentum);

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				for (int l = 0; l < Nmomentum; ++l) {
					if (i == 0) {
						a[k][j][i][l] = 0;
						b[k][j][i][l] = 1.0;
						c[k][j][i][l] = 0;
					}
					else if (i == Nx - 1) {
						a[k][j][i][l] = 0;
						b[k][j][i][l] = 1.0;
						c[k][j][i][l] = 0;
					}
					else {
						a[k][j][i][l] = -dt / (dx * dx) - dt * u[k][j][i - 1] / dx;
						b[k][j][i][l] = 1 + 2 * dt / (dx * dx) + dt * u[k][j][i] / dx;
						c[k][j][i][l] = -dt / (dx * dx);
					}
				}
			}
		}
	}

	sequentialThreeDiagonalSolverX(F, rightPart, a, b, c);

	for (int k = 0; k < Nz; ++k) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				double divu = (u[i] - u[i - 1]) / dx;
				for (int l = 0; l < Nmomentum; ++l) {
					rightPart[k][j][i][l] = F[k][j][i][l];
					double dp = 0;
					if (l == 0) {
						double dp = pgrid[1] - pgrid[0];
						a[k][j][i][l] = 0;
						b[k][j][i][l] = 1.0 - dt * pgrid[0] * divu / dp;
						c[k][j][i][l] = dt * pgrid[0] * divu / dp;
					}
					else if (l == Nx - 1) {
						a[k][j][i][l] = 0;
						b[k][j][i][l] = 1.0;
						c[k][j][i][l] = 0;
					}
					else {
						dp = pgrid[l] - pgrid[l - 1];
						a[k][j][i][l] = -dt * pgrid[l] * divu / dp;
						b[k][j][i][l] = 1 + dt * pgrid[l] * divu / dp;
						c[k][j][i][l] = 0;
					}
				}
			}
		}
	}

	sequentialThreeDiagonalSolverP(F, rightPart, a, b, c);
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
		factor = pow(pmax / pmin, 1.0 / (Nmomentum - 1.0));
	}
	double* pgrid = new double[Nmomentum];
	pgrid[0] = pmin;
	for (int l = 1; l < Nmomentum; ++l) {
		pgrid[l] = pgrid[l - 1] * factor;
	}

	FILE* xfile = fopen("xgrid.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		fprintf(xfile, "%g\n", xgrid[i]);
	}
	fclose(xfile);
	FILE* yfile = fopen("ygrid.dat", "w");
	for (int i = 0; i < Ny; ++i) {
		fprintf(yfile, "%g\n", ygrid[i]);
	}
	fclose(yfile);
	FILE* zfile = fopen("zgrid.dat", "w");
	for (int i = 0; i < Nz; ++i) {
		fprintf(zfile, "%g\n", zgrid[i]);
	}
	fclose(zfile);
	FILE* pfile = fopen("pgrid.dat", "w");
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
				else {
					u[k][j][i] = 0.25;
				}
			}
		}
	}

	//testSequentialThreeDiagonalSolver();


	int Nt = 1000000;
	double dt = 0.1;
	int writeN = 100;
	int writeCount = 0;
	for (int m = 0; m < Nt; ++m) {
		printf("timestep %d\n", m);
		//injection
		for (int k = 0; k < Nz; ++k) {
			for (int j = 0; j < Ny; ++j) {
				F[k][j][Nx / 2][0] += 1.0*dt;
			}
		}

		advanceDiffusionStep(F, rightPart, xgrid, ygrid, zgrid, pgrid, u, dt);

		if (m % writeN == 0) {
			std::string fileName = std::string("F") + std::to_string(writeCount) + std::string(".dat");
			FILE* outFile = fopen(fileName.c_str(), "w");
			for (int k = 0; k < Nz; ++k) {
				for (int j = 0; j < Ny; ++j) {
					for (int i = 0; i < Nx; ++i) {
						for (int l = 0; l < Nmomentum; ++l) {
							fprintf(outFile, "%g\n", F[k][j][i][l]);
						}
					}
				}
			}
			fclose(outFile);
			writeCount = writeCount + 1;
		}
	}
}