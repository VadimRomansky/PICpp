#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "util.h"

double max(const double& a, const double& b) {
	if (a >= b) {
		return a;
	}
	return b;
}

double max3(const double& a, const double& b, const double& c) {
	if (a > b) {
		return max(a, c);
	}
	return max(b, c);
}

double min(const double& a, const double& b) {
	if (a >= b) {
		return b;
	}
	return a;
}

double min3(const double& a, const double& b, const double& c) {
	if (a > b) {
		return min(b, c);
	}
	return min(a, c);
}

double*** create3dArray(const int Nx, const int Ny, const int Nz, const double& value)
{
    double*** a = new double** [Nz];
    for (int i = 0; i < Nz; ++i) {
        a[i] = new double* [Ny];
        for (int j = 0; j < Ny; ++j) {
            a[i][j] = new double[Nx];
            for (int k = 0; k < Nx; ++k) {
                a[i][j][k] = value;
            }
        }
    }
    return a;
}

void delete3dArray(double*** a, const int Nx, const int Ny, const int Nz)
{
    for (int i = 0; i < Nz; ++i) {
        for (int j = 0; j < Ny; ++j) {
            delete[] a[i][j];
        }
        delete[] a[i];
    }
    delete[] a;
}

double**** create4dArray(const int Nx, const int Ny, const int Nz, const int Nmomentum, const double& value)
{
    double**** a = new double*** [Nz];
    for (int i = 0; i < Nz; ++i) {
        a[i] = new double** [Ny];
        for (int j = 0; j < Ny; ++j) {
            a[i][j] = new double* [Nx];
            for (int k = 0; k < Nx; ++k) {
                a[i][j][k] = new double[Nmomentum];
                for (int l = 0; l < Nmomentum; ++l) {
                    a[i][j][k][l] = value;
                }
            }
        }
    }
    return a;
}

void delete4dArray(double**** a, const int Nx, const int Ny, const int Nz, const int Nmomentum)
{
    for (int i = 0; i < Nz; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nx; ++k) {
                delete[] a[i][j][k];
            }
            delete[] a[i][j];
        }
        delete[] a[i];
    }
    delete[] a;
}


