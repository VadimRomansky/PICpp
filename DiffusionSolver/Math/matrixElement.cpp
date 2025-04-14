//#include <crtdbg.h>

#include <stdio.h>
#include <stdlib.h>

//#include "memory_debug.h"
//#include "../util.h"
#include "matrixElement.h"

MatrixElement::MatrixElement() {
	value = 0;
	i = 0;
	j = 0;
	k = 0;
	l = 0;
}

MatrixElement::MatrixElement(const double& v, int iv, int jv, int kv, int lv) {
	value = v;
	i = iv;
	j = jv;
	k = kv;
	l = lv;
}

bool MatrixElement::equalsIndex(const MatrixElement& element) {
	if (i != element.i) return false;
	if (j != element.j) return false;
	if (k != element.k) return false;
	if (l != element.l) return false;

	return true;
}

int convert4dIndexTo1d(int i, int j, int k, int l, int Nx, int Ny, int Nz, int Nl)
{
	if (i >= Nx) {
		printf("i >= Nx in convert4dIndex\n");
		exit(0);
	}
	if (j >= Ny) {
		printf("j >= Ny in convert4dIndex\n");
		exit(0);
	}
	if (k >= Nz) {
		printf("k >= Nz in convert4dIndex\n");
		exit(0);
	}
	if (l >= Nl) {
		printf("l >= Nl in convert4dIndex\n");
		exit(0);
	}

	return i * Ny * Nz * Nl + j * Nz * Nl + k * Nl + l;
}

int convert3dIndexTo1d(int i, int j, int k, int Nx, int Ny, int Nz)
{
	if (i >= Nx) {
		printf("i >= Nx in convert3dIndex\n");
		exit(0);
	}
	if (j >= Ny) {
		printf("j >= Ny in convert3dIndex\n");
		exit(0);
	}
	if (k >= Nz) {
		printf("k >= Nz in convert3dIndex\n");
		exit(0);
	}
	return i*Ny*Nz + j*Nz + k;
}

int convert2dIndexTo1d(int i, int j, int Nx, int Ny)
{
	if (i >= Nx) {
		printf("i >= Nx in convert2dIndex\n");
		exit(0);
	}
	if (j >= Ny) {
		printf("j >= Ny in convert2dIndex\n");
		exit(0);
	}
	return i*Ny + j;
}

void convert1dIndexTo4d(int index, int& i, int& j, int& k, int& l, int Nx, int Ny, int Nz, int Nl)
{
	if (index >= Nx * Ny * Nz * Nl) {
		printf("index >= N in convertIndexTo4d\n");
		exit(0);
	}

	i = index / (Ny * Nz * Nl);
	int mod = index % (Ny * Nz * Nl);

	j = mod / (Nz * Nl);

	int mod1 = mod % (Nz * Nl);

	k = mod1 / Nl;

	l = mod1 % Nl;
}

void convert1dIndexTo3d(int index, int& i, int& j, int& k, int Nx, int Ny, int Nz)
{
	if (index >= Nx * Ny * Nz) {
		printf("index >= N in convertIndexTo3d\n");
		exit(0);
	}

	i = index / (Ny * Nz);

	int mod = index % (Ny * Nz);

	j = mod / Nz;

	k = mod % Nz;
}

void convert1dIndexTo2d(int index, int& i, int& j, int Nx, int Ny)
{
	if (index >= Nx * Ny) {
		printf("index >= N in convertIndexTo2d\n");
		exit(0);
	}

	i = index / Ny;

	j = index % Ny;
}
