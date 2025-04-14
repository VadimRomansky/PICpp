#ifndef _MATRIX_ELEMENT_H_
#define _MATRIX_ELEMENT_H_

class MatrixElement
{
public:
	double value;
	int i;
	int j;
	int k;
	int l;

	MatrixElement();
	MatrixElement(const double& v, int iv, int jv, int kv, int lv);

	bool equalsIndex(const MatrixElement& element);
};

int convert4dIndexTo1d(int i, int j, int k, int l, int Nx, int Ny, int Nz, int Nl);
int convert3dIndexTo1d(int i, int j, int k, int Nx, int Ny, int Nz);
int convert2dIndexTo1d(int i, int j, int Nx, int Ny);

void convert1dIndexTo4d(int index, int& i, int& j, int& k, int& l, int Nx, int Ny, int Nz, int Nl);
void convert1dIndexTo3d(int index, int& i, int& j, int& k, int Nx, int Ny, int Nz);
void convert1dIndexTo2d(int index, int& i, int& j, int Nx, int Ny);

#endif