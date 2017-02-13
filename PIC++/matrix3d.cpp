#include "math.h"
#include "stdio.h"
#include "cmath"
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "matrix3d.h"
#include "vector3d.h"
#include "constants.h"

Matrix3d::Matrix3d() {
	//matrix = new double* [3];
	for (int i = 0; i < 3; ++i) {
		//matrix[i] = new double [3];
		for (int j = 0; j < 3; ++j) {
			if (i == j) {
				matrix[i][j] = 1;
			} else {
				matrix[i][j] = 0;
			}
		}
	}
}

Matrix3d::Matrix3d(const Matrix3d& m) {
	int i, j;
	//matrix = new double* [3];
	for (i = 0; i < 3; ++i) {
		//matrix[i] = new double [3];
		for (j = 0; j < 3; ++j) {
			matrix[i][j] = m.matrix[i][j];
		}
	}
}

Matrix3d::Matrix3d(double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33) {
	/*matrix = new double* [3];
	for (int i = 0; i < 3; ++i) {
		matrix[i] = new double [3];
	}*/
	matrix[0][0] = m11;
	matrix[0][1] = m12;
	matrix[0][2] = m13;
	matrix[1][0] = m21;
	matrix[1][1] = m22;
	matrix[1][2] = m23;
	matrix[2][0] = m31;
	matrix[2][1] = m32;
	matrix[2][2] = m33;
}

Matrix3d::~Matrix3d() {
	/*int i;
	for (i = 0; i < 3; ++i) {
		delete[] matrix[i];
	}
	delete[] matrix;*/
}

Matrix3d* Matrix3d::Inverse() {
	Matrix3d* inv = new Matrix3d();
	double det;
	det = determinant();
	if (fabs(det) > epsilon) {
		inv->matrix[0][0] = (1 / det) * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]);
		inv->matrix[0][1] = -(1 / det) * (matrix[0][1] * matrix[2][2] - matrix[0][2] * matrix[2][1]);
		inv->matrix[0][2] = (1 / det) * (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]);
		inv->matrix[1][0] = -(1 / det) * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]);
		inv->matrix[1][1] = (1 / det) * (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]);
		inv->matrix[1][2] = -(1 / det) * (matrix[0][0] * matrix[1][2] - matrix[1][0] * matrix[0][2]);
		inv->matrix[2][0] = (1 / det) * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
		inv->matrix[2][1] = -(1 / det) * (matrix[0][0] * matrix[2][1] - matrix[0][1] * matrix[2][0]);
		inv->matrix[2][2] = (1 / det) * (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]);
	} else {
		printf("determinant = 0\n");
	}
	return inv;
}

double Matrix3d::determinant() {
	double temp;
	temp = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
		- matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
		+ matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
	return temp;
}

double Matrix3d::getvalue(int i, int j) {
	return matrix[i][j];
}

void Matrix3d::setvalue(int i, int j, double value) {
	matrix[i][j] = value;
}

Matrix3d& Matrix3d::operator=(const Matrix3d& matr) {
	int i;
	int j;
	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			matrix[i][j] = matr.matrix[i][j];
		}
	}
	return *this;
}

Vector3d Matrix3d::operator*(const Vector3d& v) {
	double x = matrix[0][0] * v.x + matrix[0][1] * v.y + matrix[0][2] * v.z;
	double y = matrix[1][0] * v.x + matrix[1][1] * v.y + matrix[1][2] * v.z;
	double z = matrix[2][0] * v.x + matrix[2][1] * v.y + matrix[2][2] * v.z;
	return Vector3d(x, y, z);
}

Matrix3d Matrix3d::operator*(const Matrix3d& matr) {
	Matrix3d newMatrix;
	int i;
	int j;
	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			newMatrix.matrix[i][j] = 0;
			for (int k = 0; k < 3; ++k) {
				newMatrix.matrix[i][j] += matrix[i][k] * matr.matrix[k][j];
			}
		}
	}
	return newMatrix;
}

Matrix3d Matrix3d::operator*(const double& v) {
	Matrix3d newMatrix;
	int i;
	int j;
	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			newMatrix.matrix[i][j] = matrix[i][j] * v;
		}
	}
	return newMatrix;
}

Matrix3d Matrix3d::operator/(const double& v) {
	Matrix3d newMatrix;
	int i;
	int j;
	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			newMatrix.matrix[i][j] = matrix[i][j] / v;
		}
	}
	return newMatrix;
}

Matrix3d Matrix3d::operator+(const Matrix3d& matr) {
	Matrix3d newMatrix;
	int i;
	int j;
	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			newMatrix.matrix[i][j] = matr.matrix[i][j] + matrix[i][j];
		}
	}
	return newMatrix;
}

Matrix3d Matrix3d::operator-(const Matrix3d& matr) {
	Matrix3d newMatrix;
	int i;
	int j;
	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			newMatrix.matrix[i][j] = matrix[i][j] - matr.matrix[i][j];
		}
	}
	return newMatrix;
}

Matrix3d& Matrix3d::operator+=(const Matrix3d& matr) {
	int i;
	int j;
	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			matrix[i][j] = matr.matrix[i][j] + matrix[i][j];
		}
	}
	return *this;
}

Matrix3d& Matrix3d::operator-=(const Matrix3d& matr) {
	int i;
	int j;
	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			matrix[i][j] = matrix[i][j] - matr.matrix[i][j];
		}
	}
	return *this;
}

Matrix3d* Matrix3d::createBasisByOneVector(const Vector3d& v) {
	if ((v.z * v.z + v.y * v.y + v.x * v.x) < epsilon) {
		return new Matrix3d(1, 0, 0, 0, 1, 0, 0, 0, 1);
	}
	double cosTheta = v.z / sqrt(v.z * v.z + v.y * v.y + v.x * v.x);
	double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
	double cosPhi = v.x / sqrt(v.x * v.x + v.y * v.y);
	double sinPhi = v.y / sqrt(v.x * v.x + v.y * v.y);
	return new Matrix3d(sinPhi, cosTheta * cosPhi, sinTheta * cosPhi, -cosPhi, cosTheta * sinPhi, sinTheta * sinPhi, 0, -sinTheta, cosTheta);
}
