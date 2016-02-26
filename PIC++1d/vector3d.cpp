#include "math.h"
#include "stdlib.h"
#include "stdio.h"

#include "vector3d.h"
#include "matrix3d.h"

Vector3d::Vector3d() {
	x = 0;
	y = 0;
	z = 0;
}

Vector3d::Vector3d(double vx, double vy, double vz) {
	x = vx;
	y = vy;
	z = vz;
}

Vector3d::Vector3d(Vector3d const& vector) {
	x = vector.x;
	y = vector.y;
	z = vector.z;
}

Vector3d& Vector3d::operator=(Vector3d const& vector) {
	x = vector.x;
	y = vector.y;
	z = vector.z;

	return *this;
}

double Vector3d::norm() {
	return sqrt(x * x + y * y + z * z);
}

Vector3d Vector3d::operator+(const Vector3d& vector) {
	return Vector3d(x + vector.x, y + vector.y, z + vector.z);
}

Vector3d& Vector3d::operator+=(Vector3d const& vector) {
	x += vector.x;
	y += vector.y;
	z += vector.z;

	return *this;
}

Vector3d& Vector3d::operator-=(Vector3d const& vector) {
	x -= vector.x;
	y -= vector.y;
	z -= vector.z;

	return *this;
}

Vector3d Vector3d::operator-(const Vector3d& vector) {
	return Vector3d(x - vector.x, y - vector.y, z - vector.z);
}

Vector3d Vector3d::operator*(const double& value) {
	return Vector3d(x * value, y * value, z * value);
}

Vector3d Vector3d::operator/(const double& value) {
	return Vector3d(x / value, y / value, z / value);
}

double Vector3d::scalarMult(const Vector3d& vector) {
	return x * vector.x + y * vector.y + z * vector.z;
}

Vector3d Vector3d::vectorMult(const Vector3d& vector) {
	double vx = y * vector.z - z * vector.y;
	double vy = -x * vector.z + z * vector.x;
	double vz = x * vector.y - y * vector.x;

	return Vector3d(vx, vy, vz);
}

Matrix3d Vector3d::tensorMult(const Vector3d& vector) {
	Matrix3d result;
	result.matrix[0][0] = x * vector.x;
	result.matrix[0][1] = x * vector.y;
	result.matrix[0][2] = x * vector.z;
	result.matrix[1][0] = y * vector.x;
	result.matrix[1][1] = y * vector.y;
	result.matrix[1][2] = y * vector.z;
	result.matrix[2][0] = z * vector.x;
	result.matrix[2][1] = z * vector.y;
	result.matrix[2][2] = z * vector.z;

	return result;
}

double& Vector3d::operator[](int i) {
	switch (i) {
	case 0:
		return x;
	case 1:
		return y;
	case 2:
		return z;
	default:
		printf("i must be 0 < i < 3\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "i must be 0 < i < 3\n");
		fclose(errorLogFile);
		exit(0);
	}
}
