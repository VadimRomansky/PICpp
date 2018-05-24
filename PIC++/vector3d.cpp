#include <mpi.h>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "vector3d.h"
#include "matrix3d.h"

Vector3d::Vector3d() {
	x = 0;
	y = 0;
	z = 0;
}

Vector3d::Vector3d(const double& vx, const double& vy, const double& vz) {
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

double Vector3d::norm() const {
	return sqrt(x * x + y * y + z * z);
}

double Vector3d::norm2() const {
	return x * x + y * y + z * z;
}

Vector3d Vector3d::operator+(const Vector3d& vector) const {
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

Vector3d Vector3d::operator-(const Vector3d& vector) const {
	return Vector3d(x - vector.x, y - vector.y, z - vector.z);
}

Vector3d Vector3d::operator*(const double& value) const {
	return Vector3d(x * value, y * value, z * value);
}

Vector3d Vector3d::operator/(const double& value) const {
	return Vector3d(x / value, y / value, z / value);
}

double Vector3d::scalarMult(const Vector3d& vector) const {
	return x * vector.x + y * vector.y + z * vector.z;
}

Vector3d Vector3d::vectorMult(const Vector3d& vector) const {
	double vx = y * vector.z - z * vector.y;
	double vy = -x * vector.z + z * vector.x;
	double vz = x * vector.y - y * vector.x;

	return Vector3d(vx, vy, vz);
}

Matrix3d Vector3d::tensorMult(const Vector3d& vector) const {
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

Matrix3d Vector3d::selfTensorMult() const {
	Matrix3d result;
	result.matrix[0][0] = x * x;
	result.matrix[0][1] = x * y;
	result.matrix[0][2] = x * z;
	result.matrix[1][0] = result.matrix[0][1];
	result.matrix[1][1] = y * y;
	result.matrix[1][2] = y * z;
	result.matrix[2][0] = result.matrix[0][2];
	result.matrix[2][1] = result.matrix[1][2];
	result.matrix[2][2] = z * z;

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
		MPI_Finalize();
		exit(0);
	}
}
