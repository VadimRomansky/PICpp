#ifndef VECTOR3D_H
#define VECTOR3D_H

#include "stdio.h"

#include "matrix3d.h"

class Matrix3d;

class Vector3d{
public:
	double x;
	double y;
	double z;
	
	Vector3d();
	Vector3d(double vx, double vy, double vz);
	Vector3d(const Vector3d& vector);

	Vector3d& operator=(const Vector3d& vector);

	double norm();
	Vector3d operator-(const Vector3d& vector);
	Vector3d operator+(const Vector3d& vector);
	Vector3d& operator+=(const Vector3d& vector);
	Vector3d& operator-=(const Vector3d& vector);
	Vector3d operator*(const double& value);
	Vector3d operator/(const double& value);
	double scalarMult(const Vector3d& vector);
	Vector3d vectorMult(const Vector3d& vector);
	Matrix3d tensorMult(const Vector3d& vector);

	double& operator[](int i);
};

#endif