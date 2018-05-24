#ifndef VECTOR3D_H
#define VECTOR3D_H

#include "stdio.h"

class Matrix3d;

class Vector3d{
public:
	double x;
	double y;
	double z;
	
	Vector3d();
	Vector3d(const double& vx, const double& vy, const double& vz);
	Vector3d(const Vector3d& vector);

	Vector3d& operator=(const Vector3d& vector);

	double norm() const;
	double norm2() const;
	Vector3d operator-(const Vector3d& vector) const;
	Vector3d operator+(const Vector3d& vector) const;
	Vector3d& operator+=(const Vector3d& vector);
	Vector3d& operator-=(const Vector3d& vector);
	Vector3d operator*(const double& value) const;
	Vector3d operator/(const double& value) const;
	double scalarMult(const Vector3d& vector) const;
	Vector3d vectorMult(const Vector3d& vector) const;
	Matrix3d tensorMult(const Vector3d& vector) const;
	Matrix3d selfTensorMult() const;

	double& operator[](int i);
};


#endif