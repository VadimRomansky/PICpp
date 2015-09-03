#ifndef MATRIX3D_H
#define MATRIX3D_H

#include "vector3d.h"

class Vector3d;

class Matrix3d{
public:
  double** matrix;
  Matrix3d();
  Matrix3d(double m11,double m12,double m13,double m21,double m22,double m23,double m31,double m32,double m33);
  Matrix3d(const Matrix3d& matrix);

  ~Matrix3d();
  Matrix3d* Inverse();
  double determinant();
  double getvalue(int i,int j);
  void setvalue(int i,int j, double value);
  Matrix3d& operator=(const Matrix3d& matr);
  Matrix3d operator*(const double v);
  Matrix3d operator/(const double v);
  Vector3d operator*(const Vector3d& v);
  Matrix3d operator*(const Matrix3d& m);
  Matrix3d operator+(const Matrix3d& m);
  Matrix3d operator-(const Matrix3d& m);
  Matrix3d& operator+=(const Matrix3d& matr);
  Matrix3d& operator-=(const Matrix3d& matr);
  static Matrix3d* createBasisByOneVector(const Vector3d& v);
};



#endif