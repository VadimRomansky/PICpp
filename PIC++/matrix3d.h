#ifndef MATRIX3D_H
#define MATRIX3D_H

class Vector3d;

class Matrix3d{
public:
  double matrix[3][3];
  //double** matrix;
  Matrix3d();
  Matrix3d(const double& m11,const double& m12,const double& m13,const double& m21,const double& m22,const double& m23,const double& m31,const double& m32,const double& m33);
  Matrix3d(const Matrix3d& matrix);

  ~Matrix3d();
  Matrix3d* Inverse();
  double determinant();
  double getvalue(int i,int j);
  void setvalue(int i,int j, double value);
  Matrix3d& operator=(const Matrix3d& matr);
  Matrix3d operator*(const double& v);
  Matrix3d operator/(const double& v);
  Vector3d operator*(const Vector3d& v);
  Matrix3d operator*(const Matrix3d& m);
  Matrix3d operator+(const Matrix3d& m);
  Matrix3d operator-(const Matrix3d& m);
  Matrix3d& operator+=(const Matrix3d& matr);
  Matrix3d& operator-=(const Matrix3d& matr);
  static Matrix3d* createBasisByOneVector(const Vector3d& v);
};


#endif