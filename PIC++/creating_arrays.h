#ifndef CREATING_ARRAYS_H_
#define CREATING_ARRAYS_H_
#include "vector3d.h"
#include "matrix3d.h"

double*** create3array(int xnumber, int ynumber, int znumber);
double**** create4array(int xnumber, int ynumber, int znumber, int lnumber);
Vector3d*** create3vectorArray(int xnumber, int ynumber, int znumber);
Vector3d**** create4vectorArray(int xnumber, int ynumber, int znumber, int lnumber);
Matrix3d*** create3matrixArray(int xnumber, int ynumber, int znumber);
Matrix3d**** create4matrixArray(int xnumber, int ynumber, int znumber, int lnumber);

void delete3array(double*** array, int xnumber, int ynumber, int znumber);
void delete4array(double**** array, int xnumber, int ynumber, int znumber, int lnumber);
void delete3vectorArray(Vector3d*** array, int xnumber, int ynumber, int znumber);
void delete4vectorArray(Vector3d**** array, int xnumber, int ynumber, int znumber, int lnumber);
void delete3matrixArray(Matrix3d*** array, int xnumber, int ynumber, int znumber);
void delete4matrixArray(Matrix3d**** array, int xnumber, int ynumber, int znumber, int lnumber);

#endif