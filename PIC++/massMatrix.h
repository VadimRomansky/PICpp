#ifndef MASS_MATRIX_H
#define MASS_MATRIX_H

#include "constants.h"
#include "matrix3d.h"

struct MassMatrix{
	int xindex[2*splineOrder+3];
	int yindex[2*splineOrder+3];
	int zindex[2*splineOrder+3];
	Matrix3d matrix[2*splineOrder+3][2*splineOrder+3][2*splineOrder+3];
};

#endif