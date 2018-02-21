#include "massMatrix.h"
#include "constants.h"

MassMatrix::MassMatrix() {
	for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
		for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
			for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
				for (int curI = 0; curI < 3; ++curI) {
					for (int curJ = 0; curJ < 3; ++curJ) {
						matrix[tempI][tempJ][tempK].matrix[curI][curJ] = 0;
					}
				}
				xindex[tempI] = tempI - splineOrder - 1;
				yindex[tempJ] = tempJ - splineOrder - 1;
				zindex[tempK] = tempK - splineOrder - 1;
			}
		}
	}
}

MassMatrix::MassMatrix(int i, int j, int k) {
	for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
		for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
			for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
				for (int curI = 0; curI < 3; ++curI) {
					for (int curJ = 0; curJ < 3; ++curJ) {
						matrix[tempI][tempJ][tempK].matrix[curI][curJ] = 0;
					}
				}
				xindex[tempI] = i + tempI - splineOrder - 1;
				yindex[tempJ] = j + tempJ - splineOrder - 1;
				zindex[tempK] = k + tempK - splineOrder - 1;
			}
		}
	}
}
