#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "specialmath.h"


void Simulation::cleanupDivergence() {
	printf("cleaning up divergence\n");


	int matrixDimension = xnumber;

	double fullDensity = 0;
	/*for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				fullDensity += chargeDensity[i][j][k] * volumeB(i, j, k);
			}
		}
	}*/
	fullDensity /= xsize;

	/*for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				chargeDensity[i][j][k] -= fullDensity;
			}
		}
	}*/


	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				divergenceCleanUpMatrix[i][j][k][0].clear();
				divergenceCleanUpMatrix[i][j][k][1].clear();
				divergenceCleanUpMatrix[i][j][k][2].clear();
				divergenceCleanUpRightPart[i][j][k][0] = 0;
				divergenceCleanUpRightPart[i][j][k][1] = 0;
				divergenceCleanUpRightPart[i][j][k][2] = 0;
				if (((i > 0) && (i < xnumber)) || (boundaryConditionType == PERIODIC)) {
					createDivergenceCleanupInternalEquation(i, j, k);
				} else if (i == 0) {
					createDivergenceCleanupLeftEquation(j, k);
				} else {
					createDivergenceCleanupRightEquation(j, k);
				}
			}
		}
	}


	if (debugMode) {
		checkEquationMatrix(divergenceCleanUpMatrix, 3);
	}

	/*double summRightPart = 0;
	for(int i = 0; i < xnumber; ++i){
		summRightPart += divergenceCleanUpRightPart[i][0];
	}*/
	//double rightPart = cleanUpRightPart(1);
	generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumber, ynumber, znumber, 1);
	//generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningField, xnumber, ynumber, znumber, 3);


	//double** leftPart = multiplySpecialMatrixVector(divergenceCleanUpMatrix, divergenceCleaningPotential, xnumber, 1);

	updateFieldByCleaning();
	//double div = evaluateDivCleaningE(1);

	updateBoundariesNewField();
}

void Simulation::updateFieldByCleaning() {
	evaluateDivergenceCleaningField();
	/*Vector3d meanField = Vector3d(0, 0, 0);
	for(int i = 0; i < xnumber; ++i){
		meanField.x += divergenceCleaningField[i][0];
		meanField.y += divergenceCleaningField[i][1];
		meanField.z += divergenceCleaningField[i][2];
	}
	meanField.x /= xnumber;
	meanField.y /= xnumber;
	meanField.z /= xnumber;*/
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				newEfield[i][j][k].x += divergenceCleaningField[i][j][k][0];
				newEfield[i][j][k].y += divergenceCleaningField[i][j][k][1];
				newEfield[i][j][k].z += divergenceCleaningField[i][j][k][2];
				//newEfield[i] -= meanField;
			}
		}
	}
	if (boundaryConditionType == PERIODIC) {
		for(int j = 0; j < ynumber+1; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				newEfield[xnumber][j][k] = newEfield[0][j][k];
			}
		}
	}
}

void Simulation::evaluateDivergenceCleaningField() {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				int prevI = i - 1;
				if (prevI < 0) {
					if(boundaryConditionType == PERIODIC){
						prevI = xnumber - 1;
					} else {
						prevI = 0;
					}
				}

				int prevJ = j - 1;
				if (prevJ < 0) {
					prevJ = ynumber - 1;
				}

				int prevK = k - 1;
				if (prevK < 0) {
					prevK = znumber - 1;
				}

				divergenceCleaningField[i][j][k][0] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i][prevJ][k][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[i][prevJ][prevK][0]
					- divergenceCleaningPotential[prevI][j][k][0] - divergenceCleaningPotential[prevI][prevJ][k][0] - divergenceCleaningPotential[prevI][j][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0]) / (4 * deltaX);

				divergenceCleaningField[i][j][k][1] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[prevI][j][k][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[prevI][j][prevK][0]
					- divergenceCleaningPotential[i][prevJ][k][0] - divergenceCleaningPotential[prevI][prevJ][k][0] - divergenceCleaningPotential[i][prevJ][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0]) / (4 * deltaY);

				divergenceCleaningField[i][j][k][2] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i][prevJ][k][0] + divergenceCleaningPotential[prevI][j][k][0] + divergenceCleaningPotential[prevI][prevJ][k][0]
					- divergenceCleaningPotential[i][j][prevK][0] - divergenceCleaningPotential[i][prevJ][prevK][0] - divergenceCleaningPotential[prevI][j][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0]) / (4 * deltaZ);
			}
		}
	}

	Vector3d constantField;
	constantField.x = divergenceCleaningField[xnumber - 1][0][0][0];
	constantField.y = divergenceCleaningField[xnumber - 1][0][0][1];
	constantField.z = divergenceCleaningField[xnumber - 1][0][0][2];

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				divergenceCleaningField[i][j][k][0] -= constantField.x;
				divergenceCleaningField[i][j][k][1] -= constantField.y;
				divergenceCleaningField[i][j][k][2] -= constantField.z;
			}
		}
	}
}

void Simulation::createDivergenceCleanupInternalEquation(int i, int j, int k) {

	if (i == 0 && j == 0 && k == 0) {
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
		divergenceCleanUpRightPart[i][j][k][0] = 0;
		return;
	}

	if (i == xnumber - 1 && j == ynumber - 1 && k == znumber - 1) {
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
		divergenceCleanUpRightPart[i][j][k][0] = 0;
		return;
	}

	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;
	}
	int nextI = i + 1;
	if (nextI >= xnumber) {
		nextI = 0;
	}

	int prevJ = j - 1;
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int nextJ = j + 1;
	if (nextJ >= ynumber) {
		nextJ = 0;
	}

	int prevK = k - 1;
	if (prevK < 0) {
		prevK = ynumber - 1;
	}
	int nextK = k + 1;
	if (nextK >= ynumber) {
		nextK = 0;
	}


	divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k);

	double element = -1 / (2 * deltaX2) - 1 / (2 * deltaY2) - 1 / (2 * deltaZ2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

	element = 1 / (16 * deltaX2) + 1 / (16 * deltaY2) + 1 / (16 * deltaZ2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, nextK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, prevK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, nextK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, prevK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, nextK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, prevK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, nextK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, prevK, 0));

	element = 1 / (8 * deltaX2) + 1 / (8 * deltaY2) - 1 / (8 * deltaZ2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 0));

	element = 1 / (8 * deltaX2) - 1 / (8 * deltaY2) + 1 / (8 * deltaZ2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 0));

	element = - 1 / (8 * deltaX2) + 1 / (8 * deltaY2) + 1 / (8 * deltaZ2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, nextK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, prevK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, nextK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, prevK, 0));

	element = 1 / (4 * deltaX2) - 1 / (4 * deltaY2) - 1 / (4 * deltaZ2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

	element = - 1 / (4 * deltaX2) + 1 / (4 * deltaY2) - 1 / (4 * deltaZ2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

	element = - 1 / (4 * deltaX2) - 1 / (4 * deltaY2) + 1 / (4 * deltaZ2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
}

void Simulation::createDivergenceCleanupLeftEquation(int j, int k) {
	int i = 0;

	int nextI = i + 1;

	int prevJ = j - 1;
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int nextJ = j + 1;
	if (nextJ >= ynumber) {
		nextJ = 0;
	}

	int prevK = k - 1;
	if (prevK < 0) {
		prevK = ynumber - 1;
	}
	int nextK = k + 1;
	if (nextK >= ynumber) {
		nextK = 0;
	}

	//div for x

	divergenceCleanUpRightPart[i][0] = 0;
}

void Simulation::createDivergenceCleanupRightEquation(int j, int k) {
	int i = xnumber - 1;

	int prevJ = j - 1;
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int nextJ = j + 1;
	if (nextJ >= ynumber) {
		nextJ = 0;
	}

	int prevK = k - 1;
	if (prevK < 0) {
		prevK = ynumber - 1;
	}
	int nextK = k + 1;
	if (nextK >= ynumber) {
		nextK = 0;
	}

	//div for x


	divergenceCleanUpRightPart[i][2] = 0;
}

double Simulation::cleanUpRightPart(int i, int j, int k) {
	double div = evaluateDivNewE(i, j, k);

	return 4 * pi * chargeDensity[i][j][k] - div;
}
