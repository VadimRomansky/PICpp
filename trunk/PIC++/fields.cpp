#include "stdlib.h"
#include "stdio.h"
#include <omp.h>

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"
#include "specialmath.h"

void Simulation::evaluateFields() {
	printf("evaluating fields\n");

	updateElectroMagneticParameters();

 	evaluateMaxwellEquationMatrix();

	double**** gmresOutput = new double***[xnumber];
#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		gmresOutput[i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			gmresOutput[i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				gmresOutput[i][j][k] = new double[3];
			}
		}
	}

	generalizedMinimalResidualMethod(maxwellEquationMatrix, maxwellEquationRightPart, gmresOutput, xnumber, ynumber, znumber, 3);
#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempEfield[i][j][k][l] = gmresOutput[i][j][k][l];
				}
				delete[] gmresOutput[i][j][k];
			}
			delete[] gmresOutput[i][j];
		}
		delete[] gmresOutput[i];
	}
	delete[] gmresOutput;

	updateBoundaries();

	//evaluateMagneticField();

#pragma omp parallel for
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				newEfield[i][j][k] = (tempEfield[i][j][k] - Efield[i][j][k] * (1 - theta)) / theta;
			}
		}
	}

	updateBoundariesNewField();
}

void Simulation::updateEfield()
{
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k] = newEfield[i][j][k];
			}
		}
	}
}

void Simulation::updateBfield()
{
	for (int i = 0; i < xnumber ; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber ; ++k) {
				Bfield[i][j][k] = newBfield[i][j][k];
			}
		}
	}
}

void Simulation::updateFields() {
	updateEfield();
	updateBfield();
}

void Simulation::evaluateMaxwellEquationMatrix() {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					maxwellEquationMatrix[i][j][k][l].clear();
				}
			}
		}
	}

#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				double rightPartVector[3];
				if (boundaryConditionType == SUPERCONDUCTERLEFT) {
					if (i == 0) {
						createPerfectConductaryBoundaryCondition(j, k);
					} else {
						createInternalEquation(i, j, k);
					}
				}
				if (boundaryConditionType == PERIODIC) {
					createInternalEquation(i, j, k);
				}
			}
		}
	}

	if (debugMode) {
		checkEquationMatrix(maxwellEquationMatrix, 3);
	}
}

void Simulation::checkEquationMatrix(std::vector<MatrixElement>**** matrix, int lnumber) {
#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];
						if (element.i < 0) {
							printf("element i < 0\n");
							exit(0);
						}
						if (element.i >= xnumber) {
							printf("element i >= xnumber");
							exit(0);
						}
						if (element.j < 0) {
							printf("element j < 0\n");
							exit(0);
						}
						if (element.j >= ynumber) {
							printf("element j >= ynumber\n");
							exit(0);
						}
						if (element.k < 0) {
							printf("element k < 0\n");
							exit(0);
						}
						if (element.k >= znumber) {
							printf("eement k >= znumber\n");
							exit(0);
						}
						for (int n = m + 1; n < matrix[i][j][k][l].size(); ++n) {
							MatrixElement tempElement = matrix[i][j][k][l][n];

							if (element.equalsIndex(tempElement)) {
								printf("equals indexes\n");
								printf("current = %d %d %d %d\n", i, j, k, l);
								printf("temp = %d %d %d %d\n", element.i, element.j, element.k, element.l);
								exit(0);
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::createPerfectConductaryBoundaryCondition(int j, int k) {
	int i = 0;
	int nextJ = j + 1;
	int nextK = k + 1;
	if (nextJ >= ynumber) {
		nextJ = 0;
	}
	if (nextK >= znumber) {
		nextK = 0;
	}

	/*maxwellEquationRightPart[i][j][k][0] = 4 * pi * electricDensity[0][j][k]*deltaX*deltaY*deltaZ;

	double element = -deltaZ * deltaY * (1 + dielectricTensor[0][j][k].matrix[0][0]) / 4;
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, 0, j, k, 0));

	element = -deltaZ * deltaY * (1 + dielectricTensor[0][nextJ][k].matrix[0][0]) / 4;
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, 0, nextJ, k, 0));

	element = -deltaZ * deltaY * (1 + dielectricTensor[0][j][nextK].matrix[0][0]) / 4;
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, 0, j, nextK, 0));

	element = -deltaZ * deltaY * (1 + dielectricTensor[0][nextJ][nextK].matrix[0][0]) / 4;
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, 0, nextJ, nextK, 0));

	element = (deltaZ * deltaY * (1 + dielectricTensor[i + 1][j][k].matrix[0][0]) / 4) - (deltaY * deltaX * (dielectricTensor[i + 1][j][k].matrix[2][0]) / 4) - (deltaX * deltaZ * (dielectricTensor[i + 1][j][k].matrix[1][0]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, k, 0));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][j][k].matrix[0][1] / 4) - (deltaY * deltaX * dielectricTensor[i + 1][j][k].matrix[2][1] / 4) - (deltaX * deltaZ * (1 + dielectricTensor[i + 1][j][k].matrix[1][1]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, k, 1));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][j][k].matrix[0][2] / 4) - (deltaY * deltaX * (1 + dielectricTensor[i + 1][j][k].matrix[2][2]) / 4) - (deltaX * deltaZ * dielectricTensor[i + 1][j][k].matrix[1][2] / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, k, 2));

	element = (deltaZ * deltaY * (1 + dielectricTensor[i + 1][nextJ][k].matrix[0][0]) / 4) - (deltaY * deltaX * (dielectricTensor[i + 1][nextJ][k].matrix[2][0]) / 4) + (deltaX * deltaZ * (dielectricTensor[i + 1][nextJ][k].matrix[1][0]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, k, 0));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][nextJ][k].matrix[0][1] / 4) - (deltaY * deltaX * dielectricTensor[i + 1][nextJ][k].matrix[2][1] / 4) + (deltaX * deltaZ * (1 + dielectricTensor[i + 1][nextJ][k].matrix[1][1]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, k, 1));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][nextJ][k].matrix[0][2] / 4) - (deltaY * deltaX * (1 + dielectricTensor[i + 1][nextJ][k].matrix[2][2]) / 4) + (deltaX * deltaZ * dielectricTensor[i + 1][nextJ][k].matrix[1][2] / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, k, 2));

	element = (deltaZ * deltaY * (1 + dielectricTensor[i + 1][j][nextK].matrix[0][0]) / 4) + (deltaY * deltaX * (dielectricTensor[i + 1][j][nextK].matrix[2][0]) / 4) - (deltaX * deltaZ * (dielectricTensor[i + 1][j][nextK].matrix[1][0]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, nextK, 0));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][j][nextK].matrix[0][1] / 4) + (deltaY * deltaX * dielectricTensor[i + 1][j][nextK].matrix[2][1] / 4) - (deltaX * deltaZ * (1 + dielectricTensor[i + 1][j][nextK].matrix[1][1]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, nextK, 1));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][j][nextK].matrix[0][2] / 4) + (deltaY * deltaX * (1 + dielectricTensor[i + 1][j][nextK].matrix[2][2]) / 4) - (deltaX * deltaZ * dielectricTensor[i + 1][j][nextK].matrix[1][2] / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, nextK, 2));

	element = (deltaZ * deltaY * (1 + dielectricTensor[i + 1][nextJ][nextK].matrix[0][0]) / 4) + (deltaY * deltaX * (dielectricTensor[i + 1][nextJ][nextK].matrix[2][0]) / 4) + (deltaX * deltaZ * (dielectricTensor[i + 1][nextJ][nextK].matrix[1][0]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, nextK, 0));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][nextJ][nextK].matrix[0][1] / 4) + (deltaY * deltaX * dielectricTensor[i + 1][nextJ][nextK].matrix[2][1] / 4) + (deltaX * deltaZ * (1 + dielectricTensor[i + 1][nextJ][nextK].matrix[1][1]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, nextK, 1));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][nextJ][nextK].matrix[0][2] / 4) + (deltaY * deltaX * (1 + dielectricTensor[i + 1][nextJ][nextK].matrix[2][2]) / 4) + (deltaX * deltaZ * dielectricTensor[i + 1][nextJ][nextK].matrix[1][2] / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, nextK, 2));*/

	/*double element = -0.25/deltaX;
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, nextK, 0));

	element = 0.25/deltaX;
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i+1, j, k, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i+1, nextJ, k, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i+1, j, nextK, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i+1, nextJ, nextK, 0));*/

	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
	maxwellEquationRightPart[i][j][k][0] = 0;
	//Ey and Ez = 0
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(1.0, i, j, k, 1));
	maxwellEquationRightPart[i][j][k][1] = 0;
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(1.0, i, j, k, 2));
	maxwellEquationRightPart[i][j][k][2] = 0;
}

void Simulation::createInternalEquation(int i, int j, int k) {
	Vector3d rightPart = Efield[i][j][k];

	//rightPart = rightPart + (evaluateRotB(i, j, k) - electricFlux[i][j][k] * 4 * pi / speed_of_light_normalized) * speed_of_light_normalized * theta * deltaT;
	rightPart = rightPart + (evaluateRotB(i, j, k) - electricFlux[i][j][k] * 4 * pi / speed_of_light_normalized) * speed_of_light_normalized * theta * deltaT;
	Vector3d a = (evaluateRotB(i, j, k) - electricFlux[i][j][k] * 4 * pi / speed_of_light_normalized) * speed_of_light_normalized * theta * deltaT;
	rightPart = rightPart - evaluateGradDensity(i, j, k) * 4 * pi * sqr(speed_of_light_normalized * theta * deltaT);
	Vector3d b = evaluateGradDensity(i, j, k) * 4 * pi * sqr(speed_of_light_normalized * theta * deltaT);

	alertNaNOrInfinity(rightPart.x, "right part x = NaN");
	alertNaNOrInfinity(rightPart.y, "right part y = NaN");
	alertNaNOrInfinity(rightPart.z, "right part z = NaN");

	createInternalEquationX(i, j, k, rightPart);
	createInternalEquationY(i, j, k, rightPart);
	createInternalEquationZ(i, j, k, rightPart);

	maxwellEquationRightPart[i][j][k][0] = rightPart.x;
	maxwellEquationRightPart[i][j][k][1] = rightPart.y;
	maxwellEquationRightPart[i][j][k][2] = rightPart.z;
}

void Simulation::createInternalEquationX(int i, int j, int k, Vector3d& rightPart) {
	int prevJ = j - 1;
	int prevK = k - 1;
	int nextJ = j + 1;
	int nextK = k + 1;
	int prevI = i - 1;
	int nextI = i + 1;

	if (prevI < 0) {
		prevI = xnumber - 1;
	}
	if (nextI >= xnumber) {
		nextI = 0;
	}
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	if (prevK < 0) {
		prevK = znumber - 1;
	}
	if (nextJ >= ynumber) {
		nextJ = 0;
	}
	if (nextK >= znumber) {
		nextK = 0;
	}

	double cthetadt2 = sqr(speed_of_light_normalized * theta * deltaT);

	double elementX;
	double elementY;
	double elementZ;

	//E i j k
	elementX = 1 + dielectricTensor[i][j][k].matrix[0][0] + cthetadt2 * ((2 / (deltaX * deltaX)) + (2 / (deltaY * deltaY)) + (2 / (deltaZ * deltaZ)) + (0.5 * dielectricTensor[i][j][k].matrix[0][0] / (deltaX * deltaX)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i, j, k, 0));
	elementY = dielectricTensor[i][j][k].matrix[0][1] + cthetadt2 * (0.5 * dielectricTensor[i][j][k].matrix[0][1] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i, j, k, 1));
	elementZ = dielectricTensor[i][j][k].matrix[0][2] + cthetadt2 * (0.5 * dielectricTensor[i][j][k].matrix[0][2] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i, j, k, 2));

	//E i+1 j k

	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (-(1 / (deltaX * deltaX))  - (0.25 * dielectricTensor[i + 1][j][k].matrix[0][0] / (deltaX * deltaX)));
		elementY = cthetadt2 * (-0.25 * dielectricTensor[i + 1][j][k].matrix[0][1] / (deltaX * deltaX));
		elementZ = cthetadt2 * (-0.25 * dielectricTensor[i + 1][j][k].matrix[0][2] / (deltaX * deltaX));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][0].push_back((MatrixElement(elementX, nextI, j, k, 0)));
			maxwellEquationMatrix[i][j][k][0].push_back((MatrixElement(elementY, nextI, j, k, 1)));
			maxwellEquationMatrix[i][j][k][0].push_back((MatrixElement(elementZ, nextI, j, k, 2)));
		} else {
			rightPart.x -= elementX * E0.x;
			rightPart.x -= elementY * E0.y;
			rightPart.x -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (-(1 / (deltaX * deltaX)) - (0.25 * dielectricTensor[nextI][j][k].matrix[0][0] / (deltaX * deltaX)));
		elementY = cthetadt2 * (-0.25 * dielectricTensor[nextI][j][k].matrix[0][1] / (deltaX * deltaX));
		elementZ = cthetadt2 * (-0.25 * dielectricTensor[nextI][j][k].matrix[0][2] / (deltaX * deltaX));
		maxwellEquationMatrix[i][j][k][0].push_back((MatrixElement(elementX, nextI, j, k, 0)));
		maxwellEquationMatrix[i][j][k][0].push_back((MatrixElement(elementY, nextI, j, k, 1)));
		maxwellEquationMatrix[i][j][k][0].push_back((MatrixElement(elementZ, nextI, j, k, 2)));
	}

	//E i-1 j k
	elementX = cthetadt2 * ((1 / (deltaX * deltaX)) - (0.25 * dielectricTensor[prevI][j][k].matrix[0][0] / (deltaX * deltaX)));
	maxwellEquationMatrix[i][j][k][0].push_back((MatrixElement(elementX, prevI, j, k, 0)));
	elementY = cthetadt2 * (-0.25 * dielectricTensor[prevI][j][k].matrix[0][1] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back((MatrixElement(elementY, prevI, j, k, 1)));
	elementZ = cthetadt2 * (-0.25 * dielectricTensor[prevI][j][k].matrix[0][2] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back((MatrixElement(elementZ, prevI, j, k, 2)));

	//E i j+1 k
	elementX = cthetadt2 * ( - (1 / (deltaY * deltaY)) + (0.25 * dielectricTensor[i][nextJ][k].matrix[0][0] / (deltaX * deltaX)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i, nextJ, k, 0));
	elementY = cthetadt2 * ((0.25 * dielectricTensor[i][nextJ][k].matrix[0][1]) / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i, nextJ, k, 1));
	elementZ = cthetadt2 * ((0.25 * dielectricTensor[i][nextJ][k].matrix[0][2]) / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i, nextJ, k, 2));

	//E i j-1 k
	elementX = cthetadt2 * ( - (1 / (deltaY * deltaY)) + (0.25 * dielectricTensor[i][prevJ][k].matrix[0][0] / (deltaX * deltaX)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i, prevJ, k, 0));
	elementY = cthetadt2 * ((0.25 * dielectricTensor[i][prevJ][k].matrix[0][1]) / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i, prevJ, k, 1));
	elementZ = cthetadt2 * ((0.25 * dielectricTensor[i][prevJ][k].matrix[0][2]) / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i, prevJ, k, 2));

	//E i j k+1
	elementX = cthetadt2 * ( - (1 / (deltaZ * deltaZ)) + (0.25 * dielectricTensor[i][j][nextK].matrix[0][0] / (deltaX * deltaX)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i, j, nextK, 0));
	elementY = cthetadt2 * ((0.25 * dielectricTensor[i][j][nextK].matrix[0][1]) / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i, j, nextK, 1));
	elementZ = cthetadt2 * ((0.25 * dielectricTensor[i][j][nextK].matrix[0][2]) / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i, j, nextK, 2));

	//E i j k-1
	elementX = cthetadt2 * ( - (1 / (deltaZ * deltaZ)) + (0.25 * dielectricTensor[i][j][nextK].matrix[0][0] / (deltaX * deltaX)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i, j, prevK, 0));
	elementY = cthetadt2 * ((0.25 * dielectricTensor[i][j][prevK].matrix[0][1]) / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i, j, prevK, 1));
	elementZ = cthetadt2 * ((0.25 * dielectricTensor[i][j][prevK].matrix[0][2]) / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i, j, prevK, 2));

	//E i+1 j+1 k
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ( - (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[0][0] / (deltaX * deltaX)) - (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[1][0] / (deltaX * deltaY)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[0][1] / (deltaX * deltaX)) - (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[1][1] / (deltaX * deltaY)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[0][2] / (deltaX * deltaX)) - (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[1][2] / (deltaX * deltaY)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i + 1, nextJ, k, 0));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i + 1, nextJ, k, 1));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i + 1, nextJ, k, 2));
		} else {
			rightPart.x -= elementX * E0.x;
			rightPart.x -= elementY * E0.y;
			rightPart.x -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ( - (0.125 * dielectricTensor[nextI][nextJ][k].matrix[0][0] / (deltaX * deltaX)) - (0.125 * dielectricTensor[nextI][nextJ][k].matrix[1][0] / (deltaX * deltaY)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[nextI][nextJ][k].matrix[0][1] / (deltaX * deltaX)) - (0.125 * dielectricTensor[nextI][nextJ][k].matrix[1][1] / (deltaX * deltaY)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[nextI][nextJ][k].matrix[0][2] / (deltaX * deltaX)) - (0.125 * dielectricTensor[nextI][nextJ][k].matrix[1][2] / (deltaX * deltaY)));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, nextI, nextJ, k, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, nextI, nextJ, k, 1));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, nextI, nextJ, k, 2));
	}

	//E i+1 j-1 k
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ( - (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[0][0] / (deltaX * deltaX)) + (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[1][0] / (deltaX * deltaY)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[0][1] / (deltaX * deltaX)) + (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[1][1] / (deltaX * deltaY)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[0][2] / (deltaX * deltaX)) + (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[1][2] / (deltaX * deltaY)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i + 1, prevJ, k, 0));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i + 1, prevJ, k, 1));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i + 1, prevJ, k, 2));
		} else {
			rightPart.x -= elementX * E0.x;
			rightPart.x -= elementY * E0.y;
			rightPart.x -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ( - (0.125 * dielectricTensor[nextI][prevJ][k].matrix[0][0] / (deltaX * deltaX)) + (0.125 * dielectricTensor[nextI][prevJ][k].matrix[1][0] / (deltaX * deltaY)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[nextI][prevJ][k].matrix[0][1] / (deltaX * deltaX)) + (0.125 * dielectricTensor[nextI][prevJ][k].matrix[1][1] / (deltaX * deltaY)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[nextI][prevJ][k].matrix[0][2] / (deltaX * deltaX)) + (0.125 * dielectricTensor[nextI][prevJ][k].matrix[1][2] / (deltaX * deltaY)));

		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, nextI, prevJ, k, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, nextI, prevJ, k, 1));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, nextI, prevJ, k, 2));
	}

	//E i-1 j+1 k
	elementX = cthetadt2 * ( - (0.125 * dielectricTensor[prevI][nextJ][k].matrix[0][0] / (deltaX * deltaX)) + (0.125 * dielectricTensor[prevI][nextJ][k].matrix[1][0] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, prevI, nextJ, k, 0));
	elementY = cthetadt2 * (- (0.125 * dielectricTensor[prevI][nextJ][k].matrix[0][1] / (deltaX * deltaX)) + (0.125 * dielectricTensor[prevI][nextJ][k].matrix[1][1] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, prevI, nextJ, k, 1));
	elementZ = cthetadt2 * (- (0.125 * dielectricTensor[prevI][nextJ][k].matrix[0][2] / (deltaX * deltaX)) + (0.125 * dielectricTensor[prevI][nextJ][k].matrix[1][2] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, prevI, nextJ, k, 2));

	//E i-1 j-1 k
	elementX = cthetadt2 * ( - (0.125 * dielectricTensor[prevI][prevJ][k].matrix[0][0] / (deltaX * deltaX)) - (0.125 * dielectricTensor[prevI][prevJ][k].matrix[1][0] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, prevI, prevJ, k, 0));
	elementY = cthetadt2 * (- (0.125 * dielectricTensor[prevI][prevJ][k].matrix[0][1] / (deltaX * deltaX)) - (0.125 * dielectricTensor[prevI][prevJ][k].matrix[1][1] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, prevI, prevJ, k, 1));
	elementZ = cthetadt2 * (- (0.125 * dielectricTensor[prevI][prevJ][k].matrix[0][2] / (deltaX * deltaX)) - (0.125 * dielectricTensor[prevI][prevJ][k].matrix[1][2] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, prevI, prevJ, k, 2));

	//E i+1 j k+1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ( - (0.125 * dielectricTensor[i + 1][j][nextK].matrix[0][0] / (deltaX * deltaX)) - (0.125 * dielectricTensor[i + 1][j][nextK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][j][nextK].matrix[0][1] / (deltaX * deltaX)) - (0.125 * dielectricTensor[i + 1][j][nextK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][j][nextK].matrix[0][2] / (deltaX * deltaX)) - (0.125 * dielectricTensor[i + 1][j][nextK].matrix[2][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i + 1, j, nextK, 0));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i + 1, j, nextK, 1));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i + 1, j, nextK, 2));
		} else {
			rightPart.x -= elementX * E0.x;
			rightPart.x -= elementY * E0.y;
			rightPart.x -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ( - (0.125 * dielectricTensor[nextI][j][nextK].matrix[0][0] / (deltaX * deltaX)) - (0.125 * dielectricTensor[nextI][j][nextK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[nextI][j][nextK].matrix[0][1] / (deltaX * deltaX)) - (0.125 * dielectricTensor[nextI][j][nextK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[nextI][j][nextK].matrix[0][2] / (deltaX * deltaX)) - (0.125 * dielectricTensor[nextI][j][nextK].matrix[2][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, nextI, j, nextK, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, nextI, j, nextK, 1));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, nextI, j, nextK, 2));
	}

	//E i+1 j k-1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ( - (0.125 * dielectricTensor[i + 1][j][prevK].matrix[0][0] / (deltaX * deltaX)) + (0.125 * dielectricTensor[i + 1][j][prevK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][j][prevK].matrix[0][1] / (deltaX * deltaX)) + (0.125 * dielectricTensor[i + 1][j][prevK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][j][prevK].matrix[0][2] / (deltaX * deltaX)) + (0.125 * dielectricTensor[i + 1][j][prevK].matrix[2][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i + 1, j, prevK, 0));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i + 1, j, prevK, 1));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i + 1, j, prevK, 2));
		} else {
			rightPart.x -= elementX * E0.x;
			rightPart.x -= elementY * E0.y;
			rightPart.x -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ( - (0.125 * dielectricTensor[nextI][j][prevK].matrix[0][0] / (deltaX * deltaX)) + (0.125 * dielectricTensor[nextI][j][prevK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[nextI][j][prevK].matrix[0][1] / (deltaX * deltaX)) + (0.125 * dielectricTensor[nextI][j][prevK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[nextI][j][prevK].matrix[0][2] / (deltaX * deltaX)) + (0.125 * dielectricTensor[nextI][j][prevK].matrix[2][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, nextI, j, prevK, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, nextI, j, prevK, 1));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, nextI, j, prevK, 2));
	}

	//E i-1 j k+1
	elementX = cthetadt2 * ( - (0.125 * dielectricTensor[prevI][j][nextK].matrix[0][0] / (deltaX * deltaX)) + (0.125 * dielectricTensor[prevI][j][nextK].matrix[2][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, prevI, j, nextK, 0));
	elementY = cthetadt2 * (- (0.125 * dielectricTensor[prevI][j][nextK].matrix[0][1] / (deltaX * deltaX)) + (0.125 * dielectricTensor[prevI][j][nextK].matrix[2][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, prevI, j, nextK, 1));
	elementZ = cthetadt2 * (- (0.125 * dielectricTensor[prevI][j][nextK].matrix[0][2] / (deltaX * deltaX)) + (0.125 * dielectricTensor[prevI][j][nextK].matrix[2][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, prevI, j, nextK, 2));

	//E i-1 j k-1
	elementX = cthetadt2 * ( - (0.125 * dielectricTensor[prevI][j][prevK].matrix[0][0] / (deltaX * deltaX)) - (0.125 * dielectricTensor[prevI][j][prevK].matrix[2][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, prevI, j, prevK, 0));
	elementY = cthetadt2 * (- (0.125 * dielectricTensor[prevI][j][prevK].matrix[0][1] / (deltaX * deltaX)) - (0.125 * dielectricTensor[prevI][j][prevK].matrix[2][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, prevI, j, prevK, 1));
	elementZ = cthetadt2 * (- (0.125 * dielectricTensor[prevI][j][prevK].matrix[0][2] / (deltaX * deltaX)) - (0.125 * dielectricTensor[prevI][j][prevK].matrix[2][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, prevI, j, prevK, 2));

	//E i j+1 k+1
	elementX = cthetadt2 * ( + (0.125 * dielectricTensor[i][nextJ][nextK].matrix[0][0] / (deltaX * deltaX)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i, nextJ, nextK, 0));
	elementY = cthetadt2 * (0.125 * dielectricTensor[i][nextJ][nextK].matrix[0][1] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i, nextJ, nextK, 1));
	elementZ = cthetadt2 * (0.125 * dielectricTensor[i][nextJ][nextK].matrix[0][2] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i, nextJ, nextK, 2));

	//E i j+1 k-1
	elementX = cthetadt2 * ( + (0.125 * dielectricTensor[i][nextJ][prevK].matrix[0][0] / (deltaX * deltaX)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i, nextJ, prevK, 0));
	elementY = cthetadt2 * (0.125 * dielectricTensor[i][nextJ][prevK].matrix[0][1] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i, nextJ, prevK, 1));
	elementZ = cthetadt2 * (0.125 * dielectricTensor[i][nextJ][prevK].matrix[0][2] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i, nextJ, prevK, 2));

	//E i j-1 k+1
	elementX = cthetadt2 * ( + (0.125 * dielectricTensor[i][prevJ][nextK].matrix[0][0] / (deltaX * deltaX)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i, prevJ, nextK, 0));
	elementY = cthetadt2 * (0.125 * dielectricTensor[i][prevJ][nextK].matrix[0][1] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i, prevJ, nextK, 1));
	elementZ = cthetadt2 * (0.125 * dielectricTensor[i][prevJ][nextK].matrix[0][2] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i, prevJ, nextK, 2));

	//E i j-1 k-1
	elementX = cthetadt2 * ( + (0.125 * dielectricTensor[i][prevJ][prevK].matrix[0][0] / (deltaX * deltaX)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i, prevJ, prevK, 0));
	elementY = cthetadt2 * (0.125 * dielectricTensor[i][prevJ][prevK].matrix[0][1] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i, prevJ, prevK, 1));
	elementZ = cthetadt2 * (0.125 * dielectricTensor[i][prevJ][prevK].matrix[0][2] / (deltaX * deltaX));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i, prevJ, prevK, 2));

	//E i+1 j+1 k+1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[0][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[1][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[0][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[1][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[0][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[1][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[2][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i + 1, nextJ, nextK, 0));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i + 1, nextJ, nextK, 1));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i + 1, nextJ, nextK, 2));
		} else {
			rightPart.x -= elementX * E0.x;
			rightPart.x -= elementY * E0.y;
			rightPart.x -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[0][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[1][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[0][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[1][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[0][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[1][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[2][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, nextI, nextJ, nextK, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, nextI, nextJ, nextK, 1));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, nextI, nextJ, nextK, 2));
	}

	//E i+1 j+1 k-1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[0][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[1][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[0][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[1][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[0][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[1][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[2][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i + 1, nextJ, prevK, 0));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i + 1, nextJ, prevK, 1));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i + 1, nextJ, prevK, 2));
		} else {
			rightPart.x -= elementX * E0.x;
			rightPart.x -= elementY * E0.y;
			rightPart.x -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[0][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[1][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[0][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[1][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[0][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[1][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[2][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, nextI, nextJ, prevK, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, nextI, nextJ, prevK, 1));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, nextI, nextJ, prevK, 2));
	}

	//E i+1 j-1 k+1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[0][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[1][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[0][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[1][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[0][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[1][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[2][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i + 1, prevJ, nextK, 0));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i + 1, prevJ, nextK, 1));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i + 1, prevJ, nextK, 2));
		} else {
			rightPart.x -= elementX * E0.x;
			rightPart.x -= elementY * E0.y;
			rightPart.x -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[0][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[1][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[0][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[1][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[0][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[1][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[2][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, nextI, prevJ, nextK, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, nextI, prevJ, nextK, 1));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, nextI, prevJ, nextK, 2));
	}

	//E i+1 j-1 k-1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[0][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[1][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[0][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[1][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[0][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[1][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[2][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, i + 1, prevJ, prevK, 0));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, i + 1, prevJ, prevK, 1));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, i + 1, prevJ, prevK, 2));
		} else {
			rightPart.x -= elementX * E0.x;
			rightPart.x -= elementY * E0.y;
			rightPart.x -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[0][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[1][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[2][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[0][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[1][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[2][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[0][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[1][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[2][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, nextI, prevJ, prevK, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, nextI, prevJ, prevK, 1));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, nextI, prevJ, prevK, 2));
	}

	//E i-1 j+1 k+1
	elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[0][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[1][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[2][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, prevI, nextJ, nextK, 0));
	elementY = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[0][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[1][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[2][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, prevI, nextJ, nextK, 1));
	elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[0][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[1][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[2][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, prevI, nextJ, nextK, 2));

	//E i-1 j+1 k-1
	elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[0][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[1][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[2][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, prevI, nextJ, prevK, 0));
	elementY = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[0][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[1][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[2][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, prevI, nextJ, prevK, 1));
	elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[0][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[1][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[2][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, prevI, nextJ, prevK, 2));

	//E i-1 j-1 k+1
	elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[0][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[1][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[2][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, prevI, prevJ, nextK, 0));
	elementY = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[0][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[1][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[2][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, prevI, prevJ, nextK, 1));
	elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[0][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[1][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[2][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, prevI, prevJ, nextK, 2));

	//E i-1 j-1 k-1
	elementX = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[0][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[1][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[2][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementX, prevI, prevJ, prevK, 0));
	elementY = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[0][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[1][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[2][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementY, prevI, prevJ, prevK, 1));
	elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[0][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[1][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[2][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(elementZ, prevI, prevJ, prevK, 2));

}

void Simulation::createInternalEquationY(int i, int j, int k, Vector3d& rightPart) {
	int prevJ = j - 1;
	int prevK = k - 1;
	int nextJ = j + 1;
	int nextK = k + 1;
	int prevI = i - 1;
	int nextI = i + 1;

	if (prevI < 0) {
		prevI = xnumber - 1;
	}
	if (nextI >= xnumber) {
		nextI = 0;
	}
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	if (prevK < 0) {
		prevK = znumber - 1;
	}
	if (nextJ >= ynumber) {
		nextJ = 0;
	}
	if (nextK >= znumber) {
		nextK = 0;
	}
	double cthetadt2 = sqr(speed_of_light_normalized * theta * deltaT);

	double elementX;
	double elementY;
	double elementZ;

	//E i j k
	elementX = dielectricTensor[i][j][k].matrix[1][0] + cthetadt2 * (0.5 * dielectricTensor[i][j][k].matrix[1][0] / (deltaY * deltaY));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i, j, k, 0));
	elementY = 1 + dielectricTensor[i][j][k].matrix[1][1] + cthetadt2 * ((2 / (deltaX * deltaX)) + (2 / (deltaY * deltaY)) + (2 / (deltaZ * deltaZ)) + (0.5 * dielectricTensor[i][j][k].matrix[1][1] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i, j, k, 1));
	elementZ = dielectricTensor[i][j][k].matrix[1][2] + cthetadt2 * (0.5 * dielectricTensor[i][j][k].matrix[1][2] / (deltaY * deltaY));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i, j, k, 2));

	//E i j+1 k
	elementX = cthetadt2 * (- (0.25 * dielectricTensor[i][nextJ][k].matrix[1][0] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i, nextJ, k, 0));
	elementY = cthetadt2 * ( - (1 / (deltaY * deltaY)) - (0.25 * dielectricTensor[i][nextJ][k].matrix[1][1] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i, nextJ, k, 1));
	elementZ = cthetadt2 * (- (0.25 * dielectricTensor[i][nextJ][k].matrix[1][2] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i, nextJ, k, 2));

	//E i j-1 k
	elementX = cthetadt2 * (- (0.25 * dielectricTensor[i][prevJ][k].matrix[1][0] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i, prevJ, k, 0));
	elementY = cthetadt2 * ( - (1 / (deltaY * deltaY))  - (0.25 * dielectricTensor[i][prevJ][k].matrix[1][1] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i, prevJ, k, 1));
	elementZ = cthetadt2 * (- (0.25 * dielectricTensor[i][prevJ][k].matrix[1][2] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i, prevJ, k, 2));

	//E i+1 j k
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ((0.25 * dielectricTensor[i + 1][j][k].matrix[1][0] / (deltaY * deltaY)));
		elementY = cthetadt2 * ((-1 / (deltaX * deltaX)) + (0.25 * dielectricTensor[i + 1][j][k].matrix[1][1] / (deltaY * deltaY)));
		elementZ = cthetadt2 * ((0.25 * dielectricTensor[i + 1][j][k].matrix[1][2] / (deltaY * deltaY)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i + 1, j, k, 0));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i + 1, j, k, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i + 1, j, k, 2));
		} else {
			rightPart.y -= elementX * E0.x;
			rightPart.y -= elementY * E0.y;
			rightPart.y -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ((0.25 * dielectricTensor[nextI][j][k].matrix[1][0] / (deltaY * deltaY)));
		elementY = cthetadt2 * ((-1 / (deltaX * deltaX)) + (0.25 * dielectricTensor[nextI][j][k].matrix[1][1] / (deltaY * deltaY)));
		elementZ = cthetadt2 * ((0.25 * dielectricTensor[nextI][j][k].matrix[1][2] / (deltaY * deltaY)));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, nextI, j, k, 0));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, nextI, j, k, 1));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, nextI, j, k, 2));
	}

	//E i-1 j k
	elementX = cthetadt2 * ((0.25 * dielectricTensor[prevI][j][k].matrix[1][0] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, prevI, j, k, 0));
	elementY = cthetadt2 * ((-1 / (deltaX * deltaX)) + (0.25 * dielectricTensor[prevI][j][k].matrix[1][1] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, prevI, j, k, 1));
	elementZ = cthetadt2 * ((0.25 * dielectricTensor[prevI][j][k].matrix[1][2] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, prevI, j, k, 2));

	//E i j k+1
	elementX = cthetadt2 * ((0.25 * dielectricTensor[i][j][nextK].matrix[1][0] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i, j, nextK, 0));
	elementY = cthetadt2 * ( - (1 / (deltaZ * deltaZ)) + (0.25 * dielectricTensor[i][j][nextK].matrix[1][1] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i, j, nextK, 1));
	elementZ = cthetadt2 * ((0.25 * dielectricTensor[i][j][nextK].matrix[1][2] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i, j, nextK, 2));

	//E i j k-1
	elementX = cthetadt2 * ((0.25 * dielectricTensor[i][j][prevK].matrix[1][0] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i, j, prevK, 0));
	elementY = cthetadt2 * ( - (1 / (deltaZ * deltaZ)) + (0.25 * dielectricTensor[i][j][prevK].matrix[1][1] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i, j, prevK, 1));
	elementZ = cthetadt2 * ((0.25 * dielectricTensor[i][j][prevK].matrix[1][2] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i, j, prevK, 2));

	//E i+1 j+1 k
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[1][0] / (deltaY * deltaY)) - (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[0][0] / (deltaX * deltaY)));
		elementY = cthetadt2 * ( - (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[1][1] / (deltaY * deltaY)) - (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[0][1] / (deltaX * deltaY)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[1][2] / (deltaY * deltaY)) - (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[0][2] / (deltaX * deltaY)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i + 1, nextJ, k, 0));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i + 1, nextJ, k, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i + 1, nextJ, k, 2));
		} else {
			rightPart.y -= elementX * E0.x;
			rightPart.y -= elementY * E0.y;
			rightPart.y -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (- (0.125 * dielectricTensor[nextI][nextJ][k].matrix[1][0] / (deltaY * deltaY)) - (0.125 * dielectricTensor[nextI][nextJ][k].matrix[0][0] / (deltaX * deltaY)));
		elementY = cthetadt2 * ( - (0.125 * dielectricTensor[nextI][nextJ][k].matrix[1][1] / (deltaY * deltaY)) - (0.125 * dielectricTensor[nextI][nextJ][k].matrix[0][1] / (deltaX * deltaY)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[nextI][nextJ][k].matrix[1][2] / (deltaY * deltaY)) - (0.125 * dielectricTensor[nextI][nextJ][k].matrix[0][2] / (deltaX * deltaY)));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, nextI, nextJ, k, 0));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, nextI, nextJ, k, 1));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, nextI, nextJ, k, 2));
	}
	//E i+1 j-1 k
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[1][0] / (deltaY * deltaY)) + (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[0][0] / (deltaX * deltaY)));
		elementY = cthetadt2 * ( - (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[1][1] / (deltaY * deltaY)) + (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[0][1] / (deltaX * deltaY)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[1][2] / (deltaY * deltaY)) + (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[0][2] / (deltaX * deltaY)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i + 1, prevJ, k, 0));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i + 1, prevJ, k, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i + 1, prevJ, k, 2));
		} else {
			rightPart.y -= elementX * E0.x;
			rightPart.y -= elementY * E0.y;
			rightPart.y -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (- (0.125 * dielectricTensor[nextI][prevJ][k].matrix[1][0] / (deltaY * deltaY)) + (0.125 * dielectricTensor[nextI][prevJ][k].matrix[0][0] / (deltaX * deltaY)));
		elementY = cthetadt2 * ( - (0.125 * dielectricTensor[nextI][prevJ][k].matrix[1][1] / (deltaY * deltaY)) + (0.125 * dielectricTensor[nextI][prevJ][k].matrix[0][1] / (deltaX * deltaY)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[nextI][prevJ][k].matrix[1][2] / (deltaY * deltaY)) + (0.125 * dielectricTensor[nextI][prevJ][k].matrix[0][2] / (deltaX * deltaY)));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, nextI, prevJ, k, 0));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, nextI, prevJ, k, 1));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, nextI, prevJ, k, 2));
	}

	//E i-1 j+1 k
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[prevI][nextJ][k].matrix[1][0] / (deltaY * deltaY)) + (0.125 * dielectricTensor[prevI][nextJ][k].matrix[0][0] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, prevI, nextJ, k, 0));
	elementY = cthetadt2 * ( - (0.125 * dielectricTensor[prevI][nextJ][k].matrix[1][1] / (deltaY * deltaY)) + (0.125 * dielectricTensor[prevI][nextJ][k].matrix[0][1] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, prevI, nextJ, k, 1));
	elementZ = cthetadt2 * (- (0.125 * dielectricTensor[prevI][nextJ][k].matrix[1][2] / (deltaY * deltaY)) + (0.125 * dielectricTensor[prevI][nextJ][k].matrix[0][2] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, prevI, nextJ, k, 2));

	//E i-1 j-1 k
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[prevI][prevJ][k].matrix[1][0] / (deltaY * deltaY)) - (0.125 * dielectricTensor[prevI][prevJ][k].matrix[0][0] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, prevI, prevJ, k, 0));
	elementY = cthetadt2 * ( - (0.125 * dielectricTensor[prevI][prevJ][k].matrix[1][1] / (deltaY * deltaY)) - (0.125 * dielectricTensor[prevI][prevJ][k].matrix[0][1] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, prevI, prevJ, k, 1));
	elementZ = cthetadt2 * (- (0.125 * dielectricTensor[prevI][prevJ][k].matrix[1][2] / (deltaY * deltaY)) - (0.125 * dielectricTensor[prevI][prevJ][k].matrix[0][2] / (deltaX * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, prevI, prevJ, k, 2));

	//E i j+1 k+1
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[i][nextJ][nextK].matrix[1][0] / (deltaY * deltaY)) - (0.125 * dielectricTensor[i][nextJ][nextK].matrix[2][0] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i, nextJ, nextK, 0));
	elementY = cthetadt2 * ( - (0.125 * dielectricTensor[i][nextJ][nextK].matrix[1][1] / (deltaY * deltaY)) - (0.125 * dielectricTensor[i][nextJ][nextK].matrix[2][1] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i, nextJ, nextK, 1));
	elementZ = cthetadt2 * (- (0.125 * dielectricTensor[i][nextJ][nextK].matrix[1][2] / (deltaY * deltaY)) - (0.125 * dielectricTensor[i][nextJ][nextK].matrix[2][2] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i, nextJ, nextK, 2));

	//E i j+1 k-1
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[i][nextJ][prevK].matrix[1][0] / (deltaY * deltaY)) + (0.125 * dielectricTensor[i][nextJ][prevK].matrix[2][0] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i, nextJ, prevK, 0));
	elementY = cthetadt2 * ( - (0.125 * dielectricTensor[i][nextJ][prevK].matrix[1][1] / (deltaY * deltaY)) + (0.125 * dielectricTensor[i][nextJ][prevK].matrix[2][1] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i, nextJ, prevK, 1));
	elementZ = cthetadt2 * (- (0.125 * dielectricTensor[i][nextJ][prevK].matrix[1][2] / (deltaY * deltaY)) + (0.125 * dielectricTensor[i][nextJ][prevK].matrix[2][2] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i, nextJ, prevK, 2));

	//E i j-1 k+1
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[i][prevJ][nextK].matrix[1][0] / (deltaY * deltaY)) + (0.125 * dielectricTensor[i][nextJ][nextK].matrix[2][0] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i, prevJ, nextK, 0));
	elementY = cthetadt2 * ( - (0.125 * dielectricTensor[i][prevJ][nextK].matrix[1][1] / (deltaY * deltaY)) + (0.125 * dielectricTensor[i][prevJ][nextK].matrix[2][1] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i, prevJ, nextK, 1));
	elementZ = cthetadt2 * (- (0.125 * dielectricTensor[i][nextJ][nextK].matrix[1][2] / (deltaY * deltaY)) + (0.125 * dielectricTensor[i][prevJ][nextK].matrix[2][2] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i, prevJ, nextK, 2));

	//E i j-1 k-1
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[i][prevJ][prevK].matrix[1][0] / (deltaY * deltaY)) - (0.125 * dielectricTensor[i][nextJ][prevK].matrix[2][0] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i, prevJ, prevK, 0));
	elementY = cthetadt2 * ( - (0.125 * dielectricTensor[i][prevJ][prevK].matrix[1][1] / (deltaY * deltaY)) - (0.125 * dielectricTensor[i][prevJ][prevK].matrix[2][1] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i, prevJ, prevK, 1));
	elementZ = cthetadt2 * (- (0.125 * dielectricTensor[i][nextJ][prevK].matrix[1][2] / (deltaY * deltaY)) - (0.125 * dielectricTensor[i][prevJ][prevK].matrix[2][2] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i, prevJ, prevK, 2));

	//E i+1 j k+1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ((0.125 * dielectricTensor[i + 1][j][nextK].matrix[1][0] / (deltaY * deltaY)));
		elementY = cthetadt2 * ( + (0.125 * dielectricTensor[i + 1][j][nextK].matrix[1][1] / (deltaY * deltaY)));
		elementZ = cthetadt2 * ((0.125 * dielectricTensor[i + 1][j][nextK].matrix[1][2] / (deltaY * deltaY)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i + 1, j, nextK, 0));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i + 1, j, nextK, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i + 1, j, nextK, 2));
		} else {
			rightPart.y -= elementX * E0.x;
			rightPart.y -= elementY * E0.y;
			rightPart.y -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ((0.125 * dielectricTensor[nextI][j][nextK].matrix[1][0] / (deltaY * deltaY)));
		elementY = cthetadt2 * ( + (0.125 * dielectricTensor[nextI][j][nextK].matrix[1][1] / (deltaY * deltaY)));
		elementZ = cthetadt2 * ((0.125 * dielectricTensor[nextI][j][nextK].matrix[1][2] / (deltaY * deltaY)));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, nextI, j, nextK, 0));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, nextI, j, nextK, 1));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, nextI, j, nextK, 2));
	}

	//E i+1 j k-1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ((0.125 * dielectricTensor[i + 1][j][prevK].matrix[1][0] / (deltaY * deltaY)));
		elementY = cthetadt2 * ( + (0.125 * dielectricTensor[i + 1][j][prevK].matrix[1][1] / (deltaY * deltaY)));
		elementZ = cthetadt2 * ((0.125 * dielectricTensor[i + 1][j][prevK].matrix[1][2] / (deltaY * deltaY)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i + 1, j, prevK, 0));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i + 1, j, prevK, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i + 1, j, prevK, 2));
		} else {
			rightPart.y -= elementX * E0.x;
			rightPart.y -= elementY * E0.y;
			rightPart.y -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ((0.125 * dielectricTensor[nextI][j][prevK].matrix[1][0] / (deltaY * deltaY)));
		elementY = cthetadt2 * ( + (0.125 * dielectricTensor[nextI][j][prevK].matrix[1][1] / (deltaY * deltaY)));
		elementZ = cthetadt2 * ((0.125 * dielectricTensor[nextI][j][prevK].matrix[1][2] / (deltaY * deltaY)));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, nextI, j, prevK, 0));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, nextI, j, prevK, 1));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, nextI, j, prevK, 2));
	}

	//E i-1 j k+1
	elementX = cthetadt2 * ((0.125 * dielectricTensor[prevI][j][nextK].matrix[1][0] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, prevI, j, nextK, 0));
	elementY = cthetadt2 * ( + (0.125 * dielectricTensor[prevI][j][nextK].matrix[1][1] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, prevI, j, nextK, 1));
	elementZ = cthetadt2 * ((0.125 * dielectricTensor[prevI][j][nextK].matrix[1][2] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, prevI, j, nextK, 2));

	//E i-1 j k-1
	elementX = cthetadt2 * ((0.125 * dielectricTensor[prevI][j][prevK].matrix[1][0] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, prevI, j, prevK, 0));
	elementY = cthetadt2 * ( + (0.125 * dielectricTensor[prevI][j][prevK].matrix[1][1] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, prevI, j, prevK, 1));
	elementZ = cthetadt2 * ((0.125 * dielectricTensor[prevI][j][prevK].matrix[1][2] / (deltaY * deltaY)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, prevI, j, prevK, 2));

	//E i+1 j+1 k+1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[1][0] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[0][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[2][0] / (deltaY * deltaZ)));
		elementY = cthetadt2 * ( - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[1][1] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[0][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[2][1] / (deltaY * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[1][2] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[0][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[2][2] / (deltaY * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i + 1, nextJ, nextK, 0));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i + 1, nextJ, nextK, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i + 1, nextJ, nextK, 2));
		} else {
			rightPart.y -= elementX * E0.x;
			rightPart.y -= elementY * E0.y;
			rightPart.y -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[1][0] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[0][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[2][0] / (deltaY * deltaZ)));
		elementY = cthetadt2 * ( - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[1][1] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[0][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[2][1] / (deltaY * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[1][2] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[0][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[2][2] / (deltaY * deltaZ)));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, nextI, nextJ, nextK, 0));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, nextI, nextJ, nextK, 1));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, nextI, nextJ, nextK, 2));
	}

	//E i+1 j+1 k-1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[1][0] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[0][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[2][0] / (deltaY * deltaZ)));
		elementY = cthetadt2 * ( - (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[1][1] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[0][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[2][1] / (deltaY * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[1][2] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[0][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[2][2] / (deltaY * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i + 1, nextJ, prevK, 0));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i + 1, nextJ, prevK, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i + 1, nextJ, prevK, 2));
		} else {
			rightPart.y -= elementX * E0.x;
			rightPart.y -= elementY * E0.y;
			rightPart.y -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[1][0] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[0][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[2][0] / (deltaY * deltaZ)));
		elementY = cthetadt2 * ( - (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[1][1] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[0][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[2][1] / (deltaY * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[1][2] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[0][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[2][2] / (deltaY * deltaZ)));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, nextI, nextJ, prevK, 0));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, nextI, nextJ, prevK, 1));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, nextI, nextJ, prevK, 2));
	}

	//E i+1 j-1 k+1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[1][0] / (deltaY * deltaY)) + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[0][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[2][0] / (deltaY * deltaZ)));
		elementY = cthetadt2 * ( + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[1][1] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[0][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[2][1] / (deltaY * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[1][2] / (deltaY * deltaY)) + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[0][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[2][2] / (deltaY * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i + 1, prevJ, nextK, 0));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i + 1, prevJ, nextK, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i + 1, prevJ, nextK, 2));
		} else {
			rightPart.y -= elementX * E0.x;
			rightPart.y -= elementY * E0.y;
			rightPart.y -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[1][0] / (deltaY * deltaY)) + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[0][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[2][0] / (deltaY * deltaZ)));
		elementY = cthetadt2 * ( + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[1][1] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[0][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[2][1] / (deltaY * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[1][2] / (deltaY * deltaY)) + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[0][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[2][2] / (deltaY * deltaZ)));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, nextI, prevJ, nextK, 0));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, nextI, prevJ, nextK, 1));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, nextI, prevJ, nextK, 2));
	}

	//E i+1 j-1 k-1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[1][0] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[0][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[2][0] / (deltaY * deltaZ)));
		elementY = cthetadt2 * ( - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[1][1] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[0][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[2][1] / (deltaY * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[1][2] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[0][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[2][2] / (deltaY * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, i + 1, prevJ, prevK, 0));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, i + 1, prevJ, prevK, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, i + 1, prevJ, prevK, 2));
		} else {
			rightPart.x -= elementX * E0.x;
			rightPart.x -= elementY * E0.y;
			rightPart.x -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[1][0] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[0][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[2][0] / (deltaY * deltaZ)));
		elementY = cthetadt2 * ( - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[1][1] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[0][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[2][1] / (deltaY * deltaZ)));
		elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[1][2] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[0][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[2][2] / (deltaY * deltaZ)));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, nextI, prevJ, prevK, 0));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, nextI, prevJ, prevK, 1));
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, nextI, prevJ, prevK, 2));
	}

	//E i-1 j+1 k+1
	elementX = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[1][0] / (deltaY * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[0][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[2][0] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, prevI, nextJ, nextK, 0));
	elementY = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[1][1] / (deltaY * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[0][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[2][1] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, prevI, nextJ, nextK, 1));
	elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[1][2] / (deltaY * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[0][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[2][2] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, prevI, nextJ, nextK, 2));

	//E i-1 j+1 k-1
	elementX = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[1][0] / (deltaY * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[0][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[2][0] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, prevI, nextJ, prevK, 0));
	elementY = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[1][1] / (deltaY * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[0][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[2][1] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, prevI, nextJ, prevK, 1));
	elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[1][2] / (deltaY * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[0][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[2][2] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, prevI, nextJ, prevK, 2));

	//E i-1 j-1 k+1
	elementX = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[1][0] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[0][0] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[2][0] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, prevI, prevJ, nextK, 0));
	elementY = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[1][1] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[0][1] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[2][1] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, prevI, prevJ, nextK, 1));
	elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[1][2] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[0][2] / (deltaX * deltaY)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[2][2] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, prevI, prevJ, nextK, 2));

	//E i-1 j-1 k-1
	elementX = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[1][0] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[0][0] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[2][0] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementX, prevI, prevJ, prevK, 0));
	elementY = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[1][1] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[0][1] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[2][1] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementY, prevI, prevJ, prevK, 1));
	elementZ = cthetadt2 * (-(0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[1][2] / (deltaY * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[0][2] / (deltaX * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[2][2] / (deltaY * deltaZ)));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(elementZ, prevI, prevJ, prevK, 2));

}

void Simulation::createInternalEquationZ(int i, int j, int k, Vector3d& rightPart) {
	int prevJ = j - 1;
	int prevK = k - 1;
	int nextJ = j + 1;
	int nextK = k + 1;
	int prevI = i - 1;
	int nextI = i + 1;

	if (prevI < 0) {
		prevI = xnumber - 1;
	}
	if (nextI >= xnumber) {
		nextI = 0;
	}
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	if (prevK < 0) {
		prevK = znumber - 1;
	}
	if (nextJ >= ynumber) {
		nextJ = 0;
	}
	if (nextK >= znumber) {
		nextK = 0;
	}
	double cthetadt2 = sqr(speed_of_light_normalized * theta * deltaT);

	double elementX;
	double elementY;
	double elementZ;

	//E i j k
	elementX = dielectricTensor[i][j][k].matrix[2][0] + cthetadt2 * (0.5 * dielectricTensor[i][j][k].matrix[2][0] / (deltaZ * deltaZ));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i, j, k, 0));
	elementY = dielectricTensor[i][j][k].matrix[2][1] + cthetadt2 * (0.5 * dielectricTensor[i][j][k].matrix[2][1] / (deltaZ * deltaZ));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i, j, k, 1));
	elementZ = 1 + dielectricTensor[i][j][k].matrix[2][2] + cthetadt2 * ((2 / (deltaX * deltaX) + (2 / (deltaY * deltaY))) + (2 / (deltaZ * deltaZ)) + (0.5 * dielectricTensor[i][j][k].matrix[2][2] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i, j, k, 2));

	//E i j k+1
	elementX = cthetadt2 * (- (0.25 * dielectricTensor[i][j][nextK].matrix[2][0] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i, j, nextK, 0));
	elementY = cthetadt2 * (- (0.25 * dielectricTensor[i][j][nextK].matrix[2][1] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i, j, nextK, 1));
	elementZ = cthetadt2 * ((1 / (deltaZ * deltaZ)) - (0.25 * dielectricTensor[i][j][nextK].matrix[2][2] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i, j, nextK, 2));

	//E i j k-1
	elementX = cthetadt2 * (- (0.25 * dielectricTensor[i][j][prevK].matrix[2][0] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i, j, prevK, 0));
	elementY = cthetadt2 * (- (0.25 * dielectricTensor[i][j][prevK].matrix[2][1] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i, j, prevK, 1));
	elementZ = cthetadt2 * ((1 / (deltaZ * deltaZ)) - (0.25 * dielectricTensor[i][j][prevK].matrix[2][2] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i, j, prevK, 2));

	//E i+1 j k
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ((0.25 * dielectricTensor[i + 1][j][k].matrix[2][0] / (deltaZ * deltaZ)));
		elementY = cthetadt2 * ((0.25 * dielectricTensor[i + 1][j][k].matrix[2][1] / (deltaZ * deltaZ)));
		elementZ = cthetadt2 * ((-1 / (deltaX * deltaX)) + (0.25 * dielectricTensor[i + 1][j][k].matrix[2][2] / (deltaZ * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i + 1, j, k, 0));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i + 1, j, k, 1));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i + 1, j, k, 2));
		} else {
			rightPart.z -= elementX * E0.x;
			rightPart.z -= elementY * E0.y;
			rightPart.z -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ((0.25 * dielectricTensor[nextI][j][k].matrix[2][0] / (deltaZ * deltaZ)));
		elementY = cthetadt2 * ((0.25 * dielectricTensor[nextI][j][k].matrix[2][1] / (deltaZ * deltaZ)));
		elementZ = cthetadt2 * ((-1 / (deltaX * deltaX)) + (0.25 * dielectricTensor[nextI][j][k].matrix[2][2] / (deltaZ * deltaZ)));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, nextI, j, k, 0));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, nextI, j, k, 1));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, nextI, j, k, 2));
	}

	//E i-1 j k
	elementX = cthetadt2 * ((0.25 * dielectricTensor[prevI][j][k].matrix[2][0] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, prevI, j, k, 0));
	elementY = cthetadt2 * ((0.25 * dielectricTensor[prevI][j][k].matrix[2][1] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, prevI, j, k, 1));
	elementZ = cthetadt2 * ((-1 / (deltaX * deltaX)) + (0.25 * dielectricTensor[prevI][j][k].matrix[2][2] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, prevI, j, k, 2));

	//E i j+1 k
	elementX = cthetadt2 * ((0.25 * dielectricTensor[i][nextJ][k].matrix[2][0]) / (deltaZ * deltaZ));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i, nextJ, k, 0));
	elementY = cthetadt2 * ((0.25 * dielectricTensor[i][nextJ][k].matrix[2][1]) / (deltaZ * deltaZ));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i, nextJ, k, 1));
	elementZ = cthetadt2 * ( - (1 / (deltaY * deltaY)) + (0.25 * dielectricTensor[i][nextJ][k].matrix[2][2] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i, nextJ, k, 2));

	//E i j-1 k
	elementX = cthetadt2 * ((0.25 * dielectricTensor[i][prevJ][k].matrix[2][0]) / (deltaZ * deltaZ));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i, prevJ, k, 0));
	elementY = cthetadt2 * ((0.25 * dielectricTensor[i][prevJ][k].matrix[2][1]) / (deltaZ * deltaZ));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i, prevJ, k, 1));
	elementZ = cthetadt2 * ( - (1 / (deltaY * deltaY)) + (0.25 * dielectricTensor[i][prevJ][k].matrix[2][2] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i, prevJ, k, 2));

	//E i+1 j k+1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][j][nextK].matrix[2][0] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[i + 1][j][nextK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][j][nextK].matrix[2][1] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[i + 1][j][nextK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * ( - (0.125 * dielectricTensor[i + 1][j][nextK].matrix[2][2] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[i + 1][j][nextK].matrix[0][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i + 1, j, nextK, 0));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i + 1, j, nextK, 1));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i + 1, j, nextK, 2));
		} else {
			rightPart.z -= elementX * E0.x;
			rightPart.z -= elementY * E0.y;
			rightPart.z -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (- (0.125 * dielectricTensor[nextI][j][nextK].matrix[2][0] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[nextI][j][nextK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[nextI][j][nextK].matrix[2][1] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[nextI][j][nextK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * (- (0.125 * dielectricTensor[nextI][j][nextK].matrix[2][2] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[nextI][j][nextK].matrix[0][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, nextI, j, nextK, 0));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, nextI, j, nextK, 1));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, nextI, j, nextK, 2));
	}

	//E i+1 j k-1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][j][prevK].matrix[2][0] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[i + 1][j][prevK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[i + 1][j][prevK].matrix[2][1] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[i + 1][j][prevK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * ( - (0.125 * dielectricTensor[i + 1][j][prevK].matrix[2][2] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[i + 1][j][prevK].matrix[0][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i + 1, j, prevK, 0));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i + 1, j, prevK, 1));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i + 1, j, prevK, 2));
		} else {
			rightPart.z -= elementX * E0.x;
			rightPart.z -= elementY * E0.y;
			rightPart.z -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (- (0.125 * dielectricTensor[nextI][j][prevK].matrix[2][0] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[nextI][j][prevK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.125 * dielectricTensor[nextI][j][prevK].matrix[2][1] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[nextI][j][prevK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * ( - (0.125 * dielectricTensor[nextI][j][prevK].matrix[2][2] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[nextI][j][prevK].matrix[0][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, nextI, j, prevK, 0));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, nextI, j, prevK, 1));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, nextI, j, prevK, 2));
	}

	//E i-1 j k+1
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[prevI][j][nextK].matrix[2][0] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[prevI][j][nextK].matrix[0][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, prevI, j, nextK, 0));
	elementY = cthetadt2 * (- (0.125 * dielectricTensor[prevI][j][nextK].matrix[2][1] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[prevI][j][nextK].matrix[0][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, prevI, j, nextK, 1));
	elementZ = cthetadt2 * ( - (0.125 * dielectricTensor[prevI][j][nextK].matrix[2][2] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[prevI][j][nextK].matrix[0][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, prevI, j, nextK, 2));

	//E i-1 j k-1
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[prevI][j][prevK].matrix[2][0] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[prevI][j][prevK].matrix[0][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, prevI, j, prevK, 0));
	elementY = cthetadt2 * (- (0.125 * dielectricTensor[prevI][j][prevK].matrix[2][1] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[prevI][j][prevK].matrix[0][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, prevI, j, prevK, 1));
	elementZ = cthetadt2 * ( - (0.125 * dielectricTensor[prevI][j][prevK].matrix[2][2] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[prevI][j][prevK].matrix[0][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, prevI, j, prevK, 2));

	//E i j+1 k+1
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[i][nextJ][nextK].matrix[2][0] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[i][nextJ][nextK].matrix[1][0] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i, nextJ, nextK, 0));
	elementY = cthetadt2 * (- (0.125 * dielectricTensor[i][nextJ][nextK].matrix[2][1] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[i][nextJ][nextK].matrix[1][1] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i, nextJ, nextK, 1));
	elementZ = cthetadt2 * ( - (0.125 * dielectricTensor[i][nextJ][nextK].matrix[2][2] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[i][nextJ][nextK].matrix[1][2] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i, nextJ, nextK, 2));

	//E i j+1 k-1
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[i][nextJ][prevK].matrix[2][0] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[i][nextJ][prevK].matrix[1][0] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i, nextJ, prevK, 0));
	elementY = cthetadt2 * (- (0.125 * dielectricTensor[i][nextJ][prevK].matrix[2][1] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[i][nextJ][prevK].matrix[1][1] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i, nextJ, prevK, 1));
	elementZ = cthetadt2 * ( - (0.125 * dielectricTensor[i][nextJ][prevK].matrix[2][2] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[i][nextJ][prevK].matrix[1][2] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i, nextJ, prevK, 2));

	//E i j-1 k+1
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[i][prevJ][nextK].matrix[2][0] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[i][prevJ][nextK].matrix[1][0] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i, prevJ, nextK, 0));
	elementY = cthetadt2 * (- (0.125 * dielectricTensor[i][prevJ][nextK].matrix[2][1] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[i][prevJ][nextK].matrix[1][1] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i, prevJ, nextK, 1));
	elementZ = cthetadt2 * ( - (0.125 * dielectricTensor[i][prevJ][nextK].matrix[2][2] / (deltaZ * deltaZ)) + (0.125 * dielectricTensor[i][prevJ][nextK].matrix[1][2] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i, prevJ, nextK, 2));

	//E i j-1 k-1
	elementX = cthetadt2 * (- (0.125 * dielectricTensor[i][prevJ][prevK].matrix[2][0] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[i][prevJ][prevK].matrix[1][0] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i, prevJ, prevK, 0));
	elementY = cthetadt2 * (- (0.125 * dielectricTensor[i][prevJ][prevK].matrix[2][1] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[i][prevJ][prevK].matrix[1][1] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i, prevJ, prevK, 1));
	elementZ = cthetadt2 * ( - (0.125 * dielectricTensor[i][prevJ][prevK].matrix[2][2] / (deltaZ * deltaZ)) - (0.125 * dielectricTensor[i][prevJ][prevK].matrix[1][2] / (deltaZ * deltaY)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i, prevJ, prevK, 2));

	//E i+1 j+1 k
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ((0.125 * dielectricTensor[i + 1][nextJ][k].matrix[2][0] / (deltaZ * deltaZ)));
		elementY = cthetadt2 * ((0.125 * dielectricTensor[i + 1][nextJ][k].matrix[2][1] / (deltaZ * deltaZ)));
		elementZ = cthetadt2 * ( + (0.125 * dielectricTensor[i + 1][nextJ][k].matrix[2][2] / (deltaZ * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i + 1, nextJ, k, 0));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i + 1, nextJ, k, 1));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i + 1, nextJ, k, 2));
		} else {
			rightPart.z -= elementX * E0.x;
			rightPart.z -= elementY * E0.y;
			rightPart.z -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ((0.125 * dielectricTensor[nextI][nextJ][k].matrix[2][0] / (deltaZ * deltaZ)));
		elementY = cthetadt2 * ((0.125 * dielectricTensor[nextI][nextJ][k].matrix[2][1] / (deltaZ * deltaZ)));
		elementZ = cthetadt2 * ( + (0.125 * dielectricTensor[nextI][nextJ][k].matrix[2][2] / (deltaZ * deltaZ)));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, nextI, nextJ, k, 0));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, nextI, nextJ, k, 1));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, nextI, nextJ, k, 2));
	}

	//E i+1 j-1 k
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * ((0.125 * dielectricTensor[i + 1][prevJ][k].matrix[2][0] / (deltaZ * deltaZ)));
		elementY = cthetadt2 * ((0.125 * dielectricTensor[i + 1][prevJ][k].matrix[2][1] / (deltaZ * deltaZ)));
		elementZ = cthetadt2 * ( + (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[2][2] / (deltaZ * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i + 1, prevJ, k, 0));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i + 1, prevJ, k, 1));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i + 1, prevJ, k, 2));
		} else {
			rightPart.z -= elementX * E0.x;
			rightPart.z -= elementY * E0.y;
			rightPart.z -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * ((0.125 * dielectricTensor[i + 1][prevJ][k].matrix[2][0] / (deltaZ * deltaZ)));
		elementY = cthetadt2 * ((0.125 * dielectricTensor[i + 1][prevJ][k].matrix[2][1] / (deltaZ * deltaZ)));
		elementZ = cthetadt2 * ( + (0.125 * dielectricTensor[i + 1][prevJ][k].matrix[2][2] / (deltaZ * deltaZ)));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, nextI, prevJ, k, 0));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, nextI, prevJ, k, 1));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, nextI, prevJ, k, 2));
	}

	//E i-1 j+1 k
	elementX = cthetadt2 * ((0.125 * dielectricTensor[prevI][nextJ][k].matrix[2][0] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, prevI, nextJ, k, 0));
	elementY = cthetadt2 * ((0.125 * dielectricTensor[prevI][nextJ][k].matrix[2][1] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, prevI, nextJ, k, 1));
	elementZ = cthetadt2 * ( + (0.125 * dielectricTensor[prevI][nextJ][k].matrix[2][2] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, prevI, nextJ, k, 2));

	//E i-1 j-1 k
	elementX = cthetadt2 * ((0.125 * dielectricTensor[prevI][prevJ][k].matrix[2][0] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, prevI, prevJ, k, 0));
	elementY = cthetadt2 * ((0.125 * dielectricTensor[prevI][prevJ][k].matrix[2][1] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, prevI, prevJ, k, 1));
	elementZ = cthetadt2 * ( + (0.125 * dielectricTensor[prevI][prevJ][k].matrix[2][2] / (deltaZ * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, prevI, prevJ, k, 2));

	//E i+1 j+1 k+1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (- (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[2][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[1][0] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[2][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[1][1] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[2][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[1][2] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[i + 1][nextJ][nextK].matrix[0][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i + 1, nextJ, nextK, 0));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i + 1, nextJ, nextK, 1));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i + 1, nextJ, nextK, 2));
		} else {
			rightPart.z -= elementX * E0.x;
			rightPart.z -= elementY * E0.y;
			rightPart.z -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (- (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[2][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[1][0] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[2][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[1][1] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[2][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[1][2] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[nextI][nextJ][nextK].matrix[0][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, nextI, nextJ, nextK, 0));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, nextI, nextJ, nextK, 1));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, nextI, nextJ, nextK, 2));
	}

	//E i+1 j+1 k-1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (- (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[2][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[1][0] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[2][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[1][1] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[2][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[1][2] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[i + 1][nextJ][prevK].matrix[0][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i + 1, nextJ, prevK, 0));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i + 1, nextJ, prevK, 1));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i + 1, nextJ, prevK, 2));
		} else {
			rightPart.z -= elementX * E0.x;
			rightPart.z -= elementY * E0.y;
			rightPart.z -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (- (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[2][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[1][0] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[2][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[1][1] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[2][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[1][2] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[nextI][nextJ][prevK].matrix[0][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, nextI, nextJ, prevK, 0));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, nextI, nextJ, prevK, 1));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, nextI, nextJ, prevK, 2));
	}

	//E i+1 j-1 k+1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (- (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[2][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[1][0] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[2][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[1][1] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[2][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[1][2] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[i + 1][prevJ][nextK].matrix[0][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i + 1, prevJ, nextK, 0));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i + 1, prevJ, nextK, 1));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i + 1, prevJ, nextK, 2));
		} else {
			rightPart.z -= elementX * E0.x;
			rightPart.z -= elementY * E0.y;
			rightPart.z -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (- (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[2][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[1][0] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[2][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[1][1] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[2][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[1][2] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[nextI][prevJ][nextK].matrix[0][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, nextI, prevJ, nextK, 0));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, nextI, prevJ, nextK, 1));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, nextI, prevJ, nextK, 2));
	}

	//E i+1 j-1 k-1
	if (boundaryConditionType == SUPERCONDUCTERLEFT) {
		elementX = cthetadt2 * (- (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[2][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[1][0] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[2][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[1][1] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[2][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[1][2] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[i + 1][prevJ][prevK].matrix[0][2] / (deltaX * deltaZ)));
		if (i < xnumber - 1) {
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, i + 1, prevJ, prevK, 0));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, i + 1, prevJ, prevK, 1));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, i + 1, prevJ, prevK, 2));
		} else {
			rightPart.z -= elementX * E0.x;
			rightPart.z -= elementY * E0.y;
			rightPart.z -= elementZ * E0.z;
		}
	}
	if (boundaryConditionType == PERIODIC) {
		elementX = cthetadt2 * (- (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[2][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[1][0] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[0][0] / (deltaX * deltaZ)));
		elementY = cthetadt2 * (- (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[2][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[1][1] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[0][1] / (deltaX * deltaZ)));
		elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[2][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[1][2] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[nextI][prevJ][prevK].matrix[0][2] / (deltaX * deltaZ)));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, nextI, prevJ, prevK, 0));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, nextI, prevJ, prevK, 1));
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, nextI, prevJ, prevK, 2));
	}

	//E i-1 j+1 k+1
	elementX = cthetadt2 * (- (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[2][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[1][0] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[0][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, prevI, nextJ, nextK, 0));
	elementY = cthetadt2 * (- (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[2][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[1][1] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[0][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, prevI, nextJ, nextK, 1));
	elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[2][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[1][2] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[prevI][nextJ][nextK].matrix[0][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, prevI, nextJ, nextK, 2));

	//E i-1 j+1 k-1
	elementX = cthetadt2 * (- (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[2][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[1][0] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[0][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, prevI, nextJ, prevK, 0));
	elementY = cthetadt2 * (- (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[2][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[1][1] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[0][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, prevI, nextJ, prevK, 1));
	elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[2][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[1][2] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[prevI][nextJ][prevK].matrix[0][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, prevI, nextJ, prevK, 2));

	//E i-1 j-1 k+1
	elementX = cthetadt2 * (- (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[2][0] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[1][0] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[0][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, prevI, prevJ, nextK, 0));
	elementY = cthetadt2 * (- (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[2][1] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[1][1] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[0][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, prevI, prevJ, nextK, 1));
	elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[2][2] / (deltaX * deltaX)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[1][2] / (deltaZ * deltaY)) + (0.0625 * dielectricTensor[prevI][prevJ][nextK].matrix[0][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, prevI, prevJ, nextK, 2));

	//E i-1 j-1 k-1
	elementX = cthetadt2 * (- (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[2][0] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[1][0] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[0][0] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementX, prevI, prevJ, prevK, 0));
	elementY = cthetadt2 * (- (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[2][1] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[1][1] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[0][1] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementY, prevI, prevJ, prevK, 1));
	elementZ = cthetadt2 * ( - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[2][2] / (deltaX * deltaX)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[1][2] / (deltaZ * deltaY)) - (0.0625 * dielectricTensor[prevI][prevJ][prevK].matrix[0][2] / (deltaX * deltaZ)));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(elementZ, prevI, prevJ, prevK, 2));
}

void Simulation::evaluateMagneticField() {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Vector3d rotEold = evaluateRotE(i, j, k);
				Vector3d rotEnew = evaluateRotNewE(i, j, k);
				newBfield[i][j][k] = Bfield[i][j][k] - (rotEold*(1 - theta) + rotEnew*theta) * speed_of_light_normalized * deltaT;
			}
		}
	}
}

void Simulation::updateBoundaries() {
	for (int i = 0; i < xnumber; ++i) {
		//periodic by Y
		for (int k = 0; k < znumber; ++k) {
			tempEfield[i][ynumber][k] = tempEfield[i][0][k];
		}
		//periodic by Z
		//note j <= number because corner point

		for (int j = 0; j <= ynumber; ++j) {
			tempEfield[i][j][znumber] = tempEfield[i][j][0];
		}
	}

	if(boundaryConditionType == PERIODIC) {
		for(int j = 0; j <= ynumber; ++j) {
			for(int k = 0; k <= znumber; ++k) {
				tempEfield[xnumber][j][k] = tempEfield[0][j][k];
			}
		}
	}
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		for (int j = 0; j <= ynumber; ++j) {
			for (int k = 0; k <= znumber; ++k) {
				tempEfield[xnumber][j][k] = E0;
			}
		}
	}
}

void Simulation::updateBoundariesOldField() {
	for (int i = 0; i < xnumber; ++i) {
		//periodic by Y
		for (int k = 0; k < znumber; ++k) {
			Efield[i][ynumber][k] = Efield[i][0][k];
		}
		//periodic by Z
		//note j <= number because corner point

		for (int j = 0; j <= ynumber; ++j) {
			Efield[i][j][znumber] = Efield[i][j][0];
		}
	}

	if(boundaryConditionType == PERIODIC) {
		for(int j = 0; j <= ynumber; ++j) {
			for(int k = 0; k <= znumber; ++k) {
				Efield[xnumber][j][k] = Efield[0][j][k];
			}
		}
	}
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		for (int j = 0; j <= ynumber; ++j) {
			for (int k = 0; k <= znumber; ++k) {
				Efield[xnumber][j][k] = E0;
			}
		}
	}
}

void Simulation::updateBoundariesNewField() {
	for (int i = 0; i < xnumber; ++i) {
		//periodic by Y
		for (int k = 0; k < znumber; ++k) {
			newEfield[i][ynumber][k] = newEfield[i][0][k];
		}
		//periodic by Z
		//note j <= number because corner point

		for (int j = 0; j <= ynumber; ++j) {
			newEfield[i][j][znumber] = newEfield[i][j][0];
		}
	}

	if(boundaryConditionType == PERIODIC) {
		for(int j = 0; j <= ynumber; ++j) {
			for(int k = 0; k <= znumber; ++k) {
				newEfield[xnumber][j][k] = newEfield[0][j][k];
			}
		}
	}
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		for (int j = 0; j <= ynumber; ++j) {
			for (int k = 0; k <= znumber; ++k) {
				newEfield[xnumber][j][k] = E0;
			}
		}
	}
}

double Simulation::evaluateDivFlux(int i, int j, int k) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			exit(0);
		}

		if (j < 0) {
			printf("j < 0\n");
			exit(0);
		}

		if (j >= ynumber) {
			printf("j >= ynumber\n");
			exit(0);
		}

		if (k < 0) {
			printf("k < 0\n");
			exit(0);
		}

		if (k >= znumber) {
			printf("k >= znumber\n");
			exit(0);
		}
	}


	double rfluxX = (electricFlux[i + 1][j][k].x + electricFlux[i + 1][j + 1][k].x + electricFlux[i + 1][j][k + 1].x + electricFlux[i + 1][j + 1][k + 1].x) / 4;
	double lfluxX = (electricFlux[i][j][k].x + electricFlux[i][j + 1][k].x + electricFlux[i][j][k + 1].x + electricFlux[i][j + 1][k + 1].x) / 4;

	double rfluxY = (electricFlux[i][j + 1][k].y + electricFlux[i + 1][j + 1][k].y + electricFlux[i][j + 1][k + 1].y + electricFlux[i + 1][j + 1][k + 1].y) / 4;
	double lfluxY = (electricFlux[i][j][k].y + electricFlux[i + 1][j][k].y + electricFlux[i][j][k + 1].y + electricFlux[i + 1][j][k + 1].y) / 4;

	double rfluxZ = (electricFlux[i][j][k + 1].z + electricFlux[i + 1][j][k + 1].z + electricFlux[i][j + 1][k + 1].z + electricFlux[i + 1][j + 1][k + 1].z) / 4;
	double lfluxZ = (electricFlux[i][j][k].z + electricFlux[i + 1][j][k].z + electricFlux[i][j + 1][k].z + electricFlux[i + 1][j + 1][k].z) / 4;

	return ((rfluxX - lfluxX) / deltaX) + ((rfluxY - lfluxY) / deltaY) + ((rfluxZ - lfluxZ) / deltaZ);
}

Vector3d Simulation::evaluateRotB(int i, int j, int k) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i > xnumber) {
			printf("x > xnumber\n");
			exit(0);
		}

		if (j < 0) {
			printf("j < 0\n");
			exit(0);
		}

		if (j > ynumber) {
			printf("j > ynumber\n");
			exit(0);
		}

		if (k < 0) {
			printf("k < 0\n");
			exit(0);
		}

		if (k > znumber) {
			printf("k > znumber\n");
			exit(0);
		}
	}

	Vector3d BrightX;
	Vector3d BleftX;
	Vector3d BrightY;
	Vector3d BleftY;
	Vector3d BrightZ;
	Vector3d BleftZ;

	if (i == 0) {
		if(boundaryConditionType == PERIODIC) {
			int prevI = i - 1;
			if(prevI < 0) {
				prevI = xnumber - 1;
			}
			BrightX = (getBfield(i, j, k) + getBfield(i, j - 1, k) + getBfield(i, j - 1, k - 1) + getBfield(i, j - 1, k - 1)) / 4;
			BleftX = (getBfield(prevI, j, k) + getBfield(prevI, j - 1, k) + getBfield(prevI, j - 1, k - 1) + getBfield(prevI, j - 1, k - 1)) / 4;

			BrightY = (getBfield(i, j, k) + getBfield(prevI, j, k) + getBfield(i, j, k - 1) + getBfield(prevI, j, k - 1)) / 4;
			BleftY = (getBfield(i, j - 1, k) + getBfield(prevI, j - 1, k) + getBfield(i, j - 1, k - 1) + getBfield(prevI, j - 1, k - 1)) / 4;

			BrightZ = (getBfield(i, j, k) + getBfield(prevI, j, k) + getBfield(i, j - 1, k) + getBfield(prevI, j - 1, k)) / 4;
			BleftZ = (getBfield(i, j, k - 1) + getBfield(prevI, j, k - 1) + getBfield(i, j - 1, k - 1) + getBfield(prevI, j - 1, k - 1)) / 4;
		}
	} else {
		BrightX = (getBfield(i, j, k) + getBfield(i, j - 1, k) + getBfield(i, j - 1, k - 1) + getBfield(i, j - 1, k - 1)) / 4;
		BleftX = (getBfield(i - 1, j, k) + getBfield(i - 1, j - 1, k) + getBfield(i - 1, j - 1, k - 1) + getBfield(i - 1, j - 1, k - 1)) / 4;

		BrightY = (getBfield(i, j, k) + getBfield(i - 1, j, k) + getBfield(i, j, k - 1) + getBfield(i - 1, j, k - 1)) / 4;
		BleftY = (getBfield(i, j - 1, k) + getBfield(i - 1, j - 1, k) + getBfield(i, j - 1, k - 1) + getBfield(i - 1, j - 1, k - 1)) / 4;

		BrightZ = (getBfield(i, j, k) + getBfield(i - 1, j, k) + getBfield(i, j - 1, k) + getBfield(i - 1, j - 1, k)) / 4;
		BleftZ = (getBfield(i, j, k - 1) + getBfield(i - 1, j, k - 1) + getBfield(i, j - 1, k - 1) + getBfield(i - 1, j - 1, k - 1)) / 4;
	}


	double x = 0;
	double y = 0;
	double z = 0;

	x = ((BrightY.z - BleftY.z) / deltaY) - ((BrightZ.y - BleftZ.y) / deltaZ);
	y = ((BrightZ.x - BleftZ.x) / deltaZ) - ((BrightX.z - BleftX.z) / deltaX);
	z = ((BrightX.y - BleftX.y) / deltaX) - ((BrightY.x - BleftY.x) / deltaY);

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotTempE(int i, int j, int k) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			exit(0);
		}

		if (j < 0) {
			printf("j < 0\n");
			exit(0);
		}

		if (j >= ynumber) {
			printf("j >= ynumber\n");
			exit(0);
		}

		if (k < 0) {
			printf("k < 0\n");
			exit(0);
		}

		if (k >= znumber) {
			printf("k >= znumber\n");
			exit(0);
		}
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = (tempEfield[i + 1][j][k] + tempEfield[i + 1][j + 1][k] + tempEfield[i + 1][j][k + 1] + tempEfield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftX = (tempEfield[i][j][k] + tempEfield[i][j + 1][k] + tempEfield[i][j][k + 1] + tempEfield[i][j + 1][k + 1]) / 4;

	Vector3d ErightY = (tempEfield[i][j + 1][k] + tempEfield[i + 1][j + 1][k] + tempEfield[i][j + 1][k + 1] + tempEfield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftY = (tempEfield[i][j][k] + tempEfield[i + 1][j][k] + tempEfield[i][j][k + 1] + tempEfield[i + 1][j][k + 1]) / 4;

	Vector3d ErightZ = (tempEfield[i][j][k + 1] + tempEfield[i + 1][j][k + 1] + tempEfield[i][j + 1][k + 1] + tempEfield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftZ = (tempEfield[i][j][k] + tempEfield[i + 1][j][k] + tempEfield[i][j + 1][k] + tempEfield[i + 1][j + 1][k]) / 4;

	x = ((ErightY.z - EleftY.z) / deltaY) - ((ErightZ.y - EleftZ.y) / deltaZ);
	y = ((ErightZ.x - EleftZ.x) / deltaZ) - ((ErightX.z - EleftX.z) / deltaX);
	z = ((ErightX.y - EleftX.y) / deltaX) - ((ErightY.x - EleftY.x) / deltaY);

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotE(int i, int j, int k) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			exit(0);
		}

		if (j < 0) {
			printf("j < 0\n");
			exit(0);
		}

		if (j >= ynumber) {
			printf("j >= ynumber\n");
			exit(0);
		}

		if (k < 0) {
			printf("k < 0\n");
			exit(0);
		}

		if (k >= znumber) {
			printf("k >= znumber\n");
			exit(0);
		}
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = (Efield[i + 1][j][k] + Efield[i + 1][j + 1][k] + Efield[i + 1][j][k + 1] + Efield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftX = (Efield[i][j][k] + Efield[i][j + 1][k] + Efield[i][j][k + 1] + Efield[i][j + 1][k + 1]) / 4;

	Vector3d ErightY = (Efield[i][j + 1][k] + Efield[i + 1][j + 1][k] + Efield[i][j + 1][k + 1] + Efield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftY = (Efield[i][j][k] + Efield[i + 1][j][k] + Efield[i][j][k + 1] + Efield[i + 1][j][k + 1]) / 4;

	Vector3d ErightZ = (Efield[i][j][k + 1] + Efield[i + 1][j][k + 1] + Efield[i][j + 1][k + 1] + Efield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftZ = (Efield[i][j][k] + Efield[i + 1][j][k] + Efield[i][j + 1][k] + Efield[i + 1][j + 1][k]) / 4;

	x = ((ErightY.z - EleftY.z) / deltaY) - ((ErightZ.y - EleftZ.y) / deltaZ);
	y = ((ErightZ.x - EleftZ.x) / deltaZ) - ((ErightX.z - EleftX.z) / deltaX);
	z = ((ErightX.y - EleftX.y) / deltaX) - ((ErightY.x - EleftY.x) / deltaY);

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotNewE(int i, int j, int k) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			exit(0);
		}

		if (j < 0) {
			printf("j < 0\n");
			exit(0);
		}

		if (j >= ynumber) {
			printf("j >= ynumber\n");
			exit(0);
		}

		if (k < 0) {
			printf("k < 0\n");
			exit(0);
		}

		if (k >= znumber) {
			printf("k >= znumber\n");
			exit(0);
		}
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = (newEfield[i + 1][j][k] + newEfield[i + 1][j + 1][k] + newEfield[i + 1][j][k + 1] + newEfield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftX = (newEfield[i][j][k] + newEfield[i][j + 1][k] + newEfield[i][j][k + 1] + newEfield[i][j + 1][k + 1]) / 4;

	Vector3d ErightY = (newEfield[i][j + 1][k] + newEfield[i + 1][j + 1][k] + newEfield[i][j + 1][k + 1] + newEfield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftY = (newEfield[i][j][k] + newEfield[i + 1][j][k] + newEfield[i][j][k + 1] + newEfield[i + 1][j][k + 1]) / 4;

	Vector3d ErightZ = (newEfield[i][j][k + 1] + newEfield[i + 1][j][k + 1] + newEfield[i][j + 1][k + 1] + newEfield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftZ = (newEfield[i][j][k] + newEfield[i + 1][j][k] + newEfield[i][j + 1][k] + newEfield[i + 1][j + 1][k]) / 4;

	x = ((ErightY.z - EleftY.z) / deltaY) - ((ErightZ.y - EleftZ.y) / deltaZ);
	y = ((ErightZ.x - EleftZ.x) / deltaZ) - ((ErightX.z - EleftX.z) / deltaX);
	z = ((ErightX.y - EleftX.y) / deltaX) - ((ErightY.x - EleftY.x) / deltaY);

	return Vector3d(x, y, z);
}

double Simulation::evaluateDivE(int i, int j, int k) {
	double ErightX = (Efield[i + 1][j][k].x + Efield[i + 1][j + 1][k].x + Efield[i + 1][j][k + 1].x + Efield[i + 1][j + 1][k + 1].x) / 4;
	double EleftX = (Efield[i][j][k].x + Efield[i][j + 1][k].x + Efield[i][j][k + 1].x + Efield[i][j + 1][k + 1].x) / 4;

	double ErightY = (Efield[i][j + 1][k].y + Efield[i + 1][j + 1][k].y + Efield[i][j + 1][k + 1].y + Efield[i + 1][j + 1][k + 1].y) / 4;
	double EleftY = (Efield[i][j][k].y + Efield[i + 1][j][k].y + Efield[i][j][k + 1].y + Efield[i + 1][j][k + 1].y) / 4;

	double ErightZ = (Efield[i][j][k + 1].z + Efield[i + 1][j][k + 1].z + Efield[i][j + 1][k + 1].z + Efield[i + 1][j + 1][k + 1].z) / 4;
	double EleftZ = (Efield[i][j][k].z + Efield[i + 1][j][k].z + Efield[i][j + 1][k].z + Efield[i + 1][j + 1][k].z) / 4;

	return ((ErightX - EleftX) / deltaX) + ((ErightY - EleftY) / deltaY) + ((ErightZ - EleftZ) / deltaZ);
}

double Simulation::evaluateDivCleaningE(int i, int j, int k) {
	double ErightX = (divergenceCleaningField[i + 1][j][k][0] + divergenceCleaningField[i + 1][j + 1][k][0] + divergenceCleaningField[i + 1][j][k + 1][0] + divergenceCleaningField[i + 1][j + 1][k + 1][0]) / 4;
	double EleftX = (divergenceCleaningField[i][j][k][0] + divergenceCleaningField[i][j + 1][k][0] + divergenceCleaningField[i][j][k + 1][0] + divergenceCleaningField[i][j + 1][k + 1][0]) / 4;

	double ErightY = (divergenceCleaningField[i][j + 1][k][1] + divergenceCleaningField[i + 1][j + 1][k][1] + divergenceCleaningField[i][j + 1][k + 1][1] + divergenceCleaningField[i + 1][j + 1][k + 1][1]) / 4;
	double EleftY = (divergenceCleaningField[i][j][k][1] + divergenceCleaningField[i + 1][j][k][1] + divergenceCleaningField[i][j][k + 1][1] + divergenceCleaningField[i + 1][j][k + 1][1]) / 4;

	double ErightZ = (divergenceCleaningField[i][j][k + 1][2] + divergenceCleaningField[i + 1][j][k + 1][2] + divergenceCleaningField[i][j + 1][k + 1][2] + divergenceCleaningField[i + 1][j + 1][k + 1][2]) / 4;
	double EleftZ = (divergenceCleaningField[i][j][k][2] + divergenceCleaningField[i + 1][j][k][2] + divergenceCleaningField[i][j + 1][k][2] + divergenceCleaningField[i + 1][j + 1][k][2]) / 4;

	return ((ErightX - EleftX) / deltaX) + ((ErightY - EleftY) / deltaY) + ((ErightZ - EleftZ) / deltaZ);
}

double Simulation::evaluateDivTempE(int i, int j, int k) {
	double ErightX = (tempEfield[i + 1][j][k].x + tempEfield[i + 1][j + 1][k].x + tempEfield[i + 1][j][k + 1].x + tempEfield[i + 1][j + 1][k + 1].x) / 4;
	double EleftX = (tempEfield[i][j][k].x + tempEfield[i][j + 1][k].x + tempEfield[i][j][k + 1].x + tempEfield[i][j + 1][k + 1].x) / 4;

	double ErightY = (tempEfield[i][j + 1][k].y + tempEfield[i + 1][j + 1][k].y + tempEfield[i][j + 1][k + 1].y + tempEfield[i + 1][j + 1][k + 1].y) / 4;
	double EleftY = (tempEfield[i][j][k].y + tempEfield[i + 1][j][k].y + tempEfield[i][j][k + 1].y + tempEfield[i + 1][j][k + 1].y) / 4;

	double ErightZ = (tempEfield[i][j][k + 1].z + tempEfield[i + 1][j][k + 1].z + tempEfield[i][j + 1][k + 1].z + tempEfield[i + 1][j + 1][k + 1].z) / 4;
	double EleftZ = (tempEfield[i][j][k].z + tempEfield[i + 1][j][k].z + tempEfield[i][j + 1][k].z + tempEfield[i + 1][j + 1][k].z) / 4;

	return ((ErightX - EleftX) / deltaX) + ((ErightY - EleftY) / deltaY) + ((ErightZ - EleftZ) / deltaZ);
}

double Simulation::evaluateDivNewE(int i, int j, int k) {
	double ErightX = (newEfield[i + 1][j][k].x + newEfield[i + 1][j + 1][k].x + newEfield[i + 1][j][k + 1].x + newEfield[i + 1][j + 1][k + 1].x) / 4;
	double EleftX = (newEfield[i][j][k].x + newEfield[i][j + 1][k].x + newEfield[i][j][k + 1].x + newEfield[i][j + 1][k + 1].x) / 4;

	double ErightY = (newEfield[i][j + 1][k].y + newEfield[i + 1][j + 1][k].y + newEfield[i][j + 1][k + 1].y + newEfield[i + 1][j + 1][k + 1].y) / 4;
	double EleftY = (newEfield[i][j][k].y + newEfield[i + 1][j][k].y + newEfield[i][j][k + 1].y + newEfield[i + 1][j][k + 1].y) / 4;

	double ErightZ = (newEfield[i][j][k + 1].z + newEfield[i + 1][j][k + 1].z + newEfield[i][j + 1][k + 1].z + newEfield[i + 1][j + 1][k + 1].z) / 4;
	double EleftZ = (newEfield[i][j][k].z + newEfield[i + 1][j][k].z + newEfield[i][j + 1][k].z + newEfield[i + 1][j + 1][k].z) / 4;

	return ((ErightX - EleftX) / deltaX) + ((ErightY - EleftY) / deltaY) + ((ErightZ - EleftZ) / deltaZ);
}

Vector3d Simulation::evaluateDivPressureTensor(int i, int j, int k) {
	Vector3d result = Vector3d(0, 0, 0);

	Matrix3d tensorDerX = (getPressureTensor(i, j, k) + getPressureTensor(i, j - 1, k) + getPressureTensor(i, j, k - 1) + getPressureTensor(i, j - 1, k - 1) - getPressureTensor(i - 1, j, k) - getPressureTensor(i - 1, j - 1, k) - getPressureTensor(i - 1, j, k - 1) - getPressureTensor(i - 1, j - 1, k - 1)) / (4 * deltaX);
	Matrix3d tensorDerY = (getPressureTensor(i, j, k) + getPressureTensor(i - 1, j, k) + getPressureTensor(i, j, k - 1) + getPressureTensor(i - 1, j, k - 1) - getPressureTensor(i, j - 1, k) - getPressureTensor(i - 1, j - 1, k) - getPressureTensor(i, j - 1, k - 1) - getPressureTensor(i - 1, j - 1, k - 1)) / (4 * deltaY);
	Matrix3d tensorDerZ = (getPressureTensor(i, j, k) + getPressureTensor(i - 1, j, k) + getPressureTensor(i, j - 1, k) + getPressureTensor(i - 1, j - 1, k) - getPressureTensor(i, j, k - 1) - getPressureTensor(i - 1, j, k - 1) - getPressureTensor(i, j - 1, k - 1) - getPressureTensor(i - 1, j - 1, k - 1)) / (4 * deltaZ);

	result.x = tensorDerX.matrix[0][0] + tensorDerY.matrix[0][1] + tensorDerZ.matrix[0][2];
	result.y = tensorDerX.matrix[1][0] + tensorDerY.matrix[1][1] + tensorDerZ.matrix[1][2];
	result.z = tensorDerX.matrix[2][0] + tensorDerY.matrix[2][1] + tensorDerZ.matrix[2][2];

	return result;
}

Vector3d Simulation::evaluateGradDensity(int i, int j, int k) {
	int prevI = i - 1;
	if(prevI < 0) {
		prevI = xnumber - 1;
	}

	double densityRightX = (getDensity(i, j, k) + getDensity(i, j + 1, k) + getDensity(i, j, k + 1) + getDensity(i, j + 1, k + 1)) / 4;
	double densityLeftX = (getDensity(prevI, j, k) + getDensity(prevI, j + 1, k) + getDensity(prevI, j, k + 1) + getDensity(prevI, j + 1, k + 1)) / 4;

	double densityRightY = (getDensity(prevI, j + 1, k) + getDensity(i, j + 1, k) + getDensity(prevI, j + 1, k + 1) + getDensity(i, j + 1, k + 1)) / 4;
	double densityLeftY = (getDensity(prevI, j, k) + getDensity(i, j, k) + getDensity(prevI, j, k + 1) + getDensity(i, j, k + 1)) / 4;

	double densityRightZ = (getDensity(prevI, j, k + 1) + getDensity(i, j, k + 1) + getDensity(prevI, j + 1, k + 1) + getDensity(i, j + 1, k + 1)) / 4;
	double densityLeftZ = (getDensity(prevI, j, k) + getDensity(i, j, k) + getDensity(prevI, j + 1, k) + getDensity(i, j + 1, k)) / 4;

	double x = (densityRightX - densityLeftX) / deltaX;
	double y = (densityRightY - densityLeftY) / deltaY;
	double z = (densityRightZ - densityLeftZ) / deltaZ;

	return Vector3d(x, y, z);
}

/*Vector3d Simulation::evaluateGradPotential(int i, int j, int k) {
	int prevJ = j - 1;
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int prevK = k - 1;
	if (prevK < 0) {
		prevK = znumber - 1;
	}

	Vector3d gradient;

	double phiRightX = (divergenceCleaningPotential[i + 1][j][k][0] + divergenceCleaningPotential[i + 1][prevJ][k][0] + divergenceCleaningPotential[i + 1][j][prevK][0] + divergenceCleaningPotential[i + 1][prevJ][prevK][0]) / 4;
	double phiLeftX = (divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i][prevJ][k][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[i][prevJ][prevK][0]) / 4;
	double phiRightY = (divergenceCleaningPotential[i + 1][j][k][0] + divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i + 1][j][prevK][0] + divergenceCleaningPotential[i][j][prevK][0]) / 4;
	double phiLeftY = (divergenceCleaningPotential[i + 1][prevJ][k][0] + divergenceCleaningPotential[i][prevJ][k][0] + divergenceCleaningPotential[i + 1][prevJ][prevK][0] + divergenceCleaningPotential[i][prevJ][prevK][0]) / 4;
	double phiRightZ = (divergenceCleaningPotential[i + 1][j][k][0] + divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i + 1][prevJ][k][0] + divergenceCleaningPotential[i][prevJ][k][0]) / 4;
	double phiLeftZ = (divergenceCleaningPotential[i + 1][j][prevK][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[i + 1][prevJ][prevK][0] + divergenceCleaningPotential[i][prevJ][prevK][0]) / 4;


	if (i == 0) {
		gradient.x = 2 * (phiRightX - phiLeftX) / deltaX;
		gradient.y = 0;
		gradient.z = 0;
	} else {
		gradient.x = (phiRightX - phiLeftX) / deltaX;
		gradient.y = (phiRightY - phiLeftY) / deltaY;
		gradient.z = (phiRightZ - phiLeftZ) / deltaZ;
	}

	return gradient;
}*/