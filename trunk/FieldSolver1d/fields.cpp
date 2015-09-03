#include "stdlib.h"
#include "stdio.h"
#include <omp.h>

#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "matrix3d.h"
#include "specialmath.h"

void Simulation::evaluateFields() {
	printf("evaluating fields\n");

	evaluateMaxwellEquationMatrix();

	double** gmresOutput = new double*[xnumber];
#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		gmresOutput[i] = new double[3];
	}

	generalizedMinimalResidualMethod(maxwellEquationMatrix, maxwellEquationRightPart, gmresOutput, xnumber, 3);
#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		for (int l = 0; l < 3; ++l) {
			tempEfield[i][l] = gmresOutput[i][l];
		}
		delete[] gmresOutput[i];
	}
	delete[] gmresOutput;

	updateBoundaries();

	//evaluateMagneticField();

	for(int i = 0; i < xnumber; ++i){
		implicitField[i] += evaluateRotB(i)*speed_of_light_normalized*deltaT;
	}
	implicitField[xnumber] = implicitField[0];

	for (int i = 0; i < xnumber + 1; ++i) {
		newEfield[i] = (tempEfield[i] - Efield[i] * (1 - theta)) / theta;
		//newEfield[i] = implicitField[i];
	}
}

void Simulation::updateEfield() {
	for (int i = 0; i < xnumber + 1; ++i) {
		Efield[i] = newEfield[i];
	}
}

void Simulation::updateBfield() {
	for (int i = 0; i < xnumber; ++i) {
		Bfield[i] = newBfield[i];
	}
}

void Simulation::updateFields() {
	updateEfield();
	updateBfield();
}

void Simulation::evaluateMaxwellEquationMatrix() {
	for (int i = 0; i < xnumber; ++i) {
		for (int l = 0; l < 3; ++l) {
			maxwellEquationMatrix[i][l].clear();
		}
	}

#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		createInternalEquation(i);

	}

	if (debugMode) {
		checkEquationMatrix(maxwellEquationMatrix, 3);
	}
}

void Simulation::checkEquationMatrix(std::vector<MatrixElement>** matrix, int lnumber) {
#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		for (int l = 0; l < lnumber; ++l) {
			for (int m = 0; m < matrix[i][l].size(); ++m) {
				MatrixElement element = matrix[i][l][m];
				if (element.i < 0) {
					printf("element i < 0\n");
					exit(0);
				}
				if (element.i >= xnumber) {
					printf("element i >= xnumber");
					exit(0);
				}
				for (int n = m + 1; n < matrix[i][l].size(); ++n) {
					MatrixElement tempElement = matrix[i][l][n];

					if (element.equalsIndex(tempElement)) {
						printf("equals indexes\n");
						printf("current = %d %d\n", i, l);
						printf("temp = %d %d\n", element.i, element.l);
						exit(0);
					}
				}
			}
		}
	}
}

void Simulation::createPerfectConductaryBoundaryCondition() {
	int i = 0;

	maxwellEquationMatrix[i][0].push_back(MatrixElement(1.0, i, 0));
	maxwellEquationRightPart[i][0] = 0;
	maxwellEquationMatrix[i][1].push_back(MatrixElement(1.0, i, 1));
	maxwellEquationRightPart[i][1] = 0;
	maxwellEquationMatrix[i][2].push_back(MatrixElement(1.0, i, 2));
	maxwellEquationRightPart[i][2] = 0;
}

void Simulation::createInternalEquation(int i) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);
	Vector3d rightPart = Efield[i];

	//rightPart = rightPart + (evaluateRotB(i, j, k) - electricFlux[i][j][k] * 4 * pi / speed_of_light_normalized) * speed_of_light_normalized * theta * deltaT;
	rightPart = rightPart + (evaluateRotB(i)) * speed_of_light_normalized * theta * deltaT;

	double element = 1 + c_theta_deltaT2 * (2 / deltaX2);
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 0));
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, i, 1));
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, i, 2));

	int nextI = i + 1;
	if (nextI >= xnumber) {
		nextI = 0;
	}
	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;
	}

	element = -c_theta_deltaT2 / deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, nextI, 0));
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, nextI, 1));
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, nextI, 2));

	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 0));
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, prevI, 1));
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, prevI, 2));

	alertNaNOrInfinity(rightPart.x, "right part x = NaN");
	alertNaNOrInfinity(rightPart.y, "right part y = NaN");
	alertNaNOrInfinity(rightPart.z, "right part z = NaN");


	maxwellEquationRightPart[i][0] = rightPart.x;
	maxwellEquationRightPart[i][1] = rightPart.y;
	maxwellEquationRightPart[i][2] = rightPart.z;
}

void Simulation::evaluateMagneticField() {
	for (int i = 0; i < xnumber; ++i) {
		Vector3d rotEold = evaluateRotE(i);
		Vector3d rotEnew = evaluateRotNewE(i);
		newBfield[i] = Bfield[i] - (rotEold * (1 - theta) + rotEnew * theta) * speed_of_light_normalized * deltaT;
	}
}

void Simulation::updateBoundaries() {
	tempEfield[xnumber] = tempEfield[0];
}

void Simulation::updateBoundariesOldField() {
	Efield[xnumber] = Efield[0];
}

void Simulation::updateBoundariesNewField() {
	newEfield[xnumber] = newEfield[0];
}


Vector3d Simulation::evaluateRotB(int i) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i > xnumber) {
			printf("x > xnumber\n");
			exit(0);
		}
	}

	Vector3d BrightX;
	Vector3d BleftX;
	Vector3d BrightY;
	Vector3d BleftY;
	Vector3d BrightZ;
	Vector3d BleftZ;

	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;
	}
	BrightX = (getBfield(i) + getBfield(i) + getBfield(i) + getBfield(i)) / 4;
	BleftX = (getBfield(prevI) + getBfield(prevI) + getBfield(prevI) + getBfield(prevI)) / 4;

	BrightY = (getBfield(i) + getBfield(prevI) + getBfield(i) + getBfield(prevI)) / 4;
	BleftY = (getBfield(i) + getBfield(prevI) + getBfield(i) + getBfield(prevI)) / 4;

	BrightZ = (getBfield(i) + getBfield(prevI) + getBfield(i) + getBfield(prevI)) / 4;
	BleftZ = (getBfield(i) + getBfield(prevI) + getBfield(i) + getBfield(prevI)) / 4;


	double x = 0;
	double y = 0;
	double z = 0;

	x = 0;
	y = - ((BrightX.z - BleftX.z) / deltaX);
	z = (BrightX.y - BleftX.y) / deltaX;

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotTempE(int i) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			exit(0);
		}
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = (tempEfield[i + 1] + tempEfield[i + 1] + tempEfield[i + 1] + tempEfield[i + 1]) / 4;
	Vector3d EleftX = (tempEfield[i] + tempEfield[i] + tempEfield[i] + tempEfield[i]) / 4;

	Vector3d ErightY = (tempEfield[i] + tempEfield[i + 1] + tempEfield[i] + tempEfield[i + 1]) / 4;
	Vector3d EleftY = (tempEfield[i] + tempEfield[i + 1] + tempEfield[i] + tempEfield[i + 1]) / 4;

	Vector3d ErightZ = (tempEfield[i] + tempEfield[i + 1] + tempEfield[i] + tempEfield[i + 1]) / 4;
	Vector3d EleftZ = (tempEfield[i] + tempEfield[i + 1] + tempEfield[i] + tempEfield[i + 1]) / 4;

	x = 0;
	y = - ((ErightX.z - EleftX.z) / deltaX);
	z = (ErightX.y - EleftX.y) / deltaX;

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotE(int i) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			exit(0);
		}
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = (Efield[i + 1] + Efield[i + 1] + Efield[i + 1] + Efield[i + 1]) / 4;
	Vector3d EleftX = (Efield[i] + Efield[i] + Efield[i] + Efield[i]) / 4;

	Vector3d ErightY = (Efield[i] + Efield[i + 1] + Efield[i] + Efield[i + 1]) / 4;
	Vector3d EleftY = (Efield[i] + Efield[i + 1] + Efield[i] + Efield[i + 1]) / 4;

	Vector3d ErightZ = (Efield[i] + Efield[i + 1] + Efield[i] + Efield[i + 1]) / 4;
	Vector3d EleftZ = (Efield[i] + Efield[i + 1] + Efield[i] + Efield[i + 1]) / 4;

	x = 0;
	y = - ((ErightX.z - EleftX.z) / deltaX);
	z = (ErightX.y - EleftX.y) / deltaX;

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotNewE(int i) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			exit(0);
		}
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = (newEfield[i + 1] + newEfield[i + 1] + newEfield[i + 1] + newEfield[i + 1]) / 4;
	Vector3d EleftX = (newEfield[i] + newEfield[i] + newEfield[i] + newEfield[i]) / 4;

	Vector3d ErightY = (newEfield[i] + newEfield[i + 1] + newEfield[i] + newEfield[i + 1]) / 4;
	Vector3d EleftY = (newEfield[i] + newEfield[i + 1] + newEfield[i] + newEfield[i + 1]) / 4;

	Vector3d ErightZ = (newEfield[i] + newEfield[i + 1] + newEfield[i] + newEfield[i + 1]) / 4;
	Vector3d EleftZ = (newEfield[i] + newEfield[i + 1] + newEfield[i] + newEfield[i + 1]) / 4;

	x = 0;
	y = - ((ErightX.z - EleftX.z) / deltaX);
	z = (ErightX.y - EleftX.y) / deltaX;

	return Vector3d(x, y, z);
}

double Simulation::evaluateDivE(int i) {
	double ErightX = (Efield[i + 1].x + Efield[i + 1].x + Efield[i + 1].x + Efield[i + 1].x) / 4;
	double EleftX = (Efield[i].x + Efield[i].x + Efield[i].x + Efield[i].x) / 4;

	double ErightY = (Efield[i].y + Efield[i + 1].y + Efield[i].y + Efield[i + 1].y) / 4;
	double EleftY = (Efield[i].y + Efield[i + 1].y + Efield[i].y + Efield[i + 1].y) / 4;

	double ErightZ = (Efield[i].z + Efield[i + 1].z + Efield[i].z + Efield[i + 1].z) / 4;
	double EleftZ = (Efield[i].z + Efield[i + 1].z + Efield[i].z + Efield[i + 1].z) / 4;

	return (ErightX - EleftX) / deltaX;
}

double Simulation::evaluateDivCleaningE(int i) {
	double ErightX = (divergenceCleaningField[i + 1][0] + divergenceCleaningField[i + 1][0] + divergenceCleaningField[i + 1][0] + divergenceCleaningField[i + 1][0]) / 4;
	double EleftX = (divergenceCleaningField[i][0] + divergenceCleaningField[i][0] + divergenceCleaningField[i][0] + divergenceCleaningField[i][0]) / 4;

	double ErightY = (divergenceCleaningField[i][1] + divergenceCleaningField[i + 1][1] + divergenceCleaningField[i][1] + divergenceCleaningField[i + 1][1]) / 4;
	double EleftY = (divergenceCleaningField[i][1] + divergenceCleaningField[i + 1][1] + divergenceCleaningField[i][1] + divergenceCleaningField[i + 1][1]) / 4;

	double ErightZ = (divergenceCleaningField[i][2] + divergenceCleaningField[i + 1][2] + divergenceCleaningField[i][2] + divergenceCleaningField[i + 1][2]) / 4;
	double EleftZ = (divergenceCleaningField[i][2] + divergenceCleaningField[i + 1][2] + divergenceCleaningField[i][2] + divergenceCleaningField[i + 1][2]) / 4;

	return (ErightX - EleftX) / deltaX;
}

double Simulation::evaluateDivTempE(int i) {
	double ErightX = (tempEfield[i + 1].x + tempEfield[i + 1].x + tempEfield[i + 1].x + tempEfield[i + 1].x) / 4;
	double EleftX = (tempEfield[i].x + tempEfield[i].x + tempEfield[i].x + tempEfield[i].x) / 4;

	double ErightY = (tempEfield[i].y + tempEfield[i + 1].y + tempEfield[i].y + tempEfield[i + 1].y) / 4;
	double EleftY = (tempEfield[i].y + tempEfield[i + 1].y + tempEfield[i].y + tempEfield[i + 1].y) / 4;

	double ErightZ = (tempEfield[i].z + tempEfield[i + 1].z + tempEfield[i].z + tempEfield[i + 1].z) / 4;
	double EleftZ = (tempEfield[i].z + tempEfield[i + 1].z + tempEfield[i].z + tempEfield[i + 1].z) / 4;

	return (ErightX - EleftX) / deltaX;
}

double Simulation::evaluateDivNewE(int i) {
	double ErightX = (newEfield[i + 1].x + newEfield[i + 1].x + newEfield[i + 1].x + newEfield[i + 1].x) / 4;
	double EleftX = (newEfield[i].x + newEfield[i].x + newEfield[i].x + newEfield[i].x) / 4;

	double ErightY = (newEfield[i].y + newEfield[i + 1].y + newEfield[i].y + newEfield[i + 1].y) / 4;
	double EleftY = (newEfield[i].y + newEfield[i + 1].y + newEfield[i].y + newEfield[i + 1].y) / 4;

	double ErightZ = (newEfield[i].z + newEfield[i + 1].z + newEfield[i].z + newEfield[i + 1].z) / 4;
	double EleftZ = (newEfield[i].z + newEfield[i + 1].z + newEfield[i].z + newEfield[i + 1].z) / 4;

	return (ErightX - EleftX) / deltaX;
}