#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>

#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "matrix3d.h"
#include "specialmath.h"

void Simulation::evaluateFields() {
	printf("evaluating fields\n");

	if (solverType == IMPLICIT) {

		evaluateMaxwellEquationMatrix();

		double** gmresOutput = new double*[xnumber];
		//#pragma omp parallel for
		for (int i = 0; i < xnumber; ++i) {
			gmresOutput[i] = new double[maxwellEquationMatrixSize];
		}

		generalizedMinimalResidualMethod(maxwellEquationMatrix, maxwellEquationRightPart, gmresOutput, xnumber, maxwellEquationMatrixSize);
		printf("end of GMRES\n");
		//#pragma omp parallel for
		for (int i = 0; i < xnumber; ++i) {
			for (int l = 0; l < 3; ++l) {
				tempEfield[i][l] = gmresOutput[i][l];
			}
			/*for (int l = 3; l <maxwellEquationMatrixSize; ++l){
				newBfield[i][l-3] = gmresOutput[i][l];
			}*/
			delete[] gmresOutput[i];
		}
		delete[] gmresOutput;

		//updateBoundaries();

		double alfvenV = B0.norm() / sqrt(4 * pi * density);
		double k = 2 * pi / xsize;

		evaluateExplicitDerivative();
		//smoothEderivative();
		for (int i = 0; i < xnumber; ++i) {
			explicitEfield[i] += Ederivative[i] * deltaT;
		}

		//evaluateMagneticField();
		if (boundaryConditionType == PERIODIC) {
			tempEfield[xnumber] = tempEfield[0];
			explicitEfield[xnumber] = explicitEfield[0];
		} else {
			tempEfield[xnumber] = currentRightField;

			explicitEfield[xnumber] = currentRightField;

			//tempEfield[xnumber] = Efield[xnumber] + evaluateRotB(xnumber - 1)*speed_of_light_normalized*deltaT*theta;

			//tempEfield[xnumber] = tempEfield[xnumber - 1];

			//explicitEfield[xnumber] = tempEfield[xnumber];
		}
		//smoothEfield();
		for (int i = 0; i < xnumber + 1; ++i) {
			newEfield[i] = (tempEfield[i] - Efield[i] * (1 - theta)) / theta;
			//newEfield[i].x= 0;
		}
	}

	if (solverType == EXPLICIT) {
		evaluateExplicitDerivative();
		//smoothEderivative();
		for (int i = 0; i < xnumber; ++i) {
			explicitEfield[i] += Ederivative[i] * deltaT;
		}
		if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
			explicitEfield[xnumber] = explicitEfield[xnumber - 1];
			explicitEfield[0].y = 0;
			explicitEfield[0].z = 0;
			explicitEfield[0].x = explicitEfield[1].x - 4 * pi * electricDensity[0] * deltaX;
		}
		if (boundaryConditionType == PERIODIC) {
			explicitEfield[xnumber] = explicitEfield[0];
		}
		for (int i = 0; i < xnumber + 1; ++i) {
			newEfield[i] = explicitEfield[i];
			tempEfield[i] = newEfield[i] * theta + Efield[i] * (1 - theta);
			//newEfield[i].x = 0;
		}
	}
}

void Simulation::evaluateExplicitDerivative() {
	int maxX = xnumber + 1;
	if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
		maxX = xnumber;
		Ederivative[xnumber] = Vector3d(0, 0, 0);
	}
	int minX = 0;
	if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
		minX = 1;
		Ederivative[0] = Vector3d(0, 0, 0);
	}
	for (int i = minX; i < maxX; ++i) {
		rotB[i] = evaluateRotB(i) * speed_of_light_normalized;
		Ederivative[i] = (evaluateRotB(i) * speed_of_light_normalized - (electricFlux[i] * 4 * pi / fieldScale));
		if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
			if (i == 0) {
				Ederivative[i].y = 0;
				Ederivative[i].z = 0;
			}
			if (i == xnumber) {
				Ederivative[i] = Vector3d(0, 0, 0);
			}
		}
	}
}

/*void Simulation::smoothEfield() {
	for(int i = 1; i < xnumber-1; ++i) {
		Efield[i] = (newEfield[i-1] + newEfield[i]*2 + newEfield[i+1])/4.0;
		//Efield[i].x = 0;
	}

	Efield[0] = newEfield[0];
	Efield[xnumber] = newEfield[xnumber];
	//Efield[0] = (newEfield[0]*2 + newEfield[1] + newEfield[xnumber-1])/4.0;
	//Efield[0].x = 0;
	//Efield[xnumber-1] = (newEfield[0] + newEfield[xnumber - 1]*2 + newEfield[xnumber - 2])/4.0;
	//Efield[xnumber-1].x = 0;
	//Efield[xnumber] = Efield[0];

	for(int i = 0; i < xnumber + 1; ++i){
		newEfield[i] = Efield[i];
	}
}*/

void Simulation::smoothEfield() {
	double x = 0.001;
	for (int i = 1; i < xnumber; ++i) {
		newEfield[i] = tempEfield[i] * (1 - x) + ((tempEfield[i - 1] + tempEfield[i] * 2.0 + tempEfield[i + 1]) * x / 4.0);
		//Efield[i].x = 0;
	}

	newEfield[0] = tempEfield[0];
	newEfield[xnumber] = tempEfield[xnumber];
	//Efield[0] = (newEfield[0]*2 + newEfield[1] + newEfield[xnumber-1])/4.0;
	//Efield[0].x = 0;
	//Efield[xnumber-1] = (newEfield[0] + newEfield[xnumber - 1]*2 + newEfield[xnumber - 2])/4.0;
	//Efield[xnumber-1].x = 0;
	//Efield[xnumber] = Efield[0];

	for (int i = 0; i < xnumber + 1; ++i) {
		tempEfield[i] = newEfield[i];
	}
}

void Simulation::updateEfield() {
	for (int i = 0; i < xnumber + 1; ++i) {
		Efield[i] = newEfield[i];
		//Efield[i].y = newEfield[i].y;
		//Efield[i].z = newEfield[i].z;
	}
	if (boundaryConditionType == PERIODIC) {
		Efield[xnumber] = Efield[0];
	}
	//smoothEfield();
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
		for (int l = 0; l < maxwellEquationMatrixSize; ++l) {
			maxwellEquationMatrix[i][l].clear();
		}
	}

#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		if ((i > 0 && i < xnumber - 1) || boundaryConditionType == PERIODIC) {
			createInternalEquation(i);
		} else {
			if (i == 0) {
				if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
					createSuperConductorLeftEquation();
				} else {
					createFreeLeftEquation();
				}
			} else {
				createFreeRightEquation();
			}
		}

	}

	if (debugMode) {
		checkEquationMatrix(maxwellEquationMatrix, maxwellEquationMatrixSize);
	}
}

void Simulation::checkEquationMatrix(std::vector<MatrixElement>** matrix, int lnumber) {
	//#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		for (int l = 0; l < lnumber; ++l) {
			for (int m = 0; m < matrix[i][l].size(); ++m) {
				MatrixElement element = matrix[i][l][m];
				if (element.i < 0) {
					printf("element i < 0\n");
					errorLogFile = fopen("./output/errorLog.dat", "w");
					fprintf(errorLogFile, "element i = %d < 0\n", element.i);
					fclose(errorLogFile);
					exit(0);
				}
				if (element.i >= xnumber) {
					printf("element i >= xnumber");
					errorLogFile = fopen("./output/errorLog.dat", "w");
					fprintf(errorLogFile, "element i = %d >= xnumber = %d\n", element.i, xnumber);
					fclose(errorLogFile);
					exit(0);
				}
				for (int n = m + 1; n < matrix[i][l].size(); ++n) {
					MatrixElement tempElement = matrix[i][l][n];

					if (element.equalsIndex(tempElement)) {
						printf("equals indexes\n");
						printf("current = %d %d\n", i, l);
						printf("temp = %d %d\n", element.i, element.l);
						errorLogFile = fopen("./output/errorLog.dat", "w");
						fprintf(errorLogFile, "equal indexes current = %d %d temp = %d %d\n", i, l, element.i, element.l);
						fclose(errorLogFile);
						exit(0);
					}
				}
			}
		}
	}
}

void Simulation::createSuperConductorLeftEquation() {
	int i = 0;
	//todo!!!
	maxwellEquationMatrix[i][0].push_back(MatrixElement(1.0, i, 0));
	maxwellEquationMatrix[i][0].push_back(MatrixElement(-1.0, i + 1, 0));
	maxwellEquationRightPart[i][0] = -4 * pi * electricDensity[0] * deltaX / fieldScale;
	maxwellEquationMatrix[i][1].push_back(MatrixElement(1.0, i, 1));
	maxwellEquationRightPart[i][1] = 0;
	maxwellEquationMatrix[i][2].push_back(MatrixElement(1.0, i, 2));
	maxwellEquationRightPart[i][2] = 0;
}

void Simulation::createFreeRightEquation() {
	//Vector3d rightPart = Vector3d(0, 0, 0);
	Vector3d rightPart = currentRightField;

	maxwellEquationMatrix[xnumber - 1][0].push_back(MatrixElement(1.0, xnumber - 1, 0));
	maxwellEquationMatrix[xnumber - 1][1].push_back(MatrixElement(1.0, xnumber - 1, 1));
	maxwellEquationMatrix[xnumber - 1][2].push_back(MatrixElement(1.0, xnumber - 1, 2));

	//maxwellEquationMatrix[xnumber - 1][0].push_back(MatrixElement(-1.0, xnumber - 2, 0));
	//maxwellEquationMatrix[xnumber - 1][1].push_back(MatrixElement(-1.0, xnumber - 2, 1));
	//maxwellEquationMatrix[xnumber - 1][2].push_back(MatrixElement(-1.0, xnumber - 2, 2));


	alertNaNOrInfinity(rightPart.x, "right part x = NaN");
	alertNaNOrInfinity(rightPart.y, "right part y = NaN");
	alertNaNOrInfinity(rightPart.z, "right part z = NaN");

	maxwellEquationRightPart[xnumber - 1][0] = rightPart.x;
	maxwellEquationRightPart[xnumber - 1][1] = rightPart.y;
	maxwellEquationRightPart[xnumber - 1][2] = rightPart.z;
}

void Simulation::createFreeLeftEquation() {
	Vector3d rightPart = Vector3d(0, 0, 0);

	maxwellEquationMatrix[0][0].push_back(MatrixElement(1.0, 0, 0));
	maxwellEquationMatrix[0][1].push_back(MatrixElement(1.0, 0, 1));
	maxwellEquationMatrix[0][2].push_back(MatrixElement(1.0, 0, 2));

	maxwellEquationMatrix[0][0].push_back(MatrixElement(-1.0, 1, 0));
	maxwellEquationMatrix[0][1].push_back(MatrixElement(-1.0, 1, 1));
	maxwellEquationMatrix[0][2].push_back(MatrixElement(-1.0, 1, 2));


	alertNaNOrInfinity(rightPart.x, "right part x = NaN");
	alertNaNOrInfinity(rightPart.y, "right part y = NaN");
	alertNaNOrInfinity(rightPart.z, "right part z = NaN");

	maxwellEquationRightPart[0][0] = rightPart.x;
	maxwellEquationRightPart[0][1] = rightPart.y;
	maxwellEquationRightPart[0][2] = rightPart.z;
}

/*void Simulation::createFreeRightEquation(){
	Vector3d rightPart = Efield[xnumber - 1];
	//rightPart = rightPart + (evaluateRotB(xnumber - 1)* speed_of_light_normalized - (electricFlux[xnumber - 1]*4*pi/fieldScale)) * (theta * deltaT) - (evaluateGradDensity(xnumber - 1)*speed_of_light_normalized_sqr*theta*theta*deltaT*deltaT*4*pi/fieldScale);
	rightPart = rightPart + (evaluateRotB(xnumber - 1)* speed_of_light_normalized - (electricFlux[xnumber - 1]*4*pi/fieldScale)) * (theta * deltaT);

	createFreeRightEquationX(rightPart);
	createFreeRightEquationY(rightPart);
	createFreeRightEquationZ(rightPart);

	alertNaNOrInfinity(rightPart.x, "right part x = NaN");
	alertNaNOrInfinity(rightPart.y, "right part y = NaN");
	alertNaNOrInfinity(rightPart.z, "right part z = NaN");

	maxwellEquationRightPart[xnumber - 1][0] = rightPart.x;
	maxwellEquationRightPart[xnumber - 1][1] = rightPart.y;
	maxwellEquationRightPart[xnumber - 1][2] = rightPart.z;
}

void Simulation::createFreeRightEquationX(Vector3d& rightPart) {
	int i = xnumber - 1;

	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);
	double element = 1.0 - dielectricTensor[i].matrix[0][0] + c_theta_deltaT2*(2.0 - 2*dielectricTensor[i].matrix[0][0])/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 0));

	element = -dielectricTensor[i].matrix[0][1] - c_theta_deltaT2*2*dielectricTensor[i].matrix[0][1]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 1));

	element = -dielectricTensor[i].matrix[0][2] - c_theta_deltaT2*2*dielectricTensor[i].matrix[0][2]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 2));

	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;                   
	}

	element = -c_theta_deltaT2*(1.0  - dielectricTensor[xnumber-1].matrix[0][0])/ deltaX2;
	rightPart.x -= element*Efield[xnumber-1].x;
	element = c_theta_deltaT2*dielectricTensor[xnumber-1].matrix[0][1]/deltaX2;
	rightPart.x -= element*Efield[xnumber-1].y;
	element = c_theta_deltaT2*dielectricTensor[xnumber-1].matrix[0][2]/deltaX2;
	rightPart.x -= element*Efield[xnumber-1].z;

	element = -c_theta_deltaT2*(1.0  - dielectricTensor[prevI].matrix[0][0])/ deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 0));
	element = c_theta_deltaT2*dielectricTensor[prevI].matrix[0][1]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 1));
	element = c_theta_deltaT2*dielectricTensor[prevI].matrix[0][2]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 2));
}

void Simulation::createFreeRightEquationY(Vector3d& rightPart) {
	int i = xnumber - 1;

	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);

	double element = 1.0 - dielectricTensor[i].matrix[1][1] + c_theta_deltaT2*2/deltaX2;
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, i, 1));

	element = -dielectricTensor[i].matrix[1][0];
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, i, 0));
	
	element = -dielectricTensor[i].matrix[1][2];
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, i, 2));

	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;
	}
	element = -c_theta_deltaT2/deltaX2;
	rightPart.y -= element*Efield[xnumber-1].y;
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, prevI, 1));
}

void Simulation::createFreeRightEquationZ(Vector3d& rightPart) {
	int i = xnumber - 1;
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);

	double element = 1.0 - dielectricTensor[i].matrix[2][2] + c_theta_deltaT2 * 2 /deltaX2;
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, i, 2));

	element = -dielectricTensor[i].matrix[2][0];
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, i, 0));

	element = -dielectricTensor[i].matrix[2][1];
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, i, 1));

	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;
	}
	element = -c_theta_deltaT2/deltaX2;
	rightPart.z -= element*Efield[xnumber-1].z;
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, prevI, 2));
}*/

void Simulation::createInternalEquationX(int i) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);
	double element = 1.0 - dielectricTensor[i].matrix[0][0] + c_theta_deltaT2 * (2.0 - 2 * dielectricTensor[i].matrix[0][0]) / deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 0));

	element = -dielectricTensor[i].matrix[0][1] - c_theta_deltaT2 * 2 * dielectricTensor[i].matrix[0][1] / deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 1));

	element = -dielectricTensor[i].matrix[0][2] - c_theta_deltaT2 * 2 * dielectricTensor[i].matrix[0][2] / deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 2));

	int nextI = i + 1;
	if (nextI >= xnumber) {
		nextI = 0;
	}
	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;
	}

	element = c_theta_deltaT2 * (-1.0 + dielectricTensor[nextI].matrix[0][0]) / deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, nextI, 0));
	element = c_theta_deltaT2 * dielectricTensor[nextI].matrix[0][1] / deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, nextI, 1));
	element = c_theta_deltaT2 * dielectricTensor[nextI].matrix[0][2] / deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, nextI, 2));

	element = c_theta_deltaT2 * (-1.0 + dielectricTensor[prevI].matrix[0][0]) / deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 0));
	element = c_theta_deltaT2 * dielectricTensor[prevI].matrix[0][1] / deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 1));
	element = c_theta_deltaT2 * dielectricTensor[prevI].matrix[0][2] / deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 2));
}

void Simulation::createInternalEquationY(int i) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);

	double element = 1.0 - dielectricTensor[i].matrix[1][1] + c_theta_deltaT2 * 2.0 / deltaX2;
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, i, 1));

	element = -dielectricTensor[i].matrix[1][0];
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, i, 0));

	element = -dielectricTensor[i].matrix[1][2];
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, i, 2));

	int nextI = i + 1;
	if (nextI >= xnumber) {
		nextI = 0;
	}
	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;
	}
	element = -c_theta_deltaT2 / deltaX2;
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, nextI, 1));
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, prevI, 1));
}

void Simulation::createInternalEquationZ(int i) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);

	double element = 1.0 - dielectricTensor[i].matrix[2][2] + c_theta_deltaT2 * 2.0 / deltaX2;
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, i, 2));

	element = -dielectricTensor[i].matrix[2][0];
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, i, 0));

	element = -dielectricTensor[i].matrix[2][1];
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, i, 1));

	int nextI = i + 1;
	if (nextI >= xnumber) {
		nextI = 0;
	}
	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;
	}
	element = -c_theta_deltaT2 / deltaX2;
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, nextI, 2));
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, prevI, 2));
}

void Simulation::createInternalEquation(int i) {
	Vector3d rightPart = Efield[i];
	Vector3d rightPart2 = Bfield[i];

	//rightPart = rightPart + (evaluateRotB(i)* speed_of_light_normalized - electricFlux[i]*4*pi/fieldScale) * (theta * deltaT);
	rightPart = rightPart + (evaluateRotB(i) * speed_of_light_normalized - (electricFlux[i] * 4 * pi / fieldScale)) * (theta * deltaT) - (evaluateGradDensity(i) * speed_of_light_normalized_sqr * theta * theta * deltaT * deltaT * 4 * pi / fieldScale);
	//rightPart = rightPart + evaluateRotB(i)*speed_of_light_normalized*theta*deltaT - electricFlux[i]*4*pi*theta*deltaT/fieldScale;
	createInternalEquationX(i);
	createInternalEquationY(i);
	createInternalEquationZ(i);

	alertNaNOrInfinity(rightPart.x, "right part x = NaN");
	alertNaNOrInfinity(rightPart.y, "right part y = NaN");
	alertNaNOrInfinity(rightPart.z, "right part z = NaN");

	alertNaNOrInfinity(rightPart2.x, "right part 2 x = NaN");
	alertNaNOrInfinity(rightPart2.y, "right part 2 y = NaN");
	alertNaNOrInfinity(rightPart2.z, "right part 2 z = NaN");


	maxwellEquationRightPart[i][0] = rightPart.x;
	maxwellEquationRightPart[i][1] = rightPart.y;
	maxwellEquationRightPart[i][2] = rightPart.z;

	//maxwellEquationRightPart[i][3] = rightPart2.x;
	//maxwellEquationRightPart[i][4] = rightPart2.y;
	//maxwellEquationRightPart[i][5] = rightPart2.z;
}

void Simulation::evaluateMagneticField() {
	printf("updating magnetic field\n");
	for (int i = 0; i < xnumber; ++i) {
		Vector3d rotEold = evaluateRotE(i);
		Vector3d rotEnew = evaluateRotNewE(i);
		Vector3d fakeRotE;
		if (i >= 1 && i < xnumber - 1) {
			fakeRotE = (evaluateRotNewE(i - 1) + evaluateRotNewE(i) * 2 + evaluateRotNewE(i + 1)) / 4.0;
		}
		newBfield[i] = Bfield[i] - (rotEold * (1 - theta) + rotEnew * theta) * (speed_of_light_normalized * deltaT);
		//newBfield[i] = Bfield[i] - rotEold*speed_of_light_normalized*deltaT;
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
		if ((i < 0) || ((i == 0) && (boundaryConditionType == SUPER_CONDUCTOR_LEFT))) {
			printf("i < 0\n");
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "i = %d < 0 in evaluateRotB\n", i);
			fclose(errorLogFile);
			exit(0);
		}

		if ((i > xnumber) || ((i == xnumber) && (boundaryConditionType == SUPER_CONDUCTOR_LEFT))) {
			printf("i >= xnumber\n");
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "i = %d > xnumber = %d in evaluzteRotB\n", i, xnumber);
			fclose(errorLogFile);
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
	BrightX = getBfield(i);
	BleftX = getBfield(prevI);

	BrightY = (getBfield(prevI) + getBfield(i)) / 2;
	BleftY = (getBfield(i) + getBfield(prevI)) / 2;

	BrightZ = (getBfield(i) + getBfield(prevI)) / 2;
	BleftZ = (getBfield(i) + getBfield(prevI)) / 2;


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
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "i = %d < 0 in evaluateRotTempE\n", i);
			fclose(errorLogFile);
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "i = %d >= xnumber = %d in evaluateRotTempE\n", i, xnumber);
			fclose(errorLogFile);
			exit(0);
		}
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = tempEfield[i + 1];
	Vector3d EleftX = tempEfield[i];

	x = 0;
	y = - ((ErightX.z - EleftX.z) / deltaX);
	z = (ErightX.y - EleftX.y) / deltaX;

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotE(int i) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "i = %d < 0 in evaluateRotE\n", i);
			fclose(errorLogFile);
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "i = %d >= xnumber = %d in evaluateRotE\n", i, xnumber);
			fclose(errorLogFile);
			exit(0);
		}
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = Efield[i + 1];
	Vector3d EleftX = Efield[i];

	x = 0;
	y = - ((ErightX.z - EleftX.z) / deltaX);
	z = (ErightX.y - EleftX.y) / deltaX;

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotNewE(int i) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "i = %d < 0 in evaluateRotNewE\n", i);
			fclose(errorLogFile);
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "i = %d >= xnumber = %d in evaluateRotNewE\n", i, xnumber);
			fclose(errorLogFile);
			exit(0);
		}
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = newEfield[i + 1];
	Vector3d EleftX = newEfield[i];


	x = 0;
	y = - ((ErightX.z - EleftX.z) / deltaX);
	z = (ErightX.y - EleftX.y) / deltaX;

	return Vector3d(x, y, z);
}

double Simulation::evaluateDivE(int i) {
	double ErightX = Efield[i + 1].x;
	double EleftX = Efield[i].x;


	return (ErightX - EleftX) / deltaX;
}

double Simulation::evaluateDivCleaningE(int i) {
	double ErightX = divergenceCleaningField[i + 1][0];
	double EleftX = divergenceCleaningField[i][0];

	return (ErightX - EleftX) / deltaX;
}

double Simulation::evaluateDivTempE(int i) {
	double ErightX = tempEfield[i + 1].x;
	double EleftX = tempEfield[i].x;

	return (ErightX - EleftX) / deltaX;
}

double Simulation::evaluateDivNewE(int i) {
	double ErightX = newEfield[i + 1].x;
	double EleftX = newEfield[i].x;

	return (ErightX - EleftX) / deltaX;
}

double Simulation::evaluateDivFlux(int i) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "i = %d < 0 in evaluateDivFlux\n", i);
			fclose(errorLogFile);
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "i = %d >= xnumber = %d in evaluateDivFlux\n", i, xnumber);
			fclose(errorLogFile);
			exit(0);
		}
	}


	double rfluxX = electricFlux[i + 1].x;
	double lfluxX = electricFlux[i].x;

	return (rfluxX - lfluxX) / deltaX;
}

Vector3d Simulation::evaluateDivPressureTensor(int i) {
	Vector3d result = Vector3d(0, 0, 0);

	Matrix3d tensorDerX = (getPressureTensor(i) - getPressureTensor(i - 1)) / (deltaX);

	result.x = tensorDerX.matrix[0][0];
	result.y = tensorDerX.matrix[0][1];
	result.z = tensorDerX.matrix[0][2];

	return result;
}

Vector3d Simulation::evaluateGradDensity(int i) {
	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;
	}

	double densityRightX = getDensity(i);
	double densityLeftX = getDensity(prevI);

	double x = (densityRightX - densityLeftX) / deltaX;

	return Vector3d(x, 0, 0);
}

void Simulation::fourierFilter() {
	printf("fourier filtering\n");
	fourierFilter(0, shockWavePoint);
	//fourierFilter(shockWavePoint, xnumber);
}


void Simulation::fourierFilter(int startPoint, int endPoint) {

	double length = xgrid[endPoint] - xgrid[startPoint];

	double minWaveLength = 2 * deltaX;
	double maxWaveLength = length / 4;


	int maxHarmonicNumber = length / minWaveLength;
	int minHarmonicNumber = length / maxWaveLength;

	for (int harmCounter = minHarmonicNumber; harmCounter <= maxHarmonicNumber; ++harmCounter) {
		Vector3d Ecos = Vector3d(0, 0, 0);
		Vector3d Esin = Vector3d(0, 0, 0);
		Vector3d Bcos = Vector3d(0, 0, 0);
		Vector3d Bsin = Vector3d(0, 0, 0);

		double k = 2 * pi * harmCounter / length;

		//for(int i = 0; i < xnumber; ++i
		for (int i = startPoint; i < endPoint; ++i) {
			Bcos = Bcos + Bfield[i] * cos(k * middleXgrid[i]) * deltaX * 2 / length;
			Bsin = Bsin + Bfield[i] * sin(k * middleXgrid[i]) * deltaX * 2 / length;
		}

		Ecos = Ecos + Efield[startPoint] * cos(k * xgrid[startPoint]) * (deltaX / 2) * 2 / length;
		Esin = Esin + Efield[startPoint] * sin(k * xgrid[startPoint]) * (deltaX / 2) * 2 / length;
		for (int i = startPoint + 1; i < endPoint; ++i) {
			Ecos = Ecos + Efield[i] * cos(k * xgrid[i]) * deltaX * 2 / length;
			Esin = Esin + Efield[i] * sin(k * xgrid[i]) * deltaX * 2 / length;
		}
		Ecos = Ecos + Efield[endPoint] * cos(k * xgrid[endPoint]) * (deltaX / 2) * 2 / length;
		Esin = Esin + Efield[endPoint] * sin(k * xgrid[endPoint]) * (deltaX / 2) * 2 / length;

		Bcos.x = 0;
		Bsin.x = 0;

		for (int i = startPoint; i < endPoint; ++i) {
			//for(int i = 0; i < xnumber; ++i){
			Bfield[i] = Bfield[i] - Bcos * cos(k * middleXgrid[i]);
			Bfield[i] = Bfield[i] - Bsin * sin(k * middleXgrid[i]);
		}

		for (int i = startPoint; i < endPoint + 1; ++i) {
			//for(int i = 0; i < xnumber+1; ++i){
			Efield[i] = Efield[i] - Ecos * cos(k * xgrid[i]);
			Efield[i] = Efield[i] - Esin * sin(k * xgrid[i]);
		}
	}

}

void Simulation::updateRightFields() {
	int n = 20;
	double halfLength = xsize / (2 * n);
	double halfPeriod = halfLength / fabs(V0.x);
	int count = time / halfPeriod;
	if ((count % 2) == 0) {
		currentRightField = E0;
	} else {
		currentRightField = E0 * (-1.0);
	}
}
