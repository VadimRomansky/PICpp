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

	if(solverType == IMPLICIT){

		evaluateMaxwellEquationMatrix();

		double**** gmresOutput = new double***[xnumber];
		//#pragma omp parallel for
		for (int i = 0; i < xnumber; ++i) {
			gmresOutput[i] = new double**[ynumber];
			for(int j = 0; j < ynumber; ++j){
				gmresOutput[i][j] = new double*[znumber];
				for(int k = 0; k < znumber; ++k){
					gmresOutput[i][j][k] = new double[maxwellEquationMatrixSize];
				}
			}
		}

		generalizedMinimalResidualMethod(maxwellEquationMatrix, maxwellEquationRightPart, gmresOutput, xnumber, ynumber, znumber, maxwellEquationMatrixSize);
		//#pragma omp parallel for
		for (int i = 0; i < xnumber; ++i) {
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
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

		//updateBoundaries();

		double alfvenV = B0.norm()/sqrt(4*pi*density);
		double k = 2*pi/xsize;

		evaluateExplicitDerivative();
		//smoothEderivative();
		for(int i = 0; i < xnumber; ++i){
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					explicitEfield[i][j][k] += Ederivative[i][j][k]*deltaT;
				}
			}
		}

		//evaluateMagneticField();
		if(boundaryConditionType == PERIODIC){
			tempEfield[xnumber] = tempEfield[0];
			explicitEfield[xnumber] = explicitEfield[0];
		} else {
			//tempEfield[xnumber] = E0;

			//tempEfield[xnumber] = Efield[xnumber] + evaluateRotB(xnumber - 1)*speed_of_light_normalized*deltaT*theta;

			tempEfield[xnumber] = tempEfield[xnumber - 1];

			explicitEfield[xnumber] = tempEfield[xnumber];
		}
		for (int i = 0; i < xnumber + 1; ++i) {
			for(int j = 0; j < ynumber + 1; ++j){
				for(int k = 0; k < znumber + 1; ++k){
					newEfield[i][j][k] = (tempEfield[i][j][k] - Efield[i][j][k] * (1 - theta)) / theta;
					newEfield[i][j][k].x= 0;
				}
			}
		}
	}

	if(solverType == EXPLICIT){
		evaluateExplicitDerivative();
		//smoothEderivative();
		for(int i = 0; i < xnumber; ++i){
			for(int j = 0; j < ynumber; ++j){
				for(int k = 0; k < znumber; ++k){
					explicitEfield[i][j][k] += Ederivative[i][j][k]*deltaT;
				}
			}
		}
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
			}
		}

		for(int i = 0; i < xnumber + 1; ++i){
			for(int k = 0; k < znumber; ++k){
				explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
			}
		}

		for(int i = 0; i < xnumber + 1; ++i){
			for(int j = 0; j < ynumber + 1; ++j){
				explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
			}
		}

		for (int i = 0; i < xnumber + 1; ++i) {
			for(int j = 0; j < ynumber + 1; ++j){
				for(int k = 0; k < znumber + 1; ++k){
					newEfield[i][j][k] = explicitEfield[i][j][k];
					//newEfield[i].x = 0;
				}
			}
		}
	}
}

void Simulation::evaluateExplicitDerivative(){
	int maxX = xnumber+1;
	if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
		maxX = xnumber;
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				Ederivative[xnumber][j][k] = Vector3d(0, 0, 0);
			}
		}
	}
	int minX = 0;
	if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
		minX = 1;
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				Ederivative[0][j][k] = Vector3d(0, 0, 0);
			}
		}
	}
	for(int i = minX; i < maxX; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				rotB[i][j][k] = evaluateRotB(i, j, k)*speed_of_light_normalized;
				Ederivative[i][j][k] = (evaluateRotB(i, j, k)*speed_of_light_normalized - (electricFlux[i][j][k]*4*pi/fieldScale));
				if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
					if(i == 0){
						Ederivative[i][j][k].y = 0;
						Ederivative[i][j][k].z = 0;
					}
					if(i == xnumber){
						Ederivative[i][j][k] = Vector3d(0, 0, 0);
					}
				}
			}
		}
	}
}

void Simulation::updateEfield() {
	for (int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				Efield[i][j][k] = newEfield[i][j][k];
				//Efield[i].y = newEfield[i].y;
				//Efield[i].z = newEfield[i].z;
			}
		}
	}
	if(boundaryConditionType == PERIODIC){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				Efield[xnumber][j][k] = Efield[0][j][k];
			}
		}
	}
}

void Simulation::updateBfield() {
	for (int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
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
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				for (int l = 0; l < maxwellEquationMatrixSize; ++l) {
					maxwellEquationMatrix[i][j][k][l].clear();
				}
			}
		}
	}

	#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				if((i > 0 && i < xnumber - 1) || boundaryConditionType == PERIODIC){
					createInternalEquation(i, j, k);
				} else {
					if(i == 0){
						createSuperConductorLeftEquation(j, k);
					} else {
						createFreeRightEquation(j, k);
					}
				}
			}
		}

	}

	if (debugMode) {
		checkEquationMatrix(maxwellEquationMatrix, maxwellEquationMatrixSize);
	}
}

void Simulation::checkEquationMatrix(std::vector<MatrixElement>**** matrix, int lnumber) {
//#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		for(int j = ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				for (int l = 0; l < lnumber; ++l) {
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];
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

						if (element.j < 0) {
							printf("element j < 0\n");
							errorLogFile = fopen("./output/errorLog.dat", "w");
							fprintf(errorLogFile, "element j = %d < 0\n", element.j);
							fclose(errorLogFile);
							exit(0);
						}
						if (element.j >= ynumber) {
							printf("element j >= ynumber");
							errorLogFile = fopen("./output/errorLog.dat", "w");
							fprintf(errorLogFile, "element j = %d >= ynumber = %d\n", element.j, ynumber);
							fclose(errorLogFile);
							exit(0);
						}

						if (element.k < 0) {
							printf("element k < 0\n");
							errorLogFile = fopen("./output/errorLog.dat", "w");
							fprintf(errorLogFile, "element k = %d < 0\n", element.k);
							fclose(errorLogFile);
							exit(0);
						}
						if (element.k >= znumber) {
							printf("element k >= znumber");
							errorLogFile = fopen("./output/errorLog.dat", "w");
							fprintf(errorLogFile, "element k = %d >= xnumber = %d\n", element.k, znumber);
							fclose(errorLogFile);
							exit(0);
						}
						for (int n = m + 1; n < matrix[i][j][k][l].size(); ++n) {
							MatrixElement tempElement = matrix[i][j][k][l][n];
	
							if (element.equalsIndex(tempElement)) {
								printf("equals indexes\n");
								printf("current = %d %d %d %d\n", i, l);
								printf("temp = %d %d %d %d\n", element.i, element.j, element.k, element.l);
								errorLogFile = fopen("./output/errorLog.dat", "w");
								fprintf(errorLogFile, "equal indexes current = %d %d %d %d temp = %d %d %d %d\n", i, j, k, l, element.i, element.j, element.k, element.l);
								fclose(errorLogFile);
								exit(0);
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::createSuperConductorLeftEquation(int j, int k) {
	int i = 0;

	int nextJ = j+1;
	if(nextJ >= ynumber){
		nextJ = 0;
	}
	int nextK = k+1;
	if(nextK >= znumber){
		nextK = 0;
	}
	//todo!!!
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, nextJ, k, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, nextK, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, nextJ, nextK, 0));

	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-1.0, i+1, j, k, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-1.0, i+1, nextJ, k, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-1.0, i+1, j, nextK, 0));
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-1.0, i+1, nextJ, nextK, 0));

	maxwellEquationRightPart[i][j][k][0] = -4*pi*electricDensity[0][j][k]*deltaX/fieldScale;
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(1.0, i, j, k, 1));
	maxwellEquationRightPart[i][j][k][1] = 0;
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(1.0, i, j, k, 2));
	maxwellEquationRightPart[i][j][k][2] = 0;
}

void Simulation::createFreeRightEquation(int j, int k){
	Vector3d rightPart = Vector3d(0, 0, 0);

	maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(1.0, xnumber - 1, j, k, 0));
	maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(1.0, xnumber - 1, j, k, 1));
	maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(1.0, xnumber - 1, j, k, 2));

	maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(-1.0, xnumber - 2, j, k, 0));
	maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(-1.0, xnumber - 2, j, k, 1));
	maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(-1.0, xnumber - 2, j, k, 2));
	

	alertNaNOrInfinity(rightPart.x, "right part x = NaN");
	alertNaNOrInfinity(rightPart.y, "right part y = NaN");
	alertNaNOrInfinity(rightPart.z, "right part z = NaN");

	maxwellEquationRightPart[xnumber - 1][j][k][0] = rightPart.x;
	maxwellEquationRightPart[xnumber - 1][j][k][1] = rightPart.y;
	maxwellEquationRightPart[xnumber - 1][j][k][2] = rightPart.z;
}

void Simulation::createInternalEquationX(int i, int j, int k) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);
	double element = 1.0 - dielectricTensor[i].matrix[0][0] + c_theta_deltaT2*(2.0 - 2*dielectricTensor[i].matrix[0][0])/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 0));

	element = -dielectricTensor[i].matrix[0][1] - c_theta_deltaT2*2*dielectricTensor[i].matrix[0][1]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 1));

	element = -dielectricTensor[i].matrix[0][2] - c_theta_deltaT2*2*dielectricTensor[i].matrix[0][2]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 2));

	int nextI = i + 1;
	if (nextI >= xnumber) {
		nextI = 0;
	}
	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;                   
	}

	element = -c_theta_deltaT2*(1.0  - dielectricTensor[nextI].matrix[0][0])/ deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, nextI, 0));
	element = c_theta_deltaT2*dielectricTensor[nextI].matrix[0][1]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, nextI, 1));
	element = c_theta_deltaT2*dielectricTensor[nextI].matrix[0][2]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, nextI, 2));

	element = -c_theta_deltaT2*(1.0  - dielectricTensor[prevI].matrix[0][0])/ deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 0));
	element = c_theta_deltaT2*dielectricTensor[prevI].matrix[0][1]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 1));
	element = c_theta_deltaT2*dielectricTensor[prevI].matrix[0][2]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 2));
}

void Simulation::createInternalEquationY(int i) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);

	double element = 1.0 - dielectricTensor[i].matrix[1][1] + c_theta_deltaT2*2/deltaX2;
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
	element = -c_theta_deltaT2/deltaX2;
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, nextI, 1));
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, prevI, 1));
}

void Simulation::createInternalEquationZ(int i) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);

	double element = 1.0 - dielectricTensor[i].matrix[2][2] + c_theta_deltaT2 * 2 /deltaX2;
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
	element = -c_theta_deltaT2/deltaX2;
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, nextI, 2));
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, prevI, 2));
}

void Simulation::createInternalEquation(int i) {
	Vector3d rightPart = Efield[i];
	Vector3d rightPart2 = Bfield[i];

	//rightPart = rightPart + (evaluateRotB(i)* speed_of_light_normalized - electricFlux[i]*4*pi/fieldScale) * (theta * deltaT);
	rightPart = rightPart + (evaluateRotB(i)* speed_of_light_normalized - (electricFlux[i]*4*pi/fieldScale)) * (theta * deltaT) - (evaluateGradDensity(i)*speed_of_light_normalized_sqr*theta*theta*deltaT*deltaT*4*pi/fieldScale);
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
	for (int i = 0; i < xnumber; ++i) {
		Vector3d rotEold = evaluateRotE(i);
		Vector3d rotEnew = evaluateRotNewE(i);
		Vector3d fakeRotE;
		if(i >= 1 && i < xnumber - 1) {
			fakeRotE = (evaluateRotNewE(i-1) + evaluateRotNewE(i)*2 + evaluateRotNewE(i+1))/4.0;
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

	// why /4 ?
	Matrix3d tensorDerX = (getPressureTensor(i)  - getPressureTensor(i - 1)) / (deltaX);

	result.x = tensorDerX.matrix[0][0];
	result.y = tensorDerX.matrix[0][1];
	result.z = tensorDerX.matrix[0][2];

	return result;
}

Vector3d Simulation::evaluateGradDensity(int i) {
	int prevI = i - 1;
	if(prevI < 0) {
		prevI = xnumber - 1;
	}

	double densityRightX = getDensity(i);
	double densityLeftX = getDensity(prevI);

	double x = (densityRightX - densityLeftX) / deltaX;

	return Vector3d(x, 0, 0);
}