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

	evaluateMaxwellEquationMatrix();

	double** gmresOutput = new double*[xnumber];
//#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		gmresOutput[i] = new double[3];
	}

	generalizedMinimalResidualMethod(maxwellEquationMatrix, maxwellEquationRightPart, gmresOutput, xnumber, 3);
//#pragma omp parallel for
	for (int i = 0; i < xnumber; ++i) {
		for (int l = 0; l < 3; ++l) {
			tempEfield[i][l] = gmresOutput[i][l];
		}
		delete[] gmresOutput[i];
	}
	delete[] gmresOutput;

	updateBoundaries();

	double alfvenV = B0.norm()/sqrt(4*pi*density);
	double k = 2*pi/xsize;

	evaluateExplicitDerivative();
	//smoothEderivative();
	for(int i = 0; i < xnumber; ++i){
		implicitEfield[i] += Ederivative[i]*deltaT;
	}
	implicitEfield[xnumber] = implicitEfield[0];

	//evaluateMagneticField();

	tempEfield[xnumber] = tempEfield[0];

	for (int i = 0; i < xnumber+1; ++i) {
		newEfield[i] = (tempEfield[i] - Efield[i] * (1 - theta)) / theta;
		if(solverType == EXPLICIT){
			newEfield[i] = implicitEfield[i];
		}
	    newEfield[i].x= 0;
	}
}

void Simulation::evaluateExplicitDerivative(){
	for(int i = 0; i < xnumber + 1; ++i) {
		rotB[i] = evaluateRotB(i)*speed_of_light_normalized;
		Ederivative[i] = (evaluateRotB(i)*speed_of_light_normalized - (electricFlux[i]*4*pi/fieldScale));
	}
}

void Simulation::smoothEfield() {
	for(int i = 1; i < xnumber-1; ++i) {
		Efield[i] = (newEfield[i-1] + newEfield[i]*2 + newEfield[i+1])/4.0;
		//Efield[i].x = 0;
	}

	Efield[0] = (newEfield[0]*2 + newEfield[1] + newEfield[xnumber-1])/4.0;
	//Efield[0].x = 0;
	Efield[xnumber-1] = (newEfield[0] + newEfield[xnumber - 1]*2 + newEfield[xnumber - 2])/4.0;
	//Efield[xnumber-1].x = 0;
	Efield[xnumber] = Efield[0];
}

void Simulation::updateEfield() {
	for (int i = 0; i < xnumber + 1; ++i) {
		Efield[i] = newEfield[i];
		//Efield[i].y = newEfield[i].y;
		//Efield[i].z = newEfield[i].z;
	}
	Efield[xnumber] = Efield[0];
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

void Simulation::createPerfectConductaryBoundaryCondition() {
	int i = 0;

	maxwellEquationMatrix[i][0].push_back(MatrixElement(1.0, i, 0));
	maxwellEquationRightPart[i][0] = 0;
	maxwellEquationMatrix[i][1].push_back(MatrixElement(1.0, i, 1));
	maxwellEquationRightPart[i][1] = 0;
	maxwellEquationMatrix[i][2].push_back(MatrixElement(1.0, i, 2));
	maxwellEquationRightPart[i][2] = 0;
}

void Simulation::createInternalEquationX(int i) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);
	double element = 1 + dielectricTensor[i].matrix[0][0] + c_theta_deltaT2 * (2 + 2*dielectricTensor[i].matrix[0][0])/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 0));

	element = dielectricTensor[i].matrix[0][1] + c_theta_deltaT2*(2*dielectricTensor[i].matrix[0][1]/deltaX2);
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 1));

	element = dielectricTensor[i].matrix[0][2] + c_theta_deltaT2*(2*dielectricTensor[i].matrix[0][2]/deltaX2);
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, i, 2));

	int nextI = i + 1;
	if (nextI >= xnumber) {
		nextI = 0;
	}
	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;
	}

	element = -c_theta_deltaT2*(1  + dielectricTensor[nextI].matrix[0][0])/ deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, nextI, 0));
	element = -c_theta_deltaT2*dielectricTensor[nextI].matrix[0][1]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, nextI, 1));
	element = -c_theta_deltaT2*dielectricTensor[nextI].matrix[0][2]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, nextI, 2));

	element = -c_theta_deltaT2*(1  + dielectricTensor[prevI].matrix[0][0])/ deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 0));
	element = -c_theta_deltaT2*dielectricTensor[prevI].matrix[0][1]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 1));
	element = -c_theta_deltaT2*dielectricTensor[prevI].matrix[0][2]/deltaX2;
	maxwellEquationMatrix[i][0].push_back(MatrixElement(element, prevI, 2));

}

void Simulation::createInternalEquationY(int i) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);

	double element = 1 + dielectricTensor[i].matrix[1][1] + c_theta_deltaT2 * (2  + 2*dielectricTensor[i].matrix[1][1])/deltaX2;
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, i, 1));

	element = dielectricTensor[i].matrix[1][0];
	maxwellEquationMatrix[i][1].push_back(MatrixElement(element, i, 0));

	element = dielectricTensor[i].matrix[1][2];
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

	double element = 1 + dielectricTensor[i].matrix[2][2] + c_theta_deltaT2 * (2 + 2*dielectricTensor[i].matrix[2][2])/deltaX2;
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, i, 2));

	element = dielectricTensor[i].matrix[2][0];
	maxwellEquationMatrix[i][2].push_back(MatrixElement(element, i, 0));

	element = dielectricTensor[i].matrix[2][1];
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

	//rightPart = rightPart - (evaluateRotB(i)* speed_of_light_normalized - electricFlux[i]*4*pi) * (theta * deltaT);
	rightPart = rightPart + (evaluateRotB(i)* speed_of_light_normalized - (electricFlux[i]*4*pi/fieldScale)) * (theta * deltaT) - (evaluateGradDensity(i)*speed_of_light_normalized_sqr*theta*theta*deltaT*deltaT*4*pi/fieldScale);
	createInternalEquationX(i);
	createInternalEquationY(i);
	createInternalEquationZ(i);

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
		if (i < 0) {
			printf("i < 0\n");
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "i = %d < 0 in evaluateRotB\n", i);
			fclose(errorLogFile);
			exit(0);
		}

		if (i > xnumber) {
			printf("i > xnumber\n");
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