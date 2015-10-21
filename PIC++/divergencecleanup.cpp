#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "specialmath.h"



void Simulation::cleanupDivergence() {
	printf("cleaning up divergence\n");


	int matrixDimension = xnumber;

	double fullDensity = 0;
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				fullDensity += chargeDensity[i][j][k]*volumeB(i, j, k);
			}
		}
	}
	fullDensity /= xsize;

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				chargeDensity[i][j][k] -= fullDensity;
			}
		}
	}



	for (int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				divergenceCleanUpMatrix[i][j][k][0].clear();
				divergenceCleanUpMatrix[i][j][k][1].clear();
				divergenceCleanUpMatrix[i][j][k][2].clear();
				divergenceCleanUpRightPart[i][j][k][0] = 0;
				divergenceCleanUpRightPart[i][j][k][1] = 0;
				divergenceCleanUpRightPart[i][j][k][2] = 0;
				createDivergenceCleanupInternalEquation(i, j, k);
			}
		}
	}


	if(debugMode) {
		checkEquationMatrix(divergenceCleanUpMatrix, 3);
	}

	/*double summRightPart = 0;
	for(int i = 0; i < xnumber; ++i){
		summRightPart += divergenceCleanUpRightPart[i][0];
	}*/
	//double rightPart = cleanUpRightPart(1);
	//generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumber, 1);
	generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningField, xnumber, ynumber, znumber, 3);


	//double** leftPart = multiplySpecialMatrixVector(divergenceCleanUpMatrix, divergenceCleaningPotential, xnumber, 1);
	
	updateFieldByCleaning();
	//double div = evaluateDivCleaningE(1);

	updateBoundariesNewField();
}

void Simulation::updateFieldByCleaning() {
	//evaluateDivergenceCleaningField();
	/*Vector3d meanField = Vector3d(0, 0, 0);
	for(int i = 0; i < xnumber; ++i){
		meanField.x += divergenceCleaningField[i][0];
		meanField.y += divergenceCleaningField[i][1];
		meanField.z += divergenceCleaningField[i][2];
	}
	meanField.x /= xnumber;
	meanField.y /= xnumber;
	meanField.z /= xnumber;*/
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				newEfield[i][j][k].x += divergenceCleaningField[i][j][k][0];
				newEfield[i][j][k].y += divergenceCleaningField[i][j][k][1];
				newEfield[i][j][k].z += divergenceCleaningField[i][j][k][2];
				//newEfield[i] -= meanField;
			}
		}
	}
	if(boundaryConditionType == PERIODIC){
		newEfield[xnumber] = newEfield[0];
	}
}

void Simulation::evaluateDivergenceCleaningField(){
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
					int prevI = i - 1;
					if(prevI < 0){
						prevI = xnumber - 1;
					}

					int prevJ = j - 1;
					if(prevJ < 0){
						prevJ = ynumber - 1;
					}

					int prevK = k - 1;
					if(prevK < 0){
						prevK = znumber - 1;
					}

					divergenceCleaningField[i][j][k][0] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i][prevJ][k][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[i][prevJ][prevK][0]
														  - divergenceCleaningPotential[prevI][j][k][0] - divergenceCleaningPotential[prevI][prevJ][k][0] - divergenceCleaningPotential[prevI][j][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0])/(4*deltaX);

					divergenceCleaningField[i][j][k][1] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[prevI][j][k][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[prevI][j][prevK][0]
														  - divergenceCleaningPotential[i][prevJ][k][0] - divergenceCleaningPotential[prevI][prevJ][k][0] - divergenceCleaningPotential[i][prevJ][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0])/(4*deltaY);

					divergenceCleaningField[i][j][k][2] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i][prevJ][k][0] + divergenceCleaningPotential[prevI][j][k][0] + divergenceCleaningPotential[prevI][prevJ][k][0]
														  - divergenceCleaningPotential[i][j][prevK][0] - divergenceCleaningPotential[i][prevJ][prevK][0] - divergenceCleaningPotential[prevI][j][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0])/(4*deltaZ);
			}
		}
	}
}

void Simulation::createDivergenceCleanupInternalEquation(int i, int j, int k) {

	/*if(i == 0){
		divergenceCleanUpMatrix[i][0].push_back(MatrixElement(1, i, 0));
		divergenceCleanUpMatrix[i][0].push_back(MatrixElement(-1, i+1, 0));
		divergenceCleanUpRightPart[i][0] = 0;
		return;
	}*/
	
	int prevI = i - 1;
	if(prevI < 0) {
		prevI = xnumber - 1;
	}
	int nextI = i + 1;
	if(nextI >= xnumber) {
		nextI = 0;
	}

	//div for x
	/*double element = -2/deltaX2;
	divergenceCleanUpMatrix[i][0].push_back(MatrixElement(element, i, 0));

	element = 1/deltaX2;
	divergenceCleanUpMatrix[i][0].push_back(MatrixElement(element, nextI, 0));
	divergenceCleanUpMatrix[i][0].push_back(MatrixElement(element, prevI, 0));

	divergenceCleanUpRightPart[i][0] = -cleanUpRightPart(i);*/

	//divergenceCleanUpMatrix[i][j][k][1].push_back(MatrixElement(1, i,  1));
	//divergenceCleanUpMatrix[i][j][k][2].push_back(MatrixElement(1, i, 2));
	divergenceCleanUpRightPart[i][j][k][1] = 0;
	divergenceCleanUpRightPart[i][j][k][2] = 0;

	divergenceCleanUpRightPart[i][j][k][0] = cleanUpRightPart(i, j, k);

	//divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(-1.0/deltaX, i, 0));
	//divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1.0/deltaX, nextI, 0));

}

void Simulation::createDivergenceCleanupLeftEquation() {
	int i = 0;

	int nextI = i+1;

	//div for x

	divergenceCleanUpRightPart[i][0] = 0;
}

void Simulation::createDivergenceCleanupRightEquation() {
	int i = xnumber - 1;

	//div for x
	

	divergenceCleanUpRightPart[i][2] = 0;
}

double Simulation::cleanUpRightPart(int i, int j, int k) {
	double div = evaluateDivNewE(i, j, k);

	return  4*pi*chargeDensity[i][j][k] - div;
}