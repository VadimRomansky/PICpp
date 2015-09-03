#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "specialmath.h"


void Simulation::cleanupDivergence() {
	printf("cleaning up divergence\n");


	int matrixDimension = xnumber;



	for (int i = 0; i < xnumber; ++i) {
				divergenceCleanUpMatrix[i][0].clear();
				divergenceCleanUpMatrix[i][1].clear();
				divergenceCleanUpMatrix[i][2].clear();
				divergenceCleanUpRightPart[i][0] = 0;
				divergenceCleanUpRightPart[i][1] = 0;
				divergenceCleanUpRightPart[i][2] = 0;
				createDivergenceCleanupInternalEquation(i);
	}


	if(debugMode) {
		checkEquationMatrix(divergenceCleanUpMatrix, 1);
	}
	double rightPart = cleanUpRightPart(1);
	generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumber, 1);

	double** leftPart = multiplySpecialMatrixVector(divergenceCleanUpMatrix, divergenceCleaningPotential, xnumber, 1);
	
	updateFieldByCleaning();
	double div = evaluateDivCleaningE(1);

	updateBoundariesNewField();
}

void Simulation::updateFieldByCleaning() {
	evaluateDivergenceCleaningField();
	for(int i = 0; i < xnumber; ++i) {
				newEfield[i].x += divergenceCleaningField[i][0];
				newEfield[i].y += divergenceCleaningField[i][1];
				newEfield[i].z += divergenceCleaningField[i][2];
	}
}

void Simulation::evaluateDivergenceCleaningField(){
	for(int i = 0; i < xnumber; ++i){
				int prevI = i - 1;
				if(prevI < 0){
					prevI = xnumber - 1;
				}
				divergenceCleaningField[i][0] = -0.25*(divergenceCleaningPotential[i][0] + divergenceCleaningPotential[i][0] + divergenceCleaningPotential[i][0] + divergenceCleaningPotential[i][0] -
														     divergenceCleaningPotential[prevI][0] - divergenceCleaningPotential[prevI][0] - divergenceCleaningPotential[prevI][0] - divergenceCleaningPotential[prevI][0])/deltaX;

				divergenceCleaningField[i][1] = 0;

				divergenceCleaningField[i][2] = 0;
	}
}

void Simulation::createDivergenceCleanupInternalEquation(int i) {

	if(i == 0){
		divergenceCleanUpMatrix[i][0].push_back(MatrixElement(1, i, 0));
		divergenceCleanUpRightPart[i][0] = 0;
		return;
	}

	int prevI = i - 1;
	if(prevI < 0) {
		prevI = xnumber - 1;
	}
	int nextI = i + 1;
	if(nextI >= xnumber) {
		nextI = 0;
	}

	//div for x
	double element = -2/deltaX2;
	divergenceCleanUpMatrix[i][0].push_back(MatrixElement(element, i, 0));

	element = 1/deltaX2;
	divergenceCleanUpMatrix[i][0].push_back(MatrixElement(element, nextI, 0));
	divergenceCleanUpMatrix[i][0].push_back(MatrixElement(element, prevI, 0));

	divergenceCleanUpRightPart[i][0] = -cleanUpRightPart(i);

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

double Simulation::cleanUpRightPart(int i) {
	double div = evaluateDivNewE(i);

	return - div;
}