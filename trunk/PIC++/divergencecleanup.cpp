#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "specialmath.h"


void Simulation::cleanupDivergence() {
	printf("cleaning up divergence\n");

	double fullDensity = 0;
	double fullDiv = 0;
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fullDensity += chargeDensity[i][j][k]*volume(i, j, k);
				fullDiv += evaluateDivNewE(i, j, k);
			}
		}
	}
	fullDensity /= (xsize*ysize*zsize);

	int matrixDimension = xnumber*ynumber*znumber;
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				//chargeDensity[i][j][k] -= (fullDensity);
			}
		}
	}
	
	fullDensity = 0;
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fullDensity += chargeDensity[i][j][k]*volume(i, j, k);
				fullDiv += evaluateDivNewE(i, j, k);
			}
		}
	}
	fullDensity /= (xsize*ysize*zsize);

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				divergenceCleanUpMatrix[i][j][k][0].clear();
				divergenceCleanUpMatrix[i][j][k][1].clear();
				divergenceCleanUpMatrix[i][j][k][2].clear();
				divergenceCleanUpRightPart[i][j][k][0] = 0;
				divergenceCleanUpRightPart[i][j][k][1] = 0;
				divergenceCleanUpRightPart[i][j][k][2] = 0;
				if(i == 0 && boundaryConditionType == SUPERCONDUCTERLEFT){
					createDivergenceCleanupLeftEquation(j, k);
				} else if((i == xnumber - 1) && (boundaryConditionType == SUPERCONDUCTERLEFT)) {
					createDivergenceCleanupRightEquation(j, k);
				} else {
					createDivergenceCleanupInternalEquation(i, j, k);
				}
			}
		}
	}


	if(debugMode) {
		checkEquationMatrix(divergenceCleanUpMatrix, 1);
	}
	double rightPart = cleanUpRightPart(1,1,1);
	generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumber, ynumber, znumber, 1);

	double**** leftPart = multiplySpecialMatrixVector(divergenceCleanUpMatrix, divergenceCleaningPotential, xnumber, ynumber, znumber, 1);
	
	updateFieldByCleaning();
	double div = evaluateDivCleaningE(1, 1, 1);

	updateBoundariesNewField();
}

void Simulation::updateFieldByCleaning() {
	evaluateDivergenceCleaningField();
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				newEfield[i][j][k].x += divergenceCleaningField[i][j][k][0];
				newEfield[i][j][k].y += divergenceCleaningField[i][j][k][1];
				newEfield[i][j][k].z += divergenceCleaningField[i][j][k][2];
			}
		}
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

				divergenceCleaningField[i][j][k][0] = -0.25*(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i][prevJ][k][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[i][prevJ][prevK][0] -
														     divergenceCleaningPotential[prevI][j][k][0] - divergenceCleaningPotential[prevI][prevJ][k][0] - divergenceCleaningPotential[prevI][j][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0])/deltaX;

				divergenceCleaningField[i][j][k][1] = -0.25*(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[prevI][j][k][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[prevI][j][prevK][0] -
														     divergenceCleaningPotential[i][prevJ][k][0] - divergenceCleaningPotential[prevI][prevJ][k][0] - divergenceCleaningPotential[i][prevJ][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0])/deltaY;

				divergenceCleaningField[i][j][k][2] = -0.25*(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i][prevJ][k][0] + divergenceCleaningPotential[prevI][j][k][0] + divergenceCleaningPotential[prevI][prevJ][k][0] -
														     divergenceCleaningPotential[i][j][prevK][0] - divergenceCleaningPotential[i][prevJ][prevK][0] - divergenceCleaningPotential[prevI][j][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0])/deltaZ;
			}
		}
	}
}

void Simulation::createDivergenceCleanupInternalEquation(int i, int j, int k) {

	if(i == 0 && j == 0 && k == 0){
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
		divergenceCleanUpRightPart[i][j][k][0] = 0;
		return;
	}

	int prevJ = j - 1;
	if(prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int prevK = k - 1;
	if(prevK < 0) {
		prevK = znumber - 1;
	}
	int prevI = i - 1;
	if(prevI < 0) {
		prevI = xnumber - 1;
	}
	int nextJ = j + 1;
	if(nextJ >= ynumber) {
		nextJ = 0;
	}
	int nextK = k + 1;
	if(nextK >= znumber) {
		nextK = 0;
	}

	int nextI = i + 1;
	if(nextI >= xnumber) {
		nextI = 0;
	}

	//div for x
	MatrixElement element = MatrixElement( -0.5/deltaX2 - 0.5/deltaY2 - 0.5/deltaZ2, i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaX2 - 0.25/deltaY2 - 0.25/deltaZ2, nextI, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaX2 - 0.25/deltaY2 - 0.25/deltaZ2, prevI, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaX2 + 0.25/deltaY2 - 0.25/deltaZ2, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX2 + 0.25/deltaY2 - 0.25/deltaZ2, i, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaX2 - 0.25/deltaY2 + 0.25/deltaZ2, i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX2 - 0.25/deltaY2 + 0.25/deltaZ2, i, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.125/deltaX2 + 0.125/deltaY2 - 0.125/deltaZ2, nextI, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.125/deltaX2 + 0.125/deltaY2 - 0.125/deltaZ2, nextI, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.125/deltaX2 + 0.125/deltaY2 - 0.125/deltaZ2, prevI, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.125/deltaX2 + 0.125/deltaY2 - 0.125/deltaZ2, prevI, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.125/deltaX2 - 0.125/deltaY2 + 0.125/deltaZ2, nextI, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.125/deltaX2 - 0.125/deltaY2 + 0.125/deltaZ2, nextI, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.125/deltaX2 - 0.125/deltaY2 + 0.125/deltaZ2, prevI, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.125/deltaX2 - 0.125/deltaY2 + 0.125/deltaZ2, prevI, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.125/deltaX2 + 0.125/deltaY2 + 0.125/deltaZ2, i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.125/deltaX2 + 0.125/deltaY2 + 0.125/deltaZ2, i, prevJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.125/deltaX2 + 0.125/deltaY2 + 0.125/deltaZ2, i, nextJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.125/deltaX2 + 0.125/deltaY2 + 0.125/deltaZ2, i, prevJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.0625/deltaX2 + 0.0625/deltaY2 + 0.0625/deltaZ2, nextI, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.0625/deltaX2 + 0.0625/deltaY2 + 0.0625/deltaZ2, nextI, nextJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.0625/deltaX2 + 0.0625/deltaY2 + 0.0625/deltaZ2, nextI, prevJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.0625/deltaX2 + 0.0625/deltaY2 + 0.0625/deltaZ2, nextI, prevJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.0625/deltaX2 + 0.0625/deltaY2 + 0.0625/deltaZ2, prevI, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.0625/deltaX2 + 0.0625/deltaY2 + 0.0625/deltaZ2, prevI, nextJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.0625/deltaX2 + 0.0625/deltaY2 + 0.0625/deltaZ2, prevI, prevJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.0625/deltaX2 + 0.0625/deltaY2 + 0.0625/deltaZ2, prevI, prevJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	

	divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k);

}

void Simulation::createDivergenceCleanupLeftEquation(int j, int k) {
	int i = 0;

	int prevJ = j - 1;
	if(prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int prevK = k - 1;
	if(prevK < 0) {
		prevK = znumber - 1;
	}
	int nextJ = j + 1;
	if(nextJ >= ynumber) {
		nextJ = 0;
	}
	int nextK = k + 1;
	if(nextK >= znumber) {
		nextK = 0;
	}

	int nextI = i+1;

	//div for x

	MatrixElement element = MatrixElement( -0.25/deltaX, i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaX, nextI, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaY, i, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaY, i, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaY, nextI, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaY, nextI, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaY, i, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaY, i, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaY, nextI, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaY, nextI, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaZ, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaZ, i, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaZ, nextI, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaZ, nextI, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaZ, i, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaZ, i, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaZ, nextI, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaZ, nextI, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	divergenceCleanUpRightPart[i][j][k][0] = cleanUpRightPart(i, j, k);

	//for y

	element = MatrixElement(1, i, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	divergenceCleanUpRightPart[i][j][k][1] = 0;

	//for z;

	element = MatrixElement(1, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	divergenceCleanUpRightPart[i][j][k][2] = 0;
}

void Simulation::createDivergenceCleanupRightEquation(int j, int k) {
	int i = xnumber - 1;

	int prevJ = j - 1;
	if(prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int prevK = k - 1;
	if(prevK < 0) {
		prevK = znumber - 1;
	}
	int nextJ = j + 1;
	if(nextJ >= ynumber) {
		nextJ = 0;
	}
	int nextK = k + 1;
	if(nextK >= znumber) {
		nextK = 0;
	}

	//div for x
	MatrixElement element = MatrixElement( -0.25/deltaX, i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaY, i, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaY, i, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaY, i, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaY, i, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaZ, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaZ, i, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaZ, i, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaZ, i, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	divergenceCleanUpRightPart[i][j][k][0] = cleanUpRightPart(i, j, k);

	//rot x for y

	element = MatrixElement(-0.25/deltaY, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-0.25/deltaY, i, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	element = MatrixElement(0.25/deltaY, i, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(0.25/deltaY, i, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	element = MatrixElement(-0.25/deltaZ, i, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-0.25/deltaZ, i, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	element = MatrixElement(0.25/deltaZ, i, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(0.25/deltaZ, i, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	divergenceCleanUpRightPart[i][j][k][1] = 0;

	//rot y for z;

	element = MatrixElement(-0.25/deltaZ, i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-0.25/deltaZ, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	element = MatrixElement(0.25/deltaZ, i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaZ, i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	element = MatrixElement(0.25/deltaX, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaX, i, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaX, i, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaX, i, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	divergenceCleanUpRightPart[i][j][k][2] = 0;
}

double Simulation::cleanUpRightPart(int i, int j, int k) {
	double div = evaluateDivNewE(i, j, k);

	return 4*pi*chargeDensity[i][j][k] - div;
}