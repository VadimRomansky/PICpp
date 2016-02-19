#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "specialmath.h"


void Simulation::cleanupDivergence() {
	printf("cleaning up divergence\n");


	int matrixDimension = xnumber;

	double fullDensity = 0;
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				fullDensity += chargeDensity[i][j][k] * volumeB(i, j, k);
			}
		}
	}
	fullDensity /= (xsize*ysize*zsize);

	if(boundaryConditionType == PERIODIC){
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					chargeDensity[i][j][k] -= fullDensity;
				}
			}
		}
		/*for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					fullDensity += chargeDensity[i][j][k] * volumeB(i, j, k);
				}
			}
		}
		fullDensity /= (xsize*ysize*zsize);*/
	} else {
		double Elinear = -4*pi*fullDensity*xsize + newEfield[xnumber][0][0].x - newEfield[0][0][0].x;
		for(int i = 0; i < xnumber + 1; ++i){
			for(int j = 0; j< ynumber + 1; ++j){
				for(int k = 0; k < znumber + 1; ++k){
					double factor = (xgrid[xnumber] - xgrid[i])/xsize;
					newEfield[i][j][k].x = newEfield[i][j][k].x + Elinear*factor;
				}
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				chargeDensity[i][j][k] -= fullDensity;
			}
		}
	}

	/*for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				chargeDensity[i][j][k] -= fullDensity;
			}
		}
	}*/


	/*for (int i = 0; i < xnumber; ++i) {
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
	}*/


	if (debugMode) {
		//checkEquationMatrix(divergenceCleanUpMatrix, 3);
	}

	/*double summRightPart = 0;
	for(int i = 0; i < xnumber; ++i){
		summRightPart += divergenceCleanUpRightPart[i][0];
	}*/
	//double rightPart = cleanUpRightPart(1);
	//generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumber, ynumber, znumber, 1, 1E-100, 40000);
	//generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningField, xnumber, ynumber, znumber, 3);


	//double** leftPart = multiplySpecialMatrixVector(divergenceCleanUpMatrix, divergenceCleaningPotential, xnumber, 1);

	double*** rightPart = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		rightPart[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			rightPart[i][j] = new double[znumber];
			for(int k = 0; k < znumber; ++k){
				rightPart[i][j][k] = -cleanUpRightPart(i, j, k);
			}
		}
	}

	Complex*** rightPartFourier = evaluateFourierTranslation(rightPart);

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				if(i == 0 && j == 0 && k == 0){
					rightPartFourier[i][j][k] = Complex(0, 0);
				} else {
					double kx = i;
					if(i >= xnumber/2.0){
						kx = xnumber - i - 2;
					} 
					double ky = j;
					if(j >= ynumber/2.0){
						ky = ynumber - j - 2;
					}
					double kz = k;
					if(k >= znumber/2.0){
						kz = znumber - k - 2;
					}
					rightPartFourier[i][j][k] = rightPartFourier[i][j][k]*(-4*pi*pi*((kx*kx/(xsize*xsize))  + (ky*ky/(ysize*ysize)) + (kz*kz/(zsize*zsize))));
				}
			}
		}
	}

	double*** potential = evaluateReverceFourierTranslation(rightPartFourier);

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				divergenceCleaningPotential[i][j][k][0] = potential[i][j][k];
			}
			delete[] potential[i][j];
			delete[] rightPartFourier[i][j];
			delete[] rightPart[i][j];
		}
		delete[] potential[i];
		delete[] rightPartFourier[i];
		delete[] rightPart[i];
	}
	delete[] potential;
	delete[] rightPartFourier;
	delete[] rightPart;

	updateFieldByCleaning();
	//double div = evaluateDivCleaningE(1);

	updateBoundariesNewField();
}

void Simulation::updateFieldByCleaning() {
	evaluateDivergenceCleaningField();
	/*if((ynumber == 1) && (znumber == 1)){
		divergenceCleaningField[xnumber - 1][0][0][0] = 0;
		divergenceCleaningField[xnumber - 1][0][0][1] = 0;
		divergenceCleaningField[xnumber - 1][0][0][2] = 0;

		for (int i = xnumber - 2; i >= 0; --i) {
			divergenceCleaningField[i][0][0][0] = divergenceCleaningField[i + 1][0][0][0] - cleanUpRightPart(i, 0, 0) * deltaX;
			divergenceCleaningField[i][0][0][1] = 0;
			divergenceCleaningField[i][0][0][2] = 0;
		}
	}*/
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
					//if(boundaryConditionType == PERIODIC){
						prevI = xnumber - 1;
					//} else {
						//prevI = 0;
					//}
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

	/*Vector3d constantField;
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
	}*/
}

void Simulation::createDivergenceCleanupInternalEquation(int i, int j, int k) {

	if (i == 0 && j == 0 && k == 0) {
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
		divergenceCleanUpRightPart[i][j][k][0] = 1;
		return;
	}

	if (i == xnumber - 1 && j == 0 && k == 0) {
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(-1, i-1, j, k, 0));
		divergenceCleanUpRightPart[i][j][k][0] = 0;
		return;
	}

	/*if (i == xnumber - 1 && j == ynumber - 1 && k == znumber - 1) {
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
		divergenceCleanUpRightPart[i][j][k][0] = 0;
		return;
	}*/

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


	//divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k);
	divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k)*deltaX2;

	if((ynumber > 1) && (znumber > 1)){
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
	} else if(ynumber > 1){
	} else if(znumber > 1){
	} else {
		//double element = -2 / (deltaX2);
		double element = -2;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		//element = 1/deltaX2;
		element = 1;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));
	}
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

Complex*** Simulation::evaluateFourierTranslation(double*** a){
	Complex*** result = new Complex**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		result[i] = new Complex*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			result[i][j] = new Complex[znumber];
			for(int k = 0; k < znumber; ++k){
				result[i][j][k] = Complex(0, 0);

				for(int tempi = 0; tempi < xnumber; ++tempi){
					for(int tempj = 0; tempj < ynumber; ++tempj){
						for(int tempk = 0; tempk < znumber; ++tempk){
							result[i][j][k] += complexExp(-2*pi*((i*tempi*1.0/xnumber) + (j*tempj*1.0/ynumber) + (k*tempk*1.0/znumber)))*a[tempi][tempj][tempk];
						}
					}
				}

				result[i][j][k] = result[i][j][k]/(xnumber*ynumber*znumber);
			}
		}
	}

	return result;
}

double*** Simulation::evaluateReverceFourierTranslation(Complex*** a){
	double*** result = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		result[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			result[i][j] = new double[znumber];
			for(int k = 0; k < znumber; ++k){
				result[i][j][k] = 0;

				for(int tempi = 0; tempi < xnumber; ++tempi){
					for(int tempj = 0; tempj < ynumber; ++tempj){
						for(int tempk = 0; tempk < znumber; ++tempk){
							double kx = 2*pi*tempi;
							if(tempi >= xnumber/2.0){
								kx = -2*pi*(tempi - xnumber + 2);
							}

							double ky = 2*pi*tempj;
							if(tempj >= ynumber/2.0){
								ky = -2*pi*(tempj - ynumber + 2);
							}

							double kz = 2*pi*tempk;
							if(tempk >= znumber/2.0){
								kz = -2*pi*(tempk - znumber + 2);
							}
							result[i][j][k] += (complexExp(((i*kx/xnumber) + (j*ky/ynumber) + (k*kz/znumber)))*a[tempi][tempj][tempk]).re;
						}
					}
				}
			}
		}
	}

	return result;
}
