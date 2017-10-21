#include <mpi.h>
#include <time.h>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "constants.h"
#include "specialmath.h"
#include "util.h"
#include "matrixElement.h"
#include "vector3d.h"
#include "matrix3d.h"
#include "particle.h"
#include "complex.h"
#include "simulation.h"
#include "mpi_util.h"
#include "paths.h"

void Simulation::cleanupDivergenceBuneman() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) printf("cleaning up divergence\n");
	fflush(stdout);

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				divergenceCleanUpMatrix[i][j][k][0].clear();
				divergenceCleanUpMatrix[i][j][k][1].clear();
				divergenceCleanUpMatrix[i][j][k][2].clear();
				divergenceCleanUpRightPart[i][j][k][0] = 0;
				divergenceCleanUpRightPart[i][j][k][1] = 0;
				divergenceCleanUpRightPart[i][j][k][2] = 0;
				divergenceCleaningPotential[i][j][k][0] = 0;
				divergenceCleaningField[i][j][k][0] = 0;
				divergenceCleaningField[i][j][k][1] = 0;
				divergenceCleaningField[i][j][k][2] = 0;
			}
		}
	}

	for(int i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				bunemanDivergenceCleaningPotential[i][j][k][0] = 0;
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < maxwellEquationMatrixSize; ++l) {
					gmresOutput[i][j][k][l] = 0;
				}
			}
		}
	}

	for(int  i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				bunemanDivCleaningEx[i][j][k] = 0;
			}
		}
	}

	for(int  i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				bunemanDivCleaningEy[i][j][k] = 0;
			}
		}
	}

	for(int  i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				bunemanDivCleaningEz[i][j][k] = 0;
			}
		}
	}

	if (cartDim[0] > 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					if ((j <= additionalBinNumber) || (j >= ynumberAdded - 1 - additionalBinNumber) || (k <= additionalBinNumber) || (k >= znumberAdded - 1 - additionalBinNumber)) {
						createDivergenceFakeEquation(i, j, k);
					} else {
						if (cartCoord[0] == 0) {
							if (i <= additionalBinNumber) {
								createDivergenceFakeEquation(i, j, k);
							} else if (i == 1 + additionalBinNumber) {
								if (boundaryConditionType == PERIODIC) {
									createDivergenceCleanupInternalEquationBuneman(i, j, k);
								} else {
									createDivergenceZeroEquation(i, j, k);
								}
							} else if (i < xnumberAdded - 1 - additionalBinNumber) {
								createDivergenceCleanupInternalEquationBuneman(i, j, k);
							} else {
								createDivergenceFakeEquation(i, j, k);
							}
						} else if (cartCoord[0] == cartDim[0] - 1) {
							if (i <= additionalBinNumber) {
								createDivergenceFakeEquation(i, j, k);
							} else if (i < xnumberAdded - 1 - additionalBinNumber) {
								createDivergenceCleanupInternalEquationBuneman(i, j, k);
							} else {
								if (boundaryConditionType == PERIODIC) {
									createDivergenceFakeEquation(i, j, k);
								} else {
									//createDivergenceFakeEquation(i, j, k);
									createDivergenceZeroEquation(i, j, k);
								}
							}
						} else {
							if (i <= additionalBinNumber) {
								createDivergenceFakeEquation(i, j, k);
							} else if (i < xnumberAdded - 1 - additionalBinNumber) {
								createDivergenceCleanupInternalEquationBuneman(i, j, k);
							} else {
								createDivergenceFakeEquation(i, j, k);
							}
						}
					}
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					if ((j <= additionalBinNumber) || (j >= ynumberAdded - 1 - additionalBinNumber) || (k <= additionalBinNumber) || (k >= znumberAdded - 1 - additionalBinNumber)) {
						createDivergenceFakeEquation(i, j, k);
					} else {
						if (i <= additionalBinNumber) {
							createDivergenceFakeEquation(i, j, k);
						} else if (i == 1 + additionalBinNumber) {
							if (boundaryConditionType == PERIODIC) {
								createDivergenceCleanupInternalEquationBuneman(i, j, k);
							} else {
								createDivergenceZeroEquation(i, j, k);
							}
						} else if (i < xnumberAdded - 1 - additionalBinNumber) {
							createDivergenceCleanupInternalEquationBuneman(i, j, k);
						} else {
							if (boundaryConditionType == PERIODIC) {
								createDivergenceFakeEquation(i, j, k);
							} else {
								createDivergenceZeroEquation(i, j, k);
							}
						}
					}
				}
			}
		}
	}
	//if (debugMode) {
		checkEquationMatrix(divergenceCleanUpMatrix, 1);
	//}
		bool converges;
	biconjugateStabilizedGradientMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, gmresOutput, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, 1, xnumberGeneral, ynumberGeneral, znumberGeneral, maxCleanupErrorLevel, maxDivergenceCleanupIterations, boundaryConditionType == PERIODIC, verbosity, cartComm, cartCoord, cartDim, residualBiconjugateDivE, firstResidualBiconjugateDivE, vBiconjugateDivE, pBiconjugateDivE, sBiconjugateDivE, tBiconjugateDivE, leftOutDivergenceBuffer, rightOutDivergenceBuffer, leftInDivergenceBuffer, rightInDivergenceBuffer, frontOutDivergenceBuffer, backOutDivergenceBuffer, frontInDivergenceBuffer, backInDivergenceBuffer, bottomOutDivergenceBuffer, topOutDivergenceBuffer, bottomInDivergenceBuffer, topInDivergenceBuffer, converges);
	//generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, gmresOutput, xnumberAdded, ynumberAdded,
		                                 //znumberAdded, 1, xnumberGeneral, znumberGeneral, ynumberGeneral, additionalBinNumber, maxErrorLevel, maxDivergenceCleanupIterations, boundaryConditionType == PERIODIC, verbosity, leftOutDivergenceBuffer, rightOutDivergenceBuffer, leftInDivergenceBuffer, rightInDivergenceBuffer, frontOutDivergenceBuffer, backOutDivergenceBuffer, frontInDivergenceBuffer, backInDivergenceBuffer, bottomOutDivergenceBuffer, topOutDivergenceBuffer, bottomInDivergenceBuffer, topInDivergenceBuffer, gmresCleanupBasis, cartComm, cartCoord, cartDim);

	for(int i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				bunemanDivergenceCleaningPotential[i][j][k][0] = gmresOutput[i][j][k][0];
			}
		}
	}

	if(cartCoord[0] == cartDim[0] - 1 && boundaryConditionType != PERIODIC){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				bunemanDivergenceCleaningPotential[xnumberAdded][j][k][0] = 0;
			}
		}
	}

	exchangeGeneralScalarNodeField(bunemanDivergenceCleaningPotential);

	evaluateDivergenceCleaningFieldBuneman();
	//substractMeanEfield();
	updateFieldByCleaningBuneman();

	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("cleaning divergence time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
	//double div = evaluateDivCleaningE(1);

	//updateBoundariesNewField();
}

void Simulation::createDivergenceCleanupInternalEquationBuneman(int i, int j, int k) {
	/*if(cartCoord[0] == 0 && cartCoord[1] == 0 & cartCoord[2] == 0){
	if(i == 1 + additionalBinNumber && j == 1 + additionalBinNumber && k == 1 + additionalBinNumber){
		createDivergenceZeroEquation(i, j, k);
		return;
	}

	if(i == 2 + additionalBinNumber && j == 1 + additionalBinNumber && k == 1 + additionalBinNumber){
		createDivergenceZeroEquation(i, j, k);
		return;
	}

	if(cartDim[1] > 1){
		if(i == 1 + additionalBinNumber && j == 2 + additionalBinNumber && k == 1 + additionalBinNumber){
		createDivergenceZeroEquation(i, j, k);
		return;
	}
	}
	if(cartDim[2] > 1){
		if(i == 1 + additionalBinNumber && j == 1 + additionalBinNumber && k == 2 + additionalBinNumber){
		createDivergenceZeroEquation(i, j, k);
		return;
	}
	}
	}*/
	int prevI = i - 1;

	int nextI = i + 1;


	int prevJ = j - 1;

	int nextJ = j + 1;


	int prevK = k - 1;

	int nextK = k + 1;


	divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPartBuneman(i, j, k);
	//divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k)*deltaX2;

	if ((ynumberGeneral > 1) && (znumberGeneral > 1)) {
		double element = -2 / (deltaX2) - 2 / (deltaY2) - 2 / (deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / (deltaX2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

		element = 1 / (deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

		element = 1 / (deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));

		
	} else if (ynumberGeneral > 1) {
		double element = -2 / (deltaX2) - 2 / (deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / (deltaX2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

		element = 1 / (deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

	} else if (znumberGeneral > 1) {
		double element = -2 / (deltaX2) - 2 / (deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / (deltaX2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

		element = 1 / (deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
	} else {
		double element = -2 / (deltaX2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / deltaX2;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));
	}
}

void Simulation::evaluateDivergenceCleaningFieldBuneman() {
	for(int  i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				bunemanDivCleaningEx[i][j][k] = 0;
			}
		}
	}

	for(int  i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				bunemanDivCleaningEy[i][j][k] = 0;
			}
		}
	}

	for(int  i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				bunemanDivCleaningEz[i][j][k] = 0;
			}
		}
	}

	for (int i = 1; i < xnumberAdded; ++i) {
		for (int j = 1; j < ynumberAdded; ++j) {
			for (int k = 1; k < znumberAdded; ++k) {
				int prevI = i - 1;
				int nextI = i + 1;

				int prevJ = j - 1;
				int nextJ = j + 1;

				int prevK = k - 1;
				int nextK = k + 1;

				bunemanDivCleaningEx[i][j][k] = - (bunemanDivergenceCleaningPotential[nextI][j][k][0] - bunemanDivergenceCleaningPotential[i][j][k][0])/deltaX;
				bunemanDivCleaningEy[i][j][k] = - (bunemanDivergenceCleaningPotential[i][nextJ][k][0] - bunemanDivergenceCleaningPotential[i][j][k][0])/deltaY;
				bunemanDivCleaningEz[i][j][k] = - (bunemanDivergenceCleaningPotential[i][j][nextK][0] - bunemanDivergenceCleaningPotential[i][j][k][0])/deltaZ;
			}
		}
	}

	if (cartCoord[0] == 0 && boundaryConditionType != PERIODIC) {
		for (int i = 0; i <= 1 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					if(i < 1 + additionalBinNumber){
						bunemanDivCleaningEx[i][j][k] = 0;
					}
					bunemanDivCleaningEy[i][j][k] = 0;
					bunemanDivCleaningEz[i][j][k] = 0;
				}
			}
		}
	}

	if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionType != PERIODIC) {
		for (int i = xnumberAdded - 1 - additionalBinNumber; i <= xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					if(i < xnumberAdded){
						bunemanDivCleaningEx[i][j][k] = 0;
					}
					bunemanDivCleaningEy[i][j][k] = 0;
					bunemanDivCleaningEz[i][j][k] = 0;
				}
			}
		}
	}
}

void Simulation::updateFieldByCleaningBuneman() {
	exchangeBunemanEfield(bunemanDivCleaningEx, bunemanDivCleaningEy, bunemanDivCleaningEz);
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanNewEx[i][j][k] += bunemanDivCleaningEx[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanNewEy[i][j][k] += bunemanDivCleaningEy[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				bunemanNewEz[i][j][k] += bunemanDivCleaningEz[i][j][k];
			}
		}
	}

	exchangeBunemanEfield(bunemanNewEx, bunemanNewEy, bunemanNewEz);
}

double Simulation::cleanUpRightPartBuneman(int i, int j, int k) {
	double div = evaluateDivBunemanNewE(i, j, k);

	return 4 * pi * bunemanChargeDensity[i][j][k] - div;
	//return -4 * pi * bunemanChargeDensity[i][j][k] + div;
}

void Simulation::cleanupDivergenceBunemanMagnetic() {
	if(ynumberGeneral == 1 && znumberGeneral == 1){
		return;
	}
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) printf("cleaning up divergence\n");
	fflush(stdout);

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				divergenceCleanUpMatrix[i][j][k][0].clear();
				divergenceCleanUpMatrix[i][j][k][1].clear();
				divergenceCleanUpMatrix[i][j][k][2].clear();
				divergenceCleanUpRightPart[i][j][k][0] = 0;
				divergenceCleanUpRightPart[i][j][k][1] = 0;
				divergenceCleanUpRightPart[i][j][k][2] = 0;
				divergenceCleaningPotential[i][j][k][0] = 0;
				divergenceCleaningField[i][j][k][0] = 0;
				divergenceCleaningField[i][j][k][1] = 0;
				divergenceCleaningField[i][j][k][2] = 0;
			}
		}
	}

	for(int i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				divergenceCleaningPotential[i][j][k][0] = 0;
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < maxwellEquationMatrixSize; ++l) {
					gmresOutput[i][j][k][l] = 0;
				}
			}
		}
	}

	for(int  i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				bunemanDivCleaningBx[i][j][k] = 0;
			}
		}
	}

	for(int  i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				bunemanDivCleaningBy[i][j][k] = 0;
			}
		}
	}

	for(int  i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				bunemanDivCleaningBz[i][j][k] = 0;
			}
		}
	}

	if (cartDim[0] > 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					if ((j <= additionalBinNumber) || (j >= ynumberAdded - 1 - additionalBinNumber) || (k <= additionalBinNumber) || (k >= znumberAdded - 1 - additionalBinNumber)) {
						createDivergenceFakeEquation(i, j, k);
					} else {
						if (cartCoord[0] == 0) {
							if (i <= additionalBinNumber) {
								createDivergenceFakeEquation(i, j, k);
							} else if (i == 1 + additionalBinNumber) {
								if (boundaryConditionType == PERIODIC) {
									createDivergenceCleanupInternalEquationBunemanMagnetic(i, j, k);
								} else {
									createDivergenceZeroEquation(i, j, k);
								}
							} else if (i < xnumberAdded - 1 - additionalBinNumber) {
								createDivergenceCleanupInternalEquationBunemanMagnetic(i, j, k);
							} else {
								createDivergenceFakeEquation(i, j, k);
							}
						} else if (cartCoord[0] == cartDim[0] - 1) {
							if (i <= additionalBinNumber) {
								createDivergenceFakeEquation(i, j, k);
							} else if (i < xnumberAdded - 1 - additionalBinNumber) {
								createDivergenceCleanupInternalEquationBunemanMagnetic(i, j, k);
							} else {
								if (boundaryConditionType == PERIODIC) {
									createDivergenceFakeEquation(i, j, k);
								} else {
									//createDivergenceFakeEquation(i, j, k);
									createDivergenceZeroEquation(i, j, k);
								}
							}
						} else {
							if (i <= additionalBinNumber) {
								createDivergenceFakeEquation(i, j, k);
							} else if (i < xnumberAdded - 1 - additionalBinNumber) {
								createDivergenceCleanupInternalEquationBunemanMagnetic(i, j, k);
							} else {
								createDivergenceFakeEquation(i, j, k);
							}
						}
					}
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					if ((j <= additionalBinNumber) || (j >= ynumberAdded - 1 - additionalBinNumber) || (k <= additionalBinNumber) || (k >= znumberAdded - 1 - additionalBinNumber)) {
						createDivergenceFakeEquation(i, j, k);
					} else {
						if (i <= additionalBinNumber) {
							createDivergenceFakeEquation(i, j, k);
						} else if (i == 1 + additionalBinNumber) {
							if (boundaryConditionType == PERIODIC) {
								createDivergenceCleanupInternalEquationBunemanMagnetic(i, j, k);
							} else {
								createDivergenceZeroEquation(i, j, k);
							}
						} else if (i < xnumberAdded - 1 - additionalBinNumber) {
							createDivergenceCleanupInternalEquationBunemanMagnetic(i, j, k);
						} else {
							if (boundaryConditionType == PERIODIC) {
								createDivergenceFakeEquation(i, j, k);
							} else {
								createDivergenceZeroEquation(i, j, k);
							}
						}
					}
				}
			}
		}
	}
	//if (debugMode) {
		checkEquationMatrix(divergenceCleanUpMatrix, 1);
	//}
		bool converges;
	biconjugateStabilizedGradientMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, gmresOutput, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, 1, xnumberGeneral, ynumberGeneral, znumberGeneral, maxCleanupErrorLevel, maxDivergenceCleanupIterations, boundaryConditionType == PERIODIC, verbosity, cartComm, cartCoord, cartDim, residualBiconjugateDivE, firstResidualBiconjugateDivE, vBiconjugateDivE, pBiconjugateDivE, sBiconjugateDivE, tBiconjugateDivE, leftOutDivergenceBuffer, rightOutDivergenceBuffer, leftInDivergenceBuffer, rightInDivergenceBuffer, frontOutDivergenceBuffer, backOutDivergenceBuffer, frontInDivergenceBuffer, backInDivergenceBuffer, bottomOutDivergenceBuffer, topOutDivergenceBuffer, bottomInDivergenceBuffer, topInDivergenceBuffer, converges);
	//generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, gmresOutput, xnumberAdded, ynumberAdded,
		                                 //znumberAdded, 1, xnumberGeneral, znumberGeneral, ynumberGeneral, additionalBinNumber, maxErrorLevel, maxDivergenceCleanupIterations, boundaryConditionType == PERIODIC, verbosity, leftOutDivergenceBuffer, rightOutDivergenceBuffer, leftInDivergenceBuffer, rightInDivergenceBuffer, frontOutDivergenceBuffer, backOutDivergenceBuffer, frontInDivergenceBuffer, backInDivergenceBuffer, bottomOutDivergenceBuffer, topOutDivergenceBuffer, bottomInDivergenceBuffer, topInDivergenceBuffer, gmresCleanupBasis, cartComm, cartCoord, cartDim);

	for(int i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				divergenceCleaningPotential[i][j][k][0] = gmresOutput[i][j][k][0];
			}
		}
	}

	if(cartCoord[0] == cartDim[0] - 1 && boundaryConditionType != PERIODIC){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				divergenceCleaningPotential[xnumberAdded][j][k][0] = 0;
			}
		}
	}

	exchangeGeneralScalarCellField(divergenceCleaningPotential);

	evaluateDivergenceCleaningFieldBunemanMagnetic();
	//substractMeanEfield();
	updateFieldByCleaningBunemanMagnetic();

	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("cleaning divergence time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
	//double div = evaluateDivCleaningE(1);

	//updateBoundariesNewField();
}

void Simulation::createDivergenceCleanupInternalEquationBunemanMagnetic(int i, int j, int k) {
	/*if(cartCoord[0] == 0 && cartCoord[1] == 0 & cartCoord[2] == 0){
	if(i == 1 + additionalBinNumber && j == 1 + additionalBinNumber && k == 1 + additionalBinNumber){
		createDivergenceZeroEquation(i, j, k);
		return;
	}

	if(i == 2 + additionalBinNumber && j == 1 + additionalBinNumber && k == 1 + additionalBinNumber){
		createDivergenceZeroEquation(i, j, k);
		return;
	}

	if(cartDim[1] > 1){
		if(i == 1 + additionalBinNumber && j == 2 + additionalBinNumber && k == 1 + additionalBinNumber){
		createDivergenceZeroEquation(i, j, k);
		return;
	}
	}
	if(cartDim[2] > 1){
		if(i == 1 + additionalBinNumber && j == 1 + additionalBinNumber && k == 2 + additionalBinNumber){
		createDivergenceZeroEquation(i, j, k);
		return;
	}
	}
	}*/
	int prevI = i - 1;

	int nextI = i + 1;


	int prevJ = j - 1;

	int nextJ = j + 1;


	int prevK = k - 1;

	int nextK = k + 1;


	divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPartBunemanMagnetic(i, j, k);
	//divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k)*deltaX2;

	if ((ynumberGeneral > 1) && (znumberGeneral > 1)) {
		double element = -2 / (deltaX2) - 2 / (deltaY2) - 2 / (deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / (deltaX2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

		element = 1 / (deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

		element = 1 / (deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));

		
	} else if (ynumberGeneral > 1) {
		double element = -2 / (deltaX2) - 2 / (deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / (deltaX2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

		element = 1 / (deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

	} else if (znumberGeneral > 1) {
		double element = -2 / (deltaX2) - 2 / (deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / (deltaX2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

		element = 1 / (deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
	} else {
		double element = -2 / (deltaX2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / deltaX2;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));
	}
}

void Simulation::evaluateDivergenceCleaningFieldBunemanMagnetic() {
	for(int  i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				bunemanDivCleaningBx[i][j][k] = 0;
			}
		}
	}

	for(int  i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				bunemanDivCleaningBy[i][j][k] = 0;
			}
		}
	}

	for(int  i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				bunemanDivCleaningBz[i][j][k] = 0;
			}
		}
	}

	for (int i = 1; i < xnumberAdded; ++i) {
		for (int j = 1; j < ynumberAdded; ++j) {
			for (int k = 1; k < znumberAdded; ++k) {
				int prevI = i - 1;
				int nextI = i + 1;

				int prevJ = j - 1;
				int nextJ = j + 1;

				int prevK = k - 1;
				int nextK = k + 1;

				bunemanDivCleaningBx[i][j][k] = - (divergenceCleaningPotential[i][j][k][0] - divergenceCleaningPotential[prevI][j][k][0])/deltaX;
				bunemanDivCleaningBy[i][j][k] = - (divergenceCleaningPotential[i][j][k][0] - divergenceCleaningPotential[i][prevJ][k][0])/deltaY;
				bunemanDivCleaningBz[i][j][k] = - (divergenceCleaningPotential[i][j][k][0] - divergenceCleaningPotential[i][j][prevK][0])/deltaZ;
			}
		}
	}

	if (cartCoord[0] == 0 && boundaryConditionType != PERIODIC) {
		for (int i = 0; i <= 1 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanDivCleaningBx[i][j][k] = 0;
					if(i < 1 + additionalBinNumber){
						bunemanDivCleaningBy[i][j][k] = 0;
						bunemanDivCleaningBz[i][j][k] = 0;
					}
				}
			}
		}
	}

	if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionType != PERIODIC) {
		for (int i = xnumberAdded - 1 - additionalBinNumber; i <= xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanDivCleaningBx[i][j][k] = 0;
					if(i < xnumberAdded){
						bunemanDivCleaningBy[i][j][k] = 0;
						bunemanDivCleaningBz[i][j][k] = 0;
					}
				}
			}
		}
	}
}

void Simulation::updateFieldByCleaningBunemanMagnetic() {
	exchangeBunemanBfield(bunemanDivCleaningBx, bunemanDivCleaningBy, bunemanDivCleaningBz);
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				bunemanNewBx[i][j][k] += bunemanDivCleaningBx[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				bunemanNewBy[i][j][k] += bunemanDivCleaningBy[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanNewBz[i][j][k] += bunemanDivCleaningBz[i][j][k];
			}
		}
	}

	exchangeBunemanBfield(bunemanNewBx, bunemanNewBy, bunemanNewBz);
}

double Simulation::cleanUpRightPartBunemanMagnetic(int i, int j, int k) {
	double div = evaluateDivBunemanNewB(i, j, k);

	return - div;
}