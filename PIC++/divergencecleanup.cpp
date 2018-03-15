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
#include "fourier.h"
#include "paths.h"


void Simulation::cleanupDivergence(Vector3d*** field, double*** density) {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) printf("cleaning up divergence\n");
	fflush(stdout);

	if (ynumberGeneral == 1 && znumberGeneral == 1) {
		cleanupDivergence1d(field, density);
		MPI_Barrier(cartComm);
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("cleaning divergence time = %g sec\n", procTime / CLOCKS_PER_SEC);
		}
		return;
	}

	//substractMeanChargeDensity();
	//return;

	bool fourier = false;
	bool iterative = false;

	if (!fourier) {
		if(iterative) {
			double**** tempVector1 = new double***[xnumberAdded];
			double**** tempVector = new double***[xnumberAdded];
			for (int i = 0; i < xnumberAdded; ++i) {
				tempVector1[i] = new double**[ynumberAdded];
				tempVector[i] = new double**[ynumberAdded];
				for (int j = 0; j < ynumberAdded; ++j) {
					tempVector1[i][j] = new double* [znumberAdded];
					tempVector[i][j] = new double* [znumberAdded];
					for (int k = 0; k < znumberAdded; ++k) {
						tempVector1[i][j][k] = new double[1];
						tempVector[i][j][k] = new double[1];
						for (int l = 0; l < 1; ++l) {
							tempVector1[i][j][k][l] = 0;
							tempVector[i][j][k][l] = 0;
						}
					}
				}
			}

			RightPartIterativeEvaluator* rightPartEvaluator = new PoissonRightPartEvaluator(chargeDensity, tempVector1, deltaX, deltaY, deltaZ, xnumberAdded, ynumberAdded, znumberAdded, boundaryConditionTypeX == PERIODIC, boundaryConditionTypeY == PERIODIC, boundaryConditionTypeZ == PERIODIC, rank, nprocs, cartComm, cartCoord, cartDim);

			simpleIterationSolver(divergenceCleaningPotential, tempVector, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, 1, rank, nprocs, xnumberGeneral,
	                      ynumberGeneral, znumberGeneral, maxErrorLevel, maxSimpleIterationSolverIterations,
	                      boundaryConditionTypeX == PERIODIC, boundaryConditionTypeY == PERIODIC,
	                      boundaryConditionTypeZ == PERIODIC, verbosity, leftOutGmresBuffer, rightOutGmresBuffer,
	                      leftInGmresBuffer, rightInGmresBuffer, frontOutGmresBuffer, backOutGmresBuffer,
	                      frontInGmresBuffer, backInGmresBuffer, bottomOutGmresBuffer, topOutGmresBuffer,
	                      bottomInGmresBuffer, topInGmresBuffer, rightPartEvaluator, cartComm, cartCoord, cartDim);

			delete rightPartEvaluator;

			for (int i = 0; i < xnumberAdded; ++i) {;
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						delete[] tempVector1[i][j][k];
						delete[] tempVector[i][j][k];
					}
					delete[] tempVector1[i][j];
					delete[] tempVector[i][j];
				}
				delete[] tempVector1[i];
				delete[] tempVector[i];
			}
			delete[] tempVector1;
			delete[] tempVector;
		} else {
		double foolRightPart = 0;
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

		if (cartDim[0] > 1) {
			for (int i = 0; i < xnumberAdded; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						if ((j <= additionalBinNumber) || (j >= ynumberAdded - 1 - additionalBinNumber) || (k <= additionalBinNumber) || (
							k >= znumberAdded - 1 - additionalBinNumber)) {
							createDivergenceFakeEquation(i, j, k);
						} else {
							if (cartCoord[0] == 0) {
								if (i <= additionalBinNumber) {
									createDivergenceFakeEquation(i, j, k);
								} else if (i == 1 + additionalBinNumber) {
									if (boundaryConditionTypeX == PERIODIC) {
										createDivergenceCleanupInternalEquation(i, j, k, field, density);
									} else {
										createDivergenceZeroEquation(i, j, k);
									}
								} else if (i < xnumberAdded - 1 - additionalBinNumber) {
									createDivergenceCleanupInternalEquation(i, j, k, field, density);
								} else {
									createDivergenceFakeEquation(i, j, k);
								}
							} else if (cartCoord[0] == cartDim[0] - 1) {
								if (i <= additionalBinNumber) {
									createDivergenceFakeEquation(i, j, k);
								} else if (i < xnumberAdded - 1 - additionalBinNumber) {
									createDivergenceCleanupInternalEquation(i, j, k, field, density);
								} else {
									if (boundaryConditionTypeX == PERIODIC) {
										createDivergenceFakeEquation(i, j, k);
									} else {
										createDivergenceFakeEquation(i, j, k);
										//createDivergenceZeroEquation(i, j, k);
									}
								}
							} else {
								if (i <= additionalBinNumber) {
									createDivergenceFakeEquation(i, j, k);
								} else if (i < xnumberAdded - 1 - additionalBinNumber) {
									createDivergenceCleanupInternalEquation(i, j, k, field, density);
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
						if ((j <= additionalBinNumber) || (j >= ynumberAdded - 1 - additionalBinNumber) || (k <= additionalBinNumber) || (
							k >= znumberAdded - 1 - additionalBinNumber)) {
							createDivergenceFakeEquation(i, j, k);
						} else {
							if (i <= additionalBinNumber) {
								createDivergenceFakeEquation(i, j, k);
							} else if (i == 1 + additionalBinNumber) {
								if (boundaryConditionTypeX == PERIODIC) {
									createDivergenceCleanupInternalEquation(i, j, k, field, density);
								} else {
									createDivergenceZeroEquation(i, j, k);
								}
							} else if (i < xnumberAdded - 1 - additionalBinNumber) {
								createDivergenceCleanupInternalEquation(i, j, k, field, density);
							} else {
								if (boundaryConditionTypeX == PERIODIC) {
									createDivergenceFakeEquation(i, j, k);
								} else {
									createDivergenceZeroEquation(i, j, k);
								}
							}
						}
					}
				}
			}
			double fullRightPart = 0;
			double fullDensity = 0;
			for (int i = 1 + additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
				for (int j = 1 + additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
					for (int k = 1 + additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
						fullRightPart += divergenceCleanUpRightPart[i][j][k][0];
						fullDensity += chargeDensity[i][j][k];
					}
				}
			}
			fullRightPart = fullRightPart;
		}
		if (debugMode) {
			//checkEquationMatrix(divergenceCleanUpMatrix, 1);
		}
		bool converges = true;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < 1; ++l) {
						//divergenceCleaningPotential[i][j][k][l] = chargeDensity[i][j][k]*deltaX2;
					}
				}
			}
		}
		//exchangeLargeVector(divergenceCleaningPotential, xnumberAdded, ynumberAdded, znumberAdded, 1, additionalBinNumber, boundaryConditionTypeX == PERIODIC, cartComm, cartCoord, cartDim, leftOutDivergenceBuffer, rightOutDivergenceBuffer, leftInDivergenceBuffer, rightInDivergenceBuffer, frontOutDivergenceBuffer, backOutDivergenceBuffer, frontInDivergenceBuffer, backInDivergenceBuffer, bottomOutDivergenceBuffer, topOutDivergenceBuffer, bottomInDivergenceBuffer, topInDivergenceBuffer);
		/*biconjugateStabilizedGradientMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential,
		                                    xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, 1, xnumberGeneral,
		                                    ynumberGeneral, znumberGeneral, maxCleanupErrorLevel,
		                                    maxDivergenceCleanupIterations, boundaryConditionTypeX == PERIODIC,
		                                    boundaryConditionTypeY == PERIODIC, boundaryConditionTypeZ == PERIODIC, verbosity,
		                                    cartComm, cartCoord, cartDim, residualBiconjugateDivE,
		                                    firstResidualBiconjugateDivE, vBiconjugateDivE, pBiconjugateDivE,
		                                    sBiconjugateDivE, tBiconjugateDivE, leftOutDivergenceBuffer,
		                                    rightOutDivergenceBuffer, leftInDivergenceBuffer, rightInDivergenceBuffer,
		                                    frontOutDivergenceBuffer, backOutDivergenceBuffer, frontInDivergenceBuffer,
		                                    backInDivergenceBuffer, bottomOutDivergenceBuffer, topOutDivergenceBuffer,
		                                    bottomInDivergenceBuffer, topInDivergenceBuffer, converges);
		if (!converges) {
			if (rank == 0) printf("conjugate gradients did not converge\n");

			for (int i = 0; i < xnumberAdded; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						divergenceCleaningPotential[i][j][k][0] = 0;
					}
				}
			}
		}*/
		generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumberAdded, ynumberAdded,
		znumberAdded, 1, xnumberGeneral, znumberGeneral, ynumberGeneral, additionalBinNumber, maxErrorLevel, maxDivergenceCleanupIterations, boundaryConditionTypeX == PERIODIC, boundaryConditionTypeY == PERIODIC, boundaryConditionTypeZ == PERIODIC, verbosity, leftOutDivergenceBuffer, rightOutDivergenceBuffer, leftInDivergenceBuffer, rightInDivergenceBuffer, frontOutDivergenceBuffer, backOutDivergenceBuffer, frontInDivergenceBuffer, backInDivergenceBuffer, bottomOutDivergenceBuffer, topOutDivergenceBuffer, bottomInDivergenceBuffer, topInDivergenceBuffer, gmresCleanupBasis, cartComm, cartCoord, cartDim);
		}
	} else {
		//test
		/*for(int i = 0; i < xnumberAdded; ++i) {
			for(int j = 0; j < ynumberAdded; ++j) {
				for(int k = 0; k < znumberAdded; ++k) {
					chargeDensity[i][j][k] = 0;
				}
			}
		}
		for(int i = 0; i < xnumberAdded + 1; ++i) {
			for(int j = 0; j < ynumberAdded + 1; ++j) {
				for(int k = 0; k < znumberAdded + 1; ++k) {
					newEfield[i][j][k].x = sin(2*pi*xgrid[i]/xsizeGeneral);
					newEfield[i][j][k].y = 0;
					newEfield[i][j][k].z = 0;
				}
			}
		}*/
		int* xabsoluteIndex = new int[xnumber];
		int* yabsoluteIndex = new int[ynumber];
		int* zabsoluteIndex = new int[znumber];

		for (int i = 0; i < xnumber; ++i) {
			xabsoluteIndex[i] = firstAbsoluteXindex + i;
		}
		for (int j = 0; j < ynumber; ++j) {
			yabsoluteIndex[j] = firstAbsoluteYindex + j;
		}
		for (int k = 0; k < znumber; ++k) {
			zabsoluteIndex[k] = firstAbsoluteZindex + k;
		}

		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					fourierInput[i][j][k].re = -cleanUpRightPart(1 + additionalBinNumber + i, 1 + additionalBinNumber + j,
					                                             1 + additionalBinNumber + k, field, density);
					fourierInput[i][j][k].im = 0;
				}
			}
		}

		fourierTranslation(fourierInput, fourierImage, fourierScalarTempOutput, fourierScalarTempOutput1, true, xnumber,
		                   ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY,
		                   cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);

		int halfx = xnumberGeneral / 2;
		int halfy = ynumberGeneral / 2;
		int halfz = znumberGeneral / 2;

		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					double kx = xabsoluteIndex[i];
					if (kx > halfx) {
						kx = xnumberGeneral - kx;
					}
					kx *= 2 * pi / xsizeGeneral;
					double ky = yabsoluteIndex[j];
					if (ky > halfy) {
						ky = ynumberGeneral - ky;
					}
					ky *= 2 * pi / ysizeGeneral;
					double kz = zabsoluteIndex[k];
					if (kz > halfz) {
						kz = znumberGeneral - kz;
					}
					kz *= 2 * pi / zsizeGeneral;
					if (xabsoluteIndex[i] != 0 || yabsoluteIndex[j] != 0 || zabsoluteIndex[k] != 0) {
						fourierImage[i][j][k] = fourierImage[i][j][k] / (-kx * kx - ky * ky - kz * kz);
					} else {
						fourierImage[i][j][k] = Complex(0, 0);
					}
				}
			}
		}
		fourierTranslation(fourierImage, fourierOutput, fourierScalarTempOutput, fourierScalarTempOutput1, false, xnumber,
		                   ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY,
		                   cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);

		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					divergenceCleaningPotential[1 + additionalBinNumber + i][1 + additionalBinNumber + j][1 + additionalBinNumber + k][
						0] = fourierOutput[i][j][k].re;
				}
			}
		}

		double meanPotential[1];
		meanPotential[0] = 0;
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					meanPotential[0] += divergenceCleaningPotential[1 + additionalBinNumber + i][1 + additionalBinNumber + j][1 +
						additionalBinNumber + k][0];
				}
			}
		}
		double temp[1];
		MPI_Allreduce(meanPotential, temp, 1, MPI_DOUBLE, MPI_SUM, cartComm);
		meanPotential[0] = temp[0] / (xnumberGeneral * ynumberGeneral * znumberGeneral);
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					divergenceCleaningPotential[1 + additionalBinNumber + i][1 + additionalBinNumber + j][1 + additionalBinNumber + k][
						0] -= meanPotential[0];
				}
			}
		}

		if (boundaryConditionTypeX != PERIODIC) {
			if (cartCoord[0] == 0) {
				for (int i = 0; i <= additionalBinNumber; ++i) {
					for (int j = 0; j < ynumberAdded; ++j) {
						for (int k = 0; k < znumberAdded; ++k) {
							divergenceCleaningPotential[i][j][k][0] = 0;
						}
					}
				}
			}
			if (cartCoord[0] == cartDim[0] - 1) {
				for (int i = xnumberAdded - additionalBinNumber - 1; i < xnumberAdded; ++i) {
					for (int j = 0; j < ynumberAdded; ++j) {
						for (int k = 0; k < znumberAdded; ++k) {
							divergenceCleaningPotential[i][j][k][0] = 0;
						}
					}
				}
			}
		}
		exchangeGeneralScalarCellField(divergenceCleaningPotential);

		delete[] xabsoluteIndex;
		delete[] yabsoluteIndex;
		delete[] zabsoluteIndex;
	}
	double a = cleanUpRightPart(4, 2, 2, field, density);
	Vector3d E1 = newEfield[4][2][2];
	evaluateDivergenceCleaningField();
	//substractMeanEfield();
	updateFieldByCleaning(field);
	double b = cleanUpRightPart(4, 2, 2, field, density);
	Vector3d E2 = newEfield[4][2][2];
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("cleaning divergence time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
	//double div = evaluateDivCleaningE(1);

	//updateBoundariesNewField();
}

void Simulation::cleanupDivergence1d(Vector3d*** field, double*** density) {
	MPI_Barrier(cartComm);
	double rightField[3];
	if (cartCoord[0] != cartDim[0] - 1) {
		MPI_Status status;
		MPI_Recv(rightField, 3, MPI_DOUBLE, rightRank, MPI_SEND_VECTOR_LEFT, cartComm, &status);
	} else {
		rightField[0] = 0;
		rightField[1] = 0;
		rightField[2] = 0;
	}

	for (int i = 0; i < 3; ++i) {
		divergenceCleaningField[xnumberAdded][1 + additionalBinNumber][1 + additionalBinNumber][i] = rightField[i];
	}

	for (int i = xnumberAdded - 1; i >= 0; --i) {
		double clean_up_right_part = cleanUpRightPart(i, 1 + additionalBinNumber, 1 + additionalBinNumber, field, density);
		divergenceCleaningField[i][1 + additionalBinNumber][1 + additionalBinNumber][0] = divergenceCleaningField[i + 1][1 +
			additionalBinNumber][1 + additionalBinNumber][0] - clean_up_right_part * deltaX;
		divergenceCleaningField[i][1 + additionalBinNumber][1 + additionalBinNumber][1] = 0;
		divergenceCleaningField[i][1 + additionalBinNumber][1 + additionalBinNumber][2] = 0;
		if ((cartCoord[0] == cartDim[0] - 1) && (boundaryConditionTypeX != PERIODIC)) {
			if (i >= xnumberAdded - 1 - additionalBinNumber) {
				//if(i >= xnumber){
				divergenceCleaningField[i][1 + additionalBinNumber][1 + additionalBinNumber][0] = 0;
			}
		}
	}

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					divergenceCleaningField[i][j][k][l] = divergenceCleaningField[i][1 + additionalBinNumber][1 + additionalBinNumber][
						l];
				}
			}
		}
	}

	double leftField[3];

	for (int i = 0; i < 3; ++i) {
		leftField[i] = divergenceCleaningField[2 + 2 * additionalBinNumber][1 + additionalBinNumber][1 + additionalBinNumber][
			i];
	}

	if (cartCoord[0] != 0) {
		MPI_Send(leftField, 3, MPI_DOUBLE, leftRank, MPI_SEND_VECTOR_LEFT, cartComm);
	}

	MPI_Barrier(cartComm);

	if (boundaryConditionTypeX == PERIODIC) {
		//substractMeanEfield(divergenceCleaningField);
	}

	updateFieldByCleaning(field);

	if (cartCoord[0] == 0 && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
		for (int i = 0; i <= 1 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					newEfield[i][j][k].y = 0;
					newEfield[i][j][k].z = 0;
				}
			}

		}
	}
	for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
				newEfield[xnumber][j][k] = rightBoundaryFieldEvaluator->evaluateEfield(time+deltaT, j, k);
				newEfield[xnumber + 1][j][k] = newEfield[xnumber][j][k];
			}
		}
	}

}

void Simulation::substractMeanEfield(double**** field) {
	Vector3d meanField = Vector3d(0, 0, 0);
	double meanFieldSquare = 0;

	/*int maxI = xnumber;
	if(rank == nprocs-1 && boundaryConditionTypeX != PERIODIC){
		maxI = xnumber+1;
	}*/

	int minI = 1 + additionalBinNumber;
	if (cartCoord[0] == 0 && boundaryConditionTypeX != PERIODIC) {
		minI = 0;
	}
	int maxI = xnumberAdded - 1 - additionalBinNumber;
	if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
		maxI = xnumberAdded - 1;
	}
	int minJ = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - 1 - additionalBinNumber;
	int minK = 1 + additionalBinNumber;
	int maxK = znumberAdded - 1 - additionalBinNumber;

	double meanEfield[3];

	meanEfield[0] = 0;
	meanEfield[1] = 0;
	meanEfield[2] = 0;
	for (int i = minI; i < maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				meanEfield[0] += field[i][j][k][0];
				meanEfield[1] += field[i][j][k][1];
				meanEfield[2] += field[i][j][k][2];
			}
		}
	}

	double tempE[3];

	MPI_Allreduce(meanEfield, tempE, 3, MPI_DOUBLE, MPI_SUM, cartComm);

	meanField.x = tempE[0] / (xnumberGeneral * ynumberGeneral * znumberGeneral);
	meanField.y = tempE[1] / (xnumberGeneral * ynumberGeneral * znumberGeneral);
	meanField.z = tempE[2] / (xnumberGeneral * ynumberGeneral * znumberGeneral);

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				//if (i < xnumber - additionalBinNumber ||  cartCoord[0] < cartDim[0] - 1 || boundaryConditionTypeX == PERIODIC) {
				field[i][j][k][0] -= meanField.x;
				//}
			}
		}
	}
}

void Simulation::updateFieldByCleaning(Vector3d*** field) {

	Vector3d E1 = newEfield[3][1][1];
	for (int i = 1; i < xnumberAdded; ++i) {
		for (int j = 1; j < ynumberAdded; ++j) {
			for (int k = 1; k < znumberAdded; ++k) {
				Vector3d E2 = field[i][j][k];
				field[i][j][k].x += divergenceCleaningField[i][j][k][0];
				field[i][j][k].y += divergenceCleaningField[i][j][k][1];
				field[i][j][k].z += divergenceCleaningField[i][j][k][2];
				Vector3d E3 = field[i][j][k];
				Vector3d E6 = field[i][j][k];
			}
		}
	}
	Vector3d E4 = newEfield[3][1][1];
	MPI_Barrier(cartComm);

	exchangeGeneralEfield(field);
	Vector3d E5 = newEfield[3][1][1];
}

void Simulation::evaluateDivergenceCleaningField() {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				divergenceCleaningField[i][j][k][0] = 0;
				divergenceCleaningField[i][j][k][1] = 0;
				divergenceCleaningField[i][j][k][2] = 0;
			}
		}
	}

	for (int i = 1; i < xnumberAdded; ++i) {
		for (int j = 1; j < ynumberAdded; ++j) {
			for (int k = 1; k < znumberAdded; ++k) {
				int prevI = i - 1;

				int prevJ = j - 1;

				int prevK = k - 1;

				divergenceCleaningField[i][j][k][0] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i][
						prevJ][k][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[i][prevJ][prevK][0]
					- divergenceCleaningPotential[prevI][j][k][0] - divergenceCleaningPotential[prevI][prevJ][k][0] -
					divergenceCleaningPotential[prevI][j][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0]) / (4 *
					deltaX);

				divergenceCleaningField[i][j][k][1] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[prevI]
					[j][k][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[prevI][j][prevK][0]
					- divergenceCleaningPotential[i][prevJ][k][0] - divergenceCleaningPotential[prevI][prevJ][k][0] -
					divergenceCleaningPotential[i][prevJ][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0]) / (4 *
					deltaY);

				divergenceCleaningField[i][j][k][2] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i][
						prevJ][k][0] + divergenceCleaningPotential[prevI][j][k][0] + divergenceCleaningPotential[prevI][prevJ][k][0]
					- divergenceCleaningPotential[i][j][prevK][0] - divergenceCleaningPotential[i][prevJ][prevK][0] -
					divergenceCleaningPotential[prevI][j][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0]) / (4 *
					deltaZ);
			}
		}
	}

	if (cartCoord[0] == 0 && boundaryConditionTypeX != PERIODIC) {
		for (int i = 0; i <= 1 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					divergenceCleaningField[i][j][k][0] = 0;
					divergenceCleaningField[i][j][k][1] = 0;
					divergenceCleaningField[i][j][k][2] = 0;
				}
			}
		}
	}

	if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
		for (int i = xnumberAdded - 1 - additionalBinNumber; i <= xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					divergenceCleaningField[i][j][k][0] = 0;
					divergenceCleaningField[i][j][k][1] = 0;
					divergenceCleaningField[i][j][k][2] = 0;
				}
			}

		}
	}
}

void Simulation::createDivergenceCleanupInternalEquation(int i, int j, int k, Vector3d*** field, double*** density) {
	if (cartCoord[0] == 0 && cartCoord[1] == 0 & cartCoord[2] == 0) {
		/*if (i == 1 + additionalBinNumber && j == 1 + additionalBinNumber && k == 1 + additionalBinNumber) {
			createDivergenceFixEquation(i, j, k, density);
			return;
		}*/

		/*if(i == xnumberAdded/2 && j == ynumberAdded/2 && k == znumberAdded/2){
			createDivergenceZeroEquation(i, j, k);
			return;
		}*/

		/*if(i == 2 + additionalBinNumber && j == 1 + additionalBinNumber && k == 1 + additionalBinNumber){
			createDivergenceZeroEquation(i, j, k);
			return;
		}*/

		/*if(ynumberGeneral > 1){
			if(i == 1 + additionalBinNumber && j == 2 + additionalBinNumber && k == 1 + additionalBinNumber){
			createDivergenceZeroEquation(i, j, k);
			return;
		}
		}
		if(znumberGeneral > 1){
			if(i == 1 + additionalBinNumber && j == 1 + additionalBinNumber && k == 2 + additionalBinNumber){
			createDivergenceZeroEquation(i, j, k);
			return;
		}
		}*/
	}
	/*if(cartCoord[1] == 0 && cartCoord[2] == 0 && ynumberGeneral > 1 && znumberGeneral > 1){
		if(i == 1 + additionalBinNumber && j == 1 + additionalBinNumber){
			createDivergenceZeroEquation(i, j, k);
			return;
		}
	}*/
	int prevI = i - 1;

	int nextI = i + 1;


	int prevJ = j - 1;

	int nextJ = j + 1;


	int prevK = k - 1;

	int nextK = k + 1;


	divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k, field, density);
	//divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k)*deltaX2;

	bool _27points = true;

	if ((ynumberGeneral > 1) && (znumberGeneral > 1)) {
		if (! _27points) {
			double element = -2 / deltaX2 - 2 / deltaY2 - 2 / deltaZ2;
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

			element = 1 / deltaX2;
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

			element = 1 / deltaY2;
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

			element = 1 / deltaZ2;
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
		} else {
			double element = -1 / (2 * deltaX2) - 1 / (2 * deltaY2) - 1 / (2 * deltaZ2);
			//double element = -44/(13*deltaX2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

			element = 1 / (16 * deltaX2) + 1 / (16 * deltaY2) + 1 / (16 * deltaZ2);
			//element = 1/(13*deltaX2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, nextK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, prevK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, nextK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, prevK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, nextK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, prevK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, nextK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, prevK, 0));

			element = 1 / (8 * deltaX2) + 1 / (8 * deltaY2) - 1 / (8 * deltaZ2);
			//element = 3/(26*deltaX2);
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

			//element = 3/(13*deltaX2);
			element = 1 / (4 * deltaX2) - 1 / (4 * deltaY2) - 1 / (4 * deltaZ2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

			element = - 1 / (4 * deltaX2) + 1 / (4 * deltaY2) - 1 / (4 * deltaZ2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

			element = - 1 / (4 * deltaX2) - 1 / (4 * deltaY2) + 1 / (4 * deltaZ2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
		}
	} else if (ynumberGeneral > 1) {
		if (!_27points) {
			double element = -2 / deltaX2 - 2 / deltaY2;
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

			element = 1 / deltaX2;
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

			element = 1 / deltaY2;
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));
		} else {
			double element = -1.0 / (deltaX2) - 1.0 / (deltaY2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

			element = 1.0 / (4.0 * deltaX2) + 1 / (4.0 * deltaY2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 0));

			element = 1.0 / (2.0 * deltaX2) - 1.0 / (2.0 * deltaY2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

			element = - 1.0 / (2.0 * deltaX2) + 1.0 / (2.0 * deltaY2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));
		}

	} else if (znumberGeneral > 1) {
		if (!_27points) {
			double element = -2 / deltaX2 - 2 / deltaZ2;
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

			element = 1 / deltaX2;
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

			element = 1 / deltaZ2;
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
		} else {
			double element = -1 / (deltaX2) - 1 / (deltaZ2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));


			element = 1 / (4 * deltaX2) + 1 / (4 * deltaZ2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 0));


			element = 1 / (2 * deltaX2) - 1 / (2 * deltaZ2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

			element = - 1 / (2 * deltaX2) + 1 / (2 * deltaZ2);
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
			divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
		}
	} else {
		double element = -2 / (deltaX2);
		//double element = -2;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / deltaX2;
		//element = 1;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));
	}
}

void Simulation::createDivergenceFakeEquation(int i, int j, int k) {
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
	divergenceCleanUpRightPart[i][j][k][0] = 0;
}

void Simulation::createDivergenceZeroEquation(int i, int j, int k) {
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
	divergenceCleanUpRightPart[i][j][k][0] = 0;
}

void Simulation::createDivergenceFixEquation(int i, int j, int k, double*** density) {
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
	//divergenceCleanUpRightPart[i][j][k][0] = density[i][j][k]*deltaX2;
	divergenceCleanUpRightPart[i][j][k][0] = 0;
	//divergenceCleanUpRightPart[i][j][k][0] = 1;
}

double Simulation::cleanUpRightPart(int i, int j, int k, Vector3d*** field, double*** density) {
	double div = evaluateDivEgeneral(field, i, j, k);

	return 4 * pi * density[i][j][k] - div;
	//return -4 * pi * chargeDensity[i][j][k] + div;
}

double Simulation::evaluateMeanChargeDensity() {
	double fullCharge[1];
	double temp[1];
	fullCharge[0] = 0;
	for (int i = 1 + additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
				fullCharge[0] += chargeDensity[i][j][k] * volumeB();
			}
		}
	}
	MPI_Allreduce(fullCharge, temp, 1, MPI_DOUBLE, MPI_SUM, cartComm);
	return temp[0] / (xsizeGeneral * ysizeGeneral * zsizeGeneral);
}

void Simulation::substractMeanChargeDensity() {
	double meanChargeDensity = evaluateMeanChargeDensity();
	int minI = 0;
	if (cartCoord[0] == 0 && boundaryConditionTypeX != PERIODIC) {
		minI = 1 + additionalBinNumber;
	}
	int maxI = xnumberAdded - 1;
	if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
		maxI = xnumber - 2 - additionalBinNumber;
	}
	for (int i = 0; i < znumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				chargeDensity[i][j][k] -= meanChargeDensity;
			}
		}
	}
	meanChargeDensity = evaluateMeanChargeDensity();
}

Complex*** Simulation::evaluateFourierTranslation(double*** a) {
	Complex*** result = new Complex**[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		result[i] = new Complex*[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			result[i][j] = new Complex[znumber];
			for (int k = 0; k < znumber; ++k) {
				result[i][j][k] = Complex(0, 0);

				for (int tempi = 0; tempi < xnumber; ++tempi) {
					for (int tempj = 0; tempj < ynumber; ++tempj) {
						for (int tempk = 0; tempk < znumber; ++tempk) {
							result[i][j][k] += complexExp(
								-2 * pi * ((i * tempi * 1.0 / xnumber) + (j * tempj * 1.0 / ynumber) + (k * tempk * 1.0 / znumber))) * a[tempi][
								tempj][tempk];
						}
					}
				}

				result[i][j][k] = result[i][j][k] / (xnumber * ynumber * znumber);
			}
		}
	}

	return result;
}

double*** Simulation::evaluateReverceFourierTranslation(Complex*** a) {
	double*** result = new double**[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		result[i] = new double*[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			result[i][j] = new double[znumber];
			for (int k = 0; k < znumber; ++k) {
				result[i][j][k] = 0;

				for (int tempi = 0; tempi < xnumber; ++tempi) {
					for (int tempj = 0; tempj < ynumber; ++tempj) {
						for (int tempk = 0; tempk < znumber; ++tempk) {
							double kx = 2 * pi * tempi;
							if (tempi >= xnumber / 2.0) {
								kx = -2 * pi * (tempi - xnumber + 2);
							}

							double ky = 2 * pi * tempj;
							if (tempj >= ynumber / 2.0) {
								ky = -2 * pi * (tempj - ynumber + 2);
							}

							double kz = 2 * pi * tempk;
							if (tempk >= znumber / 2.0) {
								kz = -2 * pi * (tempk - znumber + 2);
							}
							result[i][j][k] += (complexExp(((i * kx / xnumber) + (j * ky / ynumber) + (k * kz / znumber))) * a[tempi][tempj][
								tempk]).re;
						}
					}
				}
			}
		}
	}

	return result;
}

void Simulation::cleanupDivergenceMagnetic() {
	if (ynumberGeneral == 1 && znumberGeneral == 1) {
		return;
	}
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) printf("cleaning up divergence\n");
	fflush(stdout);

	double foolRightPart = 0;
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

	if (cartDim[0] > 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					if ((j <= additionalBinNumber) || (j >= ynumberAdded - 1 - additionalBinNumber) || (k <= additionalBinNumber) || (k
						>= znumberAdded - 1 - additionalBinNumber)) {
						createDivergenceFakeEquation(i, j, k);
					} else {
						if (cartCoord[0] == 0) {
							if (i <= additionalBinNumber) {
								createDivergenceFakeEquation(i, j, k);
							} else if (i == 1 + additionalBinNumber) {
								if (boundaryConditionTypeX == PERIODIC) {
									createDivergenceCleanupInternalEquationMagnetic(i, j, k);
								} else {
									createDivergenceZeroEquation(i, j, k);
								}
							} else if (i < xnumberAdded - 1 - additionalBinNumber) {
								createDivergenceCleanupInternalEquationMagnetic(i, j, k);
							} else {
								createDivergenceFakeEquation(i, j, k);
							}
						} else if (cartCoord[0] == cartDim[0] - 1) {
							if (i <= additionalBinNumber) {
								createDivergenceFakeEquation(i, j, k);
							} else if (i < xnumberAdded - 1 - additionalBinNumber) {
								createDivergenceCleanupInternalEquationMagnetic(i, j, k);
							} else {
								if (boundaryConditionTypeX == PERIODIC) {
									createDivergenceFakeEquation(i, j, k);
								} else {
									createDivergenceFakeEquation(i, j, k);
									//createDivergenceZeroEquation(i, j, k);
								}
							}
						} else {
							if (i <= additionalBinNumber) {
								createDivergenceFakeEquation(i, j, k);
							} else if (i < xnumberAdded - 1 - additionalBinNumber) {
								createDivergenceCleanupInternalEquationMagnetic(i, j, k);
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
					if ((j <= additionalBinNumber) || (j >= ynumberAdded - 1 - additionalBinNumber) || (k <= additionalBinNumber) || (k
						>= znumberAdded - 1 - additionalBinNumber)) {
						createDivergenceFakeEquation(i, j, k);
					} else {
						if (i <= additionalBinNumber) {
							createDivergenceFakeEquation(i, j, k);
						} else if (i == 1 + additionalBinNumber) {
							if (boundaryConditionTypeX == PERIODIC) {
								createDivergenceCleanupInternalEquationMagnetic(i, j, k);
							} else {
								createDivergenceZeroEquation(i, j, k);
							}
						} else if (i < xnumberAdded - 1 - additionalBinNumber) {
							createDivergenceCleanupInternalEquationMagnetic(i, j, k);
						} else {
							if (boundaryConditionTypeX == PERIODIC) {
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
	biconjugateStabilizedGradientMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential,
	                                    xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, 1, xnumberGeneral,
	                                    ynumberGeneral, znumberGeneral, maxCleanupErrorLevel,
	                                    maxDivergenceCleanupIterations, boundaryConditionTypeX == PERIODIC,
	                                    boundaryConditionTypeY == PERIODIC, boundaryConditionTypeZ == PERIODIC, verbosity,
	                                    cartComm, cartCoord, cartDim, residualBiconjugateDivE,
	                                    firstResidualBiconjugateDivE, vBiconjugateDivE, pBiconjugateDivE, sBiconjugateDivE,
	                                    tBiconjugateDivE, leftOutDivergenceBuffer, rightOutDivergenceBuffer,
	                                    leftInDivergenceBuffer, rightInDivergenceBuffer, frontOutDivergenceBuffer,
	                                    backOutDivergenceBuffer, frontInDivergenceBuffer, backInDivergenceBuffer,
	                                    bottomOutDivergenceBuffer, topOutDivergenceBuffer, bottomInDivergenceBuffer,
	                                    topInDivergenceBuffer, converges);
	//generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumberAdded, ynumberAdded,
	//znumberAdded, 1, xnumberGeneral, znumberGeneral, ynumberGeneral, additionalBinNumber, maxErrorLevel, maxDivergenceCleanupIterations, boundaryConditionTypeX == PERIODIC, verbosity, leftOutDivergenceBuffer, rightOutDivergenceBuffer, leftInDivergenceBuffer, rightInDivergenceBuffer, frontOutDivergenceBuffer, backOutDivergenceBuffer, frontInDivergenceBuffer, backInDivergenceBuffer, bottomOutDivergenceBuffer, topOutDivergenceBuffer, bottomInDivergenceBuffer, topInDivergenceBuffer, gmresCleanupBasis, cartComm, cartCoord, cartDim);

	exchangeGeneralScalarCellField(divergenceCleaningPotential);
	evaluateDivergenceCleaningFieldMagnetic();
	updateFieldByCleaningMagnetic();

	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("cleaning divergence time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::updateFieldByCleaningMagnetic() {

	for (int i = 0; i < xnumberAdded - additionalBinNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				newBfield[i][j][k].x += divergenceCleaningField[i][j][k][0];
				newBfield[i][j][k].y += divergenceCleaningField[i][j][k][1];
				newBfield[i][j][k].z += divergenceCleaningField[i][j][k][2];
			}
		}
	}

	exchangeGeneralBfield(newBfield);
}

void Simulation::evaluateDivergenceCleaningFieldMagnetic() {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				divergenceCleaningField[i][j][k][0] = 0;
				divergenceCleaningField[i][j][k][1] = 0;
				divergenceCleaningField[i][j][k][2] = 0;
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				int prevI = i - 1;
				int nextI = i + 1;

				int prevJ = j - 1;
				int nextJ = j + 1;

				int prevK = k - 1;
				int nextK = k + 1;

				divergenceCleaningField[i][j][k][0] = -(divergenceCleaningPotential[nextI][j][k][0] + divergenceCleaningPotential[
						nextI][nextJ][k][0] + divergenceCleaningPotential[nextI][j][nextK][0] + divergenceCleaningPotential[nextI][nextJ][
						nextK][0]
					- divergenceCleaningPotential[i][j][k][0] - divergenceCleaningPotential[i][nextJ][k][0] -
					divergenceCleaningPotential[i][j][nextK][0] - divergenceCleaningPotential[i][nextJ][nextK][0]) / (4 * deltaX);

				divergenceCleaningField[i][j][k][1] = -(divergenceCleaningPotential[i][nextJ][k][0] + divergenceCleaningPotential[
						nextI][nextJ][k][0] + divergenceCleaningPotential[i][nextJ][nextK][0] + divergenceCleaningPotential[nextI][nextJ][
						nextK][0]
					- divergenceCleaningPotential[i][j][k][0] - divergenceCleaningPotential[nextI][j][k][0] -
					divergenceCleaningPotential[i][j][nextK][0] - divergenceCleaningPotential[nextI][j][nextK][0]) / (4 * deltaY);

				divergenceCleaningField[i][j][k][2] = -(divergenceCleaningPotential[i][j][nextK][0] + divergenceCleaningPotential[i]
					[nextJ][nextK][0] + divergenceCleaningPotential[nextI][j][nextK][0] + divergenceCleaningPotential[nextI][nextJ][
						nextK][0]
					- divergenceCleaningPotential[i][j][k][0] - divergenceCleaningPotential[i][nextJ][k][0] -
					divergenceCleaningPotential[nextI][j][k][0] - divergenceCleaningPotential[nextI][nextJ][k][0]) / (4 * deltaZ);
			}
		}
	}

	if (cartCoord[0] == 0 && boundaryConditionTypeX != PERIODIC) {
		for (int i = 0; i < 1 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					divergenceCleaningField[i][j][k][0] = 0;
					divergenceCleaningField[i][j][k][1] = 0;
					divergenceCleaningField[i][j][k][2] = 0;
				}
			}
		}
	}

	if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
		for (int i = xnumberAdded - 1 - additionalBinNumber; i <= xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					divergenceCleaningField[i][j][k][0] = 0;
					divergenceCleaningField[i][j][k][1] = 0;
					divergenceCleaningField[i][j][k][2] = 0;
				}
			}

		}
	}
}

void Simulation::createDivergenceCleanupInternalEquationMagnetic(int i, int j, int k) {
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
	//if(prevI < 1 +additionalBinNumber){
	//	prevI = xnumberAdded - 2 - additionalBinNumber;
	//}

	int nextI = i + 1;
	//if(nextI > xnumberAdded - 2 - additionalBinNumber){
	//	nextI = 1 + additionalBinNumber;
	//}


	int prevJ = j - 1;

	int nextJ = j + 1;


	int prevK = k - 1;

	int nextK = k + 1;


	divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPartMagnetic(i, j, k);
	//divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k)*deltaX2;

	if ((ynumberGeneral > 1) && (znumberGeneral > 1)) {
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
	} else if (ynumberGeneral > 1) {
		double element = -1.0 / (deltaX2) - 1.0 / (deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1.0 / (4.0 * deltaX2) + 1 / (4.0 * deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 0));

		element = 1.0 / (2.0 * deltaX2) - 1.0 / (2.0 * deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

		element = - 1.0 / (2.0 * deltaX2) + 1.0 / (2.0 * deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

	} else if (znumberGeneral > 1) {
		double element = -1 / (deltaX2) - 1 / (deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));


		element = 1 / (4 * deltaX2) + 1 / (4 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 0));


		element = 1 / (2 * deltaX2) - 1 / (2 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

		element = - 1 / (2 * deltaX2) + 1 / (2 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
	} else {
		double element = -2 / (deltaX2);
		//double element = -2;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / deltaX2;
		//element = 1;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));
	}
}

double Simulation::cleanUpRightPartMagnetic(int i, int j, int k) {
	double div = evaluateDivNewB(i, j, k);

	return - div;
}
