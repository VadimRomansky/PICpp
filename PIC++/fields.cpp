#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <mpi.h>
#include <time.h>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "specialmath.h"
#include "util.h"
#include "matrixElement.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "particle.h"
#include "simulation.h"
#include "mpi_util.h"
#include "output.h"
#include "paths.h"

void Simulation::evaluateElectricField() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) printf("evaluating fields\n");
	fflush(stdout);
	if ((rank == 0) && (verbosity > 0)) printLog("evaluating fields\n");
	if (solverType == BUNEMAN) {
		evaluateBunemanElectricField();
	} else {
		//fopen("./output/outputEverythingFile.dat","a");

		if (solverType == IMPLICIT || solverType == IMPLICIT_EC) {
			evaluateMaxwellEquationMatrix();
			//outputVectorCellArray((outputDir + "rightPart.dat").c_str(), maxwellEquationRightPart, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, cartCoord, cartDim);
			//maxwellMatrixFile = fopen("./output/maxwellMatrixFile.dat", "w");
			//outputMaxwellEquationMatrixFull(maxwellMatrixFile, maxwellEquationMatrix, xnumber, ynumber, znumber, maxwellEquationMatrixSize);
			//fclose(maxwellMatrixFile);
			//outputMaxwellEquationMatrixSimple(maxwellEquationMatrix, xnumber, ynumber, znumber, maxwellEquationMatrixSize);

			//FILE* rightPartFile = fopen("./output/rightPartFile.dat", "w");
			//for(int i = 0; i < xnumber; ++i){
			//fprintf(rightPartFile, "%28.22g %28.22g %28.22g\n", maxwellEquationRightPart[i][0][0][0], maxwellEquationRightPart[i][0][0][1], maxwellEquationRightPart[i][0][0][2]);
			//	}
			//fclose(rightPartFile);
			for (int i = 0; i < xnumberAdded; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						for (int l = 0; l < maxwellEquationMatrixSize; ++l) {
							gmresOutput[i][j][k][l] = 0;
						}
					}
				}
			}


			bool periodicX = (boundaryConditionTypeX == PERIODIC);
			bool periodicY = (boundaryConditionTypeY == PERIODIC);
			bool periodicZ = (boundaryConditionTypeZ == PERIODIC);
			generalizedMinimalResidualMethod(maxwellEquationMatrix, maxwellEquationRightPart, gmresOutput, xnumberAdded,
			                                 ynumberAdded,
			                                 znumberAdded, maxwellEquationMatrixSize, xnumberGeneral, ynumberGeneral,
			                                 znumberGeneral, additionalBinNumber, maxErrorLevel, maxGMRESIterations, periodicX,
			                                 periodicY, periodicZ, verbosity, leftOutGmresBuffer, rightOutGmresBuffer,
			                                 leftInGmresBuffer, rightInGmresBuffer, frontOutGmresBuffer, backOutGmresBuffer,
			                                 frontInGmresBuffer, backInGmresBuffer, bottomOutGmresBuffer, topOutGmresBuffer,
			                                 bottomInGmresBuffer, topInGmresBuffer, gmresMaxwellBasis, cartComm, cartCoord,
			                                 cartDim);
			//biconjugateStabilizedGradientMethod(maxwellEquationMatrix, maxwellEquationRightPart, gmresOutput, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, maxwellEquationMatrixSize, maxErrorLevel, maxGMRESIterations, boundaryConditionTypeX == PERIODIC, verbosity, cartComm, cartCoord, cartDim, residualBiconjugateMaxwell, firstResidualBiconjugateMaxwell, vBiconjugateMaxwell, pBiconjugateMaxwell, sBiconjugateMaxwell, tBiconjugateMaxwell, leftOutGmresBuffer, rightOutGmresBuffer, leftInGmresBuffer, rightInGmresBuffer, frontOutGmresBuffer, backOutGmresBuffer, frontInGmresBuffer, backInGmresBuffer, bottomOutGmresBuffer, topOutGmresBuffer, bottomInGmresBuffer, topInGmresBuffer);


			MPI_Barrier(cartComm);

			if (periodicX || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
				sendLargeVectorToLeftReceiveFromRight(gmresOutput, leftOutGmresBuffer, rightInGmresBuffer, xnumberAdded,
				                                      ynumberAdded, znumberAdded, 3, additionalBinNumber, cartComm);
			} else if (cartCoord[0] == 0) {
				receiveLargeVectorFromRight(gmresOutput, rightInGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, 3,
				                            additionalBinNumber, cartComm);
			} else if (cartCoord[0] == cartDim[0] - 1) {
				sendLargeVectorToLeft(gmresOutput, leftOutGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, 3,
				                      additionalBinNumber, cartComm);
			}

			if (periodicX || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
				sendLargeVectorToRightReceiveFromLeft(gmresOutput, rightOutGmresBuffer, leftInGmresBuffer, xnumberAdded,
				                                      ynumberAdded, znumberAdded, 3, additionalBinNumber, cartComm);
			} else if (cartCoord[0] == 0) {
				sendLargeVectorToRight(gmresOutput, rightOutGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, 3,
				                       additionalBinNumber, cartComm);
			} else if (cartCoord[0] == cartDim[0] - 1) {
				receiveLargeVectorFromLeft(gmresOutput, leftInGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, 3,
				                           additionalBinNumber, cartComm);
			}

			sendLargeVectorToFrontReceiveFromBack(gmresOutput, frontOutGmresBuffer, backInGmresBuffer, xnumberAdded,
			                                      ynumberAdded, znumberAdded, 3, additionalBinNumber, cartComm);
			sendLargeVectorToBackReceiveFromFront(gmresOutput, backOutGmresBuffer, frontInGmresBuffer, xnumberAdded,
			                                      ynumberAdded, znumberAdded, 3, additionalBinNumber, cartComm);

			sendLargeVectorToBottomReceiveFromTop(gmresOutput, bottomOutGmresBuffer, topInGmresBuffer, xnumberAdded,
			                                      ynumberAdded, znumberAdded, 3, additionalBinNumber, cartComm);
			sendLargeVectorToTopReceiveFromBottom(gmresOutput, topOutGmresBuffer, bottomInGmresBuffer, xnumberAdded,
			                                      ynumberAdded, znumberAdded, 3, additionalBinNumber, cartComm);

			MPI_Barrier(cartComm);

			double omega2 = 0;
			for (int i = 0; i < typesNumber; ++i) {
				omega2 += 4 * pi * types[i].concentration * types[i].charge * types[i].charge / types[i].mass;
			}

			double gamma = 1.0 / sqrt(1 - V0.scalarMult(V0) / speed_of_light_normalized_sqr);
			double omega = sqrt(omega2 / (gamma * gamma * gamma));
			double a = cos(omega * (time + theta * deltaT));
			for (int i = 0; i < xnumberAdded; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						for (int l = 0; l < 3; ++l) {
							tempEfield[i][j][k][l] = gmresOutput[i][j][k][l];
						}
						//tempEfield[i][j][k].x = 0;
						//tempEfield[i][j][k].y = 0;
						//tempEfield[i][j][k].z = 0;
					}
				}
			}

			//updateBoundaries();

			double alfvenV = B0.norm() / sqrt(4 * pi * density);
			//double k = 2 * pi / xsize;

			evaluateExplicitDerivative();
			//smoothEderivative();
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						explicitEfield[i][j][k] += Ederivative[i][j][k] * deltaT;
					}
				}
			}

			//evaluateMagneticField();
			/*if (boundaryConditionTypeX == PERIODIC) {
			    for (int j = 0; j < ynumber; ++j) {
			        for (int k = 0; k < znumber; ++k) {
			            tempEfield[xnumber][j][k] = tempEfield[0][j][k];
			            explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
			        }
			    }
			} else {
			    for (int j = 0; j < ynumber; ++j) {
			        for (int k = 0; k < znumber; ++k) {
			            tempEfield[xnumber][j][k] = tempEfield[xnumber - 1][j][k];
	
			            explicitEfield[xnumber][j][k] = tempEfield[xnumber][j][k];
			        }
			    }
			}*/

			/*for (int i = 0; i < xnumber + 1; ++i) {
				for (int j = 0; j < ynumber + 1; ++j) {
					tempEfield[i][j][znumber] = tempEfield[i][j][0];
					explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
				}
			}
	
			for (int i = 0; i < xnumber + 1; ++i) {
				for (int k = 0; k < znumber + 1; ++k) {
					tempEfield[i][ynumber][k] = tempEfield[i][0][k];
					explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
				}
			}*/

			if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						tempEfield[xnumberAdded][j][k] = E0;
					}
				}
			}

			//smoothTempEfield();

			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						newEfield[i][j][k] = (tempEfield[i][j][k] - Efield[i][j][k] * (1 - theta)) / theta;
						//newEfield[i][j][k].x= 0;
						//newEfield[i][j][k].x = E0.x*cos(omega*(time + deltaT));
						//newEfield[i][j][k].y= 0;
						//newEfield[i][j][k].z= 0;
						if (newEfield[i][j][k].x > Efield[i][j][k].x) {
							//printf("aaa\n");
						}
					}
				}
			}

			if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						newEfield[xnumberAdded][j][k] = E0;
					}
				}
			}
		}

		if (solverType == EXPLICIT) {
			evaluateExplicitDerivative();
			//smoothEderivative();
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						explicitEfield[i][j][k] += Ederivative[i][j][k] * deltaT;
					}
				}
			}

			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						newEfield[i][j][k] = explicitEfield[i][j][k];
						tempEfield[i][j][k] = Efield[i][j][k] * (1 - theta) + newEfield[i][j][k] * theta;
					}
				}
			}
		}
	}
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("evaluating electric field time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
	//fclose(outputEverythingFile);
}

void Simulation::evaluateExplicitDerivative() {
	for (int i = 1; i < xnumberAdded; ++i) {
		for (int j = 1; j < ynumberAdded; ++j) {
			for (int k = 1; k < znumberAdded; ++k) {
				rotB[i][j][k] = evaluateRotB(i, j, k) * speed_of_light_normalized;
				Ederivative[i][j][k] = rotB[i][j][k] -
					(electricFlux[i][j][k] * 4 * pi);
				if (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
					if (i <= 1 + additionalBinNumber) {
						Ederivative[i][j][k].y = 0;
						Ederivative[i][j][k].z = 0;
					}
					if (i == xnumberAdded) {
						Ederivative[i][j][k] = Vector3d(0, 0, 0);
					}
				}
			}
		}
	}
}

void Simulation::updateEfield() {
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				Efield[i][j][k] = newEfield[i][j][k];
				//Efield[i].y = newEfield[i].y;
				//Efield[i].z = newEfield[i].z;
			}
		}
	}
}

void Simulation::updateBfield() {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k] = newBfield[i][j][k];
			}
		}
	}
}

void Simulation::updateFields() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if (solverType == BUNEMAN) {
		updateBunemanElectricField();
		updateBunemanMagneticField();
	} else {
		updateEfield();
		updateBfield();
	}
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating fields time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::evaluateMaxwellEquationMatrix() {
	int resistiveLayerWidth = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < maxwellEquationMatrixSize; ++l) {
					maxwellEquationMatrix[i][j][k][l].clear();
				}
			}
		}
	}

	if (cartDim[0] > 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					if ((j <= additionalBinNumber) || (j >= ynumberAdded - 1 - additionalBinNumber) || (k <= additionalBinNumber) || (k
						>= znumberAdded - 1 - additionalBinNumber)) {
						createFakeEquation(i, j, k);
					} else {
						if (cartCoord[0] == 0) {
							if (i <= additionalBinNumber) {
								if (boundaryConditionTypeX == FREE_BOTH) {
									createFreeLeftEquation(i, j, k);
								} else {
									createFakeEquation(i, j, k);
								}
							} else if (i == 1 + additionalBinNumber) {
								if (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
									createSuperConductorLeftEquation(i, j, k);
								} else {
									createInternalEquation(i, j, k);
								}
							} else if (i < xnumberAdded - 1) {
								createInternalEquation(i, j, k);
							} else {
								createFakeEquation(i, j, k);
							}
						} else if (cartCoord[0] == cartDim[0] - 1) {
							if (i <= additionalBinNumber) {
								createFakeEquation(i, j, k);
								//} else if (i < xnumberAdded - 1 - additionalBinNumber - resistiveLayerWidth){
								//	createInternalEquation(i, j, k);
							} else if (i < xnumberAdded - 1 - additionalBinNumber) {
								createInternalEquation(i, j, k);
								//createResistiveEquation(i, j, k);
							} else if (i < xnumberAdded - 1) {
								if (boundaryConditionTypeX != PERIODIC) {
									createFreeRightEquation(i, j, k);
								} else {
									createInternalEquation(i, j, k);
								}
							} else {
								if (boundaryConditionTypeX != PERIODIC) {
									createFreeRightEquation(i, j, k);
								} else {
									createFakeEquation(i, j, k);
								}
							}
						} else {
							if ((i > additionalBinNumber && i < xnumberAdded - additionalBinNumber - 1)) {
								createInternalEquation(i, j, k);
							} else {
								createFakeEquation(i, j, k);
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
						createFakeEquation(i, j, k);
					} else {
						if (i <= additionalBinNumber) {
							if (boundaryConditionTypeX == FREE_BOTH) {
								createFreeLeftEquation(i, j, k);
							} else {
								createFakeEquation(i, j, k);
							}
						} else if (i == 1 + additionalBinNumber) {
							if (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
								createSuperConductorLeftEquation(i, j, k);
							} else {
								createInternalEquation(i, j, k);
							}
							//} else if (i < xnumberAdded - 1 - additionalBinNumber - resistiveLayerWidth) {
							//	createInternalEquation(i, j, k);
						} else if (i < xnumberAdded - 1 - additionalBinNumber) {
							createInternalEquation(i, j, k);
							//createResistiveEquation(i, j, k);
						} else if (i < xnumberAdded - 1) {
							if (boundaryConditionTypeX != PERIODIC) {
								createFreeRightEquation(i, j, k);
							} else {
								createInternalEquation(i, j, k);
							}
						} else {
							if (boundaryConditionTypeX != PERIODIC) {
								createFreeRightEquation(i, j, k);
							} else {
								createFakeEquation(i, j, k);
							}
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				if (fabs(maxwellEquationRightPart[i][j][k][0]) > 1E100) {
					printf("too large maxweell rightPart x\n");
					MPI_Finalize();
					exit(0);
				}
				if (fabs(maxwellEquationRightPart[i][j][k][1]) > 1E100) {
					printf("too large maxweell rightPart y\n");
					MPI_Finalize();
					exit(0);
				}
				if (fabs(maxwellEquationRightPart[i][j][k][2]) > 1E100) {
					printf("too large maxweell rightPart z\n");
					MPI_Finalize();
					exit(0);
				}
			}
		}
	}

	if (debugMode) {
		checkEquationMatrix(maxwellEquationMatrix, maxwellEquationMatrixSize);
	}
}

void Simulation::checkEquationMatrix(std::vector < MatrixElement >**** matrix, int lnumber) {
	//#pragma omp parallel for
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];
						if (element.i < 0) {
							if (rank == 0) printf("element i < 0\n");
							fflush(stdout);
							errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
							fprintf(errorLogFile, "element i = %d < 0\n", element.i);
							fclose(errorLogFile);
							MPI_Finalize();
							exit(0);
						}
						if (element.i > xnumberAdded - 1) {
							if (rank == 0) printf("element i > xnumber");
							fflush(stdout);
							errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
							fprintf(errorLogFile, "element i = %d > xnumber = %d\n", element.i, xnumberAdded - 1);
							fclose(errorLogFile);
							MPI_Finalize();
							exit(0);
						}

						if (element.j < 0) {
							if (rank == 0) printf("element j < 0\n");
							fflush(stdout);
							errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
							fprintf(errorLogFile, "element j = %d < 0\n", element.j);
							fclose(errorLogFile);
							MPI_Finalize();
							exit(0);
						}
						if (element.j >= ynumberAdded) {
							if (rank == 0) printf("element j >= ynumber");
							fflush(stdout);
							errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
							fprintf(errorLogFile, "element j = %d >= ynumber = %d\n", element.j, ynumberAdded);
							fclose(errorLogFile);
							MPI_Finalize();
							exit(0);
						}

						if (element.k < 0) {
							if (rank == 0) printf("element k < 0\n");
							fflush(stdout);
							errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
							fprintf(errorLogFile, "element k = %d < 0\n", element.k);
							fclose(errorLogFile);
							MPI_Finalize();
							exit(0);
						}
						if (element.k >= znumberAdded) {
							if (rank == 0) printf("element k >= znumber");
							fflush(stdout);
							errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
							fprintf(errorLogFile, "element k = %d >= xnumber = %d\n", element.k, znumber);
							fclose(errorLogFile);
							MPI_Finalize();
							exit(0);
						}
						for (int n = m + 1; n < matrix[i][j][k][l].size(); ++n) {
							MatrixElement tempElement = matrix[i][j][k][l][n];

							if (element.equalsIndex(tempElement)) {
								printf("equals indexes\n");
								printf("current = %d %d %d %d\n", i, j, k, l);
								printf("temp = %d %d %d %d\n", tempElement.i, tempElement.j, tempElement.k,
								       tempElement.l);
								fflush(stdout);
								errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
								fprintf(errorLogFile, "equal indexes current = %d %d %d %d temp = %d %d %d %d\n", i, j,
								        k, l, element.i, element.j, element.k, element.l);
								fclose(errorLogFile);
								MPI_Finalize();
								exit(0);
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::createSuperConductorLeftEquation(int i, int j, int k) {
	//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
	//maxwellEquationRightPart[i][j][k][0] = 0;
	//maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(1.0, i, j, k, 1));
	//maxwellEquationRightPart[i][j][k][1] = 0;
	//maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(1.0, i, j, k, 2));
	//maxwellEquationRightPart[i][j][k][2] = 0;
	//return;

	int nextJ = j + 1;

	int nextK = k + 1;


	if ((ynumberGeneral) > 1 && (znumberGeneral > 1)) {
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-1.0, i + 1, j, k, 0));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaX, i, j, k, 0));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaX, i + 1, j, k, 0));

		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaX, i, nextJ, k, 0));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaX, i + 1, nextJ, k, 0));

		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaX, i, j, nextK, 0));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaX, i + 1, j, nextK, 0));

		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaX, i, nextJ, nextK, 0));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaX, i + 1, nextJ, nextK, 0));

		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaY, i + 1, j, k, 1));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaY, i + 1, nextJ, k, 1));

		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaY, i + 1, j, nextK, 1));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaY, i + 1, nextJ, nextK, 1));

		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaZ, i + 1, j, k, 2));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaZ, i + 1, j, nextK, 2));

		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaZ, i + 1, nextJ, k, 2));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaZ, i + 1, nextJ, nextK, 2));
	} else if (ynumberGeneral > 1) {
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-1.0, i + 1, j, k, 0));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaX, i, j, k, 0));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaX, i + 1, j, k, 0));

		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaX, i, nextJ, k, 0));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaX, i + 1, nextJ, k, 0));

		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaY, i + 1, j, k, 1));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaY, i + 1, nextJ, k, 1));
	} else if (znumberGeneral > 1) {
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-1.0, i + 1, j, k, 0));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaX, i, j, k, 0));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaX, i + 1, j, k, 0));

		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaX, i, j, nextK, 0));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaX, i + 1, j, nextK, 0));

		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaZ, i + 1, j, k, 2));
		//maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaZ, i + 1, j, nextK, 2));
	} else {
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0 / deltaX, i, j, k, 0));
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-1.0 / deltaX, i + 1, j, k, 0));
	}

	maxwellEquationRightPart[i][j][k][0] = -pi * (chargeDensityHat[i][j][k] + chargeDensityHat[i][j - 1][k] +
		chargeDensityHat[i][j][k - 1] + chargeDensityHat[i][j - 1][k - 1]);
	//maxwellEquationRightPart[i][j][k][0] = -4 * pi * chargeDensity[i][j][k];

	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(1.0, i, j, k, 1));
	maxwellEquationRightPart[i][j][k][1] = 0;
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(1.0, i, j, k, 2));
	maxwellEquationRightPart[i][j][k][2] = 0;
}

void Simulation::createFreeRightEquation(int i, int j, int k) {
	Vector3d rightPart = E0;

	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(1.0, i, j, k, 1));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(1.0, i, j, k, 2));

	alertNaNOrInfinity(rightPart.x, "right part x = NaN in create free right");
	alertNaNOrInfinity(rightPart.y, "right part y = NaN");
	alertNaNOrInfinity(rightPart.z, "right part z = NaN");

	maxwellEquationRightPart[i][j][k][0] = rightPart.x;
	maxwellEquationRightPart[i][j][k][1] = rightPart.y;
	maxwellEquationRightPart[i][j][k][2] = rightPart.z;
}

void Simulation::createFreeLeftEquation(int i, int j, int k) {
	Vector3d rightPart = Efield[i][j][k];

	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(1.0, i, j, k, 1));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(1.0, i, j, k, 2));

	alertNaNOrInfinity(rightPart.x, "right part x = NaN in create free right");
	alertNaNOrInfinity(rightPart.y, "right part y = NaN");
	alertNaNOrInfinity(rightPart.z, "right part z = NaN");

	maxwellEquationRightPart[i][j][k][0] = rightPart.x;
	maxwellEquationRightPart[i][j][k][1] = rightPart.y;
	maxwellEquationRightPart[i][j][k][2] = rightPart.z;
}

void Simulation::createInternalECEquationX(int i, int j, int k) {
	double element = 1.0 + massMatrix[i][j][k].matrix[splineOrder + 1][splineOrder + 1][splineOrder + 1].matrix[0][0];
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

	element = massMatrix[i][j][k].matrix[splineOrder + 1][splineOrder + 1][splineOrder + 1].matrix[0][1];
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 1));

	element = massMatrix[i][j][k].matrix[splineOrder + 1][splineOrder + 1][splineOrder + 1].matrix[0][2];
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 2));

	for (int tempI = 0; tempI < 2 * splineOrder + 3; ++ tempI) {
		for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
			for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
				if (tempI != splineOrder + 1 && tempJ != splineOrder + 1 && tempK != splineOrder + 1) {
					int xindex = massMatrix[i][j][k].xindex[tempI];
					int yindex = massMatrix[i][j][k].yindex[tempJ];
					int zindex = massMatrix[i][j][k].zindex[tempK];
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[0][0];
					maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, xindex, yindex, zindex, 0));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[0][1];
					maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, xindex, yindex, zindex, 1));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[0][2];
					maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, xindex, yindex, zindex, 2));
				}
			}
		}
	}
}

void Simulation::createInternalECEquationY(int i, int j, int k) {
	double element = 1.0 + massMatrix[i][j][k].matrix[splineOrder + 1][splineOrder + 1][splineOrder + 1].matrix[1][1];
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 1));

	element = massMatrix[i][j][k].matrix[splineOrder + 1][splineOrder + 1][splineOrder + 1].matrix[1][0];
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 0));

	element = massMatrix[i][j][k].matrix[splineOrder + 1][splineOrder + 1][splineOrder + 1].matrix[1][2];
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 2));

	for (int tempI = 0; tempI < 2 * splineOrder + 3; ++ tempI) {
		for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
			for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
				if (tempI != splineOrder + 1 && tempJ != splineOrder + 1 && tempK != splineOrder + 1) {
					int xindex = massMatrix[i][j][k].xindex[tempI];
					int yindex = massMatrix[i][j][k].yindex[tempJ];
					int zindex = massMatrix[i][j][k].zindex[tempK];
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[1][0];
					maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, xindex, yindex, zindex, 0));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[1][1];
					maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, xindex, yindex, zindex, 1));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[1][2];
					maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, xindex, yindex, zindex, 2));
				}
			}
		}
	}
}

void Simulation::createInternalECEquationZ(int i, int j, int k) {
	double element = 1.0 + massMatrix[i][j][k].matrix[splineOrder + 1][splineOrder + 1][splineOrder + 1].matrix[2][2];
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 2));

	element = massMatrix[i][j][k].matrix[splineOrder + 1][splineOrder + 1][splineOrder + 1].matrix[2][1];
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 1));

	element = massMatrix[i][j][k].matrix[splineOrder + 1][splineOrder + 1][splineOrder + 1].matrix[2][0];
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 0));

	for (int tempI = 0; tempI < 2 * splineOrder + 3; ++ tempI) {
		for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
			for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
				if (tempI != splineOrder + 1 && tempJ != splineOrder + 1 && tempK != splineOrder + 1) {
					int xindex = massMatrix[i][j][k].xindex[tempI];
					int yindex = massMatrix[i][j][k].yindex[tempJ];
					int zindex = massMatrix[i][j][k].zindex[tempK];
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[2][0];
					maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, xindex, yindex, zindex, 0));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[2][1];
					maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, xindex, yindex, zindex, 1));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[2][2];
					maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, xindex, yindex, zindex, 2));
				}
			}
		}
	}
}

void Simulation::createInternalEquationX(int i, int j, int k) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT) * speed_of_light_correction_sqr;
	//c_theta_deltaT2 = 0;
	double element = 1.0 - dielectricTensor[i][j][k].matrix[0][0];
	if (isInResistiveLayer(i, j, k)) {
		element += fakeCondactivity * theta * deltaT;
	}

	bool _27_points = false;
	//bool _27_points = true;
	int nextI = i + 1;
	int prevI = i - 1;
	int nextJ = j + 1;
	int prevJ = j - 1;
	int nextK = k + 1;
	int prevK = k - 1;

	if (_27_points) {
		if (ynumberGeneral > 1 && znumberGeneral > 1) {
			element += -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][k].matrix[0][0]) / (2 * deltaX2) - 1.0 / (2 * deltaY2)
				- 1.0 / (2 * deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][k].matrix[0][0]) / (4 * deltaX2) - 1.0 / (4 * deltaY2
			) - 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][k].matrix[0][0]) / (4 * deltaX2) - 1.0 / (4 * deltaY2
			) - 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

			element = -c_theta_deltaT2 * (- (1.0 - dielectricTensor[i][nextJ][k].matrix[0][0]) / (4 * deltaX2) + 1.0 / (4 *
				deltaY2) - 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
			element = -c_theta_deltaT2 * (- (1.0 - dielectricTensor[i][prevJ][k].matrix[0][0]) / (4 * deltaX2) + 1.0 / (4 *
				deltaY2) - 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][nextK].matrix[0][0]) / (4 * deltaX2) - 1.0 / (4 *
				deltaY2) + 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][prevK].matrix[0][0]) / (4 * deltaX2) - 1.0 / (4 *
				deltaY2) + 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));

			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][nextJ][nextK].matrix[0][0]) / (8 * deltaX2) + 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, nextK, 0));
			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][nextJ][prevK].matrix[0][0]) / (8 * deltaX2) + 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, prevK, 0));
			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][prevJ][nextK].matrix[0][0]) / (8 * deltaX2) + 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, nextK, 0));
			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][prevJ][prevK].matrix[0][0]) / (8 * deltaX2) + 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, prevK, 0));

			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][nextK].matrix[0][0]) / (8 * deltaX2) - 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[2][0] / (8 * deltaX *
				deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 0));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][prevK].matrix[0][0]) / (8 * deltaX2) - 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[2][0] / (8 * deltaX *
				deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 0));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][nextK].matrix[0][0]) / (8 * deltaX2) - 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[2][0] / (8 * deltaX *
				deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 0));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][prevK].matrix[0][0]) / (8 * deltaX2) - 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[2][0] / (8 * deltaX *
				deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 0));

			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][nextJ][k].matrix[0][0]) / (8 * deltaX2) + 1.0 / (8 *
				deltaY2) - 1.0 / (8 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[1][0] / (8 * deltaX *
				deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 0));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][prevJ][k].matrix[0][0]) / (8 * deltaX2) + 1.0 / (8 *
				deltaY2) - 1.0 / (8 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[1][0] / (8 * deltaX *
				deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 0));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][prevJ][k].matrix[0][0]) / (8 * deltaX2) + 1.0 / (8 *
				deltaY2) - 1.0 / (8 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[1][0] / (8 * deltaX *
				deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 0));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][nextJ][k].matrix[0][0]) / (8 * deltaX2) + 1.0 / (8 *
				deltaY2) - 1.0 / (8 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[1][0] / (8 * deltaX *
				deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 0));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][nextJ][nextK].matrix[0][0]) / (16 * deltaX2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][nextK].matrix[1][0] / (16 * deltaX
				* deltaY) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][nextK].matrix[2][0] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, nextK, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][nextJ][prevK].matrix[0][0]) / (16 * deltaX2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][prevK].matrix[1][0] / (16 * deltaX
				* deltaY) - c_theta_deltaT2 * dielectricTensor[nextI][nextJ][prevK].matrix[2][0] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, prevK, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][prevJ][nextK].matrix[0][0]) / (16 * deltaX2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][nextK].matrix[1][0] / (16 * deltaX
				* deltaY) + c_theta_deltaT2 * dielectricTensor[nextI][prevJ][nextK].matrix[2][0] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, nextK, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][prevJ][prevK].matrix[0][0]) / (16 * deltaX2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][prevK].matrix[1][0] / (16 * deltaX
				* deltaY) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][prevK].matrix[2][0] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, prevK, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][nextJ][nextK].matrix[0][0]) / (16 * deltaX2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][nextK].matrix[1][0] / (16 * deltaX
				* deltaY) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][nextK].matrix[2][0] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, nextK, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][nextJ][prevK].matrix[0][0]) / (16 * deltaX2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][prevK].matrix[1][0] / (16 * deltaX
				* deltaY) + c_theta_deltaT2 * dielectricTensor[prevI][nextJ][prevK].matrix[2][0] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, prevK, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][prevJ][nextK].matrix[0][0]) / (16 * deltaX2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][nextK].matrix[1][0] / (16 * deltaX
				* deltaY) - c_theta_deltaT2 * dielectricTensor[prevI][prevJ][nextK].matrix[2][0] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, nextK, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][prevJ][prevK].matrix[0][0]) / (16 * deltaX2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][prevK].matrix[1][0] / (16 * deltaX
				* deltaY) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][prevK].matrix[2][0] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, prevK, 0));

			///y
			element = -dielectricTensor[i][j][k].matrix[0][1] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[0][1]) / (2 *
				deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 1));


			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[0][1]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[0][1]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 1));

			element = -c_theta_deltaT2 * (dielectricTensor[i][nextJ][k].matrix[0][1]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 1));
			element = -c_theta_deltaT2 * (dielectricTensor[i][prevJ][k].matrix[0][1]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 1));

			element = -c_theta_deltaT2 * (dielectricTensor[i][j][nextK].matrix[0][1]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 1));
			element = -c_theta_deltaT2 * (dielectricTensor[i][j][prevK].matrix[0][1]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 1));

			element = - c_theta_deltaT2 * (dielectricTensor[i][nextJ][nextK].matrix[0][1]) / (8 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, nextK, 1));
			element = - c_theta_deltaT2 * (dielectricTensor[i][nextJ][prevK].matrix[0][1]) / (8 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, prevK, 1));
			element = - c_theta_deltaT2 * (dielectricTensor[i][prevJ][nextK].matrix[0][1]) / (8 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, nextK, 1));
			element = - c_theta_deltaT2 * (dielectricTensor[i][prevJ][prevK].matrix[0][1]) / (8 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, prevK, 1));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][nextK].matrix[0][1]) / (8 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][j][nextK].matrix[2][1] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][prevK].matrix[0][1]) / (8 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][j][prevK].matrix[2][1] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][nextK].matrix[0][1]) / (8 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][j][nextK].matrix[2][1] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][prevK].matrix[0][1]) / (8 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][j][prevK].matrix[2][1] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 1));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][k].matrix[0][1]) / (8 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][k].matrix[1][1] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][k].matrix[0][1]) / (8 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][k].matrix[1][1] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][k].matrix[0][1]) / (8 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][k].matrix[1][1] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][k].matrix[0][1]) / (8 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][k].matrix[1][1] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 1));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][nextK].matrix[0][1]) / (16 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][nextK].matrix[1][1] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][nextK].matrix[2][1] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][prevK].matrix[0][1]) / (16 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][prevK].matrix[1][1] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][prevK].matrix[2][1] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, prevK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][nextK].matrix[0][1]) / (16 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][nextK].matrix[1][1] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][nextK].matrix[2][1] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][prevK].matrix[0][1]) / (16 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][prevK].matrix[1][1] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][prevK].matrix[2][1] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, prevK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][nextK].matrix[0][1]) / (16 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][nextK].matrix[1][1] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][nextK].matrix[2][1] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][prevK].matrix[0][1]) / (16 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][prevK].matrix[1][1] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][prevK].matrix[2][1] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, prevK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][nextK].matrix[0][1]) / (16 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][nextK].matrix[1][1] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][nextK].matrix[2][1] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][prevK].matrix[0][1]) / (16 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][prevK].matrix[1][1] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][prevK].matrix[2][1] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, prevK, 1));

			///z
			element = - dielectricTensor[i][j][k].matrix[0][2] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[0][2]) / (2
				* deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 2));


			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[0][2]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[0][2]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 2));

			element = -c_theta_deltaT2 * (dielectricTensor[i][nextJ][k].matrix[0][2]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 2));
			element = -c_theta_deltaT2 * (dielectricTensor[i][prevJ][k].matrix[0][2]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 2));

			element = -c_theta_deltaT2 * (dielectricTensor[i][j][nextK].matrix[0][2]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 2));
			element = -c_theta_deltaT2 * (dielectricTensor[i][j][prevK].matrix[0][2]) / (4 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 2));

			element = - c_theta_deltaT2 * (dielectricTensor[i][nextJ][nextK].matrix[0][2]) / (8 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, nextK, 2));
			element = - c_theta_deltaT2 * (dielectricTensor[i][nextJ][prevK].matrix[0][2]) / (8 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, prevK, 2));
			element = - c_theta_deltaT2 * (dielectricTensor[i][prevJ][nextK].matrix[0][2]) / (8 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, nextK, 2));
			element = - c_theta_deltaT2 * (dielectricTensor[i][prevJ][prevK].matrix[0][2]) / (8 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, prevK, 2));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][nextK].matrix[0][2]) / (8 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][j][nextK].matrix[2][2] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][prevK].matrix[0][2]) / (8 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][j][prevK].matrix[2][2] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][nextK].matrix[0][2]) / (8 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][j][nextK].matrix[2][2] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][prevK].matrix[0][2]) / (8 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][j][prevK].matrix[2][2] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 2));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][k].matrix[0][2]) / (8 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][k].matrix[1][2] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][k].matrix[0][2]) / (8 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][k].matrix[1][2] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][k].matrix[0][2]) / (8 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][k].matrix[1][2] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][k].matrix[0][2]) / (8 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][k].matrix[1][2] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 2));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][nextK].matrix[0][2]) / (16 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][nextK].matrix[1][2] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][nextK].matrix[2][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][prevK].matrix[0][2]) / (16 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][prevK].matrix[1][2] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][prevK].matrix[2][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, prevK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][nextK].matrix[0][2]) / (16 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][nextK].matrix[1][2] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][nextK].matrix[2][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][prevK].matrix[0][2]) / (16 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][prevK].matrix[1][2] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][prevK].matrix[2][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, prevK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][nextK].matrix[0][2]) / (16 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][nextK].matrix[1][2] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][nextK].matrix[2][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][prevK].matrix[0][2]) / (16 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][prevK].matrix[1][2] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][prevK].matrix[2][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, prevK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][nextK].matrix[0][2]) / (16 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][nextK].matrix[1][2] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][nextK].matrix[2][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][prevK].matrix[0][2]) / (16 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][prevK].matrix[1][2] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][prevK].matrix[2][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, prevK, 2));
		} else if (ynumberGeneral > 1) {
			element += -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][k].matrix[0][0]) / (deltaX2) - 1.0 / (deltaY2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][k].matrix[0][0]) / (2 * deltaX2) - 1.0 / (2 * deltaY2
			));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][k].matrix[0][0]) / (2 * deltaX2) - 1.0 / (2 * deltaY2
			));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][nextJ][k].matrix[0][0]) / (2 * deltaX2) + 1.0 / (2 *
				deltaY2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][prevJ][k].matrix[0][0]) / (2 * deltaX2) + 1.0 / (2 *
				deltaY2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][nextJ][k].matrix[0][0]) / (4 * deltaX2) + 1.0 / (4 *
				deltaY2)) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[1][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][prevJ][k].matrix[0][0]) / (4 * deltaX2) + 1.0 / (4 *
				deltaY2)) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[1][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][nextJ][k].matrix[0][0]) / (4 * deltaX2) + 1.0 / (4 *
				deltaY2)) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[1][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][prevJ][k].matrix[0][0]) / (4 * deltaX2) + 1.0 / (4 *
				deltaY2)) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[1][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 0));

			element = -dielectricTensor[i][j][k].matrix[0][1] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[0][1]) / (
				deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[0][1]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[0][1]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 1));

			element = -c_theta_deltaT2 * (dielectricTensor[i][nextJ][k].matrix[0][1]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 1));
			element = -c_theta_deltaT2 * (dielectricTensor[i][prevJ][k].matrix[0][1]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 1));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][k].matrix[0][1]) / (4 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][k].matrix[1][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][k].matrix[0][1]) / (4 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][k].matrix[1][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][k].matrix[0][1]) / (4 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][k].matrix[1][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][k].matrix[0][1]) / (4 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][k].matrix[1][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 1));

			element = -dielectricTensor[i][j][k].matrix[0][2] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[0][2]) / (
				deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[0][2]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[0][2]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 2));

			element = -c_theta_deltaT2 * (dielectricTensor[i][nextJ][k].matrix[0][2]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 2));
			element = -c_theta_deltaT2 * (dielectricTensor[i][prevJ][k].matrix[0][2]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 2));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][k].matrix[0][2]) / (4 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][k].matrix[1][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][k].matrix[0][2]) / (4 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][k].matrix[1][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][k].matrix[0][2]) / (4 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][k].matrix[1][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][k].matrix[0][2]) / (4 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][k].matrix[1][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 2));
		} else if (znumberGeneral > 1) {
			element += -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][k].matrix[0][0]) / (deltaX2) - 1.0 / (deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][k].matrix[0][0]) / (2 * deltaX2) - 1.0 / (2 * deltaZ2
			));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][k].matrix[0][0]) / (2 * deltaX2) - 1.0 / (2 * deltaZ2
			));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][nextK].matrix[0][0]) / (2 * deltaX2) + 1.0 / (2 *
				deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][prevK].matrix[0][0]) / (2 * deltaX2) + 1.0 / (2 *
				deltaZ2));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][nextK].matrix[0][0]) / (4 * deltaX2) + 1.0 / (4 *
				deltaZ2)) + c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[2][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][prevK].matrix[0][0]) / (4 * deltaX2) + 1.0 / (4 *
				deltaZ2)) - c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[2][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][nextK].matrix[0][0]) / (4 * deltaX2) + 1.0 / (4 *
				deltaZ2)) - c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[2][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 0));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][prevK].matrix[0][0]) / (4 * deltaX2) + 1.0 / (4 *
				deltaZ2)) + c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[2][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 0));

			element = -dielectricTensor[i][j][k].matrix[0][1] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[0][1]) / (
				deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[0][1]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[0][1]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 1));

			element = -c_theta_deltaT2 * (dielectricTensor[i][j][nextK].matrix[0][1]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 1));
			element = -c_theta_deltaT2 * (dielectricTensor[i][j][prevK].matrix[0][1]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 1));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][nextK].matrix[0][1]) / (4 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][j][nextK].matrix[2][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][prevK].matrix[0][1]) / (4 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][j][prevK].matrix[2][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][nextK].matrix[0][1]) / (4 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][j][nextK].matrix[2][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][prevK].matrix[0][1]) / (4 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][j][prevK].matrix[2][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 1));

			element = -dielectricTensor[i][j][k].matrix[0][2] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[0][2]) / (
				deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[0][2]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[0][2]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 2));

			element = -c_theta_deltaT2 * (dielectricTensor[i][j][nextK].matrix[0][2]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 2));
			element = -c_theta_deltaT2 * (dielectricTensor[i][j][prevK].matrix[0][2]) / (2 * deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 2));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][nextK].matrix[0][2]) / (4 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[nextI][j][nextK].matrix[1][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][prevK].matrix[0][2]) / (4 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[nextI][j][prevK].matrix[1][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][nextK].matrix[0][2]) / (4 * deltaX2) - c_theta_deltaT2 *
				dielectricTensor[prevI][j][nextK].matrix[1][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][prevK].matrix[0][2]) / (4 * deltaX2) + c_theta_deltaT2 *
				dielectricTensor[prevI][j][prevK].matrix[1][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 2));
		} else {
			element += c_theta_deltaT2 * ((2.0 - 2 * dielectricTensor[i][j][k].matrix[0][0]) / deltaX2);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

			element = -dielectricTensor[i][j][k].matrix[0][1];
			element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[0][1] / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 1));

			element = -dielectricTensor[i][j][k].matrix[0][2];
			element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[0][2] / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 2));

			element = c_theta_deltaT2 * (-1.0 + dielectricTensor[nextI][j][k].matrix[0][0]) / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
			element = c_theta_deltaT2 * dielectricTensor[nextI][j][k].matrix[0][1] / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 1));
			element = c_theta_deltaT2 * dielectricTensor[nextI][j][k].matrix[0][2] / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 2));

			element = c_theta_deltaT2 * (-1.0 + dielectricTensor[prevI][j][k].matrix[0][0]) / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));
			element = c_theta_deltaT2 * dielectricTensor[prevI][j][k].matrix[0][1] / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 1));
			element = c_theta_deltaT2 * dielectricTensor[prevI][j][k].matrix[0][2] / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 2));
		}

		///
		/*for(int tempI = 0; tempI < 2*splineOrder + 3; ++ tempI){
		for(int tempJ = 0; tempJ < 2*splineOrder + 3; ++tempJ){
			for(int tempK = 0; tempK < 2*splineOrder + 3; ++tempK){
					int xindex = massMatrix[i][j][k].xindex[tempI];
					int yindex = massMatrix[i][j][k].yindex[tempJ];
					int zindex = massMatrix[i][j][k].zindex[tempK];
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[0][0];
					maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, xindex, yindex, zindex, 0));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[0][1];
					maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, xindex, yindex, zindex, 1));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[0][2];
					maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, xindex, yindex, zindex, 2));
			}
			}
		}*/
		///
	} else {
		if (xnumberGeneral > 1) {
			element += c_theta_deltaT2 * ((2.0 - 2 * dielectricTensor[i][j][k].matrix[0][0]) / deltaX2);
		}
		if (ynumberGeneral > 1) {
			element += c_theta_deltaT2 * 2.0 / deltaY2;
		}
		if (znumberGeneral > 1) {
			element += c_theta_deltaT2 * 2.0 / deltaZ2;
		}

		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = -dielectricTensor[i][j][k].matrix[0][1];
		if (xnumberGeneral > 1) {
			element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[0][1] / deltaX2;
		}
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 1));

		element = -dielectricTensor[i][j][k].matrix[0][2];
		if (xnumberGeneral > 1) {
			element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[0][2] / deltaX2;
		}
		maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 2));


		if (xnumberGeneral > 1) {
			element = c_theta_deltaT2 * (-1.0 + dielectricTensor[nextI][j][k].matrix[0][0]) / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
			element = c_theta_deltaT2 * dielectricTensor[nextI][j][k].matrix[0][1] / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 1));
			element = c_theta_deltaT2 * dielectricTensor[nextI][j][k].matrix[0][2] / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 2));

			element = c_theta_deltaT2 * (-1.0 + dielectricTensor[prevI][j][k].matrix[0][0]) / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));
			element = c_theta_deltaT2 * dielectricTensor[prevI][j][k].matrix[0][1] / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 1));
			element = c_theta_deltaT2 * dielectricTensor[prevI][j][k].matrix[0][2] / deltaX2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 2));
		}

		if (ynumberGeneral > 1) {
			element = -c_theta_deltaT2 / deltaY2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));
		}

		if (znumberGeneral > 1) {
			element = -c_theta_deltaT2 / deltaZ2;
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
		}

		if ((ynumberGeneral > 1) && (xnumberGeneral > 1)) {
			element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[1][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 0));
			element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[1][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 1));
			element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[1][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 2));

			element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[1][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 0));
			element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[1][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 1));
			element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[1][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 2));

			element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[1][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 0));
			element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[1][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 1));
			element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[1][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 2));

			element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[1][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 0));
			element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[1][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 1));
			element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[1][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 2));
		}

		if ((znumberGeneral > 1) && (xnumberGeneral > 1)) {
			element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[2][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 0));
			element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[2][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 1));
			element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[2][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 2));

			element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[2][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 0));
			element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[2][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 1));
			element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[2][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 2));

			element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[2][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 0));
			element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[2][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 1));
			element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[2][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 2));

			element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[2][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 0));
			element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[2][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 1));
			element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[2][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 2));
		}
	}
}

void Simulation::createInternalEquationY(int i, int j, int k) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT) * speed_of_light_correction_sqr;
	//c_theta_deltaT2 = 0;
	double element = 1.0 - dielectricTensor[i][j][k].matrix[1][1];
	if (isInResistiveLayer(i, j, k)) {
		element += fakeCondactivity * theta * deltaT;
	}


	int nextI = i + 1;
	int prevI = i - 1;
	int nextJ = j + 1;
	int prevJ = j - 1;
	int nextK = k + 1;
	int prevK = k - 1;
	bool _27_points = false;
	//bool _27_points = true;

	if (_27_points) {
		if (ynumberGeneral > 1 && znumberGeneral > 1) {
			element += -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][k].matrix[1][1]) / (2 * deltaY2) - 1.0 / (2 * deltaX2)
				- 1.0 / (2 * deltaZ2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 1));

			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[nextI][j][k].matrix[1][1]) / (4 * deltaY2) + 1.0 / (4 *
				deltaX2) - 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, k, 1));
			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[prevI][j][k].matrix[1][1]) / (4 * deltaY2) + 1.0 / (4 *
				deltaX2) - 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, k, 1));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[i][nextJ][k].matrix[1][1]) / (4 * deltaY2) - 1.0 / (4 * deltaX2
			) - 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[i][prevJ][k].matrix[1][1]) / (4 * deltaY2) - 1.0 / (4 * deltaX2
			) - 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 1));

			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][nextK].matrix[1][1]) / (4 * deltaY2) - 1.0 / (4 *
				deltaX2) + 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, nextK, 1));
			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][prevK].matrix[1][1]) / (4 * deltaY2) - 1.0 / (4 *
				deltaX2) + 1.0 / (4 * deltaZ2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, prevK, 1));

			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[i][nextJ][nextK].matrix[1][1]) / (8 * deltaY2) - 1.0 / (8 *
				deltaX2) + 1.0 / (8 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[2][1] / (8 * deltaY *
				deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, nextK, 1));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[i][nextJ][prevK].matrix[1][1]) / (8 * deltaY2) - 1.0 / (8 *
				deltaX2) + 1.0 / (8 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[2][1] / (8 * deltaY *
				deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, prevK, 1));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[i][prevJ][nextK].matrix[1][1]) / (8 * deltaY2) - 1.0 / (8 *
				deltaX2) + 1.0 / (8 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[2][1] / (8 * deltaY *
				deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, nextK, 1));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[i][prevJ][prevK].matrix[1][1]) / (8 * deltaY2) - 1.0 / (8 *
				deltaX2) + 1.0 / (8 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[2][1] / (8 * deltaY *
				deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, prevK, 1));

			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[nextI][j][nextK].matrix[1][1]) / (8 * deltaY2) + 1.0 / (8 *
				deltaX2) + 1.0 / (8 * deltaZ2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, nextK, 1));
			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[nextI][j][prevK].matrix[1][1]) / (8 * deltaY2) + 1.0 / (8 *
				deltaX2) + 1.0 / (8 * deltaZ2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, prevK, 1));
			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[prevI][j][nextK].matrix[1][1]) / (8 * deltaY2) + 1.0 / (8 *
				deltaX2) + 1.0 / (8 * deltaZ2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, nextK, 1));
			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[prevI][j][prevK].matrix[1][1]) / (8 * deltaY2) + 1.0 / (8 *
				deltaX2) + 1.0 / (8 * deltaZ2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, prevK, 1));

			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][nextJ][k].matrix[1][1]) / (8 * deltaY2) + 1.0 / (8 *
				deltaX2) - 1.0 / (8 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[0][1] / (8 * deltaX *
				deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 1));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][prevJ][k].matrix[1][1]) / (8 * deltaY2) + 1.0 / (8 *
				deltaX2) - 1.0 / (8 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[0][1] / (8 * deltaX *
				deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 1));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][prevJ][k].matrix[1][1]) / (8 * deltaY2) + 1.0 / (8 *
				deltaX2) - 1.0 / (8 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[0][1] / (8 * deltaX *
				deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 1));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][nextJ][k].matrix[1][1]) / (8 * deltaY2) + 1.0 / (8 *
				deltaX2) - 1.0 / (8 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[0][1] / (8 * deltaX *
				deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 1));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][nextJ][nextK].matrix[1][1]) / (16 * deltaY2) + 1 / (16 *
				deltaX2) + 1 / (16 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][nextK].matrix[0][1] / (16 * deltaX
				* deltaY) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][nextK].matrix[2][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, nextK, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][nextJ][prevK].matrix[1][1]) / (16 * deltaY2) + 1 / (16 *
				deltaX2) + 1 / (16 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][prevK].matrix[0][1] / (16 * deltaX
				* deltaY) - c_theta_deltaT2 * dielectricTensor[nextI][nextJ][prevK].matrix[2][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, prevK, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][prevJ][nextK].matrix[1][1]) / (16 * deltaY2) + 1 / (16 *
				deltaX2) + 1 / (16 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][nextK].matrix[0][1] / (16 * deltaX
				* deltaY) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][nextK].matrix[2][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, nextK, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][prevJ][prevK].matrix[1][1]) / (16 * deltaY2) + 1 / (16 *
				deltaX2) + 1 / (16 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][prevK].matrix[0][1] / (16 * deltaX
				* deltaY) + c_theta_deltaT2 * dielectricTensor[nextI][prevJ][prevK].matrix[2][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, prevK, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][nextJ][nextK].matrix[1][1]) / (16 * deltaY2) + 1 / (16 *
				deltaX2) + 1 / (16 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][nextK].matrix[0][1] / (16 * deltaX
				* deltaY) + c_theta_deltaT2 * dielectricTensor[prevI][nextJ][nextK].matrix[2][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, nextK, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][nextJ][prevK].matrix[1][1]) / (16 * deltaY2) + 1 / (16 *
				deltaX2) + 1 / (16 * deltaZ2)) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][prevK].matrix[0][1] / (16 * deltaX
				* deltaY) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][prevK].matrix[2][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, prevK, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][prevJ][nextK].matrix[1][1]) / (16 * deltaY2) + 1 / (16 *
				deltaX2) + 1 / (16 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][nextK].matrix[0][1] / (16 * deltaX
				* deltaY) - c_theta_deltaT2 * dielectricTensor[prevI][prevJ][nextK].matrix[2][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, nextK, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][prevJ][prevK].matrix[1][1]) / (16 * deltaY2) + 1 / (16 *
				deltaX2) + 1 / (16 * deltaZ2)) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][prevK].matrix[0][1] / (16 * deltaX
				* deltaY) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][prevK].matrix[2][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, prevK, 1));

			//// x
			element = -dielectricTensor[i][j][k].matrix[1][0] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[1][0]) / (2 *
				deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 0));

			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[1][0]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, k, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[1][0]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, k, 0));

			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][k].matrix[1][0]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][k].matrix[1][0]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 0));

			element = -c_theta_deltaT2 * (dielectricTensor[i][j][nextK].matrix[1][0]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, nextK, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[i][j][prevK].matrix[1][0]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, prevK, 0));

			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][nextK].matrix[1][0]) / (8 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[i][nextJ][nextK].matrix[2][0] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][prevK].matrix[1][0]) / (8 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[i][nextJ][prevK].matrix[2][0] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, prevK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][nextK].matrix[1][0]) / (8 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[i][prevJ][nextK].matrix[2][0] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][prevK].matrix[1][0]) / (8 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[i][prevJ][prevK].matrix[2][0] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, prevK, 0));

			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][nextK].matrix[1][0]) / (8 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, nextK, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][prevK].matrix[1][0]) / (8 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, prevK, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][nextK].matrix[1][0]) / (8 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, nextK, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][prevK].matrix[1][0]) / (8 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, prevK, 0));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][k].matrix[1][0]) / (8 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][k].matrix[0][0] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][k].matrix[1][0]) / (8 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][k].matrix[0][0] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][k].matrix[1][0]) / (8 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][k].matrix[0][0] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][k].matrix[1][0]) / (8 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][k].matrix[0][0] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 0));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][nextK].matrix[1][0]) / (16 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][nextK].matrix[0][0] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][nextK].matrix[2][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][prevK].matrix[1][0]) / (16 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][prevK].matrix[0][0] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][prevK].matrix[2][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, prevK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][nextK].matrix[1][0]) / (16 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][nextK].matrix[0][0] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][nextK].matrix[2][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][prevK].matrix[1][0]) / (16 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][prevK].matrix[0][0] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][prevK].matrix[2][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, prevK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][nextK].matrix[1][0]) / (16 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][nextK].matrix[0][0] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][nextK].matrix[2][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][prevK].matrix[1][0]) / (16 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][prevK].matrix[0][0] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][prevK].matrix[2][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, prevK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][nextK].matrix[1][0]) / (16 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][nextK].matrix[0][0] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][nextK].matrix[2][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][prevK].matrix[1][0]) / (16 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][prevK].matrix[0][0] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][prevK].matrix[2][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, prevK, 0));

			//// z
			element = -dielectricTensor[i][j][k].matrix[1][2] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[1][2]) / (2 *
				deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 2));

			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[1][2]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, k, 2));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[1][2]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, k, 2));

			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][k].matrix[1][2]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][k].matrix[1][2]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 2));

			element = -c_theta_deltaT2 * (dielectricTensor[i][j][nextK].matrix[1][2]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, nextK, 2));
			element = -c_theta_deltaT2 * (dielectricTensor[i][j][prevK].matrix[1][2]) / (4 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, prevK, 2));

			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][nextK].matrix[1][2]) / (8 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[i][nextJ][nextK].matrix[2][2] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][prevK].matrix[1][2]) / (8 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[i][nextJ][prevK].matrix[2][2] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, prevK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][nextK].matrix[1][2]) / (8 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[i][prevJ][nextK].matrix[2][2] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][prevK].matrix[1][2]) / (8 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[i][prevJ][prevK].matrix[2][2] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, prevK, 2));

			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][nextK].matrix[1][2]) / (8 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, nextK, 2));
			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][prevK].matrix[1][2]) / (8 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, prevK, 2));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][nextK].matrix[1][2]) / (8 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, nextK, 2));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][prevK].matrix[1][2]) / (8 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, prevK, 2));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][k].matrix[1][2]) / (8 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][k].matrix[0][2] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][k].matrix[1][2]) / (8 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][k].matrix[0][2] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][k].matrix[1][2]) / (8 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][k].matrix[0][2] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][k].matrix[1][2]) / (8 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][k].matrix[0][2] / (8 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 2));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][nextK].matrix[1][2]) / (16 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][nextK].matrix[0][2] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][nextK].matrix[2][2] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][prevK].matrix[1][2]) / (16 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][prevK].matrix[0][2] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][prevK].matrix[2][2] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, prevK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][nextK].matrix[1][2]) / (16 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][nextK].matrix[0][2] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][nextK].matrix[2][2] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][prevK].matrix[1][2]) / (16 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][prevK].matrix[0][2] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][prevK].matrix[2][2] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, prevK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][nextK].matrix[1][2]) / (16 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][nextK].matrix[0][2] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][nextK].matrix[2][2] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][prevK].matrix[1][2]) / (16 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][prevK].matrix[0][2] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][prevK].matrix[2][2] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, prevK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][nextK].matrix[1][2]) / (16 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][nextK].matrix[0][2] / (16 * deltaX * deltaY) - c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][nextK].matrix[2][2] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, nextK, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][prevK].matrix[1][2]) / (16 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][prevK].matrix[0][2] / (16 * deltaX * deltaY) + c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][prevK].matrix[2][2] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, prevK, 2));

		} else if (ynumberGeneral > 1) {
			element += -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][k].matrix[1][1]) / (deltaY2) - 1.0 / (deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 1));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][k].matrix[1][1]) / (2 * deltaY2) - 1.0 / (2 * deltaX2
			));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, k, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][k].matrix[1][1]) / (2 * deltaY2) - 1.0 / (2 * deltaX2
			));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, k, 1));

			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][nextJ][k].matrix[1][1]) / (2 * deltaY2) + 1.0 / (2 *
				deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 1));
			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][prevJ][k].matrix[1][1]) / (2 * deltaY2) + 1.0 / (2 *
				deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 1));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][nextJ][k].matrix[1][1]) / (4 * deltaY2) + 1.0 / (4 *
				deltaX2)) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[0][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][prevJ][k].matrix[1][1]) / (4 * deltaY2) + 1.0 / (4 *
				deltaX2)) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[0][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][nextJ][k].matrix[1][1]) / (4 * deltaY2) + 1.0 / (4 *
				deltaX2)) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[0][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 1));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][prevJ][k].matrix[1][1]) / (4 * deltaY2) + 1.0 / (4 *
				deltaX2)) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[0][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 1));

			/////x

			element = -dielectricTensor[i][j][k].matrix[1][0] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[1][0]) / (
				deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[1][0]) / (2 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, k, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[1][0]) / (2 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, k, 0));

			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][k].matrix[1][0]) / (2 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][k].matrix[1][0]) / (2 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 0));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][k].matrix[1][0]) / (4 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][k].matrix[0][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][k].matrix[1][0]) / (4 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][k].matrix[0][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][k].matrix[1][0]) / (4 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][k].matrix[0][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][k].matrix[1][0]) / (4 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][k].matrix[0][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 0));

			////z
			element = -dielectricTensor[i][j][k].matrix[1][2] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[1][2]) / (
				deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 2));
			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[1][2]) / (2 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, k, 2));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[1][2]) / (2 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, k, 2));

			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][k].matrix[1][2]) / (2 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][k].matrix[1][2]) / (2 * deltaY2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 2));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][k].matrix[1][2]) / (4 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][k].matrix[0][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][k].matrix[1][2]) / (4 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][k].matrix[0][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][k].matrix[1][2]) / (4 * deltaY2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][k].matrix[0][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 2));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][k].matrix[1][2]) / (4 * deltaY2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][k].matrix[0][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 2));
		} else if (znumberGeneral > 1) {
			element += -c_theta_deltaT2 * (-1.0 / (deltaZ2) - 1.0 / (deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 1));

			element = -c_theta_deltaT2 * (-1.0 / (2 * deltaZ2) + 1.0 / (2 * deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, k, 1));
			element = -c_theta_deltaT2 * (-1.0 / (2 * deltaZ2) + 1.0 / (2 * deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, k, 1));

			element = -c_theta_deltaT2 * (1.0 / (2 * deltaZ2) - 1.0 / (2 * deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, nextK, 1));
			element = -c_theta_deltaT2 * (1.0 / (2 * deltaZ2) - 1.0 / (2 * deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, prevK, 1));

			element = -c_theta_deltaT2 * (1.0 / (4 * deltaZ2) + 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, nextK, 1));
			element = -c_theta_deltaT2 * (1.0 / (4 * deltaZ2) + 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, prevK, 1));
			element = -c_theta_deltaT2 * (1.0 / (4 * deltaZ2) + 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, nextK, 1));
			element = -c_theta_deltaT2 * (1. / (4 * deltaZ2) + 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, prevK, 1));

			////x
			element = -dielectricTensor[i][j][k].matrix[1][0];
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 0));

			////z
			element = -dielectricTensor[i][j][k].matrix[1][2];
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 2));

		} else {
			element += -c_theta_deltaT2 * (-2.0 / (deltaX2));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 1));

			element = -dielectricTensor[i][j][k].matrix[1][0];
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 0));

			element = -dielectricTensor[i][j][k].matrix[1][2];
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 2));

			element = -c_theta_deltaT2 * (1.0 / deltaX2);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, k, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, k, 1));
		}
		///
		/*for(int tempI = 0; tempI < 2*splineOrder + 3; ++ tempI){
		for(int tempJ = 0; tempJ < 2*splineOrder + 3; ++tempJ){
			for(int tempK = 0; tempK < 2*splineOrder + 3; ++tempK){
					int xindex = massMatrix[i][j][k].xindex[tempI];
					int yindex = massMatrix[i][j][k].yindex[tempJ];
					int zindex = massMatrix[i][j][k].zindex[tempK];
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[1][0];
					maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, xindex, yindex, zindex, 0));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[1][1];
					maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, xindex, yindex, zindex, 1));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[1][2];
					maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, xindex, yindex, zindex, 2));
			}
			}
		}*/
		///
	} else {
		if (xnumberGeneral > 1) {
			element += c_theta_deltaT2 * (2.0 / deltaX2);
		}
		if (ynumberGeneral > 1) {
			element += c_theta_deltaT2 * (2.0 - 2 * dielectricTensor[i][j][k].matrix[1][1]) / deltaY2;
		}
		if (znumberGeneral > 1) {
			element += c_theta_deltaT2 * 2.0 / deltaZ2;
		}
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 1));

		element = -dielectricTensor[i][j][k].matrix[1][0];
		if (ynumberGeneral > 1) {
			element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[1][0] / deltaY2;
		}
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 0));

		element = -dielectricTensor[i][j][k].matrix[1][2];
		if (ynumberGeneral > 1) {
			element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[1][2] / deltaY2;
		}
		maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 2));

		if (ynumberGeneral > 1) {
			element = c_theta_deltaT2 * (-1.0 + dielectricTensor[i][nextJ][k].matrix[1][1]) / deltaY2;
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 1));
			element = c_theta_deltaT2 * dielectricTensor[i][nextJ][k].matrix[1][0] / deltaY2;
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 0));
			element = c_theta_deltaT2 * dielectricTensor[i][nextJ][k].matrix[1][2] / deltaY2;
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 2));

			element = c_theta_deltaT2 * (-1.0 + dielectricTensor[i][prevJ][k].matrix[1][1]) / deltaY2;
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 1));
			element = c_theta_deltaT2 * dielectricTensor[i][prevJ][k].matrix[1][0] / deltaY2;
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 0));
			element = c_theta_deltaT2 * dielectricTensor[i][prevJ][k].matrix[1][2] / deltaY2;
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 2));
		}

		if (xnumberGeneral > 1) {
			element = -c_theta_deltaT2 / deltaX2;
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, k, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, k, 1));
		}

		if (znumberGeneral > 1) {
			element = -c_theta_deltaT2 / deltaZ2;
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, nextK, 1));
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, prevK, 1));
		}

		if ((ynumberGeneral > 1) && (xnumberGeneral > 1)) {
			element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[0][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 0));
			element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[0][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 1));
			element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[0][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 2));

			element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[0][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 0));
			element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[0][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 1));
			element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[0][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 2));

			element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[0][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 0));
			element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[0][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 1));
			element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[0][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 2));

			element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[0][0] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 0));
			element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[0][1] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 1));
			element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[0][2] / (4 * deltaX * deltaY);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 2));
		}

		if ((znumberGeneral > 1) && (ynumberGeneral > 1)) {
			element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[2][0] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, nextK, 0));
			element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[2][1] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, nextK, 1));
			element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[2][2] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, nextK, 2));

			element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[2][0] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, prevK, 0));
			element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[2][1] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, prevK, 1));
			element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[2][2] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, prevK, 2));

			element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[2][0] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, nextK, 0));
			element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[2][1] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, nextK, 1));
			element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[2][2] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, nextK, 2));

			element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[2][0] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, prevK, 0));
			element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[2][1] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, prevK, 1));
			element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[2][2] / (4 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, prevK, 2));
		}
	}
}

void Simulation::createInternalEquationZ(int i, int j, int k) {
	double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT) * speed_of_light_correction_sqr;
	double element = 1.0 - dielectricTensor[i][j][k].matrix[2][2];
	//c_theta_deltaT2 = 0;
	if (isInResistiveLayer(i, j, k)) {
		element += fakeCondactivity * theta * deltaT;
	}
	int nextI = i + 1;
	int prevI = i - 1;
	int nextJ = j + 1;
	int prevJ = j - 1;
	int nextK = k + 1;
	int prevK = k - 1;
	bool _27_points = false;
	//bool _27_points = true;

	if (_27_points) {
		if (ynumberGeneral > 1 && znumberGeneral > 1) {
			element += -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][k].matrix[2][2]) / (2 * deltaZ2) - 1.0 / (2 * deltaY2)
				- 1.0 / (2 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 2));

			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[nextI][j][k].matrix[2][2]) / (4 * deltaZ2) - 1.0 / (4 *
				deltaY2) + 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, k, 2));
			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[prevI][j][k].matrix[2][2]) / (4 * deltaZ2) - 1.0 / (4 *
				deltaY2) + 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, k, 2));

			element = -c_theta_deltaT2 * (- (1.0 - dielectricTensor[i][nextJ][k].matrix[2][2]) / (4 * deltaZ2) + 1.0 / (4 *
				deltaY2) - 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, k, 2));
			element = -c_theta_deltaT2 * (- (1.0 - dielectricTensor[i][prevJ][k].matrix[2][2]) / (4 * deltaZ2) + 1.0 / (4 *
				deltaY2) - 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, k, 2));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[i][j][nextK].matrix[2][2]) / (4 * deltaZ2) - 1.0 / (4 * deltaY2
			) - 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[i][j][prevK].matrix[2][2]) / (4 * deltaZ2) - 1.0 / (4 * deltaY2
			) - 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 2));

			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[i][nextJ][nextK].matrix[2][2]) / (8 * deltaZ2) + 1.0 / (8 *
				deltaY2) - 1.0 / (8 * deltaX2)) + c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[1][2] / (8 * deltaZ *
				deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, nextK, 2));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[i][nextJ][prevK].matrix[2][2]) / (8 * deltaZ2) + 1.0 / (8 *
				deltaY2) - 1.0 / (8 * deltaX2)) - c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[1][2] / (8 * deltaZ *
				deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, prevK, 2));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[i][prevJ][nextK].matrix[2][2]) / (8 * deltaZ2) + 1.0 / (8 *
				deltaY2) - 1.0 / (8 * deltaX2)) - c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[1][2] / (8 * deltaZ *
				deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, nextK, 2));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[i][prevJ][prevK].matrix[2][2]) / (8 * deltaZ2) + 1.0 / (8 *
				deltaY2) - 1.0 / (8 * deltaX2)) + c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[1][2] / (8 * deltaZ *
				deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, prevK, 2));

			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][nextK].matrix[2][2]) / (8 * deltaZ2) - 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaX2)) + c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[0][2] / (8 * deltaX *
				deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 2));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][prevK].matrix[2][2]) / (8 * deltaZ2) - 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaX2)) - c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[0][2] / (8 * deltaX *
				deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 2));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][nextK].matrix[2][2]) / (8 * deltaZ2) - 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaX2)) - c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[0][2] / (8 * deltaX *
				deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 2));
			element = - c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][prevK].matrix[2][2]) / (8 * deltaZ2) - 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaX2)) + c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[0][2] / (8 * deltaX *
				deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 2));

			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[nextI][nextJ][k].matrix[2][2]) / (8 * deltaZ2) + 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, nextJ, k, 2));
			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[nextI][prevJ][k].matrix[2][2]) / (8 * deltaZ2) + 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, prevJ, k, 2));
			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[prevI][prevJ][k].matrix[2][2]) / (8 * deltaZ2) + 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, prevJ, k, 2));
			element = - c_theta_deltaT2 * (-(1.0 - dielectricTensor[prevI][nextJ][k].matrix[2][2]) / (8 * deltaZ2) + 1.0 / (8 *
				deltaY2) + 1.0 / (8 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, nextJ, k, 2));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][nextJ][nextK].matrix[2][2]) / (16 * deltaZ2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaX2)) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][nextK].matrix[1][2] / (16 * deltaZ
				* deltaY) + c_theta_deltaT2 * dielectricTensor[nextI][nextJ][nextK].matrix[1][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, nextJ, nextK, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][nextJ][prevK].matrix[2][2]) / (16 * deltaZ2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaX2)) - c_theta_deltaT2 * dielectricTensor[nextI][nextJ][prevK].matrix[1][2] / (16 * deltaZ
				* deltaY) - c_theta_deltaT2 * dielectricTensor[nextI][nextJ][prevK].matrix[1][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, nextJ, prevK, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][prevJ][nextK].matrix[2][2]) / (16 * deltaZ2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaX2)) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][nextK].matrix[1][2] / (16 * deltaZ
				* deltaY) + c_theta_deltaT2 * dielectricTensor[nextI][prevJ][nextK].matrix[1][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, prevJ, nextK, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][prevJ][prevK].matrix[2][2]) / (16 * deltaZ2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaX2)) + c_theta_deltaT2 * dielectricTensor[nextI][prevJ][prevK].matrix[1][2] / (16 * deltaZ
				* deltaY) - c_theta_deltaT2 * dielectricTensor[nextI][prevJ][prevK].matrix[1][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, prevJ, prevK, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][nextJ][nextK].matrix[2][2]) / (16 * deltaZ2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaX2)) + c_theta_deltaT2 * dielectricTensor[prevI][nextJ][nextK].matrix[1][2] / (16 * deltaZ
				* deltaY) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][nextK].matrix[1][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, nextJ, nextK, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][nextJ][prevK].matrix[2][2]) / (16 * deltaZ2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaX2)) - c_theta_deltaT2 * dielectricTensor[prevI][nextJ][prevK].matrix[1][2] / (16 * deltaZ
				* deltaY) + c_theta_deltaT2 * dielectricTensor[prevI][nextJ][prevK].matrix[1][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, nextJ, prevK, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][prevJ][nextK].matrix[2][2]) / (16 * deltaZ2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaX2)) - c_theta_deltaT2 * dielectricTensor[prevI][prevJ][nextK].matrix[1][2] / (16 * deltaZ
				* deltaY) - c_theta_deltaT2 * dielectricTensor[prevI][prevJ][nextK].matrix[1][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, prevJ, nextK, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][prevJ][prevK].matrix[2][2]) / (16 * deltaZ2) + 1 / (16 *
				deltaY2) + 1 / (16 * deltaX2)) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][prevK].matrix[1][2] / (16 * deltaX
				* deltaY) + c_theta_deltaT2 * dielectricTensor[prevI][prevJ][prevK].matrix[1][2] / (16 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, prevJ, prevK, 2));

			////x
			element = -dielectricTensor[i][j][k].matrix[2][0] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[2][0]) / (2 *
				deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 0));

			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[2][0]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, k, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[2][0]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, k, 0));

			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][k].matrix[2][0]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][k].matrix[2][0]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, k, 0));

			element = -c_theta_deltaT2 * (dielectricTensor[i][j][nextK].matrix[2][0]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[i][j][prevK].matrix[2][0]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 0));

			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][nextK].matrix[2][0]) / (8 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[i][nextJ][nextK].matrix[1][0] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][prevK].matrix[2][0]) / (8 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[i][nextJ][prevK].matrix[1][0] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, prevK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][nextK].matrix[2][0]) / (8 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[i][prevJ][nextK].matrix[1][0] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][prevK].matrix[2][0]) / (8 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[i][prevJ][prevK].matrix[1][0] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, prevK, 0));

			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][nextK].matrix[2][0]) / (8 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[nextI][j][nextK].matrix[0][0] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][prevK].matrix[2][0]) / (8 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[nextI][j][prevK].matrix[0][0] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][nextK].matrix[2][0]) / (8 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[prevI][j][nextK].matrix[0][0] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][prevK].matrix[2][0]) / (8 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[prevI][j][prevK].matrix[0][0] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 0));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][k].matrix[2][0]) / (8 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, nextJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][k].matrix[2][0]) / (8 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, prevJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][k].matrix[2][0]) / (8 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, prevJ, k, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][k].matrix[2][0]) / (8 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, nextJ, k, 0));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][nextK].matrix[2][0]) / (16 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][nextK].matrix[0][0] / (16 * deltaX * deltaZ) + c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][nextK].matrix[1][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, nextJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][prevK].matrix[2][0]) / (16 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][prevK].matrix[0][0] / (16 * deltaX * deltaZ) - c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][prevK].matrix[1][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, nextJ, prevK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][nextK].matrix[2][0]) / (16 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][nextK].matrix[0][0] / (16 * deltaX * deltaZ) - c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][nextK].matrix[1][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, prevJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][prevK].matrix[2][0]) / (16 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][prevK].matrix[0][0] / (16 * deltaX * deltaZ) + c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][prevK].matrix[1][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, prevJ, prevK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][nextK].matrix[2][0]) / (16 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][nextK].matrix[0][0] / (16 * deltaX * deltaZ) + c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][nextK].matrix[1][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, nextJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][prevK].matrix[2][0]) / (16 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][prevK].matrix[0][0] / (16 * deltaX * deltaZ) - c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][prevK].matrix[1][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, nextJ, prevK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][nextK].matrix[2][0]) / (16 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][nextK].matrix[0][0] / (16 * deltaX * deltaZ) - c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][nextK].matrix[1][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, prevJ, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][prevK].matrix[2][0]) / (16 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][prevK].matrix[0][0] / (16 * deltaX * deltaZ) + c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][prevK].matrix[1][0] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, prevJ, prevK, 0));

			////y

			element = -dielectricTensor[i][j][k].matrix[2][1] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[2][1]) / (2 *
				deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 1));

			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[2][1]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, k, 1));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[2][1]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, k, 1));

			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][k].matrix[2][1]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][k].matrix[2][1]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, k, 1));

			element = -c_theta_deltaT2 * (dielectricTensor[i][j][nextK].matrix[2][1]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 1));
			element = -c_theta_deltaT2 * (dielectricTensor[i][j][prevK].matrix[2][1]) / (4 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 1));

			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][nextK].matrix[2][1]) / (8 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[i][nextJ][nextK].matrix[1][1] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[i][nextJ][prevK].matrix[2][1]) / (8 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[i][nextJ][prevK].matrix[1][1] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, prevK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][nextK].matrix[2][1]) / (8 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[i][prevJ][nextK].matrix[1][1] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[i][prevJ][prevK].matrix[2][1]) / (8 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[i][prevJ][prevK].matrix[1][1] / (8 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, prevK, 1));

			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][nextK].matrix[2][1]) / (8 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[nextI][j][nextK].matrix[0][1] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 1));
			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][prevK].matrix[2][1]) / (8 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[nextI][j][prevK].matrix[0][1] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 1));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][nextK].matrix[2][1]) / (8 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[prevI][j][nextK].matrix[0][1] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 1));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][prevK].matrix[2][1]) / (8 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[prevI][j][prevK].matrix[0][1] / (8 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 1));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][k].matrix[2][1]) / (8 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, nextJ, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][k].matrix[2][1]) / (8 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, prevJ, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][k].matrix[2][1]) / (8 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, prevJ, k, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][k].matrix[2][1]) / (8 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, nextJ, k, 1));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][nextK].matrix[2][1]) / (16 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][nextK].matrix[0][1] / (16 * deltaX * deltaZ) + c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][nextK].matrix[1][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, nextJ, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][nextJ][prevK].matrix[2][1]) / (16 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[nextI][nextJ][prevK].matrix[0][1] / (16 * deltaX * deltaZ) - c_theta_deltaT2 * dielectricTensor[
					nextI][nextJ][prevK].matrix[1][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, nextJ, prevK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][nextK].matrix[2][1]) / (16 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][nextK].matrix[0][1] / (16 * deltaX * deltaZ) - c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][nextK].matrix[1][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, prevJ, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][prevJ][prevK].matrix[2][1]) / (16 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[nextI][prevJ][prevK].matrix[0][1] / (16 * deltaX * deltaZ) + c_theta_deltaT2 * dielectricTensor[
					nextI][prevJ][prevK].matrix[1][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, prevJ, prevK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][nextK].matrix[2][1]) / (16 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][nextK].matrix[0][1] / (16 * deltaX * deltaZ) + c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][nextK].matrix[1][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, nextJ, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][nextJ][prevK].matrix[2][1]) / (16 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[prevI][nextJ][prevK].matrix[0][1] / (16 * deltaX * deltaZ) - c_theta_deltaT2 * dielectricTensor[
					prevI][nextJ][prevK].matrix[1][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, nextJ, prevK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][nextK].matrix[2][1]) / (16 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][nextK].matrix[0][1] / (16 * deltaX * deltaZ) - c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][nextK].matrix[1][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, prevJ, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][prevJ][prevK].matrix[2][1]) / (16 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[prevI][prevJ][prevK].matrix[0][1] / (16 * deltaX * deltaZ) + c_theta_deltaT2 * dielectricTensor[
					prevI][prevJ][prevK].matrix[1][1] / (16 * deltaY * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, prevJ, prevK, 1));
		} else if (ynumberGeneral > 1) {
			element += -c_theta_deltaT2 * (-1.0 / (deltaY2) - 1.0 / (deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 2));

			element = -c_theta_deltaT2 * (-1.0 / (2 * deltaY2) + 1.0 / (2 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, k, 2));
			element = -c_theta_deltaT2 * (-1.0 / (2 * deltaY2) + 1.0 / (2 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, k, 2));

			element = -c_theta_deltaT2 * (1.0 / (2 * deltaY2) - 1.0 / (2 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, k, 2));
			element = -c_theta_deltaT2 * (1.0 / (2 * deltaY2) - 1.0 / (2 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, k, 2));

			element = -c_theta_deltaT2 * (1.0 / (4 * deltaY2) + 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, nextJ, k, 2));
			element = -c_theta_deltaT2 * (1.0 / (4 * deltaY2) + 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, prevJ, k, 2));
			element = -c_theta_deltaT2 * (1.0 / (4 * deltaY2) + 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, nextJ, k, 2));
			element = -c_theta_deltaT2 * (1.0 / (4 * deltaY2) + 1.0 / (4 * deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, prevJ, k, 2));

			////x
			element = -dielectricTensor[i][j][k].matrix[2][0];
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 0));

			////z
			element = -dielectricTensor[i][j][k].matrix[2][1];
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 1));
		} else if (znumberGeneral > 1) {
			element += -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][k].matrix[2][2]) / (deltaZ2) - 1.0 / (deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 2));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][k].matrix[2][2]) / (2 * deltaZ2) - 1.0 / (2 * deltaX2
			));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, k, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][k].matrix[2][2]) / (2 * deltaZ2) - 1.0 / (2 * deltaX2
			));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, k, 2));

			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][nextK].matrix[2][2]) / (2 * deltaZ2) + 1.0 / (2 *
				deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 2));
			element = -c_theta_deltaT2 * (-(1.0 - dielectricTensor[i][j][prevK].matrix[2][2]) / (2 * deltaZ2) + 1.0 / (2 *
				deltaX2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 2));

			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][nextK].matrix[2][2]) / (4 * deltaZ2) + 1.0 / (4 *
				deltaX2)) + c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[0][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[nextI][j][prevK].matrix[2][2]) / (4 * deltaZ2) + 1.0 / (4 *
				deltaX2)) - c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[0][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][nextK].matrix[2][2]) / (4 * deltaZ2) + 1.0 / (4 *
				deltaX2)) - c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[0][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 2));
			element = -c_theta_deltaT2 * ((1.0 - dielectricTensor[prevI][j][prevK].matrix[2][2]) / (4 * deltaZ2) + 1.0 / (4 *
				deltaX2)) + c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[0][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 2));

			/////x

			element = -dielectricTensor[i][j][k].matrix[2][0] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[2][0]) / (
				deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[2][0]) / (2 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, k, 0));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[2][0]) / (2 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, k, 0));

			element = c_theta_deltaT2 * (dielectricTensor[i][j][nextK].matrix[2][0]) / (2 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[i][j][prevK].matrix[2][0]) / (2 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 0));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][nextK].matrix[2][0]) / (4 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[nextI][j][nextK].matrix[0][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][prevK].matrix[2][0]) / (4 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[nextI][j][prevK].matrix[0][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][nextK].matrix[2][0]) / (4 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[prevI][j][nextK].matrix[0][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 0));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][prevK].matrix[2][0]) / (4 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[prevI][j][prevK].matrix[0][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 0));

			////y
			element = -dielectricTensor[i][j][k].matrix[2][1] - c_theta_deltaT2 * (dielectricTensor[i][j][k].matrix[2][1]) / (
				deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 1));
			element = -c_theta_deltaT2 * (dielectricTensor[nextI][j][k].matrix[2][1]) / (2 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, k, 1));
			element = -c_theta_deltaT2 * (dielectricTensor[prevI][j][k].matrix[2][1]) / (2 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, k, 1));

			element = c_theta_deltaT2 * (dielectricTensor[i][j][nextK].matrix[2][1]) / (2 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[i][j][prevK].matrix[2][1]) / (2 * deltaZ2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 1));

			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][nextK].matrix[2][1]) / (4 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[nextI][j][nextK].matrix[0][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[nextI][j][prevK].matrix[2][1]) / (4 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[nextI][j][prevK].matrix[0][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][nextK].matrix[2][1]) / (4 * deltaZ2) - c_theta_deltaT2 *
				dielectricTensor[prevI][j][nextK].matrix[0][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 1));
			element = c_theta_deltaT2 * (dielectricTensor[prevI][j][prevK].matrix[2][1]) / (4 * deltaZ2) + c_theta_deltaT2 *
				dielectricTensor[prevI][j][prevK].matrix[0][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, k, 1));
		} else {
			element += c_theta_deltaT2 * (2.0 / deltaX2);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 2));

			element = -dielectricTensor[i][j][k].matrix[2][0];
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 0));

			element = -dielectricTensor[i][j][k].matrix[2][1];
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 1));

			element = -c_theta_deltaT2 / deltaX2;
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, k, 2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, k, 2));
		}
		///
		/*for(int tempI = 0; tempI < 2*splineOrder + 3; ++ tempI){
		for(int tempJ = 0; tempJ < 2*splineOrder + 3; ++tempJ){
			for(int tempK = 0; tempK < 2*splineOrder + 3; ++tempK){
					int xindex = massMatrix[i][j][k].xindex[tempI];
					int yindex = massMatrix[i][j][k].yindex[tempJ];
					int zindex = massMatrix[i][j][k].zindex[tempK];
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[2][0];
					maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, xindex, yindex, zindex, 0));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[2][1];
					maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, xindex, yindex, zindex, 1));
					element = massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[2][2];
					maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, xindex, yindex, zindex, 2));
			}
			}
		}*/
		///
	} else {
		if (xnumberGeneral > 1) {
			element += c_theta_deltaT2 * (2.0 / deltaX2);
		}
		if (ynumberGeneral > 1) {
			element += c_theta_deltaT2 * 2.0 / deltaY2;
		}
		if (znumberGeneral > 1) {
			element += c_theta_deltaT2 * (2.0 - 2 * dielectricTensor[i][j][k].matrix[2][2]) / deltaZ2;
		}
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 2));

		element = -dielectricTensor[i][j][k].matrix[2][0];
		if (znumberGeneral > 1) {
			element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[2][0] / deltaZ2;
		}
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 0));

		element = -dielectricTensor[i][j][k].matrix[2][1];
		if (znumberGeneral > 1) {
			element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[2][1] / deltaZ2;
		}
		maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 1));

		if (znumberGeneral > 1) {
			element = c_theta_deltaT2 * (-1.0 + dielectricTensor[i][j][nextK].matrix[2][2]) / deltaZ2;
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 2));
			element = c_theta_deltaT2 * dielectricTensor[i][j][nextK].matrix[2][1] / deltaZ2;
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 1));
			element = c_theta_deltaT2 * dielectricTensor[i][j][nextK].matrix[2][0] / deltaZ2;
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 0));

			element = c_theta_deltaT2 * (-1.0 + dielectricTensor[i][j][prevK].matrix[2][2]) / deltaZ2;
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 2));
			element = c_theta_deltaT2 * dielectricTensor[i][j][prevK].matrix[2][1] / deltaZ2;
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 1));
			element = c_theta_deltaT2 * dielectricTensor[i][j][prevK].matrix[2][0] / deltaZ2;
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 0));
		}

		if (ynumberGeneral > 1) {
			element = -c_theta_deltaT2 / deltaY2;
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, k, 2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, k, 2));
		}

		if (xnumberGeneral > 1) {
			element = -c_theta_deltaT2 / deltaX2;
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, k, 2));
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, k, 2));
		}

		if ((ynumberGeneral > 1) && (znumberGeneral > 1)) {
			element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[1][0] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, nextK, 0));
			element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[1][1] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, nextK, 1));
			element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[1][2] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, nextK, 2));

			element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[1][0] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, prevK, 0));
			element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[1][1] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, prevK, 1));
			element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[1][2] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, prevK, 2));

			element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[1][0] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, nextK, 0));
			element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[1][1] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, nextK, 1));
			element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[1][2] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, nextK, 2));

			element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[1][0] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, prevK, 0));
			element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[1][1] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, prevK, 1));
			element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[1][2] / (4 * deltaZ * deltaY);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, prevK, 2));
		}

		if ((znumberGeneral > 1) && (xnumberGeneral > 1)) {
			element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[0][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 0));
			element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[0][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 1));
			element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[0][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 2));

			element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[0][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 0));
			element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[0][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 1));
			element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[0][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 2));

			element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[0][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 0));
			element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[0][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 1));
			element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[0][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 2));

			element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[0][0] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 0));
			element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[0][1] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 1));
			element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[0][2] / (4 * deltaX * deltaZ);
			maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 2));
		}
	}
}

void Simulation::createInternalEquation(int i, int j, int k) {
	Vector3d rightPart = Efield[i][j][k];
	Vector3d rightPart2 = Bfield[i][j][k];

	//rightPart = rightPart + (evaluateRotB(i)* speed_of_light_normalized - electricFlux[i]*4*pi/fieldScale) * (theta * deltaT);
	Vector3d rotorB = evaluateRotB(i, j, k);
	//alertNaNOrInfinity(rotorB.x, "rotorB x = NaN in create internal\n");
	//alertNaNOrInfinity(rotorB.y, "rotorB y = NaN in create internal\n");
	//alertNaNOrInfinity(rotorB.z, "rotorB z = NaN in create internal\n");

	//alertNaNOrInfinity(electricFlux[i][j][k].x, "electricFlux[i][j][k] x = NaN in create internal\n");
	//alertNaNOrInfinity(electricFlux[i][j][k].y, "electricFlux[i][j][k] y = NaN in create internal\n");
	//alertNaNOrInfinity(electricFlux[i][j][k].z, "electricFlux[i][j][k] z = NaN in create internal\n");

	Vector3d gradDensity = evaluateGradDensity(i, j, k);

	//alertNaNOrInfinity(gradDensity.x, "gradDensity x = NaN in create internal\n");
	//alertNaNOrInfinity(gradDensity.y, "gradDensity y = NaN in create internal\n");
	//alertNaNOrInfinity(gradDensity.z, "gradDensity z = NaN in create internal\n");
	if (solverType == IMPLICIT) {
		rightPart = rightPart +
			(rotorB * speed_of_light_normalized * speed_of_light_correction - (electricFlux[i][j][k] * 4 * pi)) * (theta * deltaT
			) -
			(gradDensity * speed_of_light_normalized_sqr * speed_of_light_correction_sqr * theta * theta * deltaT * deltaT * 4 *
				pi);
	} else if (solverType == IMPLICIT_EC) {
		rightPart = rightPart +
			(rotorB * speed_of_light_normalized - (electricFlux[i][j][k] * 4 * pi)) * (theta * deltaT);
	} else {
		printf("wrong soler type in create internal equation\n");
		MPI_Finalize();
		exit(0);
	}
	if (isInResistiveLayer(i, j, k)) {
		rightPart = rightPart + E0 * (theta * deltaT * fakeCondactivity);
	}

	//rightPart = rightPart + evaluateRotB(i)*speed_of_light_normalized*theta*deltaT - electricFlux[i]*4*pi*theta*deltaT/fieldScale;
	if (solverType == IMPLICIT) {
		createInternalEquationX(i, j, k);
		createInternalEquationY(i, j, k);
		createInternalEquationZ(i, j, k);
	} else if (solverType == IMPLICIT_EC) {
		createInternalECEquationX(i, j, k);
		createInternalECEquationY(i, j, k);
		createInternalECEquationZ(i, j, k);
	} else {
		printf("wrong solver type in create internal equation\n");
		MPI_Finalize();
		exit(0);
	}

	if (fabs(rightPart.x) > 1E100) {
		printf("too large maxweell rightPart x\n");
		MPI_Finalize();
		exit(0);
	}

	if (fabs(rightPart.y) > 1E100) {
		printf("too large maxweell rightPart y\n");
		MPI_Finalize();
		exit(0);
	}

	if (fabs(rightPart.z) > 1E100) {
		printf("too large maxweell rightPart z\n");
		MPI_Finalize();
		exit(0);
	}

	alertNaNOrInfinity(rightPart.x, "right part x = NaN in create internal");
	alertNaNOrInfinity(rightPart.y, "right part y = NaN");
	alertNaNOrInfinity(rightPart.z, "right part z = NaN");

	//alertNaNOrInfinity(rightPart2.x, "right part 2 x = NaN");
	//alertNaNOrInfinity(rightPart2.y, "right part 2 y = NaN");
	//alertNaNOrInfinity(rightPart2.z, "right part 2 z = NaN");


	maxwellEquationRightPart[i][j][k][0] = rightPart.x;
	maxwellEquationRightPart[i][j][k][1] = rightPart.y;
	maxwellEquationRightPart[i][j][k][2] = rightPart.z;

	//printf("right part = %g %g %g\n", rightPart.x, rightPart.y, rightPart.z);

	//maxwellEquationRightPart[i][3] = rightPart2.x;
	//maxwellEquationRightPart[i][4] = rightPart2.y;
	//maxwellEquationRightPart[i][5] = rightPart2.z;
}

void Simulation::evaluateMagneticField() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if (solverType == BUNEMAN) {
		evaluateBunemanMagneticField();
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					rotE[i][j][k] = evaluateRotTempE(i, j, k);
					if (isInResistiveLayer(i, j, k)) {
						newBfield[i][j][k] = (Bfield[i][j][k] + B0 * (fakeCondactivity * deltaT) - (rotE[i][j][k]) * (
							speed_of_light_normalized * speed_of_light_correction * deltaT)) / (1 + fakeCondactivity * deltaT);
					} else {
						newBfield[i][j][k] = Bfield[i][j][k] - (rotE[i][j][k]) * (speed_of_light_normalized * speed_of_light_correction *
							deltaT);
					}
					if ((boundaryConditionTypeX != PERIODIC) && (cartCoord[0] == cartDim[0] - 1) && (i >= xnumberAdded - 1 -
						additionalBinNumber)) {
						newBfield[i][j][k] = B0;
					}
				}
			}
		}

		if ((boundaryConditionTypeX != PERIODIC) && (cartCoord[0] == 0)) {
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						newBfield[i][j][k] = newBfield[additionalBinNumber + 1][j][k];
					}
				}
			}
		}
		/*double kw = (2 * pi / xsizeGeneral);
		for (int i = 0; i <= xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					newBfield[i][j][k].z = B0.norm() * sin(kw * middleXgrid[i] - kw*speed_of_light_normalized*(time + deltaT));
				}
			}
		}*/
	}
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("evaluating magnetic field time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

bool Simulation::isInResistiveLayer(int i, int j, int k) {
	//return false;
	if (boundaryConditionTypeX == PERIODIC) {
		return false;
	}
	if (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) {
		return false;
	}

	if (i >= xnumberAdded - 1 - additionalBinNumber - resistiveLayerWidth && i < xnumberAdded - 1 - additionalBinNumber) {
		return cartCoord[0] == cartDim[0] - 1;
	}

	if (i > 1 + additionalBinNumber && i <= 1 + additionalBinNumber + resistiveLayerWidth) {
		if (boundaryConditionTypeX != FREE_BOTH) {
			return false;
		}
		return cartCoord[0] == 0;
	}

	return false;
}


Vector3d Simulation::evaluateRotB(int i, int j, int k) {
	Vector3d BrightX;
	Vector3d BleftX;
	Vector3d BrightY;
	Vector3d BleftY;
	Vector3d BrightZ;
	Vector3d BleftZ;

	int curI = i;

	int prevI = i - 1;

	int curJ = j;

	int prevJ = j - 1;

	int prevK = k - 1;

	int curK = k;

	//printf("BrightX\n");
	if (curI == 0) {
		printf("curI = 0 in evaluateRotB\n");
		fflush(stdout);
		MPI_Finalize();
		exit(0);
	}
	BrightX = (Bfield[curI][curJ][curK] + Bfield[curI][prevJ][curK] + Bfield[curI][curJ][prevK] + Bfield[curI][prevJ][prevK
	]) / 4.0;
	BleftX = (Bfield[prevI][curJ][curK] + Bfield[prevI][prevJ][curK] + Bfield[prevI][curJ][prevK] + Bfield[prevI][prevJ][
		prevK]) / 4.0;

	BrightY = (Bfield[curI][curJ][curK] + Bfield[prevI][curJ][curK] + Bfield[curI][curJ][prevK] + Bfield[prevI][curJ][prevK
	]) / 4.0;
	BleftY = (Bfield[curI][prevJ][curK] + Bfield[prevI][prevJ][curK] + Bfield[curI][prevJ][prevK] + Bfield[prevI][prevJ][
		prevK]) / 4.0;

	BrightZ = (Bfield[curI][curJ][curK] + Bfield[prevI][curJ][curK] + Bfield[curI][prevJ][curK] + Bfield[prevI][prevJ][curK
	]) / 4.0;
	BleftZ = (Bfield[curI][curJ][prevK] + Bfield[prevI][curJ][prevK] + Bfield[curI][prevJ][prevK] + Bfield[prevI][prevJ][
		prevK]) / 4.0;

	double x = 0;
	double y = 0;
	double z = 0;

	x = (BrightY.z - BleftY.z) / deltaY - (BrightZ.y - BleftZ.y) / deltaZ;
	y = (BrightZ.x - BleftZ.x) / deltaZ - (BrightX.z - BleftX.z) / deltaX;
	z = (BrightX.y - BleftX.y) / deltaX - (BrightY.x - BleftY.x) / deltaY;

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotEgeneral(Vector3d*** E, int i, int j, int k) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "i = %d < 0 in evaluateRotTempE\n", i);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}

		if (i >= xnumberAdded) {
			printf("x >= xnumberAdded\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "i = %d >= xnumber = %d in evaluateRotTempE\n", i, xnumberAdded);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}

		if (j < 0) {
			printf("j < 0\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "j = %d < 0 in evaluateRotTempE\n", j);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}

		if (j >= ynumberAdded) {
			printf("y >= ynumberAdded\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "j = %d >= ynumber = %d in evaluateRotTempE\n", j, ynumberAdded);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}

		if (k < 0) {
			printf("k < 0\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "k = %d < 0 in evaluateRotTempE\n", k);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}

		if (k >= znumberAdded) {
			printf("z >= znumberAdded\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "k = %d >= znumber = %d in evaluateRotTempE\n", k, znumber);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = (E[i + 1][j][k] + E[i + 1][j + 1][k] + E[i + 1][j][k + 1] + E[i + 1][j + 1][k + 1]) / 4.0;
	Vector3d EleftX = (E[i][j][k] + E[i][j + 1][k] + E[i][j][k + 1] + E[i][j + 1][k + 1]) / 4.0;

	Vector3d ErightY = (E[i][j + 1][k] + E[i + 1][j + 1][k] + E[i][j + 1][k + 1] + E[i + 1][j + 1][k + 1]) / 4.0;
	Vector3d EleftY = (E[i][j][k] + E[i + 1][j][k] + E[i][j][k + 1] + E[i + 1][j][k + 1]) / 4.0;
	if (ynumberGeneral == 1) {
		ErightY = Vector3d(0, 0, 0);
		EleftY = Vector3d(0, 0, 0);
	}

	Vector3d ErightZ = (E[i][j][k + 1] + E[i + 1][j][k + 1] + E[i][j + 1][k + 1] + E[i + 1][j + 1][k + 1]) / 4.0;
	Vector3d EleftZ = (E[i][j][k] + E[i + 1][j][k] + E[i][j + 1][k] + E[i + 1][j + 1][k]) / 4.0;
	if (znumberGeneral == 1) {
		ErightZ = Vector3d(0, 0, 0);
		EleftZ = Vector3d(0, 0, 0);
	}

	x = (ErightY.z - EleftY.z) / deltaY - (ErightZ.y - EleftZ.y) / deltaZ;
	//x = - (ErightZ.y - EleftZ.y) / deltaZ;
	y = (ErightZ.x - EleftZ.x) / deltaZ - (ErightX.z - EleftX.z) / deltaX;
	z = (ErightX.y - EleftX.y) / deltaX - (ErightY.x - EleftY.x) / deltaY;
	//z = (ErightX.y - EleftX.y) / deltaX;

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotTempE(int i, int j, int k) {
	return evaluateRotEgeneral(tempEfield, i, j, k);
}

Vector3d Simulation::evaluateRotE(int i, int j, int k) {
	return evaluateRotEgeneral(Efield, i, j, k);
}

Vector3d Simulation::evaluateRotNewE(int i, int j, int k) {
	return evaluateRotEgeneral(newEfield, i, j, k);
}

double Simulation::evaluateDivEgeneral(Vector3d*** E, int i, int j, int k) {
	double ErightX = (E[i + 1][j][k].x + E[i + 1][j + 1][k].x + E[i + 1][j][k + 1].x + E[i + 1][j + 1][k + 1].x) / 4.0;
	double EleftX = (E[i][j][k].x + E[i][j + 1][k].x + E[i][j][k + 1].x + E[i][j + 1][k + 1].x) / 4.0;

	double ErightY = (E[i][j + 1][k].y + E[i + 1][j + 1][k].y + E[i][j + 1][k + 1].y + E[i + 1][j + 1][k + 1].y) / 4.0;
	double EleftY = (E[i][j][k].y + E[i + 1][j][k].y + E[i][j][k + 1].y + E[i + 1][j][k + 1].y) / 4.0;

	double ErightZ = (E[i][j][k + 1].z + E[i + 1][j][k + 1].z + E[i][j + 1][k + 1].z + E[i + 1][j + 1][k + 1].z) / 4.0;
	double EleftZ = (E[i][j][k].z + E[i + 1][j][k].z + E[i][j + 1][k].z + E[i + 1][j + 1][k].z) / 4.0;

	return ((ErightX - EleftX) / deltaX) + ((ErightY - EleftY) / deltaY) + ((ErightZ - EleftZ) / deltaZ);
}

double Simulation::evaluateDivE(int i, int j, int k) {
	return evaluateDivEgeneral(Efield, i, j, k);
}

double Simulation::evaluateDivCleaningE(int i, int j, int k) {
	double ErightX = (divergenceCleaningField[i + 1][j][k][0] + divergenceCleaningField[i + 1][j + 1][k][0] +
		divergenceCleaningField[i + 1][j][k + 1][0] + divergenceCleaningField[i + 1][j + 1][k + 1][0]) / 4.0;
	double EleftX = (divergenceCleaningField[i][j][k][0] + divergenceCleaningField[i][j + 1][k][0] +
		divergenceCleaningField[i][j][k + 1][0] + divergenceCleaningField[i][j + 1][k + 1][0]) / 4.0;

	double ErightY = (divergenceCleaningField[i][j + 1][k][1] + divergenceCleaningField[i + 1][j + 1][k][1] +
		divergenceCleaningField[i][j + 1][k + 1][1] + divergenceCleaningField[i + 1][j + 1][k + 1][1]) / 4.0;
	double EleftY = (divergenceCleaningField[i][j][k][1] + divergenceCleaningField[i + 1][j][k][1] +
		divergenceCleaningField[i][j][k + 1][1] + divergenceCleaningField[i + 1][j][k + 1][1]) / 4.0;

	double ErightZ = (divergenceCleaningField[i][j][k + 1][2] + divergenceCleaningField[i + 1][j][k + 1][2] +
		divergenceCleaningField[i][j + 1][k + 1][2] + divergenceCleaningField[i + 1][j + 1][k + 1][2]) / 4.0;
	double EleftZ = (divergenceCleaningField[i][j][k][2] + divergenceCleaningField[i + 1][j][k][2] +
		divergenceCleaningField[i][j + 1][k][2] + divergenceCleaningField[i + 1][j + 1][k][2]) / 4.0;

	return ((ErightX - EleftX) / deltaX) + ((ErightY - EleftY) / deltaY) + ((ErightZ - EleftZ) / deltaZ);
}

double Simulation::evaluateDivTempE(int i, int j, int k) {
	return evaluateDivEgeneral(tempEfield, i, j, k);
}

double Simulation::evaluateDivNewE(int i, int j, int k) {
	return evaluateDivEgeneral(newEfield, i, j, k);
}

double Simulation::evaluateDivB(int i, int j, int k) {
	return evaluateDivBgeneral(Bfield, i, j, k);
}

double Simulation::evaluateDivNewB(int i, int j, int k) {
	return evaluateDivBgeneral(newBfield, i, j, k);
}

double Simulation::evaluateDivBgeneral(Vector3d*** B, int i, int j, int k) {
	double BrightX = (B[i][j - 1][k - 1].x + B[i][j][k - 1].x + B[i][j - 1][k].x + B[i][j][k].x) / 4.0;
	double BleftX = (B[i - 1][j - 1][k - 1].x + B[i - 1][j][k - 1].x + B[i - 1][j - 1][k].x + B[i - 1][j][k].x) / 4.0;

	double BrightY = (B[i - 1][j][k - 1].y + B[i][j][k - 1].y + B[i - 1][j][k].y + B[i][j][k].y) / 4.0;
	double BleftY = (B[i - 1][j - 1][k - 1].y + B[i][j - 1][k - 1].y + B[i - 1][j - 1][k].y + B[i][j - 1][k].y) / 4.0;

	double BrightZ = (B[i - 1][j - 1][k].z + B[i][j - 1][k].z + B[i - 1][j][k].z + B[i][j][k].z) / 4.0;
	double BleftZ = (B[i - 1][j - 1][k - 1].z + B[i][j - 1][k - 1].z + B[i - 1][j][k - 1].z + B[i][j][k - 1].z) / 4.0;

	return ((BrightX - BleftX) / deltaX) + ((BrightY - BleftY) / deltaY) + ((BrightZ - BleftZ) / deltaZ);
}

double Simulation::
evaluateDivBunemanEgeneral(double*** fieldX, double*** fieldY, double*** fieldZ, int i, int j, int k) {
	return ((fieldX[i][j][k] - fieldX[i - 1][j][k]) / deltaX) + ((fieldY[i][j][k] - fieldY[i][j - 1][k]) / deltaY) + ((
		fieldZ[i][j][k] - fieldZ[i][j][k - 1]) / deltaZ);
}

double Simulation::evaluateDivBunemanNewE(int i, int j, int k) {
	return evaluateDivBunemanEgeneral(bunemanNewEx, bunemanNewEy, bunemanNewEz, i, j, k);
}

double Simulation::evaluateDivBunemanE(int i, int j, int k) {
	return evaluateDivBunemanEgeneral(bunemanEx, bunemanEy, bunemanEz, i, j, k);
}

double Simulation::
evaluateDivBunemanBgeneral(double*** fieldX, double*** fieldY, double*** fieldZ, int i, int j, int k) {
	return ((fieldX[i + 1][j][k] - fieldX[i][j][k]) / deltaX) + ((fieldY[i][j + 1][k] - fieldY[i][j][k]) / deltaY) + ((
		fieldZ[i][j][k + 1] - fieldZ[i][j][k]) / deltaZ);
}

double Simulation::evaluateDivBunemanNewB(int i, int j, int k) {
	return evaluateDivBunemanBgeneral(bunemanNewBx, bunemanNewBy, bunemanNewBz, i, j, k);
}

double Simulation::evaluateDivBunemanB(int i, int j, int k) {
	return evaluateDivBunemanBgeneral(bunemanBx, bunemanBy, bunemanBz, i, j, k);
}

double Simulation::evaluateDivFlux(int i, int j, int k) {
	if (debugMode) {
		if (i < -additionalBinNumber) {
			printf("i < 0\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "i = %d < 0 in evaluateDivFlux\n", i);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}

		if (i > xnumberAdded - 1) {
			printf("x >= xnumber\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "i = %d >= xnumber = %d in evaluateDivFlux\n", i, xnumberAdded - 1);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}

		if (j < 0) {
			printf("i < 0\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "j = %d < 0 in evaluateDivFlux\n", j);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}

		if (j >= ynumberAdded) {
			printf("y >= ynumberAdded\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "j = %d >= ynumber = %d in evaluateDivFlux\n", j, ynumberAdded);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}

		if (k < 0) {
			printf("k < 0\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "k = %d < 0 in evaluateDivFlux\n", k);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}

		if (k >= znumberAdded) {
			printf("z >= znumberAdded\n");
			fflush(stdout);
			errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "k = %d >= znumber = %d in evaluateDivFlux\n", k, znumber);
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}
	}


	double rightFluxX =
	(getElectricFlux(i + 1, j, k).x + getElectricFlux(i + 1, j + 1, k).x + getElectricFlux(i + 1, j, k + 1).x +
		getElectricFlux(i + 1, j + 1, k + 1).x) * 0.25;
	double leftFluxX = (getElectricFlux(i, j, k).x + getElectricFlux(i, j + 1, k).x + getElectricFlux(i, j, k + 1).x +
		getElectricFlux(i, j + 1, k + 1).x) * 0.25;

	double rightFluxY = 0;
	if (ynumberGeneral > 1) {
		rightFluxY = (getElectricFlux(i, j + 1, k).y + getElectricFlux(i + 1, j + 1, k).y + getElectricFlux(i, j + 1, k + 1).y
			+
			getElectricFlux(i + 1, j + 1, k + 1).y) * 0.25;
	}
	double leftFluxY = 0;
	if (ynumberGeneral > 1) {
		leftFluxY = (getElectricFlux(i, j, k).y + getElectricFlux(i + 1, j, k).y + getElectricFlux(i, j, k + 1).y +
			getElectricFlux(i + 1, j, k + 1).y) * 0.25;
	}

	double rightFluxZ = 0;
	if (znumberGeneral > 1) {
		rightFluxZ = (getElectricFlux(i + 1, j, k + 1).z + getElectricFlux(i, j, k + 1).z + getElectricFlux(i, j + 1, k + 1).z
			+
			getElectricFlux(i + 1, j + 1, k + 1).z) * 0.25;
	}
	double leftFluxZ = 0;
	if (znumberGeneral > 1) {
		leftFluxZ = (getElectricFlux(i + 1, j, k).z + getElectricFlux(i, j, k).z + getElectricFlux(i, j + 1, k).z +
			getElectricFlux(i + 1, j + 1, k).z) * 0.25;
	}

	return ((rightFluxX - leftFluxX) / deltaX) + ((rightFluxY - leftFluxY) / deltaY) +
		((rightFluxZ - leftFluxZ) / deltaZ);
}

Vector3d Simulation::getElectricFlux(int i, int j, int k) {
	return electricFlux[i][j][k];
}

Vector3d Simulation::evaluateGradDensity(int i, int j, int k) {
	int prevI = i - 1;

	int prevJ = j - 1;

	int prevK = k - 1;

	double densityRightX = (chargeDensityHat[i][j][k] + chargeDensityHat[i][prevJ][k] + chargeDensityHat[i][j][prevK] +
		chargeDensityHat[i][prevJ][prevK]) * 0.25;
	double densityLeftX =
	(chargeDensityHat[prevI][j][k] + chargeDensityHat[prevI][prevJ][k] + chargeDensityHat[prevI][j][prevK] +
		chargeDensityHat[prevI][prevJ][prevK]) * 0.25;

	double densityRightY = (chargeDensityHat[i][j][k] + chargeDensityHat[prevI][j][k] + chargeDensityHat[i][j][prevK] +
		chargeDensityHat[prevI][j][prevK]) * 0.25;
	double densityLeftY =
	(chargeDensityHat[i][prevJ][k] + chargeDensityHat[prevI][prevJ][k] + chargeDensityHat[i][prevJ][prevK] +
		chargeDensityHat[prevI][prevJ][prevK]) * 0.25;

	double densityRightZ = (chargeDensityHat[i][j][k] + chargeDensityHat[i][prevJ][k] + chargeDensityHat[prevI][j][k] +
		chargeDensityHat[prevI][prevJ][k]) * 0.25;
	double densityLeftZ =
	(chargeDensityHat[i][j][prevK] + chargeDensityHat[i][prevJ][prevK] + chargeDensityHat[prevI][j][prevK] +
		chargeDensityHat[prevI][prevJ][prevK]) * 0.25;

	double x = (densityRightX - densityLeftX) / deltaX;
	double y = (densityRightY - densityLeftY) / deltaY;
	double z = (densityRightZ - densityLeftZ) / deltaZ;

	return Vector3d(x, y, z);
}

void Simulation::createFakeEquation(int i, int j, int k) {
	Vector3d rightPart = Vector3d(0, 0, 0);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(1.0, i, j, k, 1));
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(1.0, i, j, k, 2));


	maxwellEquationRightPart[i][j][k][0] = rightPart.x;
	maxwellEquationRightPart[i][j][k][1] = rightPart.y;
	maxwellEquationRightPart[i][j][k][2] = rightPart.z;
}

double Simulation::evaluateBunemanRotEx(int i, int j, int k) {
	int nextJ = j + 1;
	int nextK = k + 1;
	double result = ((bunemanEz[i][nextJ][k] - bunemanEz[i][j][k]) / deltaY) - ((bunemanEy[i][j][nextK] - bunemanEy[i][j][k
	]) / deltaZ);
	return result;
}

double Simulation::evaluateBunemanRotEy(int i, int j, int k) {
	int nextI = i + 1;
	int nextK = k + 1;
	double result = ((bunemanEx[i][j][nextK] - bunemanEx[i][j][k]) / deltaZ) - ((bunemanEz[nextI][j][k] - bunemanEz[i][j][k
	]) / deltaZ);
	return result;
}

double Simulation::evaluateBunemanRotEz(int i, int j, int k) {
	int nextI = i + 1;
	int nextJ = j + 1;
	double result = ((bunemanEy[nextI][j][k] - bunemanEy[i][j][k]) / deltaX) - ((bunemanEx[i][nextJ][k] - bunemanEx[i][j][k
	]) / deltaY);
	return result;
}

double Simulation::evaluateBunemanRotBx(int i, int j, int k) {
	int prevJ = j - 1;
	int prevK = k - 1;
	double result = ((bunemanBz[i][j][k] - bunemanBz[i][prevJ][k]) / deltaY) - ((bunemanBy[i][j][k] - bunemanBy[i][j][prevK
	]) / deltaZ);
	return result;
}

double Simulation::evaluateBunemanRotBy(int i, int j, int k) {
	int prevI = i - 1;
	int prevK = k - 1;
	double result = ((bunemanBx[i][j][k] - bunemanBx[i][j][prevK]) / deltaZ) - ((bunemanBz[i][j][k] - bunemanBz[prevI][j][k
	]) / deltaX);
	return result;
}

double Simulation::evaluateBunemanRotBz(int i, int j, int k) {
	int prevJ = j - 1;
	int prevI = i - 1;
	double result = ((bunemanBy[i][j][k] - bunemanBy[prevI][j][k]) / deltaX) - ((bunemanBx[i][j][k] - bunemanBx[i][prevJ][k
	]) / deltaY);
	return result;
}

void Simulation::evaluateBunemanElectricField() {
	evaluateBunemanElectricFieldX();
	evaluateBunemanElectricFieldY();
	evaluateBunemanElectricFieldZ();
}

void Simulation::evaluateBunemanElectricFieldX() {
	double omega2 = 0;
	for (int i = 0; i < typesNumber; ++i) {
		omega2 += 4 * pi * types[i].concentration * types[i].charge * types[i].charge / types[i].mass;
	}

	double gamma = 1.0 / sqrt(1 - V0.scalarMult(V0) / speed_of_light_normalized_sqr);
	double omega = sqrt(omega2 / (gamma) * gamma * gamma);

	for (int i = 1 + additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
		for (int j = 1 + additionalBinNumber; j <= ynumberAdded - 1 - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k <= znumberAdded - 1 - additionalBinNumber; ++k) {
				bunemanNewEx[i][j][k] = bunemanEx[i][j][k] + (speed_of_light_normalized * evaluateBunemanRotBx(i, j, k) - 4 * pi *
					bunemanJx[i][j][k]) * deltaT;
				//bunemanNewEx[i][j][k] = E0.x*cos(omega*(time + deltaT));
				//bunemanNewEx[i][j][k] = 0;
			}
		}
	}
}

void Simulation::evaluateBunemanElectricFieldY() {
	for (int i = 1 + additionalBinNumber; i <= xnumberAdded - 1 - additionalBinNumber; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k <= znumberAdded - 1 - additionalBinNumber; ++k) {
				bunemanNewEy[i][j][k] = bunemanEy[i][j][k] + (speed_of_light_normalized * evaluateBunemanRotBy(i, j, k) - 4 * pi *
					bunemanJy[i][j][k]) * deltaT;
				if (cartCoord[0] == 0 && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && i <= 1 + additionalBinNumber) {
					bunemanNewEy[i][j][k] = 0;
				}
				//bunemanNewEy[i][j][k] = 0;
			}
		}
	}
}

void Simulation::evaluateBunemanElectricFieldZ() {
	for (int i = 1 + additionalBinNumber; i <= xnumberAdded - 1 - additionalBinNumber; ++i) {
		for (int j = 1 + additionalBinNumber; j <= ynumberAdded - 1 - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
				bunemanNewEz[i][j][k] = bunemanEz[i][j][k] + (speed_of_light_normalized * evaluateBunemanRotBz(i, j, k) - 4 * pi *
					bunemanJz[i][j][k]) * deltaT;
				if (cartCoord[0] == 0 && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && i <= 1 + additionalBinNumber) {
					bunemanNewEz[i][j][k] = 0;
				}
				//bunemanNewEz[i][j][k] = 0;
			}
		}
	}
}

void Simulation::evaluateBunemanMagneticField() {
	evaluateBunemanMagneticFieldX();
	evaluateBunemanMagneticFieldY();
	evaluateBunemanMagneticFieldZ();
}

void Simulation::evaluateBunemanMagneticFieldX() {
	for (int i = 1 + additionalBinNumber; i <= xnumberAdded - 1 - additionalBinNumber; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
				bunemanNewBx[i][j][k] = bunemanBx[i][j][k] - speed_of_light_normalized * evaluateBunemanRotEx(i, j, k) * deltaT;
			}
		}
	}
}

void Simulation::evaluateBunemanMagneticFieldY() {
	for (int i = 1 + additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
		for (int j = 1 + additionalBinNumber; j <= ynumberAdded - 1 - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
				bunemanNewBy[i][j][k] = bunemanBy[i][j][k] - speed_of_light_normalized * evaluateBunemanRotEy(i, j, k) * deltaT;
			}
		}
	}
}

void Simulation::evaluateBunemanMagneticFieldZ() {
	for (int i = 1 + additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k <= znumberAdded - 1 - additionalBinNumber; ++k) {
				bunemanNewBz[i][j][k] = bunemanBz[i][j][k] - speed_of_light_normalized * evaluateBunemanRotEz(i, j, k) * deltaT;
			}
		}
	}
}

void Simulation::updateBunemanFields() {
	updateBunemanElectricField();
	updateBunemanMagneticField();
	/*for(int i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				Bfield[i][j][k] = getBunemanMagneticField(i, j, k);
				Efield[i][j][k] = getBunemanElectricField(i, j, k);
			}
		}
	}*/

	/*for(int i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				Efield[i][j][k] = getBunemanElectricField(i, j, k);
			}
		}
	}*/
}

void Simulation::updateBunemanElectricField() {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanEx[i][j][k] = bunemanNewEx[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanEy[i][j][k] = bunemanNewEy[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				bunemanEz[i][j][k] = bunemanNewEz[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Efield[i][j][k] = getBunemanElectricField(i, j, k);
			}
		}
	}
}

void Simulation::updateBunemanMagneticField() {
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				bunemanBx[i][j][k] = bunemanNewBx[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				bunemanBy[i][j][k] = bunemanNewBy[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanBz[i][j][k] = bunemanNewBz[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k] = getBunemanMagneticField(i, j, k);
			}
		}
	}
}

Vector3d Simulation::getBunemanElectricField(int i, int j, int k) {
	return Vector3d(bunemanEx[i][j][k], bunemanEy[i][j][k], bunemanEz[i][j][k]);
}

Vector3d Simulation::getBunemanMagneticField(int i, int j, int k) {
	return Vector3d(bunemanBx[i][j][k], bunemanBy[i][j][k], bunemanBz[i][j][k]);
}

Vector3d Simulation::averageFieldXY(Vector3d*** field, int k) {
	double tempField[3];
	double sumField[3];
	tempField[0] = 0;
	tempField[1] = 0;
	tempField[2] = 0;

	for (int i = 1 + additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
			tempField[0] += field[i][j][k].x;
			tempField[1] += field[i][j][k].y;
			tempField[2] += field[i][j][k].z;
		}
	}

	MPI_Allreduce(tempField, sumField, 3, MPI_DOUBLE, MPI_SUM, cartCommXY);

	Vector3d result;

	result.x = sumField[0] / (xnumberGeneral * ynumberGeneral);
	result.y = sumField[1] / (xnumberGeneral * ynumberGeneral);
	result.z = sumField[2] / (xnumberGeneral * ynumberGeneral);

	return result;
}

Vector3d Simulation::averageFieldXZ(Vector3d*** field, int j) {
	double tempField[3];
	double sumField[3];
	tempField[0] = 0;
	tempField[1] = 0;
	tempField[2] = 0;

	for (int i = 1 + additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
		for (int k = 1 + additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
			tempField[0] += field[i][j][k].x;
			tempField[1] += field[i][j][k].y;
			tempField[2] += field[i][j][k].z;
		}
	}

	MPI_Allreduce(tempField, sumField, 3, MPI_DOUBLE, MPI_SUM, cartCommXZ);

	Vector3d result;

	result.x = sumField[0] / (xnumberGeneral * znumberGeneral);
	result.y = sumField[1] / (xnumberGeneral * znumberGeneral);
	result.z = sumField[2] / (xnumberGeneral * znumberGeneral);

	return result;
}

Vector3d Simulation::averageFieldYZ(Vector3d*** field, int i) {
	double tempField[3];
	double sumField[3];
	tempField[0] = 0;
	tempField[1] = 0;
	tempField[2] = 0;

	for (int j = 1 + additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
		for (int k = 1 + additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
			tempField[0] += field[i][j][k].x;
			tempField[1] += field[i][j][k].y;
			tempField[2] += field[i][j][k].z;
		}
	}

	MPI_Allreduce(tempField, sumField, 3, MPI_DOUBLE, MPI_SUM, cartCommYZ);

	Vector3d result;

	result.x = sumField[0] / (ynumberGeneral * znumberGeneral);
	result.y = sumField[1] / (ynumberGeneral * znumberGeneral);
	result.z = sumField[2] / (ynumberGeneral * znumberGeneral);

	return result;
}

double Simulation::averageConcentrationYZ(double*** concentration, int i) {
	double tempField[1];
	double sumField[1];
	tempField[0] = 0;

	for (int j = 1 + additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
		for (int k = 1 + additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
			tempField[0] += concentration[i][j][k];
		}
	}

	MPI_Allreduce(tempField, sumField, 1, MPI_DOUBLE, MPI_SUM, cartCommYZ);

	double result;

	result = sumField[0] / (ynumberGeneral * znumberGeneral);

	return result;
}
