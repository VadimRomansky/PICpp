#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include <time.h>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "mpi_util.h"
#include "util.h"
#include "output.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "particle.h"
#include "random.h"
#include "simulation.h"
#include "paths.h"

void Simulation::sumNodeMassMatrixParametersX() {
	if (cartDim[0] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum node vector parameters x\n");
		int n = 9 * (2 * splineOrder + 3) * (2 * splineOrder + 3) * (2 * splineOrder + 3);
		double* inBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * n];
		double* outBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * n];
		double* inBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * n];
		double* outBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * n];

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("sending left flux sum node vector parameters x rank = %d\n", rank);
		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendNodeMassMatrixParametersToLeftReceiveFromRight(massMatrix, outBufferLeft, tempNodeMassMatrixParameterRight,
			                                                   inBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
			                                                   additionalBinNumber, cartComm, rank, leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			receiveNodeMassMatrixParametersRight(tempNodeMassMatrixParameterRight, inBufferRight, xnumberAdded, ynumberAdded,
			                                     znumberAdded, additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			sendNodeMassMatrixParametersLeft(massMatrix, outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
			                                 additionalBinNumber, cartComm, rank, leftRank);
		}

		//MPI_Barrier(cartComm);

		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendNodeMassMatrixParametersToRightReceiveFromLeft(massMatrix, outBufferRight, tempNodeMassMatrixParameterLeft,
			                                                   inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
			                                                   additionalBinNumber, cartComm, rank, leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			sendNodeMassMatrixParametersRight(massMatrix, outBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
			                                  additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			receiveNodeMassMatrixParametersLeft(tempNodeMassMatrixParameterLeft, inBufferLeft, xnumberAdded, ynumberAdded,
			                                    znumberAdded, additionalBinNumber, cartComm, rank, leftRank);
		}


		sumTempNodeMassMatrixParametersX(massMatrix);


		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters x rank = %d\n", rank);
		delete[] inBufferLeft;
		delete[] inBufferRight;
		delete[] outBufferLeft;
		delete[] outBufferRight;

	} else {
		if (boundaryConditionTypeX == PERIODIC) {
			for (int j = 0; j <= ynumberAdded; ++j) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int i = 0; i <= 2 * additionalBinNumber + 2; ++i) {
						for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
							for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
								for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
									for (int l = 0; l < 3; ++l) {
										for (int m = 0; m < 3; ++m) {
											massMatrix[xnumberAdded - i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] += massMatrix[2 + 2 *
												additionalBinNumber - i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
											massMatrix[2 + 2 * additionalBinNumber - i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = massMatrix[
												xnumberAdded - i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumNodeMassMatrixParametersY() {
	if (cartDim[1] > 1) {
		int n = 9 * (2 * splineOrder + 3) * (2 * splineOrder + 3) * (2 * splineOrder + 3);
		double* inBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * n];
		double* outBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * n];
		double* inBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * n];
		double* outBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * n];

		if ((verbosity > 2)) printf("sending left flux sum node vector parameters y rank = %d\n", rank);

		if ((boundaryConditionTypeY == PERIODIC) || (cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1)) {
			sendNodeMassMatrixParametersToFrontReceiveFromBack(massMatrix, outBufferFront, tempNodeMassMatrixParameterBack,
			                                                   inBufferBack, xnumberAdded, ynumberAdded, znumberAdded,
			                                                   additionalBinNumber, cartComm, rank, frontRank, backRank);
		} else if (cartCoord[1] == 0) {
			receiveNodeMassMatrixParametersBack(tempNodeMassMatrixParameterBack, inBufferBack, xnumberAdded, ynumberAdded,
			                                    znumberAdded, additionalBinNumber, cartComm, rank, backRank);
		} else if (cartCoord[1] == cartDim[1] - 1) {
			sendNodeMassMatrixParametersFront(massMatrix, outBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
			                                  additionalBinNumber, cartComm, rank, frontRank);
		}

		//MPI_Barrier(cartComm);

		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1)) {
			sendNodeMassMatrixParametersToBackReceiveFromFront(massMatrix, outBufferBack, tempNodeMassMatrixParameterFront,
			                                                   inBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
			                                                   additionalBinNumber, cartComm, rank, frontRank, backRank);
		} else if (cartCoord[1] == 0) {
			sendNodeMassMatrixParametersBack(massMatrix, outBufferBack, xnumberAdded, ynumberAdded, znumberAdded,
			                                 additionalBinNumber, cartComm, rank, backRank);
		} else if (cartCoord[1] == cartDim[1] - 1) {
			receiveNodeMassMatrixParametersFront(tempNodeMassMatrixParameterFront, inBufferFront, xnumberAdded, ynumberAdded,
			                                     znumberAdded, additionalBinNumber, cartComm, rank, frontRank);
		}

		sumTempNodeMassMatrixParametersY(massMatrix);

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters y rank = %d\n", rank);

		delete[] inBufferFront;
		delete[] inBufferBack;
		delete[] outBufferBack;
		delete[] outBufferFront;
	} else {
		if (boundaryConditionTypeY == PERIODIC) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= 2 * additionalBinNumber + 2; ++j) {
						for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
							for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
								for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
									for (int l = 0; l < 3; ++l) {
										for (int m = 0; m < 3; ++m) {
											massMatrix[i][ynumberAdded - j][k].matrix[tempI][tempJ][tempK].matrix[l][m] += massMatrix[i][2 + 2 *
												additionalBinNumber - j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
											massMatrix[i][2 + 2 * additionalBinNumber - j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = massMatrix[i][
												ynumberAdded - j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
										}
									}
								}
							}
						}
					}
				}
			}
		}

		if (ynumberGeneral == 1) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= ynumberAdded; ++j) {
						for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
							for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
								for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
									for (int l = 0; l < 3; ++l) {
										for (int m = 0; m < 3; ++m) {
											massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = massMatrix[i][1 + additionalBinNumber][k].
												matrix[tempI][tempJ][tempK].matrix[l][m];
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumNodeMassMatrixParametersZ() {
	if (cartDim[2] > 1) {
		int n = 9 * (2 * splineOrder + 3) * (2 * splineOrder + 3) * (2 * splineOrder + 3);
		double* inBufferTop = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1) * n];
		double* outBufferTop = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1) * n];
		double* inBufferBottom = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1) * n];
		double* outBufferBottom = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1) * n];

		if ((verbosity > 2)) printf("sending left flux sum node vector parameters z rank = %d\n", rank);

		if ((boundaryConditionTypeZ == PERIODIC) || (cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1)) {
			sendNodeMassMatrixParametersToBottomReceiveFromTop(massMatrix, outBufferBottom, tempNodeMassMatrixParameterTop,
			                                                   inBufferTop, xnumberAdded, ynumberAdded, znumberAdded,
			                                                   additionalBinNumber, cartComm, rank, bottomRank, topRank);
		} else if (cartCoord[2] == 0) {
			receiveNodeMassMatrixParametersTop(tempNodeMassMatrixParameterTop, inBufferTop, xnumberAdded, ynumberAdded,
			                                   znumberAdded, additionalBinNumber, cartComm, rank, topRank);
		} else if (cartCoord[2] == cartDim[2] - 1) {
			sendNodeMassMatrixParametersBottom(massMatrix, outBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
			                                   additionalBinNumber, cartComm, rank, bottomRank);
		}

		//MPI_Barrier(cartComm);

		if ((boundaryConditionTypeZ == PERIODIC) || (cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1)) {
			sendNodeMassMatrixParametersToTopReceiveFromBottom(massMatrix, outBufferTop, tempNodeMassMatrixParameterBottom,
			                                                   inBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
			                                                   additionalBinNumber, cartComm, rank, bottomRank, topRank);
		} else if (cartCoord[2] == 0) {
			sendNodeMassMatrixParametersTop(massMatrix, outBufferTop, xnumberAdded, ynumberAdded, znumberAdded,
			                                additionalBinNumber, cartComm, rank, topRank);
		} else if (cartCoord[2] == cartDim[2] - 1) {
			receiveNodeMassMatrixParametersBottom(tempNodeMassMatrixParameterBottom, inBufferBottom, xnumberAdded, ynumberAdded,
			                                      znumberAdded, additionalBinNumber, cartComm, rank, bottomRank);
		}

		sumTempNodeMassMatrixParametersZ(massMatrix);


		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters z rank = %d\n", rank);

		delete[] inBufferTop;
		delete[] inBufferBottom;
		delete[] outBufferBottom;
		delete[] outBufferTop;

	} else {
		if (boundaryConditionTypeZ == PERIODIC) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int j = 0; j <= ynumberAdded; ++j) {
					for (int k = 0; k <= 2 * additionalBinNumber + 2; ++k) {
						for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
							for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
								for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
									for (int l = 0; l < 3; ++l) {
										for (int m = 0; m < 3; ++m) {
											massMatrix[i][j][znumberAdded - k].matrix[tempI][tempJ][tempK].matrix[l][m] += massMatrix[i][j][2 + 2 *
												additionalBinNumber - k].matrix[tempI][tempJ][tempK].matrix[l][m];
											massMatrix[i][j][2 + 2 * additionalBinNumber - k].matrix[tempI][tempJ][tempK].matrix[l][m] = massMatrix[i][j]
											[
												znumberAdded - k].matrix[tempI][tempJ][tempK].matrix[l][m];
										}
									}
								}
							}
						}
					}
				}
			}
		}
		if (znumberGeneral == 1) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= ynumberAdded; ++j) {
						for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
							for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
								for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
									for (int l = 0; l < 3; ++l) {
										for (int m = 0; m < 3; ++m) {
											massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = massMatrix[i][j][1 + additionalBinNumber].
												matrix[tempI][tempJ][tempK].matrix[l][m];
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumTempNodeMassMatrixParametersX(MassMatrix*** array) {
	if (cartCoord[0] > 0 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int i = 0; i < 2 * additionalBinNumber + 3; ++i) {
					for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
						for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
							for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
								for (int l = 0; l < 3; ++l) {
									for (int m = 0; m < 3; ++m) {
										array[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] += tempNodeMassMatrixParameterLeft[i][j][k].matrix[
											tempI][tempJ][tempK].matrix[l][m];
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//todo
	if (cartCoord[0] < cartDim[0] - 1 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int i = 0; i < 3 + 2 * additionalBinNumber; ++i) {
					for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
						for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
							for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
								for (int l = 0; l < 3; ++l) {
									for (int m = 0; m < 3; ++m) {
										array[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] +=
											tempNodeMassMatrixParameterRight[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumTempNodeMassMatrixParametersY(MassMatrix*** array) {
	if (cartCoord[1] > 0 || boundaryConditionTypeY == PERIODIC) {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int j = 0; j < 2 * additionalBinNumber + 3; ++j) {
					for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
						for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
							for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
								for (int l = 0; l < 3; ++l) {
									for (int m = 0; m < 3; ++m) {
										array[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] += tempNodeMassMatrixParameterFront[i][j][k].matrix[
											tempI][tempJ][tempK].matrix[l][m];
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if (cartCoord[1] < cartDim[1] - 1 || boundaryConditionTypeY == PERIODIC) {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int j = 0; j < 2 * additionalBinNumber + 3; ++j) {
					for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
						for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
							for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
								for (int l = 0; l < 3; ++l) {
									for (int m = 0; m < 3; ++m) {
										array[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k].matrix[tempI][tempJ][tempK].matrix[l][m] +=
											tempNodeMassMatrixParameterBack[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumTempNodeMassMatrixParametersZ(MassMatrix*** array) {
	if (cartCoord[2] > 0 || boundaryConditionTypeZ == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int k = 0; k < 2 * additionalBinNumber + 3; ++k) {
					for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
						for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
							for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
								for (int l = 0; l < 3; ++l) {
									for (int m = 0; m < 3; ++m) {
										array[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] += tempNodeMassMatrixParameterBottom[i][j][k].matrix[
											tempI][tempJ][tempK].matrix[l][m];
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if (cartCoord[1] < cartDim[1] - 1 || boundaryConditionTypeY == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int k = 0; k < 2 * additionalBinNumber + 3; ++k) {
					for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
						for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
							for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
								for (int l = 0; l < 3; ++l) {
									for (int m = 0; m < 3; ++m) {
										array[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k].matrix[tempI][tempJ][tempK].matrix[l][m] +=
											tempNodeMassMatrixParameterTop[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


void Simulation::sumNodeVectorParametersGeneralX(Vector3d*** vector, double* inBufferRight, double* outBufferRight,
                                                 double* inBufferLeft, double* outBufferLeft) {
	if ((verbosity > 2)) printf("sending left flux sum node vector parameters x rank = %d\n", rank);
	if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
		sendNodeVectorParametersToLeftReceiveFromRight(vector, outBufferLeft, tempNodeVectorParameterRight,
		                                               inBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
		                                               additionalBinNumber, cartComm, rank, leftRank, rightRank);
	} else if (cartCoord[0] == 0) {
		receiveNodeVectorParametersRight(tempNodeVectorParameterRight, inBufferRight, xnumberAdded, ynumberAdded,
		                                 znumberAdded, additionalBinNumber, cartComm, rank, rightRank);
	} else if (cartCoord[0] == cartDim[0] - 1) {
		sendNodeVectorParametersLeft(vector, outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
		                             additionalBinNumber, cartComm, rank, leftRank);
	}

	//MPI_Barrier(cartComm);

	if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
		sendNodeVectorParametersToRightReceiveFromLeft(vector, outBufferRight, tempNodeVectorParameterLeft,
		                                               inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
		                                               additionalBinNumber, cartComm, rank, leftRank, rightRank);
	} else if (cartCoord[0] == 0) {
		sendNodeVectorParametersRight(vector, outBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
		                              additionalBinNumber, cartComm, rank, rightRank);
	} else if (cartCoord[0] == cartDim[0] - 1) {
		receiveNodeVectorParametersLeft(tempNodeVectorParameterLeft, inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
		                                additionalBinNumber, cartComm, rank, leftRank);
	}


	sumTempNodeVectorParametersX(vector);
}

void Simulation::sumNodeVectorParametersGeneralY(Vector3d*** vector, double* inBufferBack, double* outBufferBack,
                                                 double* inBufferFront, double* outBufferFront) {
	if ((verbosity > 2)) printf("sending left flux sum node vector parameters x rank = %d\n", rank);
	if ((boundaryConditionTypeY == PERIODIC) || (cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1)) {
		sendNodeVectorParametersToFrontReceiveFromBack(vector, outBufferFront, tempNodeVectorParameterBack,
		                                               inBufferBack, xnumberAdded, ynumberAdded, znumberAdded,
		                                               additionalBinNumber, cartComm, rank, frontRank, backRank);
	} else if (cartCoord[1] == 0) {
		receiveNodeVectorParametersBack(tempNodeVectorParameterBack, inBufferBack, xnumberAdded, ynumberAdded,
		                                znumberAdded, additionalBinNumber, cartComm, rank, backRank);
	} else if (cartCoord[1] == cartDim[1] - 1) {
		sendNodeVectorParametersFront(vector, outBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
		                              additionalBinNumber, cartComm, rank, frontRank);
	}

	//MPI_Barrier(cartComm);

	if ((boundaryConditionTypeY == PERIODIC) || (cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1)) {
		sendNodeVectorParametersToBackReceiveFromFront(vector, outBufferBack, tempNodeVectorParameterFront,
		                                               inBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
		                                               additionalBinNumber, cartComm, rank, frontRank, backRank);
	} else if (cartCoord[1] == 0) {
		sendNodeVectorParametersBack(vector, outBufferBack, xnumberAdded, ynumberAdded, znumberAdded,
		                             additionalBinNumber, cartComm, rank, backRank);
	} else if (cartCoord[1] == cartDim[1] - 1) {
		receiveNodeVectorParametersFront(tempNodeVectorParameterFront, inBufferFront, xnumberAdded, ynumberAdded,
		                                 znumberAdded,
		                                 additionalBinNumber, cartComm, rank, frontRank);
	}


	sumTempNodeVectorParametersY(vector);
}

void Simulation::sumNodeVectorParametersGeneralZ(Vector3d*** vector, double* inBufferTop, double* outBufferTop,
                                                 double* inBufferBottom, double* outBufferBottom) {
	if ((verbosity > 2)) printf("sending left flux sum node vector parameters x rank = %d\n", rank);
	if ((boundaryConditionTypeZ == PERIODIC) || (cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1)) {
		sendNodeVectorParametersToBottomReceiveFromTop(vector, outBufferBottom, tempNodeVectorParameterTop,
		                                               inBufferTop, xnumberAdded, ynumberAdded, znumberAdded,
		                                               additionalBinNumber, cartComm, rank, bottomRank, topRank);
	} else if (cartCoord[2] == 0) {
		receiveNodeVectorParametersTop(tempNodeVectorParameterTop, inBufferTop, xnumberAdded, ynumberAdded,
		                               znumberAdded, additionalBinNumber, cartComm, rank, topRank);
	} else if (cartCoord[2] == cartDim[2] - 1) {
		sendNodeVectorParametersBottom(vector, outBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
		                               additionalBinNumber, cartComm, rank, bottomRank);
	}

	//MPI_Barrier(cartComm);

	if ((boundaryConditionTypeZ == PERIODIC) || (cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1)) {
		sendNodeVectorParametersToTopReceiveFromBottom(vector, outBufferTop, tempNodeVectorParameterBottom,
		                                               inBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
		                                               additionalBinNumber, cartComm, rank, bottomRank, topRank);
	} else if (cartCoord[2] == 0) {
		sendNodeVectorParametersTop(vector, outBufferTop, xnumberAdded, ynumberAdded, znumberAdded,
		                            additionalBinNumber, cartComm, rank, topRank);
	} else if (cartCoord[2] == cartDim[2] - 1) {
		receiveNodeVectorParametersBottom(tempNodeVectorParameterBottom, inBufferBottom, xnumberAdded, ynumberAdded,
		                                  znumberAdded,
		                                  additionalBinNumber, cartComm, rank, bottomRank);
	}


	sumTempNodeVectorParametersZ(vector);
}

void Simulation::sumNodeVectorParametersX() {
	if (cartDim[0] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum node vector parameters x\n");
		double* inBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * 3];
		double* outBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * 3];
		double* inBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * 3];
		double* outBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * 3];

		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralX(electricFlux, inBufferRight, outBufferRight, inBufferLeft, outBufferLeft);
		//////
		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralX(electricFluxMinus, inBufferRight, outBufferRight, inBufferLeft, outBufferLeft);

		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralX(divPressureTensor, inBufferRight, outBufferRight, inBufferLeft, outBufferLeft);

		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralX(divPressureTensorMinus, inBufferRight, outBufferRight, inBufferLeft, outBufferLeft);

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters x rank = %d\n", rank);
		delete[] inBufferLeft;
		delete[] inBufferRight;
		delete[] outBufferLeft;
		delete[] outBufferRight;

		if (debugMode) {
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].x, "electricFlux[i][j][k].x = NaN");
						if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].y, "electricFlux[i][j][k].y = NaN");
						if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].z, "electricFlux[i][j][k].z = NaN");
						alertLargeModlue(electricFlux[i][j][k].x, 1E100, "electric Flux x is too large\n");
						alertLargeModlue(electricFlux[i][j][k].y, 1E100, "electric Flux y is too large\n");
						alertLargeModlue(electricFlux[i][j][k].z, 1E100, "electric Flux z is too large\n");
					}
				}
			}
		}
	} else {
		if (boundaryConditionTypeX == PERIODIC) {
			for (int j = 0; j <= ynumberAdded; ++j) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int i = 0; i <= 2 * additionalBinNumber + 2; ++i) {
						for (int l = 0; l < 3; ++l) {
							electricFlux[xnumberAdded - i][j][k][l] += electricFlux[2 + 2 * additionalBinNumber - i][j][k][l];
							electricFlux[2 + 2 * additionalBinNumber - i][j][k][l] = electricFlux[xnumberAdded - i][j][k][l];

							electricFluxMinus[xnumberAdded - i][j][k][l] += electricFluxMinus[2 + 2 * additionalBinNumber - i][j][k][l];
							electricFluxMinus[2 + 2 * additionalBinNumber - i][j][k][l] = electricFluxMinus[xnumberAdded - i][j][k][l];

							divPressureTensor[xnumberAdded - i][j][k][l] += divPressureTensor[2 + 2 * additionalBinNumber - i][j][k][l];
							divPressureTensor[2 + 2 * additionalBinNumber - i][j][k][l] = divPressureTensor[xnumberAdded - i][j][k][l];

							divPressureTensorMinus[xnumberAdded - i][j][k][l] += divPressureTensorMinus[2 + 2 * additionalBinNumber - i][j][k][l];
							divPressureTensorMinus[2 + 2 * additionalBinNumber - i][j][k][l] = divPressureTensorMinus[xnumberAdded - i][j][k][l];
						}
					}
				}
			}
		}
	}
}


void Simulation::sumNodeVectorParametersY() {
	if (cartDim[1] > 1) {
		double* inBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * 3];
		double* outBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * 3];
		double* inBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * 3];
		double* outBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * 3];

		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralY(electricFlux, inBufferBack, outBufferBack, inBufferFront, outBufferFront);
		//////
		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralY(electricFluxMinus, inBufferBack, outBufferBack, inBufferFront, outBufferFront);

		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralY(divPressureTensor, inBufferBack, outBufferBack, inBufferFront, outBufferFront);

		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralY(divPressureTensorMinus, inBufferBack, outBufferBack, inBufferFront, outBufferFront);

		delete[] inBufferFront;
		delete[] inBufferBack;
		delete[] outBufferBack;
		delete[] outBufferFront;

		if (debugMode) {
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].x, "electricFlux[i][j][k].x = NaN");
						if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].y, "electricFlux[i][j][k].y = NaN");
						if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].z, "electricFlux[i][j][k].z = NaN");
						//alertLargeModlue(electricFlux[i][j][k].x, 1E100, "electric Flux x is too large\n");
						//alertLargeModlue(electricFlux[i][j][k].y, 1E100, "electric Flux y is too large\n");
						//alertLargeModlue(electricFlux[i][j][k].z, 1E100, "electric Flux z is too large\n");
					}
				}
			}
		}
	} else {
		if (boundaryConditionTypeY == PERIODIC) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= 2 * additionalBinNumber + 2; ++j) {
						electricFlux[i][ynumberAdded - j][k] += electricFlux[i][2 + 2 * additionalBinNumber - j][k];
						electricFlux[i][2 + 2 * additionalBinNumber - j][k] = electricFlux[i][ynumberAdded - j][k];

						electricFluxMinus[i][ynumberAdded - j][k] += electricFluxMinus[i][2 + 2 * additionalBinNumber - j][k];
						electricFluxMinus[i][2 + 2 * additionalBinNumber - j][k] = electricFluxMinus[i][ynumberAdded - j][k];

						divPressureTensor[i][ynumberAdded - j][k] += divPressureTensor[i][2 + 2 * additionalBinNumber - j][k];
						divPressureTensor[i][2 + 2 * additionalBinNumber - j][k] = divPressureTensor[i][ynumberAdded - j][k];

						divPressureTensorMinus[i][ynumberAdded - j][k] += divPressureTensorMinus[i][2 + 2 * additionalBinNumber - j][k];
						divPressureTensorMinus[i][2 + 2 * additionalBinNumber - j][k] = divPressureTensorMinus[i][ynumberAdded - j][k];
					}
				}
			}
		}

		if (ynumberGeneral == 1) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= ynumberAdded; ++j) {
						electricFlux[i][j][k] = electricFlux[i][1 + additionalBinNumber][k];
						electricFluxMinus[i][j][k] = electricFluxMinus[i][1 + additionalBinNumber][k];
						divPressureTensor[i][j][k] = divPressureTensor[i][1 + additionalBinNumber][k];
						divPressureTensorMinus[i][j][k] = divPressureTensorMinus[i][1 + additionalBinNumber][k];
					}
				}
			}
		}
	}
}

void Simulation::sumNodeVectorParametersZ() {
	if (cartDim[2] > 1) {
		double* inBufferTop = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1) * 3];
		double* outBufferTop = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1) * 3];
		double* inBufferBottom = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1) * 3];
		double* outBufferBottom = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1) * 3];

		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralZ(electricFlux, inBufferTop, outBufferTop, inBufferBottom, outBufferBottom);
		//////
		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralZ(electricFluxMinus, inBufferTop, outBufferTop, inBufferBottom, outBufferBottom);

		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralZ(divPressureTensor, inBufferTop, outBufferTop, inBufferBottom, outBufferBottom);

		//MPI_Barrier(cartComm);
		sumNodeVectorParametersGeneralZ(divPressureTensorMinus, inBufferTop, outBufferTop, inBufferBottom, outBufferBottom);

		delete[] inBufferTop;
		delete[] inBufferBottom;
		delete[] outBufferBottom;
		delete[] outBufferTop;

		if (debugMode) {
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].x, "electricFlux[i][j][k].x = NaN");
						if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].y, "electricFlux[i][j][k].y = NaN");
						if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].z, "electricFlux[i][j][k].z = NaN");
						//alertLargeModlue(electricFlux[i][j][k].x, 1E100, "electric Flux x is too large\n");
						//alertLargeModlue(electricFlux[i][j][k].y, 1E100, "electric Flux y is too large\n");
						//alertLargeModlue(electricFlux[i][j][k].z, 1E100, "electric Flux z is too large\n");
					}
				}
			}
		}
	} else {
		if (boundaryConditionTypeZ == PERIODIC) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int j = 0; j <= ynumberAdded; ++j) {
					for (int k = 0; k <= 2 * additionalBinNumber + 2; ++k) {
						for (int l = 0; l < 3; ++l) {
							electricFlux[i][j][znumberAdded - k][l] += electricFlux[i][j][2 + 2 * additionalBinNumber - k][l];
							electricFlux[i][j][2 + 2 * additionalBinNumber - k][l] = electricFlux[i][j][znumberAdded - k][l];

							electricFluxMinus[i][j][znumberAdded - k][l] += electricFluxMinus[i][j][2 + 2 * additionalBinNumber - k][l];
							electricFluxMinus[i][j][2 + 2 * additionalBinNumber - k][l] = electricFluxMinus[i][j][znumberAdded - k][l];

							divPressureTensor[i][j][znumberAdded - k][l] += divPressureTensor[i][j][2 + 2 * additionalBinNumber - k][l];
							divPressureTensor[i][j][2 + 2 * additionalBinNumber - k][l] = divPressureTensor[i][j][znumberAdded - k][l];

							divPressureTensorMinus[i][j][znumberAdded - k][l] += divPressureTensorMinus[i][j][2 + 2 * additionalBinNumber - k][l];
							divPressureTensorMinus[i][j][2 + 2 * additionalBinNumber - k][l] = divPressureTensorMinus[i][j][znumberAdded - k][l];
						}
					}
				}
			}
		}
		if (znumberGeneral == 1) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= ynumberAdded; ++j) {
						electricFlux[i][j][k] = electricFlux[i][j][1 + additionalBinNumber];
						electricFluxMinus[i][j][k] = electricFluxMinus[i][j][1 + additionalBinNumber];
						divPressureTensor[i][j][k] = divPressureTensor[i][j][1 + additionalBinNumber];
						divPressureTensorMinus[i][j][k] = divPressureTensorMinus[i][j][1 + additionalBinNumber];
					}
				}
			}
		}
	}
}

void Simulation::sumTempNodeVectorParametersX(Vector3d*** array) {
	if (cartCoord[0] > 0 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int i = 0; i < 2 * additionalBinNumber + 3; ++i) {
					array[i][j][k] += tempNodeVectorParameterLeft[i][j][k];
				}
			}
		}
	}
	//todo
	if (cartCoord[0] < cartDim[0] - 1 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int i = 0; i < 3 + 2 * additionalBinNumber; ++i) {
					array[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k] += tempNodeVectorParameterRight[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumTempNodeVectorParametersY(Vector3d*** array) {
	if (cartCoord[1] > 0 || boundaryConditionTypeY == PERIODIC) {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int j = 0; j < 2 * additionalBinNumber + 3; ++j) {
					array[i][j][k] += tempNodeVectorParameterFront[i][j][k];
				}
			}
		}
	}
	if (cartCoord[1] < cartDim[1] - 1 || boundaryConditionTypeY == PERIODIC) {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int j = 0; j < 2 * additionalBinNumber + 3; ++j) {
					array[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k] += tempNodeVectorParameterBack[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumTempNodeVectorParametersZ(Vector3d*** array) {
	if (cartCoord[2] > 0 || boundaryConditionTypeZ == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int k = 0; k < 2 * additionalBinNumber + 3; ++k) {
					array[i][j][k] += tempNodeVectorParameterBottom[i][j][k];
				}
			}
		}
	}
	if (cartCoord[2] < cartDim[2] - 1 || boundaryConditionTypeZ == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int k = 0; k < 2 * additionalBinNumber + 3; ++k) {
					array[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k] += tempNodeVectorParameterTop[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumNodeMatrixParametersX() {
	if (cartDim[0] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum node matrix parameters x rank = %d\n", rank);
		double* inBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * 9];
		double* outBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * 9];
		double* inBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * 9];
		double* outBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1) * 9];

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("send left dielectric tensor sum node matrix parameters x rank = %d\n", rank);
		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendNodeMatrixParametersToLeftReceiveFromRight(dielectricTensor, outBufferLeft, tempNodeMatrixParameterRight,
			                                               inBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
			                                               additionalBinNumber, cartComm, rank, leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			receiveNodeMatrixParametersRight(tempNodeMatrixParameterRight, inBufferRight, xnumberAdded, ynumberAdded,
			                                 znumberAdded, additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			sendNodeMatrixParametersLeft(dielectricTensor, outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
			                             additionalBinNumber, cartComm, rank, leftRank);
		}

		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendNodeMatrixParametersToRightReceiveFromLeft(dielectricTensor, outBufferRight, tempNodeMatrixParameterLeft,
			                                               inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
			                                               additionalBinNumber, cartComm, rank, leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			sendNodeMatrixParametersRight(dielectricTensor, outBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
			                              additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			receiveNodeMatrixParametersLeft(tempNodeMatrixParameterLeft, inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
			                                additionalBinNumber, cartComm, rank, leftRank);
		}

		sumTempNodeMatrixParametersX(dielectricTensor);
		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node matrix parameters x rank = %d\n", rank);
		delete[] inBufferLeft;
		delete[] inBufferRight;
		delete[] outBufferLeft;
		delete[] outBufferRight;
	} else {
		if (boundaryConditionTypeX == PERIODIC) {
			for (int j = 0; j <= ynumberAdded; ++j) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int l = 0; l < 3; ++l) {
						for (int m = 0; m < 3; ++m) {
							for (int i = 0; i <= 2 * additionalBinNumber + 2; ++i) {
								dielectricTensor[xnumberAdded - i][j][k].matrix[l][m] += dielectricTensor[2 + 2 * additionalBinNumber - i][j][k]
									.matrix[l][m];
								dielectricTensor[2 + 2 * additionalBinNumber - i][j][k].matrix[l][m] = dielectricTensor[xnumberAdded - i][j][k].
									matrix[l][m];
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumNodeMatrixParametersY() {
	if (cartDim[1] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum node matrix parameters y rank = %d\n", rank);
		double* inBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * 9];
		double* outBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * 9];
		double* inBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * 9];
		double* outBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1) * 9];

		if ((boundaryConditionTypeY == PERIODIC) || (cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1)) {
			sendNodeMatrixParametersToFrontReceiveFromBack(dielectricTensor, outBufferFront, tempNodeMatrixParameterBack,
			                                               inBufferBack, xnumberAdded, ynumberAdded, znumberAdded,
			                                               additionalBinNumber, cartComm, rank, frontRank, backRank);
		} else if (cartCoord[1] == 0) {
			receiveNodeMatrixParametersBack(tempNodeMatrixParameterBack, inBufferBack, xnumberAdded, ynumberAdded,
			                                znumberAdded, additionalBinNumber, cartComm, rank, backRank);
		} else if (cartCoord[1] == cartDim[1] - 1) {
			sendNodeMatrixParametersFront(dielectricTensor, outBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
			                              additionalBinNumber, cartComm, rank, frontRank);
		}

		if ((boundaryConditionTypeY == PERIODIC) || (cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1)) {
			sendNodeMatrixParametersToBackReceiveFromFront(dielectricTensor, outBufferBack, tempNodeMatrixParameterFront,
			                                               inBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
			                                               additionalBinNumber, cartComm, rank, frontRank, backRank);
		} else if (cartCoord[1] == 0) {
			sendNodeMatrixParametersBack(dielectricTensor, outBufferBack, xnumberAdded, ynumberAdded, znumberAdded,
			                             additionalBinNumber, cartComm, rank, backRank);
		} else if (cartCoord[1] == cartDim[1] - 1) {
			receiveNodeMatrixParametersFront(tempNodeMatrixParameterFront, inBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
			                                 additionalBinNumber, cartComm, rank, frontRank);
		}


		sumTempNodeMatrixParametersY(dielectricTensor);
		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node matrix parameters y rank = %d\n", rank);
		delete[] inBufferFront;
		delete[] inBufferBack;
		delete[] outBufferFront;
		delete[] outBufferBack;
	} else {
		if (boundaryConditionTypeY == PERIODIC) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= 2 * additionalBinNumber + 2; ++j) {
						dielectricTensor[i][ynumberAdded - j][k] = dielectricTensor[i][ynumberAdded - j][k] + dielectricTensor[i][2 + 2 *
							additionalBinNumber - j][k];
						dielectricTensor[i][2 + 2 * additionalBinNumber - j][k] = dielectricTensor[i][ynumberAdded - j][k];
					}
				}
			}
		}
		if (ynumberGeneral == 1) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= ynumberAdded; ++j) {
						dielectricTensor[i][j][k] = dielectricTensor[i][1 + additionalBinNumber][k];
					}
				}
			}
		}
	}
}

void Simulation::sumNodeMatrixParametersZ() {
	if (cartDim[2] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum node matrix parameters z rank = %d\n", rank);
		double* inBufferBottom = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded + 1) * 9];
		double* outBufferBottom = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded + 1) * 9];
		double* inBufferTop = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded + 1) * 9];
		double* outBufferTop = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded + 1) * 9];

		if ((boundaryConditionTypeZ == PERIODIC) || (cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1)) {
			sendNodeMatrixParametersToBottomReceiveFromTop(dielectricTensor, outBufferBottom, tempNodeMatrixParameterTop,
			                                               inBufferTop, xnumberAdded, ynumberAdded, znumberAdded,
			                                               additionalBinNumber, cartComm, rank, bottomRank, topRank);
		} else if (cartCoord[2] == 0) {
			receiveNodeMatrixParametersTop(tempNodeMatrixParameterTop, inBufferTop, xnumberAdded, ynumberAdded,
			                               znumberAdded, additionalBinNumber, cartComm, rank, topRank);
		} else if (cartCoord[2] == cartDim[2] - 1) {
			sendNodeMatrixParametersBottom(dielectricTensor, outBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
			                               additionalBinNumber, cartComm, rank, bottomRank);
		}

		if ((boundaryConditionTypeZ == PERIODIC) || (cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1)) {
			sendNodeMatrixParametersToTopReceiveFromBottom(dielectricTensor, outBufferTop, tempNodeMatrixParameterBottom,
			                                               inBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
			                                               additionalBinNumber, cartComm, rank, bottomRank, topRank);
		} else if (cartCoord[2] == 0) {
			sendNodeMatrixParametersTop(dielectricTensor, outBufferTop, xnumberAdded, ynumberAdded, znumberAdded,
			                            additionalBinNumber, cartComm, rank, topRank);
		} else if (cartCoord[2] == cartDim[2] - 1) {
			receiveNodeMatrixParametersBottom(tempNodeMatrixParameterBottom, inBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
			                                  additionalBinNumber, cartComm, rank, bottomRank);
		}

		sumTempNodeMatrixParametersZ(dielectricTensor);
		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node matrix parameters z rank = %d\n", rank);
		delete[] inBufferBottom;
		delete[] inBufferTop;
		delete[] outBufferBottom;
		delete[] outBufferTop;
	} else {
		if (boundaryConditionTypeZ == PERIODIC) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int j = 0; j <= ynumberAdded; ++j) {
					for (int k = 0; k <= 2 * additionalBinNumber + 2; ++k) {
						dielectricTensor[i][j][znumberAdded - k] = dielectricTensor[i][j][znumberAdded - k] + dielectricTensor[i][j][2 + 2
							* additionalBinNumber - k];
						dielectricTensor[i][j][2 + 2 * additionalBinNumber - k] = dielectricTensor[i][j][znumberAdded - k];
					}
				}
			}
		}
		if (znumberGeneral == 1) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int j = 0; j <= ynumberAdded; ++j) {
					for (int k = 0; k <= additionalBinNumber + 2; ++k) {
						dielectricTensor[i][j][k] = dielectricTensor[i][j][1 + additionalBinNumber];
					}
				}
			}
		}
	}
}

void Simulation::sumTempNodeMatrixParametersX(Matrix3d*** array) {
	if (cartCoord[0] > 0 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int i = 0; i < 2 * additionalBinNumber + 3; ++i) {
					array[i][j][k] += tempNodeMatrixParameterLeft[i][j][k];
				}
			}
		}
	}

	if (cartCoord[0] < cartDim[0] - 1 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int i = 0; i < 3 + 2 * additionalBinNumber; ++i) {
					array[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k] += tempNodeMatrixParameterRight[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumTempNodeMatrixParametersY(Matrix3d*** array) {
	if (cartCoord[1] > 0 || boundaryConditionTypeY == PERIODIC) {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int j = 0; j < 2 * additionalBinNumber + 3; ++j) {
					array[i][j][k] += tempNodeMatrixParameterFront[i][j][k];
				}
			}
		}
	}
	if (cartCoord[1] < cartDim[1] - 1 || boundaryConditionTypeY == PERIODIC) {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int j = 0; j < 2 * additionalBinNumber + 3; ++j) {
					array[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k] += tempNodeMatrixParameterBack[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumTempNodeMatrixParametersZ(Matrix3d*** array) {
	if (cartCoord[2] > 0 || boundaryConditionTypeZ == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int k = 0; k < 2 * additionalBinNumber + 3; ++k) {
					array[i][j][k] += tempNodeMatrixParameterBottom[i][j][k];
				}
			}
		}
	}
	if (cartCoord[2] < cartDim[2] - 1 || boundaryConditionTypeZ == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int k = 0; k < 2 * additionalBinNumber + 3; ++k) {
					array[i][j][znumberAdded - 2 - 2 * additionalBinNumber - k] += tempNodeMatrixParameterTop[i][j][k];
				}
			}
		}
	}
}


void Simulation::sumTempNodeParametersX(double*** array) {
	if (cartCoord[0] > 0 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int i = 0; i < 2 * additionalBinNumber + 3; ++i) {
					array[i][j][k] += tempNodeParameterLeft[i][j][k];
				}
			}
		}
	}

	if (cartCoord[0] < cartDim[0] - 1 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int i = 0; i < 3 + 2 * additionalBinNumber; ++i) {
					array[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k] += tempNodeParameterRight[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumTempNodeParametersY(double*** array) {
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int k = 0; k < znumberAdded + 1; ++k) {
			for (int j = 0; j < 2 * additionalBinNumber + 3; ++j) {
				array[i][j][k] += tempNodeParameterFront[i][j][k];
				array[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k] += tempNodeParameterBack[i][j][k];
			}
		}
	}
}

void Simulation::sumTempNodeParametersZ(double*** array) {
	for (int j = 0; j < ynumberAdded + 1; ++j) {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < 2 * additionalBinNumber + 3; ++k) {
				array[i][j][k] += tempNodeParameterBottom[i][j][k];
				array[i][j][znumberAdded - 2 - 2 * additionalBinNumber - k] += tempNodeParameterTop[i][j][k];
			}
		}
	}
}

void Simulation::sumCellParametersGeneralX(double*** array, double* inBufferRight, double* outBufferRight, double* inBufferLeft, double* outBufferLeft) {
	if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
		sendCellParametersToLeftReceiveFromRight(array, outBufferLeft, tempCellParameterRight, inBufferRight,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, leftRank, rightRank);
	} else if (cartCoord[0] == 0) {
		receiveCellParametersRight(tempCellParameterRight, inBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
		                           additionalBinNumber, cartComm, rank, rightRank);
	} else if (cartCoord[0] == cartDim[0] - 1) {
		sendCellParametersLeft(array, outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
		                       additionalBinNumber, cartComm, rank, leftRank);
	}

	if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
		sendCellParametersToRightReceiveFromLeft(array, outBufferRight, tempCellParameterLeft, inBufferLeft,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, leftRank, rightRank);
	} else if (cartCoord[0] == 0) {
		sendCellParametersRight(array, outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
		                        additionalBinNumber, cartComm, rank, rightRank);
	} else if (cartCoord[0] == cartDim[0] - 1) {
		receiveCellParametersLeft(tempCellParameterLeft, inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
		                          additionalBinNumber, cartComm, rank, leftRank);
	}

	sumCellTempParametersX(array);
}

void Simulation::sumCellParametersGeneralY(double*** array, double* inBufferBack, double* outBufferBack, double* inBufferFront, double* outBufferFront) {
	if (boundaryConditionTypeY == PERIODIC || (cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1)) {
		sendCellParametersToFrontReceiveFromBack(array, outBufferFront, tempCellParameterBack, inBufferBack,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, frontRank, backRank);
	} else if (cartCoord[1] == 0) {
		receiveCellParametersBack(tempCellParameterBack, inBufferBack, xnumberAdded, ynumberAdded, znumberAdded,
		                           additionalBinNumber, cartComm, rank, backRank);
	} else if (cartCoord[1] == cartDim[1] - 1) {
		sendCellParametersFront(array, outBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
		                       additionalBinNumber, cartComm, rank, frontRank);
	}

	if (boundaryConditionTypeY == PERIODIC || (cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1)) {
		sendCellParametersToBackReceiveFromFront(array, outBufferBack, tempCellParameterFront, inBufferFront,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, frontRank, backRank);
	} else if (cartCoord[1] == 0) {
		sendCellParametersBack(array, outBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
		                        additionalBinNumber, cartComm, rank, backRank);
	} else if (cartCoord[1] == cartDim[1] - 1) {
		receiveCellParametersFront(tempCellParameterFront, inBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
		                          additionalBinNumber, cartComm, rank, frontRank);
	}

	sumCellTempParametersY(array);
}

void Simulation::sumCellParametersGeneralZ(double*** array, double* inBufferTop, double* outBufferTop, double* inBufferBottom, double* outBufferBottom) {
	if (boundaryConditionTypeZ == PERIODIC || (cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1)) {
		sendCellParametersToBottomReceiveFromTop(array, outBufferBottom, tempCellParameterTop, inBufferTop,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, bottomRank, topRank);
	} else if (cartCoord[2] == 0) {
		receiveCellParametersTop(tempCellParameterTop, inBufferTop, xnumberAdded, ynumberAdded, znumberAdded,
		                           additionalBinNumber, cartComm, rank, topRank);
	} else if (cartCoord[2] == cartDim[2] - 1) {
		sendCellParametersBottom(array, outBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
		                       additionalBinNumber, cartComm, rank, bottomRank);
	}

	if (boundaryConditionTypeZ == PERIODIC || (cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1)) {
		sendCellParametersToTopReceiveFromBottom(array, outBufferTop, tempCellParameterBottom, inBufferBottom,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, bottomRank, topRank);
	} else if (cartCoord[2] == 0) {
		sendCellParametersTop(array, outBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
		                        additionalBinNumber, cartComm, rank, topRank);
	} else if (cartCoord[2] == cartDim[2] - 1) {
		receiveCellParametersBottom(tempCellParameterBottom, inBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
		                          additionalBinNumber, cartComm, rank, bottomRank);
	}

	sumCellTempParametersZ(array);
}

void Simulation::sumChargeDensityHatX() {
	if (cartDim[0] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum charge density hat rank = %d\n", rank);
		double* inBufferRight = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * znumberAdded];
		double* outBufferRight = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * znumberAdded];
		double* inBufferLeft = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * znumberAdded];
		double* outBufferLeft = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * znumberAdded];
		MPI_Barrier(cartComm);
		//MPI_Barrier(cartComm);
		if ((rank == 0) && (verbosity > 2)) printf("sending left sum charge density hat x rank = %d\n", rank);

		sumCellParametersGeneralX(chargeDensityHat, inBufferRight, outBufferRight, inBufferLeft, outBufferLeft);
		/////
		
		sumCellParametersGeneralX(chargeDensityMinus, inBufferRight, outBufferRight, inBufferLeft, outBufferLeft);
		

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum charge density hat x rank = %d\n", rank);
		delete[] inBufferRight;
		delete[] outBufferRight;
		delete[] inBufferLeft;
		delete[] outBufferLeft;
	} else {
		if (boundaryConditionTypeX == PERIODIC) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i < 2 * additionalBinNumber + 2; ++i) {
						chargeDensityHat[xnumberAdded - 1 - i][j][k] += chargeDensityHat[1 + 2 * additionalBinNumber - i][j][k];
						chargeDensityHat[1 + 2 * additionalBinNumber - i][j][k] = chargeDensityHat[xnumberAdded - 1 - i][j][k];

						chargeDensityMinus[xnumberAdded - 1 - i][j][k] += chargeDensityMinus[1 + 2 * additionalBinNumber - i][j][k];
						chargeDensityMinus[1 + 2 * additionalBinNumber - i][j][k] = chargeDensityMinus[xnumberAdded - 1 - i][j][k];
					}
				}
			}
		}
	}
}

void Simulation::sumChargeDensityHatY() {
	if (cartDim[1] > 1) {
		double* inBufferBack = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * znumberAdded];
		double* outBufferBack = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * znumberAdded];
		double* inBufferFront = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * znumberAdded];
		double* outBufferFront = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * znumberAdded];

		sumCellParametersGeneralY(chargeDensityHat, inBufferBack, outBufferBack, inBufferFront, outBufferFront);

		sumCellParametersGeneralY(chargeDensityMinus, inBufferBack, outBufferBack, inBufferFront, outBufferFront);


		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum charge density hat y rank = %d\n", rank);
		delete[] inBufferFront;
		delete[] outBufferFront;
		delete[] inBufferBack;
		delete[] outBufferBack;
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				if(boundaryConditionTypeY == PERIODIC){
					for (int j = 0; j < 2 * additionalBinNumber + 2; ++j) {
						chargeDensityHat[i][ynumberAdded - 1 - j][k] += chargeDensityHat[i][1 + 2 * additionalBinNumber - j][k];
						chargeDensityHat[i][1 + 2 * additionalBinNumber - j][k] = chargeDensityHat[i][ynumberAdded - 1 - j][k];

						chargeDensityMinus[i][ynumberAdded - 1 - j][k] += chargeDensityMinus[i][1 + 2 * additionalBinNumber - j][k];
						chargeDensityMinus[i][1 + 2 * additionalBinNumber - j][k] = chargeDensityMinus[i][ynumberAdded - 1 - j][k];
					}
				}
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded; ++j) {
						chargeDensityHat[i][j][k] = chargeDensityHat[i][1 + additionalBinNumber][k];
						chargeDensityMinus[i][j][k] = chargeDensityMinus[i][1 + additionalBinNumber][k];
					}
				}
			}
		}
	}
}

void Simulation::sumChargeDensityHatZ() {
	if (cartDim[2] > 1) {
		double* inBufferTop = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * ynumberAdded];
		double* outBufferTop = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * ynumberAdded];
		double* inBufferBottom = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * ynumberAdded];
		double* outBufferBottom = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * ynumberAdded];

		sumCellParametersGeneralZ(chargeDensityHat, inBufferTop, outBufferTop, inBufferBottom, outBufferBottom);

		sumCellParametersGeneralZ(chargeDensityMinus, inBufferTop, outBufferTop, inBufferBottom, outBufferBottom);

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum charge density hat z rank = %d\n", rank);
		delete[] inBufferBottom;
		delete[] outBufferBottom;
		delete[] inBufferTop;
		delete[] outBufferTop;
	} else {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int i = 0; i < xnumberAdded; ++i) {
				if(boundaryConditionTypeZ == PERIODIC){
					for (int k = 0; k < 2 * additionalBinNumber + 2; ++k) {
						chargeDensityHat[i][j][znumberAdded - 1 - k] += chargeDensityHat[i][j][1 + 2 * additionalBinNumber - k];
						chargeDensityHat[i][j][1 + 2 * additionalBinNumber - k] = chargeDensityHat[i][j][znumberAdded - 1 - k];
						chargeDensityMinus[i][j][znumberAdded - 1 - k] += chargeDensityMinus[i][j][1 + 2 * additionalBinNumber - k];
						chargeDensityMinus[i][j][1 + 2 * additionalBinNumber - k] = chargeDensityMinus[i][j][znumberAdded - 1 - k];
					}
				}
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded; ++k) {
						chargeDensityHat[i][j][k] = chargeDensityHat[i][j][1 + additionalBinNumber];
						chargeDensityMinus[i][j][k] = chargeDensityMinus[i][j][1 + additionalBinNumber];
					}
				}
			}
		}
	}
}

void Simulation::sumCellParametersX() {
	if (cartDim[0] > 1) {
		double* inBufferRight = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * znumberAdded];
		double* outBufferRight = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * znumberAdded];
		double* inBufferLeft = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * znumberAdded];
		double* outBufferLeft = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * znumberAdded];

		for (int t = 0; t < typesNumber; ++t) {
			//MPI_Barrier(cartComm);
			if (types[t].particlesPerBin > 0) {
				if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
					sendCellParametersToLeftReceiveFromRight(particleConcentrations[t], outBufferLeft, tempCellParameterRight,
					                                         inBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
					                                         additionalBinNumber, cartComm, rank, leftRank, rightRank);
				} else if (cartCoord[0] == 0) {
					receiveCellParametersRight(tempCellParameterRight, inBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
					                           additionalBinNumber, cartComm, rank, rightRank);
				} else if (cartCoord[0] == cartDim[0] - 1) {
					sendCellParametersLeft(particleConcentrations[t], outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
					                       additionalBinNumber, cartComm, rank, leftRank);
				}

				if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
					sendCellParametersToRightReceiveFromLeft(particleConcentrations[t], outBufferRight, tempCellParameterLeft,
					                                         inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
					                                         additionalBinNumber, cartComm, rank, leftRank, rightRank);
				} else if (cartCoord[0] == 0) {
					sendCellParametersRight(particleConcentrations[t], outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
					                        additionalBinNumber, cartComm, rank, rightRank);
				} else if (cartCoord[0] == cartDim[0] - 1) {
					receiveCellParametersLeft(tempCellParameterLeft, inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
					                          additionalBinNumber, cartComm, rank, leftRank);
				}

				sumCellTempParametersX(particleConcentrations[t]);

			}
		}

		for (int t = 0; t < typesNumber; ++t) {
			//MPI_Barrier(cartComm);
			if (types[t].particlesPerBin > 0) {
				if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
					sendCellParametersToLeftReceiveFromRight(particleEnergies[t], outBufferLeft, tempCellParameterRight, inBufferRight,
					                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
					                                         rank, leftRank, rightRank);
				} else if (cartCoord[0] == 0) {
					receiveCellParametersRight(tempCellParameterRight, inBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
					                           additionalBinNumber, cartComm, rank, rightRank);
				} else if (cartCoord[0] == cartDim[0] - 1) {
					sendCellParametersLeft(particleEnergies[t], outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
					                       additionalBinNumber, cartComm, rank, leftRank);
				}

				if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
					sendCellParametersToRightReceiveFromLeft(particleEnergies[t], outBufferRight, tempCellParameterLeft, inBufferLeft,
					                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
					                                         rank, leftRank, rightRank);
				} else if (cartCoord[0] == 0) {
					sendCellParametersRight(particleEnergies[t], outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
					                        additionalBinNumber, cartComm, rank, rightRank);
				} else if (cartCoord[0] == cartDim[0] - 1) {
					receiveCellParametersLeft(tempCellParameterLeft, inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
					                          additionalBinNumber, cartComm, rank, leftRank);
				}

				sumCellTempParametersX(particleEnergies[t]);

			}
		}
		//MPI_Barrier(cartComm);


		if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendCellParametersToLeftReceiveFromRight(chargeDensity, outBufferLeft, tempCellParameterRight, inBufferRight,
			                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
			                                         rank, leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			receiveCellParametersRight(tempCellParameterRight, inBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
			                           additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			sendCellParametersLeft(chargeDensity, outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber,
			                       cartComm, rank, leftRank);
		}

		if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendCellParametersToRightReceiveFromLeft(chargeDensity, outBufferRight, tempCellParameterLeft, inBufferLeft,
			                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
			                                         rank, leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			sendCellParametersRight(chargeDensity, outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber,
			                        cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			receiveCellParametersLeft(tempCellParameterLeft, inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
			                          additionalBinNumber, cartComm, rank, leftRank);
		}

		sumCellTempParametersX(chargeDensity);

		//MPI_Barrier(cartComm);

		if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
			//todo!
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						chargeDensity[xnumberAdded - 2 - additionalBinNumber][j][k] += chargeDensity[xnumberAdded - 1 -
							additionalBinNumber + i][j][k];
						chargeDensity[xnumberAdded - 1 - additionalBinNumber + i][j][k] = 0;
						for (int t = 0; t < typesNumber; ++t) {
							particleConcentrations[t][xnumberAdded - 2 - additionalBinNumber][j][k] += particleConcentrations[t][xnumberAdded
								- 1 - additionalBinNumber + i][j][k];
							particleConcentrations[t][xnumberAdded - 1 - additionalBinNumber + i][j][k] = 0;
							particleEnergies[t][xnumberAdded - 2 - additionalBinNumber][j][k] += particleEnergies[t][xnumberAdded - 1 -
								additionalBinNumber + i][j][k];
							particleEnergies[t][xnumberAdded - 1 - additionalBinNumber + i][j][k] = 0;
						}
					}
				}
			}
		}

		delete[] inBufferRight;
		delete[] outBufferRight;
		delete[] inBufferLeft;
		delete[] outBufferLeft;
	} else {
		if (boundaryConditionTypeX == PERIODIC) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i < 2 + 2 * additionalBinNumber; ++i) {
						for (int t = 0; t < typesNumber; ++t) {
							particleConcentrations[t][xnumberAdded - 1 - i][j][k] += particleConcentrations[t][2 * additionalBinNumber + 1 -
								i][j][k];
							particleConcentrations[t][2 * additionalBinNumber + 1 - i][j][k] = particleConcentrations[t][xnumberAdded - 1 - i
							][j][k];

							particleEnergies[t][xnumberAdded - 1 - i][j][k] += particleEnergies[t][2 * additionalBinNumber + 1 - i][j][k];
							particleEnergies[t][2 * additionalBinNumber + 1 - i][j][k] = particleEnergies[t][xnumberAdded - 1 - i][j][k];
						}
						chargeDensity[xnumberAdded - 1 - i][j][k] += chargeDensity[2 * additionalBinNumber + 1 - i][j][k];
						chargeDensity[2 * additionalBinNumber + 1 - i][j][k] = chargeDensity[xnumberAdded - 1 - i][j][k];
					}
				}
			}
		} else {
			//manualy set zero to the last bin //todo whaaat?
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						chargeDensity[xnumberAdded - 2 - additionalBinNumber][j][k] += chargeDensity[xnumberAdded - 1 -
							additionalBinNumber + i][j][k];
						chargeDensity[xnumberAdded - 1 - additionalBinNumber + i][j][k] = 0;
						for (int t = 0; t < typesNumber; ++t) {
							particleConcentrations[t][xnumberAdded - 2 - additionalBinNumber][j][k] += particleConcentrations[t][xnumberAdded
								- 1 - additionalBinNumber + i][j][k];
							particleConcentrations[t][xnumberAdded - 1 - additionalBinNumber + i][j][k] = 0;

							particleEnergies[t][xnumberAdded - 2 - additionalBinNumber][j][k] += particleEnergies[t][xnumberAdded - 1 -
								additionalBinNumber + i][j][k];
							particleEnergies[t][xnumberAdded - 1 - additionalBinNumber + i][j][k] = 0;
						}
					}
				}
			}

			if (boundaryConditionTypeX == FREE_BOTH) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						for (int i = 0; i <= additionalBinNumber; ++i) {
							chargeDensity[1 + additionalBinNumber][j][k] += chargeDensity[additionalBinNumber - i][j][k];
							chargeDensity[xnumberAdded - 1 - additionalBinNumber + i][j][k] = 0;
							for (int t = 0; t < typesNumber; ++t) {
								particleConcentrations[t][xnumberAdded - 2 - additionalBinNumber][j][k] += particleConcentrations[t][
									xnumberAdded - 1 - additionalBinNumber + i][j][k];
								particleConcentrations[t][xnumberAdded - 1 - additionalBinNumber + i][j][k] = 0;

								particleEnergies[t][xnumberAdded - 2 - additionalBinNumber][j][k] += particleEnergies[t][xnumberAdded - 1 -
									additionalBinNumber + i][j][k];
								particleEnergies[t][xnumberAdded - 1 - additionalBinNumber + i][j][k] = 0;
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumCellParametersY() {
	if (cartDim[1] > 1) {
		double* inBufferFront = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * znumberAdded];
		double* outBufferFront = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * znumberAdded];
		double* inBufferBack = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * znumberAdded];
		double* outBufferBack = new double[(2 + 2 * additionalBinNumber) * xnumberAdded * znumberAdded];

		for (int t = 0; t < typesNumber; ++t) {
			if (types[t].particlesPerBin > 0) {

				sendCellParametersToBackReceiveFromFront(particleConcentrations[t], outBufferBack, tempCellParameterFront,
				                                         inBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
				                                         additionalBinNumber, cartComm, rank, backRank, frontRank);
				sendCellParametersToFrontReceiveFromBack(particleConcentrations[t], outBufferFront, tempCellParameterBack,
				                                         inBufferBack, xnumberAdded, ynumberAdded, znumberAdded,
				                                         additionalBinNumber, cartComm, rank, backRank, frontRank);

				/*sendCellParametersBack(particleConcentrations[t], outBufferBack, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, backRank);
				receiveCellParametersFront(tempCellParameterFront, inBufferFront, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, frontRank);
				sendCellParametersFront(particleConcentrations[t], outBufferFront, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, frontRank);
				receiveCellParametersBack(tempCellParameterBack, inBufferBack, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, backRank);*/

				sumCellTempParametersY(particleConcentrations[t]);
			}
		}

		for (int t = 0; t < typesNumber; ++t) {
			if (types[t].particlesPerBin > 0) {

				sendCellParametersToBackReceiveFromFront(particleEnergies[t], outBufferBack, tempCellParameterFront, inBufferFront,
				                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
				                                         rank, backRank, frontRank);
				sendCellParametersToFrontReceiveFromBack(particleEnergies[t], outBufferFront, tempCellParameterBack, inBufferBack,
				                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
				                                         rank, backRank, frontRank);

				/*sendCellParametersBack(particleConcentrations[t], outBufferBack, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, backRank);
				receiveCellParametersFront(tempCellParameterFront, inBufferFront, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, frontRank);
				sendCellParametersFront(particleConcentrations[t], outBufferFront, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, frontRank);
				receiveCellParametersBack(tempCellParameterBack, inBufferBack, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, backRank);*/

				sumCellTempParametersY(particleEnergies[t]);
			}
		}

		sendCellParametersToBackReceiveFromFront(chargeDensity, outBufferBack, tempCellParameterFront, inBufferFront,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, backRank, frontRank);
		sendCellParametersToFrontReceiveFromBack(chargeDensity, outBufferFront, tempCellParameterBack, inBufferBack,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, backRank, frontRank);

		/*sendCellParametersBack(chargeDensity, outBufferBack, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, backRank);
		receiveCellParametersFront(tempCellParameterFront, inBufferFront, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, frontRank);
		sendCellParametersFront(chargeDensity, outBufferFront, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, frontRank);
		receiveCellParametersBack(tempCellParameterBack, inBufferBack, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, backRank);*/

		sumCellTempParametersY(chargeDensity);

		delete[] inBufferFront;
		delete[] outBufferFront;
		delete[] inBufferBack;
		delete[] outBufferBack;
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int j = 0; j < 2 * additionalBinNumber + 2; ++j) {
					chargeDensity[i][ynumberAdded - 1 - j][k] += chargeDensity[i][1 + 2 * additionalBinNumber - j][k];
					chargeDensity[i][1 + 2 * additionalBinNumber - j][k] = chargeDensity[i][ynumberAdded - 1 - j][k];
					for (int t = 0; t < typesNumber; ++t) {
						if (types[t].particlesPerBin > 0) {
							particleConcentrations[t][i][ynumberAdded - 1 - j][k] += particleConcentrations[t][i][1 + 2 * additionalBinNumber
								- j][k];
							particleConcentrations[t][i][1 + 2 * additionalBinNumber - j][k] = particleConcentrations[t][i][ynumberAdded - 1
								- j][k];

							particleEnergies[t][i][ynumberAdded - 1 - j][k] += particleEnergies[t][i][1 + 2 * additionalBinNumber - j][k];
							particleEnergies[t][i][1 + 2 * additionalBinNumber - j][k] = particleEnergies[t][i][ynumberAdded - 1 - j][k];
						}
					}
				}
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded; ++j) {
						chargeDensity[i][j][k] = chargeDensity[i][1 + additionalBinNumber][k];
						for (int t = 0; t < typesNumber; ++t) {
							if (types[t].particlesPerBin > 0) {
								particleConcentrations[t][i][j][k] = particleConcentrations[t][i][1 + additionalBinNumber][k];
								particleEnergies[t][i][j][k] = particleEnergies[t][i][1 + additionalBinNumber][k];
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumCellParametersZ() {
	if (cartDim[2] > 1) {
		double* inBufferTop = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * xnumberAdded];
		double* outBufferTop = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * xnumberAdded];
		double* inBufferBottom = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * xnumberAdded];
		double* outBufferBottom = new double[(2 + 2 * additionalBinNumber) * ynumberAdded * xnumberAdded];

		for (int t = 0; t < typesNumber; ++t) {
			if (types[t].particlesPerBin > 0) {
				sendCellParametersToBottomReceiveFromTop(particleConcentrations[t], outBufferBottom, tempCellParameterTop,
				                                         inBufferTop, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber,
				                                         cartComm, rank, bottomRank, topRank);
				sendCellParametersToTopReceiveFromBottom(particleConcentrations[t], outBufferTop, tempCellParameterBottom,
				                                         inBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
				                                         additionalBinNumber, cartComm, rank, bottomRank, topRank);

				/*sendCellParametersBottom(particleConcentrations[t], outBufferBottom, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, bottomRank);
				receiveCellParametersTop(tempCellParameterTop, inBufferTop, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, topRank);
				sendCellParametersTop(particleConcentrations[t], outBufferTop, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, topRank);
				receiveCellParametersBottom(tempCellParameterBottom, inBufferBottom, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, bottomRank);*/

				sumCellTempParametersZ(particleConcentrations[t]);
			}
		}

		for (int t = 0; t < typesNumber; ++t) {
			if (types[t].particlesPerBin > 0) {
				sendCellParametersToBottomReceiveFromTop(particleEnergies[t], outBufferBottom, tempCellParameterTop, inBufferTop,
				                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
				                                         rank, bottomRank, topRank);
				sendCellParametersToTopReceiveFromBottom(particleEnergies[t], outBufferTop, tempCellParameterBottom, inBufferBottom,
				                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
				                                         rank, bottomRank, topRank);

				/*sendCellParametersBottom(particleConcentrations[t], outBufferBottom, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, bottomRank);
				receiveCellParametersTop(tempCellParameterTop, inBufferTop, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, topRank);
				sendCellParametersTop(particleConcentrations[t], outBufferTop, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, topRank);
				receiveCellParametersBottom(tempCellParameterBottom, inBufferBottom, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, bottomRank);*/

				sumCellTempParametersZ(particleEnergies[t]);
			}
		}

		sendCellParametersToBottomReceiveFromTop(chargeDensity, outBufferBottom, tempCellParameterTop, inBufferTop,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, bottomRank, topRank);
		sendCellParametersToTopReceiveFromBottom(chargeDensity, outBufferTop, tempCellParameterBottom, inBufferBottom,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, bottomRank, topRank);

		/*sendCellParametersBottom(chargeDensity, outBufferBottom, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, bottomRank);
		receiveCellParametersTop(tempCellParameterTop, inBufferTop, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, topRank);
		sendCellParametersTop(chargeDensity, outBufferTop, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, topRank);
		receiveCellParametersBottom(tempCellParameterBottom, inBufferBottom, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, bottomRank);*/

		sumCellTempParametersZ(chargeDensity);

		delete[] inBufferTop;
		delete[] outBufferTop;
		delete[] inBufferBottom;
		delete[] outBufferBottom;
	} else {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int i = 0; i < xnumberAdded; ++i) {
				for (int k = 0; k < 2 * additionalBinNumber + 2; ++k) {
					chargeDensity[i][j][znumberAdded - 1 - k] += chargeDensity[i][j][1 + 2 * additionalBinNumber - k];
					chargeDensity[i][j][1 + 2 * additionalBinNumber - k] = chargeDensity[i][j][znumberAdded - 1 - k];
					for (int t = 0; t < typesNumber; ++t) {
						particleConcentrations[t][i][j][znumberAdded - 1 - k] += particleConcentrations[t][i][j][1 + 2 *
							additionalBinNumber - k];
						particleConcentrations[t][i][j][1 + 2 * additionalBinNumber - k] = particleConcentrations[t][i][j][znumberAdded -
							1 - k];

						particleEnergies[t][i][j][znumberAdded - 1 - k] += particleEnergies[t][i][j][1 + 2 * additionalBinNumber - k];
						particleEnergies[t][i][j][1 + 2 * additionalBinNumber - k] = particleEnergies[t][i][j][znumberAdded - 1 - k];
					}

				}
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded; ++k) {
						chargeDensity[i][j][k] = chargeDensity[i][j][1 + additionalBinNumber];
						for (int t = 0; t < typesNumber; ++t) {
							if (types[t].particlesPerBin > 0) {
								particleConcentrations[t][i][j][k] = particleConcentrations[t][i][j][1 + additionalBinNumber];
								particleEnergies[t][i][j][k] = particleEnergies[t][i][j][1 + additionalBinNumber];
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumCellTempParametersX(double*** array) {
	if (cartCoord[0] > 0 || boundaryConditionTypeX == PERIODIC) {
		for (int i = 0; i < 2 * additionalBinNumber + 2; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					array[i][j][k] += tempCellParameterLeft[i][j][k];
				}
			}
		}
	}

	if (cartCoord[0] < cartDim[0] - 1 || boundaryConditionTypeX == PERIODIC) {
		for (int i = 0; i < 2 * additionalBinNumber + 2; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					array[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k] += tempCellParameterRight[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumCellTempParametersY(double*** array) {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 * additionalBinNumber + 2; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				array[i][j][k] += tempCellParameterFront[i][j][k];
				array[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k] += tempCellParameterBack[i][j][k];
			}
		}
	}
}

void Simulation::sumCellTempParametersZ(double*** array) {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 * additionalBinNumber + 2; ++k) {
				array[i][j][k] += tempCellParameterBottom[i][j][k];
				array[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k] += tempCellParameterTop[i][j][k];
			}
		}
	}
}

void Simulation::sumCellVectorParametersX() {
	if (cartDim[0] > 1) {
		double* inBufferRight = new double[(2 + 2 * additionalBinNumber) * 3 * ynumberAdded * znumberAdded];
		double* outBufferRight = new double[(2 + 2 * additionalBinNumber) * 3 * ynumberAdded * znumberAdded];
		double* inBufferLeft = new double[(2 + 2 * additionalBinNumber) * 3 * ynumberAdded * znumberAdded];
		double* outBufferLeft = new double[(2 + 2 * additionalBinNumber) * 3 * ynumberAdded * znumberAdded];

		for (int t = 0; t < typesNumber; ++t) {
			//MPI_Barrier(cartComm);
			if (types[t].particlesPerBin > 0) {
				if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
					sendCellVectorParametersToLeftReceiveFromRight(particleBulkVelocities[t], outBufferLeft,
					                                               tempCellVectorParameterRight, inBufferRight, xnumberAdded,
					                                               ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank,
					                                               leftRank, rightRank);
				} else if (cartCoord[0] == 0) {
					receiveCellVectorParametersRight(tempCellVectorParameterRight, inBufferRight, xnumberAdded, ynumberAdded,
					                                 znumberAdded, additionalBinNumber, cartComm, rank, rightRank);
				} else if (cartCoord[0] == cartDim[0] - 1) {
					sendCellVectorParametersLeft(particleBulkVelocities[t], outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
					                             additionalBinNumber, cartComm, rank, leftRank);
				}

				if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
					sendCellVectorParametersToRightReceiveFromLeft(particleBulkVelocities[t], outBufferRight,
					                                               tempCellVectorParameterLeft, inBufferLeft, xnumberAdded,
					                                               ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank,
					                                               leftRank, rightRank);
				} else if (cartCoord[0] == 0) {
					sendCellVectorParametersRight(particleBulkVelocities[t], outBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
					                              additionalBinNumber, cartComm, rank, rightRank);
				} else if (cartCoord[0] == cartDim[0] - 1) {
					receiveCellVectorParametersLeft(tempCellVectorParameterLeft, inBufferLeft, xnumberAdded, ynumberAdded,
					                                znumberAdded, additionalBinNumber, cartComm, rank, leftRank);
				}

				sumCellVectorTempParametersX(particleBulkVelocities[t]);
			}
		}

		if (cartCoord[0] == cartDim[0] - 1 && boundaryConditionTypeX != PERIODIC) {
			//todo!
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						for (int t = 0; t < typesNumber; ++t) {
							particleBulkVelocities[t][xnumberAdded - 2 - additionalBinNumber][j][k] += particleBulkVelocities[t][xnumberAdded
								- 1 - additionalBinNumber + i][j][k];
							particleBulkVelocities[t][xnumberAdded - 1 - additionalBinNumber + i][j][k] = V0;
						}
					}
				}
			}
		}

		delete[] inBufferRight;
		delete[] outBufferRight;
		delete[] inBufferLeft;
		delete[] outBufferLeft;

	} else {
		if (boundaryConditionTypeX == PERIODIC) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i < 2 + 2 * additionalBinNumber; ++i) {
						for (int t = 0; t < typesNumber; ++t) {
							particleBulkVelocities[t][xnumberAdded - 1 - i][j][k] += particleBulkVelocities[t][2 * additionalBinNumber + 1 -
								i][j][k];
							particleBulkVelocities[t][2 * additionalBinNumber + 1 - i][j][k] = particleBulkVelocities[t][xnumberAdded - 1 - i
							][j][k];
						}
					}
				}
			}
		} else {
			//manualy set zero to the last bin
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						for (int t = 0; t < typesNumber; ++t) {
							particleBulkVelocities[t][xnumberAdded - 2 - additionalBinNumber][j][k] += particleBulkVelocities[t][xnumberAdded
								- 1 - additionalBinNumber + i][j][k];
							//particleBulkVelocities[t][xnumberAdded - 1 - additionalBinNumber + i][j][k] = V0;
						}
					}
				}
			}

			if (boundaryConditionTypeX == FREE_BOTH) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						for (int i = 0; i <= additionalBinNumber; ++i) {
							for (int t = 0; t < typesNumber; ++t) {
								particleBulkVelocities[t][xnumberAdded - 2 - additionalBinNumber][j][k] += particleBulkVelocities[t][
									xnumberAdded - 1 - additionalBinNumber + i][j][k];
								particleBulkVelocities[t][xnumberAdded - 1 - additionalBinNumber + i][j][k] = V0;
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumCellVectorParametersY() {
	if (cartDim[1] > 1) {
		double* inBufferFront = new double[(2 + 2 * additionalBinNumber) * 3 * xnumberAdded * znumberAdded];
		double* outBufferFront = new double[(2 + 2 * additionalBinNumber) * 3 * xnumberAdded * znumberAdded];
		double* inBufferBack = new double[(2 + 2 * additionalBinNumber) * 3 * xnumberAdded * znumberAdded];
		double* outBufferBack = new double[(2 + 2 * additionalBinNumber) * 3 * xnumberAdded * znumberAdded];

		for (int t = 0; t < typesNumber; ++t) {
			if (types[t].particlesPerBin > 0) {
				sendCellVectorParametersToBackReceiveFromFront(particleBulkVelocities[t], outBufferBack,
				                                               tempCellVectorParameterFront, inBufferFront, xnumberAdded,
				                                               ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank,
				                                               frontRank, backRank);
				sendCellVectorParametersToFrontReceiveFromBack(particleBulkVelocities[t], outBufferFront,
				                                               tempCellVectorParameterBack, inBufferBack, xnumberAdded,
				                                               ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank,
				                                               frontRank, backRank);

				sumCellVectorTempParametersY(particleBulkVelocities[t]);
			}
		}


		delete[] inBufferFront;
		delete[] outBufferFront;
		delete[] inBufferBack;
		delete[] outBufferBack;
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int j = 0; j < 2 * additionalBinNumber + 2; ++j) {
					for (int t = 0; t < typesNumber; ++t) {
						if (types[t].particlesPerBin > 0) {
							particleBulkVelocities[t][i][ynumberAdded - 1 - j][k] += particleBulkVelocities[t][i][1 + 2 * additionalBinNumber
								- j][k];
							particleBulkVelocities[t][i][1 + 2 * additionalBinNumber - j][k] = particleBulkVelocities[t][i][ynumberAdded - 1
								- j][k];
						}
					}
				}
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded; ++j) {
						for (int t = 0; t < typesNumber; ++t) {
							if (types[t].particlesPerBin > 0) {
								particleBulkVelocities[t][i][j][k] = particleBulkVelocities[t][i][1 + additionalBinNumber][k];
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumCellVectorParametersZ() {
	if (cartDim[2] > 1) {
		double* inBufferTop = new double[(2 + 2 * additionalBinNumber) * 3 * ynumberAdded * xnumberAdded];
		double* outBufferTop = new double[(2 + 2 * additionalBinNumber) * 3 * ynumberAdded * xnumberAdded];
		double* inBufferBottom = new double[(2 + 2 * additionalBinNumber) * 3 * ynumberAdded * xnumberAdded];
		double* outBufferBottom = new double[(2 + 2 * additionalBinNumber) * 3 * ynumberAdded * xnumberAdded];

		for (int t = 0; t < typesNumber; ++t) {
			if (types[t].particlesPerBin > 0) {
				sendCellVectorParametersToBottomReceiveFromTop(particleBulkVelocities[t], outBufferBottom,
				                                               tempCellVectorParameterTop, inBufferTop, xnumberAdded, ynumberAdded,
				                                               znumberAdded, additionalBinNumber, cartComm, rank, bottomRank,
				                                               topRank);
				sendCellVectorParametersToTopReceiveFromBottom(particleBulkVelocities[t], outBufferTop,
				                                               tempCellVectorParameterBottom, inBufferBottom, xnumberAdded,
				                                               ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank,
				                                               bottomRank, topRank);


				sumCellVectorTempParametersZ(particleBulkVelocities[t]);
			}
		}

		delete[] inBufferTop;
		delete[] outBufferTop;
		delete[] inBufferBottom;
		delete[] outBufferBottom;
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 * additionalBinNumber + 2; ++k) {
					for (int t = 0; t < typesNumber; ++t) {
						particleBulkVelocities[t][i][j][znumberAdded - 1 - k] += particleBulkVelocities[t][i][j][1 + 2 *
							additionalBinNumber - k];
						particleBulkVelocities[t][i][j][1 + 2 * additionalBinNumber - k] = particleBulkVelocities[t][i][j][znumberAdded -
							1 - k];
					}
				}
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded; ++k) {
						for (int t = 0; t < typesNumber; ++t) {
							if (types[t].particlesPerBin > 0) {
								particleBulkVelocities[t][i][j][k] = particleBulkVelocities[t][i][j][1 + additionalBinNumber];
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::sumCellVectorTempParametersX(Vector3d*** array) {
	if (cartCoord[0] > 0 || boundaryConditionTypeX == PERIODIC) {
		for (int i = 0; i < 2 * additionalBinNumber + 2; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					array[i][j][k] += tempCellVectorParameterLeft[i][j][k];
				}
			}
		}
	}

	if (cartCoord[0] < cartDim[0] - 1 || boundaryConditionTypeX == PERIODIC) {
		for (int i = 0; i < 2 * additionalBinNumber + 2; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					array[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k] += tempCellVectorParameterRight[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumCellVectorTempParametersY(Vector3d*** array) {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 * additionalBinNumber + 2; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				array[i][j][k] += tempCellVectorParameterFront[i][j][k];
				array[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k] += tempCellVectorParameterBack[i][j][k];
			}
		}
	}
}

void Simulation::sumCellVectorTempParametersZ(Vector3d*** array) {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 * additionalBinNumber + 2; ++k) {
				array[i][j][k] += tempCellVectorParameterBottom[i][j][k];
				array[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k] += tempCellVectorParameterTop[i][j][k];
			}
		}
	}
}

void Simulation::sumCellMatrixParametersX() {
	if (cartDim[0] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum cell matrix parameters rank = %d\n", rank);
		double* inBufferRight = new double[(2 + 2 * additionalBinNumber) * 9 * ynumberAdded * znumberAdded];
		double* outBufferRight = new double[(2 + 2 * additionalBinNumber) * 9 * ynumberAdded * znumberAdded];
		double* inBufferLeft = new double[(2 + 2 * additionalBinNumber) * 9 * ynumberAdded * znumberAdded];
		double* outBufferLeft = new double[(2 + 2 * additionalBinNumber) * 9 * ynumberAdded * znumberAdded];

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("sending left  pressure tensor sum cell matrix parameters rank = %d\n", rank);

		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendCellMatrixParametersToLeftReceiveFromRight(pressureTensor, outBufferLeft, tempCellMatrixParameterRight,
			                                               inBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
			                                               additionalBinNumber, cartComm, rank, leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			receiveCellMatrixParametersRight(tempCellMatrixParameterRight, inBufferRight, xnumberAdded, ynumberAdded,
			                                 znumberAdded, additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			sendCellMatrixParametersLeft(pressureTensor, outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
			                             additionalBinNumber, cartComm, rank, leftRank);
		}

		if ((verbosity > 2)) printf("sending right pressure tensor sum cell matrix parameters rank = %d\n", rank);
		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendCellMatrixParametersToRightReceiveFromLeft(pressureTensor, outBufferRight, tempCellMatrixParameterLeft,
			                                               inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
			                                               additionalBinNumber, cartComm, rank, leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			sendCellMatrixParametersRight(pressureTensor, outBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
			                              additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			receiveCellMatrixParametersLeft(tempCellMatrixParameterLeft, inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
			                                additionalBinNumber, cartComm, rank, leftRank);
		}

		sumCellTempMatrixParametersX(pressureTensor);
		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum cell matrix parameters rank = %d\n", rank);
		delete[] inBufferRight;
		delete[] outBufferRight;
		delete[] inBufferLeft;
		delete[] outBufferLeft;
	} else {
		if (boundaryConditionTypeX == PERIODIC) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i < additionalBinNumber + 2; ++i) {
						pressureTensor[xnumberAdded - 1 - i][j][k] = pressureTensor[xnumberAdded - 1 - i][j][k] + pressureTensor[
							additionalBinNumber + 1 - i][j][k];
						pressureTensor[additionalBinNumber + 1 - i][j][k] = pressureTensor[xnumberAdded - 1 - i][j][k];
					}
				}
			}
		}
	}
}

void Simulation::sumCellMatrixParametersY() {
	if (cartDim[1] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum cell matrix parameters rank = %d\n", rank);
		double* inBufferFront = new double[(2 + 2 * additionalBinNumber) * 9 * xnumberAdded * znumberAdded];
		double* outBufferFront = new double[(2 + 2 * additionalBinNumber) * 9 * xnumberAdded * znumberAdded];
		double* inBufferBack = new double[(2 + 2 * additionalBinNumber) * 9 * xnumberAdded * znumberAdded];
		double* outBufferBack = new double[(2 + 2 * additionalBinNumber) * 9 * xnumberAdded * znumberAdded];

		sendCellMatrixParametersToBackReceiveFromFront(pressureTensor, outBufferBack, tempCellMatrixParameterFront,
		                                               inBufferFront, xnumberAdded, ynumberAdded, znumberAdded,
		                                               additionalBinNumber, cartComm, rank, frontRank, backRank);
		sendCellMatrixParametersToFrontReceiveFromBack(pressureTensor, outBufferFront, tempCellMatrixParameterBack,
		                                               inBufferBack, xnumberAdded, ynumberAdded, znumberAdded,
		                                               additionalBinNumber, cartComm, rank, frontRank, backRank);

		sumCellTempMatrixParametersY(pressureTensor);
		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum cell matrix parameters rank = %d\n", rank);
		delete[] inBufferFront;
		delete[] outBufferFront;
		delete[] inBufferBack;
		delete[] outBufferBack;
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int j = 0; j < additionalBinNumber + 2; ++j) {
					pressureTensor[i][ynumberAdded - 1 - j][k] = pressureTensor[i][ynumberAdded - 1 - j][k] + pressureTensor[i][
						additionalBinNumber + 1 - j][k];
					pressureTensor[i][additionalBinNumber + 1 - j][k] = pressureTensor[i][ynumberAdded - 1 - j][k];
				}
			}
		}
		if (ynumberGeneral == 1) {
			for (int i = 0; i < xnumberAdded; ++i) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int j = 0; j < ynumberAdded; ++j) {
						pressureTensor[i][j][k] = pressureTensor[i][1 + additionalBinNumber][k];
					}
				}
			}
		}
	}
}

void Simulation::sumCellMatrixParametersZ() {
	if (cartDim[2] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum cell matrix parameters rank = %d\n", rank);
		double* inBufferBottom = new double[(2 + 2 * additionalBinNumber) * 9 * xnumberAdded * ynumberAdded];
		double* outBufferBottom = new double[(2 + 2 * additionalBinNumber) * 9 * xnumberAdded * ynumberAdded];
		double* inBufferTop = new double[(2 + 2 * additionalBinNumber) * 9 * xnumberAdded * ynumberAdded];
		double* outBufferTop = new double[(2 + 2 * additionalBinNumber) * 9 * xnumberAdded * ynumberAdded];

		sendCellMatrixParametersToTopReceiveFromBottom(pressureTensor, outBufferTop, tempCellMatrixParameterBottom,
		                                               inBufferBottom, xnumberAdded, ynumberAdded, znumberAdded,
		                                               additionalBinNumber, cartComm, rank, bottomRank, topRank);
		sendCellMatrixParametersToBottomReceiveFromTop(pressureTensor, outBufferBottom, tempCellMatrixParameterTop,
		                                               inBufferTop, xnumberAdded, ynumberAdded, znumberAdded,
		                                               additionalBinNumber, cartComm, rank, bottomRank, topRank);

		sumCellTempMatrixParametersZ(pressureTensor);
		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum cell matrix parameters rank = %d\n", rank);
		delete[] inBufferBottom;
		delete[] outBufferBottom;
		delete[] inBufferTop;
		delete[] outBufferTop;
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < additionalBinNumber + 2; ++k) {
					pressureTensor[i][j][znumberAdded - 1 - k] = pressureTensor[i][j][znumberAdded - 1 - k] + pressureTensor[i][j][
						additionalBinNumber + 1 - k];
					pressureTensor[i][j][additionalBinNumber + 1 - k] = pressureTensor[i][j][znumberAdded - 1 - k];
				}
			}
		}
		if (znumberGeneral == 1) {
			for (int i = 0; i < xnumberAdded; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						pressureTensor[i][j][k] = pressureTensor[i][j][1 + additionalBinNumber];
					}
				}
			}
		}

	}
}

void Simulation::sumCellTempMatrixParametersX(Matrix3d*** array) {
	if (cartCoord[0] > 0 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int i = 0; i < 2 * additionalBinNumber + 2; ++i) {
					array[i][j][k] += tempCellMatrixParameterLeft[i][j][k];
				}
			}
		}
	}

	if (cartCoord[0] < cartDim[0] - 1 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int i = 0; i < 2 * additionalBinNumber + 2; ++i) {
					array[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k] += tempCellMatrixParameterRight[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumCellTempMatrixParametersY(Matrix3d*** array) {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int k = 0; k < znumberAdded; ++k) {
			for (int j = 0; j < 2 * additionalBinNumber + 2; ++j) {
				array[i][j][k] += tempCellMatrixParameterFront[i][j][k];
				array[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k] += tempCellMatrixParameterBack[i][j][k];
			}
		}
	}
}

void Simulation::sumCellTempMatrixParametersZ(Matrix3d*** array) {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 * additionalBinNumber + 2; ++k) {
				array[i][j][k] += tempCellMatrixParameterBottom[i][j][k];
				array[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k] += tempCellMatrixParameterTop[i][j][k];
			}
		}
	}
}

void Simulation::exchangeBunemanFlux() {
	sumBunemanJxAlongX();
	sumBunemanJxAlongY();
	sumBunemanJxAlongZ();

	sumBunemanJyAlongX();
	sumBunemanJyAlongY();
	sumBunemanJyAlongZ();

	sumBunemanJzAlongX();
	sumBunemanJzAlongY();
	sumBunemanJzAlongZ();
}

void Simulation::sumBunemanJxAlongX() {
	if (cartDim[0] > 1) {
		double* inBufferRight = new double[(2 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1)];
		double* outBufferRight = new double[(2 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1)];
		double* inBufferLeft = new double[(2 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1)];
		double* outBufferLeft = new double[(2 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1)];

		if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendCellParametersToLeftReceiveFromRight(bunemanJx, outBufferLeft, tempBunemanJxRight, inBufferRight, xnumberAdded,
			                                         ynumberAdded + 1, znumberAdded + 1, additionalBinNumber, cartComm, rank,
			                                         leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			receiveCellParametersRight(tempBunemanJxRight, inBufferRight, xnumberAdded, ynumberAdded + 1, znumberAdded + 1,
			                           additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			sendCellParametersLeft(bunemanJx, outBufferLeft, xnumberAdded, ynumberAdded + 1, znumberAdded + 1,
			                       additionalBinNumber, cartComm, rank, leftRank);
		}

		if (boundaryConditionTypeX == PERIODIC || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendCellParametersToRightReceiveFromLeft(bunemanJx, outBufferRight, tempBunemanJxLeft, inBufferLeft, xnumberAdded,
			                                         ynumberAdded + 1, znumberAdded + 1, additionalBinNumber, cartComm, rank,
			                                         leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			sendCellParametersRight(bunemanJx, outBufferRight, xnumberAdded, ynumberAdded + 1, znumberAdded + 1,
			                        additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			receiveCellParametersLeft(tempBunemanJxLeft, inBufferLeft, xnumberAdded, ynumberAdded + 1, znumberAdded + 1,
			                          additionalBinNumber, cartComm, rank, leftRank);
		}

		sumTempBunemanJxAlongX();


		delete[] inBufferRight;
		delete[] outBufferRight;
		delete[] inBufferLeft;
		delete[] outBufferLeft;

	} else {
		if (boundaryConditionTypeX == PERIODIC) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i < 2 + 2 * additionalBinNumber; ++i) {
						bunemanJx[xnumberAdded - 1 - i][j][k] += bunemanJx[2 * additionalBinNumber + 1 - i][j][k];
						bunemanJx[2 * additionalBinNumber + 1 - i][j][k] = bunemanJx[xnumberAdded - 1 - i][j][k];
					}
				}
			}
		} else {
			//manualy set zero to the last bin
			/*for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						//todo!!!
						bunemanJx[xnumberAdded - 2 - additionalBinNumber][j][k] += bunemanJx[xnumberAdded - 1 - additionalBinNumber + i][j][k];
						bunemanJx[xnumberAdded - 1 - additionalBinNumber + i][j][k] = 0;
					}
				}
			}

			if (boundaryConditionTypeX == FREE_BOTH) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						for (int i = 0; i <= additionalBinNumber; ++i) {
							for (int t = 0; t < typesNumber; ++t) {
								bunemanJx[xnumberAdded - 2 - additionalBinNumber][j][k] += bunemanJx[xnumberAdded - 1 - additionalBinNumber + i][j][k];
								bunemanJx[xnumberAdded - 1 - additionalBinNumber + i][j][k] = 0;
							}
						}
					}
				}
			}*/
		}
	}
}

void Simulation::sumBunemanJxAlongY() {
	if (cartDim[1] > 1) {
		double* inBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded) * (znumberAdded + 1)];
		double* outBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded) * (znumberAdded + 1)];
		double* inBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded) * (znumberAdded + 1)];
		double* outBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded) * (znumberAdded + 1)];

		if ((verbosity > 2)) printf("sending left flux sum node vector parameters y rank = %d\n", rank);

		sendNodeParametersToFrontReceiveFromBack(bunemanJx, outBufferFront, tempBunemanJxBack, inBufferBack, xnumberAdded - 1,
		                                         ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, frontRank,
		                                         backRank);
		sendNodeParametersToBackReceiveFromFront(bunemanJx, outBufferBack, tempBunemanJxFront, inBufferFront,
		                                         xnumberAdded - 1, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, frontRank, backRank);

		sumTempBunemanJxAlongY();


		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters y rank = %d\n", rank);

		delete[] inBufferFront;
		delete[] inBufferBack;
		delete[] outBufferBack;
		delete[] outBufferFront;

	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k <= znumberAdded; ++k) {
				for (int j = 0; j <= 2 * additionalBinNumber + 2; ++j) {
					bunemanJx[i][ynumberAdded - j][k] += bunemanJx[i][2 + 2 * additionalBinNumber - j][k];
					bunemanJx[i][2 + 2 * additionalBinNumber - j][k] = bunemanJx[i][ynumberAdded - j][k];
				}
			}
		}

		if (ynumberGeneral == 1) {
			for (int i = 0; i < xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= ynumberAdded; ++j) {
						bunemanJx[i][j][k] = bunemanJx[i][1 + additionalBinNumber][k];
					}
				}
			}
		}
	}
}

void Simulation::sumBunemanJxAlongZ() {
	if (cartDim[2] > 1) {
		double* inBufferTop = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded)];
		double* outBufferTop = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded)];
		double* inBufferBottom = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded)];
		double* outBufferBottom = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded)];

		if ((verbosity > 2)) printf("sending left flux sum node vector parameters z rank = %d\n", rank);

		sendNodeParametersToBottomReceiveFromTop(bunemanJx, outBufferBottom, tempBunemanJxTop, inBufferTop, xnumberAdded - 1,
		                                         ynumberAdded, znumberAdded, additionalBinNumber, cartComm, rank, bottomRank,
		                                         topRank);
		sendNodeParametersToTopReceiveFromBottom(bunemanJx, outBufferTop, tempBunemanJxBottom, inBufferBottom,
		                                         xnumberAdded - 1, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, bottomRank, topRank);

		sumTempBunemanJxAlongZ();

		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters z rank = %d\n", rank);

		delete[] inBufferTop;
		delete[] inBufferBottom;
		delete[] outBufferBottom;
		delete[] outBufferTop;

	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= ynumberAdded; ++j) {
				for (int k = 0; k <= 2 * additionalBinNumber + 2; ++k) {
					bunemanJx[i][j][znumberAdded - k] += bunemanJx[i][j][2 + 2 * additionalBinNumber - k];
					bunemanJx[i][j][2 + 2 * additionalBinNumber - k] = bunemanJx[i][j][znumberAdded - k];
				}
			}
		}
		if (znumberGeneral == 1) {
			for (int i = 0; i < xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= ynumberAdded; ++j) {
						bunemanJx[i][j][k] = bunemanJx[i][j][1 + additionalBinNumber];
					}
				}
			}
		}
	}
}

void Simulation::sumBunemanJyAlongX() {
	if (cartDim[0] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum node vector parameters x\n");
		double* inBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded) * (znumberAdded + 1)];
		double* outBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded) * (znumberAdded + 1)];
		double* inBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded) * (znumberAdded + 1)];
		double* outBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded) * (znumberAdded + 1)];

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("sending left flux sum node vector parameters x rank = %d\n", rank);
		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendNodeParametersToLeftReceiveFromRight(bunemanJy, outBufferLeft, tempBunemanJyRight, inBufferRight, xnumberAdded,
			                                         ynumberAdded - 1, znumberAdded, additionalBinNumber, cartComm, rank,
			                                         leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			receiveNodeParametersRight(tempBunemanJyRight, inBufferRight, xnumberAdded, ynumberAdded - 1, znumberAdded,
			                           additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			sendNodeParametersLeft(bunemanJy, outBufferLeft, xnumberAdded, ynumberAdded - 1, znumberAdded, additionalBinNumber,
			                       cartComm, rank, leftRank);
		}

		//MPI_Barrier(cartComm);

		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendNodeParametersToRightReceiveFromLeft(bunemanJy, outBufferRight, tempBunemanJyLeft, inBufferLeft, xnumberAdded,
			                                         ynumberAdded - 1, znumberAdded, additionalBinNumber, cartComm, rank,
			                                         leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			sendNodeParametersRight(bunemanJy, outBufferRight, xnumberAdded, ynumberAdded - 1, znumberAdded, additionalBinNumber,
			                        cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			receiveNodeParametersLeft(tempBunemanJyLeft, inBufferLeft, xnumberAdded, ynumberAdded - 1, znumberAdded,
			                          additionalBinNumber, cartComm, rank, leftRank);
		}


		sumTempBunemanJyAlongX();


		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters x rank = %d\n", rank);
		delete[] inBufferLeft;
		delete[] inBufferRight;
		delete[] outBufferLeft;
		delete[] outBufferRight;
	} else {
		if (boundaryConditionTypeX == PERIODIC) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int i = 0; i <= 2 * additionalBinNumber + 2; ++i) {
						bunemanJy[xnumberAdded - i][j][k] += bunemanJy[2 + 2 * additionalBinNumber - i][j][k];
						bunemanJy[2 + 2 * additionalBinNumber - i][j][k] = bunemanJy[xnumberAdded - i][j][k];
					}
				}
			}
		}
	}
}

void Simulation::sumBunemanJyAlongY() {
	if (cartDim[1] > 1) {
		double* inBufferBack = new double[(2 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1)];
		double* outBufferBack = new double[(2 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1)];
		double* inBufferFront = new double[(2 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1)];
		double* outBufferFront = new double[(2 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1)];

		sendCellParametersToFrontReceiveFromBack(bunemanJy, outBufferFront, tempBunemanJyBack, inBufferBack, xnumberAdded + 1,
		                                         ynumberAdded, znumberAdded + 1, additionalBinNumber, cartComm, rank,
		                                         frontRank, backRank);
		sendCellParametersToBackReceiveFromFront(bunemanJy, outBufferBack, tempBunemanJyFront, inBufferFront,
		                                         xnumberAdded + 1, ynumberAdded, znumberAdded + 1, additionalBinNumber,
		                                         cartComm, rank, frontRank, backRank);


		sumTempBunemanJyAlongY();

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum charge density hat y rank = %d\n", rank);
		delete[] inBufferFront;
		delete[] outBufferFront;
		delete[] inBufferBack;
		delete[] outBufferBack;
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int j = 0; j < 2 * additionalBinNumber + 2; ++j) {
					bunemanJy[i][ynumberAdded - 1 - j][k] += bunemanJy[i][1 + 2 * additionalBinNumber - j][k];
					bunemanJy[i][1 + 2 * additionalBinNumber - j][k] = bunemanJy[i][ynumberAdded - 1 - j][k];
				}
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded; ++j) {
						bunemanJy[i][j][k] = bunemanJy[i][1 + additionalBinNumber][k];
					}
				}
			}
		}
	}
}

void Simulation::sumBunemanJyAlongZ() {
	if (cartDim[2] > 1) {
		double* inBufferTop = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded) * (xnumberAdded + 1)];
		double* outBufferTop = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded) * (xnumberAdded + 1)];
		double* inBufferBottom = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded) * (xnumberAdded + 1)];
		double* outBufferBottom = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded) * (xnumberAdded + 1)];

		if ((verbosity > 2)) printf("sending left flux sum node vector parameters z rank = %d\n", rank);

		sendNodeParametersToBottomReceiveFromTop(bunemanJy, outBufferBottom, tempBunemanJyTop, inBufferTop, xnumberAdded,
		                                         ynumberAdded - 1, znumberAdded, additionalBinNumber, cartComm, rank,
		                                         bottomRank, topRank);
		sendNodeParametersToTopReceiveFromBottom(bunemanJy, outBufferTop, tempBunemanJyBottom, inBufferBottom, xnumberAdded,
		                                         ynumberAdded - 1, znumberAdded, additionalBinNumber, cartComm, rank,
		                                         bottomRank, topRank);

		sumTempBunemanJyAlongZ();

		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters z rank = %d\n", rank);

		delete[] inBufferTop;
		delete[] inBufferBottom;
		delete[] outBufferBottom;
		delete[] outBufferTop;

	} else {
		for (int i = 0; i <= xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= 2 * additionalBinNumber + 2; ++k) {
					bunemanJy[i][j][znumberAdded - k] += bunemanJy[i][j][2 + 2 * additionalBinNumber - k];
					bunemanJy[i][j][2 + 2 * additionalBinNumber - k] = bunemanJy[i][j][znumberAdded - k];
				}
			}
		}
		if (znumberGeneral == 1) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j < ynumberAdded; ++j) {
						bunemanJy[i][j][k] = bunemanJy[i][j][1 + additionalBinNumber];
					}
				}
			}
		}
	}
}

void Simulation::sumBunemanJzAlongX() {
	if (cartDim[0] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum node vector parameters x\n");
		double* inBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded)];
		double* outBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded)];
		double* inBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded)];
		double* outBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded)];

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("sending left flux sum node vector parameters x rank = %d\n", rank);
		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendNodeParametersToLeftReceiveFromRight(bunemanJz, outBufferLeft, tempBunemanJzRight, inBufferRight, xnumberAdded,
			                                         ynumberAdded, znumberAdded - 1, additionalBinNumber, cartComm, rank,
			                                         leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			receiveNodeParametersRight(tempBunemanJzRight, inBufferRight, xnumberAdded, ynumberAdded, znumberAdded - 1,
			                           additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			sendNodeParametersLeft(bunemanJz, outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded - 1, additionalBinNumber,
			                       cartComm, rank, leftRank);
		}

		//MPI_Barrier(cartComm);

		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendNodeParametersToRightReceiveFromLeft(bunemanJz, outBufferRight, tempBunemanJzLeft, inBufferLeft, xnumberAdded,
			                                         ynumberAdded, znumberAdded - 1, additionalBinNumber, cartComm, rank,
			                                         leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			sendNodeParametersRight(bunemanJz, outBufferRight, xnumberAdded, ynumberAdded, znumberAdded - 1, additionalBinNumber,
			                        cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			receiveNodeParametersLeft(tempBunemanJzLeft, inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded - 1,
			                          additionalBinNumber, cartComm, rank, leftRank);
		}


		sumTempBunemanJzAlongX();
		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters x rank = %d\n", rank);
		delete[] inBufferLeft;
		delete[] inBufferRight;
		delete[] outBufferLeft;
		delete[] outBufferRight;
	} else {
		if (boundaryConditionTypeX == PERIODIC) {
			for (int j = 0; j <= ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= 2 * additionalBinNumber + 2; ++i) {
						bunemanJz[xnumberAdded - i][j][k] += bunemanJz[2 + 2 * additionalBinNumber - i][j][k];
						bunemanJz[2 + 2 * additionalBinNumber - i][j][k] = bunemanJz[xnumberAdded - i][j][k];
					}
				}
			}
		}
	}
}

void Simulation::sumBunemanJzAlongY() {
	if (cartDim[1] > 1) {
		double* inBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded) ];
		double* outBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded)];
		double* inBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded)];
		double* outBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded)];

		if ((verbosity > 2)) printf("sending left flux sum node vector parameters y rank = %d\n", rank);

		sendNodeParametersToFrontReceiveFromBack(bunemanJz, outBufferFront, tempBunemanJzBack, inBufferBack, xnumberAdded,
		                                         ynumberAdded, znumberAdded - 1, additionalBinNumber, cartComm, rank,
		                                         frontRank, backRank);
		sendNodeParametersToBackReceiveFromFront(bunemanJz, outBufferBack, tempBunemanJzFront, inBufferFront, xnumberAdded,
		                                         ynumberAdded, znumberAdded - 1, additionalBinNumber, cartComm, rank,
		                                         frontRank, backRank);

		sumTempBunemanJzAlongY();

		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters y rank = %d\n", rank);

		delete[] inBufferFront;
		delete[] inBufferBack;
		delete[] outBufferBack;
		delete[] outBufferFront;
	} else {
		for (int i = 0; i <= xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int j = 0; j <= 2 * additionalBinNumber + 2; ++j) {
					bunemanJz[i][ynumberAdded - j][k] += bunemanJz[i][2 + 2 * additionalBinNumber - j][k];
					bunemanJz[i][2 + 2 * additionalBinNumber - j][k] = bunemanJz[i][ynumberAdded - j][k];
				}
			}
		}

		if (ynumberGeneral == 1) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int j = 0; j <= ynumberAdded; ++j) {
						bunemanJz[i][j][k] = bunemanJz[i][1 + additionalBinNumber][k];
					}
				}
			}
		}
	}
}

void Simulation::sumBunemanJzAlongZ() {
	if (cartDim[2] > 1) {
		double* inBufferTop = new double[(2 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded + 1)];
		double* outBufferTop = new double[(2 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded + 1)];
		double* inBufferBottom = new double[(2 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded + 1)];
		double* outBufferBottom = new double[(2 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded + 1)];

		sendCellParametersToBottomReceiveFromTop(bunemanJz, outBufferBottom, tempBunemanJzTop, inBufferTop, xnumberAdded + 1,
		                                         ynumberAdded + 1, znumberAdded, additionalBinNumber, cartComm, rank,
		                                         bottomRank, topRank);
		sendCellParametersToTopReceiveFromBottom(bunemanJz, outBufferTop, tempBunemanJzBottom, inBufferBottom,
		                                         xnumberAdded + 1, ynumberAdded + 1, znumberAdded, additionalBinNumber,
		                                         cartComm, rank, bottomRank, topRank);

		sumTempBunemanJzAlongZ();

		if ((verbosity > 2)) printf("deleting buffer in sum charge density hat z rank = %d\n", rank);
		delete[] inBufferBottom;
		delete[] outBufferBottom;
		delete[] inBufferTop;
		delete[] outBufferTop;
	} else {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int k = 0; k < 2 * additionalBinNumber + 2; ++k) {
					bunemanJz[i][j][znumberAdded - 1 - k] += bunemanJz[i][j][1 + 2 * additionalBinNumber - k];
					bunemanJz[i][j][1 + 2 * additionalBinNumber - k] = bunemanJz[i][j][znumberAdded - 1 - k];
				}
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded; ++k) {
						bunemanJz[i][j][k] = bunemanJz[i][j][1 + additionalBinNumber];
					}
				}
			}
		}
	}
}

///

void Simulation::sumTempBunemanJxAlongX() {
	if (cartCoord[0] > 0 || boundaryConditionTypeX == PERIODIC) {
		for (int i = 0; i < 2 * additionalBinNumber + 2; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanJx[i][j][k] += tempBunemanJxLeft[i][j][k];
				}
			}
		}
	}

	if (cartCoord[0] < cartDim[0] - 1 || boundaryConditionTypeX == PERIODIC) {
		for (int i = 0; i < 2 * additionalBinNumber + 2; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanJx[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k] += tempBunemanJxRight[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumTempBunemanJxAlongY() {
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int k = 0; k < znumberAdded + 1; ++k) {
			for (int j = 0; j < 2 * additionalBinNumber + 3; ++j) {
				bunemanJx[i][j][k] += tempBunemanJxFront[i][j][k];
				bunemanJx[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k] += tempBunemanJxBack[i][j][k];
			}
		}
	}
}

void Simulation::sumTempBunemanJxAlongZ() {
	for (int j = 0; j < ynumberAdded + 1; ++j) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < 2 * additionalBinNumber + 3; ++k) {
				bunemanJx[i][j][k] += tempBunemanJxBottom[i][j][k];
				bunemanJx[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k] += tempBunemanJxTop[i][j][k];
			}
		}
	}
}

void Simulation::sumTempBunemanJyAlongX() {
	if (cartCoord[0] > 0 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int i = 0; i < 2 * additionalBinNumber + 3; ++i) {
					bunemanJy[i][j][k] += tempBunemanJyLeft[i][j][k];
				}
			}
		}
	}
	//todo
	if (cartCoord[0] < cartDim[0] - 1 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int i = 0; i < 3 + 2 * additionalBinNumber; ++i) {
					bunemanJy[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k] += tempBunemanJyRight[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumTempBunemanJyAlongY() {
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 2 * additionalBinNumber + 2; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanJy[i][j][k] += tempBunemanJyFront[i][j][k];
				bunemanJy[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k] += tempBunemanJyBack[i][j][k];
			}
		}
	}
}

void Simulation::sumTempBunemanJyAlongZ() {
	for (int j = 0; j < ynumberAdded; ++j) {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < 2 * additionalBinNumber + 3; ++k) {
				bunemanJy[i][j][k] += tempBunemanJyBottom[i][j][k];
				bunemanJy[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k] += tempBunemanJyTop[i][j][k];
			}
		}
	}
}

void Simulation::sumTempBunemanJzAlongX() {
	if (cartCoord[0] > 0 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int i = 0; i < 2 * additionalBinNumber + 3; ++i) {
					bunemanJz[i][j][k] += tempBunemanJzLeft[i][j][k];
				}
			}
		}
	}
	//todo
	if (cartCoord[0] < cartDim[0] - 1 || boundaryConditionTypeX == PERIODIC) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int i = 0; i < 3 + 2 * additionalBinNumber; ++i) {
					bunemanJz[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k] += tempBunemanJzRight[i][j][k];
				}
			}
		}
	}
}

void Simulation::sumTempBunemanJzAlongY() {
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int k = 0; k < znumberAdded; ++k) {
			for (int j = 0; j < 2 * additionalBinNumber + 3; ++j) {
				bunemanJz[i][j][k] += tempBunemanJzFront[i][j][k];
				bunemanJz[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k] += tempBunemanJzBack[i][j][k];
			}
		}
	}
}

void Simulation::sumTempBunemanJzAlongZ() {
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 2 * additionalBinNumber + 2; ++k) {
				bunemanJz[i][j][k] += tempBunemanJzBottom[i][j][k];
				bunemanJz[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k] += tempBunemanJzTop[i][j][k];
			}
		}
	}
}

void Simulation::sumNodeParametersX() {
	if (cartDim[0] > 1) {
		if ((verbosity > 2)) printf("crating buffer in sum node vector parameters x\n");
		double* inBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1)];
		double* outBufferRight = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1)];
		double* inBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1)];
		double* outBufferLeft = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1)];

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("sending left flux sum node vector parameters x rank = %d\n", rank);
		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendNodeParametersToLeftReceiveFromRight(bunemanChargeDensity, outBufferLeft, tempNodeParameterRight, inBufferRight,
			                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
			                                         rank, leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			receiveNodeParametersRight(tempNodeParameterRight, inBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
			                           additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			sendNodeParametersLeft(bunemanChargeDensity, outBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
			                       additionalBinNumber, cartComm, rank, leftRank);
		}

		//MPI_Barrier(cartComm);

		if ((boundaryConditionTypeX == PERIODIC) || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
			sendNodeParametersToRightReceiveFromLeft(bunemanChargeDensity, outBufferRight, tempNodeParameterLeft, inBufferLeft,
			                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
			                                         rank, leftRank, rightRank);
		} else if (cartCoord[0] == 0) {
			sendNodeParametersRight(bunemanChargeDensity, outBufferRight, xnumberAdded, ynumberAdded, znumberAdded,
			                        additionalBinNumber, cartComm, rank, rightRank);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			receiveNodeParametersLeft(tempNodeParameterLeft, inBufferLeft, xnumberAdded, ynumberAdded, znumberAdded,
			                          additionalBinNumber, cartComm, rank, leftRank);
		}


		sumTempNodeParametersX(bunemanChargeDensity);

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters x rank = %d\n", rank);
		delete[] inBufferLeft;
		delete[] inBufferRight;
		delete[] outBufferLeft;
		delete[] outBufferRight;
	} else {
		if (boundaryConditionTypeX == PERIODIC) {
			for (int j = 0; j <= ynumberAdded; ++j) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int i = 0; i <= 2 * additionalBinNumber + 2; ++i) {
						bunemanChargeDensity[xnumberAdded - i][j][k] += bunemanChargeDensity[2 + 2 * additionalBinNumber - i][j][k];
						bunemanChargeDensity[2 + 2 * additionalBinNumber - i][j][k] = bunemanChargeDensity[xnumberAdded - i][j][k];
					}
				}
			}
		}
	}
}

void Simulation::sumNodeParametersY() {
	if (cartDim[1] > 1) {
		double* inBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1)];
		double* outBufferFront = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1)];
		double* inBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1)];
		double* outBufferBack = new double[(3 + 2 * additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1)];

		if ((verbosity > 2)) printf("sending left flux sum node vector parameters y rank = %d\n", rank);

		sendNodeParametersToFrontReceiveFromBack(bunemanChargeDensity, outBufferFront, tempNodeParameterBack, inBufferBack,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, frontRank, backRank);
		sendNodeParametersToBackReceiveFromFront(bunemanChargeDensity, outBufferBack, tempNodeParameterFront, inBufferFront,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, frontRank, backRank);

		sumTempNodeParametersY(bunemanChargeDensity);

		//MPI_Barrier(cartComm);
		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters y rank = %d\n", rank);

		delete[] inBufferFront;
		delete[] inBufferBack;
		delete[] outBufferBack;
		delete[] outBufferFront;

	} else {
		for (int i = 0; i <= xnumberAdded; ++i) {
			for (int k = 0; k <= znumberAdded; ++k) {
				for (int j = 0; j <= 2 * additionalBinNumber + 2; ++j) {
					bunemanChargeDensity[i][ynumberAdded - j][k] += bunemanChargeDensity[i][2 + 2 * additionalBinNumber - j][k];
					bunemanChargeDensity[i][2 + 2 * additionalBinNumber - j][k] = bunemanChargeDensity[i][ynumberAdded - j][k];
				}
			}
		}

		if (ynumberGeneral == 1) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= ynumberAdded; ++j) {
						bunemanChargeDensity[i][j][k] = bunemanChargeDensity[i][1 + additionalBinNumber][k];
					}
				}
			}
		}
	}
}

void Simulation::sumNodeParametersZ() {
	if (cartDim[2] > 1) {
		double* inBufferTop = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1)];
		double* outBufferTop = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1)];
		double* inBufferBottom = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1)];
		double* outBufferBottom = new double[(3 + 2 * additionalBinNumber) * (ynumberAdded + 1) * (xnumberAdded + 1)];

		if ((verbosity > 2)) printf("sending left flux sum node vector parameters z rank = %d\n", rank);

		sendNodeParametersToBottomReceiveFromTop(bunemanChargeDensity, outBufferBottom, tempNodeParameterTop, inBufferTop,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, bottomRank, topRank);
		sendNodeParametersToTopReceiveFromBottom(bunemanChargeDensity, outBufferTop, tempNodeParameterBottom, inBufferBottom,
		                                         xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber, cartComm,
		                                         rank, bottomRank, topRank);

		sumTempNodeParametersZ(bunemanChargeDensity);

		if ((verbosity > 2)) printf("deleting buffer in sum node vector parameters z rank = %d\n", rank);

		delete[] inBufferTop;
		delete[] inBufferBottom;
		delete[] outBufferBottom;
		delete[] outBufferTop;
	} else {
		for (int i = 0; i <= xnumberAdded; ++i) {
			for (int j = 0; j <= ynumberAdded; ++j) {
				for (int k = 0; k <= 2 * additionalBinNumber + 2; ++k) {
					bunemanChargeDensity[i][j][znumberAdded - k] += bunemanChargeDensity[i][j][2 + 2 * additionalBinNumber - k];
					bunemanChargeDensity[i][j][2 + 2 * additionalBinNumber - k] = bunemanChargeDensity[i][j][znumberAdded - k];
				}
			}
		}
		if (znumberGeneral == 1) {
			for (int i = 0; i <= xnumberAdded; ++i) {
				for (int k = 0; k <= znumberAdded; ++k) {
					for (int j = 0; j <= ynumberAdded; ++j) {
						bunemanChargeDensity[i][j][k] = bunemanChargeDensity[i][j][1 + additionalBinNumber];
					}
				}
			}
		}
	}
}
