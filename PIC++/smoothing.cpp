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
#include "paths.h"

void Simulation::smoothCellParameter(double*** array) {
	int minI = 1 + additionalBinNumber;
	/*if (cartCoord[0] == 0 && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
		minI = 2 + additionalBinNumber;
	}*/
	int maxI = xnumberAdded - additionalBinNumber;
	double smoothingmatrix[3][3][3];
	createSmoothingMatrix(smoothingmatrix);

	for (int i = minI; i < maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				if (boundaryConditionTypeX == PERIODIC || i > 1 + additionalBinNumber && i < xnumberAdded - additionalBinNumber - 1) {
					tempCellParameter[i][j][k] = 0;
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -1; tempJ <= 1; ++tempJ) {
							for (int tempK = -1; tempK <= 1; ++tempK) {
								tempCellParameter[i][j][k] += array[i + tempI][j + tempJ][k + tempK] *smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					tempCellParameter[i][j][k] = array[i][j][k];
					/*for (int tempI = 0; tempI <= 1; ++tempI) {
						for (int tempJ = -1; tempJ <= 1; ++tempJ) {
							for (int tempK = -1; tempK <= 1; ++tempK) {
								tempCellParameter[i][j][k] += smoothingParameter * array[i + tempI][j + tempJ][k + tempK] / 27.0;
							}
						}
					}*/
				} else if (i == xnumberAdded - additionalBinNumber - 1) {
					tempCellParameter[i][j][k] = array[i][j][k];
					/*for (int tempI = -1; tempI <= 0; ++tempI) {
						for (int tempJ = -1; tempJ <= 1; ++tempJ) {
							for (int tempK = -1; tempK <= 1; ++tempK) {
								tempCellParameter[i][j][k] += smoothingParameter * array[i + tempI][j + tempJ][k + tempK] / 27.0;
							}
						}
					}*/
				}
			}
		}
	}

	for (int i = minI; i < maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				array[i][j][k] = tempCellParameter[i][j][k];
			}
		}
	}
}

void Simulation::smoothChargeDensity() {
	smoothCellParameter(chargeDensity);
}

void Simulation::smoothChargeDensityHat() {
	smoothCellParameter(chargeDensityHat);
}

void Simulation::createSmoothingMatrix(double smoothingmatrix[3][3][3]) {
	for(int i = 0; i < 3; ++i) {
		double factorX = smoothingParameter/2;
		if(i == 1) {
			factorX = 1 - smoothingParameter;
		}
		for(int j = 0; j < 3; ++j) {
			double factorY = smoothingParameter/2;
			if(j == 1) {
				factorY = 1 - smoothingParameter;
			}
			for(int k = 0; k < 3; ++k) {
				double factorZ = smoothingParameter/2;
				if(k == 1) {
					factorZ = 1 - smoothingParameter;
				}
				smoothingmatrix[i][j][k] = factorX*factorY*factorZ;
			}
		}
	}
}

void Simulation::smoothVectorNodeParameter(Vector3d*** E) {
	int minI = 1 + additionalBinNumber;
	/*if (cartCoord[0] == 0 && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
		minI = 2 + additionalBinNumber;
	}*/
	int maxI = xnumberAdded - 1 - additionalBinNumber;

	double smoothingmatrix[3][3][3];
	createSmoothingMatrix(smoothingmatrix);

	for (int i = minI; i <= maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber; ++k) {
				if (i > 1 + additionalBinNumber && i < xnumberAdded - 1 - additionalBinNumber) {
					tempNodeVectorParameter[i][j][k] = Vector3d(0, 0, 0);
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -1; tempJ <= 1; ++tempJ) {
							for (int tempK = -1; tempK <= 1; ++tempK) {
								tempNodeVectorParameter[i][j][k] += E[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1+tempI][1+tempJ][1+tempK];
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempNodeVectorParameter[i][j][k] = E[i][j][k];
						/*for (int tempI = 0; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempNodeVectorParameter[i][j][k] += E[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}*/
					} else {
						tempNodeVectorParameter[i][j][k] = Vector3d(0, 0, 0);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempNodeVectorParameter[i][j][k] += E[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1+tempI][1+tempJ][1+tempK];
								}
							}
						}
					}
				} else if (i == xnumberAdded - 1 - additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempNodeVectorParameter[i][j][k] = E[i][j][k];
						/*for (int tempI = -1; tempI <= 0; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempNodeVectorParameter[i][j][k] += E[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}*/
					} else {
						tempNodeVectorParameter[i][j][k] = Vector3d(0, 0, 0);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempNodeVectorParameter[i][j][k] += E[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1+tempI][1+tempJ][1+tempK];
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i <= maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber; ++k) {
				E[i][j][k] = tempNodeVectorParameter[i][j][k];
			}
		}
	}
}

void Simulation::smoothTempEfield() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}

	smoothVectorNodeParameter(tempEfield);

	exchangeGeneralEfield(tempEfield);
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("smoothing new Efield time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::smoothNewEfield() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}

	smoothVectorNodeParameter(newEfield);

	exchangeGeneralEfield(newEfield);
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("smoothing new Efield time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::smoothFlux() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}

	smoothVectorNodeParameter(electricFlux);

	exchangeGeneralEfield(electricFlux);
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("smoothing new Efield time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::smoothVectorCellParameter(Vector3d*** B) {
	int minI = 1 + additionalBinNumber;
	/*if (cartCoord[0] == 0 && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
		minI = 2 + additionalBinNumber;
	}*/
	int maxI = xnumberAdded - 1 - additionalBinNumber;
	double smoothingmatrix[3][3][3];
	createSmoothingMatrix(smoothingmatrix);

	for (int i = minI; i < maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				if (boundaryConditionTypeX == PERIODIC || i > 1 + additionalBinNumber && i < xnumberAdded - additionalBinNumber - 2) {
					tempCellVectorParameter[i][j][k] = Vector3d(0, 0, 0);
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -1; tempJ <= 1; ++tempJ) {
							for (int tempK = -1; tempK <= 1; ++tempK) {
								tempCellVectorParameter[i][j][k] += B[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempCellVectorParameter[i][j][k] = B[i][j][k];
						/*for (int tempI = 0; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempCellVectorParameter[i][j][k] += B[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}*/
					} else {
						tempCellVectorParameter[i][j][k] = Vector3d(0, 0, 0);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempCellVectorParameter[i][j][k] += B[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				} else if (i == xnumberAdded - additionalBinNumber - 2) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempCellVectorParameter[i][j][k] = B[i][j][k];
						/*for (int tempI = -1; tempI <= 0; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempCellVectorParameter[i][j][k] += B[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}*/
					} else {
						tempCellVectorParameter[i][j][k] = Vector3d(0, 0, 0);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempCellVectorParameter[i][j][k] += B[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i < maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				B[i][j][k] = tempCellVectorParameter[i][j][k];
			}
		}
	}
}

void Simulation::smoothBfield() {
	smoothVectorCellParameter(Bfield);
	exchangeGeneralBfield(Bfield);
}

void Simulation::smoothNewBfield() {
	smoothVectorCellParameter(newBfield);
	exchangeGeneralBfield(newBfield);
}

void Simulation::smoothBunemanEfieldGeneral(double*** fieldX, double*** fieldY, double*** fieldZ) {
	int minI = 1 + additionalBinNumber;
	int maxI = xnumberAdded - 1 - additionalBinNumber;

	////Ex
	for (int i = minI; i < maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber; ++k) {
				if (boundaryConditionTypeX == PERIODIC || i > 1 + additionalBinNumber && i < xnumberAdded - additionalBinNumber - 1
				) {
					tempBunemanExParameter[i][j][k] = fieldX[i][j][k] * (1 - smoothingParameter);
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -1; tempJ <= 1; ++tempJ) {
							for (int tempK = -1; tempK <= 1; ++tempK) {
								tempBunemanExParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempBunemanExParameter[i][j][k] = fieldX[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = 0; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanExParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanExParameter[i][j][k] = fieldX[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanExParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				} else if (i == xnumberAdded - additionalBinNumber - 2) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempBunemanExParameter[i][j][k] = fieldX[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = -1; tempI <= 0; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanExParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanExParameter[i][j][k] = fieldX[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanExParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i < maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber; ++k) {
				fieldX[i][j][k] = tempBunemanExParameter[i][j][k];
			}
		}
	}

	///Ey
	minI = 1 + additionalBinNumber;
	maxI = xnumberAdded - 1 - additionalBinNumber;


	for (int i = minI; i <= maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber; ++k) {
				if (i > 1 + additionalBinNumber && i < xnumberAdded - 1 - additionalBinNumber) {
					tempBunemanEyParameter[i][j][k] = fieldY[i][j][k] * (1 - smoothingParameter);
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -1; tempJ <= 1; ++tempJ) {
							for (int tempK = -1; tempK <= 1; ++tempK) {
								tempBunemanEyParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempBunemanEyParameter[i][j][k] = fieldY[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = 0; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanEyParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanEyParameter[i][j][k] = fieldY[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanEyParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				} else if (i == xnumberAdded - 1 - additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempBunemanEyParameter[i][j][k] = fieldY[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = -1; tempI <= 0; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanEyParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanEyParameter[i][j][k] = fieldY[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanEyParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i <= maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber; ++k) {
				fieldY[i][j][k] = tempBunemanEyParameter[i][j][k];
			}
		}
	}

	////Ez
	minI = 1 + additionalBinNumber;
	maxI = xnumberAdded - 1 - additionalBinNumber;


	for (int i = minI; i <= maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				if (i > 1 + additionalBinNumber && i < xnumberAdded - 1 - additionalBinNumber) {
					tempBunemanEzParameter[i][j][k] = fieldZ[i][j][k] * (1 - smoothingParameter);
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -1; tempJ <= 1; ++tempJ) {
							for (int tempK = -1; tempK <= 1; ++tempK) {
								tempBunemanEzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempBunemanEzParameter[i][j][k] = fieldZ[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = 0; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanEzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanEzParameter[i][j][k] = fieldZ[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanEzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				} else if (i == xnumberAdded - 1 - additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempBunemanEzParameter[i][j][k] = fieldZ[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = -1; tempI <= 0; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanEzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanEzParameter[i][j][k] = fieldZ[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanEzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i <= maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				fieldZ[i][j][k] = tempBunemanEzParameter[i][j][k];
			}
		}
	}

	exchangeBunemanEfield(fieldX, fieldY, fieldZ);
}

void Simulation::smoothBunemanBfieldGeneral(double*** fieldX, double*** fieldY, double*** fieldZ) {
	int minI = 1 + additionalBinNumber;
	int maxI = xnumberAdded - 1 - additionalBinNumber;


	for (int i = minI; i <= maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				if (i > 1 + additionalBinNumber && i < xnumberAdded - 1 - additionalBinNumber) {
					tempBunemanBxParameter[i][j][k] = fieldX[i][j][k] * (1 - smoothingParameter);
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -1; tempJ <= 1; ++tempJ) {
							for (int tempK = -1; tempK <= 1; ++tempK) {
								tempBunemanBxParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempBunemanBxParameter[i][j][k] = fieldX[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = 0; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanBxParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanBxParameter[i][j][k] = fieldX[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanBxParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				} else if (i == xnumberAdded - 1 - additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempBunemanBxParameter[i][j][k] = fieldX[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = -1; tempI <= 0; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanBxParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanBxParameter[i][j][k] = fieldX[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanBxParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i <= maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				fieldX[i][j][k] = tempBunemanBxParameter[i][j][k];
			}
		}
	}

	////By
	minI = 1 + additionalBinNumber;
	maxI = xnumberAdded - 1 - additionalBinNumber;


	for (int i = minI; i < maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				if (boundaryConditionTypeX == PERIODIC || i > 1 + additionalBinNumber && i < xnumberAdded - additionalBinNumber - 2
				) {
					tempBunemanByParameter[i][j][k] = fieldY[i][j][k] * (1 - smoothingParameter);
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -1; tempJ <= 1; ++tempJ) {
							for (int tempK = -1; tempK <= 1; ++tempK) {
								tempBunemanByParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempBunemanByParameter[i][j][k] = fieldY[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = 0; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanByParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanByParameter[i][j][k] = fieldY[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanByParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				} else if (i == xnumberAdded - additionalBinNumber - 2) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempBunemanByParameter[i][j][k] = fieldY[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = -1; tempI <= 0; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanByParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanByParameter[i][j][k] = fieldY[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanByParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i < maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				fieldY[i][j][k] = tempBunemanByParameter[i][j][k];
			}
		}
	}

	///Bz
	minI = 1 + additionalBinNumber;
	maxI = xnumberAdded - 1 - additionalBinNumber;


	for (int i = minI; i < maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber; ++k) {
				if (boundaryConditionTypeX == PERIODIC || i > 1 + additionalBinNumber && i < xnumberAdded - additionalBinNumber - 2
				) {
					tempBunemanBzParameter[i][j][k] = fieldZ[i][j][k] * (1 - smoothingParameter);
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -1; tempJ <= 1; ++tempJ) {
							for (int tempK = -1; tempK <= 1; ++tempK) {
								tempBunemanBzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempBunemanBzParameter[i][j][k] = fieldZ[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = 0; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanBzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanBzParameter[i][j][k] = fieldZ[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanBzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				} else if (i == xnumberAdded - additionalBinNumber - 2) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempBunemanBzParameter[i][j][k] = fieldZ[i][j][k] * (1 - smoothingParameter * 2.0 / 3.0);
						for (int tempI = -1; tempI <= 0; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanBzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					} else {
						tempBunemanBzParameter[i][j][k] = fieldZ[i][j][k] * (1 - smoothingParameter);
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempBunemanBzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i < maxI; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber; ++k) {
				fieldZ[i][j][k] = tempBunemanBzParameter[i][j][k];
			}
		}
	}

	exchangeBunemanBfield(fieldX, fieldY, fieldZ);
}
