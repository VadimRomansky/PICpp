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
	int minJ = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - additionalBinNumber - 1;
	int minK = 1 + additionalBinNumber;
	int maxK = znumberAdded - additionalBinNumber - 1;
	if(ynumberGeneral == 1) {
		minJ = 0;
		maxJ = 1;
	}
	if(znumberGeneral == 1) {
		minK = 0;
		maxK = 1;
	}


	double smoothingmatrix[3][3][3];
	createSmoothingMatrix(smoothingmatrix);

	int jwidth = 1;
	if(ynumberGeneral == 1) {
		jwidth = 0;
	}
	int kwidth = 1;
		if(znumberGeneral == 1) {
			kwidth = 0;
	}

	for (int i = minI; i < maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minJ; k < maxK; ++k) {
				if (boundaryConditionTypeX == PERIODIC || ((i > 1 + additionalBinNumber) && (i < xnumberAdded - additionalBinNumber - 1))) {
					tempCellParameter[i][j][k] = 0;
				
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
							for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
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
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				array[i][j][k] = tempCellParameter[i][j][k];
			}
		}
	}
}

void Simulation::smoothChargeDensity() {
	smoothCellParameter(chargeDensity);
	exchangeGeneralScalarCellField(chargeDensity);
}

void Simulation::smoothChargeDensityHat() {
	smoothCellParameter(chargeDensityHat);
	exchangeGeneralScalarCellField(chargeDensityHat);
}

void Simulation::createSmoothingMatrix(double smoothingmatrix[3][3][3]) {
	for(int i = 0; i < 3; ++i) {
		double factorX = smoothingParameter/2;
		if(i == 1) {
			factorX = 1 - smoothingParameter;
		}
		for(int j = 0; j < 3; ++j) {
			double factorY = smoothingParameter/2;
			if(ynumberGeneral == 1){
				factorY = 0;
			}
			if(j == 1) {
				factorY = 1 - smoothingParameter;
				if(ynumberGeneral == 1){
					factorY = 1;
				}
			}
			for(int k = 0; k < 3; ++k) {
				double factorZ = smoothingParameter/2;
				if(znumberGeneral == 1){
					factorZ = 0;
				}
				if(k == 1) {
					factorZ = 1 - smoothingParameter;
					if(znumberGeneral == 1){
						factorZ = 1;
					}
				}
				smoothingmatrix[i][j][k] = factorX*factorY*factorZ;
			}
		}
	}
}

void Simulation::smoothVectorNodeParameter(Vector3d*** E) {
	int minI = 1 + additionalBinNumber;
	if (cartCoord[0] == 0 && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
		minI = 2 + additionalBinNumber;
	}
	int maxI = xnumberAdded - 1 - additionalBinNumber;
	int minJ = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - additionalBinNumber;
	int minK = 1 + additionalBinNumber;
	int maxK = znumberAdded - additionalBinNumber;
	if(ynumberGeneral == 1) {
		minJ = 0;
		maxJ = 1;
	}
	if(znumberGeneral == 1) {
		minK = 0;
		maxK = 1;
	}

	double smoothingmatrix[3][3][3];
	createSmoothingMatrix(smoothingmatrix);

	for (int i = minI; i <= maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				if (i > 1 + additionalBinNumber && i < xnumberAdded - 1 - additionalBinNumber) {
					tempNodeVectorParameter[i][j][k] = Vector3d(0, 0, 0);
					int jwidth = 1;
					if(ynumberGeneral == 1) {
						jwidth = 0;
					}
					int kwidth = 1;
					if(znumberGeneral == 1) {
						kwidth = 0;
					}
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
							for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
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
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
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
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
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
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				E[i][j][k] = tempNodeVectorParameter[i][j][k];
			}
		}
	}
}

void Simulation::smoothMatrixNodeParameter(Matrix3d*** E) {
	int minI = 1 + additionalBinNumber;
	if (cartCoord[0] == 0 && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
		minI = 2 + additionalBinNumber;
	}
	int maxI = xnumberAdded - 1 - additionalBinNumber;
	int minJ = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - additionalBinNumber;
	int minK = 1 + additionalBinNumber;
	int maxK = znumberAdded - additionalBinNumber;
	if(ynumberGeneral == 1) {
		minJ = 0;
		maxJ = 1;
	}
	if(znumberGeneral == 1) {
		minK = 0;
		maxK = 1;
	}
	double smoothingmatrix[3][3][3];
	createSmoothingMatrix(smoothingmatrix);

	for (int i = minI; i <= maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				if (i > 1 + additionalBinNumber && i < xnumberAdded - 1 - additionalBinNumber) {
					tempNodeMatrixParameter[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
					int jwidth = 1;
					if(ynumberGeneral == 1) {
						jwidth = 0;
					}
					int kwidth = 1;
					if(znumberGeneral == 1) {
						kwidth = 0;
					}
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
							for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
								tempNodeMatrixParameter[i][j][k] += E[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1+tempI][1+tempJ][1+tempK];
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempNodeMatrixParameter[i][j][k] = E[i][j][k];
						/*for (int tempI = 0; tempI <= 1; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempNodeVectorParameter[i][j][k] += E[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}*/
					} else {
						tempNodeMatrixParameter[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempNodeMatrixParameter[i][j][k] += E[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1+tempI][1+tempJ][1+tempK];
								}
							}
						}
					}
				} else if (i == xnumberAdded - 1 - additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempNodeMatrixParameter[i][j][k] = E[i][j][k];
						/*for (int tempI = -1; tempI <= 0; ++tempI) {
							for (int tempJ = -1; tempJ <= 1; ++tempJ) {
								for (int tempK = -1; tempK <= 1; ++tempK) {
									tempNodeVectorParameter[i][j][k] += E[i + tempI][j + tempJ][k + tempK] * smoothingParameter / 27.0;
								}
							}
						}*/
					} else {
						tempNodeMatrixParameter[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempNodeMatrixParameter[i][j][k] += E[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1+tempI][1+tempJ][1+tempK];
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i <= maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				for(int l = 0; l < 3; ++l) {
					for(int m = 0; m < 3; ++m){			
						E[i][j][k].matrix[l][m] = tempNodeMatrixParameter[i][j][k].matrix[l][m];
					}
				}
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
	//MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("smoothing new Efield time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::smoothNewEfield() {
	/*double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}*/

	smoothVectorNodeParameter(newEfield);

	exchangeGeneralEfield(newEfield);
	/*MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("smoothing new Efield time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}*/
}

void Simulation::smoothFlux() {
	/*double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}*/

	smoothVectorNodeParameter(electricFlux);

	exchangeGeneralEfield(electricFlux);
	/*MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("smoothing new Efield time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}*/
}

void Simulation::smoothVectorCellParameter(Vector3d*** B) {
	int minI = 1 + additionalBinNumber;
	/*if (cartCoord[0] == 0 && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
		minI = 2 + additionalBinNumber;
	}*/
	int maxI = xnumberAdded - 1 - additionalBinNumber;
	int minJ = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - additionalBinNumber - 1;
	int minK = 1 + additionalBinNumber;
	int maxK = znumberAdded - additionalBinNumber - 1;
	if(ynumberGeneral == 1) {
		minJ = 0;
		maxJ = 1;
	}
	if(znumberGeneral == 1) {
		minK = 0;
		maxK = 1;
	}
	double smoothingmatrix[3][3][3];
	createSmoothingMatrix(smoothingmatrix);

	for (int i = minI; i < maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minJ; k < maxJ; ++k) {
				if (boundaryConditionTypeX == PERIODIC || i > 1 + additionalBinNumber && i < xnumberAdded - additionalBinNumber - 2) {
					tempCellVectorParameter[i][j][k] = Vector3d(0, 0, 0);
					int jwidth = 1;
					if(ynumberGeneral == 1) {
						jwidth = 0;
					}
					int kwidth = 1;
					if(znumberGeneral == 1) {
						kwidth = 0;
					}
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
							for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
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
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
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
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
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
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
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
	int minJ = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - additionalBinNumber - 1 ;
	int minK = 1 + additionalBinNumber;
	int maxK = znumberAdded - additionalBinNumber - 1;
	if(ynumberGeneral == 1) {
		minJ = 0;
		maxJ = 0;
	}
	if(znumberGeneral == 1) {
		minK = 0;
		maxK = 0;
	}

	double smoothingmatrix[3][3][3];
	createSmoothingMatrix(smoothingmatrix);

	int jwidth = 1;
	if(ynumberGeneral == 1) {
		jwidth = 0;
	}
	int kwidth = 1;
	if(znumberGeneral == 1) {
		kwidth = 0;
	}
	////Ex
	for (int i = minI; i < maxI; ++i) {
		for (int j = minJ; j <= maxJ; ++j) {
			for (int k = minK; k <= maxK; ++k) {
				if (boundaryConditionTypeX == PERIODIC || i > 1 + additionalBinNumber && i < xnumberAdded - additionalBinNumber - 1) {
					tempBunemanExParameter[i][j][k] = 0;					
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
							for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
								tempBunemanExParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (cartCoord[0] == 0) {
						tempBunemanExParameter[i][j][k] = fieldX[i][j][k];
					} else {
						tempBunemanExParameter[i][j][k] = 0;
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanExParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				} else if (i == xnumberAdded - additionalBinNumber - 2) {
					if (cartCoord[0] == cartDim[0] - 1) {
						tempBunemanExParameter[i][j][k] = fieldX[i][j][k];
					} else {
						tempBunemanExParameter[i][j][k] = 0;
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanExParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i < maxI; ++i) {
		for (int j = minJ; j <= maxJ; ++j) {
			for (int k = minK; k <= maxK; ++k) {
				fieldX[i][j][k] = tempBunemanExParameter[i][j][k];
			}
		}
	}

	///Ey
	minI = 1 + additionalBinNumber;
	maxI = xnumberAdded - 1 - additionalBinNumber;

	maxJ = ynumberAdded - additionalBinNumber - 1;
	if(ynumberGeneral == 1) {
		maxJ = 1;
	}
	
	for (int i = minI; i <= maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k <= maxK; ++k) {
				if (i > 1 + additionalBinNumber && i < xnumberAdded - 1 - additionalBinNumber) {
					tempBunemanEyParameter[i][j][k] = 0;
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
							for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
								tempBunemanEyParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempBunemanEyParameter[i][j][k] = fieldY[i][j][k];
					} else {
						tempBunemanEyParameter[i][j][k] = 0;
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanEyParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				} else if (i == xnumberAdded - 1 - additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempBunemanEyParameter[i][j][k] = fieldY[i][j][k];
					} else {
						tempBunemanEyParameter[i][j][k] = 0;
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanEyParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i <= maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k <= maxK; ++k) {
				fieldY[i][j][k] = tempBunemanEyParameter[i][j][k];
			}
		}
	}

	////Ez
	minI = 1 + additionalBinNumber;
	maxI = xnumberAdded - 1 - additionalBinNumber;
	maxJ = ynumberAdded - additionalBinNumber - 1;
	if(ynumberGeneral == 1) {
		maxJ = 0;
	}
	maxK = znumberAdded - additionalBinNumber - 1;
	if(znumberGeneral == 1) {
		maxK = 1;
	}

	for (int i = minI; i <= maxI; ++i) {
		for (int j = minJ; j <= maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				if (i > 1 + additionalBinNumber && i < xnumberAdded - 1 - additionalBinNumber) {
					tempBunemanEzParameter[i][j][k] = 0;
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
							for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
								tempBunemanEzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempBunemanEzParameter[i][j][k] = fieldZ[i][j][k];
					} else {
						tempBunemanEzParameter[i][j][k] = 0;
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanEzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				} else if (i == xnumberAdded - 1 - additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempBunemanEzParameter[i][j][k] = fieldZ[i][j][k];
					} else {
						tempBunemanEzParameter[i][j][k] = 0;
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanEzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i <= maxI; ++i) {
		for (int j = minJ; j <= maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				fieldZ[i][j][k] = tempBunemanEzParameter[i][j][k];
			}
		}
	}

	exchangeBunemanEfield(fieldX, fieldY, fieldZ);
}

void Simulation::smoothBunemanBfieldGeneral(double*** fieldX, double*** fieldY, double*** fieldZ) {
	int minI = 1 + additionalBinNumber;
	int maxI = xnumberAdded - 1 - additionalBinNumber;
	int minJ = 1 + additionalBinNumber;
	int maxJ = ynumberAdded - additionalBinNumber - 1;
	int minK = 1 + additionalBinNumber;
	int maxK = znumberAdded - additionalBinNumber - 1;
	if(ynumberGeneral == 1) {
		minJ = 0;
		maxJ = 0;
	}
	if(znumberGeneral == 1) {
		minK = 0;
		maxK = 1;
	}

	double smoothingmatrix[3][3][3];
	createSmoothingMatrix(smoothingmatrix);

	for (int i = minI; i <= maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				if (i > 1 + additionalBinNumber && i < xnumberAdded - 1 - additionalBinNumber) {
					tempBunemanBxParameter[i][j][k] = 0;
					int jwidth = 1;
					if(ynumberGeneral == 1) {
						jwidth = 0;
					}
					int kwidth = 1;
					if(znumberGeneral == 1) {
						kwidth = 0;
					}
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
							for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
								tempBunemanBxParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempBunemanBxParameter[i][j][k] = fieldX[i][j][k];
					} else {
						tempBunemanBxParameter[i][j][k] = 0;
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanBxParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				} else if (i == xnumberAdded - 1 - additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempBunemanBxParameter[i][j][k] = fieldX[i][j][k];
					} else {
						tempBunemanBxParameter[i][j][k] = 0;
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanBxParameter[i][j][k] += fieldX[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i <= maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				fieldX[i][j][k] = tempBunemanBxParameter[i][j][k];
			}
		}
	}

	////By
	minI = 1 + additionalBinNumber;
	maxI = xnumberAdded - 1 - additionalBinNumber;
	
	maxJ = ynumberAdded - additionalBinNumber;
	if(ynumberGeneral == 1) {
		maxJ = 1;
	}


	for (int i = minI; i < maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				if (boundaryConditionTypeX == PERIODIC || i > 1 + additionalBinNumber && i < xnumberAdded - additionalBinNumber - 2
				) {
					tempBunemanByParameter[i][j][k] = 0;
					int jwidth = 1;
					if(ynumberGeneral == 1) {
						jwidth = 0;
					}
					int kwidth = 1;
					if(znumberGeneral == 1) {
						kwidth = 0;
					}
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
							for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
								tempBunemanByParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempBunemanByParameter[i][j][k] = fieldY[i][j][k];
					} else {
						tempBunemanByParameter[i][j][k] = 0;
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanByParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				} else if (i == xnumberAdded - additionalBinNumber - 2) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempBunemanByParameter[i][j][k] = fieldY[i][j][k];
					} else {
						tempBunemanByParameter[i][j][k] = 0;
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanByParameter[i][j][k] += fieldY[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i < maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				fieldY[i][j][k] = tempBunemanByParameter[i][j][k];
			}
		}
	}

	///Bz
	minI = 1 + additionalBinNumber;
	maxI = xnumberAdded - 1 - additionalBinNumber;
	maxJ = ynumberAdded - additionalBinNumber - 1;
	if(ynumberGeneral == 1) {
		maxJ = 1;
	}
	maxK = znumberAdded - additionalBinNumber;
	if(znumberGeneral == 1) {
		maxK = 1;
	}


	for (int i = minI; i < maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				if (boundaryConditionTypeX == PERIODIC || i > 1 + additionalBinNumber && i < xnumberAdded - additionalBinNumber - 2
				) {
					tempBunemanBzParameter[i][j][k] = 0;
					int jwidth = 1;
					if(ynumberGeneral == 1) {
						jwidth = 0;
					}
					int kwidth = 1;
					if(znumberGeneral == 1) {
						kwidth = 0;
					}
					for (int tempI = -1; tempI <= 1; ++tempI) {
						for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
							for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
								tempBunemanBzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
							}
						}
					}
				} else if (i == 1 + additionalBinNumber) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == 0) {
						tempBunemanBzParameter[i][j][k] = fieldZ[i][j][k];
					} else {
						tempBunemanBzParameter[i][j][k] = 0;
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanBzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				} else if (i == xnumberAdded - additionalBinNumber - 2) {
					if (boundaryConditionTypeX != PERIODIC && cartCoord[0] == cartDim[0] - 1) {
						tempBunemanBzParameter[i][j][k] = fieldZ[i][j][k];
					} else {
						tempBunemanBzParameter[i][j][k] = 0;
						int jwidth = 1;
						if(ynumberGeneral == 1) {
							jwidth = 0;
						}
						int kwidth = 1;
						if(znumberGeneral == 1) {
							kwidth = 0;
						}
						for (int tempI = -1; tempI <= 1; ++tempI) {
							for (int tempJ = -jwidth; tempJ <= jwidth; ++tempJ) {
								for (int tempK = -kwidth; tempK <= kwidth; ++tempK) {
									tempBunemanBzParameter[i][j][k] += fieldZ[i + tempI][j + tempJ][k + tempK] * smoothingmatrix[1 + tempI][1 + tempJ][1 + tempK];
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = minI; i < maxI; ++i) {
		for (int j = minJ; j < maxJ; ++j) {
			for (int k = minK; k < maxK; ++k) {
				fieldZ[i][j][k] = tempBunemanBzParameter[i][j][k];
			}
		}
	}

	exchangeBunemanBfield(fieldX, fieldY, fieldZ);
}

void Simulation::fuckingStrangeSmoothingBunemanFields(double*** oldEx, double*** oldEy, double*** oldEz, double*** oldBx, double*** oldBy,
	double*** oldBz, double*** newEx, double*** newEy, double*** newEz, double*** newBx, double*** newBy, double*** newBz) {
	for(int i = 0; i < xnumberAdded; ++i) {
		for(int j = 0; j <= ynumberAdded; ++j) {
			for(int k = 0; k <= znumberAdded; ++k) {
				newEx[i][j][k] = oldEx[i][j][k];
			}
		}
	}
	for(int i = 0; i <= xnumberAdded; ++i) {
		for(int j = 0; j < ynumberAdded; ++j) {
			for(int k = 0; k <= znumberAdded; ++k) {
				newEy[i][j][k] = oldEy[i][j][k];
			}
		}
	}
	for(int i = 0; i <= xnumberAdded; ++i) {
		for(int j = 0; j <= ynumberAdded; ++j) {
			for(int k = 0; k < znumberAdded; ++k) {
				newEz[i][j][k] = oldEz[i][j][k];
			}
		}
	}

	for(int i = 0; i <= xnumberAdded; ++i) {
		for(int j = 0; j < ynumberAdded; ++j) {
			for(int k = 0; k < znumberAdded; ++k) {
				newBx[i][j][k] = oldBx[i][j][k];
			}
		}
	}
	for(int i = 0; i < xnumberAdded; ++i) {
		for(int j = 0; j <= ynumberAdded; ++j) {
			for(int k = 0; k < znumberAdded; ++k) {
				newBy[i][j][k] = oldBy[i][j][k];
			}
		}
	}
	for(int i = 0; i < xnumberAdded; ++i) {
		for(int j = 0; j < ynumberAdded; ++j) {
			for(int k = 0; k <= znumberAdded; ++k) {
				newBz[i][j][k] = oldBz[i][j][k];
			}
		}
	}


	int minI = 1;
	int maxI = xnumberAdded - 1;
	int minJ = 1;
	int maxJ = ynumberAdded - 1;
	int minK = 1;
	int maxK = znumberAdded - 1;

	for(int i = minI; i < maxI; ++i) {
		for(int j = minJ; j < maxJ; ++j) {
			for(int k = 0; k < maxK; ++k) {
				newBx[i][j][k] = (oldEx[i][j][k] + oldEx[i][j+1][k] + oldEx[i][j][k+1] + oldEx[i][j+1][k+1] +
					oldEx[i-1][j][k] + oldEx[i-1][j+1][k] + oldEx[i-1][j][k+1] + oldEx[i-1][j+1][k+1])/8;
			}
		}
	}

	for(int i = minI; i < maxI; ++i) {
		for(int j = minJ; j < maxJ; ++j) {
			for(int k = 0; k < maxK; ++k) {
				newEx[i][j][k] = (oldBx[i][j][k] + oldBx[i][j-1][k] + oldBx[i][j][k-1] + oldBx[i][j-1][k-1] +
					oldBx[i+1][j][k] + oldBx[i+1][j-1][k] + oldBx[i+1][j][k-1] + oldBx[i+1][j-1][k-1])/8;
			}
		}
	}
}
