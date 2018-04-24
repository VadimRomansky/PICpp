//
// Created by vadim on 19.06.16.
//
#include <mpi.h>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "mpi_util.h"
#include "particle.h"
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "paths.h"

void sendInput(Simulation& simulation, int nprocs) {
	int integerData[15];
	double doubleData[26];

	if (simulation.rank == 0) {
		integerData[0] = simulation.xnumberGeneral;
		integerData[1] = simulation.ynumberGeneral;
		integerData[2] = simulation.znumberGeneral;
		integerData[3] = simulation.maxIteration;
		integerData[4] = simulation.particlesPerBin[0];
		integerData[5] = simulation.particlesPerBin[1];
		integerData[6] = simulation.particlesPerBin[2];
		integerData[7] = simulation.particlesPerBin[3];
		integerData[8] = simulation.particlesPerBin[4];
		integerData[9] = simulation.particlesPerBin[5];
		integerData[10] = simulation.particlesPerBin[6];
		integerData[11] = simulation.particlesPerBin[7];
		integerData[12] = simulation.inputType;
		integerData[13] = simulation.verbosity;
		int solverType = -1;
		if (simulation.solverType == IMPLICIT) {
			solverType = 0;
		} else if (simulation.solverType == EXPLICIT) {
			solverType = 1;
		} else if (simulation.solverType == BUNEMAN) {
			solverType = 2;
		} else if (simulation.solverType == IMPLICIT_EC) {
			solverType = 3;
		}
		integerData[14] = solverType;

		doubleData[0] = simulation.xsizeGeneral;
		doubleData[1] = simulation.ysizeGeneral;
		doubleData[2] = simulation.zsizeGeneral;
		doubleData[3] = simulation.temperature;
		doubleData[4] = simulation.concentrations[0];
		doubleData[5] = simulation.concentrations[1];
		doubleData[6] = simulation.concentrations[2];
		doubleData[7] = simulation.concentrations[3];
		doubleData[8] = simulation.concentrations[4];
		doubleData[9] = simulation.concentrations[5];
		doubleData[10] = simulation.concentrations[6];
		doubleData[11] = simulation.concentrations[7];
		doubleData[12] = simulation.V0.x;
		doubleData[13] = simulation.V0.y;
		doubleData[14] = simulation.V0.z;
		doubleData[15] = simulation.E0.x;
		doubleData[16] = simulation.E0.y;
		doubleData[17] = simulation.E0.z;
		doubleData[18] = simulation.B0.x;
		doubleData[19] = simulation.B0.y;
		doubleData[20] = simulation.B0.z;
		doubleData[21] = simulation.maxTime;
		doubleData[22] = simulation.preferedDeltaT;
		doubleData[23] = simulation.electronMassInput;
		doubleData[24] = simulation.plasma_period;
		doubleData[25] = simulation.scaleFactor;

		for (int i = 1; i < nprocs; ++i) {
			// printf("sending input to %d\n", i);
			MPI_Send(integerData, 15, MPI_INT, i, MPI_INPUT_INTEGER_TAG, simulation.cartComm);
			MPI_Send(doubleData, 26, MPI_DOUBLE, i, MPI_INPUT_DOUBLE_TAG, simulation.cartComm);
		}
	} else {
		printf("send input only with 0 rank!\n");
		MPI_Finalize();
		exit(0);
	}
}

Simulation recieveInput(MPI_Comm cartComm) {
	int integerData[15];
	double doubleData[26];

	int rank;
	MPI_Comm_rank(cartComm, &rank);
	int nprocs;
	MPI_Comm_size(cartComm, &nprocs);

	if (rank != 0) {
		//printf("receiving\n");
		MPI_Status status;
		MPI_Recv(integerData, 15, MPI_INT, 0, MPI_INPUT_INTEGER_TAG, cartComm, &status);
		MPI_Recv(doubleData, 26, MPI_DOUBLE, 0, MPI_INPUT_DOUBLE_TAG, cartComm, &status);

		int* particlesPerBin = new int[8];
		double* concentrations = new double[8];

		int xnumber = integerData[0];
		int ynumber = integerData[1];
		int znumber = integerData[2];
		int maxIteration = integerData[3];
		particlesPerBin[0] = integerData[4];
		particlesPerBin[1] = integerData[5];
		particlesPerBin[2] = integerData[6];
		particlesPerBin[3] = integerData[7];
		particlesPerBin[4] = integerData[8];
		particlesPerBin[5] = integerData[9];
		particlesPerBin[6] = integerData[10];
		particlesPerBin[7] = integerData[11];
		int inputType = integerData[12];
		int verbosity = integerData[13];
		int solverType = integerData[14];
		SolverType solverTypev;
		if (solverType == 0) {
			solverTypev = IMPLICIT;
		} else if (solverType == 1) {
			solverTypev = EXPLICIT;
		} else if (solverType == 2) {
			solverTypev = BUNEMAN;
		} else if (solverType == 3) {
			solverTypev = IMPLICIT_EC;
		} else {
			printf("wrong solver type\n");
			MPI_Finalize();
			exit(0);
		}
		double xsize = doubleData[0];
		double ysize = doubleData[1];
		double zsize = doubleData[2];
		double temperature = doubleData[3];
		concentrations[0] = doubleData[4];
		concentrations[1] = doubleData[5];
		concentrations[2] = doubleData[6];
		concentrations[3] = doubleData[7];
		concentrations[4] = doubleData[8];
		concentrations[5] = doubleData[9];
		concentrations[6] = doubleData[10];
		concentrations[7] = doubleData[11];
		double V0x = doubleData[12];
		double V0y = doubleData[13];
		double V0z = doubleData[14];
		double E0x = doubleData[15];
		double E0y = doubleData[16];
		double E0z = doubleData[17];
		double B0x = doubleData[18];
		double B0y = doubleData[19];
		double B0z = doubleData[20];
		double maxTime = doubleData[21];
		double preferedDeltaT = doubleData[22];
		double electronMassInput = doubleData[23];
		double plasmaPeriod = doubleData[24];
		double scaleFactor = doubleData[25];
		return Simulation(xnumber, ynumber, znumber, xsize, ysize, zsize, temperature, V0x, V0y, V0z, E0x, E0y, E0z,
		                  B0x, B0y, B0z, maxIteration, maxTime, 8, particlesPerBin, concentrations, inputType, nprocs,
		                  verbosity, preferedDeltaT, electronMassInput, plasmaPeriod, scaleFactor, solverTypev, cartComm);
	}
	printf("recieve input only with not 0 rank!\n");
	MPI_Finalize();
	exit(0);
}

void sendLargeVectorToRightReceiveFromLeft(double**** tempVector, double* outBuffer, double* inBuffer, int xnumberAdded,
                                           int ynumberAdded, int znumberAdded,
                                           int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[0] > 1) {
		int bcount = 0;
		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outBuffer[bcount] = tempVector[xnumberAdded - 2 - additionalBinNumber - i][j][k][l];
						++bcount;
					}
				}
			}
		}
		int rightCoord[3];
		for (int i = 0; i < 3; ++i) {
			rightCoord[i] = cartCoord[i];
		}
		rightCoord[0] += 1;
		if (rightCoord[0] >= cartDim[0]) {
			rightCoord[0] = 0;
		}
		int rightRank;
		MPI_Cart_rank(cartComm, rightCoord, &rightRank);

		int leftCoord[3];
		for (int i = 0; i < 3; ++i) {
			leftCoord[i] = cartCoord[i];
		}
		leftCoord[0] -= 1;
		if (leftCoord[0] < 0) {
			leftCoord[0] = cartDim[0] - 1;
		}
		int leftRank;
		MPI_Cart_rank(cartComm, leftCoord, &leftRank);


		MPI_Status status;
		int number = ynumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber);
		MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, rightRank, MPI_TEMPVECTOR_RIGHT, inBuffer, number, MPI_DOUBLE, leftRank,
		             MPI_TEMPVECTOR_RIGHT, cartComm, &status);

		bcount = 0;
		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[additionalBinNumber - i][j][k][l] = inBuffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[additionalBinNumber - i][j][k][l] = tempVector[xnumberAdded - 2 - additionalBinNumber - i][j][k][l];
					}
				}
			}
		}
	}
}

void sendLargeVectorToLeftReceiveFromRight(double**** tempVector, double* outBuffer, double* inBuffer, int xnumberAdded,
                                           int ynumberAdded, int znumberAdded,
                                           int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	if (cartDim[0] > 1) {
		int bcount = 0;
		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outBuffer[bcount] = tempVector[1 + additionalBinNumber + i][j][k][l];
						++bcount;
					}
				}
			}
		}

		int leftCoord[3];
		for (int i = 0; i < 3; ++i) {
			leftCoord[i] = cartCoord[i];
		}
		leftCoord[0] -= 1;
		if (leftCoord[0] < 0) {
			leftCoord[0] = cartDim[0] - 1;
		}
		int leftRank;
		MPI_Cart_rank(cartComm, leftCoord, &leftRank);

		int rightCoord[3];
		for (int i = 0; i < 3; ++i) {
			rightCoord[i] = cartCoord[i];
		}
		rightCoord[0] += 1;
		if (rightCoord[0] >= cartDim[0]) {
			rightCoord[0] = 0;
		}
		int rightRank;
		MPI_Cart_rank(cartComm, rightCoord, &rightRank);

		MPI_Status status;
		int number = ynumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber);
		MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, leftRank, MPI_TEMPVECTOR_LEFT, inBuffer, number, MPI_DOUBLE, rightRank,
		             MPI_TEMPVECTOR_LEFT, cartComm, &status);


		bcount = 0;
		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[xnumberAdded - 1 - additionalBinNumber + i][j][k][l] = inBuffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[xnumberAdded - 1 - additionalBinNumber + i][j][k][l] = tempVector[1 + additionalBinNumber + i][j][k][l
						];
					}
				}
			}
		}
	}
}

void sendLargeVectorToRight(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[0] > 1) {
		int bcount = 0;
		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						buffer[bcount] = tempVector[xnumberAdded - 2 - additionalBinNumber - i][j][k][l];
						++bcount;
					}
				}
			}
		}
		int rightCoord[3];
		for (int i = 0; i < 3; ++i) {
			rightCoord[i] = cartCoord[i];
		}
		rightCoord[0] += 1;
		if (rightCoord[0] >= cartDim[0]) {
			rightCoord[0] = 0;
		}
		int rightRank;
		MPI_Cart_rank(cartComm, rightCoord, &rightRank);

		int leftCoord[3];
		for (int i = 0; i < 3; ++i) {
			leftCoord[i] = cartCoord[i];
		}
		leftCoord[0] -= 1;
		if (leftCoord[0] < 0) {
			leftCoord[0] = cartDim[0] - 1;
		}
		int leftRank;
		MPI_Cart_rank(cartComm, leftCoord, &leftRank);

		MPI_Send(buffer, ynumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, rightRank,
		         MPI_TEMPVECTOR_RIGHT, cartComm);
	}

}

void sendLargeVectorToLeft(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {

	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	if (cartDim[0] > 1) {
		int bcount = 0;
		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						buffer[bcount] = tempVector[1 + additionalBinNumber + i][j][k][l];
						++bcount;
					}
				}
			}
		}

		int leftCoord[3];
		for (int i = 0; i < 3; ++i) {
			leftCoord[i] = cartCoord[i];
		}
		leftCoord[0] -= 1;
		if (leftCoord[0] < 0) {
			leftCoord[0] = cartDim[0] - 1;
		}
		int leftRank;
		MPI_Cart_rank(cartComm, leftCoord, &leftRank);

		MPI_Send(buffer, ynumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, leftRank,
		         MPI_TEMPVECTOR_LEFT, cartComm);
	}

}

void receiveLargeVectorFromRight(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded,
                                 int znumberAdded,
                                 int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[0] > 1) {
		int rightCoord[3];
		for (int i = 0; i < 3; ++i) {
			rightCoord[i] = cartCoord[i];
		}
		rightCoord[0] += 1;
		if (rightCoord[0] >= cartDim[0]) {
			rightCoord[0] = 0;
		}
		int rightRank;
		MPI_Cart_rank(cartComm, rightCoord, &rightRank);
		MPI_Status status;
		MPI_Recv(buffer, ynumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, rightRank,
		         MPI_TEMPVECTOR_LEFT, cartComm,
		         &status);

		int bcount = 0;
		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[xnumberAdded - 1 - additionalBinNumber + i][j][k][l] = buffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		//printf("aaa\n");
		/*for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[xnumberAdded - 1 - additionalBinNumber + i][j][k][l] = tempVector[1 + additionalBinNumber + i][j][k][l];
					}
				}
			}
		}*/
	}
}

void receiveLargeVectorFromLeft(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded,
                                int znumberAdded,
                                int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[0] > 1) {
		int leftCoord[3];
		for (int i = 0; i < 3; ++i) {
			leftCoord[i] = cartCoord[i];
		}
		leftCoord[0] -= 1;
		if (leftCoord[0] < 0) {
			leftCoord[0] = cartDim[0] - 1;
		}
		int leftRank;
		MPI_Cart_rank(cartComm, leftCoord, &leftRank);
		MPI_Status status;
		MPI_Recv(buffer, ynumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, leftRank,
		         MPI_TEMPVECTOR_RIGHT, cartComm,
		         &status);

		int bcount = 0;
		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[additionalBinNumber - i][j][k][l] = buffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		//printf("aaa\n");
		/*for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[additionalBinNumber - i][j][k][l] = tempVector[xnumberAdded - 2 - additionalBinNumber - i][j][k][l];
					}
				}
			}
		}*/
	}
}

void sendLargeVectorToBackReceiveFromFront(double**** tempVector, double* outBuffer, double* inBuffer, int xnumberAdded,
                                           int ynumberAdded, int znumberAdded,
                                           int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[1] > 1) {
		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outBuffer[bcount] = tempVector[i][ynumberAdded - 2 - additionalBinNumber - j][k][l];
						++bcount;
					}
				}
			}
		}
		int backCoord[3];
		for (int i = 0; i < 3; ++i) {
			backCoord[i] = cartCoord[i];
		}
		backCoord[1] += 1;
		if (backCoord[1] >= cartDim[1]) {
			backCoord[1] = 0;
		}
		int backRank;
		MPI_Cart_rank(cartComm, backCoord, &backRank);

		int frontCoord[3];
		for (int i = 0; i < 3; ++i) {
			frontCoord[i] = cartCoord[i];
		}
		frontCoord[1] -= 1;
		if (frontCoord[1] < 0) {
			frontCoord[1] = cartDim[1] - 1;
		}
		int frontRank;
		MPI_Cart_rank(cartComm, frontCoord, &frontRank);

		MPI_Status status;
		int number = xnumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber);
		MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, backRank, MPI_TEMPVECTOR_RIGHT, inBuffer, number, MPI_DOUBLE, frontRank,
		             MPI_TEMPVECTOR_RIGHT, cartComm, &status);


		bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][additionalBinNumber - j][k][l] = inBuffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][additionalBinNumber - j][k][l] = tempVector[i][ynumberAdded - 2 - additionalBinNumber - j][k][l];
						//tempVector[i][ynumberAdded - 1 - additionalBinNumber + j][k][l] = tempVector[i][1 + additionalBinNumber + j][k][l];
					}
				}
			}
		}
	}
}

void sendLargeVectorToFrontReceiveFromBack(double**** tempVector, double* outBuffer, double* inBuffer, int xnumberAdded,
                                           int ynumberAdded, int znumberAdded,
                                           int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	if (cartDim[1] > 1) {
		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outBuffer[bcount] = tempVector[i][1 + additionalBinNumber + j][k][l];
						++bcount;
					}
				}
			}
		}

		int frontCoord[3];
		for (int i = 0; i < 3; ++i) {
			frontCoord[i] = cartCoord[i];
		}
		frontCoord[1] -= 1;
		if (frontCoord[1] < 0) {
			frontCoord[1] = cartDim[1] - 1;
		}
		int frontRank;

		MPI_Cart_rank(cartComm, frontCoord, &frontRank);

		int backCoord[3];
		for (int i = 0; i < 3; ++i) {
			backCoord[i] = cartCoord[i];
		}
		backCoord[1] += 1;
		if (backCoord[1] >= cartDim[1]) {
			backCoord[1] = 0;
		}
		int backRank;
		MPI_Cart_rank(cartComm, backCoord, &backRank);


		MPI_Status status;
		int number = xnumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber);
		MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, frontRank, MPI_TEMPVECTOR_LEFT, inBuffer, number, MPI_DOUBLE, backRank,
		             MPI_TEMPVECTOR_LEFT, cartComm, &status);


		bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][ynumberAdded - 1 - additionalBinNumber + j][k][l] = inBuffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][ynumberAdded - 1 - additionalBinNumber + j][k][l] = tempVector[i][1 + additionalBinNumber + j][k][l
						];
					}
				}
			}
		}
	}
}

void sendLargeVectorToBack(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[1] > 1) {
		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						buffer[bcount] = tempVector[i][ynumberAdded - 2 - additionalBinNumber - j][k][l];
						++bcount;
					}
				}
			}
		}
		int backCoord[3];
		for (int i = 0; i < 3; ++i) {
			backCoord[i] = cartCoord[i];
		}
		backCoord[1] += 1;
		if (backCoord[1] >= cartDim[1]) {
			backCoord[1] = 0;
		}
		int backRank;
		MPI_Cart_rank(cartComm, backCoord, &backRank);


		MPI_Send(buffer, xnumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, backRank,
		         MPI_TEMPVECTOR_RIGHT, cartComm);
	}

}

void sendLargeVectorToFront(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {

	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	if (cartDim[1] > 1) {
		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						buffer[bcount] = tempVector[i][1 + additionalBinNumber + j][k][l];
						++bcount;
					}
				}
			}
		}

		int frontCoord[3];
		for (int i = 0; i < 3; ++i) {
			frontCoord[i] = cartCoord[i];
		}
		frontCoord[1] -= 1;
		if (frontCoord[1] < 0) {
			frontCoord[1] = cartDim[1] - 1;
		}
		int frontRank;
		MPI_Cart_rank(cartComm, frontCoord, &frontRank);


		MPI_Send(buffer, xnumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, frontRank,
		         MPI_TEMPVECTOR_LEFT, cartComm);
	}

}

void receiveLargeVectorFromBack(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded,
                                int znumberAdded,
                                int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[1] > 1) {
		int backCoord[3];
		for (int i = 0; i < 3; ++i) {
			backCoord[i] = cartCoord[i];
		}
		backCoord[1] += 1;
		if (backCoord[1] >= cartDim[1]) {
			backCoord[1] = 0;
		}
		int backRank;
		MPI_Cart_rank(cartComm, backCoord, &backRank);
		MPI_Status status;
		MPI_Recv(buffer, xnumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, backRank,
		         MPI_TEMPVECTOR_LEFT, cartComm,
		         &status);

		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][ynumberAdded - 1 - additionalBinNumber + j][k][l] = buffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][ynumberAdded - 1 - additionalBinNumber + j][k][l] = tempVector[i][1 + additionalBinNumber + j][k][l
						];
					}
				}
			}
		}
	}
}

void receiveLargeVectorFromFront(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded,
                                 int znumberAdded,
                                 int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[1] > 1) {
		int frontCoord[3];
		for (int i = 0; i < 3; ++i) {
			frontCoord[i] = cartCoord[i];
		}
		frontCoord[1] -= 1;
		if (frontCoord[1] < 0) {
			frontCoord[1] = cartDim[1] - 1;
		}
		int frontRank;
		MPI_Cart_rank(cartComm, frontCoord, &frontRank);
		MPI_Status status;
		MPI_Recv(buffer, xnumberAdded * znumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, frontRank,
		         MPI_TEMPVECTOR_RIGHT, cartComm,
		         &status);

		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][additionalBinNumber - j][k][l] = buffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][additionalBinNumber - j][k][l] = tempVector[i][ynumberAdded - 2 - additionalBinNumber - j][k][l];
					}
				}
			}
		}
	}
}

void sendLargeVectorToTopReceiveFromBottom(double**** tempVector, double* outBuffer, double* inBuffer, int xnumberAdded,
                                           int ynumberAdded, int znumberAdded,
                                           int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[2] > 1) {
		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outBuffer[bcount] = tempVector[i][j][znumberAdded - 2 - additionalBinNumber - k][l];
						++bcount;
					}
				}
			}
		}
		int topCoord[3];
		for (int i = 0; i < 3; ++i) {
			topCoord[i] = cartCoord[i];
		}
		topCoord[2] += 1;
		if (topCoord[2] >= cartDim[2]) {
			topCoord[2] = 0;
		}
		int topRank;
		MPI_Cart_rank(cartComm, topCoord, &topRank);

		int bottomCoord[3];
		for (int i = 0; i < 3; ++i) {
			bottomCoord[i] = cartCoord[i];
		}
		bottomCoord[2] -= 1;
		if (bottomCoord[2] < 0) {
			bottomCoord[2] = cartDim[2] - 1;
		}
		int bottomRank;
		MPI_Cart_rank(cartComm, bottomCoord, &bottomRank);

		MPI_Status status;
		int number = xnumberAdded * ynumberAdded * lnumber * (1 + additionalBinNumber);
		MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, topRank, MPI_TEMPVECTOR_RIGHT, inBuffer, number, MPI_DOUBLE, bottomRank,
		             MPI_TEMPVECTOR_RIGHT, cartComm, &status);


		bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][j][additionalBinNumber - k][l] = inBuffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][j][additionalBinNumber - k][l] = tempVector[i][j][znumberAdded - 2 - additionalBinNumber - k][l];
					}
				}
			}
		}
	}
}

void sendLargeVectorToBottomReceiveFromTop(double**** tempVector, double* outBuffer, double* inBuffer, int xnumberAdded,
                                           int ynumberAdded, int znumberAdded,
                                           int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	if (cartDim[2] > 1) {
		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outBuffer[bcount] = tempVector[i][j][1 + additionalBinNumber + k][l];
						++bcount;
					}
				}
			}
		}

		int bottomCoord[3];
		for (int i = 0; i < 3; ++i) {
			bottomCoord[i] = cartCoord[i];
		}
		bottomCoord[2] -= 1;
		if (bottomCoord[2] < 0) {
			bottomCoord[2] = cartDim[2] - 1;
		}
		int bottomRank;
		MPI_Cart_rank(cartComm, bottomCoord, &bottomRank);

		int topCoord[3];
		for (int i = 0; i < 3; ++i) {
			topCoord[i] = cartCoord[i];
		}
		topCoord[2] += 1;
		if (topCoord[2] >= cartDim[2]) {
			topCoord[2] = 0;
		}
		int topRank;
		MPI_Cart_rank(cartComm, topCoord, &topRank);

		MPI_Status status;
		int number = xnumberAdded * ynumberAdded * lnumber * (1 + additionalBinNumber);
		MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, bottomRank, MPI_TEMPVECTOR_LEFT, inBuffer, number, MPI_DOUBLE, topRank,
		             MPI_TEMPVECTOR_LEFT, cartComm, &status);


		bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][j][znumberAdded - 1 - additionalBinNumber + k][l] = inBuffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][j][znumberAdded - 1 - additionalBinNumber + k][l] = tempVector[i][j][1 + additionalBinNumber + k][l
						];
					}
				}
			}
		}
	}
}

void sendLargeVectorToTop(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                          int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[2] > 1) {
		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						buffer[bcount] = tempVector[i][j][znumberAdded - 2 - additionalBinNumber - k][l];
						++bcount;
					}
				}
			}
		}
		int topCoord[3];
		for (int i = 0; i < 3; ++i) {
			topCoord[i] = cartCoord[i];
		}
		topCoord[2] += 1;
		if (topCoord[2] >= cartDim[2]) {
			topCoord[2] = 0;
		}
		int topRank;
		MPI_Cart_rank(cartComm, topCoord, &topRank);

		MPI_Send(buffer, xnumberAdded * ynumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, topRank,
		         MPI_TEMPVECTOR_RIGHT, cartComm);
	}

}

void sendLargeVectorToBottom(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded,
                             int znumberAdded,
                             int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {

	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	if (cartDim[2] > 1) {
		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						buffer[bcount] = tempVector[i][j][1 + additionalBinNumber + k][l];
						++bcount;
					}
				}
			}
		}

		int bottomCoord[3];
		for (int i = 0; i < 3; ++i) {
			bottomCoord[i] = cartCoord[i];
		}
		bottomCoord[2] -= 1;
		if (bottomCoord[2] < 0) {
			bottomCoord[2] = cartDim[2] - 1;
		}
		int bottomRank;
		MPI_Cart_rank(cartComm, bottomCoord, &bottomRank);

		MPI_Send(buffer, xnumberAdded * ynumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, bottomRank,
		         MPI_TEMPVECTOR_LEFT, cartComm);
	}

}

void receiveLargeVectorFromTop(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded,
                               int znumberAdded,
                               int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[2] > 1) {
		int topCoord[3];
		for (int i = 0; i < 3; ++i) {
			topCoord[i] = cartCoord[i];
		}
		topCoord[2] += 1;
		if (topCoord[2] >= cartDim[2]) {
			topCoord[2] = 0;
		}
		int topRank;
		MPI_Cart_rank(cartComm, topCoord, &topRank);
		MPI_Status status;
		MPI_Recv(buffer, xnumberAdded * ynumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, topRank,
		         MPI_TEMPVECTOR_LEFT, cartComm,
		         &status);

		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][j][znumberAdded - 1 - additionalBinNumber + k][l] = buffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][j][znumberAdded - 1 - additionalBinNumber + k][l] = tempVector[i][j][1 + additionalBinNumber + k][l
						];
					}
				}
			}
		}
	}
}

void receiveLargeVectorFromBottom(double**** tempVector, double* buffer, int xnumberAdded, int ynumberAdded,
                                  int znumberAdded,
                                  int lnumber, int additionalBinNumber, MPI_Comm& cartComm) {
	int rank;
	int size;
	MPI_Comm_size(cartComm, &size);
	MPI_Comm_rank(cartComm, &rank);
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[2] > 1) {
		int bottomCoord[3];
		for (int i = 0; i < 3; ++i) {
			bottomCoord[i] = cartCoord[i];
		}
		bottomCoord[2] -= 1;
		if (bottomCoord[2] < 0) {
			bottomCoord[2] = cartDim[2] - 1;
		}
		int bottomRank;
		MPI_Cart_rank(cartComm, bottomCoord, &bottomRank);
		MPI_Status status;
		MPI_Recv(buffer, xnumberAdded * ynumberAdded * lnumber * (1 + additionalBinNumber), MPI_DOUBLE, bottomRank,
		         MPI_TEMPVECTOR_RIGHT, cartComm,
		         &status);

		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][j][additionalBinNumber - k][l] = buffer[bcount];
						++bcount;
					}
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						tempVector[i][j][additionalBinNumber - k][l] = tempVector[i][j][znumberAdded - 2 - additionalBinNumber - k][l];
					}
				}
			}
		}
	}
}

void sendCellParametersToLeftReceiveFromRight(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank,
                                              int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < 2 * additionalNumber + 2; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				outBuffer[bcount] = array[i][j][k];
				bcount++;
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * ynumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, leftRank, MPI_CELL_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             rightRank, MPI_CELL_PARAMETERS_LEFT, cartComm, &status);

	bcount = 0;

	for (int i = 0; i < 2 * additionalNumber + 2; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendCellParametersToRightReceiveFromLeft(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank,
                                              int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k];
				bcount++;
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * ynumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, rightRank, MPI_CELL_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             leftRank, MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendCellParametersLeft(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	if (cartDim[0] > 1) {
		int bcount = 0;

		for (int i = 0; i < 2 * additionalNumber + 2; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					outBuffer[bcount] = array[i][j][k];
					bcount++;
				}
			}
		}

		MPI_Send(outBuffer, (2 + 2 * additionalNumber) * ynumberAdded * znumberAdded, MPI_DOUBLE, leftRank,
		         MPI_CELL_PARAMETERS_LEFT, cartComm);
	}

}

void receiveCellParametersRight(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	if (cartDim[0] > 1) {
		MPI_Status status;
		MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * ynumberAdded * znumberAdded, MPI_DOUBLE, rightRank,
		         MPI_CELL_PARAMETERS_LEFT, cartComm, &status);

		int bcount = 0;

		for (int i = 0; i < 2 * additionalNumber + 2; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					tempArray[i][j][k] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellParametersRight(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                             int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[0] > 1) {
		int bcount = 0;
		for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k];
					bcount++;
				}
			}
		}

		MPI_Send(outBuffer, (2 + 2 * additionalNumber) * ynumberAdded * znumberAdded, MPI_DOUBLE, rightRank,
		         MPI_CELL_PARAMETERS_RIGHT, cartComm);
	}
}

void receiveCellParametersLeft(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                               int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[0] > 1) {
		MPI_Status status;
		MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * ynumberAdded * znumberAdded, MPI_DOUBLE, leftRank,
		         MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);

		int bcount = 0;
		for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					tempArray[i][j][k] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellParametersToFrontReceiveFromBack(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank,
                                              int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	//if (cartDim[1] > 1) {
	int bcount = 0;

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 * additionalNumber + 2; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				outBuffer[bcount] = array[i][j][k];
				bcount++;
			}
		}
	}
	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * xnumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, frontRank, MPI_CELL_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             backRank, MPI_CELL_PARAMETERS_LEFT, cartComm, &status);

	bcount = 0;

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
	//}
}

void sendCellParametersToBackReceiveFromFront(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank,
                                              int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	//if(cartDim[1] > 1){
	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k];
				bcount++;
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * xnumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, backRank, MPI_CELL_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             frontRank, MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);

	bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
	//}
}

void sendCellParametersFront(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                             int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[1] > 1) {
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < 2 * additionalNumber + 2; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					outBuffer[bcount] = array[i][j][k];
					bcount++;
				}
			}
		}

		MPI_Send(outBuffer, (2 + 2 * additionalNumber) * xnumberAdded * znumberAdded, MPI_DOUBLE, frontRank,
		         MPI_CELL_PARAMETERS_LEFT, cartComm);
	}

}

void receiveCellParametersBack(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                               int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	if (cartDim[1] > 1) {
		MPI_Status status;
		MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * xnumberAdded * znumberAdded, MPI_DOUBLE, backRank,
		         MPI_CELL_PARAMETERS_LEFT, cartComm, &status);

		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					tempArray[i][j][k] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellParametersBack(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int additionalNumber, MPI_Comm& cartComm, int rank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[1] > 1) {
		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k];
					bcount++;
				}
			}
		}

		MPI_Send(outBuffer, (2 + 2 * additionalNumber) * xnumberAdded * znumberAdded, MPI_DOUBLE, backRank,
		         MPI_CELL_PARAMETERS_RIGHT, cartComm);
	}
}

void receiveCellParametersFront(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	if (cartDim[1] > 1) {
		MPI_Status status;
		MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * xnumberAdded * znumberAdded, MPI_DOUBLE, frontRank,
		         MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);

		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					tempArray[i][j][k] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellParametersToBottomReceiveFromTop(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank,
                                              int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 * additionalNumber + 2; ++k) {
				outBuffer[bcount] = array[i][j][k];
				bcount++;
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * ynumberAdded * xnumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, bottomRank, MPI_CELL_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             topRank, MPI_CELL_PARAMETERS_LEFT, cartComm, &status);

	bcount = 0;

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendCellParametersToTopReceiveFromBottom(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank,
                                              int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k];
				bcount++;
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * ynumberAdded * xnumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, topRank, MPI_CELL_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             bottomRank, MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendCellParametersBottom(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                              int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	//if (cartDim[2] > 1) {
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 * additionalNumber + 2; ++k) {
					outBuffer[bcount] = array[i][j][k];
					bcount++;
				}
			}
		}

		MPI_Send(outBuffer, (2 + 2 * additionalNumber) * ynumberAdded * xnumberAdded, MPI_DOUBLE, bottomRank,
		         MPI_CELL_PARAMETERS_LEFT, cartComm);
	//}
}

void receiveCellParametersTop(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                              int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	//if (cartDim[2] > 1) {
		MPI_Status status;
		MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * ynumberAdded * xnumberAdded, MPI_DOUBLE, topRank,
		         MPI_CELL_PARAMETERS_LEFT, cartComm, &status);

		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
					tempArray[i][j][k] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	//}
}

void sendCellParametersTop(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	//if (cartDim[2] > 1) {
		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
					outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k];
					bcount++;
				}
			}
		}

		MPI_Send(outBuffer, (2 + 2 * additionalNumber) * ynumberAdded * xnumberAdded, MPI_DOUBLE, topRank,
		         MPI_CELL_PARAMETERS_RIGHT, cartComm);
	//}
}

void receiveCellParametersBottom(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                 int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	//if (cartDim[2] > 1) {
		MPI_Status status;
		MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * ynumberAdded * xnumberAdded, MPI_DOUBLE, bottomRank,
		         MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);

		int bcount = 0;
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
					tempArray[i][j][k] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	//}
}

void sendCellVectorParametersToLeftReceiveFromRight(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 3 * ynumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, leftRank, MPI_CELL_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             rightRank, MPI_CELL_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellVectorParametersToRightReceiveFromLeft(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 3 * ynumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, rightRank, MPI_CELL_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             leftRank, MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellVectorParametersLeft(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                  int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 3 * ynumberAdded * znumberAdded, MPI_DOUBLE, leftRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm);

}

void receiveCellVectorParametersRight(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                      int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                      int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 3 * ynumberAdded * znumberAdded, MPI_DOUBLE, rightRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellVectorParametersRight(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                   int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                   int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 3 * ynumberAdded * znumberAdded, MPI_DOUBLE, rightRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm);
}

void receiveCellVectorParametersLeft(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                     int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                     int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 3 * ynumberAdded * znumberAdded, MPI_DOUBLE, leftRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellVectorParametersToFrontReceiveFromBack(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 3 * xnumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, frontRank, MPI_CELL_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             backRank, MPI_CELL_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellVectorParametersToBackReceiveFromFront(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 3 * xnumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, backRank, MPI_CELL_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             frontRank, MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellVectorParametersFront(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                   int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                   int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 3 * xnumberAdded * znumberAdded, MPI_DOUBLE, frontRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm);

}

void receiveCellVectorParametersBack(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                     int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                     int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 3 * xnumberAdded * znumberAdded, MPI_DOUBLE, backRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellVectorParametersBack(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                  int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 3 * xnumberAdded * znumberAdded, MPI_DOUBLE, backRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm);
}

void receiveCellVectorParametersFront(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                      int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                      int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 3 * xnumberAdded * znumberAdded, MPI_DOUBLE, frontRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellVectorParametersToBottomReceiveFromTop(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 3 * ynumberAdded * xnumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, bottomRank, MPI_CELL_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             topRank, MPI_CELL_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellVectorParametersToTopReceiveFromBottom(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 3 * ynumberAdded * xnumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, topRank, MPI_CELL_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             bottomRank, MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellVectorParametersBottom(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                    int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 3 * ynumberAdded * xnumberAdded, MPI_DOUBLE, bottomRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm);

}

void receiveCellVectorParametersTop(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 3 * ynumberAdded * xnumberAdded, MPI_DOUBLE, topRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellVectorParametersTop(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                 int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 3 * ynumberAdded * xnumberAdded, MPI_DOUBLE, topRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm);
}

void receiveCellVectorParametersBottom(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                       int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                       int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 3 * ynumberAdded * xnumberAdded, MPI_DOUBLE, bottomRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendCellMatrixParametersToLeftReceiveFromRight(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 9 * ynumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, leftRank, MPI_CELL_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             rightRank, MPI_CELL_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendCellMatrixParametersToRightReceiveFromLeft(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 9 * ynumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, rightRank, MPI_CELL_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             leftRank, MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendCellMatrixParametersLeft(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                  int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 9 * ynumberAdded * znumberAdded, MPI_DOUBLE, leftRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm);
}

void receiveCellMatrixParametersRight(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                      int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                      int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 9 * ynumberAdded * znumberAdded, MPI_DOUBLE, rightRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendCellMatrixParametersRight(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                   int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                   int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 9 * ynumberAdded * znumberAdded, MPI_DOUBLE, rightRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm);
}

void receiveCellMatrixParametersLeft(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                     int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                     int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 9 * ynumberAdded * znumberAdded, MPI_DOUBLE, leftRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm,
	         &status);

	int bcount = 0;

	for (int i = 0; i < 2 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendCellMatrixParametersToFrontReceiveFromBack(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 9 * xnumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, frontRank, MPI_CELL_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             backRank, MPI_CELL_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendCellMatrixParametersToBackReceiveFromFront(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 9 * xnumberAdded * znumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, backRank, MPI_CELL_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             frontRank, MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);

	bcount = 0;

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendCellMatrixParametersFront(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                   int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                   int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 9 * xnumberAdded * znumberAdded, MPI_DOUBLE, frontRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm);
}

void receiveCellMatrixParametersBack(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                     int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                     int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 9 * xnumberAdded * znumberAdded, MPI_DOUBLE, backRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm,
	         &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendCellMatrixParametersBack(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                  int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 9 * xnumberAdded * znumberAdded, MPI_DOUBLE, backRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm);
}

void receiveCellMatrixParametersFront(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                      int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                      int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 9 * xnumberAdded * znumberAdded, MPI_DOUBLE, frontRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm,
	         &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < 2 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendCellMatrixParametersToBottomReceiveFromTop(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 9 * xnumberAdded * ynumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, bottomRank, MPI_CELL_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             topRank, MPI_CELL_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendCellMatrixParametersToTopReceiveFromBottom(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (2 + 2 * additionalNumber) * 9 * xnumberAdded * ynumberAdded;
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, topRank, MPI_CELL_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             bottomRank, MPI_CELL_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendCellMatrixParametersBottom(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                    int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 9 * xnumberAdded * ynumberAdded, MPI_DOUBLE, bottomRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm);
}

void receiveCellMatrixParametersTop(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 9 * xnumberAdded * ynumberAdded, MPI_DOUBLE, topRank,
	         MPI_CELL_PARAMETERS_LEFT, cartComm,
	         &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendCellMatrixParametersTop(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                 int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Send(outBuffer, (2 + 2 * additionalNumber) * 9 * xnumberAdded * ynumberAdded, MPI_DOUBLE, topRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm);
}

void receiveCellMatrixParametersBottom(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                       int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                       int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (2 + 2 * additionalNumber) * 9 * xnumberAdded * ynumberAdded, MPI_DOUBLE, bottomRank,
	         MPI_CELL_PARAMETERS_RIGHT, cartComm,
	         &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < 2 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeParametersToLeftReceiveFromRight(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank,
                                              int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				outBuffer[bcount] = array[i][j][k];
				bcount++;
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * (ynumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, leftRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             rightRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeParametersToRightReceiveFromLeft(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank,
                                              int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k];
				bcount++;
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * (ynumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, rightRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             leftRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeParametersLeft(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				outBuffer[bcount] = array[i][j][k];
				bcount++;
			}
		}
	}


	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, leftRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm);
}

void receiveNodeParametersRight(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, rightRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeParametersRight(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                             int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k];
				bcount++;
			}
		}
	}

	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, rightRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm);
}

void receiveNodeParametersLeft(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                               int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, leftRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeParametersToFrontReceiveFromBack(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank,
                                              int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				outBuffer[bcount] = array[i][j][k];
				bcount++;
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, frontRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             backRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeParametersToBackReceiveFromFront(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank,
                                              int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k];
				bcount++;
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, backRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             frontRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeParametersFront(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                             int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				outBuffer[bcount] = array[i][j][k];
				bcount++;
			}
		}
	}


	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, frontRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm);
}

void receiveNodeParametersBack(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                               int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, backRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeParametersBack(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int additionalNumber, MPI_Comm& cartComm, int rank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k];
				bcount++;
			}
		}
	}


	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, backRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm);
}

void receiveNodeParametersFront(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, frontRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeParametersToBottomReceiveFromTop(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank,
                                              int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				outBuffer[bcount] = array[i][j][k];
				bcount++;
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (ynumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, bottomRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             topRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);

	bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeParametersToTopReceiveFromBottom(double*** array, double* outBuffer, double*** tempArray, double* inBuffer,
                                              int xnumberAdded, int ynumberAdded, int znumberAdded,
                                              int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank,
                                              int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k];
				bcount++;
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (ynumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, topRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             bottomRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeParametersBottom(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                              int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				outBuffer[bcount] = array[i][j][k];
				bcount++;
			}
		}
	}


	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, bottomRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm);
}

void receiveNodeParametersTop(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                              int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, topRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeParametersTop(double*** array, double* outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k];
				bcount++;
			}
		}
	}


	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, topRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm);
}

void receiveNodeParametersBottom(double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                 int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, bottomRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				tempArray[i][j][k] = inBuffer[bcount];
				bcount++;
			}
		}
	}
}

void sendNodeVectorParametersToLeftReceiveFromRight(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, leftRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             rightRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendNodeVectorParametersToRightReceiveFromLeft(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, rightRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             leftRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					bcount++;
				}
			}
		}
	}
}

void sendNodeVectorParametersLeft(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                  int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}


	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, leftRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm);
}

void receiveNodeVectorParametersRight(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                      int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                      int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	//printf("receive node vector right from %d to %d\n", rightRank, rank);
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, rightRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					//alertNaNOrInfinity(tempArray[xnumber - 2 - i][j][k][l], "tempArray[xnumber - 2 - i][j][k]l] = NaN in receiveNodeVectorParametersRiht\n");
					bcount++;
				}
			}
		}
	}
}

void sendNodeVectorParametersRight(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                   int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                   int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k][l];
					bcount++;
				}
			}
		}
	}


	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, rightRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm);
}

void receiveNodeVectorParametersLeft(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                     int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                     int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	//printf("receive node vector left from %d to %d\n", leftRank, rank);
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, leftRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					//alertNaNOrInfinity(tempArray[2 + i][j][k][l], "tempArray[2 + i][j][k]l] = NaN in receiveNodeVectorParametersLeft\n");
					bcount++;
				}
			}
		}
	}
}

void sendNodeVectorParametersToFrontReceiveFromBack(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}

	//printf("send node vector left from %d to %d\n", rank, leftRank);

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, frontRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             backRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);

	bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					//alertNaNOrInfinity(tempArray[xnumber - 2 - i][j][k][l], "tempArray[xnumber - 2 - i][j][k]l] = NaN in receiveNodeVectorParametersRiht\n");
					bcount++;
				}
			}
		}
	}
}

void sendNodeVectorParametersToBackReceiveFromFront(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k][l];
					bcount++;
				}
			}
		}
	}

	//printf("send node vector right from %d to %d\n", rank, rightRank);
	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, backRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             frontRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					//alertNaNOrInfinity(tempArray[xnumber - 2 - i][j][k][l], "tempArray[xnumber - 2 - i][j][k]l] = NaN in receiveNodeVectorParametersRiht\n");
					bcount++;
				}
			}
		}
	}
}

void sendNodeVectorParametersFront(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                   int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                   int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}

	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, frontRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm);
}

void receiveNodeVectorParametersBack(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                     int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                     int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	//printf("receive node vector right from %d to %d\n", rightRank, rank);
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, backRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					//alertNaNOrInfinity(tempArray[xnumber - 2 - i][j][k][l], "tempArray[xnumber - 2 - i][j][k]l] = NaN in receiveNodeVectorParametersRiht\n");
					bcount++;
				}
			}
		}
	}
}

void sendNodeVectorParametersBack(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                  int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k][l];
					bcount++;
				}
			}
		}
	}


	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, backRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm);
}

void receiveNodeVectorParametersFront(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                      int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                      int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	//printf("receive node vector left from %d to %d\n", frontRank, rank);
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, frontRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					//alertNaNOrInfinity(tempArray[2 + i][j][k][l], "tempArray[2 + i][j][k]l] = NaN in receiveNodeVectorParametersLeft\n");
					bcount++;
				}
			}
		}
	}
}

void sendNodeVectorParametersToBottomReceiveFromTop(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}

	//printf("send node vector left from %d to %d\n", rank, leftRank);
	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (ynumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, bottomRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             topRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					//alertNaNOrInfinity(tempArray[xnumber - 2 - i][j][k][l], "tempArray[xnumber - 2 - i][j][k]l] = NaN in receiveNodeVectorParametersRiht\n");
					bcount++;
				}
			}
		}
	}
}

void sendNodeVectorParametersToTopReceiveFromBottom(Vector3d*** array, double* outBuffer, Vector3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k][l];
					bcount++;
				}
			}
		}
	}

	//printf("send node vector right from %d to %d\n", rank, rightRank);
	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (ynumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, topRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             bottomRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);

	bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					//alertNaNOrInfinity(tempArray[2 + i][j][k][l], "tempArray[2 + i][j][k]l] = NaN in receiveNodeVectorParametersLeft\n");
					bcount++;
				}
			}
		}
	}
}

void sendNodeVectorParametersBottom(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                    int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][k][l];
					bcount++;
				}
			}
		}
	}


	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, bottomRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm);
}

void receiveNodeVectorParametersTop(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	//printf("receive node vector right from %d to %d\n", rightRank, rank);
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, topRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					//alertNaNOrInfinity(tempArray[xnumber - 2 - i][j][k][l], "tempArray[xnumber - 2 - i][j][k]l] = NaN in receiveNodeVectorParametersRiht\n");
					bcount++;
				}
			}
		}
	}
}

void sendNodeVectorParametersTop(Vector3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                 int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k][l];
					bcount++;
				}
			}
		}
	}


	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, topRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm);
}

void receiveNodeVectorParametersBottom(Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                       int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                       int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	//printf("receive node vector left from %d to %d\n", frontRank, rank);
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 3 * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, bottomRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempArray[i][j][k][l] = inBuffer[bcount];
					//alertNaNOrInfinity(tempArray[2 + i][j][k][l], "tempArray[2 + i][j][k]l] = NaN in receiveNodeVectorParametersLeft\n");
					bcount++;
				}
			}
		}
	}
}

void sendNodeMatrixParametersToLeftReceiveFromRight(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, leftRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             rightRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMatrixParametersToRightReceiveFromLeft(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, rightRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             leftRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMatrixParametersLeft(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                  int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, leftRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm);
}

void receiveNodeMatrixParametersRight(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                      int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                      int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, rightRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMatrixParametersRight(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                   int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                   int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, rightRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm);
}

void receiveNodeMatrixParametersLeft(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                     int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                     int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, leftRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMatrixParametersToFrontReceiveFromBack(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}


	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, frontRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             backRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMatrixParametersToBackReceiveFromFront(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, backRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             frontRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMatrixParametersFront(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                   int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                   int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, frontRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm);
}

void receiveNodeMatrixParametersBack(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                     int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                     int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, backRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMatrixParametersBack(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                  int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, backRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm);
}

void receiveNodeMatrixParametersFront(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                      int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                      int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1), MPI_DOUBLE, frontRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMatrixParametersToBottomReceiveFromTop(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, bottomRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             topRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMatrixParametersToTopReceiveFromBottom(Matrix3d*** array, double* outBuffer, Matrix3d*** tempArray,
                                                    double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                    int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, topRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             bottomRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMatrixParametersBottom(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                    int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}


	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, bottomRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm);
}

void receiveNodeMatrixParametersTop(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                    int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, topRank,
	         MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMatrixParametersTop(Matrix3d*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                 int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k].matrix[l][m];
						bcount++;
					}
				}
			}
		}
	}

	MPI_Send(outBuffer, (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, topRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm);
}

void receiveNodeMatrixParametersBottom(Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                       int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                       int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer, (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1), MPI_DOUBLE, bottomRank,
	         MPI_NODE_PARAMETERS_RIGHT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < 3; ++m) {
						tempArray[i][j][k].matrix[l][m] = inBuffer[bcount];
						bcount++;
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersToLeftReceiveFromRight(MassMatrix*** array, double* outBuffer, MassMatrix*** tempArray,
                                                        double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                        int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                        int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		splineOrder + 3) * (2 * splineOrder + 3);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, leftRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             rightRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersToRightReceiveFromLeft(MassMatrix*** array, double* outBuffer, MassMatrix*** tempArray,
                                                        double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                        int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                        int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k].matrix[tempI][tempJ][tempK].matrix
										[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		splineOrder + 3) * (2 * splineOrder + 3);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, rightRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             leftRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersLeft(MassMatrix*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                      int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                      int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}

	MPI_Send(outBuffer,
	         (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, leftRank, MPI_NODE_PARAMETERS_LEFT, cartComm);
}

void receiveNodeMassMatrixParametersRight(MassMatrix*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                          int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                          int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	MPI_Recv(inBuffer,
	         (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, rightRank, MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersRight(MassMatrix*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                       int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                       int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[xnumberAdded - 2 - 2 * additionalNumber + i][j][k].matrix[tempI][tempJ][tempK].matrix
										[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}

	MPI_Send(outBuffer,
	         (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, rightRank, MPI_NODE_PARAMETERS_RIGHT, cartComm);
}

void receiveNodeMassMatrixParametersLeft(MassMatrix*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                         int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                         int leftRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer,
	         (3 + 2 * additionalNumber) * 9 * (ynumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, leftRank, MPI_NODE_PARAMETERS_RIGHT, cartComm,
	         &status);

	int bcount = 0;

	for (int i = 0; i < 3 + 2 * additionalNumber; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersToFrontReceiveFromBack(MassMatrix*** array, double* outBuffer, MassMatrix*** tempArray,
                                                        double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                        int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                        int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}


	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		splineOrder + 3) * (2 * splineOrder + 3);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, frontRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             backRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersToBackReceiveFromFront(MassMatrix*** array, double* outBuffer, MassMatrix*** tempArray,
                                                        double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                        int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                        int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k].matrix[tempI][tempJ][tempK].matrix
										[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		splineOrder + 3) * (2 * splineOrder + 3);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, backRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             frontRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersFront(MassMatrix*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                       int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                       int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}

	MPI_Send(outBuffer,
	         (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, frontRank, MPI_NODE_PARAMETERS_LEFT, cartComm);
}

void receiveNodeMassMatrixParametersBack(MassMatrix*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                         int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                         int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	MPI_Recv(inBuffer,
	         (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, backRank, MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersBack(MassMatrix*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                      int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                      int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[i][ynumberAdded - 2 - 2 * additionalNumber + j][k].matrix[tempI][tempJ][tempK].matrix
										[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}

	MPI_Send(outBuffer,
	         (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, backRank, MPI_NODE_PARAMETERS_RIGHT, cartComm);
}

void receiveNodeMassMatrixParametersFront(MassMatrix*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                          int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                          int frontRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer,
	         (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (znumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, frontRank, MPI_NODE_PARAMETERS_RIGHT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < 3 + 2 * additionalNumber; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersToBottomReceiveFromTop(MassMatrix*** array, double* outBuffer, MassMatrix*** tempArray,
                                                        double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                        int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                        int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		splineOrder + 3) * (2 * splineOrder + 3);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, bottomRank, MPI_NODE_PARAMETERS_LEFT, inBuffer, number, MPI_DOUBLE,
	             topRank, MPI_NODE_PARAMETERS_LEFT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersToTopReceiveFromBottom(MassMatrix*** array, double* outBuffer, MassMatrix*** tempArray,
                                                        double* inBuffer, int xnumberAdded, int ynumberAdded,
                                                        int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                                        int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k].matrix[tempI][tempJ][tempK].matrix
										[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}

	MPI_Status status;
	int number = (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		splineOrder + 3) * (2 * splineOrder + 3);
	MPI_Sendrecv(outBuffer, number, MPI_DOUBLE, topRank, MPI_NODE_PARAMETERS_RIGHT, inBuffer, number, MPI_DOUBLE,
	             bottomRank, MPI_NODE_PARAMETERS_RIGHT, cartComm, &status);
	bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersBottom(MassMatrix*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                        int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                        int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}


	MPI_Send(outBuffer,
	         (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, bottomRank, MPI_NODE_PARAMETERS_LEFT, cartComm);
}

void receiveNodeMassMatrixParametersTop(MassMatrix*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded,
                                        int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                        int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);


	MPI_Status status;
	MPI_Recv(inBuffer,
	         (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, topRank, MPI_NODE_PARAMETERS_LEFT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void sendNodeMassMatrixParametersTop(MassMatrix*** array, double* outBuffer, int xnumberAdded, int ynumberAdded,
                                     int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank,
                                     int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									outBuffer[bcount] = array[i][j][znumberAdded - 2 - 2 * additionalNumber + k].matrix[tempI][tempJ][tempK].matrix
										[l][m];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}

	MPI_Send(outBuffer,
	         (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, topRank, MPI_NODE_PARAMETERS_RIGHT, cartComm);
}

void receiveNodeMassMatrixParametersBottom(MassMatrix*** tempArray, double* inBuffer, int xnumberAdded,
                                           int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm,
                                           int rank, int bottomRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);

	MPI_Status status;
	MPI_Recv(inBuffer,
	         (3 + 2 * additionalNumber) * 9 * (xnumberAdded + 1) * (ynumberAdded + 1) * (2 * splineOrder + 3) * (2 *
		         splineOrder + 3) * (2 * splineOrder + 3), MPI_DOUBLE, bottomRank, MPI_NODE_PARAMETERS_RIGHT,
	         cartComm, &status);

	int bcount = 0;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < 3 + 2 * additionalNumber; ++k) {
				for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
					for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
						for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
							for (int l = 0; l < 3; ++l) {
								for (int m = 0; m < 3; ++m) {
									tempArray[i][j][k].matrix[tempI][tempJ][tempK].matrix[l][m] = inBuffer[bcount];
									bcount++;
								}
							}
						}
					}
				}
			}
		}
	}
}

void collectOutDoubleParameters(std::vector < Particle* >& outParticles, double* outDoubleParticlesParameters) {
	int doubleBcount = 0;
	for (int i = 0; i < outParticles.size(); ++i) {
		Particle* particle = outParticles[i];

		outDoubleParticlesParameters[doubleBcount] = particle->mass;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->charge;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->weight;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->coordinates.x;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->coordinates.y;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->coordinates.z;
		doubleBcount++;
		Vector3d momentum = particle->getMomentum();
		outDoubleParticlesParameters[doubleBcount] = momentum.x;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = momentum.y;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = momentum.z;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->initialMomentum.x;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->initialMomentum.y;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->initialMomentum.z;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->prevMomentum.x;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->prevMomentum.y;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->prevMomentum.z;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->dx;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->dy;
		doubleBcount++;
		outDoubleParticlesParameters[doubleBcount] = particle->dz;
		doubleBcount++;
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				outDoubleParticlesParameters[doubleBcount] = particle->rotationTensor.matrix[j][k];
				doubleBcount++;
			}
		}
	}
}

void collectOutIntegerParameters(std::vector < Particle* >& outParticles, ParticleTypeContainer* types, int typesNumber,
                                 int* outIntegerParticlesParameters) {
	int bcount = 0;
	for (int i = 0; i < outParticles.size(); ++i) {
		Particle* particle = outParticles[i];
		outIntegerParticlesParameters[bcount] = particle->number;
		bcount++;
		outIntegerParticlesParameters[bcount] = particle->chargeCount;
		bcount++;
		int type = -1;
		for (int t = 0; t < typesNumber; ++t) {
			if (particle->type == types[t].type) {
				type = t;
				break;
			}
		}
		if (type == -1) {
			std::string outputDir = outputDirectory;
			FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "particle has no type in send particles\n");
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}
		outIntegerParticlesParameters[bcount] = type;
		bcount++;
		outIntegerParticlesParameters[bcount] = particle->crossBoundaryCount;
		bcount++;
	}
}

void addParticlesFromParameters(std::vector < Particle* >& inParticles, std::vector<Particle*>& reservedParticles, ParticleTypeContainer* types, int typesNumber,
                                int inParticlesNumber[1], int* inIntegerParticlesParameters, double* inDoubleParticlesParameters) {
	int bcount = 0;
	int doubleBcount = 0;
	for (int i = 0; i < inParticlesNumber[0]; ++i) {

		int number = inIntegerParticlesParameters[bcount];
		bcount++;
		int chargeCount = inIntegerParticlesParameters[bcount];
		bcount++;
		int type = inIntegerParticlesParameters[bcount];
		if (type < 0 || type >= typesNumber) {
			std::string outputDir = outputDirectory;
			FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			fprintf(errorLogFile, "particle has wrong type in receive particles\n");
			fclose(errorLogFile);
			MPI_Finalize();
			exit(0);
		}
		ParticleTypes particleType = types[type].type;
		bcount++;
		int crossBoundary = inIntegerParticlesParameters[bcount];
		bcount++;

		double mass = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double charge = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double weight = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double x = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double y = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double z = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double momentumX = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double momentumY = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double momentumZ = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double initialMomentumX = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double initialMomentumY = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double initialMomentumZ = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double prevMomentumX = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double prevMomentumY = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double prevMomentumZ = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double dx = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double dy = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		double dz = inDoubleParticlesParameters[doubleBcount];
		doubleBcount++;
		Particle* particle;

		if(reservedParticles.size() > 0) {
			particle = reservedParticles.back();
			reservedParticles.pop_back();
			particle->number = number;

			particle->mass = mass;
			particle->chargeCount = chargeCount;
			particle->charge = charge;
			particle->weight = weight;
			particle->type = particleType;

			particle->coordinates.x = x;
			particle->coordinates.y = y;
			particle->coordinates.z = z;

			particle->initialMomentum.x = initialMomentumX;
			particle->initialMomentum.y = initialMomentumY;
			particle->initialMomentum.z = initialMomentumZ;

			particle->dx = dx;
			particle->dy = dy;
			particle->dz = dz;

			particle->escaped = false;
			particle->crossBoundaryCount = 0;
		} else {
			particle = new Particle(number, mass, chargeCount, charge, weight, particleType, x, y, z,
		                                  initialMomentumX, initialMomentumY, initialMomentumZ, dx, dy, dz);
		}


		particle->setMomentum(momentumX, momentumY, momentumZ);
		particle->prevMomentum.x = prevMomentumX;
		particle->prevMomentum.x = prevMomentumY;
		particle->prevMomentum.x = prevMomentumZ;
		particle->crossBoundaryCount = crossBoundary;

		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				particle->rotationTensor.matrix[j][k] = inDoubleParticlesParameters[doubleBcount];
				doubleBcount++;
			}
		}


		inParticles.push_back(particle);
	}
}

void sendLeftReceiveRightParticlesGeneral(std::vector < Particle* >& outParticles, std::vector < Particle* >& inParticles,
                                          std::vector<Particle*>& reservedParticles, ParticleTypeContainer* types, int typesNumber, int verbosity, MPI_Comm& cartComm, int rank,
                                          int leftRank, int rightRank, const int numberOfIntegerParameters,
                                          const int numberOfDoubleParameters, int outParticlesNumber[1], int inParticlesNumber[1], const bool sendingLeft, const bool receivingRight) {
	int* outIntegerParticlesParameters;
	int* inIntegerParticlesParameters;
	double* outDoubleParticlesParameters;
	double* inDoubleParticlesParameters;

	if (sendingLeft) {
		outIntegerParticlesParameters = new int[numberOfIntegerParameters * outParticlesNumber[0]];
		collectOutIntegerParameters(outParticles, types, typesNumber, outIntegerParticlesParameters);
	}

	if (receivingRight) {
		inIntegerParticlesParameters = new int[numberOfIntegerParameters * inParticlesNumber[0]];
	}

	if (sendingLeft && receivingRight) {
		MPI_Status status;
		MPI_Sendrecv(outIntegerParticlesParameters, numberOfIntegerParameters * outParticlesNumber[0], MPI_INT, leftRank,
		             MPI_SEND_INTEGER_NUMBER_LEFT, inIntegerParticlesParameters,
		             numberOfIntegerParameters * inParticlesNumber[0], MPI_INT, rightRank, MPI_SEND_INTEGER_NUMBER_LEFT,
		             cartComm, &status);
	} else if (sendingLeft) {
		MPI_Send(outIntegerParticlesParameters, numberOfIntegerParameters * outParticlesNumber[0], MPI_INT, leftRank,
		         MPI_SEND_INTEGER_NUMBER_LEFT, cartComm);
	} else if (receivingRight) {
		MPI_Status status;
		MPI_Recv(inIntegerParticlesParameters, numberOfIntegerParameters * inParticlesNumber[0], MPI_INT, rightRank,
		         MPI_SEND_INTEGER_NUMBER_LEFT, cartComm, &status);
	}

	MPI_Barrier(cartComm);
	if (sendingLeft) {
		outDoubleParticlesParameters = new double[numberOfDoubleParameters * outParticlesNumber[0]];
		collectOutDoubleParameters(outParticles, outDoubleParticlesParameters);
	}

	if (receivingRight) {
		inDoubleParticlesParameters = new double[numberOfDoubleParameters * inParticlesNumber[0]];
		if (verbosity > 2) printf("receive inDoubleParameters right rank = %d\n", rank);
	}

	if (sendingLeft && receivingRight) {
		MPI_Status status;
		MPI_Sendrecv(outDoubleParticlesParameters, numberOfDoubleParameters * outParticlesNumber[0], MPI_DOUBLE, leftRank,
		             MPI_SEND_DOUBLE_NUMBER_LEFT, inDoubleParticlesParameters,
		             numberOfDoubleParameters * inParticlesNumber[0], MPI_DOUBLE, rightRank, MPI_SEND_DOUBLE_NUMBER_LEFT,
		             cartComm, &status);
	} else if (sendingLeft) {
		MPI_Send(outDoubleParticlesParameters, numberOfDoubleParameters * outParticlesNumber[0], MPI_DOUBLE, leftRank,
		         MPI_SEND_DOUBLE_NUMBER_LEFT, cartComm);
	} else if (receivingRight) {
		MPI_Status status;
		MPI_Recv(inDoubleParticlesParameters, numberOfDoubleParameters * inParticlesNumber[0], MPI_DOUBLE, rightRank,
		         MPI_SEND_DOUBLE_NUMBER_LEFT, cartComm, &status);
	}

	if (receivingRight) {
		addParticlesFromParameters(inParticles, reservedParticles, types, typesNumber, inParticlesNumber, inIntegerParticlesParameters, inDoubleParticlesParameters);
	}

	if (sendingLeft) {
		delete[] outIntegerParticlesParameters;
		delete[] outDoubleParticlesParameters;
	}
	if (receivingRight) {
		delete[] inIntegerParticlesParameters;
		delete[] inDoubleParticlesParameters;
	}
	if (verbosity > 2) printf("finish send left receive right rank = %d\n", rank);
}

void sendRightReceiveLeftParticlesGeneral(std::vector < Particle* >& outParticles, std::vector < Particle* >& inParticles,
                                          std::vector<Particle*>& reservedParticles, ParticleTypeContainer* types, int typesNumber, int verbosity, MPI_Comm& cartComm, int rank,
                                          int leftRank, int rightRank, const int numberOfIntegerParameters,
                                          const int numberOfDoubleParameters, int outParticlesNumber[1], int inParticlesNumber[1], const bool sendingRight, const bool receivingLeft) {
	int* outIntegerParticlesParameters;
	int* inIntegerParticlesParameters;
	double* outDoubleParticlesParameters;
	double* inDoubleParticlesParameters;

	if (sendingRight) {
		outIntegerParticlesParameters = new int[numberOfIntegerParameters * outParticlesNumber[0]];
		collectOutIntegerParameters(outParticles, types, typesNumber, outIntegerParticlesParameters);
	}

	if (receivingLeft) {
		inIntegerParticlesParameters = new int[numberOfIntegerParameters * inParticlesNumber[0]];
	}

	if (sendingRight && receivingLeft) {
		MPI_Status status;
		MPI_Sendrecv(outIntegerParticlesParameters, numberOfIntegerParameters * outParticlesNumber[0], MPI_INT, rightRank,
		             MPI_SEND_INTEGER_NUMBER_RIGHT, inIntegerParticlesParameters,
		             numberOfIntegerParameters * inParticlesNumber[0], MPI_INT, leftRank, MPI_SEND_INTEGER_NUMBER_RIGHT,
		             cartComm, &status);
	} else if (sendingRight) {
		MPI_Send(outIntegerParticlesParameters, numberOfIntegerParameters * outParticlesNumber[0], MPI_INT, rightRank,
		         MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm);
	} else if (receivingLeft) {
		MPI_Status status;
		MPI_Recv(inIntegerParticlesParameters, numberOfIntegerParameters * inParticlesNumber[0], MPI_INT, leftRank,
		         MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
	}

	MPI_Barrier(cartComm);

	if (sendingRight) {
		outDoubleParticlesParameters = new double[numberOfDoubleParameters * outParticlesNumber[0]];
		collectOutDoubleParameters(outParticles, outDoubleParticlesParameters);
	}

	if (receivingLeft) {
		inDoubleParticlesParameters = new double[numberOfDoubleParameters * inParticlesNumber[0]];
	}

	if (sendingRight && receivingLeft) {
		MPI_Status status;
		MPI_Sendrecv(outDoubleParticlesParameters, numberOfDoubleParameters * outParticlesNumber[0], MPI_DOUBLE, rightRank,
		             MPI_SEND_DOUBLE_NUMBER_RIGHT, inDoubleParticlesParameters,
		             numberOfDoubleParameters * inParticlesNumber[0], MPI_DOUBLE, leftRank, MPI_SEND_DOUBLE_NUMBER_RIGHT,
		             cartComm, &status);
	} else if (sendingRight) {
		MPI_Send(outDoubleParticlesParameters, numberOfDoubleParameters * outParticlesNumber[0], MPI_DOUBLE, rightRank,
		         MPI_SEND_DOUBLE_NUMBER_RIGHT, cartComm);
	} else if (receivingLeft) {
		MPI_Status status;
		MPI_Recv(inDoubleParticlesParameters, numberOfDoubleParameters * inParticlesNumber[0], MPI_DOUBLE, leftRank,
		         MPI_SEND_DOUBLE_NUMBER_RIGHT, cartComm, &status);
	}

	if (receivingLeft) {
		addParticlesFromParameters(inParticles, reservedParticles, types, typesNumber, inParticlesNumber, inIntegerParticlesParameters, inDoubleParticlesParameters);
	}
	if (sendingRight) {
		delete[] outDoubleParticlesParameters;
		delete[] outIntegerParticlesParameters;
	}
	if (receivingLeft) {
		delete[] inIntegerParticlesParameters;
		delete[] inDoubleParticlesParameters;
	}
	if (verbosity > 2) printf("finish send righ receive left rank = %d\n", rank);
}

void sendLeftReceiveRightParticles(std::vector < Particle * >& outParticles, std::vector < Particle * >& inParticles,
                                   std::vector<Particle*>& reservedParticles, ParticleTypeContainer* types, int typesNumber, bool periodic,
                                   int verbosity, MPI_Comm& cartComm, int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	const int numberOfIntegerParameters = 4;
	const int numberOfDoubleParameters = 27;

	int outParticlesNumber[1];
	int inParticlesNumber[1];

	outParticlesNumber[0] = outParticles.size();
	inParticlesNumber[0] = 0;

	if (verbosity > 2) printf("send outParticlesNumber left rank = %d number = %d\n", rank, outParticlesNumber[0]);
	if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || periodic) {
		MPI_Status status;
		MPI_Sendrecv(outParticlesNumber, 1, MPI_INT, leftRank, MPI_SEND_INTEGER_NUMBER_LEFT, inParticlesNumber, 1, MPI_INT,
		             rightRank, MPI_SEND_INTEGER_NUMBER_LEFT, cartComm, &status);
	} else if (cartCoord[0] == 0) {
		MPI_Status status;
		MPI_Recv(inParticlesNumber, 1, MPI_INT, rightRank, MPI_SEND_INTEGER_NUMBER_LEFT, cartComm, &status);
	} else if (cartCoord[0] == cartDim[0] - 1) {
		MPI_Send(outParticlesNumber, 1, MPI_INT, leftRank, MPI_SEND_INTEGER_NUMBER_LEFT, cartComm);
	}

	MPI_Barrier(cartComm);

	const bool sendingLeft = (outParticlesNumber[0] > 0) && (cartCoord[0] > 0 || periodic);
	const bool receivingRight = (inParticlesNumber[0] > 0) && (cartCoord[0] < cartDim[0] - 1 || periodic);

	sendLeftReceiveRightParticlesGeneral(outParticles, inParticles, reservedParticles, types, typesNumber, verbosity, cartComm, rank,
	                                     leftRank, rightRank, numberOfIntegerParameters, numberOfDoubleParameters,
	                                     outParticlesNumber, inParticlesNumber, sendingLeft, receivingRight);
}

void sendRightReceiveLeftParticles(std::vector < Particle * >& outParticles, std::vector < Particle * >& inParticles,
                                   std::vector<Particle*>& reservedParticles, ParticleTypeContainer* types, int typesNumber, bool periodic,
                                   int verbosity, MPI_Comm& cartComm, int rank, int leftRank, int rightRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	const int numberOfIntegerParameters = 4;
	const int numberOfDoubleParameters = 27;

	int outParticlesNumber[1];
	int inParticlesNumber[1];

	outParticlesNumber[0] = outParticles.size();
	inParticlesNumber[0] = 0;

	if (verbosity > 2) printf("send outParticlesNumber right rank = %d number = %d\n", rank, outParticlesNumber[0]);

	if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || periodic) {
		MPI_Status status;
		MPI_Sendrecv(outParticlesNumber, 1, MPI_INT, rightRank, MPI_SEND_INTEGER_NUMBER_RIGHT, inParticlesNumber, 1, MPI_INT,
		             leftRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
	} else if (cartCoord[0] == 0) {
		MPI_Send(outParticlesNumber, 1, MPI_INT, rightRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm);
	} else if (cartCoord[0] == cartDim[0] - 1) {
		MPI_Status status;
		MPI_Recv(inParticlesNumber, 1, MPI_INT, leftRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
	}

	if (verbosity > 2) printf("in particle number = %d\n", inParticlesNumber[0]);

	MPI_Barrier(cartComm);

	const bool sendingRight = (outParticlesNumber[0] > 0) && (cartCoord[0] < cartDim[0] - 1 || periodic);
	const bool receivingLeft = (inParticlesNumber[0] > 0) && (cartCoord[0] > 0 || periodic);

	sendRightReceiveLeftParticlesGeneral(outParticles, inParticles, reservedParticles, types, typesNumber, verbosity, cartComm, rank,
	                                     leftRank, rightRank, numberOfIntegerParameters, numberOfDoubleParameters,
	                                     outParticlesNumber, inParticlesNumber, sendingRight, receivingLeft);
}

void sendFrontReceiveBackParticles(std::vector < Particle * >& outParticles, std::vector < Particle * >& inParticles,
                                   std::vector<Particle*>& reservedParticles, ParticleTypeContainer* types, int typesNumber, bool periodic,
                                   int verbosity, MPI_Comm& cartComm, int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	const int numberOfIntegerParameters = 4;
	const int numberOfDoubleParameters = 27;

	int outParticlesNumber[1];
	int inParticlesNumber[1];

	outParticlesNumber[0] = outParticles.size();
	inParticlesNumber[0] = 0;

	if (verbosity > 2) printf("send outParticlesNumber left rank = %d number = %d\n", rank, outParticlesNumber[0]);
	if ((cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1) || periodic) {
		MPI_Status status;
		MPI_Sendrecv(outParticlesNumber, 1, MPI_INT, frontRank, MPI_SEND_INTEGER_NUMBER_LEFT, inParticlesNumber, 1, MPI_INT,
		             backRank, MPI_SEND_INTEGER_NUMBER_LEFT, cartComm, &status);
	} else if (cartCoord[1] == 0) {
		MPI_Status status;
		MPI_Recv(inParticlesNumber, 1, MPI_INT, backRank, MPI_SEND_INTEGER_NUMBER_LEFT, cartComm, &status);
	} else if (cartCoord[1] == cartDim[1] - 1) {
		MPI_Send(outParticlesNumber, 1, MPI_INT, frontRank, MPI_SEND_INTEGER_NUMBER_LEFT, cartComm);
	}

	MPI_Barrier(cartComm);

	const bool sendingLeft = (outParticlesNumber[0] > 0) && (cartCoord[1] > 0 || periodic);
	const bool receivingRight = (inParticlesNumber[0] > 0) && (cartCoord[1] < cartDim[1] - 1 || periodic);

	sendLeftReceiveRightParticlesGeneral(outParticles, inParticles, reservedParticles, types, typesNumber, verbosity, cartComm, rank,
	                                     frontRank, backRank, numberOfIntegerParameters, numberOfDoubleParameters,
	                                     outParticlesNumber, inParticlesNumber, sendingLeft, receivingRight);
}

void sendBackReceiveFrontParticles(std::vector < Particle * >& outParticles, std::vector < Particle * >& inParticles,
                                   std::vector<Particle*>& reservedParticles, ParticleTypeContainer* types, int typesNumber, bool periodic,
                                   int verbosity, MPI_Comm& cartComm, int rank, int frontRank, int backRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	const int numberOfIntegerParameters = 4;
	const int numberOfDoubleParameters = 27;

	int outParticlesNumber[1];
	int inParticlesNumber[1];

	outParticlesNumber[0] = outParticles.size();
	inParticlesNumber[0] = 0;

	if (verbosity > 2) printf("send outParticlesNumber right rank = %d number = %d\n", rank, outParticlesNumber[0]);

	if ((cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1) || periodic) {
		MPI_Status status;
		MPI_Sendrecv(outParticlesNumber, 1, MPI_INT, backRank, MPI_SEND_INTEGER_NUMBER_RIGHT, inParticlesNumber, 1, MPI_INT,
		             frontRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
	} else if (cartCoord[1] == 0) {
		MPI_Send(outParticlesNumber, 1, MPI_INT, backRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm);
	} else if (cartCoord[1] == cartDim[1] - 1) {
		MPI_Status status;
		MPI_Recv(inParticlesNumber, 1, MPI_INT, frontRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
	}

	if (verbosity > 2) printf("in particle number = %d\n", inParticlesNumber[0]);

	MPI_Barrier(cartComm);

	const bool sendingRight = (outParticlesNumber[0] > 0) && (cartCoord[1] < cartDim[1] - 1 || periodic);
	const bool receivingLeft = (inParticlesNumber[0] > 0) && (cartCoord[1] > 0 || periodic);

	sendRightReceiveLeftParticlesGeneral(outParticles, inParticles, reservedParticles, types, typesNumber, verbosity, cartComm, rank,
	                                     frontRank, backRank, numberOfIntegerParameters, numberOfDoubleParameters,
	                                     outParticlesNumber, inParticlesNumber, sendingRight, receivingLeft);
}

void sendBottomReceiveTopParticles(std::vector < Particle * >& outParticles, std::vector < Particle * >& inParticles,
                                   std::vector<Particle*>& reservedParticles, ParticleTypeContainer* types, int typesNumber, bool periodic,
                                   int verbosity, MPI_Comm& cartComm, int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	const int numberOfIntegerParameters = 4;
	const int numberOfDoubleParameters = 27;

	int outParticlesNumber[1];
	int inParticlesNumber[1];

	outParticlesNumber[0] = outParticles.size();
	inParticlesNumber[0] = 0;

	if (verbosity > 2) printf("send outParticlesNumber left rank = %d number = %d\n", rank, outParticlesNumber[0]);
	if ((cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1) || periodic) {
		MPI_Status status;
		MPI_Sendrecv(outParticlesNumber, 1, MPI_INT, bottomRank, MPI_SEND_INTEGER_NUMBER_LEFT, inParticlesNumber, 1, MPI_INT,
		             topRank, MPI_SEND_INTEGER_NUMBER_LEFT, cartComm, &status);
	} else if (cartCoord[2] == 0) {
		MPI_Status status;
		MPI_Recv(inParticlesNumber, 1, MPI_INT, topRank, MPI_SEND_INTEGER_NUMBER_LEFT, cartComm, &status);
	} else if (cartCoord[2] == cartDim[2] - 1) {
		MPI_Send(outParticlesNumber, 1, MPI_INT, bottomRank, MPI_SEND_INTEGER_NUMBER_LEFT, cartComm);
	}

	MPI_Barrier(cartComm);

	const bool sendingLeft = (outParticlesNumber[0] > 0) && (cartCoord[2] > 0 || periodic);
	const bool receivingRight = (inParticlesNumber[0] > 0) && (cartCoord[2] < cartDim[2] - 1 || periodic);

	sendLeftReceiveRightParticlesGeneral(outParticles, inParticles, reservedParticles, types, typesNumber, verbosity, cartComm, rank,
	                                     bottomRank, topRank, numberOfIntegerParameters, numberOfDoubleParameters,
	                                     outParticlesNumber, inParticlesNumber, sendingLeft, receivingRight);
}

void sendTopReceiveBottomParticles(std::vector < Particle * >& outParticles, std::vector < Particle * >& inParticles,
                                   std::vector<Particle*>& reservedParticles, ParticleTypeContainer* types, int typesNumber, bool periodic,
                                   int verbosity, MPI_Comm& cartComm, int rank, int bottomRank, int topRank) {
	int periods[MPI_dim];
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	const int numberOfIntegerParameters = 4;
	const int numberOfDoubleParameters = 27;

	int outParticlesNumber[1];
	int inParticlesNumber[1];

	outParticlesNumber[0] = outParticles.size();
	inParticlesNumber[0] = 0;

	if (verbosity > 2) printf("send outParticlesNumber right rank = %d number = %d\n", rank, outParticlesNumber[0]);

	if ((cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1) || periodic) {
		MPI_Status status;
		MPI_Sendrecv(outParticlesNumber, 1, MPI_INT, topRank, MPI_SEND_INTEGER_NUMBER_RIGHT, inParticlesNumber, 1, MPI_INT,
		             bottomRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
	} else if (cartCoord[2] == 0) {
		MPI_Send(outParticlesNumber, 1, MPI_INT, topRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm);
	} else if (cartCoord[2] == cartDim[2] - 1) {
		MPI_Status status;
		MPI_Recv(inParticlesNumber, 1, MPI_INT, bottomRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
	}

	if (verbosity > 2) printf("in particle number = %d\n", inParticlesNumber[0]);

	MPI_Barrier(cartComm);

	const bool sendingRight = (outParticlesNumber[0] > 0) && (cartCoord[2] < cartDim[2] - 1 || periodic);
	const bool receivingLeft = (inParticlesNumber[0] > 0) && (cartCoord[2] > 0 || periodic);

	sendRightReceiveLeftParticlesGeneral(outParticles, inParticles, reservedParticles, types, typesNumber, verbosity, cartComm, rank,
	                                     bottomRank, topRank, numberOfIntegerParameters, numberOfDoubleParameters,
	                                     outParticlesNumber, inParticlesNumber, sendingRight, receivingLeft);
}

void exchangeLargeVector(double**** vector, int xnumberAdded, int ynumberAdded, int znumberAdded, int lnumber,
                         int additionalBinNumber, bool periodicX, bool
                         periodicY, bool periodicZ, MPI_Comm& cartComm, int* cartCoord, int* cartDim,
                         double* leftOutGmresBuffer, double* rightOutGmresBuffer, double* leftInGmresBuffer,
                         double* rightInGmresBuffer, double* frontOutGmresBuffer, double* backOutGmresBuffer,
                         double* frontInGmresBuffer, double* backInGmresBuffer, double* bottomOutGmresBuffer,
                         double* topOutGmresBuffer, double* bottomInGmresBuffer, double* topInGmresBuffer) {
	if (periodicX || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
		sendLargeVectorToLeftReceiveFromRight(vector, leftOutGmresBuffer, rightInGmresBuffer, xnumberAdded, ynumberAdded,
		                                      znumberAdded, lnumber, additionalBinNumber, cartComm);
	} else if (cartCoord[0] == 0 && cartDim[0] > 1) {
		receiveLargeVectorFromRight(vector, rightInGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                            additionalBinNumber, cartComm);
	} else if (cartCoord[0] == cartDim[0] - 1 && cartDim[0] > 1) {
		sendLargeVectorToLeft(vector, leftOutGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                      additionalBinNumber, cartComm);
	}

	if (periodicX || (cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1)) {
		sendLargeVectorToRightReceiveFromLeft(vector, rightOutGmresBuffer, leftInGmresBuffer, xnumberAdded, ynumberAdded,
		                                      znumberAdded, lnumber, additionalBinNumber, cartComm);
	} else if (cartCoord[0] == 0 && cartDim[0] > 1) {
		sendLargeVectorToRight(vector, rightOutGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                       additionalBinNumber, cartComm);
	} else if (cartCoord[0] == cartDim[0] - 1 && cartDim[0] > 1) {
		receiveLargeVectorFromLeft(vector, leftInGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                           additionalBinNumber, cartComm);
	}

	if (periodicY || (cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1)) {
		sendLargeVectorToFrontReceiveFromBack(vector, frontOutGmresBuffer, backInGmresBuffer, xnumberAdded, ynumberAdded,
		                                      znumberAdded, lnumber, additionalBinNumber, cartComm);
	} else if (cartCoord[1] == 0 && cartDim[1] > 1) {
		receiveLargeVectorFromBack(vector, backInGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                           additionalBinNumber, cartComm);
	} else if (cartCoord[1] == cartDim[1] - 1 && cartDim[1] > 1) {
		sendLargeVectorToFront(vector, frontOutGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                       additionalBinNumber, cartComm);
	}

	if (periodicY || (cartCoord[1] > 0 && cartCoord[1] < cartDim[1] - 1)) {
		sendLargeVectorToBackReceiveFromFront(vector, backOutGmresBuffer, frontInGmresBuffer, xnumberAdded, ynumberAdded,
		                                      znumberAdded, lnumber, additionalBinNumber, cartComm);
	} else if (cartCoord[1] == 0 && cartDim[1] > 1) {
		sendLargeVectorToBack(vector, backOutGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                      additionalBinNumber, cartComm);
	} else if (cartCoord[1] == cartDim[1] - 1 && cartDim[1] > 1) {
		receiveLargeVectorFromFront(vector, frontInGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                            additionalBinNumber, cartComm);
	}

	if (periodicZ || (cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1)) {
		sendLargeVectorToBottomReceiveFromTop(vector, bottomOutGmresBuffer, topInGmresBuffer, xnumberAdded, ynumberAdded,
		                                      znumberAdded, lnumber, additionalBinNumber, cartComm);
	} else if (cartCoord[2] == 0 && cartDim[2] > 1) {
		receiveLargeVectorFromTop(vector, topInGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                          additionalBinNumber, cartComm);
	} else if (cartCoord[2] == cartDim[2] - 1 && cartDim[2] > 1) {
		sendLargeVectorToBottom(vector, bottomOutGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                        additionalBinNumber, cartComm);
	}

	if (periodicZ || (cartCoord[2] > 0 && cartCoord[2] < cartDim[2] - 1)) {
		sendLargeVectorToTopReceiveFromBottom(vector, topOutGmresBuffer, bottomInGmresBuffer, xnumberAdded, ynumberAdded,
		                                      znumberAdded, lnumber, additionalBinNumber, cartComm);
	} else if (cartCoord[2] == 0 && cartDim[2] > 1) {
		sendLargeVectorToTop(vector, topOutGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                     additionalBinNumber, cartComm);
	} else if (cartCoord[2] == cartDim[2] - 1 && cartDim[2] > 1) {
		receiveLargeVectorFromBottom(vector, bottomInGmresBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
		                             additionalBinNumber, cartComm);
	}
}
