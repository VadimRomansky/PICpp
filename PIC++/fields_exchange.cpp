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
void Simulation::exchangeEfield() {
	exchangeGeneralEfield(Efield);
	exchangeGeneralEfield(tempEfield);
	exchangeGeneralEfield(newEfield);
}

void Simulation::exchangeGeneralEfield(Vector3d*** field) {
	exchangeGeneralEfieldX(field);
	exchangeGeneralEfieldY(field);
	exchangeGeneralEfieldZ(field);
}

void Simulation::exchangeGeneralEfieldX(Vector3d*** field) {
	if (cartDim[0] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int l = 0; l < 3; ++l) {
						rightOutVectorNodeBuffer[bcount] = field[xnumberAdded - additionalBinNumber - 1 - i][j][k][l];
						bcount++;
					}
				}
			}
		}
		int number = (2 + additionalBinNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1);
		//int numberLeft = (1 + additionalBinNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1);
		if ((cartCoord[0] < cartDim[0] - 1 && cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, rightRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
			MPI_Status status;
			MPI_Sendrecv(rightOutVectorNodeBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_RIGHT, leftInVectorNodeBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_RIGHT, cartComm, &status);
			if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, rightRank);
		} else if (cartCoord[0] == 0) {
			MPI_Send(rightOutVectorNodeBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_RIGHT, cartComm);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Status status;
			MPI_Recv(leftInVectorNodeBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_RIGHT, cartComm, &status);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field left from %d to %d\n", leftRank, rank);
			for (int i = 0; i < 2 + additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						for (int l = 0; l < 3; ++l) {
							field[additionalBinNumber + 1 - i][j][k][l] = leftInVectorNodeBuffer[bcount];
							bcount++;
						}
					}
				}
			}
		}

		MPI_Barrier(cartComm);
		bcount = 0;
		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int l = 0; l < 3; ++l) {
						leftOutVectorNodeBuffer[bcount] = field[1 + additionalBinNumber + i][j][k][l];
						bcount++;
					}
				}
			}
		}

		if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			MPI_Status status;
			MPI_Sendrecv(leftOutVectorNodeBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_LEFT, rightInVectorNodeBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_LEFT, cartComm, &status);
			if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, leftRank);
		} else if (cartCoord[0] == 0) {
			MPI_Status status;
			MPI_Recv(rightInVectorNodeBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_LEFT, cartComm, &status);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Send(leftOutVectorNodeBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_LEFT, cartComm);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field rigth from %d to %d\n", rightRank, rank);
			for (int i = 0; i < 2 + additionalBinNumber; ++i)
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						for (int l = 0; l < 3; ++l) {
							field[xnumberAdded - additionalBinNumber - 1 + i][j][k][l] = rightInVectorNodeBuffer[bcount];
							bcount++;
						}
					}
				}
		}

		if ((cartCoord[0] == 0) && (boundaryConditionType != PERIODIC)) {
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						//field[0][j][k] = field[2][j][k];
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							field[i][j][k] = field[1 + additionalBinNumber][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							field[i][j][k] = field[additionalBinNumber][j][k];
						}
					}
				}
			}
		}

		if ((cartCoord[0] == cartDim[0] - 1) && (boundaryConditionType != PERIODIC)) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						field[xnumberAdded - 1 - i][j][k] = field[xnumberAdded - 1 - additionalBinNumber][j][k];
					}
				}
			}
		}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						field[xnumberAdded - additionalBinNumber + i][j][k] = field[2 + additionalBinNumber + i][j][k];
						field[i][j][k] = field[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k];
					}
					field[xnumberAdded - additionalBinNumber - 1][j][k] = field[1 + additionalBinNumber][j][k];
				}
			}
		} else {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					field[0][j][k] = field[2][j][k];
					for (int i = 0; i <= additionalBinNumber; ++i) {
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							field[i][j][k] = field[1 + additionalBinNumber][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							field[i][j][k] = field[additionalBinNumber][j][k];
						}
						//field[xnumberAdded - 1 - i][j][k] = field[xnumberAdded - 1 - additionalBinNumber][j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeGeneralEfieldY(Vector3d*** field) {
	if (cartDim[1] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int l = 0; l < 3; ++l) {
						backOutVectorNodeBuffer[bcount] = field[i][ynumberAdded - additionalBinNumber - 1 - j][k][l];
						bcount++;
					}
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, backRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (2 + additionalBinNumber) * 3 * (xnumberAdded + 1) * (znumberAdded + 1);
		MPI_Sendrecv(backOutVectorNodeBuffer, number, MPI_DOUBLE, backRank, MPI_EFIELD_RIGHT, frontInVectorNodeBuffer, number, MPI_DOUBLE, frontRank, MPI_EFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, backRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int l = 0; l < 3; ++l) {
						field[i][additionalBinNumber + 1 - j][k][l] = frontInVectorNodeBuffer[bcount];
						bcount++;
					}
				}
			}
		}


		MPI_Barrier(cartComm);

		bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int l = 0; l < 3; ++l) {
						//backOutVectorNodeBuffer[bcount] = field[i][ynumberAdded - additionalBinNumber - 2 - j][k][l];
						frontOutVectorNodeBuffer[bcount] = field[i][1 + additionalBinNumber + j][k][l];
						bcount++;
					}
				}
			}
		}

		MPI_Sendrecv(frontOutVectorNodeBuffer, number, MPI_DOUBLE, frontRank, MPI_EFIELD_LEFT, backInVectorNodeBuffer, number, MPI_DOUBLE, backRank, MPI_EFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, frontRank);

		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", backRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i)
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int l = 0; l < 3; ++l) {
						field[i][ynumberAdded - additionalBinNumber - 1 + j][k][l] = backInVectorNodeBuffer[bcount];
						bcount++;
					}
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded + 1; ++j) {
						field[i][j][k] = field[i][1 + additionalBinNumber][k];
					}
				} else {
					for (int j = 0; j <= additionalBinNumber; ++j) {
						field[i][ynumberAdded - additionalBinNumber + j][k] = field[i][2 + additionalBinNumber + j][k];
						field[i][j][k] = field[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k];
					}
					field[i][ynumberAdded - additionalBinNumber - 1][k] = field[i][1 + additionalBinNumber][k];
				}
			}
		}
	}
}

void Simulation::exchangeGeneralEfieldZ(Vector3d*** field) {
	if (cartDim[2] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						topOutVectorNodeBuffer[bcount] = field[i][j][znumberAdded - additionalBinNumber - 1 - k][l];
						bcount++;
					}
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, topRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (2 + additionalBinNumber) * 3 * (xnumberAdded + 1) * (ynumberAdded + 1);
		MPI_Sendrecv(topOutVectorNodeBuffer, number, MPI_DOUBLE, topRank, MPI_EFIELD_RIGHT, bottomInVectorNodeBuffer, number, MPI_DOUBLE, bottomRank, MPI_EFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, topRank);

		//MPI_Barrier(cartComm);
		bcount = 0;
		if (verbosity > 2) printf("receive general field left from %d to %d\n", bottomRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						field[i][j][additionalBinNumber + 1 - k][l] = bottomInVectorNodeBuffer[bcount];
						bcount++;
					}
				}
			}
		}


		MPI_Barrier(cartComm);

		bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						bottomOutVectorNodeBuffer[bcount] = field[i][j][1 + additionalBinNumber + k][l];
						bcount++;
					}
				}
			}
		}


		MPI_Sendrecv(bottomOutVectorNodeBuffer, number, MPI_DOUBLE, bottomRank, MPI_EFIELD_LEFT, topInVectorNodeBuffer, number, MPI_DOUBLE, topRank, MPI_EFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, bottomRank);

		//MPI_Barrier(cartComm);
		bcount = 0;
		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", topRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i)
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						field[i][j][znumberAdded - additionalBinNumber - 1 + k][l] = topInVectorNodeBuffer[bcount];
						bcount++;
					}
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						field[i][j][k] = field[i][j][1 + additionalBinNumber];
					}
				} else {
					for (int k = 0; k <= additionalBinNumber; ++k) {
						field[i][j][znumberAdded - additionalBinNumber + k] = field[i][j][2 + additionalBinNumber + k];
						field[i][j][k] = field[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k];
					}
					field[i][j][znumberAdded - additionalBinNumber - 1] = field[i][j][1 + additionalBinNumber];
				}
			}
		}
	}
}

void Simulation::exchangeGeneralBfield(Vector3d*** field) {
	exchangeGeneralBfieldX(field);
	exchangeGeneralBfieldY(field);
	exchangeGeneralBfieldZ(field);
}

void Simulation::exchangeGeneralBfieldX(Vector3d*** field) {
	if (cartDim[0] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < 3; ++l) {
						rightOutVectorCellBuffer[bcount] = field[xnumberAdded - additionalBinNumber - 2 - i][j][k][l];
						leftOutVectorCellBuffer[bcount] = field[1 + additionalBinNumber + i][j][k][l];
						bcount++;
					}
				}
			}
		}
		int number = (1 + additionalBinNumber) * 3 * (ynumberAdded) * (znumberAdded);
		if ((cartCoord[0] < cartDim[0] - 1 && cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, rightRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
			MPI_Status status;
			MPI_Sendrecv(rightOutVectorCellBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_RIGHT, leftInVectorCellBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_RIGHT, cartComm, &status);
			if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, rightRank);
		} else if (cartCoord[0] == 0) {
			MPI_Send(rightOutVectorCellBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_RIGHT, cartComm);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Status status;
			MPI_Recv(leftInVectorCellBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_RIGHT, cartComm, &status);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field left from %d to %d\n", leftRank, rank);
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						for (int l = 0; l < 3; ++l) {
							field[additionalBinNumber - i][j][k][l] = leftInVectorCellBuffer[bcount];
							bcount++;
						}
					}
				}
			}
		}

		MPI_Barrier(cartComm);

		if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			MPI_Status status;
			MPI_Sendrecv(leftOutVectorCellBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_LEFT, rightInVectorCellBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_LEFT, cartComm, &status);
			if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, leftRank);
		} else if (cartCoord[0] == 0) {
			MPI_Status status;
			MPI_Recv(rightInVectorCellBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_LEFT, cartComm, &status);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Send(leftOutVectorCellBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_LEFT, cartComm);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field rigth from %d to %d\n", rightRank, rank);
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						for (int l = 0; l < 3; ++l) {
							field[xnumberAdded - 1 - additionalBinNumber + i][j][k][l] = rightInVectorCellBuffer[bcount];
							bcount++;
						}
					}
				}
			}
		}

		if ((cartCoord[0] == 0) && (boundaryConditionType != PERIODIC)) {
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						//field[0][j][k] = field[2][j][k];
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							field[i][j][k] = field[additionalBinNumber + 1][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							field[i][j][k] = field[additionalBinNumber + 1][j][k];
						}
					}
				}
			}
		}

		if ((cartCoord[0] == cartDim[0] - 1) && (boundaryConditionType != PERIODIC)) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						field[xnumberAdded - 2 - i][j][k] = field[xnumberAdded - 2 - additionalBinNumber][j][k];
					}
				}
			}
		}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						field[xnumberAdded - 1 - additionalBinNumber + i][j][k] = field[1 + additionalBinNumber + i][j][k];
						field[i][j][k] = field[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k];
					}
				}
			}
		} else {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					//field[0][j][k] = field[2][j][k];
					for (int i = 0; i <= additionalBinNumber; ++i) {
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							field[i][j][k] = Vector3d(0, 0, 0);
						} else if (boundaryConditionType == FREE_BOTH) {
							field[i][j][k] = field[additionalBinNumber + 1][j][k];
						}
						//field[xnumberAdded - 1 - i][j][k] = field[xnumberAdded - 2 - additionalBinNumber][j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeGeneralBfieldY(Vector3d*** field) {
	if (cartDim[1] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < 3; ++l) {
						backOutVectorCellBuffer[bcount] = field[i][ynumberAdded - additionalBinNumber - 2 - j][k][l];
						frontOutVectorCellBuffer[bcount] = field[i][1 + additionalBinNumber + j][k][l];
						bcount++;
					}
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, backRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (1 + additionalBinNumber) * 3 * (xnumberAdded) * (znumberAdded);
		MPI_Sendrecv(backOutVectorCellBuffer, number, MPI_DOUBLE, backRank, MPI_BFIELD_RIGHT, frontInVectorCellBuffer, number, MPI_DOUBLE, frontRank, MPI_BFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, backRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field left from %d to %d\n", frontRank, rank);
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < 3; ++l) {
						field[i][additionalBinNumber - j][k][l] = frontInVectorCellBuffer[bcount];
						bcount++;
					}
				}
			}
		}


		MPI_Barrier(cartComm);


		MPI_Sendrecv(frontOutVectorCellBuffer, number, MPI_DOUBLE, frontRank, MPI_BFIELD_LEFT, backInVectorCellBuffer, number, MPI_DOUBLE, backRank, MPI_BFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, frontRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", backRank, rank);
		for (int i = 0; i < xnumberAdded; ++i)
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < 3; ++l) {
						field[i][ynumberAdded - additionalBinNumber - 1 + j][k][l] = backInVectorCellBuffer[bcount];
						bcount++;
					}
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded; ++j) {
						field[i][j][k] = field[i][1 + additionalBinNumber][k];
					}
				} else {
					for (int j = 0; j <= additionalBinNumber; ++j) {
						field[i][ynumberAdded - 1 - additionalBinNumber + j][k] = field[i][1 + additionalBinNumber + j][k];
						field[i][j][k] = field[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeGeneralBfieldZ(Vector3d*** field) {
	if (cartDim[2] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						topOutVectorCellBuffer[bcount] = field[i][j][znumberAdded - additionalBinNumber - 2 - k][l];
						bottomOutVectorCellBuffer[bcount] = field[i][j][1 + additionalBinNumber + k][l];
						bcount++;
					}
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, topRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (1 + additionalBinNumber) * 3 * (xnumberAdded) * (ynumberAdded);
		MPI_Sendrecv(topOutVectorCellBuffer, number, MPI_DOUBLE, topRank, MPI_BFIELD_RIGHT, bottomInVectorCellBuffer, number, MPI_DOUBLE, bottomRank, MPI_BFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, topRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field left from %d to %d\n", bottomRank, rank);
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						field[i][j][additionalBinNumber - k][l] = bottomInVectorCellBuffer[bcount];
						bcount++;
					}
				}
			}
		}


		MPI_Barrier(cartComm);


		MPI_Sendrecv(bottomOutVectorCellBuffer, number, MPI_DOUBLE, bottomRank, MPI_BFIELD_LEFT, topInVectorCellBuffer, number, MPI_DOUBLE, topRank, MPI_BFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, bottomRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", topRank, rank);
		for (int i = 0; i < xnumberAdded; ++i)
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						field[i][j][znumberAdded - 1 - additionalBinNumber + k][l] = topInVectorCellBuffer[bcount];
						bcount++;
					}
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded; ++k) {
						field[i][j][k] = field[i][j][1 + additionalBinNumber];
					}
				} else {
					for (int k = 0; k <= additionalBinNumber; ++k) {
						field[i][j][znumberAdded - 1 - additionalBinNumber + k] = field[i][j][1 + additionalBinNumber + k];
						field[i][j][k] = field[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeGeneralScalarCellField(double**** field) {
	exchangeGeneralScalarCellFieldX(field);
	exchangeGeneralScalarCellFieldY(field);
	exchangeGeneralScalarCellFieldZ(field);
}

void Simulation::exchangeGeneralScalarCellFieldX(double**** field) {
	if (cartDim[0] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						rightOutCellBuffer[bcount] = field[xnumberAdded - additionalBinNumber - 2 - i][j][k][0];
						leftOutCellBuffer[bcount] = field[1 + additionalBinNumber + i][j][k][0];
						bcount++;
				}
			}
		}
		int number = (1 + additionalBinNumber)* (ynumberAdded) * (znumberAdded);
		if ((cartCoord[0] < cartDim[0] - 1 && cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, rightRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
			MPI_Status status;
			MPI_Sendrecv(rightOutCellBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_RIGHT, leftInCellBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_RIGHT, cartComm, &status);
			if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, rightRank);
		} else if (cartCoord[0] == 0) {
			MPI_Send(rightOutCellBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_RIGHT, cartComm);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Status status;
			MPI_Recv(leftInCellBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_RIGHT, cartComm, &status);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field left from %d to %d\n", leftRank, rank);
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						field[additionalBinNumber - i][j][k][0] = leftInCellBuffer[bcount];
						bcount++;
					}
				}
			}
		}

		MPI_Barrier(cartComm);

		if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			MPI_Status status;
			MPI_Sendrecv(leftOutCellBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_LEFT, rightInCellBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_LEFT, cartComm, &status);
			if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, leftRank);
		} else if (cartCoord[0] == 0) {
			MPI_Status status;
			MPI_Recv(rightInCellBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_LEFT, cartComm, &status);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Send(leftOutCellBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_LEFT, cartComm);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field rigth from %d to %d\n", rightRank, rank);
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
							field[xnumberAdded - 1 - additionalBinNumber + i][j][k][0] = rightInVectorCellBuffer[bcount];
							bcount++;
					}
				}
			}
		}

		if ((cartCoord[0] == 0) && (boundaryConditionType != PERIODIC)) {
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						//field[0][j][k] = field[2][j][k];
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							field[i][j][k][0] = field[additionalBinNumber + 1][j][k][0];
						} else if (boundaryConditionType == FREE_BOTH) {
							field[i][j][k][0] = field[additionalBinNumber + 1][j][k][0];
						}
					}
				}
			}
		}

		if ((cartCoord[0] == cartDim[0] - 1) && (boundaryConditionType != PERIODIC)) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						field[xnumberAdded - 2 - i][j][k][0] = field[xnumberAdded - 2 - additionalBinNumber][j][k][0];
					}
				}
			}
		}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						field[xnumberAdded - 1 - additionalBinNumber + i][j][k][0] = field[1 + additionalBinNumber + i][j][k][0];
						field[i][j][k][0] = field[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k][0];
					}
				}
			}
		} else {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					//field[0][j][k] = field[2][j][k];
					for (int i = 0; i <= additionalBinNumber; ++i) {
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							field[i][j][k][0] = 0;
						} else if (boundaryConditionType == FREE_BOTH) {
							field[i][j][k][0] = field[additionalBinNumber + 1][j][k][0];
						}
						//field[xnumberAdded - 1 - i][j][k] = field[xnumberAdded - 2 - additionalBinNumber][j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeGeneralScalarCellFieldY(double**** field) {
	if (cartDim[1] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					backOutCellBuffer[bcount] = field[i][ynumberAdded - additionalBinNumber - 2 - j][k][0];
					frontOutCellBuffer[bcount] = field[i][1 + additionalBinNumber + j][k][0];
					bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, backRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (1 + additionalBinNumber) * (xnumberAdded) * (znumberAdded);
		MPI_Sendrecv(backOutCellBuffer, number, MPI_DOUBLE, backRank, MPI_BFIELD_RIGHT, frontInCellBuffer, number, MPI_DOUBLE, frontRank, MPI_BFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, backRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field left from %d to %d\n", frontRank, rank);
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					field[i][additionalBinNumber - j][k][0] = frontInCellBuffer[bcount];
					bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);


		MPI_Sendrecv(frontOutCellBuffer, number, MPI_DOUBLE, frontRank, MPI_BFIELD_LEFT, backInCellBuffer, number, MPI_DOUBLE, backRank, MPI_BFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, frontRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", backRank, rank);
		for (int i = 0; i < xnumberAdded; ++i)
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					field[i][ynumberAdded - additionalBinNumber - 1 + j][k][0] = backInVectorCellBuffer[bcount];
					bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded; ++j) {
						field[i][j][k][0] = field[i][1 + additionalBinNumber][k][0];
					}
				} else {
					for (int j = 0; j <= additionalBinNumber + 1; ++j) {
						field[i][ynumberAdded - 1 - additionalBinNumber + j][k][0] = field[i][1 + additionalBinNumber + j][k][0];
						field[i][j][k][0] = field[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k][0];
					}
				}
			}
		}
	}
}

void Simulation::exchangeGeneralScalarCellFieldZ(double**** field) {
	if (cartDim[2] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					topOutCellBuffer[bcount] = field[i][j][znumberAdded - additionalBinNumber - 2 - k][0];
					bottomOutCellBuffer[bcount] = field[i][j][1 + additionalBinNumber + k][0];
					bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, topRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (1 + additionalBinNumber) * (xnumberAdded) * (ynumberAdded);
		MPI_Sendrecv(topOutCellBuffer, number, MPI_DOUBLE, topRank, MPI_BFIELD_RIGHT, bottomInCellBuffer, number, MPI_DOUBLE, bottomRank, MPI_BFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, topRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field left from %d to %d\n", bottomRank, rank);
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					field[i][j][additionalBinNumber - k][0] = bottomInCellBuffer[bcount];
					bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);


		MPI_Sendrecv(bottomOutCellBuffer, number, MPI_DOUBLE, bottomRank, MPI_BFIELD_LEFT, topInCellBuffer, number, MPI_DOUBLE, topRank, MPI_BFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, bottomRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", topRank, rank);
		for (int i = 0; i < xnumberAdded; ++i)
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
					field[i][j][znumberAdded - 1 - additionalBinNumber + k][0] = topInCellBuffer[bcount];
					bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded; ++k) {
						field[i][j][k][0] = field[i][j][1 + additionalBinNumber][0];
					}
				} else {
					for (int k = 0; k <= additionalBinNumber + 1; ++k) {
						field[i][j][znumberAdded - 1 - additionalBinNumber + k][0] = field[i][j][1 + additionalBinNumber + k][0];
						field[i][j][k][0] = field[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k][0];
					}
				}
			}
		}
	}
}

void Simulation::exchangeGeneralScalarNodeField(double**** field){
	exchangeGeneralScalarNodeFieldX(field);
	exchangeGeneralScalarNodeFieldY(field);
	exchangeGeneralScalarNodeFieldZ(field);
}

void Simulation::exchangeGeneralScalarNodeFieldX(double**** field){
	if (cartDim[0] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						rightOutNodeBuffer[bcount] = field[xnumberAdded - additionalBinNumber - 1 - i][j][k][0];
						bcount++;
				}
			}
		}
		int number = (2 + additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1);
		//int numberLeft = (1 + additionalBinNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1);
		if ((cartCoord[0] < cartDim[0] - 1 && cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, rightRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
			MPI_Status status;
			MPI_Sendrecv(rightOutNodeBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_RIGHT, leftInNodeBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_RIGHT, cartComm, &status);
			if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, rightRank);
		} else if (cartCoord[0] == 0) {
			MPI_Send(rightOutNodeBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_RIGHT, cartComm);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Status status;
			MPI_Recv(leftInNodeBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_RIGHT, cartComm, &status);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field left from %d to %d\n", leftRank, rank);
			for (int i = 0; i < 2 + additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
							field[additionalBinNumber + 1 - i][j][k][0] = leftInNodeBuffer[bcount];
							bcount++;
					}
				}
			}
		}

		MPI_Barrier(cartComm);
		bcount = 0;
		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						leftOutNodeBuffer[bcount] = field[1 + additionalBinNumber + i][j][k][0];
						bcount++;
				}
			}
		}

		if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			MPI_Status status;
			MPI_Sendrecv(leftOutNodeBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_LEFT, rightInNodeBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_LEFT, cartComm, &status);
			if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, leftRank);
		} else if (cartCoord[0] == 0) {
			MPI_Status status;
			MPI_Recv(rightInNodeBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_LEFT, cartComm, &status);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Send(leftOutNodeBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_LEFT, cartComm);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field rigth from %d to %d\n", rightRank, rank);
			for (int i = 0; i < 2 + additionalBinNumber; ++i)
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
							field[xnumberAdded - additionalBinNumber - 1 + i][j][k][0] = rightInNodeBuffer[bcount];
							bcount++;
					}
				}
		}

		if ((cartCoord[0] == 0) && (boundaryConditionType != PERIODIC)) {
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						//field[0][j][k] = field[2][j][k];
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							field[i][j][k][0] = field[1 + additionalBinNumber][j][k][0];
						} else if (boundaryConditionType == FREE_BOTH) {
							field[i][j][k][0] = field[additionalBinNumber][j][k][0];
						}
					}
				}
			}
		}

		if ((cartCoord[0] == cartDim[0] - 1) && (boundaryConditionType != PERIODIC)) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						field[xnumberAdded - 1 - i][j][k][0] = field[xnumberAdded - 1 - additionalBinNumber][j][k][0];
					}
				}
			}
		}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						field[xnumberAdded - additionalBinNumber + i][j][k][0] = field[2 + additionalBinNumber + i][j][k][0];
						field[i][j][k][0] = field[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k][0];
					}
					field[xnumberAdded - additionalBinNumber - 1][j][k][0] = field[1 + additionalBinNumber][j][k][0];
				}
			}
		} else {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					//field[0][j][k] = field[2][j][k];
					for (int i = 0; i <= additionalBinNumber; ++i) {
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							field[i][j][k][0] = field[1 + additionalBinNumber][j][k][0];
						} else if (boundaryConditionType == FREE_BOTH) {
							field[i][j][k][0] = field[additionalBinNumber][j][k][0];
						}
						//field[xnumberAdded - 1 - i][j][k] = field[xnumberAdded - 1 - additionalBinNumber][j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeGeneralScalarNodeFieldY(double**** field){
	if (cartDim[1] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						backOutNodeBuffer[bcount] = field[i][ynumberAdded - additionalBinNumber - 1 - j][k][0];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, backRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (2 + additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1);
		MPI_Sendrecv(backOutNodeBuffer, number, MPI_DOUBLE, backRank, MPI_EFIELD_RIGHT, frontInNodeBuffer, number, MPI_DOUBLE, frontRank, MPI_EFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, backRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						field[i][additionalBinNumber + 1 - j][k][0] = frontInNodeBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);

		bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						frontOutNodeBuffer[bcount] = field[i][1 + additionalBinNumber + j][k][0];
						bcount++;
				}
			}
		}

		MPI_Sendrecv(frontOutNodeBuffer, number, MPI_DOUBLE, frontRank, MPI_EFIELD_LEFT, backInNodeBuffer, number, MPI_DOUBLE, backRank, MPI_EFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, frontRank);

		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", backRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i)
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						field[i][ynumberAdded - additionalBinNumber - 1 + j][k][0] = backInNodeBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded + 1; ++j) {
						field[i][j][k][0] = field[i][1 + additionalBinNumber][k][0];
					}
				} else {
					for (int j = 0; j <= additionalBinNumber; ++j) {
						field[i][ynumberAdded - additionalBinNumber + j][k][0] = field[i][2 + additionalBinNumber + j][k][0];
						field[i][j][k][0] = field[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k][0];
					}
					field[i][ynumberAdded - additionalBinNumber - 1][k][0] = field[i][1 + additionalBinNumber][k][0];
				}
			}
		}
	}
}

void Simulation::exchangeGeneralScalarNodeFieldZ(double**** field){
	if (cartDim[2] > 1) {
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						topOutNodeBuffer[bcount] = field[i][j][znumberAdded - additionalBinNumber - 1 - k][0];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, topRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (2 + additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded + 1);
		MPI_Sendrecv(topOutNodeBuffer, number, MPI_DOUBLE, topRank, MPI_EFIELD_RIGHT, bottomInNodeBuffer, number, MPI_DOUBLE, bottomRank, MPI_EFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, topRank);

		//MPI_Barrier(cartComm);
		bcount = 0;
		if (verbosity > 2) printf("receive general field left from %d to %d\n", bottomRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						field[i][j][additionalBinNumber + 1 - k][0] = bottomInNodeBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);

		bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						bottomOutNodeBuffer[bcount] = field[i][j][1 + additionalBinNumber + k][0];
						bcount++;
				}
			}
		}


		MPI_Sendrecv(bottomOutNodeBuffer, number, MPI_DOUBLE, bottomRank, MPI_EFIELD_LEFT, topInNodeBuffer, number, MPI_DOUBLE, topRank, MPI_EFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, bottomRank);

		//MPI_Barrier(cartComm);
		bcount = 0;
		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", topRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i)
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						field[i][j][znumberAdded - additionalBinNumber - 1 + k][0] = topInNodeBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						field[i][j][k][0] = field[i][j][1 + additionalBinNumber][0];
					}
				} else {
					for (int k = 0; k <= additionalBinNumber; ++k) {
						field[i][j][znumberAdded - additionalBinNumber + k][0] = field[i][j][2 + additionalBinNumber + k][0];
						field[i][j][k][0] = field[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k][0];
					}
					field[i][j][znumberAdded - additionalBinNumber - 1][0] = field[i][j][1 + additionalBinNumber][0];
				}
			}
		}
	}
}

void Simulation::exchangeBunemanEfield(double*** fieldX, double*** fieldY, double*** fieldZ){
	exchangeBunemanExAlongX(fieldX);
	exchangeBunemanExAlongY(fieldX);
	exchangeBunemanExAlongZ(fieldX);

	exchangeBunemanEyAlongX(fieldY);
	exchangeBunemanEyAlongY(fieldY);
	exchangeBunemanEyAlongZ(fieldY);

	exchangeBunemanEzAlongX(fieldZ);
	exchangeBunemanEzAlongY(fieldZ);
	exchangeBunemanEzAlongZ(fieldZ);
}

void Simulation::exchangeBunemanExAlongX(double*** fieldX){
	if(cartDim[0] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						rightOutBunemanExBuffer[bcount] = fieldX[xnumberAdded - additionalBinNumber - 2 - i][j][k];
						leftOutBunemanExBuffer[bcount] = fieldX[1 + additionalBinNumber + i][j][k];
						bcount++;
				}
			}
		}
		int number = (1 + additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded + 1);
		if ((cartCoord[0] < cartDim[0] - 1 && cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, rightRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
			MPI_Status status;
			MPI_Sendrecv(rightOutBunemanExBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_RIGHT, leftInBunemanExBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_RIGHT, cartComm, &status);
			if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, rightRank);
		} else if (cartCoord[0] == 0) {
			MPI_Send(rightOutBunemanExBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_RIGHT, cartComm);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Status status;
			MPI_Recv(leftInBunemanExBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_RIGHT, cartComm, &status);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field left from %d to %d\n", leftRank, rank);
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
							fieldX[additionalBinNumber - i][j][k] = leftInBunemanExBuffer[bcount];
							bcount++;
					}
				}
			}
		}

		MPI_Barrier(cartComm);

		if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			MPI_Status status;
			MPI_Sendrecv(leftOutBunemanExBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_LEFT, rightInBunemanExBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_LEFT, cartComm, &status);
			if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, leftRank);
		} else if (cartCoord[0] == 0) {
			MPI_Status status;
			MPI_Recv(rightInBunemanExBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_LEFT, cartComm, &status);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Send(leftOutBunemanExBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_LEFT, cartComm);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field rigth from %d to %d\n", rightRank, rank);
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
							fieldX[xnumberAdded - 1 - additionalBinNumber + i][j][k] = rightInBunemanExBuffer[bcount];
							bcount++;
					}
				}
			}
		}

		if ((cartCoord[0] == 0) && (boundaryConditionType != PERIODIC)) {
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						//field[0][j][k] = field[2][j][k];
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							fieldX[i][j][k] = fieldX[additionalBinNumber + 1][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldX[i][j][k] = fieldX[additionalBinNumber + 1][j][k];
						}
					}
				}
			}
		}

		if ((cartCoord[0] == cartDim[0] - 1) && (boundaryConditionType != PERIODIC)) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldX[xnumberAdded - 2 - i][j][k] = fieldX[xnumberAdded - 2 - additionalBinNumber][j][k];
					}
				}
			}
		}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldX[xnumberAdded - 1 - additionalBinNumber + i][j][k] = fieldX[1 + additionalBinNumber + i][j][k];
						fieldX[i][j][k] = fieldX[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k];
					}
				}
			}
		} else {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					//field[0][j][k] = field[2][j][k];
					for (int i = 0; i <= additionalBinNumber; ++i) {
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							fieldX[i][j][k] = 0;
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldX[i][j][k] = fieldX[additionalBinNumber + 1][j][k];
						}
						//field[xnumberAdded - 1 - i][j][k] = field[xnumberAdded - 2 - additionalBinNumber][j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanExAlongY(double*** fieldX){
	if(cartDim[1] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						backOutBunemanExBuffer[bcount] = fieldX[i][ynumberAdded - additionalBinNumber - 1 - j][k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, backRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (2 + additionalBinNumber) * (xnumberAdded) * (znumberAdded + 1);
		MPI_Sendrecv(backOutBunemanExBuffer, number, MPI_DOUBLE, backRank, MPI_EFIELD_RIGHT, frontInBunemanExBuffer, number, MPI_DOUBLE, frontRank, MPI_EFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, backRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						fieldX[i][additionalBinNumber + 1 - j][k] = frontInBunemanExBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);

		bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						frontOutBunemanExBuffer[bcount] = fieldX[i][1 + additionalBinNumber + j][k];
						bcount++;
				}
			}
		}

		MPI_Sendrecv(frontOutBunemanExBuffer, number, MPI_DOUBLE, frontRank, MPI_EFIELD_LEFT, backInBunemanExBuffer, number, MPI_DOUBLE, backRank, MPI_EFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, frontRank);

		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", backRank, rank);
		for (int i = 0; i < xnumberAdded; ++i)
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						fieldX[i][ynumberAdded - additionalBinNumber - 1 + j][k] = backInBunemanExBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded + 1; ++j) {
						fieldX[i][j][k] = fieldX[i][1 + additionalBinNumber][k];
					}
				} else {
					for (int j = 0; j <= additionalBinNumber; ++j) {
						fieldX[i][ynumberAdded - additionalBinNumber + j][k] = fieldX[i][2 + additionalBinNumber + j][k];
						fieldX[i][j][k] = fieldX[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k];
					}
					fieldX[i][ynumberAdded - additionalBinNumber - 1][k] = fieldX[i][1 + additionalBinNumber][k];
				}
			}
		}
	}
}

void Simulation::exchangeBunemanExAlongZ(double*** fieldX){
	if(cartDim[2] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						topOutBunemanExBuffer[bcount] = fieldX[i][j][znumberAdded - additionalBinNumber - 1 - k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, topRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (2 + additionalBinNumber) * (xnumberAdded) * (ynumberAdded + 1);
		MPI_Sendrecv(topOutBunemanExBuffer, number, MPI_DOUBLE, topRank, MPI_EFIELD_RIGHT, bottomInBunemanExBuffer, number, MPI_DOUBLE, bottomRank, MPI_EFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, topRank);

		//MPI_Barrier(cartComm);
		bcount = 0;
		if (verbosity > 2) printf("receive general field left from %d to %d\n", bottomRank, rank);
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						fieldX[i][j][additionalBinNumber + 1 - k] = bottomInBunemanExBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);

		bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						bottomOutBunemanExBuffer[bcount] = fieldX[i][j][1 + additionalBinNumber + k];
						bcount++;
				}
			}
		}


		MPI_Sendrecv(bottomOutBunemanExBuffer, number, MPI_DOUBLE, bottomRank, MPI_EFIELD_LEFT, topInBunemanExBuffer, number, MPI_DOUBLE, topRank, MPI_EFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, bottomRank);

		//MPI_Barrier(cartComm);
		bcount = 0;
		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", topRank, rank);
		for (int i = 0; i < xnumberAdded; ++i)
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						fieldX[i][j][znumberAdded - additionalBinNumber - 1 + k] = topInBunemanExBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						fieldX[i][j][k] = fieldX[i][j][1 + additionalBinNumber];
					}
				} else {
					for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldX[i][j][znumberAdded - additionalBinNumber + k] = fieldX[i][j][2 + additionalBinNumber + k];
						fieldX[i][j][k] = fieldX[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k];
					}
					fieldX[i][j][znumberAdded - additionalBinNumber - 1] = fieldX[i][j][1 + additionalBinNumber];
				}
			}
		}
	}
}

void Simulation::exchangeBunemanEyAlongX(double*** fieldY){
	if(cartDim[0] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						rightOutBunemanEyBuffer[bcount] = fieldY[xnumberAdded - additionalBinNumber - 1 - i][j][k];
						bcount++;
				}
			}
		}
		int number = (2 + additionalBinNumber) * (ynumberAdded) * (znumberAdded + 1);
		//int numberLeft = (1 + additionalBinNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1);
		if ((cartCoord[0] < cartDim[0] - 1 && cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, rightRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
			MPI_Status status;
			MPI_Sendrecv(rightOutBunemanEyBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_RIGHT, leftInBunemanEyBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_RIGHT, cartComm, &status);
			if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, rightRank);
		} else if (cartCoord[0] == 0) {
			MPI_Send(rightOutBunemanEyBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_RIGHT, cartComm);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Status status;
			MPI_Recv(leftInBunemanEyBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_RIGHT, cartComm, &status);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field left from %d to %d\n", leftRank, rank);
			for (int i = 0; i < 2 + additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
							fieldY[additionalBinNumber + 1 - i][j][k] = leftInBunemanEyBuffer[bcount];
							bcount++;
					}
				}
			}
		}

		MPI_Barrier(cartComm);
		bcount = 0;
		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						leftOutBunemanEyBuffer[bcount] = fieldY[1 + additionalBinNumber + i][j][k];
						bcount++;
				}
			}
		}

		if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			MPI_Status status;
			MPI_Sendrecv(leftOutBunemanEyBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_LEFT, rightInBunemanEyBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_LEFT, cartComm, &status);
			if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, leftRank);
		} else if (cartCoord[0] == 0) {
			MPI_Status status;
			MPI_Recv(rightInBunemanEyBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_LEFT, cartComm, &status);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Send(leftOutBunemanEyBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_LEFT, cartComm);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field rigth from %d to %d\n", rightRank, rank);
			for (int i = 0; i < 2 + additionalBinNumber; ++i)
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
							fieldY[xnumberAdded - additionalBinNumber - 1 + i][j][k] = rightInBunemanEyBuffer[bcount];
							bcount++;
					}
				}
		}

		if ((cartCoord[0] == 0) && (boundaryConditionType != PERIODIC)) {
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						//field[0][j][k] = field[2][j][k];
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							fieldY[i][j][k] = fieldY[1 + additionalBinNumber][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldY[i][j][k] = fieldY[additionalBinNumber][j][k];
						}
					}
				}
			}
		}

		if ((cartCoord[0] == cartDim[0] - 1) && (boundaryConditionType != PERIODIC)) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldY[xnumberAdded - 1 - i][j][k] = fieldY[xnumberAdded - 1 - additionalBinNumber][j][k];
					}
				}
			}
		}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldY[xnumberAdded - additionalBinNumber + i][j][k] = fieldY[2 + additionalBinNumber + i][j][k];
						fieldY[i][j][k] = fieldY[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k];
					}
					fieldY[xnumberAdded - additionalBinNumber - 1][j][k] = fieldY[1 + additionalBinNumber][j][k];
				}
			}
		} else {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					//bunemanEy[0][j][k] = bunemanEy[2][j][k];
					for (int i = 0; i <= additionalBinNumber; ++i) {
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							fieldY[i][j][k] = fieldY[1 + additionalBinNumber][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldY[i][j][k] = fieldY[additionalBinNumber][j][k];
						}
						//field[xnumberAdded - 1 - i][j][k] = field[xnumberAdded - 1 - additionalBinNumber][j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanEyAlongY(double*** fieldY){
	if(cartDim[1] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						backOutBunemanEyBuffer[bcount] = fieldY[i][ynumberAdded - additionalBinNumber - 2 - j][k];
						frontOutBunemanEyBuffer[bcount] = fieldY[i][1 + additionalBinNumber + j][k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, backRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (1 + additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded + 1);
		MPI_Sendrecv(backOutBunemanEyBuffer, number, MPI_DOUBLE, backRank, MPI_BFIELD_RIGHT, frontInBunemanEyBuffer, number, MPI_DOUBLE, frontRank, MPI_BFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, backRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field left from %d to %d\n", frontRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						fieldY[i][additionalBinNumber - j][k] = frontInBunemanEyBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);


		MPI_Sendrecv(frontOutBunemanEyBuffer, number, MPI_DOUBLE, frontRank, MPI_BFIELD_LEFT, backInBunemanEyBuffer, number, MPI_DOUBLE, backRank, MPI_BFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, frontRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", backRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i)
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						fieldY[i][ynumberAdded - additionalBinNumber - 1 + j][k] = backInBunemanEyBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded; ++j) {
						fieldY[i][j][k] = fieldY[i][1 + additionalBinNumber][k];
					}
				} else {
					for (int j = 0; j <= additionalBinNumber; ++j) {
						fieldY[i][ynumberAdded - 1 - additionalBinNumber + j][k] = fieldY[i][1 + additionalBinNumber + j][k];
						fieldY[i][j][k] = fieldY[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanEyAlongZ(double*** fieldY){
	if(cartDim[2] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						topOutBunemanEyBuffer[bcount] = fieldY[i][j][znumberAdded - additionalBinNumber - 1 - k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, topRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (2 + additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded);
		MPI_Sendrecv(topOutBunemanEyBuffer, number, MPI_DOUBLE, topRank, MPI_EFIELD_RIGHT, bottomInBunemanEyBuffer, number, MPI_DOUBLE, bottomRank, MPI_EFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, topRank);

		//MPI_Barrier(cartComm);
		bcount = 0;
		if (verbosity > 2) printf("receive general field left from %d to %d\n", bottomRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						fieldY[i][j][additionalBinNumber + 1 - k] = bottomInBunemanEyBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);

		bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						bottomOutBunemanEyBuffer[bcount] = fieldY[i][j][1 + additionalBinNumber + k];
						bcount++;
				}
			}
		}


		MPI_Sendrecv(bottomOutBunemanEyBuffer, number, MPI_DOUBLE, bottomRank, MPI_EFIELD_LEFT, topInBunemanEyBuffer, number, MPI_DOUBLE, topRank, MPI_EFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, bottomRank);

		//MPI_Barrier(cartComm);
		bcount = 0;
		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", topRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i)
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						fieldY[i][j][znumberAdded - additionalBinNumber - 1 + k] = topInBunemanEyBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						fieldY[i][j][k] = fieldY[i][j][1 + additionalBinNumber];
					}
				} else {
					for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldY[i][j][znumberAdded - additionalBinNumber + k] = fieldY[i][j][2 + additionalBinNumber + k];
						fieldY[i][j][k] = fieldY[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k];
					}
					fieldY[i][j][znumberAdded - additionalBinNumber - 1] = fieldY[i][j][1 + additionalBinNumber];
				}
			}
		}
	}
}

void Simulation::exchangeBunemanEzAlongX(double*** fieldZ){
	if(cartDim[0] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						rightOutBunemanEzBuffer[bcount] = fieldZ[xnumberAdded - additionalBinNumber - 1 - i][j][k];
						bcount++;
				}
			}
		}
		int number = (2 + additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded);
		//int numberLeft = (1 + additionalBinNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1);
		if ((cartCoord[0] < cartDim[0] - 1 && cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, rightRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
			MPI_Status status;
			MPI_Sendrecv(rightOutBunemanEzBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_RIGHT, leftInBunemanEzBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_RIGHT, cartComm, &status);
			if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, rightRank);
		} else if (cartCoord[0] == 0) {
			MPI_Send(rightOutBunemanEzBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_RIGHT, cartComm);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Status status;
			MPI_Recv(leftInBunemanEzBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_RIGHT, cartComm, &status);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field left from %d to %d\n", leftRank, rank);
			for (int i = 0; i < 2 + additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
							fieldZ[additionalBinNumber + 1 - i][j][k] = leftInBunemanEzBuffer[bcount];
							bcount++;
					}
				}
			}
		}

		MPI_Barrier(cartComm);
		bcount = 0;
		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						leftOutBunemanEzBuffer[bcount] = fieldZ[1 + additionalBinNumber + i][j][k];
						bcount++;
				}
			}
		}

		if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			MPI_Status status;
			MPI_Sendrecv(leftOutBunemanEzBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_LEFT, rightInBunemanEzBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_LEFT, cartComm, &status);
			if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, leftRank);
		} else if (cartCoord[0] == 0) {
			MPI_Status status;
			MPI_Recv(rightInBunemanEzBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_LEFT, cartComm, &status);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Send(leftOutBunemanEzBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_LEFT, cartComm);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field rigth from %d to %d\n", rightRank, rank);
			for (int i = 0; i < 2 + additionalBinNumber; ++i)
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
							fieldZ[xnumberAdded - additionalBinNumber - 1 + i][j][k] = rightInBunemanEzBuffer[bcount];
							bcount++;
					}
				}
		}

		if ((cartCoord[0] == 0) && (boundaryConditionType != PERIODIC)) {
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						//field[0][j][k] = field[2][j][k];
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							fieldZ[i][j][k] = fieldZ[1 + additionalBinNumber][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldZ[i][j][k] = fieldZ[additionalBinNumber][j][k];
						}
					}
				}
			}
		}

		if ((cartCoord[0] == cartDim[0] - 1) && (boundaryConditionType != PERIODIC)) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldZ[xnumberAdded - 1 - i][j][k] = fieldZ[xnumberAdded - 1 - additionalBinNumber][j][k];
					}
				}
			}
		}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldZ[xnumberAdded - additionalBinNumber + i][j][k] = fieldZ[2 + additionalBinNumber + i][j][k];
						fieldZ[i][j][k] = fieldZ[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k];
					}
					fieldZ[xnumberAdded - additionalBinNumber - 1][j][k] = fieldZ[1 + additionalBinNumber][j][k];
				}
			}
		} else {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							fieldZ[i][j][k] = fieldZ[1 + additionalBinNumber][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldZ[i][j][k] = fieldZ[additionalBinNumber][j][k];
						}
						//field[xnumberAdded - 1 - i][j][k] = field[xnumberAdded - 1 - additionalBinNumber][j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanEzAlongY(double*** fieldZ){
	if(cartDim[1] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						backOutBunemanEzBuffer[bcount] = fieldZ[i][ynumberAdded - additionalBinNumber - 1 - j][k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, backRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (2 + additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded);
		MPI_Sendrecv(backOutBunemanEzBuffer, number, MPI_DOUBLE, backRank, MPI_EFIELD_RIGHT, frontInBunemanEzBuffer, number, MPI_DOUBLE, frontRank, MPI_EFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, backRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						fieldZ[i][additionalBinNumber + 1 - j][k] = frontInBunemanEzBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);

		bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						frontOutBunemanEzBuffer[bcount] = fieldZ[i][1 + additionalBinNumber + j][k];
						bcount++;
				}
			}
		}

		MPI_Sendrecv(frontOutBunemanEzBuffer, number, MPI_DOUBLE, frontRank, MPI_EFIELD_LEFT, backInBunemanEzBuffer, number, MPI_DOUBLE, backRank, MPI_EFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, frontRank);

		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", backRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i)
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						fieldZ[i][ynumberAdded - additionalBinNumber - 1 + j][k] = backInBunemanEzBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded + 1; ++j) {
						fieldZ[i][j][k] = fieldZ[i][1 + additionalBinNumber][k];
					}
				} else {
					for (int j = 0; j <= additionalBinNumber; ++j) {
						fieldZ[i][ynumberAdded - additionalBinNumber + j][k] = fieldZ[i][2 + additionalBinNumber + j][k];
						fieldZ[i][j][k] = fieldZ[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k];
					}
					fieldZ[i][ynumberAdded - additionalBinNumber - 1][k] = fieldZ[i][1 + additionalBinNumber][k];
				}
			}
		}
	}
}

void Simulation::exchangeBunemanEzAlongZ(double*** fieldZ){
	if(cartDim[2] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
						topOutBunemanEzBuffer[bcount] = fieldZ[i][j][znumberAdded - additionalBinNumber - 2 - k];
						bottomOutBunemanEzBuffer[bcount] = fieldZ[i][j][1 + additionalBinNumber + k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, topRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (1 + additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded + 1);
		MPI_Sendrecv(topOutBunemanEzBuffer, number, MPI_DOUBLE, topRank, MPI_BFIELD_RIGHT, bottomInBunemanEzBuffer, number, MPI_DOUBLE, bottomRank, MPI_BFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, topRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field left from %d to %d\n", bottomRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldZ[i][j][additionalBinNumber - k] = bottomInBunemanEzBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);


		MPI_Sendrecv(bottomOutBunemanEzBuffer, number, MPI_DOUBLE, bottomRank, MPI_BFIELD_LEFT, topInBunemanEzBuffer, number, MPI_DOUBLE, topRank, MPI_BFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, bottomRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", topRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i)
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldZ[i][j][znumberAdded - 1 - additionalBinNumber + k] = topInBunemanEzBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded; ++k) {
						fieldZ[i][j][k] = fieldZ[i][j][1 + additionalBinNumber];
					}
				} else {
					for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldZ[i][j][znumberAdded - 1 - additionalBinNumber + k] = fieldZ[i][j][1 + additionalBinNumber + k];
						fieldZ[i][j][k] = fieldZ[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanBfield(double*** fieldX, double*** fieldY, double*** fieldZ){
	exchangeBunemanBxAlongX(fieldX);
	exchangeBunemanBxAlongY(fieldX);
	exchangeBunemanBxAlongZ(fieldX);

	exchangeBunemanByAlongX(fieldY);
	exchangeBunemanByAlongY(fieldY);
	exchangeBunemanByAlongZ(fieldY);

	exchangeBunemanBzAlongX(fieldZ);
	exchangeBunemanBzAlongY(fieldZ);
	exchangeBunemanBzAlongZ(fieldZ);
}

void Simulation::exchangeBunemanBxAlongX(double*** fieldX){
	if(cartDim[0] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						rightOutBunemanBxBuffer[bcount] = fieldX[xnumberAdded - additionalBinNumber - 1 - i][j][k];
						bcount++;
				}
			}
		}
		int number = (2 + additionalBinNumber) * (ynumberAdded) * (znumberAdded);
		//int numberLeft = (1 + additionalBinNumber) * 3 * (ynumberAdded + 1) * (znumberAdded + 1);
		if ((cartCoord[0] < cartDim[0] - 1 && cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, rightRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
			MPI_Status status;
			MPI_Sendrecv(rightOutBunemanBxBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_RIGHT, leftInBunemanBxBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_RIGHT, cartComm, &status);
			if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, rightRank);
		} else if (cartCoord[0] == 0) {
			MPI_Send(rightOutBunemanBxBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_RIGHT, cartComm);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Status status;
			MPI_Recv(leftInBunemanBxBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_RIGHT, cartComm, &status);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field left from %d to %d\n", leftRank, rank);
			for (int i = 0; i < 2 + additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
							fieldX[additionalBinNumber + 1 - i][j][k] = leftInBunemanBxBuffer[bcount];
							bcount++;
					}
				}
			}
		}

		MPI_Barrier(cartComm);
		bcount = 0;
		for (int i = 0; i < 2 + additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						leftOutBunemanBxBuffer[bcount] = fieldX[1 + additionalBinNumber + i][j][k];
						bcount++;
				}
			}
		}

		if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			MPI_Status status;
			MPI_Sendrecv(leftOutBunemanBxBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_LEFT, rightInBunemanBxBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_LEFT, cartComm, &status);
			if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, leftRank);
		} else if (cartCoord[0] == 0) {
			MPI_Status status;
			MPI_Recv(rightInBunemanBxBuffer, number, MPI_DOUBLE, rightRank, MPI_EFIELD_LEFT, cartComm, &status);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Send(leftOutBunemanBxBuffer, number, MPI_DOUBLE, leftRank, MPI_EFIELD_LEFT, cartComm);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field rigth from %d to %d\n", rightRank, rank);
			for (int i = 0; i < 2 + additionalBinNumber; ++i)
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
							fieldX[xnumberAdded - additionalBinNumber - 1 + i][j][k] = rightInBunemanBxBuffer[bcount];
							bcount++;
					}
				}
		}

		if ((cartCoord[0] == 0) && (boundaryConditionType != PERIODIC)) {
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						//field[0][j][k] = field[2][j][k];
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							fieldX[i][j][k] = fieldX[1 + additionalBinNumber][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldX[i][j][k] = fieldX[additionalBinNumber][j][k];
						}
					}
				}
			}
		}

		if ((cartCoord[0] == cartDim[0] - 1) && (boundaryConditionType != PERIODIC)) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldX[xnumberAdded - 1 - i][j][k] = fieldX[xnumberAdded - 1 - additionalBinNumber][j][k];
					}
				}
			}
		}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldX[xnumberAdded - additionalBinNumber + i][j][k] = fieldX[2 + additionalBinNumber + i][j][k];
						fieldX[i][j][k] = fieldX[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k];
					}
					fieldX[xnumberAdded - additionalBinNumber - 1][j][k] = fieldX[1 + additionalBinNumber][j][k];
				}
			}
		} else {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							fieldX[i][j][k] = fieldX[1 + additionalBinNumber][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldX[i][j][k] = fieldX[additionalBinNumber][j][k];
						}
						//field[xnumberAdded - 1 - i][j][k] = field[xnumberAdded - 1 - additionalBinNumber][j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanBxAlongY(double*** fieldX){
	if(cartDim[1] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						backOutBunemanByBuffer[bcount] = fieldX[i][ynumberAdded - additionalBinNumber - 2 - j][k];
						frontOutBunemanByBuffer[bcount] = fieldX[i][1 + additionalBinNumber + j][k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, backRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (1 + additionalBinNumber) * (xnumberAdded + 1) * (znumberAdded);
		MPI_Sendrecv(backOutBunemanByBuffer, number, MPI_DOUBLE, backRank, MPI_BFIELD_RIGHT, frontInBunemanByBuffer, number, MPI_DOUBLE, frontRank, MPI_BFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, backRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field left from %d to %d\n", frontRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						fieldX[i][additionalBinNumber - j][k] = frontInBunemanByBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);


		MPI_Sendrecv(frontOutBunemanByBuffer, number, MPI_DOUBLE, frontRank, MPI_BFIELD_LEFT, backInBunemanByBuffer, number, MPI_DOUBLE, backRank, MPI_BFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, frontRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", backRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i)
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						fieldX[i][ynumberAdded - additionalBinNumber - 1 + j][k] = backInBunemanByBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded; ++j) {
						fieldX[i][j][k] = fieldX[i][1 + additionalBinNumber][k];
					}
				} else {
					for (int j = 0; j <= additionalBinNumber; ++j) {
						fieldX[i][ynumberAdded - 1 - additionalBinNumber + j][k] = fieldX[i][1 + additionalBinNumber + j][k];
						fieldX[i][j][k] = fieldX[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanBxAlongZ(double*** fieldX){
	if(cartDim[2] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
						topOutBunemanBxBuffer[bcount] = fieldX[i][j][znumberAdded - additionalBinNumber - 2 - k];
						bottomOutBunemanBxBuffer[bcount] = fieldX[i][j][1 + additionalBinNumber + k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, topRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (1 + additionalBinNumber) * (xnumberAdded + 1) * (ynumberAdded);
		MPI_Sendrecv(topOutBunemanBxBuffer, number, MPI_DOUBLE, topRank, MPI_BFIELD_RIGHT, bottomInBunemanBxBuffer, number, MPI_DOUBLE, bottomRank, MPI_BFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, topRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field left from %d to %d\n", bottomRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldX[i][j][additionalBinNumber - k] = bottomInBunemanBxBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);


		MPI_Sendrecv(bottomOutBunemanBxBuffer, number, MPI_DOUBLE, bottomRank, MPI_BFIELD_LEFT, topInBunemanBxBuffer, number, MPI_DOUBLE, topRank, MPI_BFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, bottomRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", topRank, rank);
		for (int i = 0; i < xnumberAdded + 1; ++i){
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldX[i][j][znumberAdded - 1 - additionalBinNumber + k] = topInBunemanBxBuffer[bcount];
						bcount++;
				}
			}
		}
		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded; ++k) {
						fieldX[i][j][k] = fieldX[i][j][1 + additionalBinNumber];
					}
				} else {
					for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldX[i][j][znumberAdded - 1 - additionalBinNumber + k] = fieldX[i][j][1 + additionalBinNumber + k];
						fieldX[i][j][k] = fieldX[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanByAlongX(double*** fieldY){
	if(cartDim[0] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						rightOutBunemanByBuffer[bcount] = fieldY[xnumberAdded - additionalBinNumber - 2 - i][j][k];
						leftOutBunemanByBuffer[bcount] = fieldY[1 + additionalBinNumber + i][j][k];
						bcount++;
				}
			}
		}
		int number = (1 + additionalBinNumber) * (ynumberAdded + 1) * (znumberAdded);
		if ((cartCoord[0] < cartDim[0] - 1 && cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, rightRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
			MPI_Status status;
			MPI_Sendrecv(rightOutBunemanByBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_RIGHT, leftInBunemanByBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_RIGHT, cartComm, &status);
			if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, rightRank);
		} else if (cartCoord[0] == 0) {
			MPI_Send(rightOutBunemanByBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_RIGHT, cartComm);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Status status;
			MPI_Recv(leftInBunemanByBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_RIGHT, cartComm, &status);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field left from %d to %d\n", leftRank, rank);
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
							fieldY[additionalBinNumber - i][j][k] = leftInBunemanByBuffer[bcount];
							bcount++;
					}
				}
			}
		}

		MPI_Barrier(cartComm);

		if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			MPI_Status status;
			MPI_Sendrecv(leftOutBunemanByBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_LEFT, rightInBunemanByBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_LEFT, cartComm, &status);
			if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, leftRank);
		} else if (cartCoord[0] == 0) {
			MPI_Status status;
			MPI_Recv(rightInBunemanByBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_LEFT, cartComm, &status);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Send(leftOutBunemanByBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_LEFT, cartComm);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field rigth from %d to %d\n", rightRank, rank);
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
							fieldY[xnumberAdded - 1 - additionalBinNumber + i][j][k] = rightInBunemanByBuffer[bcount];
							bcount++;
					}
				}
			}
		}

		if ((cartCoord[0] == 0) && (boundaryConditionType != PERIODIC)) {
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						//field[0][j][k] = field[2][j][k];
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							fieldY[i][j][k] = fieldY[additionalBinNumber + 1][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldY[i][j][k] = fieldY[additionalBinNumber + 1][j][k];
						}
					}
				}
			}
		}

		if ((cartCoord[0] == cartDim[0] - 1) && (boundaryConditionType != PERIODIC)) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldY[xnumberAdded - 2 - i][j][k] = fieldY[xnumberAdded - 2 - additionalBinNumber][j][k];
					}
				}
			}
		}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldY[xnumberAdded - 1 - additionalBinNumber + i][j][k] = fieldY[1 + additionalBinNumber + i][j][k];
						fieldY[i][j][k] = fieldY[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k];
					}
				}
			}
		} else {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					//field[0][j][k] = field[2][j][k];
					for (int i = 0; i <= additionalBinNumber; ++i) {
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							fieldY[i][j][k] = 0;
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldY[i][j][k] = fieldY[additionalBinNumber + 1][j][k];
						}
						//field[xnumberAdded - 1 - i][j][k] = field[xnumberAdded - 2 - additionalBinNumber][j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanByAlongY(double*** fieldY){
	if(cartDim[1] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						backOutBunemanByBuffer[bcount] = fieldY[i][ynumberAdded - additionalBinNumber - 1 - j][k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, backRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (2 + additionalBinNumber) * (xnumberAdded) * (znumberAdded);
		MPI_Sendrecv(backOutBunemanByBuffer, number, MPI_DOUBLE, backRank, MPI_EFIELD_RIGHT, frontInBunemanByBuffer, number, MPI_DOUBLE, frontRank, MPI_EFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, backRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						fieldY[i][additionalBinNumber + 1 - j][k] = frontInBunemanByBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);

		bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						frontOutBunemanByBuffer[bcount] = fieldY[i][1 + additionalBinNumber + j][k];
						bcount++;
				}
			}
		}

		MPI_Sendrecv(frontOutBunemanByBuffer, number, MPI_DOUBLE, frontRank, MPI_EFIELD_LEFT, backInBunemanByBuffer, number, MPI_DOUBLE, backRank, MPI_EFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, frontRank);

		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", backRank, rank);
		for (int i = 0; i < xnumberAdded; ++i)
			for (int j = 0; j < 2 + additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
						fieldY[i][ynumberAdded - additionalBinNumber - 1 + j][k] = backInBunemanByBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded + 1; ++j) {
						fieldY[i][j][k] = fieldY[i][1 + additionalBinNumber][k];
					}
				} else {
					for (int j = 0; j <= additionalBinNumber; ++j) {
						fieldY[i][ynumberAdded - additionalBinNumber + j][k] = fieldY[i][2 + additionalBinNumber + j][k];
						fieldY[i][j][k] = fieldY[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k];
					}
					fieldY[i][ynumberAdded - additionalBinNumber - 1][k] = fieldY[i][1 + additionalBinNumber][k];
				}
			}
		}
	}
}

void Simulation::exchangeBunemanByAlongZ(double*** fieldY){
	if(cartDim[2] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
						topOutBunemanByBuffer[bcount] = fieldY[i][j][znumberAdded - additionalBinNumber - 2 - k];
						bottomOutBunemanByBuffer[bcount] = fieldY[i][j][1 + additionalBinNumber + k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, topRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (1 + additionalBinNumber) * (xnumberAdded) * (ynumberAdded + 1);
		MPI_Sendrecv(topOutBunemanByBuffer, number, MPI_DOUBLE, topRank, MPI_BFIELD_RIGHT, bottomInBunemanByBuffer, number, MPI_DOUBLE, bottomRank, MPI_BFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, topRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field left from %d to %d\n", bottomRank, rank);
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldY[i][j][additionalBinNumber - k] = bottomInBunemanByBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);


		MPI_Sendrecv(bottomOutBunemanByBuffer, number, MPI_DOUBLE, bottomRank, MPI_BFIELD_LEFT, topInBunemanByBuffer, number, MPI_DOUBLE, topRank, MPI_BFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, bottomRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", topRank, rank);
		for (int i = 0; i < xnumberAdded; ++i)
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldY[i][j][znumberAdded - 1 - additionalBinNumber + k] = topInBunemanByBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded; ++k) {
						fieldY[i][j][k] = fieldY[i][j][1 + additionalBinNumber];
					}
				} else {
					for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldY[i][j][znumberAdded - 1 - additionalBinNumber + k] = fieldY[i][j][1 + additionalBinNumber + k];
						fieldY[i][j][k] = fieldY[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanBzAlongX(double*** fieldZ){
	if(cartDim[0] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i <= additionalBinNumber; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						rightOutBunemanBzBuffer[bcount] = fieldZ[xnumberAdded - additionalBinNumber - 2 - i][j][k];
						leftOutBunemanBzBuffer[bcount] = fieldZ[1 + additionalBinNumber + i][j][k];
						bcount++;
				}
			}
		}
		int number = (1 + additionalBinNumber) * (ynumberAdded) * (znumberAdded + 1);
		if ((cartCoord[0] < cartDim[0] - 1 && cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, rightRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
			MPI_Status status;
			MPI_Sendrecv(rightOutBunemanBzBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_RIGHT, leftInBunemanBzBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_RIGHT, cartComm, &status);
			if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, rightRank);
		} else if (cartCoord[0] == 0) {
			MPI_Send(rightOutBunemanBzBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_RIGHT, cartComm);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Status status;
			MPI_Recv(leftInBunemanBzBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_RIGHT, cartComm, &status);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] > 0) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field left from %d to %d\n", leftRank, rank);
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
							fieldZ[additionalBinNumber - i][j][k] = leftInBunemanBzBuffer[bcount];
							bcount++;
					}
				}
			}
		}

		MPI_Barrier(cartComm);

		if ((cartCoord[0] > 0 && cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			MPI_Status status;
			MPI_Sendrecv(leftOutBunemanBzBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_LEFT, rightInBunemanBzBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_LEFT, cartComm, &status);
			if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, leftRank);
		} else if (cartCoord[0] == 0) {
			MPI_Status status;
			MPI_Recv(rightInBunemanBzBuffer, number, MPI_DOUBLE, rightRank, MPI_BFIELD_LEFT, cartComm, &status);
		} else if (cartCoord[0] == cartDim[0] - 1) {
			MPI_Send(leftOutBunemanBzBuffer, number, MPI_DOUBLE, leftRank, MPI_BFIELD_LEFT, cartComm);
		}
		//MPI_Barrier(cartComm);
		bcount = 0;
		if ((cartCoord[0] < cartDim[0] - 1) || (boundaryConditionType == PERIODIC)) {
			if (verbosity > 2) printf("receive general field rigth from %d to %d\n", rightRank, rank);
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
							fieldZ[xnumberAdded - 1 - additionalBinNumber + i][j][k] = rightInBunemanBzBuffer[bcount];
							bcount++;
					}
				}
			}
		}

		if ((cartCoord[0] == 0) && (boundaryConditionType != PERIODIC)) {
			for (int i = 0; i <= additionalBinNumber; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						//field[0][j][k] = field[2][j][k];
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							//field[i][j][k] = Vector3d(0, 0, 0);
							fieldZ[i][j][k] = fieldZ[additionalBinNumber + 1][j][k];
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldZ[i][j][k] = fieldZ[additionalBinNumber + 1][j][k];
						}
					}
				}
			}
		}

		if ((cartCoord[0] == cartDim[0] - 1) && (boundaryConditionType != PERIODIC)) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldZ[xnumberAdded - 2 - i][j][k] = fieldZ[xnumberAdded - 2 - additionalBinNumber][j][k];
					}
				}
			}
		}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber; ++i) {
						fieldZ[xnumberAdded - 1 - additionalBinNumber + i][j][k] = fieldZ[1 + additionalBinNumber + i][j][k];
						fieldZ[i][j][k] = fieldZ[xnumberAdded - 2 - 2 * additionalBinNumber + i][j][k];
					}
				}
			}
		} else {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					//field[0][j][k] = field[2][j][k];
					for (int i = 0; i <= additionalBinNumber; ++i) {
						if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
							fieldZ[i][j][k] = 0;
						} else if (boundaryConditionType == FREE_BOTH) {
							fieldZ[i][j][k] = fieldZ[additionalBinNumber + 1][j][k];
						}
						//field[xnumberAdded - 1 - i][j][k] = field[xnumberAdded - 2 - additionalBinNumber][j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanBzAlongY(double*** fieldZ){
	if(cartDim[1] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						backOutBunemanBzBuffer[bcount] = fieldZ[i][ynumberAdded - additionalBinNumber - 2 - j][k];
						frontOutBunemanBzBuffer[bcount] = fieldZ[i][1 + additionalBinNumber + j][k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, backRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (1 + additionalBinNumber) * (xnumberAdded) * (znumberAdded + 1);
		MPI_Sendrecv(backOutBunemanBzBuffer, number, MPI_DOUBLE, backRank, MPI_BFIELD_RIGHT, frontInBunemanBzBuffer, number, MPI_DOUBLE, frontRank, MPI_BFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, backRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field left from %d to %d\n", frontRank, rank);
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						fieldZ[i][additionalBinNumber - j][k] = frontInBunemanBzBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);


		MPI_Sendrecv(frontOutBunemanBzBuffer, number, MPI_DOUBLE, frontRank, MPI_BFIELD_LEFT, backInBunemanBzBuffer, number, MPI_DOUBLE, backRank, MPI_BFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, frontRank);

		//MPI_Barrier(cartComm);
		bcount = 0;

		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", backRank, rank);
		for (int i = 0; i < xnumberAdded; ++i){
			for (int j = 0; j <= additionalBinNumber; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
						fieldZ[i][ynumberAdded - additionalBinNumber - 1 + j][k] = backInBunemanBzBuffer[bcount];
						bcount++;
				}
			}
	}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				if (ynumberGeneral == 1) {
					for (int j = 0; j < ynumberAdded; ++j) {
						fieldZ[i][j][k] = fieldZ[i][1 + additionalBinNumber][k];
					}
				} else {
					for (int j = 0; j <= additionalBinNumber; ++j) {
						fieldZ[i][ynumberAdded - 1 - additionalBinNumber + j][k] = fieldZ[i][1 + additionalBinNumber + j][k];
						fieldZ[i][j][k] = fieldZ[i][ynumberAdded - 2 - 2 * additionalBinNumber + j][k];
					}
				}
			}
		}
	}
}

void Simulation::exchangeBunemanBzAlongZ(double*** fieldZ){
	if(cartDim[2] > 1){
		if (verbosity > 2) printf("start sending general field rank = %d\n", rank);
		int bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						topOutBunemanBzBuffer[bcount] = fieldZ[i][j][znumberAdded - additionalBinNumber - 1 - k];
						bcount++;
				}
			}
		}

		if (verbosity > 2) printf("before send general field right from %d to %d additionalBinNumber = %d, xnumber = %d, ynumber = %d, znumber = %d\n", rank, topRank, additionalBinNumber, xnumberAdded, ynumberAdded, znumberAdded);
		MPI_Status status;
		int number = (2 + additionalBinNumber) * (xnumberAdded) * (ynumberAdded);
		MPI_Sendrecv(topOutBunemanBzBuffer, number, MPI_DOUBLE, topRank, MPI_EFIELD_RIGHT, bottomInBunemanBzBuffer, number, MPI_DOUBLE, bottomRank, MPI_EFIELD_RIGHT, cartComm, &status);
		if (verbosity > 2) printf("after send general field right from %d to %d\n", rank, topRank);

		//MPI_Barrier(cartComm);
		bcount = 0;
		if (verbosity > 2) printf("receive general field left from %d to %d\n", bottomRank, rank);
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						fieldZ[i][j][additionalBinNumber + 1 - k] = bottomInBunemanBzBuffer[bcount];
						bcount++;
				}
			}
		}


		MPI_Barrier(cartComm);

		bcount = 0;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						bottomOutBunemanBzBuffer[bcount] = fieldZ[i][j][1 + additionalBinNumber + k];
						bcount++;
				}
			}
		}


		MPI_Sendrecv(bottomOutBunemanBzBuffer, number, MPI_DOUBLE, bottomRank, MPI_EFIELD_LEFT, topInBunemanBzBuffer, number, MPI_DOUBLE, topRank, MPI_EFIELD_LEFT, cartComm, &status);
		if (verbosity > 2) printf("send general fieldleft from %d to %d\n", rank, bottomRank);

		//MPI_Barrier(cartComm);
		bcount = 0;
		if (verbosity > 2) printf("receive general field rigth from %d to %d\n", topRank, rank);
		for (int i = 0; i < xnumberAdded; ++i)
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < 2 + additionalBinNumber; ++k) {
						fieldZ[i][j][znumberAdded - additionalBinNumber - 1 + k] = topInBunemanBzBuffer[bcount];
						bcount++;
				}
			}

		if (verbosity > 2) printf("finish exchanging E field rank = %d\n", rank);
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				if (znumberGeneral == 1) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						fieldZ[i][j][k] = fieldZ[i][j][1 + additionalBinNumber];
					}
				} else {
					for (int k = 0; k <= additionalBinNumber; ++k) {
						fieldZ[i][j][znumberAdded - additionalBinNumber + k] = fieldZ[i][j][2 + additionalBinNumber + k];
						fieldZ[i][j][k] = fieldZ[i][j][znumberAdded - 2 - 2 * additionalBinNumber + k];
					}
					fieldZ[i][j][znumberAdded - additionalBinNumber - 1] = fieldZ[i][j][1 + additionalBinNumber];
				}
			}
		}
	}
}