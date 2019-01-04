#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>
#include <mpi.h>
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "dichotomousSolver.h"
#include "specialmath.h"
#include "util.h"
#include "output.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "matrixElement.h"
#include "particle.h"
#include "random.h"
#include "simulation.h"
#include "mpi_util.h"
#include "input.h"
#include "complex.h"
#include "paths.h"
#include "creating_arrays.h"

void Simulation::initializeFluxFromRight() {
	boundaryConditionTypeX = SUPER_CONDUCTOR_LEFT;
	//boundaryConditionTypeX = PERIODIC;
	boundaryConditionTypeY = PERIODIC;
	boundaryConditionTypeZ = PERIODIC;
	createParticles();
	//E0 = E0 - V0.vectorMult(B0) / (speed_of_light_normalized * speed_of_light_correction);
	//E0 = E0 - V0.vectorMult(B0) / (speed_of_light_normalized);
	E0 = E0 - V0.vectorMult(B0);
	rightBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(E0, B0);
	//initializeAlfvenWaveY(10, 1.0E-4);
	if (solverType == BUNEMAN) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEx[i][j][k] = E0.x;
					bunemanNewEx[i][j][k] = bunemanEx[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanEy[i][j][k] = E0.y;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i < 1 + additionalBinNumber) {
							bunemanEy[i][j][k] = 0;
						}
					}
					bunemanNewEy[i][j][k] = bunemanEy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanEz[i][j][k] = E0.z;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i < 1 + additionalBinNumber) {
							bunemanEz[i][j][k] = 0;
						}
					}
					bunemanNewEz[i][j][k] = bunemanEz[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBx[i][j][k] = B0.x;
					bunemanNewBx[i][j][k] = bunemanBx[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					bunemanBy[i][j][k] = B0.y;
					bunemanNewBy[i][j][k] = bunemanBy[i][j][k];
				}
			}
		}
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					bunemanBz[i][j][k] = B0.z;
					bunemanNewBz[i][j][k] = bunemanBz[i][j][k];
				}
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					Efield[i][j][k] = E0;
					if (boundaryConditionTypeX != PERIODIC) {
						if (cartCoord[0] == 0 && i <= 1 + additionalBinNumber) {
							Efield[i][j][k].y = 0;
							Efield[i][j][k].z = 0;
						}
					}
					tempEfield[i][j][k] = Efield[i][j][k];
					newEfield[i][j][k] = Efield[i][j][k];
					explicitEfield[i][j][k] = Efield[i][j][k];
				}
			}
		}
	}

	//turbulence
	MPI_Barrier(cartComm);
	if(rank == 0) {
		printf("start initialize turbulence\n");
	}
	initializeRandomModes(5, 20, 0.5);

	//double gamma = 1.0 / sqrt(1 - V0.scalarMult(V0) / speed_of_light_normalized_sqr);
	double gamma = 1.0 / sqrt(1 - V0.scalarMult(V0));
	double p0 = gamma * massProton * V0.norm();
	//double protonGyroRadius = p0 * speed_of_light_normalized / (fabs(electron_charge_normalized) * B0.norm());
	double protonGyroRadius = p0 / (fabs(electron_charge_normalized) * B0.norm());
	int countGyroRadius = xsizeGeneral / protonGyroRadius;
	double deltaK = 2 * pi / xsizeGeneral;
	int maxCount = min2(2 * countGyroRadius, xnumberGeneral);
	int minCount = max2(1, countGyroRadius / 2);

	double magneticEnergy = B0.scalarMult(B0) / (8 * pi);
	double kineticEnergy = density * V0.scalarMult(V0) / 2;

	if (rank == 0) informationFile = fopen((outputDir + "information.dat").c_str(), "a");
	if (rank == 0) fprintf(informationFile, "magneticEnergy/kineticEnergy = %15.10g\n", magneticEnergy / kineticEnergy);
	if (rank == 0) fprintf(informationFile, "magneticEnergy/kinetikEnergy = %15.10g\n", magneticEnergy / kineticEnergy);
	if (rank == 0) fprintf(informationFile, "protonGyroRadius = %15.10g\n", protonGyroRadius);
	if (rank == 0) fprintf(informationFile, "countGyroRadius = %d\n", countGyroRadius);
	if (rank == 0) fprintf(informationFile, "minCount = %d\n", minCount);
	if (rank == 0) fprintf(informationFile, "maxCount = %d\n", maxCount);
	fflush(stdout);
	if (rank == 0) fclose(informationFile);

	checkDebyeParameter();

	MPI_Barrier(cartComm);
	if (rank == 0) {
		incrementFile = fopen((outputDir + "increment.dat").c_str(), "w");
		fprintf(incrementFile, "%g %g %g %g\n", 0.0, 0.0, 1.0, 1.0);
		fclose(incrementFile);
	}
	MPI_Barrier(cartComm);
}