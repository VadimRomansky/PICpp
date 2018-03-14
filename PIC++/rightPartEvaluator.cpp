#include "mpi.h"
#include "math.h"
#include "stdlib.h"
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "rightPartEvaluator.h"
#include "paths.h"
#include "specialmath.h"
#include "mpi_util.h"

RightPartIterativeEvaluator::RightPartIterativeEvaluator(int xn, int yn, int zn, int ln, int rankv, int nprocsv, MPI_Comm& cartCommv, int* cartCoordv, int* cartDimv) {
	rank = rankv;
	nprocs = nprocsv;
	xnumberAdded = xn;
	ynumberAdded = yn;
	znumberAdded = zn;
	lnumber = ln;
	cartComm = cartCommv;
	cartCoord = cartCoordv;
	cartDim = cartDimv;
}

double RightPartIterativeEvaluator::rightPart(double**** vector, int i, int j, int k, int l) {
	std::string outputDir = outputDirectory;
	FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
	fprintf(errorLogFile, "Using base class for rightPartEvaluator is restricted\n");
	printf("Using base class for rightPartEvaluator is restricted\n");
	fflush(stdout);
	fclose(errorLogFile);
	MPI_Finalize();
	exit(0);
}

void RightPartIterativeEvaluator::getError(double**** vector, double**** errorVector) {
	for(int i = 0; i < xnumberAdded; ++i) {
		for(int j = 0; j < ynumberAdded; ++j) {
			for(int k = 0; k < znumberAdded; ++k) {
				for(int l = 0; l < lnumber; ++l){
					errorVector[i][j][k][l] = rightPart(vector, i, j, k, l) - vector[i][j][k][l];
				}
			}
		}
	}
	
}

InitializeShockRightPartEvaluator::InitializeShockRightPartEvaluator(Vector3d v1, Vector3d v2, Vector3d B1, Vector3d B2,
                                                                     Vector3d E1, double dx, double cvalue,
                                                                     double fNorm, double vNorm,
                                                                     ParticleTypeContainer* types, int tn, int xn,
                                                                     int yn, int zn, int ln, int rank, int nprocs, MPI_Comm& cartComm, int* cartCoord, int* cartDim) : RightPartIterativeEvaluator(
	xn, yn, zn, ln, rank, nprocs, cartComm, cartCoord, cartDim) {
	downstreamV = v1;
	upstreamV = v2;
	downstreamB = B1;
	upstreamB = B2;
	E = E1;
	deltaX = dx;
	c = cvalue;
	fieldNorm = fNorm;
	velocityNorm = vNorm;
	typesNumber = tn;

	charge = new double[typesNumber];
	mass = new double[typesNumber];
	flux = new double[typesNumber];
	for (int i = 0; i < typesNumber; ++i) {
		charge[i] = types[i].charge;
		mass[i] = types[i].mass;
		flux[i] = types[i].concentration * upstreamV.x;
	}
}

InitializeShockRightPartEvaluator::~InitializeShockRightPartEvaluator() {
	delete[] charge;
	delete[] mass;
	delete[] flux;
}

double InitializeShockRightPartEvaluator::evaluateVx(double**** vector, int i, int l) {
	int t = l / 3;
	if (t == 0) {
		double rho = 0;
		for (int typeN = 1; typeN < typesNumber; ++typeN) {
			if (fabs(flux[typeN]) > 0) {
				rho += charge[typeN] * flux[typeN] / vector[i][0][0][3 * typeN];
			}
			return -charge[0] * flux[0] / rho;
		}
	} else {
		double field = (vector[i][0][0][lnumber - 3] + velocityNorm * (vector[i][0][0][l + 1] * vector[i][0][0][lnumber - 1] -
			vector[i][0][0][l + 2] * vector[i][0][0][lnumber - 2]) / c) * fieldNorm;
		return vector[i - 1][0][0][l] + deltaX * (charge[t] / mass[t]) * field / (vector[i][0][0][l] * velocityNorm *
			velocityNorm);
	}
}

double InitializeShockRightPartEvaluator::evaluateVy(double**** vector, int i, int l) {
	int t = l / 3;
	double field = (E.y + velocityNorm * (vector[i][0][0][l + 1] * upstreamB.x - vector[i][0][0][l - 1] * vector[i][0][0][
		lnumber - 1]) / c) * fieldNorm;
	return vector[i - 1][0][0][l] + deltaX * (charge[t] / mass[t]) * field / (vector[i][0][0][l - 1] * velocityNorm *
		velocityNorm);
}

double InitializeShockRightPartEvaluator::evaluateVz(double**** vector, int i, int l) {
	int t = l / 3;
	double field = (E.z + velocityNorm * (vector[i][0][0][l - 2] * vector[i][0][0][lnumber - 2] - vector[i][0][0][l - 1] *
		upstreamB.x) / c) * fieldNorm;
	return vector[i - 1][0][0][l] + deltaX * (charge[t] / mass[t]) * field / (vector[i][0][0][l - 2] * velocityNorm *
		velocityNorm);
}

double InitializeShockRightPartEvaluator::evaluateEx(double**** vector, int i) {
	/*double sum = 0;
	for(int t = 0; t < typesNumber; ++t) {
		sum += flux[t]*charge[t]/(vector[i-1][0][0][3*t]*velocityNorm);
	}
	return vector[i-1][0][0][lnumber-3] + 4*pi*sum*deltaX;*/
	return 0;
}

double InitializeShockRightPartEvaluator::evaluateBy(double**** vector, int i) {
	double sum = 0;
	for (int t = 0; t < typesNumber; ++t) {
		sum += flux[t] * charge[t] * vector[i - 1][0][0][3 * t + 2] / vector[i - 1][0][0][3 * t];
	}
	return vector[i - 1][0][0][lnumber - 2] - 4 * pi * sum * deltaX / c;
}

double InitializeShockRightPartEvaluator::evaluateBz(double**** vector, int i) {
	double sum = 0;
	for (int t = 0; t < typesNumber; ++t) {
		sum += flux[t] * charge[t] * vector[i - 1][0][0][3 * t + 1] / vector[i - 1][0][0][3 * t];
	}
	return vector[i - 1][0][0][lnumber - 1] + 4 * pi * sum * deltaX / c;
}

double InitializeShockRightPartEvaluator::rightPart(double**** vector, int i, int j, int k, int l) {
	if (i == 0) {
		if (rank == 0) {
			if (l < 3 * typesNumber) {
				if (l % 3 == 0) {
					return downstreamV.x / velocityNorm;
				}
				if (l % 3 == 1) {
					return downstreamV.y / velocityNorm;
				}
				return downstreamV.z / velocityNorm;
			}
			if (l == 3 * typesNumber) {
				return E.x / fieldNorm;
			}
			if (l == 3 * typesNumber + 1) {
				return downstreamB.y / fieldNorm;
			}
			return downstreamB.z / fieldNorm;
		}
		return 0;
	}
	if (i == xnumberAdded) {
		if (rank == nprocs - 1) {
			if (l < 3 * typesNumber) {
				if (l % 3 == 0) {
					return upstreamV.x / velocityNorm;
				}
				if (l % 3 == 1) {
					return upstreamV.y / velocityNorm;
				}
				return upstreamV.z / velocityNorm;
			}
			if (l == 3 * typesNumber) {
				return E.x / fieldNorm;
			}
			if (l == 3 * typesNumber + 1) {
				return upstreamB.y / fieldNorm;
			}
			return upstreamB.z / fieldNorm;
		}
		return 0;
	}
	if (l < 3 * typesNumber) {
		if (l % 3 == 0) {
			return evaluateVx(vector, i, l);
		}
		if (l % 3 == 1) {
			return evaluateVy(vector, i, l);
		}
		return evaluateVz(vector, i, l);
	}
	if (l == 3 * typesNumber) {
		return evaluateEx(vector, i);
	}
	if (l == 3 * typesNumber + 1) {
		return evaluateBy(vector, i);
	}
	return evaluateBz(vector, i);
}

//todo
double InitializeShockRightPartEvaluator::rightPartInitialNorm() {
	return 1.0;
}


PoissonRightPartEvaluator::PoissonRightPartEvaluator(double*** densityv, double**** tempLargeVector, double dx, double dy, double dz, int xn, int yn, int zn, bool periodicXv, bool periodicYv, bool periodicZv, int rankv, int nprocsv,
	                                         MPI_Comm& cartCommv, int* cartCoordv, int* cartDimv): RightPartIterativeEvaluator(xn, yn, zn, 1, rankv, nprocsv, cartCommv, cartCoordv, cartDimv) {

	deltaX = dx;
	deltaY = dy;
	deltaZ = dz;

	xnumberAdded = xn;
	ynumberAdded = yn;
	znumberAdded = zn;

	density = tempLargeVector;

	double factor = deltaX*deltaX*deltaY*deltaY*deltaZ*deltaZ/(2*(deltaX*deltaX*deltaY*deltaY + deltaY*deltaY*deltaZ*deltaZ + deltaZ*deltaZ*deltaX*deltaX));
	for(int i = 0; i < xnumberAdded; ++i) {
		for(int j = 0; j < ynumberAdded; ++j) {
			for(int k = 0; k < znumberAdded; ++k) {
				density[i][j][k][0] = densityv[i][j][k]*factor;
			}
		}
	}
}

double PoissonRightPartEvaluator::rightPart(double**** vector, int i, int j, int k, int l) {
	double denom = 2*(deltaX*deltaX*deltaY*deltaY + deltaY*deltaY*deltaZ*deltaZ + deltaZ*deltaZ*deltaX*deltaX);
	double factorx = deltaY*deltaY*deltaZ*deltaZ/denom;
	double factory = deltaX*deltaX*deltaZ*deltaZ/denom;
	double factorz = deltaY*deltaY*deltaX*deltaX/denom;

	if(i < 1 + additionalBinNumber || i >= xnumberAdded - 1 - additionalBinNumber || j < 1 + additionalBinNumber || j >= ynumberAdded - 1 - additionalBinNumber || k < 1 + additionalBinNumber || k >= znumberAdded - 1 - additionalBinNumber) {
		return 0;
	}

	double result = factorx*(vector[i-1][j][k][l] + vector[i+1][j][k][l]) +
		factory*(vector[i][j-1][k][l] + vector[i][j+1][k][l]) +
		factorz*(vector[i][j][k-1][l] + vector[i][j][k+1][l]) -
		density[i][j][k][l];

	return result;
}

double PoissonRightPartEvaluator::rightPartInitialNorm() {
	double density_norm = scalarMultiplyLargeVectors(density, density, xnumberAdded, ynumberAdded, znumberAdded,
	                                         additionalBinNumber, lnumber, periodicX, periodicY, periodicZ, rank, nprocs,
	                                         cartComm, cartCoord, cartDim);
	return sqrt(density_norm);
}



