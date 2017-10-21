#include "mpi.h"
#include "math.h"
#include "stdlib.h"
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "rightPartEvaluator.h"
#include "paths.h"

BaseRightPartEvaluator::BaseRightPartEvaluator(int xn, int yn, int zn, int ln) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	xnumber = xn;
	ynumber = yn;
	znumber = zn;
	lnumber = ln;
}

double BaseRightPartEvaluator::rightPart(double**** vector, int i, int j, int k, int l) {
	std::string outputDir = outputDirectory;
    FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
    fprintf(errorLogFile, "Using base class for rightPartEvaluator is restricted\n");
    printf("Using base class for rightPartEvaluator is restricted\n");
    fflush(stdout);
    fclose(errorLogFile);
    MPI_Finalize();
    exit(0);
}

InitializeShockRightPartEvaluator::InitializeShockRightPartEvaluator(Vector3d v1, Vector3d v2, Vector3d B1, Vector3d B2, Vector3d E1, double dx, double cvalue, double fNorm, double vNorm, ParticleTypeContainer* types, int tn, int xn, int yn, int zn, int ln) : BaseRightPartEvaluator(xn, yn, zn, ln) {
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
	int t = l/3;
	if(t == 0) {
		double rho = 0;
		for(int typeN = 1; typeN < typesNumber; ++typeN) {
			if(fabs(flux[typeN]) > 0) {
				rho += charge[typeN]*flux[typeN]/vector[i][0][0][3*typeN];
			}
			return -charge[0]*flux[0]/rho;
		}
	} else {
		double field = (vector[i][0][0][lnumber - 3] + velocityNorm*(vector[i][0][0][l+1]*vector[i][0][0][lnumber - 1] - vector[i][0][0][l+2]*vector[i][0][0][lnumber - 2])/c)*fieldNorm;
		return vector[i-1][0][0][l] + deltaX*(charge[t]/mass[t])*field/(vector[i][0][0][l]*velocityNorm*velocityNorm);
	}
}

double InitializeShockRightPartEvaluator::evaluateVy(double**** vector, int i, int l) {
	int t = l/3;
	double field = (E.y+ velocityNorm*(vector[i][0][0][l+1]*upstreamB.x - vector[i][0][0][l-1]*vector[i][0][0][lnumber - 1])/c)*fieldNorm;
	return vector[i-1][0][0][l] + deltaX*(charge[t]/mass[t])*field/(vector[i][0][0][l-1]*velocityNorm*velocityNorm);
}

double InitializeShockRightPartEvaluator::evaluateVz(double**** vector, int i, int l) {
	int t = l/3;
	double field = (E.z+ velocityNorm*(vector[i][0][0][l-2]*vector[i][0][0][lnumber-2] - vector[i][0][0][l-1]*upstreamB.x)/c)*fieldNorm;
	return vector[i-1][0][0][l] + deltaX*(charge[t]/mass[t])*field/(vector[i][0][0][l-2]*velocityNorm*velocityNorm);
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
	for(int t = 0; t < typesNumber; ++t) {
		sum += flux[t]*charge[t]*vector[i-1][0][0][3*t+2]/vector[i-1][0][0][3*t];
	}
	return vector[i-1][0][0][lnumber-2] - 4*pi*sum*deltaX/c;
}

double InitializeShockRightPartEvaluator::evaluateBz(double**** vector, int i) {
	double sum = 0;
	for(int t = 0; t < typesNumber; ++t) {
		sum += flux[t]*charge[t]*vector[i-1][0][0][3*t+1]/vector[i-1][0][0][3*t];
	}
	return vector[i-1][0][0][lnumber-1] + 4*pi*sum*deltaX/c;
}

double InitializeShockRightPartEvaluator::rightPart(double**** vector, int i, int j, int k, int l) {
	if (i == 0) {
		if (rank == 0) {
			if (l < 3 * typesNumber) {
				if (l % 3 == 0) {
					return downstreamV.x/velocityNorm;
				} else if (l % 3 == 1) {
					return downstreamV.y/velocityNorm;
				} else {
					return downstreamV.z/velocityNorm;
				}
			} else if (l == 3 * typesNumber) {
				return E.x/fieldNorm;
			} else if (l == 3 * typesNumber + 1) {
				return downstreamB.y/fieldNorm;
			} else {
				return downstreamB.z/fieldNorm;
			}
		} else {
			return 0;
		}
	} else if (i == xnumber) {
		if (rank == nprocs - 1) {
			if (l < 3 * typesNumber) {
				if (l % 3 == 0) {
					return upstreamV.x/velocityNorm;
				} else if (l % 3 == 1) {
					return upstreamV.y/velocityNorm;
				} else {
					return upstreamV.z/velocityNorm;
				}
			} else if (l == 3 * typesNumber) {
				return E.x/fieldNorm;
			} else if (l == 3 * typesNumber + 1) {
				return upstreamB.y/fieldNorm;
			} else {
				return upstreamB.z/fieldNorm;
			}
		} else {
			return 0;
		}
	} else {
		if (l < 3 * typesNumber) {
			if (l % 3 == 0) {
				return evaluateVx(vector, i, l);
			} else if (l % 3 == 1) {
				return evaluateVy(vector, i, l);
			} else {
				return evaluateVz(vector, i, l);
			}
		} else if (l == 3 * typesNumber) {
			return evaluateEx(vector, i);
		} else if (l == 3 * typesNumber + 1) {
			return evaluateBy(vector, i);
		} else {
			return evaluateBz(vector, i);
		}
	}
}

