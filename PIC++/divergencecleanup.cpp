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


void Simulation::cleanupDivergence() {
	double procTime = 0;
	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if((rank == 0) && (verbosity > 0)) printf("cleaning up divergence\n");
	fflush(stdout);

	//substractMeanDensity();

	if(ynumber == 1 && znumber == 1){
		cleanupDivergence1d();
		if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("cleaning divergence time = %g sec\n", procTime/CLOCKS_PER_SEC);
		}
		return;
	}


	int matrixDimension = xnumber;

	double fullDensity = 0;
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				fullDensity += chargeDensity[i][j][k] * volumeB(i, j, k);
			}
		}
	}
	fullDensity /= (xsize*ysize*zsize);
	if((rank == 0) && (verbosity > 1)) printf("full density = %22.15g\n", fullDensity);
	if((rank == 0) && (verbosity > 1)) printf("density[0][0][0] = %22.15g\n", chargeDensity[0][0][0]);
	fflush(stdout);


	/*for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				chargeDensity[i][j][k] -= fullDensity;
			}
		}
	}*/

	if(boundaryConditionType == PERIODIC){
		/*for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					chargeDensity[i][j][k] -= fullDensity;
				}
			}
		}*/
		/*for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					fullDensity += chargeDensity[i][j][k] * volumeB(i, j, k);
				}
			}
		}
		fullDensity /= (xsize*ysize*zsize);*/
	} else {
		//double Elinear = -4*pi*fullDensity*xsize + newEfield[xnumber][0][0].x - newEfield[0][0][0].x;
		/*double Elinear = 0;
		for(int i = 0; i < xnumber + 1; ++i){
			for(int j = 0; j< ynumber + 1; ++j){
				for(int k = 0; k < znumber + 1; ++k){
					double factor = (xgrid[xnumber] - xgrid[i])/xsize;
					newEfield[i][j][k].x = newEfield[i][j][k].x + Elinear*(factor - 1.0);
				}
			}
		}*/
	}

	/*for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				chargeDensity[i][j][k] -= fullDensity;
			}
		}
	}*/

	bool fourier = false;

	if(!fourier) {
		double foolRightPart = 0;
		for (int i = 0; i < xnumber+1; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
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
					/*if(rank == 0 && i == 1 && j == 0 && k == 0){
						divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
						divergenceCleanUpRightPart[i][j][k][0] = 0;
					} else if(rank == nprocs - 1 && i == xnumber-1 && j == 0 && k == 0){
						divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
						divergenceCleanUpRightPart[i][j][k][0] = 0;
					}else {*/
					if(i == 0){
                        createDivergenceCleanupLeftFakeEquation(0, j, k);
					} else if(i == 1){
						if(rank == 0 && boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
                            //todo check
							//createDivergenceCleanupSuperConductorEquation(i, j, k);
                            createDivergenceCleanupLeftFakeEquation(1, j, k);
						} else {
							createDivergenceCleanupInternalEquation(i, j, k);
						}
					} else if(i < xnumber - additionalBinNumber - 1) {
						createDivergenceCleanupInternalEquation(i, j, k);
					} else if(i < xnumber){
						if(rank == nprocs - 1 && boundaryConditionType != PERIODIC) {
							createDivergenceCleanupRightFakeEquation(i, j, k);
						} else {
							createDivergenceCleanupInternalEquation(i, j, k);
						}
					} else {
						createDivergenceCleanupRightFakeEquation(i, j, k);
					}	
					foolRightPart += divergenceCleanUpRightPart[i][j][k][0];
				}
			}
		}

        /*foolRightPart /= (xnumber+1)*ynumber*znumber;
        for(int i = 0; i < xnumber + 1; ++i){
            for(int j = 0; j < ynumber; ++j){
                for(int k = 0; k < znumber; ++k){
                    divergenceCleanUpRightPart[i][j][k][0] -= foolRightPart;
                }
            }
        }*/

        //at first we solve poison equationand then diffusion equation
        //debug only for 1d 1 thread
		/*FILE* matrixFile = fopen((outputDir + "matrix.dat").c_str(), "w");
        for(int i = 0; i < xnumber + 1; ++i){
            for(int j = 0; j < xnumber + 1; ++j){
                bool f = false;
                for(int k = 0; k < divergenceCleanUpMatrix[i][0][0][0].size(); ++k){
                    if(divergenceCleanUpMatrix[i][0][0][0][k].i == j){
                        fprintf(matrixFile, "%g ", divergenceCleanUpMatrix[i][0][0][0][k].value);
                        f = true;
                    }
                }
                if(!f){
                    fprintf(matrixFile, "%g ", 0.0);
                }
            }
            fprintf(matrixFile, "\n");
        }
		fclose(matrixFile);
        FILE* rightPartFile = fopen((outputDir + "rightPart.dat").c_str(), "w");
        for(int i = 0; i < xnumber + 1; ++i){
            fprintf(rightPartFile, "%15.10g\n", divergenceCleanUpRightPart[i][0][0][0]);
        }
        fclose(rightPartFile);*/


		/*double summRightPart = 0;
        for(int i = 0; i < xnumber; ++i){
            summRightPart += divergenceCleanUpRightPart[i][0];
        }*/
		//double rightPart = cleanUpRightPart(1);
		//conjugateGradientMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumber, ynumber, znumber, 1, 1E-10, xnumber * ynumber * znumber, false, verbosity);
		/*biconjugateStabilizedGradientMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumber, ynumber, znumber, 1, 1E-10, maxGMRESIterations, false, verbosity);*/
		/*generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart,
										 divergenceCleaningPotential, xnumber, ynumber, znumber, 1, xnumberGeneral,
										 ynumberGeneral, znumberGeneral, 1E-8, maxGMRESIterations, false, verbosity);*/

        gaussSeidelMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart,
                          divergenceCleaningPotential, xnumber, ynumber, znumber, 1, xnumberGeneral,
                          ynumberGeneral, znumberGeneral, 1E-8, 100, false, verbosity);
        //fake diffusion Coeffficient = 1;
       /* FILE* solutionFile = fopen((outputDir + "solution.dat").c_str(), "w");
        for(int i = 0; i < xnumber + 1; ++i){
            fprintf(solutionFile, "%15.10g\n", divergenceCleaningPotential[i][0][0][0]);
        }
        fclose(solutionFile);*/


        /*double fakeDeltaT = 1*deltaX2;
       // double fakeDeltaT = 0.1*deltaX2;
        for(int i = 1; i < xnumber; ++i){
            if(i < xnumber - additionalBinNumber || rank != nprocs - 1 || boundaryConditionType == PERIODIC) {
                for (int j = 0; j < ynumber; ++j) {
                    for (int k = 0; k < znumber; ++k) {
                        for (int m = 0; m < divergenceCleanUpMatrix[i][j][k][0].size(); ++m) {
                            divergenceCleanUpMatrix[i][j][k][0][m].value *= -fakeDeltaT;
                            if (indexEqual(divergenceCleanUpMatrix[i][j][k][0][m], i, j, k, 0)) {
                                divergenceCleanUpMatrix[i][j][k][0][m].value += 1.0;
                            }
                        }
                    }
                }
            }
        }
        double fakeT = 0;
        double maxFakeT = xsizeGeneral*xsizeGeneral;
        int maxStep = 2*maxFakeT/fakeDeltaT;
        for(int t = 0; t < maxStep; ++t){
            if(rank == 0 && verbosity > 1) printf("diffusion iteration number %d\n", t);
            for(int i = 0; i < xnumber + 1; ++i){
                for(int j = 0; j < ynumber; ++j){
                    for(int k = 0; k < znumber; ++k){
                        if(i == 0){
                            //todo check
                            divergenceCleanUpRightPart[i][j][k][0] = 0;
                        } else if(i < xnumber - additionalBinNumber - 1) {
                            divergenceCleanUpRightPart[i][j][k][0] = divergenceCleaningPotential[i][j][k][0] + fakeDeltaT*cleanUpRightPart(i, j, k);
                        } else if(i < xnumber){
                            if(rank == nprocs - 1 && boundaryConditionType != PERIODIC) {
                                divergenceCleanUpRightPart[i][j][k][0] = 0;
                            } else {
                                divergenceCleanUpRightPart[i][j][k][0] = divergenceCleaningPotential[i][j][k][0] + fakeDeltaT*cleanUpRightPart(i, j, k);
                            }
                        } else {
                            divergenceCleanUpRightPart[i][j][k][0] = 0;
                        }
                    }
                }
            }

            generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart,
                                             divergenceCleaningPotential, xnumber, ynumber, znumber, 1, xnumberGeneral,
                                             ynumberGeneral, znumberGeneral, 1E-8, maxGMRESIterations, false, verbosity);
        }*/
        /*FILE* solutionFile1 = fopen((outputDir + "solution1.dat").c_str(), "w");
        for(int i = 0; i < xnumber + 1; ++i){
            fprintf(solutionFile1, "%15.10g\n", divergenceCleaningPotential[i][0][0][0]);
        }
        fclose(solutionFile1);*/
	} else {
		double ***rightPart = new double **[xnumber];
		for (int i = 0; i < xnumber; ++i) {
			rightPart[i] = new double *[ynumber];
			for (int j = 0; j < ynumber; ++j) {
				rightPart[i][j] = new double[znumber];
				for (int k = 0; k < znumber; ++k) {
					rightPart[i][j][k] = -cleanUpRightPart(i, j, k);
				}
			}
		}

		//Complex*** rightPartFourier = evaluateFourierTranslation(rightPart);
		Complex ***rightPartFourier = fastFourierTransition(rightPart, xnumber, ynumber, znumber);

		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					if (i == 0 && j == 0 && k == 0) {
						rightPartFourier[i][j][k] = Complex(0, 0);
					} else {
						double kx = i;
						/*if(i >= xnumber/2.0){
                            kx = xnumber - i;// - 2;
                        }*/
						double ky = j;
						/*if(j >= ynumber/2.0){
                            ky = ynumber - j;
                        }*/
						double kz = k;
						/*if(k >= znumber/2.0){
                            kz = znumber - k;
                        }*/
						rightPartFourier[i][j][k] = rightPartFourier[i][j][k] * 2 / (-4 * pi * pi * ((kx * kx / (xsize * xsize)) + (ky * ky / (ysize * ysize)) + (kz * kz / (zsize * zsize))));
						//rightPartFourier[i][j][k] = rightPartFourier[i][j][k]/(-4*pi*pi*((kx*kx*1.0/(xnumber*xnumber))  + (ky*ky/(ysize*ysize)) + (kz*kz/(zsize*zsize))));
						alertNaNOrInfinity(rightPartFourier[i][j][k].re, "divergence cleaning right part = NaN\n");
						alertNaNOrInfinity(rightPartFourier[i][j][k].im, "divergence cleaning right part = NaN\n");
					}
				}
			}
		}

		//double*** potential = evaluateReverceFourierTranslation(rightPartFourier);
		double ***potential = fastFourierReverceTransition(rightPartFourier, xnumber, ynumber, znumber);

		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					divergenceCleaningPotential[i][j][k][0] = potential[i][j][k];
				}
				delete[] potential[i][j];
				delete[] rightPartFourier[i][j];
				delete[] rightPart[i][j];
			}
			delete[] potential[i];
			delete[] rightPartFourier[i];
			delete[] rightPart[i];
		}
		delete[] potential;
		delete[] rightPartFourier;
		delete[] rightPart;
	}
	evaluateDivergenceCleaningField();
	substractMeanEfield();
	updateFieldByCleaning();

	if(timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("cleaning divergence time = %g sec\n", procTime/CLOCKS_PER_SEC);
	}
	//double div = evaluateDivCleaningE(1);

	//updateBoundariesNewField();
}

void Simulation::substractMeanDensity() {
	double meanDensity[1];
	meanDensity[0] = 0;
	for(int i = 1; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				meanDensity[0] += chargeDensity[i][j][k];
			}
		}
	}
	double tempDensity = meanDensity[0];

	meanDensity[0] = tempDensity/(xnumberGeneral*ynumberGeneral*znumberGeneral);

	for(int i = 0; i < xnumber+1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				chargeDensity[i][j][k] -= meanDensity[0];
			}
		}
	}
}

void Simulation::cleanupDivergence1d() {
	double rightField[3];
		rightField[0] = 0;
		rightField[1] = 0;
		rightField[2] = 0;
	

	for(int i = 0; i < 3; ++i){
		divergenceCleaningField[xnumber+1][0][0][i] = rightField[i];
	}

	for (int i = xnumber; i >= 0; --i) {
		double clean_up_right_part = cleanUpRightPart(i, 0, 0);
		divergenceCleaningField[i][0][0][0] = divergenceCleaningField[i + 1][0][0][0] - clean_up_right_part * deltaX;
		divergenceCleaningField[i][0][0][1] = 0;
		divergenceCleaningField[i][0][0][2] = 0;
        if((rank == nprocs - 1) && (boundaryConditionType != PERIODIC)){
            if(i >= xnumber - additionalBinNumber){
                divergenceCleaningField[i][0][0][0]= 0;
            }
        }
	}

	double leftField[3];

	for(int i = 0; i < 3; ++i){
		leftField[i] = divergenceCleaningField[2][0][0][i];
	}


	//substractMeanEfield();

	updateFieldByCleaning();

	if(rank == 0 && boundaryConditionType == SUPER_CONDUCTOR_LEFT){
		newEfield[1][0][0].y = 0;
		newEfield[1][0][0].z = 0;
	}

	if(rank == nprocs - 1 && boundaryConditionType != PERIODIC){
		newEfield[xnumber][0][0] = E0;
		newEfield[xnumber+1][0][0] = E0;
	}

}

void Simulation::substractMeanEfield() {
	Vector3d meanField = Vector3d(0, 0, 0);
	double meanFieldSquare = 0;

	/*int maxI = xnumber;
	if(rank == nprocs-1 && boundaryConditionType != PERIODIC){
		maxI = xnumber+1;
	}*/

	for(int j = 0; j < ynumber; ++j) {
		for (int k = 0; k < znumber; ++k) {
			for (int i = 1; i < xnumber; ++i) {
				meanField.x += divergenceCleaningField[i][j][k][0];
				meanField.y += divergenceCleaningField[i][j][k][1];
				meanField.z += divergenceCleaningField[i][j][k][2];
				meanFieldSquare += divergenceCleaningField[i][j][k][0]*divergenceCleaningField[i][j][k][0] + 
					divergenceCleaningField[i][j][k][1]*divergenceCleaningField[i][j][k][1] + 
					divergenceCleaningField[i][j][k][2]*divergenceCleaningField[i][j][k][2];
			}
		}
	}

    int Nx = xnumberGeneral+1;

    if(boundaryConditionType != PERIODIC){
        Nx = Nx - additionalBinNumber;
    }

	meanField.x /= Nx*ynumber*znumber;
	meanField.y /= Nx*ynumber*znumber;
	meanField.z /= Nx*ynumber*znumber;
	meanFieldSquare /= Nx*ynumber*znumber;

	for (int i = 0; i < xnumber+1; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
                if(i < xnumber - additionalBinNumber || rank < nprocs - 1 || boundaryConditionType == PERIODIC) {
                    newEfield[i][j][k] -= meanField;
                }
			}
		}
	}
}

void Simulation::updateFieldByCleaning() {

	/*if((ynumber == 1) && (znumber == 1)){
		divergenceCleaningField[xnumber][0][0][0] = 0;
		divergenceCleaningField[xnumber][0][0][1] = 0;
		divergenceCleaningField[xnumber][0][0][2] = 0;

		for (int i = xnumber - 1; i >= 0; --i) {
			divergenceCleaningField[i][0][0][0] = divergenceCleaningField[i + 1][0][0][0] - cleanUpRightPart(i, 0, 0) * deltaX;
			divergenceCleaningField[i][0][0][1] = 0;
			divergenceCleaningField[i][0][0][2] = 0;
		}
	}*/
	for (int i = 0; i < xnumber+2; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				newEfield[i][j][k].x += divergenceCleaningField[i][j][k][0];
				newEfield[i][j][k].y += divergenceCleaningField[i][j][k][1];
				newEfield[i][j][k].z += divergenceCleaningField[i][j][k][2];
			}
		}
	}


	for(int i = 0; i < xnumber+2; ++i) {
		for(int k = 0; k < znumber + 1; ++k) {
			newEfield[i][ynumber][k] = newEfield[i][0][k];
		}
	}

	for(int i = 0; i < xnumber+2; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			newEfield[i][j][znumber] = newEfield[i][j][0];
		}
	}


		if (boundaryConditionType == PERIODIC) {
			for (int j = 0; j < ynumber + 1; ++j) {
				for (int k = 0; k < znumber + 1; ++k) {
					newEfield[xnumber-1][j][k] = newEfield[0][j][k];
					newEfield[xnumber][j][k] = newEfield[1][j][k];
					newEfield[xnumber+1][j][k] = newEfield[2][j][k];
				}
			}
		}
}

void Simulation::evaluateDivergenceCleaningField() {
	for (int j = 0; j < ynumber; ++j) {
		for (int k = 0; k < znumber; ++k) {
			divergenceCleaningField[0][j][k][0] = 0;
			divergenceCleaningField[0][j][k][1] = 0;
			divergenceCleaningField[0][j][k][2] = 0;
		}
	}
	for (int i = 1; i <= xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				int prevI = i - 1;
				if (prevI < 0) {
					//if(boundaryConditionType == PERIODIC){
						prevI = xnumber - 1;
					//} else {
						//prevI = 0;
					//}
				}

				int prevJ = j - 1;
				if (prevJ < 0) {
					prevJ = ynumber - 1;
				}

				int prevK = k - 1;
				if (prevK < 0) {
					prevK = znumber - 1;
				}

				divergenceCleaningField[i][j][k][0] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i][prevJ][k][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[i][prevJ][prevK][0]
					- divergenceCleaningPotential[prevI][j][k][0] - divergenceCleaningPotential[prevI][prevJ][k][0] - divergenceCleaningPotential[prevI][j][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0]) / (4 * deltaX);

				divergenceCleaningField[i][j][k][1] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[prevI][j][k][0] + divergenceCleaningPotential[i][j][prevK][0] + divergenceCleaningPotential[prevI][j][prevK][0]
					- divergenceCleaningPotential[i][prevJ][k][0] - divergenceCleaningPotential[prevI][prevJ][k][0] - divergenceCleaningPotential[i][prevJ][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0]) / (4 * deltaY);

				divergenceCleaningField[i][j][k][2] = -(divergenceCleaningPotential[i][j][k][0] + divergenceCleaningPotential[i][prevJ][k][0] + divergenceCleaningPotential[prevI][j][k][0] + divergenceCleaningPotential[prevI][prevJ][k][0]
					- divergenceCleaningPotential[i][j][prevK][0] - divergenceCleaningPotential[i][prevJ][prevK][0] - divergenceCleaningPotential[prevI][j][prevK][0] - divergenceCleaningPotential[prevI][prevJ][prevK][0]) / (4 * deltaZ);
			}
		}
	}
	for (int j = 0; j < ynumber; ++j) {
		for (int k = 0; k < znumber; ++k) {
			divergenceCleaningField[xnumber + 1][j][k][0] = 0;
			divergenceCleaningField[xnumber + 1][j][k][1] = 0;
			divergenceCleaningField[xnumber + 1][j][k][2] = 0;
		}
	}

	/*Vector3d constantField;
	constantField.x = divergenceCleaningField[xnumber - 1][0][0][0];
	constantField.y = divergenceCleaningField[xnumber - 1][0][0][1];
	constantField.z = divergenceCleaningField[xnumber - 1][0][0][2];

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				divergenceCleaningField[i][j][k][0] -= constantField.x;
				divergenceCleaningField[i][j][k][1] -= constantField.y;
				divergenceCleaningField[i][j][k][2] -= constantField.z;
			}
		}
	}*/
}

void Simulation::createDivergenceCleanupInternalEquation(int i, int j, int k) {

	/*if (i == 0 && j == 0 && k == 0) {
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
		divergenceCleanUpRightPart[i][j][k][0] = 0;
		return;
	}

	if (i == xnumber - 1 && j == 0 && k == 0) {
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(-1, i-1, j, k, 0));
		divergenceCleanUpRightPart[i][j][k][0] = 0;
		return;
	}*/

	/*if (i == xnumber - 1 && j == ynumber - 1 && k == znumber - 1) {
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
		divergenceCleanUpRightPart[i][j][k][0] = 0;
		return;
	}*/

	int prevI = i - 1;
	if (prevI < 0) {
		prevI = xnumber - 1;
	}
	int nextI = i + 1;
	if (nextI > xnumber) {
		nextI = 0;
	}

	int prevJ = j - 1;
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int nextJ = j + 1;
	if (nextJ >= ynumber) {
		nextJ = 0;
	}

	int prevK = k - 1;
	if (prevK < 0) {
		prevK = ynumber - 1;
	}
	int nextK = k + 1;
	if (nextK >= ynumber) {
		nextK = 0;
	}


	divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k);
	//divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k)*deltaX2;

	if((ynumber > 1) && (znumber > 1)){
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
	} else if(ynumber > 1){
		double element = -1 / (deltaX2) - 1 / (deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / (4 * deltaX2) + 1 / (4 * deltaY2) ;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 0));

		element = 1 / (2 * deltaX2) - 1 / (2 * deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));

		element = - 1 / (2 * deltaX2) + 1 / (2 * deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

	} else if(znumber > 1){
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

		element = 1/deltaX2;
		//element = 1;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));
	}
}

void Simulation::createDivergenceCleanupLeftFakeEquation(int i, int j, int k) {
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
	divergenceCleanUpRightPart[i][j][k][0] = 0;
}

void Simulation::createDivergenceCleanupRightFakeEquation(int i, int j, int k) {
	divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1, i, j, k, 0));
	divergenceCleanUpRightPart[i][j][k][0] = 0;
}

void Simulation::createDivergenceCleanupSuperConductorEquation(int i, int j, int k) {

	int nextI = i + 1;
	if (nextI > xnumber) {
		nextI = 0;
	}

	int prevJ = j - 1;
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int nextJ = j + 1;
	if (nextJ >= ynumber) {
		nextJ = 0;
	}

	int prevK = k - 1;
	if (prevK < 0) {
		prevK = ynumber - 1;
	}
	int nextK = k + 1;
	if (nextK >= ynumber) {
		nextK = 0;
	}


	divergenceCleanUpRightPart[i][j][k][0] = -cleanUpRightPart(i, j, k);

	if((ynumber > 1) && (znumber > 1)){
		double element = -3 / (4 * deltaX2) - 1 / (4 * deltaY2) - 1 / (4 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / (16 * deltaX2) + 1 / (16 * deltaY2) + 1 / (16 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, prevK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, prevK, 0));

		element = 1 / (8 * deltaX2) + 1 / (8 * deltaY2) - 1 / (8 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 0));

		element = 1 / (8 * deltaX2) - 1 / (8 * deltaY2) + 1 / (8 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 0));

		element = - 3 / (16 * deltaX2) + 1 / (16 * deltaY2) + 1 / (16 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, prevK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, prevK, 0));

		element = 1 / (4 * deltaX2) - 1 / (4 * deltaY2) - 1 / (4 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));

		element = - 3 / (8 * deltaX2) + 1 / (8 * deltaY2) - 1 / (8 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

		element = - 3 / (8 * deltaX2) - 1 / (8 * deltaY2) + 1 / (8 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
	} else if(ynumber > 1){
		double element = -1.5 / (deltaX2) - 0.5 / (deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1 / (4 * deltaX2) + 1 / (4 * deltaY2) ;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 0));

		element = 1 / (2 * deltaX2) - 1 / (2 * deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));

		element = - 3 / (4 * deltaX2) + 1 / (4 * deltaY2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));

	} else if(znumber > 1){
		double element = -1.5 / (deltaX2) - 0.5 / (deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));


		element = 1 / (4 * deltaX2) + 1 / (4 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 0));


		element = 1 / (2 * deltaX2) - 1 / (2 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));

		element = - 3 / (4 * deltaX2) + 1 / (4 * deltaZ2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
	} else {
		double element = -3 / (deltaX2);
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

		element = 1/deltaX2;
		//element = 1;
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
	}
}

double Simulation::cleanUpRightPart(int i, int j, int k) {
	double div = evaluateDivNewE(i, j, k);

	return 4 * pi * chargeDensity[i][j][k] - div;
}

Complex*** Simulation::evaluateFourierTranslation(double*** a){
	Complex*** result = new Complex**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		result[i] = new Complex*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			result[i][j] = new Complex[znumber];
			for(int k = 0; k < znumber; ++k){
				result[i][j][k] = Complex(0, 0);

				for(int tempi = 0; tempi < xnumber; ++tempi){
					for(int tempj = 0; tempj < ynumber; ++tempj){
						for(int tempk = 0; tempk < znumber; ++tempk){
							result[i][j][k] += complexExp(-2*pi*((i*tempi*1.0/xnumber) + (j*tempj*1.0/ynumber) + (k*tempk*1.0/znumber)))*a[tempi][tempj][tempk];
						}
					}
				}

				result[i][j][k] = result[i][j][k]/(xnumber*ynumber*znumber);
			}
		}
	}

	return result;
}

double*** Simulation::evaluateReverceFourierTranslation(Complex*** a){
	double*** result = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		result[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			result[i][j] = new double[znumber];
			for(int k = 0; k < znumber; ++k){
				result[i][j][k] = 0;

				for(int tempi = 0; tempi < xnumber; ++tempi){
					for(int tempj = 0; tempj < ynumber; ++tempj){
						for(int tempk = 0; tempk < znumber; ++tempk){
							double kx = 2*pi*tempi;
							if(tempi >= xnumber/2.0){
								kx = -2*pi*(tempi - xnumber + 2);
							}

							double ky = 2*pi*tempj;
							if(tempj >= ynumber/2.0){
								ky = -2*pi*(tempj - ynumber + 2);
							}

							double kz = 2*pi*tempk;
							if(tempk >= znumber/2.0){
								kz = -2*pi*(tempk - znumber + 2);
							}
							result[i][j][k] += (complexExp(((i*kx/xnumber) + (j*ky/ynumber) + (k*kz/znumber)))*a[tempi][tempj][tempk]).re;
						}
					}
				}
			}
		}
	}

	return result;
}