#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <mpi.h>
#include <time.h>

#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "simulation.h"
#include "complex.h"
#include "fourier.h"
#include "util.h"

void Simulation::filterFields(int cutWaveNumber){
	MPI_Barrier(cartComm);
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) printf("filtering fields\n");
	if ((rank == 0) && (verbosity > 0)) printLog("filtering fields\n");

	if(boundaryConditionType != PERIODIC){
		if(cartCoord[0] == 0){
			if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
				for(int i = 0; i <= 1 + additionalBinNumber; ++i){
					for(int j = 0; j < ynumberAdded + 1; ++j){
						for(int k = 0; k < znumberAdded + 1; ++k){
							newEfield[i][j][k] = newEfield[i][j][k] + E0;
						}
					}
				}
			}
		}
	}

	filterFieldGeneral(newEfield, cutWaveNumber);
	//filterFieldGeneral(newBfield, cutWaveNumber);

	if(boundaryConditionType != PERIODIC){
		if(cartCoord[0] == cartDim[0] - 1){
			for(int i = xnumberAdded - 1 - additionalBinNumber; i < xnumberAdded + 1; ++i){
				for(int j = 0; j < ynumberAdded + 1; ++j){
					for(int k = 0; k < znumberAdded + 1; ++k){
						newEfield[i][j][k] = E0;
					}
				}
			}

			/*for(int i = xnumberAdded - 1 - additionalBinNumber; i < xnumberAdded; ++i){
				for(int j = 0; j < ynumberAdded; ++j){
					for(int k = 0; k < znumberAdded; ++k){
						newBfield[i][j][k] = B0;
					}
				}
			}*/
		}
		if(cartCoord[0] == 0){
			if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
				for(int i = 0; i <= 1 + additionalBinNumber; ++i){
					for(int j = 0; j < ynumberAdded + 1; ++j){
						for(int k = 0; k < znumberAdded + 1; ++k){
							newEfield[i][j][k] = newEfield[i][j][k] - E0;
							newEfield[i][j][k].y = 0;
							newEfield[i][j][k].z = 0;
						}
					}
				}

				/*for (int i = 0; i <= additionalBinNumber; ++i) {
					for (int j = 0; j < ynumberAdded; ++j) {
						for (int k = 0; k < znumberAdded; ++k) {
							newBfield[i][j][k] = newBfield[additionalBinNumber + 1][j][k];
						}
					}
				}*/
			} else if(boundaryConditionType == FREE_BOTH){
				for(int i = 0; i < 1 + additionalBinNumber; ++i){
					for(int j = 0; j < ynumberAdded + 1; ++j){
						for(int k = 0; k < znumberAdded + 1; ++k){
							newEfield[i][j][k] = E0;
						}
					}
				}

				/*for (int i = 0; i <= additionalBinNumber; ++i) {
					for (int j = 0; j < ynumberAdded; ++j) {
						for (int k = 0; k < znumberAdded; ++k) {
							newBfield[i][j][k] = B0;
						}
					}
				}*/
			}
		}
	}

	exchangeGeneralEfield(newEfield);
	exchangeGeneralBfield(newBfield);

	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("filtering time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::filterFieldGeneral(Vector3d*** field, int cutWaveNumber){
	int* xabsoluteIndex = new int[xnumberAdded + 1];
	int* yabsoluteIndex = new int[ynumberAdded + 1];
	int* zabsoluteIndex = new int[znumberAdded + 1];

	for(int i = 0; i < xnumberAdded + 1; ++i){
		xabsoluteIndex[i] = firstAbsoluteXindex + i;
	}

	for(int j = 0; j < ynumberAdded + 1; ++j){
		yabsoluteIndex[j] = firstAbsoluteYindex + j;
	}
	
	for(int k = 0; k < znumberAdded + 1; ++k){
		zabsoluteIndex[k] = firstAbsoluteZindex + k;
	}

	for(int i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				fourierScalarInput[i][j][k] = Complex(field[i][j][k].x,0);
				fourierScalarOutput[i][j][k] = Complex(0, 0);
			}
		}
	}

	fourierTranslation(fourierScalarInput, fourierScalarOutput, fourierScalarTempOutput, fourierScalarTempOutput1, true, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY, cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);

	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				double kx = (i + firstAbsoluteXindex)*2*pi/xsizeGeneral;
				double ky = (j + firstAbsoluteYindex)*2*pi/ysizeGeneral;
				double kz = (k + firstAbsoluteZindex)*2*pi/zsizeGeneral;
				double kw = sqrt(kx*kx + ky*ky + kz*kz);
				if(kw > 2*pi/(cutWaveNumber*deltaX)){
					fourierScalarOutput[i][j][k] = Complex(0, 0);
				}
			}
		}
	}

	fourierTranslation(fourierScalarOutput, fourierScalarInput, fourierScalarTempOutput, fourierScalarTempOutput1, false, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY, cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				field[i][j][k].x = fourierScalarInput[i][j][k].re;
			}
		}
	}

	for(int i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				fourierScalarInput[i][j][k] = Complex(field[i][j][k].y,0);
				fourierScalarOutput[i][j][k] = Complex(0, 0);
			}
		}
	}

	fourierTranslation(fourierScalarInput, fourierScalarOutput, fourierScalarTempOutput, fourierScalarTempOutput1, true, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY, cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);

	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				double kx = (i + firstAbsoluteXindex)*2*pi/xsizeGeneral;
				double ky = (j + firstAbsoluteYindex)*2*pi/ysizeGeneral;
				double kz = (k + firstAbsoluteZindex)*2*pi/zsizeGeneral;
				double kw = sqrt(kx*kx + ky*ky + kz*kz);
				if(kw > 2*pi/(cutWaveNumber*deltaX)){
					fourierScalarOutput[i][j][k] = Complex(0, 0);
				}
			}
		}
	}


	fourierTranslation(fourierScalarOutput, fourierScalarInput, fourierScalarTempOutput, fourierScalarTempOutput1, false, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY, cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				field[i][j][k].y = fourierScalarInput[i][j][k].re;
			}
		}
	}


	for(int i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				fourierScalarInput[i][j][k] = Complex(field[i][j][k].z,0);
				fourierScalarOutput[i][j][k] = Complex(0, 0);
			}
		}
	}

	fourierTranslation(fourierScalarInput, fourierScalarOutput, fourierScalarTempOutput, fourierScalarTempOutput1, true, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY, cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);

	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				double kx = (i + firstAbsoluteXindex)*2*pi/xsizeGeneral;
				double ky = (j + firstAbsoluteYindex)*2*pi/ysizeGeneral;
				double kz = (k + firstAbsoluteZindex)*2*pi/zsizeGeneral;
				double kw = sqrt(kx*kx + ky*ky + kz*kz);
				if(kw > 2*pi/(cutWaveNumber*deltaX)){
					fourierScalarOutput[i][j][k] = Complex(0, 0);
				}
			}
		}
	}

	fourierTranslation(fourierScalarOutput, fourierScalarInput, fourierScalarTempOutput, fourierScalarTempOutput1, false, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY, cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				field[i][j][k].z = fourierScalarInput[i][j][k].re;
			}
		}
	}

	delete[] xabsoluteIndex;
	delete[] yabsoluteIndex;
	delete[] zabsoluteIndex;

}

///////local/////////

void Simulation::filterFieldsLocal(int cutWaveNumber){
	MPI_Barrier(cartComm);
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) printf("filtering fields\n");
	if ((rank == 0) && (verbosity > 0)) printLog("filtering fields\n");

	if ((rank == 0) && (verbosity > 0)) printf("filtering  E field\n");
	filterFieldGeneralLocal(newEfield, cutWaveNumber);
	MPI_Barrier(cartComm);
	if ((rank == 0) && (verbosity > 0)) printf("filtering B field\n");
	filterFieldGeneralLocal(newBfield, cutWaveNumber);

	exchangeGeneralEfield(newEfield);
	exchangeGeneralBfield(newBfield);

	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("filtering time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::filterFieldGeneralLocal(Vector3d*** field, int cutWaveNumber){
	if ((rank == 0) && (verbosity > 1)) printf("fourier field general\n");
	MPI_Barrier(cartComm);
	if ((rank == 0) && (verbosity > 1)) printf("xnumberAdded = %d ynumberAdded = %d znumberAdded = %d\n", xnumberAdded, ynumberAdded, znumberAdded);
	for(int i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				if ((rank == 0) && (verbosity > 1)) printf("i = %d j = %d k = %d\n", i, j, k);
				if ((rank == 0) && (verbosity > 1)) printf("field = %g\n", field[i][j][k].x);
				if ((rank == 0) && (verbosity > 1)) printf("fourierInput = %g\n", fourierScalarInput[i][j][k].re);
				if ((rank == 0) && (verbosity > 1)) printf("fourierOutput = %g\n", fourierScalarOutput[i][j][k].re);
				fourierScalarInput[i][j][k] = Complex(field[i][j][k].x,0);
				fourierScalarOutput[i][j][k] = Complex(0, 0);
			}
		}
	}
	if ((rank == 0) && (verbosity > 1)) printf("fourier field x\n");
	MPI_Barrier(cartComm);

	fourierTranslationLocal(fourierScalarInput, fourierScalarOutput, fourierScalarTempOutput, fourierScalarTempOutput1, true, xnumberAdded, ynumberAdded, znumberAdded);

	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				double kx = (i - 1 - additionalBinNumber)*2*pi/xsize;
				double ky = (j - 1 - additionalBinNumber)*2*pi/ysize;
				double kz = (k - 1 - additionalBinNumber)*2*pi/zsize;
				double kw = sqrt(kx*kx + ky*ky + kz*kz);
				if(kw > 2*pi/(cutWaveNumber*deltaX)){
					fourierScalarOutput[i][j][k] = Complex(0, 0);
				}
			}
		}
	}

	if ((rank == 0) && (verbosity > 1)) printf("inverse fourier field x\n");

	fourierTranslationLocal(fourierScalarOutput, fourierScalarInput, fourierScalarTempOutput, fourierScalarTempOutput1, false, xnumberAdded, ynumberAdded, znumberAdded);
	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				field[i][j][k].x = fourierScalarInput[i][j][k].re;
			}
		}
	}

	for(int i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				fourierScalarInput[i][j][k] = Complex(field[i][j][k].y,0);
				fourierScalarOutput[i][j][k] = Complex(0, 0);
			}
		}
	}

	if ((rank == 0) && (verbosity > 1)) printf("fourier field y\n");

	fourierTranslationLocal(fourierScalarInput, fourierScalarOutput, fourierScalarTempOutput, fourierScalarTempOutput1, true, xnumberAdded, ynumberAdded, znumberAdded);

	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				double kx = (i - 1 - additionalBinNumber)*2*pi/xsize;
				double ky = (j - 1 - additionalBinNumber)*2*pi/ysize;
				double kz = (k - 1 - additionalBinNumber)*2*pi/zsize;
				double kw = sqrt(kx*kx + ky*ky + kz*kz);
				if(kw > 2*pi/(cutWaveNumber*deltaX)){
					fourierScalarOutput[i][j][k] = Complex(0, 0);
				}
			}
		}
	}

	if ((rank == 0) && (verbosity > 1)) printf("inverse fourier field y\n");

	fourierTranslationLocal(fourierScalarOutput, fourierScalarInput, fourierScalarTempOutput, fourierScalarTempOutput1, false, xnumberAdded, ynumberAdded, znumberAdded);
	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				field[i][j][k].y = fourierScalarInput[i][j][k].re;
			}
		}
	}


	for(int i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				fourierScalarInput[i][j][k] = Complex(field[i][j][k].z,0);
				fourierScalarOutput[i][j][k] = Complex(0, 0);
			}
		}
	}

	if ((rank == 0) && (verbosity > 1)) printf("fourier field z\n");

	fourierTranslationLocal(fourierScalarInput, fourierScalarOutput, fourierScalarTempOutput, fourierScalarTempOutput1, true, xnumberAdded, ynumberAdded, znumberAdded);

	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				double kx = (i - 1 - additionalBinNumber)*2*pi/xsize;
				double ky = (j - 1 - additionalBinNumber)*2*pi/ysize;
				double kz = (k - 1 - additionalBinNumber)*2*pi/zsize;
				double kw = sqrt(kx*kx + ky*ky + kz*kz);
				if(kw > 2*pi/(cutWaveNumber*deltaX)){
					fourierScalarOutput[i][j][k] = Complex(0, 0);
				}
			}
		}
	}

	if ((rank == 0) && (verbosity > 1)) printf("inverse fourier field z\n");

	fourierTranslationLocal(fourierScalarOutput, fourierScalarInput, fourierScalarTempOutput, fourierScalarTempOutput1, false, xnumberAdded, ynumberAdded, znumberAdded);
	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				field[i][j][k].z = fourierScalarInput[i][j][k].re;
			}
		}
	}
}

void Simulation::updateMaxEderivativePoint(){
	double maxDer = 0;
	int maxDerPoint = 1+additionalBinNumber;
	for(int i = 1 + additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i){
		Vector3d curE = averageEfieldYZ(i);
		Vector3d nextE = averageEfieldYZ(i+1);

		//todo sign!
		double der = (curE - nextE).norm()/deltaX;
		if(der > maxDer){
			maxDer = der;
			maxDerPoint = i;
		}
	}

	maxDerPoint += firstAbsoluteXindex;
}

void Simulation::substractStep(Vector3d*** field, Vector3d left, Vector3d right, int sign){
	for(int i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				if((i + firstAbsoluteXindex) <= derExPoint){
					field[i][j][k] = field[i][j][k] - left*sign;
				} else {
					field[i][j][k] = field[i][j][k] - right*sign;
				}
			}
		}
	}
}