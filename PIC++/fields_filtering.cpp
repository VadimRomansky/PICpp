#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <mpi.h>

#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "simulation.h"
#include "complex.h"
#include "fourier.h"

void Simulation::filterFields(int cutWaveNumber){
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) printf("filtering fields\n");
	if ((rank == 0) && (verbosity > 0)) printLog("filtering fields\n");
	filterFieldGeneral(newEfield, cutWaveNumber);
	filterFieldGeneral(newBfield, cutWaveNumber);

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

	fourierTranslation(fourierScalarInput, fourierScalarOutput, fourierScalarTempOutput, fourierScalarTempOutput1, true, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);

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

	fourierTranslation(fourierScalarOutput, fourierScalarInput, fourierScalarTempOutput, fourierScalarTempOutput1, false, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
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

	fourierTranslation(fourierScalarInput, fourierScalarOutput, fourierScalarTempOutput, fourierScalarTempOutput1, true, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);

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


	fourierTranslation(fourierScalarOutput, fourierScalarInput, fourierScalarTempOutput, fourierScalarTempOutput1, false, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
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

	fourierTranslation(fourierScalarInput, fourierScalarOutput, fourierScalarTempOutput, fourierScalarTempOutput1, true, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);

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

	fourierTranslation(fourierScalarOutput, fourierScalarInput, fourierScalarTempOutput, fourierScalarTempOutput1, false, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
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
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((rank == 0) && (verbosity > 0)) printf("filtering fields\n");
	if ((rank == 0) && (verbosity > 0)) printLog("filtering fields\n");
	filterFieldGeneral(newEfield, cutWaveNumber);
	filterFieldGeneralLocal(newEfield, cutWaveNumber);
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
	for(int i = 0; i < xnumberAdded; ++i){
		for(int j = 0; j < ynumberAdded; ++j){
			for(int k = 0; k < znumberAdded; ++k){
				fourierScalarInput[i][j][k] = Complex(field[i][j][k].x,0);
				fourierScalarOutput[i][j][k] = Complex(0, 0);
			}
		}
	}

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

	fourierTranslationLocal(fourierScalarOutput, fourierScalarInput, fourierScalarTempOutput, fourierScalarTempOutput1, false, xnumberAdded, ynumberAdded, znumberAdded);
	for(int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k){
				field[i][j][k].z = fourierScalarInput[i][j][k].re;
			}
		}
	}
}