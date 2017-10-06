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
	filterElectricFieldGeneral(newEfield, cutWaveNumber);
	filterMagneticFieldGeneral(newBfield, cutWaveNumber);
}

void Simulation::filterElectricFieldGeneral(Vector3d*** field, int cutWaveNumber){
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

	Complex*** tempE = new Complex**[xnumberAdded + 1];
	Complex*** fourierE = new Complex**[xnumberAdded + 1];
	for(int i = 0; i < xnumberAdded + 1; ++i){
		tempE[i] = new Complex*[ynumberAdded + 1];
		fourierE[i] = new Complex*[ynumberAdded + 1];
		for(int j = 0; j < ynumberAdded + 1; ++j){
			tempE[i][j] = new Complex[znumberAdded + 1];
			fourierE[i][j] = new Complex[znumberAdded + 1];
			for(int k = 0; k < znumberAdded + 1; ++k){
				tempE[i][j][k] = Complex(field[i][j][k].x,0);
				fourierE[i][j][k] = Complex(0, 0);
			}
		}
	}

	fourierTranslation(tempE, fourierE, true, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);

	for(int i = 1 + additionalBinNumber; i < xnumberAdded + 1; ++i){
		for(int j = 1 + additionalBinNumber; j < ynumberAdded + 1; ++j){
			for(int k = 1 + additionalBinNumber; k < znumberAdded + 1; ++k){
				double kx = (i - 1 - additionalBinNumber + firstAbsoluteXindex)*2*pi/xsizeGeneral;
				double ky = (i - 1 - additionalBinNumber + firstAbsoluteYindex)*2*pi/ysizeGeneral;
				double kz = (i - 1 - additionalBinNumber + firstAbsoluteZindex)*2*pi/zsizeGeneral;
				double kw = sqrt(kx*kx + ky*ky + kz*kz);
				if(kw > 2*pi/(cutWaveNumber*deltaX)){
					fourierE[i][j][k] = Complex(0, 0);
				}
			}
		}
	}


	for(int i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				tempE[i][j][k] = Complex(field[i][j][k].y,0);
				fourierE[i][j][k] = Complex(0, 0);
			}
		}
	}


	for(int i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			for(int k = 0; k < znumberAdded + 1; ++k){
				tempE[i][j][k] = Complex(field[i][j][k].z,0);
				fourierE[i][j][k] = Complex(0, 0);
			}
		}
	}

	delete[] xabsoluteIndex;
	delete[] yabsoluteIndex;
	delete[] zabsoluteIndex;

	for(int i = 0; i < xnumberAdded + 1; ++i){
		for(int j = 0; j < ynumberAdded + 1; ++j){
			delete[] tempE[i][j];
			delete[] fourierE[i][j];
		}
		delete[] tempE[i];
		delete[] fourierE[i];
	}
	delete[] tempE;
	delete[] fourierE;
}

void Simulation::filterMagneticFieldGeneral(Vector3d*** field, int cutWaveNumber){
}