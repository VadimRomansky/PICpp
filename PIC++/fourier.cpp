#include "stdio.h"
#include <stdlib.h>
#include <string>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "fourier.h"
#include "complex.h"
#include "constants.h"
#include "output.h"
#include "paths.h"

Complex fourierTranslationXoneHarmonic(Complex*** input, bool direct, int xnumberAdded, int j, int k, int xnumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* xabsoluteIndex, Complex* localFactor, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	int maxI = xnumberAdded;
	for (int i = 1+additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i) {
		sum += input[i][j][k] * localFactor[i];
	}
	double out[2];
	double in[2];
	out[0] = sum.re;
	out[1] = sum.im;
	MPI_Allreduce(out, in, 2, MPI_DOUBLE, MPI_SUM, reducedCartComm);
	sum = Complex(in[0], in[1]);
	if (!direct) {
		sum = sum / xnumberGeneral;
	}
	return sum;
}

void fourierTranslationX(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int* xabsoluteIndex, MPI_Comm cartCommX, int* cartCoord, int* cartDim) {
	if (xnumberGeneral == 1) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[1+additionalBinNumber][j][k] = input[1+additionalBinNumber][j][k];
			}
		}
	} else {
		for (int knumber = 0; knumber < xnumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor*2*pi*knumber)/xnumberGeneral);
			Complex* localFactor = new Complex[xnumberAdded];
			localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * xabsoluteIndex[1+additionalBinNumber] * knumber) / xnumberGeneral);
			for(int i = 2 + additionalBinNumber; i < xnumberAdded; ++i){
				localFactor[i] = localFactor[i-1]*factor;
			}
			for (int j = 1+additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
				for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
					Complex f = fourierTranslationXoneHarmonic(input, direct, xnumberAdded, j, k, xnumberGeneral, cartCommX, knumber, xabsoluteIndex, localFactor, cartCoord, cartDim);
					if(knumber >= xabsoluteIndex[1+additionalBinNumber] && knumber < xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber]){
						output[knumber - xabsoluteIndex[0]][j][k] = f;
					}
				}
			}
			delete[] localFactor;
		}
	}
}

Complex fourierTranslationYoneHarmonic(Complex*** input, bool direct, int ynumberAdded, int i, int k, int ynumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* yabsoluteIndex, Complex* localFactor, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	for (int j = 1+additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
		sum += input[i][j][k] * localFactor[j];
	}
	double out[2];
	double in[2];
	out[0] = sum.re;
	out[1] = sum.im;
	MPI_Allreduce(out, in, 2, MPI_DOUBLE, MPI_SUM, reducedCartComm);
	sum = Complex(in[0], in[1]);
	if (!direct) {
		sum = sum / ynumberGeneral;
	}
	return sum;
}

void fourierTranslationY(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int ynumberGeneral, int* yabsoluteIndex, MPI_Comm cartCommY, int* cartCoord, int* cartDim) {
	if (ynumberGeneral == 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[i][1+additionalBinNumber][k] = input[i][1+additionalBinNumber][k];
			}
		}
	} else {
		for (int knumber = 0; knumber < ynumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor*2*pi*knumber)/ynumberGeneral);
			Complex* localFactor = new Complex[ynumberAdded];
			localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * yabsoluteIndex[1+additionalBinNumber] * knumber) / ynumberGeneral);
			for(int j = 2 + additionalBinNumber; j < ynumberAdded; ++j){
				localFactor[j] = localFactor[j-1]*factor;
			}
			for (int i = 1+additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
				for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
					Complex f = fourierTranslationYoneHarmonic(input, direct, ynumberAdded, i, k, ynumberGeneral, cartCommY, knumber, yabsoluteIndex, localFactor, cartCoord, cartDim);
					if(knumber >= yabsoluteIndex[1+additionalBinNumber] && knumber < yabsoluteIndex[ynumberAdded - 1 - additionalBinNumber]){
						output[i][knumber - yabsoluteIndex[0]][k] = f;
					}
				}
			}
			delete[] localFactor;
		}
	}
}

Complex fourierTranslationZoneHarmonic(Complex*** input, bool direct, int znumberAdded, int i, int j, int znumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* zabsoluteIndex, Complex* localFactor, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	for (int k = 1+additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
		sum += input[i][j][k] * localFactor[k];
	}
	double out[2];
	double in[2];
	out[0] = sum.re;
	out[1] = sum.im;
	MPI_Allreduce(out, in, 2, MPI_DOUBLE, MPI_SUM, reducedCartComm);
	sum = Complex(in[0], in[1]);
	if (!direct) {
		sum = sum / znumberGeneral;
	}
	return sum;
}

void fourierTranslationZ(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int znumberGeneral, int* zabsoluteIndex, MPI_Comm cartCommZ, int* cartCoord, int* cartDim) {
	if (znumberGeneral == 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				output[i][j][1+additionalBinNumber] = input[i][j][1+additionalBinNumber];
			}
		}
	} else {
		for (int knumber = 0; knumber < znumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor*2*pi*knumber)/znumberGeneral);
			Complex* localFactor = new Complex[znumberAdded];
			localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * zabsoluteIndex[1+additionalBinNumber] * knumber) / znumberGeneral);
			for(int k = 2 + additionalBinNumber; k < znumberAdded; ++k){
				localFactor[k] = localFactor[k-1]*factor;
			}
			for (int i = 1+additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
				for (int j = 1+additionalBinNumber; j < ynumberAdded -1 - additionalBinNumber; ++j) {
					Complex f = fourierTranslationZoneHarmonic(input, direct, znumberAdded, i, j, znumberGeneral, cartCommZ, knumber, zabsoluteIndex, localFactor, cartCoord, cartDim);
					if(knumber >= zabsoluteIndex[1+additionalBinNumber] && knumber < zabsoluteIndex[znumberAdded - 1 - additionalBinNumber]){
						output[i][j][knumber - zabsoluteIndex[0]] = f;
					}
				}
			}
		}
	}
}

void fourierTranslation(Complex*** input, Complex*** output, Complex*** tempResultX, Complex*** tempResultXY, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, MPI_Comm& cartComm, MPI_Comm& cartCommX, MPI_Comm& cartCommY, MPI_Comm& cartCommZ , int* xabsoluteIndex, int* yabsoluteIndex, int* zabsoluteIndex, int* cartCoord, int* cartDim) {

	fourierTranslationX(input, tempResultX, direct, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, xabsoluteIndex, cartCommX, cartCoord, cartDim);

	fourierTranslationY(tempResultX, tempResultXY, direct, xnumberAdded, ynumberAdded, znumberAdded, ynumberGeneral, yabsoluteIndex, cartCommY, cartCoord, cartDim);

	fourierTranslationZ(tempResultXY, output, direct, xnumberAdded, ynumberAdded, znumberAdded, znumberGeneral, zabsoluteIndex, cartCommZ, cartCoord, cartDim);
}

//////right

Complex fourierTranslationXoneHarmonicRight(Complex*** input, bool direct, int xnumberAdded, int j, int k, int xnumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* xabsoluteIndex, int startAbsoluteIndex, Complex* localFactor, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	if(knumber >= 0){
		for (int i = 1+additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i) {
			if(xabsoluteIndex[i] >= startAbsoluteIndex){
				sum += input[i][j][k] * localFactor[i];
			}
		}
	} else {
		if(knumber + startAbsoluteIndex >= xabsoluteIndex[1+additionalBinNumber] && knumber + startAbsoluteIndex < xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber]){
			sum = input[knumber + startAbsoluteIndex - xabsoluteIndex[0]][j][k];
		}
	}
	double out[2];
	double in[2];
	out[0] = sum.re;
	out[1] = sum.im;
	MPI_Allreduce(out, in, 2, MPI_DOUBLE, MPI_SUM, reducedCartComm);
	sum = Complex(in[0], in[1]);
	if (!direct && knumber >= 0) {
		sum = sum / (xnumberGeneral - startAbsoluteIndex);
	}
	return sum;
}

void fourierTranslationXRight(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int* xabsoluteIndex, int startAbsoluteIndex, Complex* localFactor, MPI_Comm cartCommX, int* cartCoord, int* cartDim) {
	if (xnumberGeneral == 1) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[1+additionalBinNumber][j][k] = input[1+additionalBinNumber][j][k];
			}
		}
	} else {
		for (int knumber =  - startAbsoluteIndex; knumber < xnumberGeneral - startAbsoluteIndex; ++knumber) {
				int phaseFactor = direct ? -1 : 1;
				Complex factor = complexExp((phaseFactor*2*pi*knumber)/(xnumberGeneral - startAbsoluteIndex));
				localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * (xabsoluteIndex[1+additionalBinNumber]-startAbsoluteIndex) * knumber) / (xnumberGeneral - startAbsoluteIndex));
				for(int i = 2 + additionalBinNumber; i < xnumberAdded; ++i){
					localFactor[i] = localFactor[i-1]*factor;
				}
				for (int j = 1+additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
					for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
						Complex f = fourierTranslationXoneHarmonicRight(input, direct, xnumberAdded, j, k, xnumberGeneral, cartCommX, knumber, xabsoluteIndex, startAbsoluteIndex, localFactor, cartCoord, cartDim);
						if(knumber >= xabsoluteIndex[1+additionalBinNumber] - startAbsoluteIndex && knumber < xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber] - startAbsoluteIndex){
							output[knumber - xabsoluteIndex[0] + startAbsoluteIndex][j][k] = f;
						}
					}
				}
			
		}
	}
}

void fourierTranslationYRight(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int ynumberGeneral, int* yabsoluteIndex, int startAbsoluteIndex, int* xabsoluteIndex, Complex* localFactor, MPI_Comm cartCommY, int* cartCoord, int* cartDim) {
	if (ynumberGeneral == 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[i][1+additionalBinNumber][k] = input[i][1+additionalBinNumber][k];
			}
		}
	} else {
		for (int knumber = 0; knumber < ynumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor*2*pi*knumber)/ynumberGeneral);
			localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * yabsoluteIndex[1+additionalBinNumber] * knumber) / ynumberGeneral);
			for(int j = 2 + additionalBinNumber; j < ynumberAdded; ++j){
				localFactor[j] = localFactor[j-1]*factor;
			}
			for (int i = 1+additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
				if(xabsoluteIndex[i] >= startAbsoluteIndex){
					for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
						Complex f = fourierTranslationYoneHarmonic(input, direct, ynumberAdded, i, k, ynumberGeneral, cartCommY, knumber, yabsoluteIndex, localFactor, cartCoord, cartDim);
						if(knumber >= yabsoluteIndex[1+additionalBinNumber] && knumber < yabsoluteIndex[ynumberAdded - 1 - additionalBinNumber]){
							output[i][knumber - yabsoluteIndex[0]][k] = f;
						}
					}
				} else {
					for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
						if(knumber >= yabsoluteIndex[1+additionalBinNumber] && knumber < yabsoluteIndex[ynumberAdded - 1 - additionalBinNumber]){
							output[i][knumber - yabsoluteIndex[0]][k] = input[i][knumber - yabsoluteIndex[0]][k];
						}
					}
				}
			}
		}
	}
}

void fourierTranslationZRight(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int znumberGeneral, int* zabsoluteIndex, int startAbsoluteIndex, int* xabsoluteIndex, Complex* localFactor, MPI_Comm cartCommZ, int* cartCoord, int* cartDim) {
	if (znumberGeneral == 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				output[i][j][1+additionalBinNumber] = input[i][j][1+additionalBinNumber];
			}
		}
	} else {
		for (int knumber = 0; knumber < znumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor*2*pi*knumber)/znumberGeneral);
			localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * zabsoluteIndex[1+additionalBinNumber] * knumber) / znumberGeneral);
			for(int k = 2 + additionalBinNumber; k < znumberAdded; ++k){
				localFactor[k] = localFactor[k-1]*factor;
			}
			for (int i = 1+additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
				if(xabsoluteIndex[i] >= startAbsoluteIndex){
					for (int j = 1+additionalBinNumber; j < ynumberAdded -1 - additionalBinNumber; ++j) {
						Complex f = fourierTranslationZoneHarmonic(input, direct, znumberAdded, i, j, znumberGeneral, cartCommZ, knumber, zabsoluteIndex, localFactor, cartCoord, cartDim);
						if(knumber >= zabsoluteIndex[1+additionalBinNumber] && knumber < zabsoluteIndex[znumberAdded - 1 - additionalBinNumber]){
							output[i][j][knumber - zabsoluteIndex[0]] = f;
						}
					}
				} else {
					for (int j = 1+additionalBinNumber; j < ynumberAdded -1 - additionalBinNumber; ++j) {
						if(knumber >= zabsoluteIndex[1+additionalBinNumber] && knumber < zabsoluteIndex[znumberAdded - 1 - additionalBinNumber]){
							output[i][j][knumber - zabsoluteIndex[0]] = input[i][j][knumber - zabsoluteIndex[0]];
						}
					}
				}
			}
		}
	}
}

void fourierTranslationRight(Complex*** input, Complex*** output, Complex*** tempResultX, Complex*** tempResultXY, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, int startAbsoluteIndex, MPI_Comm& cartComm, MPI_Comm& cartCommX, MPI_Comm& cartCommY, MPI_Comm& cartCommZ , int* xabsoluteIndex, int* yabsoluteIndex, int* zabsoluteIndex, Complex* localFactorX, Complex* localFactorY, Complex* localFactorZ, int* cartCoord, int* cartDim) {

	fourierTranslationXRight(input, tempResultX, direct, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, xabsoluteIndex, startAbsoluteIndex, localFactorX, cartCommX, cartCoord, cartDim);

	fourierTranslationYRight(tempResultX, tempResultXY, direct, xnumberAdded, ynumberAdded, znumberAdded, ynumberGeneral, yabsoluteIndex, startAbsoluteIndex, xabsoluteIndex, localFactorY, cartCommY, cartCoord, cartDim);

	fourierTranslationZRight(tempResultXY, output, direct, xnumberAdded, ynumberAdded, znumberAdded, znumberGeneral, zabsoluteIndex, startAbsoluteIndex, xabsoluteIndex, localFactorZ, cartCommZ, cartCoord, cartDim);
}

/////left

Complex fourierTranslationXoneHarmonicLeft(Complex*** input, bool direct, int xnumberAdded, int j, int k, int xnumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* xabsoluteIndex, int startAbsoluteIndex, Complex* localFactor, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	if(knumber < startAbsoluteIndex){
		for (int i = 1+additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i) {
			if(xabsoluteIndex[i] < startAbsoluteIndex){
				sum += input[i][j][k] * localFactor[i];
			}
		}
	} else {
		if(knumber >= xabsoluteIndex[1+additionalBinNumber] && knumber < xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber]){
			sum = input[knumber - xabsoluteIndex[0]][j][k];
		}
	}
	//double out[2];
	//double in[2];
	//out[0] = sum.re;
	//out[1] = sum.im;
	//MPI_Allreduce(out, in, 2, MPI_DOUBLE, MPI_SUM, reducedCartComm);
	//sum = Complex(in[0], in[1]);
	if (!direct && knumber < startAbsoluteIndex) {
		sum = sum / (startAbsoluteIndex);
	}
	return sum;
}

void fourierTranslationXLeft(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int* xabsoluteIndex, int startAbsoluteIndex, Complex* localFactor, MPI_Comm cartCommX, int* cartCoord, int* cartDim) {
	if (xnumberGeneral == 1) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[1+additionalBinNumber][j][k] = input[1+additionalBinNumber][j][k];
			}
		}
	} else {
		double* result = new double[2*xnumberGeneral];
		double* tempResult = new double[2*xnumberGeneral];
			
				
		for (int j = 1+additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
			for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
				for (int knumber =  0; knumber < xnumberGeneral; ++knumber) {

					int phaseFactor = direct ? -1 : 1;
						Complex factor = complexExp((phaseFactor*2*pi*knumber)/(startAbsoluteIndex));
						localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * (xabsoluteIndex[1+additionalBinNumber]) * knumber) / (startAbsoluteIndex));
						for(int i = 2 + additionalBinNumber; i < xnumberAdded; ++i){
							localFactor[i] = localFactor[i-1]*factor;
						}
							
					Complex f = fourierTranslationXoneHarmonicLeft(input, direct, xnumberAdded, j, k, xnumberGeneral, cartCommX, knumber, xabsoluteIndex, startAbsoluteIndex, localFactor, cartCoord, cartDim);
					result[knumber] = f.re;
					result[knumber + xnumberGeneral] = f.im;
					tempResult[knumber] = 0;
					tempResult[knumber + xnumberGeneral] = 0;
				}
				MPI_Allreduce(result, tempResult, 2*xnumberGeneral, MPI_DOUBLE, MPI_SUM, cartCommX);
				for (int knumber =  0; knumber < xnumberGeneral; ++knumber) {
					if(knumber >= xabsoluteIndex[1+additionalBinNumber] && knumber < xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber]){
						output[knumber - xabsoluteIndex[0]][j][k].re = tempResult[knumber];
						output[knumber - xabsoluteIndex[0]][j][k].im = tempResult[knumber + xnumberGeneral];
					}
				}
			}
		}

		delete[] result;
		delete[] tempResult;
	}
}

void fourierTranslationYLeft(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int ynumberGeneral, int* yabsoluteIndex, int startAbsoluteIndex, int* xabsoluteIndex, Complex* localFactor, MPI_Comm cartCommY, int* cartCoord, int* cartDim) {
	if (ynumberGeneral == 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[i][1+additionalBinNumber][k] = input[i][1+additionalBinNumber][k];
			}
		}
	} else {
		for (int knumber = 0; knumber < ynumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor*2*pi*knumber)/ynumberGeneral);
			localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * yabsoluteIndex[1+additionalBinNumber] * knumber) / ynumberGeneral);
			for(int j = 2 + additionalBinNumber; j < ynumberAdded; ++j){
				localFactor[j] = localFactor[j-1]*factor;
			}
			for (int i = 1+additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
				if(xabsoluteIndex[i] < startAbsoluteIndex){
					for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
						Complex f = fourierTranslationYoneHarmonic(input, direct, ynumberAdded, i, k, ynumberGeneral, cartCommY, knumber, yabsoluteIndex, localFactor, cartCoord, cartDim);
						if(knumber >= yabsoluteIndex[1+additionalBinNumber] && knumber < yabsoluteIndex[ynumberAdded - 1 - additionalBinNumber]){
							output[i][knumber - yabsoluteIndex[0]][k] = f;
						}
					}
				} else {
					for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
						if(knumber >= yabsoluteIndex[1+additionalBinNumber] && knumber < yabsoluteIndex[ynumberAdded - 1 - additionalBinNumber]){
							output[i][knumber - yabsoluteIndex[0]][k] = input[i][knumber - yabsoluteIndex[0]][k];
						}
					}
				}
			}
		}
	}
}

void fourierTranslationZLeft(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int znumberGeneral, int* zabsoluteIndex, int startAbsoluteIndex, int* xabsoluteIndex, Complex* localFactor, MPI_Comm cartCommZ, int* cartCoord, int* cartDim) {
	if (znumberGeneral == 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				output[i][j][1+additionalBinNumber] = input[i][j][1+additionalBinNumber];
			}
		}
	} else {
		for (int knumber = 0; knumber < znumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor*2*pi*knumber)/znumberGeneral);
			localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * zabsoluteIndex[1+additionalBinNumber] * knumber) / znumberGeneral);
			for(int k = 2 + additionalBinNumber; k < znumberAdded; ++k){
				localFactor[k] = localFactor[k-1]*factor;
			}
			for (int i = 1+additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
				if(xabsoluteIndex[i] < startAbsoluteIndex){
					for (int j = 1+additionalBinNumber; j < ynumberAdded -1 - additionalBinNumber; ++j) {
						Complex f = fourierTranslationZoneHarmonic(input, direct, znumberAdded, i, j, znumberGeneral, cartCommZ, knumber, zabsoluteIndex, localFactor, cartCoord, cartDim);
						if(knumber >= zabsoluteIndex[1+additionalBinNumber] && knumber < zabsoluteIndex[znumberAdded - 1 - additionalBinNumber]){
							output[i][j][knumber - zabsoluteIndex[0]] = f;
						}
					}
				} else {
					for (int j = 1+additionalBinNumber; j < ynumberAdded -1 - additionalBinNumber; ++j) {
						if(knumber >= zabsoluteIndex[1+additionalBinNumber] && knumber < zabsoluteIndex[znumberAdded - 1 - additionalBinNumber]){
							output[i][j][knumber - zabsoluteIndex[0]] = input[i][j][knumber - zabsoluteIndex[0]];
						}
					}
				}
			}
		}
	}
}

void fourierTranslationLeft(Complex*** input, Complex*** output, Complex*** tempResultX, Complex*** tempResultXY, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, int startAbsoluteIndex, MPI_Comm& cartComm, MPI_Comm& cartCommX, MPI_Comm& cartCommY, MPI_Comm& cartCommZ , int* xabsoluteIndex, int* yabsoluteIndex, int* zabsoluteIndex, Complex* localFactorX, Complex* localFactorY, Complex* localFactorZ, int* cartCoord, int* cartDim) {

	fourierTranslationXLeft(input, tempResultX, direct, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, xabsoluteIndex, startAbsoluteIndex, localFactorX, cartCommX, cartCoord, cartDim);

	fourierTranslationYLeft(tempResultX, tempResultXY, direct, xnumberAdded, ynumberAdded, znumberAdded, ynumberGeneral, yabsoluteIndex, startAbsoluteIndex, xabsoluteIndex, localFactorY, cartCommY, cartCoord, cartDim);

	fourierTranslationZLeft(tempResultXY, output, direct, xnumberAdded, ynumberAdded, znumberAdded, znumberGeneral, zabsoluteIndex, startAbsoluteIndex, xabsoluteIndex, localFactorZ, cartCommZ, cartCoord, cartDim);
}


//////////////local//////////////////


Complex fourierTranslationXoneHarmonicLocal(Complex*** input, bool direct, int xnumberAdded, int j, int k, int knumber) {
	int localXnumber = xnumberAdded - 2 - 2*additionalBinNumber;
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	int maxI = xnumberAdded;
	for (int i = 1+additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i) {
		sum += input[i][j][k] * complexExp((phaseFactor * 2 * pi * (i - 1 - additionalBinNumber) * knumber) / localXnumber);
	}

	if (!direct) {
		sum = sum / localXnumber;
	}
	return sum;
}

void fourierTranslationXonePointLocal(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int j, int k) {
	for (int knumber = 1 + additionalBinNumber; knumber < xnumberAdded - 1 - additionalBinNumber; ++knumber) {
		int kabsoluteNumber = knumber - 1 - additionalBinNumber;
		Complex f = fourierTranslationXoneHarmonicLocal(input, direct, xnumberAdded, j, k, kabsoluteNumber);
		output[knumber][j][k] = f;
	}
}

void fourierTranslationXLocal(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded) {
	if (xnumberAdded == 1 + 2 + 2*additionalBinNumber) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[1+additionalBinNumber][j][k] = input[1+additionalBinNumber][j][k];
			}
		}
	} else {
		for (int j = 1+additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
			for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
				fourierTranslationXonePointLocal(input, output, direct, xnumberAdded, j, k);
			}
		}
	}
}

Complex fourierTranslationYoneHarmonicLocal(Complex*** input, bool direct, int ynumberAdded, int i, int k, int knumber) {
	int localYnumber = ynumberAdded - 2 - 2*additionalBinNumber;
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	for (int j = 1+additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
		sum += input[i][j][k] * complexExp((phaseFactor * 2 * pi * (j - 1 - additionalBinNumber) * knumber) / localYnumber);
	}

	if (!direct) {
		sum = sum / localYnumber;
	}
	return sum;
}

void fourierTranslationYonePointLocal(Complex*** input, Complex*** output, bool direct, int ynumberAdded, int i, int k) {
	for (int knumber = 1 + additionalBinNumber; knumber < ynumberAdded - 1 - additionalBinNumber; ++knumber) {
		int kabsoluteNumber = knumber - 1 - additionalBinNumber;
		Complex f = fourierTranslationYoneHarmonicLocal(input, direct, ynumberAdded, i, k, kabsoluteNumber);
		output[i][knumber][k] = f;
	}
}

void fourierTranslationYLocal(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded) {
	if (ynumberAdded == 1 + 2 + 2*additionalBinNumber) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[i][1+additionalBinNumber][k] = input[i][1+additionalBinNumber][k];
			}
		}
	} else {
		for (int i = 1+additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
			for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
				fourierTranslationYonePointLocal(input, output, direct, ynumberAdded, i, k);
			}
		}
	}
}

Complex fourierTranslationZoneHarmonicLocal(Complex*** input, bool direct, int znumberAdded, int i, int j, int knumber) {
	int localZnumber = znumberAdded - 2 - 2*additionalBinNumber;
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	for (int k = 1+additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
		sum += input[i][j][k] * complexExp((phaseFactor * 2 * pi * (k - 1 - additionalBinNumber) * knumber) / localZnumber);
	}
	if (!direct) {
		sum = sum / localZnumber;
	}
	return sum;
}

void fourierTranslationZonePointLocal(Complex*** input, Complex*** output, bool direct, int znumberAdded, int i, int j) {
	for (int knumber = 1 + additionalBinNumber; knumber < znumberAdded - 1 - additionalBinNumber; ++knumber) {
		int kabsoluteNumber = knumber - 1 - additionalBinNumber;
		Complex f = fourierTranslationZoneHarmonicLocal(input, direct, znumberAdded, i, j, kabsoluteNumber);
		output[i][j][knumber] = f;
	}
}

void fourierTranslationZLocal(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded) {
	if (znumberAdded == 1 + 2 + 2*additionalBinNumber) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				output[i][j][1+additionalBinNumber] = input[i][j][1+additionalBinNumber];
			}
		}
	} else {
		for (int i = 1+additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i) {
			for (int j = 1+additionalBinNumber; j < ynumberAdded -1 - additionalBinNumber; ++j) {
				fourierTranslationZonePointLocal(input, output, direct, znumberAdded, i, j);
			}
		}
	}
}

void fourierTranslationLocal(Complex*** input, Complex*** output, Complex*** tempResultX, Complex*** tempResultXY, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded) {


	fourierTranslationXLocal(input, tempResultX, direct, xnumberAdded, ynumberAdded, znumberAdded);

	fourierTranslationYLocal(tempResultX, tempResultXY, direct, xnumberAdded, ynumberAdded, znumberAdded);

	fourierTranslationZLocal(tempResultXY, output, direct, xnumberAdded, ynumberAdded, znumberAdded);
}

///mirror

Complex fourierTranslationXoneHarmonicRightMirror(Complex*** input, bool direct, int xnumberAdded, int j, int k, int xnumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* xabsoluteIndex, int startAbsoluteIndex, Complex* localFactor, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	if(knumber >= 0 && knumber < 2*xnumberGeneral - 2*startAbsoluteIndex){
		for (int i = 1+additionalBinNumber; i < 2*xnumberAdded - 3*additionalBinNumber - 3; ++i) {
			if(xabsoluteIndex[i] >= startAbsoluteIndex && xabsoluteIndex[i] < 2*xnumberGeneral - startAbsoluteIndex){
				sum += input[i][j][k] * (localFactor[i]);
			}
		}
	} else if(knumber < 0){
		if(knumber + startAbsoluteIndex >= xabsoluteIndex[1+additionalBinNumber] && knumber + startAbsoluteIndex < xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber]){
			sum = input[knumber + startAbsoluteIndex - xabsoluteIndex[1 + additionalBinNumber]][j][k];
		}
	} else {
		if(knumber + startAbsoluteIndex >= 2*xnumberGeneral - xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber] && knumber + startAbsoluteIndex < 2*xnumberGeneral - xabsoluteIndex[1+additionalBinNumber]){
			sum = input[-(knumber + startAbsoluteIndex - 2*xnumberGeneral + xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber]) + 2*xnumberAdded - 2 - 2*additionalBinNumber - 1][j][k];
		}
	}
	double out[2];
	double in[2];
	out[0] = sum.re;
	out[1] = sum.im;
	MPI_Allreduce(out, in, 2, MPI_DOUBLE, MPI_SUM, reducedCartComm);
	sum = Complex(in[0], in[1]);
	if (!direct && knumber >= 0 && knumber < 2*xnumberGeneral - 2*startAbsoluteIndex) {
		sum = sum / (2*(xnumberGeneral - startAbsoluteIndex));
	}
	return sum;
}

void fourierTranslationXRightMirror(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int* xabsoluteIndex, int startAbsoluteIndex, Complex* localFactor, MPI_Comm cartCommX, int* cartCoord, int* cartDim) {
	if (xnumberGeneral == 1) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[1+additionalBinNumber][j][k] = input[1+additionalBinNumber][j][k];
			}
		}
	} else {
		for (int knumber = - startAbsoluteIndex; knumber < 2*(xnumberGeneral - startAbsoluteIndex) + startAbsoluteIndex; ++knumber) {
				int phaseFactor = direct ? -1 : 1;
				Complex factor = complexExp((phaseFactor*2*pi*knumber)/(2*(xnumberGeneral - startAbsoluteIndex)));
				localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * (xabsoluteIndex[1+additionalBinNumber]-startAbsoluteIndex) * knumber) / (2*(xnumberGeneral - startAbsoluteIndex)));
				for(int i = 2 + additionalBinNumber; i < xnumberAdded - 1 - additionalBinNumber; ++i){
					localFactor[i] = localFactor[i-1]*factor;
				}
				localFactor[xnumberAdded - 1 - additionalBinNumber] = complexExp((phaseFactor * 2 * pi * (2*xnumberGeneral - xabsoluteIndex[1 + additionalBinNumber]- 1 - startAbsoluteIndex) * knumber) / (2*(xnumberGeneral - startAbsoluteIndex)));
				for(int i = xnumberAdded - additionalBinNumber; i < 2*xnumberAdded; ++i){
					localFactor[i] = localFactor[i-1]/factor;
				}
				for (int j = 1+additionalBinNumber; j < ynumberAdded - 1 - additionalBinNumber; ++j) {
					for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
						Complex f = fourierTranslationXoneHarmonicRightMirror(input, direct, xnumberAdded, j, k, xnumberGeneral, cartCommX, knumber, xabsoluteIndex, startAbsoluteIndex, localFactor, cartCoord, cartDim);
						if((knumber + startAbsoluteIndex) >= xabsoluteIndex[1+additionalBinNumber] && (knumber + startAbsoluteIndex) < xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber]){
							output[knumber - xabsoluteIndex[0] + startAbsoluteIndex][j][k] = f;
						} else if(knumber + startAbsoluteIndex >= 2*xnumberGeneral - xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber] && knumber + startAbsoluteIndex < 2*xnumberGeneral - xabsoluteIndex[1+additionalBinNumber]){
							output[-(knumber + startAbsoluteIndex - 2*xnumberGeneral + xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber]) + 2*xnumberAdded - 2 - 2*additionalBinNumber - 1][j][k] = f;
						}
					}
				}
			
		}
	}
}

void fourierTranslationYRightMirror(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int ynumberGeneral, int* yabsoluteIndex, int startAbsoluteIndex, int* xabsoluteIndex, Complex* localFactor, MPI_Comm cartCommY, int* cartCoord, int* cartDim) {
	if (ynumberGeneral == 1) {
		for (int i = 0; i < 2*xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[i][1+additionalBinNumber][k] = input[i][1+additionalBinNumber][k];
			}
		}
	} else {
		for (int knumber = 0; knumber < ynumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor*2*pi*knumber)/ynumberGeneral);
			localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * yabsoluteIndex[1+additionalBinNumber] * knumber) / ynumberGeneral);
			for(int j = 2 + additionalBinNumber; j < ynumberAdded; ++j){
				localFactor[j] = localFactor[j-1]*factor;
			}
			for (int i = 1+additionalBinNumber; i < 2*xnumberAdded - 3 - 3*additionalBinNumber; ++i) {
				if(xabsoluteIndex[i] >= startAbsoluteIndex && xabsoluteIndex[i] < 2*xnumberGeneral - startAbsoluteIndex){
					for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
						Complex f = fourierTranslationYoneHarmonic(input, direct, ynumberAdded, i, k, ynumberGeneral, cartCommY, knumber, yabsoluteIndex, localFactor, cartCoord, cartDim);
						if(knumber >= yabsoluteIndex[1+additionalBinNumber] && knumber < yabsoluteIndex[ynumberAdded - 1 - additionalBinNumber]){
							output[i][knumber - yabsoluteIndex[0]][k] = f;
						}
					}
				} else {
					for (int k = 1+additionalBinNumber; k < znumberAdded - 1 - additionalBinNumber; ++k) {
						if(knumber >= yabsoluteIndex[1+additionalBinNumber] && knumber < yabsoluteIndex[ynumberAdded - 1 - additionalBinNumber]){
							output[i][knumber - yabsoluteIndex[0]][k] = input[i][knumber - yabsoluteIndex[0]][k];
						}
					}
				}
			}
		}
	}
}

void fourierTranslationZRightMirror(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int znumberGeneral, int* zabsoluteIndex, int startAbsoluteIndex, int* xabsoluteIndex, Complex* localFactor, MPI_Comm cartCommZ, int* cartCoord, int* cartDim) {
	if (znumberGeneral == 1) {
		for (int i = 0; i < 2*xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				output[i][j][1+additionalBinNumber] = input[i][j][1+additionalBinNumber];
			}
		}
	} else {
		for (int knumber = 0; knumber < znumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor*2*pi*knumber)/znumberGeneral);
			localFactor[1+additionalBinNumber] = complexExp((phaseFactor * 2 * pi * zabsoluteIndex[1+additionalBinNumber] * knumber) / znumberGeneral);
			for(int k = 2 + additionalBinNumber; k < znumberAdded; ++k){
				localFactor[k] = localFactor[k-1]*factor;
			}
			for (int i = 1+additionalBinNumber; i < 2*xnumberAdded - 3 - 3*additionalBinNumber; ++i) {
				if(xabsoluteIndex[i] >= startAbsoluteIndex && xabsoluteIndex[i] < 2*xnumberGeneral - startAbsoluteIndex){
					for (int j = 1+additionalBinNumber; j < ynumberAdded -1 - additionalBinNumber; ++j) {
						Complex f = fourierTranslationZoneHarmonic(input, direct, znumberAdded, i, j, znumberGeneral, cartCommZ, knumber, zabsoluteIndex, localFactor, cartCoord, cartDim);
						if(knumber >= zabsoluteIndex[1+additionalBinNumber] && knumber < zabsoluteIndex[znumberAdded - 1 - additionalBinNumber]){
							output[i][j][knumber - zabsoluteIndex[0]] = f;
						}
					}
				} else {
					for (int j = 1+additionalBinNumber; j < ynumberAdded -1 - additionalBinNumber; ++j) {
						if(knumber >= zabsoluteIndex[1+additionalBinNumber] && knumber < zabsoluteIndex[znumberAdded - 1 - additionalBinNumber]){
							output[i][j][knumber - zabsoluteIndex[0]] = input[i][j][knumber - zabsoluteIndex[0]];
						}
					}
				}
			}
		}
	}
}

void fourierTranslationRightMirror(Complex*** input, Complex*** output, Complex*** tempResultX, Complex*** tempResultXY, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, int startAbsoluteIndex, MPI_Comm& cartComm, MPI_Comm& cartCommX, MPI_Comm& cartCommY, MPI_Comm& cartCommZ , int* xabsoluteIndex, int* yabsoluteIndex, int* zabsoluteIndex, Complex* localFactorX, Complex* localFactorY, Complex* localFactorZ, int* cartCoord, int* cartDim) {

	fourierTranslationXRightMirror(input, tempResultX, direct, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, xabsoluteIndex, startAbsoluteIndex, localFactorX, cartCommX, cartCoord, cartDim);

	fourierTranslationYRightMirror(tempResultX, tempResultXY, direct, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, ynumberGeneral, yabsoluteIndex, startAbsoluteIndex, xabsoluteIndex, localFactorY, cartCommY, cartCoord, cartDim);

	fourierTranslationZRightMirror(tempResultXY, output, direct, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral,znumberGeneral, zabsoluteIndex, startAbsoluteIndex, xabsoluteIndex, localFactorZ, cartCommZ, cartCoord, cartDim);
}
