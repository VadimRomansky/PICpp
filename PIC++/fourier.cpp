#include "stdio.h"
#include <stdlib.h>
#include <string>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "fourier.h"
#include "complex.h"
#include "constants.h"
#include "output.h"

Complex fourierTranslationXoneHarmonic(Complex*** input, bool direct, int xnumberAdded, int j, int k, int xnumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* xabsoluteIndex, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	int maxI = xnumberAdded;
	for (int i = 1+additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i) {
		sum += input[i][j][k] * complexExp(phaseFactor * 2 * pi * xabsoluteIndex[i] * knumber / xnumberGeneral);
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

void fourierTranslationXonePoint(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int j, int k, int xnumberGeneral, MPI_Comm& reducedCartComm, int* xabsoluteIndex, int* cartCoord, int* cartDim) {
	for (int knumber = 0; knumber < xnumberGeneral; ++knumber) {
		Complex f = fourierTranslationXoneHarmonic(input, direct, xnumberAdded, j, k, xnumberGeneral, reducedCartComm, knumber, xabsoluteIndex, cartCoord, cartDim);
		if (knumber >= xabsoluteIndex[1+additionalBinNumber] && knumber < xabsoluteIndex[xnumberAdded - additionalBinNumber - 1]) {
			int i = knumber - xabsoluteIndex[0];
			output[i][j][k] = f;
		}
	}
}

void fourierTranslationX(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int* xabsoluteIndex, MPI_Comm cartCommX, int* cartCoord, int* cartDim) {
	if (xnumberGeneral == 1) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[0][j][k] = input[0][j][k];
			}
		}
	} else {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				fourierTranslationXonePoint(input, output, direct, xnumberAdded, j, k, xnumberGeneral, cartCommX, xabsoluteIndex, cartCoord, cartDim);
			}
		}
	}
}

Complex fourierTranslationYoneHarmonic(Complex*** input, bool direct, int ynumberAdded, int i, int k, int ynumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* yabsoluteIndex, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	for (int j = 1+additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
		sum += input[i][j][k] * complexExp(phaseFactor * 2 * pi * yabsoluteIndex[j] * knumber / ynumberGeneral);
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

void fourierTranslationYonePoint(Complex*** input, Complex*** output, bool direct, int ynumberAdded, int i, int k, int ynumberGeneral, MPI_Comm& reducedCartComm, int* yabsoluteIndex, int* cartCoord, int* cartDim) {
	for (int knumber = 0; knumber < ynumberGeneral; ++knumber) {
		Complex f = fourierTranslationYoneHarmonic(input, direct, ynumberAdded, i, k, ynumberGeneral, reducedCartComm, knumber, yabsoluteIndex, cartCoord, cartDim);
		if (knumber >= yabsoluteIndex[1+additionalBinNumber] && knumber < yabsoluteIndex[ynumberAdded - additionalBinNumber - 1]) {
			int j = knumber - yabsoluteIndex[0];
			output[i][j][k] = f;
		}
	}
}

void fourierTranslationY(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int ynumberGeneral, int* yabsoluteIndex, MPI_Comm cartCommY, int* cartCoord, int* cartDim) {
	if (ynumberGeneral == 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				output[i][0][k] = input[i][0][k];
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int k = 0; k < znumberAdded; ++k) {
				fourierTranslationYonePoint(input, output, direct, ynumberAdded, i, k, ynumberGeneral, cartCommY, yabsoluteIndex, cartCoord, cartDim);
			}
		}
	}
}

Complex fourierTranslationZoneHarmonic(Complex*** input, bool direct, int znumberAdded, int i, int j, int znumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* zabsoluteIndex, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	for (int k = 1+additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
		sum += input[i][j][k] * complexExp(phaseFactor * 2 * pi * zabsoluteIndex[k] * knumber / znumberGeneral);
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

void fourierTranslationZonePoint(Complex*** input, Complex*** output, bool direct, int znumberAdded, int i, int j, int znumberGeneral, MPI_Comm& reducedCartComm, int* zabsoluteIndex, int* cartCoord, int* cartDim) {
	for (int knumber = 0; knumber < znumberGeneral; ++knumber) {
		Complex f = fourierTranslationZoneHarmonic(input, direct, znumberAdded, i, j, znumberGeneral, reducedCartComm, knumber, zabsoluteIndex, cartCoord, cartDim);
		if (knumber >= zabsoluteIndex[0] && knumber < zabsoluteIndex[znumberAdded - additionalBinNumber - 1]) {
			int k = knumber - zabsoluteIndex[0];
			output[i][j][k] = f;
		}
	}
}

void fourierTranslationZ(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int znumberGeneral, int* zabsoluteIndex, MPI_Comm cartCommZ, int* cartCoord, int* cartDim) {
	if (znumberGeneral == 1) {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				output[i][j][0] = input[i][j][0];
			}
		}
	} else {
		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				fourierTranslationZonePoint(input, output, direct, znumberAdded, i, j, znumberGeneral, cartCommZ, zabsoluteIndex, cartCoord, cartDim);
			}
		}
	}
}

void fourierTranslation(Complex*** input, Complex*** output, bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, MPI_Comm& cartComm, int* xabsoluteIndex, int* yabsoluteIndex, int* zabsoluteIndex, int* cartCoord, int* cartDim) {
	MPI_Comm cartCommX;
	MPI_Comm cartCommY;
	MPI_Comm cartCommZ;

	int dims[MPI_dim];
	for (int i = 0; i < MPI_dim; ++i) {
		dims[i] = 0;
	}

	dims[0] = 1;
	MPI_Cart_sub(cartComm, dims, &cartCommX);
	dims[0] = 0;
	dims[1] = 1;
	MPI_Cart_sub(cartComm, dims, &cartCommY);
	dims[1] = 0;
	dims[2] = 1;
	MPI_Cart_sub(cartComm, dims, &cartCommZ);

	Complex*** tempResult = new Complex**[xnumberAdded];
	Complex*** tempResult1 = new Complex**[xnumberAdded];
	for (int i = 0; i < xnumberAdded; ++i) {
		tempResult[i] = new Complex*[ynumberAdded];
		tempResult1[i] = new Complex*[ynumberAdded];
		for (int j = 0; j < ynumberAdded; ++j) {
			tempResult[i][j] = new Complex[znumberAdded];
			tempResult1[i][j] = new Complex[znumberAdded];
			for (int k = 0; k < znumberAdded; ++k) {
				tempResult[i][j][k] = Complex(0, 0);
				tempResult1[i][j][k] = Complex(0, 0);
			}
		}
	}

	fourierTranslationX(input, tempResult, direct, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral, xabsoluteIndex, cartCommX, cartCoord, cartDim);

	fourierTranslationY(tempResult, tempResult1, direct, xnumberAdded, ynumberAdded, znumberAdded, ynumberGeneral, yabsoluteIndex, cartCommY, cartCoord, cartDim);

	fourierTranslationZ(tempResult1, output, direct, xnumberAdded, ynumberAdded, znumberAdded, znumberGeneral, zabsoluteIndex, cartCommZ, cartCoord, cartDim);

	/*std::string outputDir = outputDirectory;

	if(direct){
		outputArray((outputDir + "fourierX.dat").c_str(), tempResult, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
		outputArray((outputDir + "fourierY.dat").c_str(), tempResult1, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
		outputArray((outputDir + "fourierZ.dat").c_str(), output, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
	}*/

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			delete[] tempResult[i][j];
			delete[] tempResult1[i][j];
		}
		delete[] tempResult[i];
		delete[] tempResult1[i];
	}
	delete[] tempResult;
	delete[] tempResult1;
}
