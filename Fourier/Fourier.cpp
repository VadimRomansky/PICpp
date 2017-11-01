#include "stdio.h"
#include <stdlib.h>
#include <string>

#include "fourier.h"
#include "complex.h"
#include "constants.h"
#include "output.h"

Complex fourierTranslationXoneHarmonic(Complex*** input, bool direct, int xnumber, int j, int k, int xnumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* xabsoluteIndex) {
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	for (int i = 0; i < xnumber; ++i) {
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

void fourierTranslationXonePoint(Complex*** input, Complex*** output, bool direct, int xnumber, int j, int k, int xnumberGeneral, MPI_Comm& reducedCartComm, int* xabsoluteIndex) {
	for (int knumber = 0; knumber < xnumberGeneral; ++knumber) {
		Complex f = fourierTranslationXoneHarmonic(input, direct, xnumber, j, k, xnumberGeneral, reducedCartComm, knumber, xabsoluteIndex);
		if (knumber >= xabsoluteIndex[0] && knumber <= xabsoluteIndex[xnumber - 1]) {
			int i = knumber - xabsoluteIndex[0];
			output[i][j][k] = f;
		}
	}
}

void fourierTranslationX(Complex*** input, Complex*** output, bool direct, int xnumber, int ynumber, int znumber, int xnumberGeneral, int* xabsoluteIndex, MPI_Comm cartCommX) {
	if (xnumberGeneral == 1) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				output[0][j][k] = input[0][j][k];
			}
		}
	} else {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				fourierTranslationXonePoint(input, output, direct, xnumber, j, k, xnumberGeneral, cartCommX, xabsoluteIndex);
			}
		}
	}
}

Complex fourierTranslationYoneHarmonic(Complex*** input, bool direct, int ynumber, int i, int k, int ynumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* yabsoluteIndex) {
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	for (int j = 0; j < ynumber; ++j) {
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

void fourierTranslationYonePoint(Complex*** input, Complex*** output, bool direct, int ynumber, int i, int k, int ynumberGeneral, MPI_Comm& reducedCartComm, int* yabsoluteIndex) {
	for (int knumber = 0; knumber < ynumberGeneral; ++knumber) {
		Complex f = fourierTranslationYoneHarmonic(input, direct, ynumber, i, k, ynumberGeneral, reducedCartComm, knumber, yabsoluteIndex);
		if (knumber >= yabsoluteIndex[0] && knumber <= yabsoluteIndex[ynumber - 1]) {
			int j = knumber - yabsoluteIndex[0];
			output[i][j][k] = f;
		}
	}
}

void fourierTranslationY(Complex*** input, Complex*** output, bool direct, int xnumber, int ynumber, int znumber, int ynumberGeneral, int* yabsoluteIndex, MPI_Comm cartCommY) {
	if (ynumberGeneral == 1) {
		for (int i = 0; i < xnumber; ++i) {
			for (int k = 0; k < znumber; ++k) {
				output[i][0][k] = input[i][0][k];
			}
		}
	} else {
		for (int i = 0; i < xnumber; ++i) {
			for (int k = 0; k < znumber; ++k) {
				fourierTranslationYonePoint(input, output, direct, ynumber, i, k, ynumberGeneral, cartCommY, yabsoluteIndex);
			}
		}
	}
}

Complex fourierTranslationZoneHarmonic(Complex*** input, bool direct, int znumber, int i, int j, int znumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* zabsoluteIndex) {
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	for (int k = 0; k < znumber; ++k) {
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

void fourierTranslationZonePoint(Complex*** input, Complex*** output, bool direct, int znumber, int i, int j, int znumberGeneral, MPI_Comm& reducedCartComm, int* zabsoluteIndex) {
	for (int knumber = 0; knumber < znumberGeneral; ++knumber) {
		Complex f = fourierTranslationZoneHarmonic(input, direct, znumber, i, j, znumberGeneral, reducedCartComm, knumber, zabsoluteIndex);
		if (knumber >= zabsoluteIndex[0] && knumber <= zabsoluteIndex[znumber - 1]) {
			int k = knumber - zabsoluteIndex[0];
			output[i][j][k] = f;
		}
	}
}

void fourierTranslationZ(Complex*** input, Complex*** output, bool direct, int xnumber, int ynumber, int znumber, int znumberGeneral, int* zabsoluteIndex, MPI_Comm cartCommZ) {
	if (znumberGeneral == 1) {
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				output[i][j][0] = input[i][j][0];
			}
		}
	} else {
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				fourierTranslationZonePoint(input, output, direct, znumber, i, j, znumberGeneral, cartCommZ, zabsoluteIndex);
			}
		}
	}
}

void fourierTranslation(Complex*** input, Complex*** output, bool direct, int xnumber, int ynumber, int znumber, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, MPI_Comm& cartComm, int* xabsoluteIndex, int* yabsoluteIndex, int* zabsoluteIndex, int* cartCoord, int* cartDim) {
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

	Complex*** tempResult = new Complex**[xnumber];
	Complex*** tempResult1 = new Complex**[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		tempResult[i] = new Complex*[ynumber];
		tempResult1[i] = new Complex*[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			tempResult[i][j] = new Complex[znumber];
			tempResult1[i][j] = new Complex[znumber];
			for (int k = 0; k < znumber; ++k) {
				tempResult[i][j][k] = Complex(0, 0);
				tempResult1[i][j][k] = Complex(0, 0);
			}
		}
	}

	fourierTranslationX(input, tempResult, direct, xnumber, ynumber, znumber, xnumberGeneral, xabsoluteIndex, cartCommX);

	fourierTranslationY(tempResult, tempResult1, direct, xnumber, ynumber, znumber, ynumberGeneral, yabsoluteIndex, cartCommY);

	fourierTranslationZ(tempResult1, output, direct, xnumber, ynumber, znumber, znumberGeneral, zabsoluteIndex, cartCommZ);

	/*std::string outputDir = outputDirectory;

	if(direct){
		outputArray((outputDir + "fourierX.dat").c_str(), tempResult, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
		outputArray((outputDir + "fourierY.dat").c_str(), tempResult1, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
		outputArray((outputDir + "fourierZ.dat").c_str(), output, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
	}*/

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			delete[] tempResult[i][j];
			delete[] tempResult1[i][j];
		}
		delete[] tempResult[i];
		delete[] tempResult1[i];
	}
	delete[] tempResult;
	delete[] tempResult1;
}

Complex fastFourier1d(Complex* input, Complex* output, bool direct, int xnumberAdded, int N, int xnumberGeneral, MPI_Comm& comm, int* cartCoord, int* cartDim, int start, int step){
	if(N == 1){
		output[start] = input[start];
	} else {
		fastFourier1d(input, output, direct, xnumberAdded, N/2, xnumberGeneral, comm, cartCoord, cartDim, start, step*2);
		fastFourier1d(input, output, direct, xnumberAdded, N/2, xnumberGeneral, comm, cartCoord, cartDim, start + step, step*2);
		int sign = direct ? -1 : 1;
		for(int k = 0; k < N/2; ++k){
			Complex factor = complexExp((sign*2*pi*k)/N);
			Complex t = output[start + step*k];
			output[start + step*k] = t + factor*output[start + step*(k + N/2)];
			output[start + step*(k+N/2)) = t - factor*output[start + step*(k + N/2)];
		}
	}
}
