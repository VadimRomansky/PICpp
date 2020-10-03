#include "stdio.h"
#include <stdlib.h>
#include <string>
#include "math.h"

#include "fourier.h"
#include "complex.h"
#include "constants.h"
#include "output.h"

Complex fourierTranslationXoneHarmonic(Complex*** input, bool direct, int xnumberAdded, int j, int k,
                                       int xnumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* xabsoluteIndex,
                                       Complex* localFactor, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	int phaseFactor = direct ? -1 : 1;
	int maxI = xnumberAdded;
	for (int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i) {
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

void fourierTranslationX(Complex*** input, Complex*** output, bool direct, int xnumber, int ynumber,
                         int znumber, int xnumberGeneral, int* xabsoluteIndex, MPI_Comm cartCommX, int* cartCoord,
                         int* cartDim) {
	if (xnumberGeneral == 1) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				output[0][j][k] = input[0][j][k];
			}
		}
	} else {
		Complex* localFactor = new Complex[xnumber];
		double* sum = new double[2*ynumber*znumber];
		double* sumAll = new double[2*ynumber*znumber];
		for (int knumber = 0; knumber < xnumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor * 2 * pi * knumber) / xnumberGeneral);
			localFactor[0] = complexExp(
				(phaseFactor * 2 * pi * xabsoluteIndex[0] * knumber) / xnumberGeneral);
			for (int i = 1; i < xnumber; ++i) {
				localFactor[i] = localFactor[i - 1] * factor;
			}
			int l = 0;
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					sum[l] = 0;
					sum[l+1] = 0;
					Complex sum1;
					sum1 = Complex(0, 0);
					for (int i = 0; i < xnumber; ++i) {
						sum1 += input[i][j][k] * localFactor[i];
					}
					
					sum[l] = sum1.re;
					sum[l+1] = sum1.im;
					l=l+2;
				}
			}
			MPI_Allreduce(sum, sumAll, 2*ynumber*znumber, MPI_DOUBLE, MPI_SUM, cartCommX);
			l=0;
			if((knumber >= xabsoluteIndex[0]) && (knumber < xabsoluteIndex[0] + xnumber)){
				for (int j = 0; j < ynumber; ++j) {
					for (int k = 0; k < znumber; ++k) {
						output[knumber-xabsoluteIndex[0]][j][k] = Complex(sumAll[l], sumAll[l+1]);
						if (!direct) {
							output[knumber-xabsoluteIndex[0]][j][k] = output[knumber-xabsoluteIndex[0]][j][k] / xnumberGeneral;
						}
						l=l+2;
					}
				}
			}
		}
		delete[] localFactor;
		delete[] sum;
		delete[] sumAll;
	}
}

Complex fourierTranslationYoneHarmonic(Complex*** input, bool direct, int ynumberAdded, int i, int k,
                                       int ynumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* yabsoluteIndex,
                                       Complex* localFactor, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
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

void fourierTranslationY(Complex*** input, Complex*** output, bool direct, int xnumber, int ynumber,
                         int znumber, int ynumberGeneral, int* yabsoluteIndex, MPI_Comm cartCommY, int* cartCoord,
                         int* cartDim) {
	if (ynumberGeneral == 1) {
		for (int i = 0; i < xnumber; ++i) {
			for (int k = 0; k < znumber; ++k) {
				output[i][0][k] = input[i][0][k];
			}
		}
	} else {
		Complex* localFactor = new Complex[ynumber];
		double* sum = new double[2*xnumber*znumber];
		double* sumAll = new double[2*xnumber*znumber];
		for (int knumber = 0; knumber < ynumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor * 2 * pi * knumber) / ynumberGeneral);
			localFactor[0] = complexExp(
				(phaseFactor * 2 * pi * yabsoluteIndex[0] * knumber) / ynumberGeneral);
			for (int j = 1; j < ynumber; ++j) {
				localFactor[j] = localFactor[j - 1] * factor;
			}
			int l = 0;
			for (int i = 0; i < xnumber; ++i) {
				for (int k = 0; k < znumber; ++k) {
					sum[l] = 0;
					sum[l+1] = 0;
					Complex sum1;
					sum1 = Complex(0, 0);
					for (int j = 0; j < ynumber; ++j) {
						sum1 += input[i][j][k] * localFactor[j];
					}
					
					sum[l] = sum1.re;
					sum[l+1] = sum1.im;
					l=l+2;
				}
			}
			MPI_Allreduce(sum, sumAll, 2*xnumber*znumber, MPI_DOUBLE, MPI_SUM, cartCommY);
			l=0;
			if((knumber >= yabsoluteIndex[0]) && (knumber < yabsoluteIndex[0] + ynumber)){
				for (int i = 0; i < xnumber; ++i) {
					for (int k = 0; k < znumber; ++k) {
						output[i][knumber-yabsoluteIndex[0]][k] = Complex(sumAll[l], sumAll[l+1]);
						if (!direct) {
							output[i][knumber-yabsoluteIndex[0]][k] = output[i][knumber-yabsoluteIndex[0]][k] / ynumberGeneral;
						}
						l=l+2;
					}
				}
			}
		}
		delete[] localFactor;
		delete[] sum;
		delete[] sumAll;
	}
}

Complex fourierTranslationZoneHarmonic(Complex*** input, bool direct, int znumberAdded, int i, int j,
                                       int znumberGeneral, MPI_Comm& reducedCartComm, int knumber, int* zabsoluteIndex,
                                       Complex* localFactor, int* cartCoord, int* cartDim) {
	Complex sum;
	sum = Complex(0, 0);
	for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
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

void fourierTranslationZ(Complex*** input, Complex*** output, bool direct, int xnumber, int ynumber,
                         int znumber, int znumberGeneral, int* zabsoluteIndex, MPI_Comm cartCommZ, int* cartCoord,
                         int* cartDim) {
	if (znumberGeneral == 1) {
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				output[i][j][0] = input[i][j][0];
			}
		}
	} else {
		Complex* localFactor = new Complex[znumber];
		double* sum = new double[2*xnumber*ynumber];
		double* sumAll = new double[2*xnumber*ynumber];
		for (int knumber = 0; knumber < znumberGeneral; ++knumber) {
			int phaseFactor = direct ? -1 : 1;
			Complex factor = complexExp((phaseFactor * 2 * pi * knumber) / znumberGeneral);
			localFactor[0] = complexExp(
				(phaseFactor * 2 * pi * zabsoluteIndex[0] * knumber) / znumberGeneral);
			for (int k = 1; k < znumber; ++k) {
				localFactor[k] = localFactor[k - 1] * factor;
			}
			int l = 0;
			for (int i = 0; i < xnumber; ++i) {
				for (int j = 0; j < ynumber; ++j) {
					sum[l] = 0;
					sum[l+1] = 0;
					Complex sum1;
					sum1 = Complex(0, 0);
					for (int k = 0; k < znumber; ++k) {
						sum1 += input[i][j][k] * localFactor[k];
					}
					
					sum[l] = sum1.re;
					sum[l+1] = sum1.im;
					l=l+2;
				}
			}
			MPI_Allreduce(sum, sumAll, 2*xnumber*znumber, MPI_DOUBLE, MPI_SUM, cartCommZ);
			l=0;
			if((knumber >= zabsoluteIndex[0]) && (knumber < zabsoluteIndex[0] + znumber)){
				for (int i = 0; i < xnumber; ++i) {
					for (int j = 0; j < ynumber; ++j) {
						output[i][j][knumber-zabsoluteIndex[0]] = Complex(sumAll[l], sumAll[l+1]);
						if (!direct) {
							output[i][j][knumber-zabsoluteIndex[0]] = output[i][j][knumber-zabsoluteIndex[0]] / znumberGeneral;
						}
						l=l+2;
					}
				}
			}
		}
		delete[] localFactor;
		delete[] sum;
		delete[] sumAll;
	}
}

void fourierTranslation(Complex*** input, Complex*** output, Complex*** tempResultX, Complex*** tempResultXY,
                        bool direct, int xnumber, int ynumber, int znumber, int xnumberGeneral,
                        int ynumberGeneral, int znumberGeneral, MPI_Comm& cartComm, MPI_Comm& cartCommX,
                        MPI_Comm& cartCommY, MPI_Comm& cartCommZ, int* xabsoluteIndex, int* yabsoluteIndex,
                        int* zabsoluteIndex, int* cartCoord, int* cartDim) {

	fourierTranslationX(input, tempResultX, direct, xnumber, ynumber, znumber, xnumberGeneral,
	                    xabsoluteIndex, cartCommX, cartCoord, cartDim);

	fourierTranslationY(tempResultX, tempResultXY, direct, xnumber, ynumber, znumber, ynumberGeneral,
	                    yabsoluteIndex, cartCommY, cartCoord, cartDim);

	fourierTranslationZ(tempResultXY, output, direct, xnumber, ynumber, znumber, znumberGeneral,
	                    zabsoluteIndex, cartCommZ, cartCoord, cartDim);
}

void fastFourier1d(Complex* input, Complex* output, bool direct, int xnumberAdded, int N, int xnumberGeneral, MPI_Comm& comm, int* cartCoord, int* cartDim, int start, int step, int* xabsoluteIndex, int fourierStart){
	if(N == 1){
		//if(start >= xabsoluteIndex[1+additionalBinNumber] && start < xabsoluteIndex[xnumberAdded - 1 - additionalBinNumber]){
			output[fourierStart - xabsoluteIndex[0]] = input[start - xabsoluteIndex[0]];
			if(!direct){
				output[fourierStart - xabsoluteIndex[0]] = output[fourierStart - xabsoluteIndex[0]]/xnumberGeneral;
			}
		//}
	} else {
		fastFourier1d(input, output, direct, xnumberAdded, N/2, xnumberGeneral, comm, cartCoord, cartDim, start, step*2, xabsoluteIndex, fourierStart);
		fastFourier1d(input, output, direct, xnumberAdded, N/2, xnumberGeneral, comm, cartCoord, cartDim, start + step, step*2, xabsoluteIndex, fourierStart + N/2);
		int sign = direct ? -1 : 1;
		Complex factor = complexExp((sign*2*pi)/N);
		Complex currentFactor = Complex(1.0, 0);
		for(int k = 0; k < N/2; ++k){
			Complex factorOdd = currentFactor*output[fourierStart + (k + N/2) - xabsoluteIndex[0]];
			Complex t = output[fourierStart + k - xabsoluteIndex[0]];
			output[fourierStart + k - xabsoluteIndex[0]] = t + factorOdd;
			output[fourierStart + (k+N/2) - xabsoluteIndex[0]] = t - factorOdd;

			currentFactor = currentFactor*factor;
		}
	}
}

void filterHighHarmonic(Complex*** fourierImage, int xnumber, int ynumber, int znumber, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, int* xabsoluteIndex, int* yabsoluteIndex,
                        int* zabsoluteIndex, int cutNumber){
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				double middlex = xnumberGeneral/2.0;
				double middley = ynumberGeneral/2.0;
				double middlez = znumberGeneral/2.0;
				int kx = i + xabsoluteIndex[0];
				int ky = j + yabsoluteIndex[0];
				int kz = k + zabsoluteIndex[0];
				if(xnumberGeneral > 1){
					if(kx < middlex){
						if(kx > cutNumber){
							fourierImage[i][j][k] = Complex(0,0);
						}
					} else {
						if( kx < xnumberGeneral - cutNumber){
							fourierImage[i][j][k] = Complex(0, 0);
						}
					}
				}

				if(ynumberGeneral > 1){
					if(ky < middley){
						if(kx > cutNumber){
							fourierImage[i][j][k] = Complex(0,0);
						}
					} else {
						if( ky < ynumberGeneral - cutNumber){
							fourierImage[i][j][k] = Complex(0, 0);
						}
					}
				}

				if(znumberGeneral > 1){
					if(kz < middlez){
						if(kz > cutNumber){
							fourierImage[i][j][k] = Complex(0,0);
						}
					} else {
						if( kz < znumberGeneral - cutNumber){
							fourierImage[i][j][k] = Complex(0, 0);
						}
					}
				}
			}
		}
	}
}
