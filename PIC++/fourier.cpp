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
				(phaseFactor * 2 * pi * xabsoluteIndex[1 + additionalBinNumber] * knumber) / xnumberGeneral);
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
			if((knumber >= xabsoluteIndex[1+additionalBinNumber]) && (knumber < xabsoluteIndex[1+additionalBinNumber] + xnumber)){
				for (int j = 0; j < ynumber; ++j) {
					for (int k = 0; k < znumber; ++k) {
						output[knumber-xabsoluteIndex[1+additionalBinNumber]][j][k] = Complex(sumAll[l], sumAll[l+1]);
						if (!direct) {
							output[knumber-xabsoluteIndex[1+additionalBinNumber]][j][k] = output[knumber-xabsoluteIndex[1+additionalBinNumber]][j][k] / xnumberGeneral;
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
				(phaseFactor * 2 * pi * yabsoluteIndex[1 + additionalBinNumber] * knumber) / ynumberGeneral);
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
			if((knumber >= yabsoluteIndex[1+additionalBinNumber]) && (knumber < yabsoluteIndex[1+additionalBinNumber] + ynumber)){
				for (int i = 0; i < xnumber; ++i) {
					for (int k = 0; k < znumber; ++k) {
						output[i][knumber-yabsoluteIndex[1+additionalBinNumber]][k] = Complex(sumAll[l], sumAll[l+1]);
						if (!direct) {
							output[i][knumber-yabsoluteIndex[1+additionalBinNumber]][k] = output[i][knumber-yabsoluteIndex[1+additionalBinNumber]][k] / ynumberGeneral;
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
				(phaseFactor * 2 * pi * zabsoluteIndex[1 + additionalBinNumber] * knumber) / znumberGeneral);
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
			MPI_Allreduce(sum, sumAll, 2*xnumber*ynumber, MPI_DOUBLE, MPI_SUM, cartCommZ);
			l=0;
			if((knumber >= zabsoluteIndex[1+additionalBinNumber]) && (knumber < zabsoluteIndex[1+additionalBinNumber] + znumber)){
				for (int i = 0; i < xnumber; ++i) {
					for (int j = 0; j < ynumber; ++j) {
						output[i][j][knumber-zabsoluteIndex[1+additionalBinNumber]] = Complex(sumAll[l], sumAll[l+1]);
						if (!direct) {
							output[i][j][knumber-zabsoluteIndex[1+additionalBinNumber]] = output[i][j][knumber-zabsoluteIndex[1+additionalBinNumber]] / znumberGeneral;
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
                        bool direct, int xnumberAdded, int ynumberAdded, int znumberAdded, int xnumberGeneral,
                        int ynumberGeneral, int znumberGeneral, MPI_Comm& cartComm, MPI_Comm& cartCommX,
                        MPI_Comm& cartCommY, MPI_Comm& cartCommZ, int* xabsoluteIndex, int* yabsoluteIndex,
                        int* zabsoluteIndex, int* cartCoord, int* cartDim) {

	fourierTranslationX(input, tempResultX, direct, xnumberAdded, ynumberAdded, znumberAdded, xnumberGeneral,
	                    xabsoluteIndex, cartCommX, cartCoord, cartDim);

	fourierTranslationY(tempResultX, tempResultXY, direct, xnumberAdded, ynumberAdded, znumberAdded, ynumberGeneral,
	                    yabsoluteIndex, cartCommY, cartCoord, cartDim);

	fourierTranslationZ(tempResultXY, output, direct, xnumberAdded, ynumberAdded, znumberAdded, znumberGeneral,
	                    zabsoluteIndex, cartCommZ, cartCoord, cartDim);
}

void filterHighHarmonic(Complex*** fourierImage, int xnumber, int ynumber, int znumber, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, int* xabsoluteIndex, int* yabsoluteIndex,
                        int* zabsoluteIndex, int cutNumber){
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				double middlex = xnumberGeneral/2.0;
				double middley = ynumberGeneral/2.0;
				double middlez = znumberGeneral/2.0;
				int kxnumber = i + xabsoluteIndex[1+additionalBinNumber];
				int kynumber = j + yabsoluteIndex[1+additionalBinNumber];
				int kznumber = k + zabsoluteIndex[1+additionalBinNumber];
				if(xnumberGeneral > 1){
					if(kxnumber < middlex){
						if(kxnumber > xnumberGeneral/cutNumber){
							fourierImage[i][j][k] = Complex(0,0);
						}
					} else {
						if( xnumberGeneral - kxnumber > xnumberGeneral/cutNumber){
							fourierImage[i][j][k] = Complex(0, 0);
						}
					}
				}

				if(ynumberGeneral > 1){
					if(kynumber < middley){
						if(kynumber > ynumberGeneral/cutNumber){
							fourierImage[i][j][k] = Complex(0,0);
						}
					} else {
						if( ynumberGeneral - kynumber > ynumberGeneral/cutNumber){
							fourierImage[i][j][k] = Complex(0, 0);
						}
					}
				}

				if(znumberGeneral > 1){
					if(kznumber < middlez){
						if(kznumber > znumberGeneral/cutNumber){
							fourierImage[i][j][k] = Complex(0,0);
						}
					} else {
						if( znumberGeneral - kznumber > znumberGeneral/cutNumber){
							fourierImage[i][j][k] = Complex(0, 0);
						}
					}
				}
			}
		}
	}
}