#ifndef FOURIER_H
#define FOURIER_H

#include "mpi.h"

class Complex;
void fourierTranslationX(Complex*** input, Complex*** output, bool direct, int xnumber, int ynumber,
                         int znumber, int xnumberGeneral, int* xabsoluteIndex, MPI_Comm cartCommX, int* cartCoord,
                         int* cartDim);
void fourierTranslationY(Complex*** input, Complex*** output, bool direct, int xnumber, int ynumber,
                         int znumber, int ynumberGeneral, int* yabsoluteIndex, MPI_Comm cartCommY, int* cartCoord,
                         int* cartDim);
void fourierTranslationZ(Complex*** input, Complex*** output, bool direct, int xnumber, int ynumber,
                         int znumber, int znumberGeneral, int* zabsoluteIndex, MPI_Comm cartCommZ, int* cartCoord,
                         int* cartDim);

void fourierTranslation(Complex*** input, Complex*** output, Complex*** tempResultX, Complex*** tempResultXY, bool direct, int xnumber, int ynumber, int znumber, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, MPI_Comm& cartComm, MPI_Comm& cartCommX, MPI_Comm& cartCommY, MPI_Comm& cartCommZ, int* xabsoluteIndex, int* yabsoluteIndex, int* zabsoluteIndex, int* cartCoord, int* cartDim);

void filterHighHarmonic(Complex*** fourierImage, int xnumber, int ynumber, int znumber, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, int* xabsoluteIndex, int* yabsoluteIndex,
                        int* zabsoluteIndex, int cutNumber);
#endif