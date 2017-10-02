#ifndef OUTPUT_H
#define OUTPUT_H

class Complex;

void outputArray(const char* outputFileName, Complex*** array, int xnumber, int ynumber, int znumber, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, MPI_Comm& cartComm, int* cartCoord, int* cartDim);

#endif