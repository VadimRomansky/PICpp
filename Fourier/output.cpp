#include "stdio.h"
#include <stdlib.h>
#include <time.h>
//#include <crtdbg.h>
#include "mpi.h"

#include "output.h"
#include "complex.h"

void outputArray(const char* outputFileName, Complex*** array, int xnumber, int ynumber, int znumber, int xnumberGeneral, int ynumberGeneral, int znumberGeneral, MPI_Comm& cartComm, int* cartCoord, int* cartDim){
	int dims1[3];
	dims1[0] = 0;
	dims1[1] = 1;
	dims1[2] = 1;
	MPI_Comm subCommX;
	MPI_Cart_sub(cartComm, dims1, &subCommX);
	int dims2[3];
	dims2[0] = 0;
	dims2[1] = 0;
	dims2[2] = 1;
	MPI_Comm subCommY;
	MPI_Cart_sub(cartComm, dims2, &subCommY);
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			for(int i = 0; i < xnumber; ++i){
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					MPI_Barrier(subCommX);
					if (cartJ == cartCoord[1]) {
						for(int j = 0; j < ynumber; ++j){
							for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
								MPI_Barrier(subCommY);
								if (cartK == cartCoord[2]) {
									FILE* outFile = fopen(outputFileName, "a");
									for(int k = 0; k < znumber; ++k){
										fprintf(outFile, "%g %g\n", array[i][j][k].re, array[i][j][k].im);
									}
									fclose(outFile);
								}
							}
						}
					}
				}
			}
		}
	}
}