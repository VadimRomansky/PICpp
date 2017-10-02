#include "stdlib.h"
#include "stdio.h"
#include "mpi.h"
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "largeVectorBasis.h"

//remember that we need xnumber + 1 size!!!
LargeVectorBasis::LargeVectorBasis(int capacityv, int xnumberv, int ynumberv, int znumberv, int lnumberv) {
	size = 0;
	xnumber = xnumberv;
	ynumber = ynumberv;
	znumber = znumberv;
	lnumber = lnumberv;
	capacity = capacityv;

	array = new double**** [capacity];
	for(int m = 0; m < capacity; ++m) {
		array[m] = new double*** [xnumber];
		for(int i = 0; i < xnumber; ++i) {
			array[m][i] = new double** [ynumber];
			for(int j = 0; j < ynumber; ++j) {
				array[m][i][j] = new double* [znumber];
				for(int k = 0; k < znumber; ++k) {
					array[m][i][j][k] = new double[lnumber];
					for(int l = 0; l < lnumber; ++l) {
						array[m][i][j][k][l] = 0;
					}
				}
			}
		}
	}
}

LargeVectorBasis::~LargeVectorBasis() {
	for(int m = 0; m < capacity; ++m) {
		for(int i = 0; i < xnumber; ++i) {
			for(int j = 0; j < ynumber; ++j) {
				for(int k = 0; k < znumber; ++k) {
					delete[] array[m][i][j][k];
				}
				delete[] array[m][i][j];
			}
			delete[] array[m][i];
		}
		delete[] array[m];
	}	
	delete[] array;
}

void LargeVectorBasis::resize(int capacityv) {
	if(capacityv > capacity) {
		double***** newArray = new double**** [capacityv];
		for(int m = 0; m < capacity; ++m) {
			newArray[m] = array[m];
		}
		for(int m = capacity; m < capacityv; ++m) {
			newArray[m] = new double*** [xnumber];
			for(int i = 0; i < xnumber; ++i) {
				newArray[m][i] = new double** [ynumber];
				for(int j = 0; j < ynumber; ++j) {
					newArray[m][i][j] = new double* [znumber];
					for(int k = 0; k < znumber; ++k) {
						newArray[m][i][j][k] = new double[lnumber];
						for(int l = 0; l < lnumber; ++l) {
							newArray[m][i][j][k][l] = 0;
						}
					}
				}
			}
		}

		capacity = capacityv;
		delete[] array;
		array = newArray;
		return;
	}
	if(capacityv >= size) {
		double***** newArray = new double**** [capacityv];
		for(int m = 0; m < capacityv; ++m) {
			newArray[m] = array[m];
		}
		for(int m = capacityv; m < capacity; ++m) {
			for(int i = 0; i < xnumber; ++i) {
				for(int j = 0; j < ynumber; ++j) {
					for(int k = 0; k < znumber; ++k) {
						delete[] array[m][i][j][k];
					}
					delete[] array[m][i][j];
				}
				delete[] array[m][i];
			}
			delete[] array[m];
		}
		delete[] array;
		array = newArray;
		capacity = capacityv;
		return;
	}
	printf("capacity < size\n");
	MPI_Finalize();
	exit(0);
}

void LargeVectorBasis::clear() {
	size = 0;
	for(int m = 0; m < capacity; ++m) {
		for(int i = 0; i < xnumber; ++i) {
			for(int j = 0; j < ynumber; ++j) {
				for(int k = 0; k < znumber; ++k) {
					for(int l = 0; l < lnumber; ++l) {
						array[m][i][j][k][l] = 0;
					}
				}
			}
		}
	}
}
