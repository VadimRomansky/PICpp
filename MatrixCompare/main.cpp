#include "stdio.h"
#include <stdlib.h>

int main(){
	int xnumber = 25;
	int lnumber = 3;
	int number = xnumber*lnumber;

	FILE* matrixFile3d = fopen("./../PIC++/output/maxwellMatrixFile.dat","r");
	FILE* matrixFile1d = fopen("./../PIC++1d/output/maxwellMatrixFile.dat","r");
	FILE* output = fopen("./output/output.dat","w");

	double** matrix3d = new double*[number];
	double** matrix1d = new double*[number];
	for(int i = 0; i < number; ++i){
		matrix3d[i] = new double[number];
		matrix1d[i] = new double[number];
	}

	for(int i = 0; i < number; ++i){
		for(int j = 0; j < number; ++j){
			double a;
			double b;
			fscanf(matrixFile3d, "%lf", &a);
			fscanf(matrixFile1d, "%lf", &b);
			matrix3d[i][j] = a;
			matrix1d[i][j] = b;
			if(a - b > 0){
				printf("aa\n");
			}
			fprintf(output, "%15.10g", a - b);
		}
		fprintf(output, "\n");
	}

	for(int i = 0; i < number; ++i){
		delete[] matrix3d[i];
		delete[] matrix1d[i];
	}
	delete[] matrix3d;
	delete[] matrix1d;

	fclose(matrixFile3d);
	fclose(matrixFile1d);
	fclose(output);

	return 0;
}