#include "creating_arrays.h"
#include "vector3d.h"

double*** create3array(int xnumber, int ynumber, int znumber) {
	double*** array = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		array[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			array[i][j] = new double[znumber];
			for(int k = 0; k < znumber; ++k) {
				array[i][j][k] = 0;
			}
		}
	}
	return array;
}
double**** create4array(int xnumber, int ynumber, int znumber, int lnumber) {
	double**** array = new double***[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		array[i] = new double**[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			array[i][j] = new double*[znumber];
			for(int k = 0; k < znumber; ++k) {
				array[i][j][k] = new double[lnumber];
				for(int l = 0; l < lnumber; ++l) {
					array[i][j][k][l] = 0;
				}
			}
		}
	}
	return array;
}
Vector3d*** create3vectorArray(int xnumber, int ynumber, int znumber) {
	Vector3d*** array = new Vector3d**[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		array[i] = new Vector3d*[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			array[i][j] = new Vector3d[znumber];
			for(int k = 0; k < znumber; ++k) {
				array[i][j][k] = Vector3d(0, 0, 0);
			}
		}
	}
	return array;
}
Vector3d**** create4vectorArray(int xnumber, int ynumber, int znumber, int lnumber) {
	Vector3d**** array = new Vector3d***[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		array[i] = new Vector3d**[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			array[i][j] = new Vector3d*[znumber];
			for(int k = 0; k < znumber; ++k) {
				array[i][j][k] = new Vector3d[lnumber];
				for(int l = 0; l < lnumber; ++l) {
					array[i][j][k][l] = Vector3d(0, 0, 0);
				}
			}
		}
	}
	return array;
}
Matrix3d*** create3matrixArray(int xnumber, int ynumber, int znumber) {
	Matrix3d*** array = new Matrix3d**[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		array[i] = new Matrix3d*[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			array[i][j] = new Matrix3d[znumber];
			for(int k = 0; k < znumber; ++k) {
				array[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
			}
		}
	}
	return array;
}
Matrix3d**** create4matrixArray(int xnumber, int ynumber, int znumber, int lnumber) {
	Matrix3d**** array = new Matrix3d***[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		array[i] = new Matrix3d**[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			array[i][j] = new Matrix3d*[znumber];
			for(int k = 0; k < znumber; ++k) {
				array[i][j][k] = new Matrix3d[lnumber];
				for(int l = 0; l < lnumber; ++l) {
					array[i][j][k][l] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				}
			}
		}
	}
	return array;
}

void delete3array(double*** array, int xnumber, int ynumber, int znumber) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			delete[] array[i][j];
		}
		delete[] array[i];
	}
	delete[] array;
}
void delete4array(double**** array, int xnumber, int ynumber, int znumber, int lnumber) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k){
				delete[] array[i][j][k];
			}
			delete[] array[i][j];
		}
		delete[] array[i];
	}
	delete[] array;
}
void delete3vectorArray(Vector3d*** array, int xnumber, int ynumber, int znumber) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			delete[] array[i][j];
		}
		delete[] array[i];
	}
	delete[] array;
}
void delete4vectorArray(Vector3d**** array, int xnumber, int ynumber, int znumber, int lnumber) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k){
				delete[] array[i][j][k];
			}
			delete[] array[i][j];
		}
		delete[] array[i];
	}
	delete[] array;
}
void delete3matrixArray(Matrix3d*** array, int xnumber, int ynumber, int znumber) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			delete[] array[i][j];
		}
		delete[] array[i];
	}
	delete[] array;
}
void delete4matrixArray(Matrix3d**** array, int xnumber, int ynumber, int znumber, int lnumber) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k){
				delete[] array[i][j][k];
			}
			delete[] array[i][j];
		}
		delete[] array[i];
	}
	delete[] array;
}