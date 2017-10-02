#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <mpi.h>

#include "simulation.h"
#include "complex.h"
#include "output.h"
#include "fourier.h"

Simulation::Simulation(MPI_Comm& comm) {

	cartComm = comm;
	int periods[MPI_dim];
	MPI_Cart_get(cartComm, MPI_dim, cartDim, periods, cartCoord);
	MPI_Comm_rank(cartComm, &rank);
	MPI_Comm_size(cartComm, &nprocs);

	int tempCoord[3];
	for (int i = 0; i < 3; ++i) {
		tempCoord[i] = cartCoord[i];
	}
	tempCoord[0] -= 1;
	if (tempCoord[0] < 0) {
		tempCoord[0] = cartDim[0] - 1;
	}
	MPI_Cart_rank(cartComm, tempCoord, &leftRank);
	tempCoord[0] = cartCoord[0] + 1;
	if (tempCoord[0] >= cartDim[0]) {
		tempCoord[0] = 0;
	}
	MPI_Cart_rank(cartComm, tempCoord, &rightRank);
	tempCoord[0] = cartCoord[0];
	tempCoord[1] -= 1;
	if (tempCoord[1] < 0) {
		tempCoord[1] = cartDim[1] - 1;
	}
	MPI_Cart_rank(cartComm, tempCoord, &frontRank);
	tempCoord[1] = cartCoord[1] + 1;
	if (tempCoord[1] >= cartDim[1]) {
		tempCoord[1] = 0;
	}
	MPI_Cart_rank(cartComm, tempCoord, &backRank);
	tempCoord[1] = cartCoord[1];
	tempCoord[2] -= 1;
	if (tempCoord[2] < 0) {
		tempCoord[2] = cartDim[2] - 1;
	}
	MPI_Cart_rank(cartComm, tempCoord, &bottomRank);
	tempCoord[2] = cartCoord[2] + 1;
	if (tempCoord[2] >= cartDim[2]) {
		tempCoord[2] = 0;
	}
	MPI_Cart_rank(cartComm, tempCoord, &topRank);
	setSpaceForProc();
	createArrays();
	initializeArrays();
	
	if(rank == 0){
		createFiles();
	}
}

Simulation::~Simulation(){
	delete[] xgrid;
	delete[] ygrid;
	delete[] zgrid;

	delete[] middleXgrid;
	delete[] middleYgrid;
	delete[] middleZgrid;

	delete[] kxgrid;
	delete[] kygrid;
	delete[] kzgrid;

	delete[] xabsoluteIndex;
	delete[] yabsoluteIndex;
	delete[] zabsoluteIndex;

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			delete[] function[i][j];
			delete[] fourierImage[i][j];
			delete[] result[i][j];
		}
		delete[] function[i];
		delete[] fourierImage[i];
		delete[] result[i];
	}
	delete[] function;
	delete[] fourierImage;
	delete[] result;
}


void Simulation::createArrays() {
	xgrid = new double[xnumber + 1];
	ygrid = new double[ynumber + 1];
	zgrid = new double[znumber + 1];

	middleXgrid = new double[xnumber];
	middleYgrid = new double[ynumber];
	middleZgrid = new double[znumber];

	kxgrid = new double[xnumber];
	kygrid = new double[ynumber];
	kzgrid = new double[znumber];

	xabsoluteIndex = new int[xnumber];
	yabsoluteIndex = new int[ynumber];
	zabsoluteIndex = new int[znumber];

	function = new Complex**[xnumber];
	fourierImage = new Complex**[xnumber];
	result = new Complex**[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		function[i] = new Complex*[ynumber];
		fourierImage[i] = new Complex*[ynumber];
		result[i] = new Complex*[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			function[i][j] = new Complex[znumber];
			fourierImage[i][j] = new Complex[znumber];
			result[i][j] = new Complex[znumber];
			for(int k = 0; k < znumber; ++k) {
				function[i][j][k] = Complex(0, 0);
				fourierImage[i][j][k] = Complex(0, 0);
				result[i][j][k] = Complex(0, 0);
			}
		}
	}
}

void Simulation::setSpaceForProc() {
	int tempXnumber = ((xnumberGeneral) / cartDim[0]);
	int modXnumber = (xnumberGeneral) % cartDim[0];

	if (cartCoord[0] >= cartDim[0] - modXnumber) {
		xnumber = tempXnumber + 1;
		firstXabsoluteIndex = xnumberGeneral - (xnumber) * (cartDim[0] - cartCoord[0]);
	} else {
		xnumber = tempXnumber;
		firstXabsoluteIndex = (xnumber) * cartCoord[0];
	}

	int tempYnumber = ((ynumberGeneral) / cartDim[1]);
	int modYnumber = (ynumberGeneral) % cartDim[1];

	if (cartCoord[1] >= cartDim[1] - modYnumber) {
		ynumber = tempYnumber + 1;
		firstYabsoluteIndex = ynumberGeneral - (ynumber) * (cartDim[1] - cartCoord[1]);
	} else {
		ynumber = tempYnumber;
		firstYabsoluteIndex = (ynumber) * cartCoord[1];
	}

	int tempZnumber = ((znumberGeneral) / cartDim[2]);
	int modZnumber = (znumberGeneral) % cartDim[2];

	if (cartCoord[2] >= cartDim[2] - modZnumber) {
		znumber = tempZnumber + 1;
		firstZabsoluteIndex = znumberGeneral - (znumber) * (cartDim[2] - cartCoord[2]);
	} else {
		znumber = tempZnumber;
		firstZabsoluteIndex = (znumber) * cartCoord[2];
	}

	deltaX = xsizeGeneral / (xnumberGeneral);
	deltaY = ysizeGeneral / (ynumberGeneral);
	deltaZ = zsizeGeneral / (znumberGeneral);

	deltaKx = 2*pi/xsizeGeneral;
	deltaKy = 2*pi/ysizeGeneral;
	deltaKz = 2*pi/zsizeGeneral;
}

void Simulation::initializeArrays(){
	for(int i = 0; i < xnumber; ++i){
		xabsoluteIndex[i] = firstXabsoluteIndex + i;
	}

	for(int j = 0; j < ynumber; ++j){
		yabsoluteIndex[j] = firstYabsoluteIndex + j;
	}

	for(int k = 0; k < znumber; ++k){
		zabsoluteIndex[k] = firstZabsoluteIndex + k;
	}

	for(int i = 0; i < xnumber; ++i){
		xgrid[i] = xsizeGeneral + (xabsoluteIndex[i]) * deltaX;
	}
	xgrid[xnumber] = xgrid[0] + xnumber*deltaX;

	for(int i = 0; i < ynumber; ++i){
		ygrid[i] = ysizeGeneral + (yabsoluteIndex[i]) * deltaY;
	}
	ygrid[ynumber] = ygrid[0] + ynumber*deltaY;

	for(int i = 0; i < znumber; ++i){
		zgrid[i] = zsizeGeneral + (zabsoluteIndex[i]) * deltaZ;
	}
	zgrid[znumber] = zgrid[0] + znumber*deltaZ;

	for(int i = 0; i < xnumber; ++i){
		middleXgrid[i] = (xgrid[i] + xgrid[i+1])/2;
		kxgrid[i] = xabsoluteIndex[i]*deltaKx;
	}

	for(int i = 0; i < ynumber; ++i){
		middleYgrid[i] = (ygrid[i] + ygrid[i+1])/2;
		kygrid[i] = yabsoluteIndex[i]*deltaKy;
	}

	for(int i = 0; i < znumber; ++i){
		middleZgrid[i] = (zgrid[i] + zgrid[i+1])/2;
		kzgrid[i] = zabsoluteIndex[i]*deltaKz;
	}
}

void Simulation::createFiles(){
	std::string outputDir = outputDirectory;
	FILE* file = fopen((outputDir + "initial.dat").c_str(), "w");
	fclose(file);
	file = fopen((outputDir + "fourier.dat").c_str(), "w");
	fclose(file);
	file = fopen((outputDir + "final.dat").c_str(), "w");
	fclose(file);
}

void Simulation::randomSimulation(){
	srand(time(NULL) + rank);

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				function[i][j][k].re = rand()%randomSeed - (randomSeed/2);
				function[i][j][k].im = rand()%randomSeed - (randomSeed/2);
			}
		}
	}

	std::string outputDir = outputDirectory;
	outputArray((outputDir + "initial.dat").c_str(), function, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
	fourierTranslation(function, fourierImage, true, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
	outputArray((outputDir + "fourier.dat").c_str(), fourierImage, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
	fourierTranslation(fourierImage, result, false, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
	outputArray((outputDir + "final.dat").c_str(), result, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
}

void Simulation::poissonSolving() {
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				//function[i][j][k].re = sin(2*pi*xabsoluteIndex[i]/xnumberGeneral) + cos(2*2*pi*xabsoluteIndex[i]/xnumberGeneral) + sin(2*2*pi*yabsoluteIndex[j]/ynumberGeneral)*cos(2*pi*zabsoluteIndex[k]/znumberGeneral);
				function[i][j][k].re = sin(2*2*pi*yabsoluteIndex[j]/ynumberGeneral)*cos(2*pi*zabsoluteIndex[k]/znumberGeneral);
				//function[i][j][k].re = sin(2*pi*xabsoluteIndex[i]/xnumberGeneral);
				//function[i][j][k].re = sin(2*pi*yabsoluteIndex[j]/ynumberGeneral);
				//function[i][j][k].re = sin(2*pi*zabsoluteIndex[k]/znumberGeneral);
				function[i][j][k].im = 0;
			}
		}
	}

	Complex*** rightPart = create3dArray();
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				//rightPart[i][j][k].re = -4*pi*pi*((sin(2*pi*xabsoluteIndex[i]/xnumberGeneral) + 2*2*cos(2*2*pi*xabsoluteIndex[i]/xnumberGeneral))/(xnumberGeneral*xnumberGeneral) + (2*2/(ynumberGeneral*ynumberGeneral) + 1/(znumberGeneral*znumberGeneral))*sin(2*2*pi*yabsoluteIndex[j]/ynumberGeneral)*cos(2*pi*zabsoluteIndex[k]/znumberGeneral));
				rightPart[i][j][k].re = -4*pi*pi*((2.0*2.0/(ynumberGeneral*ynumberGeneral) + 1.0/(znumberGeneral*znumberGeneral))*sin(2*2*pi*yabsoluteIndex[j]/ynumberGeneral)*cos(2*pi*zabsoluteIndex[k]/znumberGeneral));
				//rightPart[i][j][k].re = -4*pi*pi*sin(2*pi*xabsoluteIndex[i]/xnumberGeneral)/(xnumberGeneral*xnumberGeneral);
				//rightPart[i][j][k].re = -4*pi*pi*sin(2*pi*yabsoluteIndex[j]/ynumberGeneral)/(ynumberGeneral*ynumberGeneral);
				//rightPart[i][j][k].re = -4*pi*pi*sin(2*pi*zabsoluteIndex[k]/znumberGeneral)/(znumberGeneral*znumberGeneral);
				rightPart[i][j][k].im = 0;
			}
		}
	}

	fourierTranslation(rightPart, fourierImage, true, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
	int halfx = xnumberGeneral/2;
	int halfy = ynumberGeneral/2;
	int halfz = znumberGeneral/2;

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				double kx = xabsoluteIndex[i];
				if(kx > halfx) {
					kx = xnumberGeneral - kx;
				}
				kx *= 2*pi/xnumberGeneral;
				double ky = yabsoluteIndex[j];
				if(ky > halfy) {
					ky = ynumberGeneral - ky;
				}
				ky *= 2*pi/ynumberGeneral;
				double kz = zabsoluteIndex[k];
				if(kz > halfz) {
					kz = znumberGeneral - kz;
				}
				kz *= 2*pi/znumberGeneral;
				if(xabsoluteIndex[i] != 0 || yabsoluteIndex[j] != 0 || zabsoluteIndex[k] != 0){
					fourierImage[i][j][k] = fourierImage[i][j][k]/(-kx*kx - ky*ky - kz*kz);
				} else {
					fourierImage[i][j][k] = Complex(0, 0);
				}
			}
		}
	}

	/*for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				double kx = xabsoluteIndex[i];
				kx *= 2*pi/xnumberGeneral;
				double ky = yabsoluteIndex[j];
				ky *= 2*pi/ynumberGeneral;
				double kz = zabsoluteIndex[k];
				kz *= 2*pi/znumberGeneral;
				if(xabsoluteIndex[i] != 0 || yabsoluteIndex[j] != 0 || zabsoluteIndex[k] != 0){
					fourierImage[i][j][k] = fourierImage[i][j][k]/(-kx*kx - ky*ky - kz*kz);
				} else {
					fourierImage[i][j][k] = Complex(0, 0);
				}
			}
		}
	}*/
	fourierTranslation(fourierImage, result, false, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);

	std::string outputDir = outputDirectory;

	outputArray((outputDir + "initial.dat").c_str(), function, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
	outputArray((outputDir + "fourier.dat").c_str(), fourierImage, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
	outputArray((outputDir + "final.dat").c_str(), result, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);

	delete3dArray(rightPart);
}

Complex*** Simulation::create3dArray() {
	Complex*** array = new Complex**[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		array[i] = new Complex*[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			array[i][j] = new Complex[znumber];
			for(int k = 0; k < znumber; ++k) {
				array[i][j][k] = Complex(0, 0);
			}
		}
	}
	return array;
}

void Simulation::delete3dArray(Complex*** array) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			delete[] array[i][j];
		}
		delete[] array[i];
	}
	delete[] array;
}
