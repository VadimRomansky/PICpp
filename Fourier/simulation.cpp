#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <mpi.h>
#include "math.h"

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

	int dimsYZ[3];
	dimsYZ[0] = 0;
	dimsYZ[1] = 1;
	dimsYZ[2] = 1;
	MPI_Cart_sub(cartComm, dimsYZ, &cartCommYZ);
	int dimsZ[3];
	dimsZ[0] = 0;
	dimsZ[1] = 0;
	dimsZ[2] = 1;
	MPI_Cart_sub(cartComm, dimsZ, &cartCommZ);
	int dimsXY[3];
	dimsXY[0] = 1;
	dimsXY[1] = 1;
	dimsXY[2] = 0;
	MPI_Cart_sub(cartComm, dimsXY, &cartCommXY);
	int dimsY[3];
	dimsY[0] = 0;
	dimsY[1] = 1;
	dimsY[2] = 0;
	MPI_Cart_sub(cartComm, dimsY, &cartCommY);
	int dimsXZ[3];
	dimsXZ[0] = 1;
	dimsXZ[1] = 0;
	dimsXZ[2] = 1;
	MPI_Cart_sub(cartComm, dimsXZ, &cartCommXZ);
	int dimsX[3];
	dimsX[0] = 1;
	dimsX[1] = 0;
	dimsX[2] = 0;
	MPI_Cart_sub(cartComm, dimsX, &cartCommX);

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
			delete[] fourierTempX[i][j];
			delete[] fourierTempXY[i][j];
			delete[] fourierImage[i][j];
			delete[] result[i][j];
		}
		delete[] function[i];
		delete[] fourierTempX[i];
		delete[] fourierTempXY[i];
		delete[] fourierImage[i];
		delete[] result[i];
	}
	delete[] function;
	delete[] fourierTempX;
	delete[] fourierTempXY;
	delete[] fourierImage;
	delete[] result;

	delete[] function1d;
	delete[] fourierImage1d;
	delete[] result1d;
}


void Simulation::createArrays() {
	xgrid = new double[xnumber];
	ygrid = new double[ynumber];
	zgrid = new double[znumber];

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
	fourierTempX = new Complex**[xnumber];
	fourierTempXY = new Complex**[xnumber];
	fourierImage = new Complex**[xnumber];
	result = new Complex**[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		function[i] = new Complex*[ynumber];
		fourierTempX[i] = new Complex*[ynumber];
		fourierTempXY[i] = new Complex*[ynumber];
		fourierImage[i] = new Complex*[ynumber];
		result[i] = new Complex*[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			function[i][j] = new Complex[znumber];
			fourierTempX[i][j] = new Complex[znumber];
			fourierTempXY[i][j] = new Complex[znumber];
			fourierImage[i][j] = new Complex[znumber];
			result[i][j] = new Complex[znumber];
			for(int k = 0; k < znumber; ++k) {
				function[i][j][k] = Complex(0, 0);
				fourierTempX[i][j][k] = Complex(0, 0);
				fourierTempXY[i][j][k] = Complex(0, 0);
				fourierImage[i][j][k] = Complex(0, 0);
				result[i][j][k] = Complex(0, 0);
			}
		}
	}

	function1d = new Complex[xnumber];
	fourierImage1d = new Complex[xnumber];
	result1d = new Complex[xnumber];
	for(int i = 0; i < xnumber; ++i){
		function1d[i] = Complex(0, 0);
		fourierImage1d[i] = Complex(0, 0);
		result1d[i] = Complex(0, 0);
	}
}

void Simulation::setSpaceForProc() {
	int tempXnumber = ((xnumberGeneral) / cartDim[0]);
	int modXnumber = (xnumberGeneral) % cartDim[0];

	if (cartCoord[0] >= cartDim[0] - modXnumber) {
		xnumber = tempXnumber + 1;
		firstAbsoluteXindex = xnumberGeneral - (xnumber) * (cartDim[0] - cartCoord[0]);
	} else {
		xnumber = tempXnumber;
		firstAbsoluteXindex = (xnumber) * cartCoord[0];
	}

	int tempYnumber = ((ynumberGeneral) / cartDim[1]);
	int modYnumber = (ynumberGeneral) % cartDim[1];

	if (cartCoord[1] >= cartDim[1] - modYnumber) {
		ynumber = tempYnumber + 1;
		firstAbsoluteYindex = ynumberGeneral - (ynumber) * (cartDim[1] - cartCoord[1]);
	} else {
		ynumber = tempYnumber;
		firstAbsoluteYindex = (ynumber) * cartCoord[1];
	}

	int tempZnumber = ((znumberGeneral) / cartDim[2]);
	int modZnumber = (znumberGeneral) % cartDim[2];

	if (cartCoord[2] >= cartDim[2] - modZnumber) {
		znumber = tempZnumber + 1;
		firstAbsoluteZindex = znumberGeneral - (znumber) * (cartDim[2] - cartCoord[2]);
	} else {
		znumber = tempZnumber;
		firstAbsoluteZindex = (znumber) * cartCoord[2];
	}

	//todo boundary conditiontype
	/*if (cartDim[0] > 5) {
		int firstSmallRegions = 2 * cartDim[0] / 3;
		int firstSmallXnumber = xnumberGeneral / 3;
		int lastLargeRegions = cartDim[0] - firstSmallRegions;
		int lastLargeXnumber = xnumberGeneral - firstSmallXnumber;
		if (cartCoord[0] < firstSmallRegions) {
			tempXnumber = firstSmallXnumber / firstSmallRegions;
			modXnumber = firstSmallXnumber % firstSmallRegions;
			if (cartCoord[0] >= firstSmallRegions - modXnumber) {
				xnumber = tempXnumber + 1;
				firstAbsoluteXindex = firstSmallXnumber - (xnumber) * (firstSmallRegions - cartCoord[0]) - additionalBinNumber - 1;
			} else {
				xnumber = tempXnumber;
				firstAbsoluteXindex = (xnumber) * cartCoord[0] - additionalBinNumber - 1;
			}
		} else {
			tempXnumber = lastLargeXnumber / lastLargeRegions;
			modXnumber = lastLargeXnumber % lastLargeRegions;
			if (cartCoord[0] >= cartDim[0] - modXnumber) {
				xnumber = tempXnumber + 1;
				firstAbsoluteXindex = xnumberGeneral - (xnumber) * (cartDim[0] - cartCoord[0]) - additionalBinNumber - 1;
			} else {
				xnumber = tempXnumber;
				firstAbsoluteXindex = firstSmallXnumber + (xnumber) * (cartCoord[0] - firstSmallRegions) - additionalBinNumber - 1;
			}
		}

	}*/
	//printf("xnumber = %d\n", xnumber);

	xsize = xnumber * xsizeGeneral / xnumberGeneral;
	ysize = ynumber * ysizeGeneral / ynumberGeneral;
	zsize = znumber * zsizeGeneral / znumberGeneral;

	//deltaX = xsize / (xnumber);
	//deltaY = ysize / (ynumber);
	//deltaZ = zsize / (znumber);

	deltaX = xsizeGeneral / (xnumberGeneral);
	deltaY = ysizeGeneral / (ynumberGeneral);
	deltaZ = zsizeGeneral / (znumberGeneral);


	leftX = xsizeGeneral + (firstAbsoluteXindex) * deltaX;
	rightX = leftX + xsize;
	leftY = ysizeGeneral + (firstAbsoluteYindex) * deltaY;
	rightY = leftY + ysize;
	leftZ = zsizeGeneral + (firstAbsoluteZindex) * deltaZ;
	rightZ = leftZ + zsize;
}

void Simulation::initializeArrays(){
	for(int i = 0; i < xnumber; ++i){
		xabsoluteIndex[i] = firstAbsoluteXindex + i;
	}

	for(int j = 0; j < ynumber; ++j){
		yabsoluteIndex[j] = firstAbsoluteYindex + j;
	}

	for(int k = 0; k < znumber; ++k){
		zabsoluteIndex[k] = firstAbsoluteZindex + k;
	}

	for(int i = 0; i < xnumber; ++i){
		xgrid[i] = xsizeGeneral + (xabsoluteIndex[i]) * deltaX;
	}
	xgrid[xnumber] = xgrid[0] + xnumber*deltaX;

	for(int i = 0; i < ynumber + 1; ++i){
		ygrid[i] = ysizeGeneral + (yabsoluteIndex[i]) * deltaY;
	}
	ygrid[ynumber] = ygrid[0] + ynumber*deltaY;

	for(int i = 0; i < znumber + 1; ++i){
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
	file = fopen((outputDir + "initial1d.dat").c_str(), "w");
	fclose(file);
	file = fopen((outputDir + "fourier1d.dat").c_str(), "w");
	fclose(file);
	file = fopen((outputDir + "final1d.dat").c_str(), "w");
	fclose(file);
}

void Simulation::randomSimulation(){
	srand(time(NULL) + rank);

	for(int i = 0; i < xnumber; ++i){
		function1d[i].re = rand()%randomSeed - (randomSeed/2);
		function1d[i].im = rand()%randomSeed - (randomSeed/2);
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				function[i][j][k].re = rand()%randomSeed - (randomSeed/2.0);
				//function[i][j][k].im = rand()%randomSeed - (randomSeed/2.0);

				//function[i][j][k].re = sin(2*pi*i/xnumberGeneral) + sin(0.35 + 2*pi*i*2/xnumberGeneral) + sin(0.1 + 2*pi*i*8/xnumberGeneral);
				//function[i][j][k].re = 1.0;
				function[i][j][k].im = 0;
			}
		}
	}

	std::string outputDir = outputDirectory;
	outputArray((outputDir + "initial.dat").c_str(), function, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
	fourierTranslation(function, fourierImage, fourierTempX, fourierTempXY, true, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY, cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
	//fourierTranslationX(function, fourierImage, true, xnumber, ynumber, znumber, xnumberGeneral, xabsoluteIndex, cartCommX, cartCoord, cartDim);
	outputArray((outputDir + "fourier.dat").c_str(), fourierImage, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
	filterHighHarmonic(fourierImage, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, xnumberGeneral/16);
	fourierTranslation(fourierImage, result, fourierTempX, fourierTempXY, false, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY, cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
	//fourierTranslationX(fourierImage, result, false, xnumber, ynumber, znumber, xnumberGeneral, xabsoluteIndex, cartCommX, cartCoord, cartDim);
	outputArray((outputDir + "final.dat").c_str(), result, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);

	/*outputArray1d((outputDir + "initial1d.dat").c_str(), function1d, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
	fastFourier1d(function1d, fourierImage1d, true, xnumber, xnumberGeneral, xnumberGeneral, cartComm, cartCoord, cartDim, 0, 1, xabsoluteIndex, 0);
	outputArray1d((outputDir + "fourier1d.dat").c_str(), fourierImage1d, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);
	fastFourier1d(fourierImage1d, result1d, false, xnumber, xnumberGeneral, xnumberGeneral, cartComm, cartCoord, cartDim, 0, 1, xabsoluteIndex, 0);
	outputArray1d((outputDir + "final1d.dat").c_str(), result1d, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCoord, cartDim);*/
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

	fourierTranslation(rightPart, fourierImage, fourierTempX, fourierTempXY, true, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY, cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);
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
	fourierTranslation(fourierImage, result, fourierTempX, fourierTempXY, false, xnumber, ynumber, znumber, xnumberGeneral, ynumberGeneral, znumberGeneral, cartComm, cartCommX, cartCommY, cartCommZ, xabsoluteIndex, yabsoluteIndex, zabsoluteIndex, cartCoord, cartDim);

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
