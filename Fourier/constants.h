#ifndef CONSTANTS_H
#define CONSTANTS_H

//for plytech hpc
//const char* const outputDirectory  = "/home/ipntsr/romansky/PIC3d1/output/";

//for MPI on mvs
//const char* const outputDirectory = "/home2/ioffe2/romansky/PIC3d1/output/";

//for MPI on h1 or h2
//const char* const outputDirectory = "/home/vadim/PIC3d1/output/";

//for MPI on alpha
//const char* const outputDirectory = "/home/uikam/romansky/PIC3d3/output/";

//for MPI on Ubuntu
//const char* const outputDirectory = "/home/vadim/PIC++/trunk/PIC++/output/";
 
//for windows 32
//char* const outputDirectory = "C:/users/Vadik/Documents/Visual Studio 2010/Projects/v4/PIC++/output/";

//for windows 64 home
//const char* const outputDirectory = "C:/users/Vadik/Documents/Visual Studio 2010/Projects/PICpp/Fourier/output/";

//for windows 64
const char* const outputDirectory = "C:/users/Vadim/Documents/Visual Studio 2010/Projects/PICpp/Fourier/output/";

const int randomSeed = 1024;
const double pi = 3.1415926535897932384626433832795028841971693993751;

const int xnumberGeneral = 8;
const int ynumberGeneral = 1;
const int znumberGeneral = 1;

const double xsizeGeneral = 8;
const double ysizeGeneral = 1;
const double zsizeGeneral = 1;

const int additionalBinNumber = 1;

const int MPI_nx = 1;
const int MPI_ny = 1;
const int MPI_nz = 1;
const int MPI_dim = 3;

#endif