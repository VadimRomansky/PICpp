#ifndef CONSTANTS_H
#define CONSTANTS_H

const int randomSeed = 1024;
const int numThreads = 4;

//actually, i don't know. sometimes clion prefer one, somitimes another
//for CLion
//const char* const outputDirectory = "../output/";
//const char* const inputDirectory = "../input/";
//const char* const backupDirectory = "../backup/";

//for visual studio and manualy build on cluster
const char* const outputDirectory = "output/";
const char* const inputDirectory = "input/";
const char* const backupDirectory = "backup/";

const double massProtonReal = 1.67262177E-24;
const double massAlphaReal = 6.644656E-24;
const double massElectronReal = 0.910938291E-27;
const double massDeuteriumReal = 3.34449696893E-24;
const double massHelium3Real = 5.00823792874E-24;
const double kBoltzman = 1.3806488E-16;
const double speed_of_light = 2.99792458E10;
const double electron_charge = 4.803529695E-10;
const double pi = 3.1415926535897932384626433832795028841971693993751;
const double epsilon = 1E-16;
const double maxErrorLevel = 1E-9;
const double timeEpsilon = 0.05;
const double relativisticPrecision = 0.001;
const double initialTheta = 0.5;
const int pnumber = 1000;
const int writeParameter = 1;
const int writeBackupParameter = 10000000;
const int maxNewtonIterations = 10;
const int maxGMRESIterations = 100;
const int splineOrder = 0;

const int debugPoint = 7;

#endif
