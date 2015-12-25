#ifndef CONSTANTS_H
#define CONSTANTS_H

const int randomSeed = 1024;
const int numThreads = 6;

const double massProtonReal = 1.67262177E-24;
const double massAlphaReal = 6.644656E-24;
const double massElectronReal = 0.910938291E-27;
//const double massElectron = massProtonReal;
const double massElectron = 100*massElectronReal;
//const double massProton = massElectronReal*200;
const double massProton = massProtonReal;
const double massAlpha = massAlphaReal;
const double kBoltzman = 1.3806488E-16;
const double speed_of_light = 2.99792458E10;
const double electron_charge = 4.803529695E-10;
const double pi = 3.1415926535897932384626433832795028841971693993751;
const double epsilon = 1E-16;
const double maxErrorLevel = 1E-8;
const double timeEpsilon = 0.05;
const double relativisticPrecision = 0.00001;
const double initialTheta = 0.6;
const int pnumber = 1000;
const int writeParameter = 500;
const int writeBackupParameter = 100000;
const int maxNewtonIterations = 100;
const int maxGMRESIterations = 200;
const int splineOrder = 2;

#endif
