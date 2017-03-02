#ifndef CONSTANTS_H
#define CONSTANTS_H

const int randomSeed = 4096;
const int numThreads = 4;

const char* const outputDirectory = "./output/";
const char* const inputDirectory = "./input/";
const char* const backupDirectory = "./backup/";

const double atomicUnionMass = 1.66053904E-24;
const double massProtonReal = 1.67262177E-24;
const double massAlphaReal = 6.644656E-24;
const double massElectronReal = 0.910938291E-27;
const double massElectronFactor = 100;
const double massDeuteriumReal = 3.34449696893E-24;
const double massHelium3Real = 5.00823792874E-24;
const double massOxygen = 26.5601801672E-24;
const double massSilicon = 46.4567787264E-24;
const double kBoltzman = 1.3806488E-16;
const double speed_of_light = 2.99792458E10;
const double electron_charge = 4.803529695E-10;
const double pi = 3.1415926535897932384626433832795028841971693993751;
const double epsilon = 1E-16;
const double maxErrorLevel = 1E-7;
const double particleVelocityErrorLevel = 1E-7;
const double timeEpsilonKourant = 0.5;
const double timeEpsilonFields = 0.25;
const double timeEpsilonPlasma = 0.25;
const double timeEpsilonCapture = 0.25;
const double relativisticPrecision = 0.03;
const double initialTheta = 0.5;
const double smoothingParameter = 0.25;
const double particleSplitLevel = 3;
const double particleSplitAngle = 0.02;
const int splitParticlesParameter = 100;
const int pnumber = 500;
const int writeParameter = 100;
const int writeTrajectoryNumber = 50;
const int writeParticleNumber = 50;
const int divergenceCleanUpParameter = 1;
const int writeBackupParameter = 10000000;
const int maxNewtonIterations = 10;
const int maxGMRESIterations = 50;
const int splineOrder = 0;

const int debugPoint = 7;
  
#endif
