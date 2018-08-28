#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "math.h"
#include "stdio.h"

const int initialRandom = 0;
const int randomParameter = 4096;

//const int MPI_nx = 1;
//const int MPI_ny = 1;
//const int MPI_nz = 1;
const int MPI_dim = 3;

const double atomicUnionMass = 1.66053904E-24;
const double massProtonReal = 1.67262177E-24;
const double massAlphaReal = 6.644656E-24;
const double massElectronReal = 0.910938291E-27;
//const double massElectronFactor = massProtonReal/(10.0*massElectronReal);
const double massElectronFactor = 100;
const double massDeuteriumReal = 3.34449696893E-24;
const double massHelium3Real = 5.00823792874E-24;
const double massOxygenReal = 26.5601801672E-24;
const double massSiliconReal = 46.4567787264E-24;
const double kBoltzman = 1.3806488E-16;
const double speed_of_light = 2.99792458E10;
const double electron_charge = 4.803529695E-10;
const double pi = 4*atan2(1.0,1.0);
const double four_pi = 4*pi;
const double epsilon = 1E-16;
const double maxErrorLevel = 1E-8;
const double maxCleanupErrorLevel = 1E-8;
const double particleVelocityErrorLevel = 1E-14;
const double particleCoordinatesErrorLevel = 1E-8;
const double timeEpsilonKourant = 0.45;
const double timeEpsilonFields = 0.5;
const double timeEpsilonPlasma = 0.5;
const double timeEpsilonCapture = 0.5;
const double relativisticPrecision = 0.0001;
const double betaGodfrey = 0.25*(1.0 - timeEpsilonKourant*timeEpsilonKourant);
const double initialTheta =  0.5;
//const double initialTheta =  sqrt(betaGodfrey/timeEpsilonKourant);
const double particleSplitLevel = 3;
const double particleSplitAngle = 0.02;
const double speed_of_light_correction = 1.05;
const double speed_of_light_correction_sqr = speed_of_light_correction*speed_of_light_correction;
const int splitParticlesParameter = 100;
const int pnumber = 500;
const int writeMemoryParameter = 20000;
const int divergenceCleanUpParameter = 1;
const int filteringParameter = 5;
const int frontHalfWidth = 100;
const int writeBackupParameter = 10000000;
const int maxNewtonIterations = 10;
const int maxGMRESIterations = 50;
const int maxDivergenceCleanupIterations = 100;
const int maxSimpleIterationSolverIterations = 10000;
const int reduceStepX = 5;
const int reduceStepY = 1;
const int reduceStepZ = 1;
const int splineOrder = 0;
const int additionalBinNumber = 1 + (splineOrder/2);
const int defaulResistiveLayerWidth = 20;

const double accelerationLevel = 5;
const int mostFastParticlesNumber = 20;


const int debugPoint = 7;
  
#endif
