#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "vector"
#include "stdio.h"
#include "particle.h"

class Simulation;
class Vector3d;
class Particle;
class Matrix3d;
class MatrixElement;

void outputDistribution(const char *outFileName, std::vector<Particle *> &particles, int particleType, double plasma_period,
                        double gyroradius, int verbosity, bool newFile);
void outputDistributionShiftedSystem(const char* outFileName, std::vector<Particle *>& particles, Vector3d& shiftV, double& speed_of_light_normalized, int particleType, double gyroradius,
                                     double plasma_period, int verbosity, bool newFile);
//void outputAnisotropy(const char *outFileName, Simulation *simulation, int particleType, double gyroradius,
//                      double plasma_period);
void outputTrajectoryByNumber(const char* outFileName, int number, const Simulation* simulation);
void outputTrajectory(const char *outFileName, Particle *particle, double time, double plasma_period, double gyroradius);
void outputParticlesTrajectories(const char *outFileName, const char* electronOutFileName, std::vector<Particle*> particles, int** numbers, int size, double time, double plasma_period, double scaleFactor, Simulation* simulation);
void outputGridSimple(const char *outFileName, double *grid, int number, bool newFile, double scale = 1.0);
void outputGridSimpleReduced(const char *outFileName, double *grid, int number, int step, bool newFile, double scale = 1.0);

void outputGridX(const char *outFileName, double *grid, int number, int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool
                 newFile, double scale = 1.0);
void outputGridY(const char *outFileName, double *grid, int number, int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool
                 newFile, double scale = 1.0);
void outputGridZ(const char *outFileName, double *grid, int number, int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool
                 newFile, double scale = 1.0);
void outputGridReducedX(const char *outFileName, double *grid, int number, int additionalBinNumber, int step, int rank, int prevRank, int nextRank, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool
                        newFile, double scale = 1.0);
void outputGridReducedY(const char *outFileName, double *grid, int number, int additionalBinNumber, int step, int rank, int prevRank, int nextRank, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool
                        newFile, double scale = 1.0);
void outputGridReducedZ(const char *outFileName, double *grid, int number, int additionalBinNumber, int step, int rank, int prevRank, int nextRank, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool
                        newFile, double scale = 1.0);
void outputFields(const char *outEfileName, const char *outBfileName, Vector3d ***Efield, Vector3d ***Bfield, int xnumberAdded,
                  int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_priod, double gyroradius, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool
                  newFile);
void outputFieldsCrossectionXY(const char *outEfileName, const char *outBfileName, Vector3d ***Efield, Vector3d ***Bfield, int xnumberAdded,
                               int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_priod, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int zindex, bool
                               newFile);
void outputFieldsCrossectionXZ(const char *outEfileName, const char *outBfileName, Vector3d ***Efield, Vector3d ***Bfield, int xnumberAdded,
                               int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_priod, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int yindex, bool
                               newFile);
void outputFieldsCrossectionYZ(const char *outEfileName, const char *outBfileName, Vector3d ***Efield, Vector3d ***Bfield, int xnumberAdded,
                               int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_priod, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int xindex, bool
                               newFile);

void outputFieldsLineX(const char* outEfileName, const char* outBfileName, Vector3d*** Efield, Vector3d*** Bfield, int xnumberAdded,
                       int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period, double gyroradius, MPI_Comm& subCommX, int* cartCoord, int* cartDim, int yindex, int zindex, bool
                       newFile);
void outputFieldsLineY(const char* outEfileName, const char* outBfileName, Vector3d*** Efield, Vector3d*** Bfield, int xnumberAdded,
                       int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period, double gyroradius, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int xindex, int zindex, bool
                       newFile);
void outputFieldsLineZ(const char* outEfileName, const char* outBfileName, Vector3d*** Efield, Vector3d*** Bfield, int xnumberAdded,
                       int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period, double gyroradius, MPI_Comm& subCommZ, int* cartCoord, int* cartDim, int xindex, int yindex, bool
                       newFile);
/*void outputFieldsReduced(const char *outEfileName, const char *outBfileName, Vector3d ***Efield, Vector3d ***Bfield, int xnumberAdded,
                         int ynumberAdded, int znumberAdded, int additionalBinNumber, int stepX, int stepY, int stepZ, double plasma_priod, double gyroradius);*/
void outputConcentrations(const char *outFileName, double ****particleConcentrations, double ***chargeDensity,
                          double ***shiftChargeDensity, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                          int typesNumber, double plasma_period, double gyroradius, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool
                          newFile);
void outputConcentrationsCrossectionXY(const char *outFileName, double ****particleConcentrations, double ***chargeDensity,
                                       double ***shiftChargeDensity, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                                       int typesNumber, double plasma_period, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int zindex, bool
                                       newFile);
void outputConcentrationsCrossectionXZ(const char *outFileName, double ****particleConcentrations, double ***chargeDensity,
                                       double ***shiftChargeDensity, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                                       int typesNumber, double plasma_period, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int yindex, bool
                                       newFile);
void outputConcentrationsCrossectionYZ(const char *outFileName, double ****particleConcentrations, double ***chargeDensity,
                                       double ***shiftChargeDensity, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                                       int typesNumber, double plasma_period, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int xindex, bool
                                       newFile);
void outputConcentrationsLineX(const char *outFileName, double ****particleConcentrations, double ***chargeDensity,
                               double ***shiftChargeDensity, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                               int typesNumber, double plasma_period, double gyroradius, MPI_Comm& subCommX, int* cartCoord, int* cartDim, int yindex, int zindex, bool
                               newFile);
void outputConcentrationsLineY(const char *outFileName, double ****particleConcentrations, double ***chargeDensity,
                               double ***shiftChargeDensity, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                               int typesNumber, double plasma_period, double gyroradius, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int xindex, int zindex, bool
                               newFile);
void outputConcentrationsLineZ(const char *outFileName, double ****particleConcentrations, double ***chargeDensity,
                               double ***shiftChargeDensity, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                               int typesNumber, double plasma_period, double gyroradius, MPI_Comm& subCommZ, int* cartCoord, int* cartDim, int xindex, int yindex, bool
                               newFile);

void outputVelocity(const char *outFileName, Vector3d ****velocity, ParticleTypeContainer* types,
                    int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int typesNumber,
                    double plasma_period, double gyroradius, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool newFile);
void outputVelocityCrossectionXY(const char *outFileName, Vector3d ****velocity, ParticleTypeContainer* types,
                                 int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int typesNumber,
                                 double plasma_period, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int zindex, bool
                                 newFile);
void outputVelocityCrossectionXZ(const char *outFileName, Vector3d ****velocity, ParticleTypeContainer* types,
                                 int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int typesNumber,
                                 double plasma_period, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int yindex, bool
                                 newFile);
void outputVelocityCrossectionYZ(const char *outFileName, Vector3d ****velocity, ParticleTypeContainer* types,
                                 int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int typesNumber,
                                 double plasma_period, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int xindex, bool
                                 newFile);
void outputVelocityLineX(const char *outFileName, Vector3d ****velocity, ParticleTypeContainer* types,
                         int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int typesNumber,
                         double plasma_period, double gyroradius, MPI_Comm& subCommX, int* cartCoord, int* cartDim, int yindex, int zindex, bool
                         newFile);
void outputVelocityLineY(const char *outFileName, Vector3d ****velocity, ParticleTypeContainer* types,
                         int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int typesNumber,
                         double plasma_period, double gyroradius, MPI_Comm& subCommX, int* cartCoord, int* cartDim, int xindex, int zindex, bool
                         newFile);
void outputVelocityLineZ(const char *outFileName, Vector3d ****velocity, ParticleTypeContainer* types,
                         int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int typesNumber,
                         double plasma_period, double gyroradius, MPI_Comm& subCommX, int* cartCoord, int* cartDim, int xindex, int yindex, bool
                         newFile);
void outputGeneral(const char *outFileName, Simulation *simulation);
void outputGeneralAnisotropy(const char *outFileName, Simulation *simulation);
void outputFlux(const char *outFileName, Vector3d ***electricFlux, Vector3d ***externalElectricFlux, int xnumberAdded,
                int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period, double gyroradius, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool
                newFile);
void outputDivergenceError(const char *outFileName, Simulation *simulation, double plasma_period, double gyroradius, bool newFile);
void outputVectorNodeArraySimple(const char *outFileName, Vector3d ***vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                                 int additionalBinNumber, bool newFile, double scale = 1.0);
void outputVectorNodeArray(const char *outFileName, Vector3d ***vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool newFile, double scale = 1.0);
/*void outputVectorNodeArrayReduced(const char *outFileName, Vector3d ***vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int additionalBinNumber, int stepX, int stepY, int stepZ, double scale = 1.0);*/
void outputVectorCellArray(const char *outFileName, Vector3d ***vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool newFile, double scale = 1.0);
void outputVectorCellArrayCrossectionXY(const char *outFileName, Vector3d ***vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                                        int additionalBinNumber, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int zindex, bool
                                        newFile, double scale = 1.0);
void outputVectorCellArrayCrossectionXZ(const char *outFileName, Vector3d ***vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                                        int additionalBinNumber, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int yindex, bool newFile, double scale = 1.0);
void outputVectorCellArrayCrossectionYZ(const char *outFileName, Vector3d ***vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                                        int additionalBinNumber, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int xindex, bool newFile, double scale = 1.0);
void outputVectorCellArrayLineX(const char *outFileName, Vector3d ***vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, int yindex, int zindex, bool newFile, double scale = 1.0);
void outputVectorCellArrayLineY(const char *outFileName, Vector3d ***vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, int xindex, int zindex, bool newFile, double scale = 1.0);
void outputVectorCellArrayLineZ(const char *outFileName, Vector3d ***vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, int xindex, int yindex, bool newFile, double scale = 1.0);
void outputVectorCellArray(const char *outFileName, double ****vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool newFile, double scale = 1.0);
void outputMatrixArray(const char *outFileName, Matrix3d ***matrix3d, int xnumberAdded, int ynumberAdded, int znumberAdded,
                       int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, bool newFile, double scale = 1.0);

void outputSimulationBackup(const char *generalFileName, const char *EfileName, const char *BfileName,
                            const char *particlesFileName, Simulation *simulation);

void outputParticles(const char* outFileName, Simulation* simulation, ParticleTypes type);
void outputAcceleratedParticlesNumbers(const char *outFileName, Simulation* simulation);


void outputBackupParticles(const char *outFileName, Simulation *simulation);
void outputBackupParticle(FILE* outFile, Particle* particle);

void outputGeneralInitialParameters(const char* outFileName, const char* outFileNameWithText, Simulation* simulation);
void outputMemory(const char* outFileName, MPI_Comm& cartComm, int* cartCoord, int* cartDim);
#endif