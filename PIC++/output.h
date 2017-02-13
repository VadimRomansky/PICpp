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
                        double gyroradius, int verbosity);
void outputDistributionShiftedSystem(const char* outFileName, std::vector<Particle *>& particles, Vector3d& shiftV, double& speed_of_light_normalized, int particleType, double gyroradius,
                        double plasma_period, int verbosity);
//void outputAnisotropy(const char *outFileName, Simulation *simulation, int particleType, double gyroradius,
//                      double plasma_period);
void outputTrajectoryByNumber(const char* outFileName, int number, const Simulation* simulation);
void outputTrajectory(const char *outFileName, Particle *particle, double time, double plasma_period, double gyroradius);
void outputGridSimple(const char *outFileName, double *grid, int number, double scale = 1.0);

void outputGrid(const char *outFileName, double *grid, int number, double scale = 1.0);
void outputFields(const char *outEfileName, const char *outBfileName, Vector3d ***Efield, Vector3d ***Bfield, int xnumber,
                  int ynumber, int znumber, double plasma_priod, double gyroradius);
void outputConcentrations(const char *outFileName, double ****particleConcentrations, double ***chargeDensity,
                          double ***shiftChargeDensity, int xnumber, int ynumber, int znumber, int typesNumber,
                          double plasma_period, double gyroradius);
void outputVelocity(const char *outFileName, Vector3d ****velocity, ParticleTypeContainer* types,
                    int xnumber, int ynumber, int znumber, int typesNumber, double plasma_period,
                    double gyroradius);
void outputGeneral(const char *outFileName, Simulation *simulation);
void outputGeneralAnisotropy(const char *outFileName, Simulation *simulation);
void outputFlux(const char *outFileName, Vector3d ***electricFlux, Vector3d ***externalElectricFlux, int xnumber,
                int ynumber, int znumber, double plasma_period, double gyroradius);
void outputDivergenceError(const char *outFileName, Simulation *simulation, double plasma_period, double gyroradius);
void outputVectorNodeArraySimple(const char *outFileName, Vector3d ***vector3d, int xnumber, int ynumber, int znumber,
                           double scale = 1.0);
void outputVectorNodeArray(const char *outFileName, Vector3d ***vector3d, int xnumber, int ynumber, int znumber,
                           double scale = 1.0);
void outputVectorCellArray(const char *outFileName, Vector3d ***vector3d, int xnumber, int ynumber, int znumber,
                           double scale = 1.0);
void outputMatrixArray(const char *outFileName, Matrix3d ***matrix3d, int xnumber, int ynumber, int znumber,
                       double scale = 1.0);


void outputSimulationBackup(const char *generalFileName, const char *EfileName, const char *BfileName,
                            const char *particlesFileName, Simulation *simulation);

void outputParticles(const char* outFileName, Simulation* simulation, ParticleTypes type);


void outputBackupParticles(const char *outFileName, Simulation *simulation);
void outputBackupParticle(FILE* outFile, Particle* particle);

void outputGeneralInitialParameters(const char* outFileName, const char* outFileNameWithText, Simulation* simulation);

#endif