#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "vector"
#include "stdio.h"

#include "vector3d.h"
#include "particle.h"
#include "simulation.h"

void outputDistribution(FILE* outFile, std::vector<Particle*> particles, int particleType);
void outputTraectory(FILE* outFile, Particle* particle, double time);
void outputGrid(FILE* outFile, double* grid, int number);
void outputFields(FILE* outEfile, FILE* outBfile, Vector3d*** Efield, Vector3d*** Bfield, int xnumber, int ynumber, int znumber, double plasma_priod, double gyroradius);
void outputConcentrations(FILE* outFile, double*** electronConcentration, double*** protonConcentration, double*** chargeDensity, double*** shiftChargeDensity, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius);
void outputVelocity(FILE* outFile, FILE* outElectronFile, Vector3d*** velocity, Vector3d*** electronVelocity, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius);
void outputFlux(FILE* outFile, Vector3d*** electricFlux, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius);
void outputArrayParameter(FILE* outFile, double*** arrayParameter, int xnumber, int ynumber, int znumber);
void outputGeneral(FILE* outFile, Simulation* simulatiom);
void outputDivergenceError(FILE* outFile, Simulation* simulation);


#endif