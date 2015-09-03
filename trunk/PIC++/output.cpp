#include "stdio.h"
#include "math.h"
#include "vector"

#include "output.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "util.h"

void outputDistribution(FILE* outFile, std::vector<Particle*> particles, int particleType){
	double minMomentum = -1; //todo something else
	double maxMomentum = 0;
	for(int i = 0; i < particles.size(); ++i){
		if(particles[i]->type == particleType){
			if(minMomentum < 0){
				minMomentum = particles[i]->momentum.norm();
			}
			if(particles[i]->momentum.norm() < minMomentum){
				minMomentum = particles[i]->momentum.norm();
			} else {
				if(particles[i]->momentum.norm() > maxMomentum){
					maxMomentum = particles[i]->momentum.norm();
				}
			}
		}
	}

	double pgrid[pnumber+1];
	double distribution[pnumber];
	double logMinMomentum = log(minMomentum);
	pgrid[0] = minMomentum;
	distribution[0] = 0;
	double deltaLogP = (log(maxMomentum) - log(minMomentum))/(pnumber);
	for(int i = 1; i < pnumber; ++i){
		distribution[i] = 0;
		pgrid[i] = exp(logMinMomentum + i*deltaLogP);
	}
	pgrid[pnumber] = maxMomentum;

	double weight = 0;

	for(int i = 0; i < particles.size(); ++i){
		if(particles[i]->type == particleType){
			int j = (log(particles[i]->momentum.norm()) - logMinMomentum)/deltaLogP;
			if( j >= 0 && j < pnumber){
				distribution[j] += particles[i]->weight;
				weight += particles[i]->weight;
			}
		}
	}

	for(int i = 0; i < pnumber; ++i){
		distribution[i] /= (weight*(pgrid[i+1] - pgrid[i]));
	}

	for(int i = 0; i < pnumber; ++i){
		fprintf(outFile, "%g %g\n", (pgrid[i] + pgrid[i+1])/2, distribution[i]);
	}
}

void outputTraectory(FILE* outFile, Particle* particle, double time){
	fprintf(outFile, "%g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g\n", time, particle->coordinates.x, particle->coordinates.y, particle->coordinates.z, particle->momentum.x, particle->momentum.y, particle->momentum.z, particle->momentum.norm());
}

void outputGrid(FILE* outFile, double* grid, int number) {
	for(int i = 0; i < number; ++i) {
		fprintf(outFile, "%15.10g %15.10g\n", grid[i], (grid[i+1] + grid[i])/2);
	}
}

void outputFields(FILE* outEfile, FILE* outBfile, Vector3d*** Efield, Vector3d*** Bfield, int xnumber, int ynumber, int znumber, double plasma_preiod, double gyroradius) {
	double scale = 1.0/(plasma_preiod*gyroradius);
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outBfile, "%15.10g %15.10g %15.10g\n", scale*Bfield[i][j][k].x, scale*Bfield[i][j][k].y, scale*Bfield[i][j][k].z);
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outEfile, "%15.10g %15.10g %15.10g\n", scale*Efield[i][j][k].x, scale*Efield[i][j][k].y, scale*Efield[i][j][k].z);
			}
		}
	}
}

void outputConcentrations(FILE* outFile, double*** electronConcentration, double*** protonConcentration, double*** chargeDensity, double*** shiftChargeDensity, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%15.10g %15.10g %15.10g %15.10g\n", electronConcentration[i][j][k]/cube(gyroradius), protonConcentration[i][j][k]/cube(gyroradius), chargeDensity[i][j][k]/(sqrt(cube(gyroradius))*plasma_period), shiftChargeDensity[i][j][k]/(sqrt(cube(gyroradius))*plasma_period));
			}
		}
	}
}

void outputVelocity(FILE* outFile, FILE* outElectronFile, Vector3d*** velocity, Vector3d*** electronVelocity, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%15.10g %15.10g %15.10g\n", velocity[i][j][k].x*gyroradius/plasma_period, velocity[i][j][k].y*gyroradius/plasma_period, velocity[i][j][k].z*gyroradius/plasma_period);
				fprintf(outElectronFile, "%15.10g %15.10g %15.10g\n", electronVelocity[i][j][k].x*gyroradius/plasma_period, electronVelocity[i][j][k].y*gyroradius/plasma_period, electronVelocity[i][j][k].z*gyroradius/plasma_period);
			}
		}
	}
}

void outputFlux(FILE* outFile, Vector3d*** electricFlux, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%15.10g %15.10g %15.10g\n", electricFlux[i][j][k].x, electricFlux[i][j][k].y, electricFlux[i][j][k].z);
			}
		}
	}
}

void outputArrayParameter(FILE* outFile, double*** arrayParameter, int xnumber, int ynumber, int znumber) {
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%g\n", arrayParameter[i][j][k]);
			}
		}
	}
}

void outputGeneral(FILE* outFile, Simulation* simulation) {
	fprintf(outFile, "%d %g %g %g %g %g %g %g %g %g\n", simulation->currentIteration, simulation->time, simulation->time*simulation->plasma_period, simulation->particleEnergy,
		simulation->electricFieldEnergy, simulation->magneticFieldEnergy, simulation->energy, simulation->momentum.x, simulation->momentum.y, simulation->momentum.z);
}

void outputDivergenceError(FILE* outFile, Simulation* simulation) {
	for(int i = 0; i < simulation->xnumber; ++i) {
		for(int j = 0; j < simulation->ynumber; ++j) {
			for(int k = 0; k < simulation->znumber; ++k) {
				double div = simulation->evaluateDivE(i, j, k);
				double div2 = simulation->evaluateDivTempE(i, j, k);
				fprintf(outFile, "%g %g %g\n", div, div - 4*pi*simulation->chargeDensity[i][j][k], div2 - 4*pi*simulation->electricDensity[i][j][k]);
			}
		}
	}
}