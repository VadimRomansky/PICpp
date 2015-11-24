#include "stdio.h"
#include "math.h"
#include "vector"

#include "output.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "util.h"

void outputDistributionUpstream(FILE* outFile, std::vector<Particle*> particles, int particleType, double shockWavePoint, double plasma_period, double gyroradius){
	double minMomentum = 0; //todo something else
	double maxMomentum = 0;
	for(int i = 0; i < particles.size(); ++i){
		if(particles[i]->x > shockWavePoint){
			if(particles[i]->type == particleType){
				if(minMomentum <= 0){
					minMomentum = particles[i]->momentum.norm()*gyroradius/plasma_period;
				}
				if(particles[i]->momentum.norm()*gyroradius/plasma_period < minMomentum){
					minMomentum = particles[i]->momentum.norm()*gyroradius/plasma_period;
				} else {
					if(particles[i]->momentum.norm()*gyroradius/plasma_period > maxMomentum){
						maxMomentum = particles[i]->momentum.norm()*gyroradius/plasma_period;
					}
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
		if(particles[i]->x > shockWavePoint){
			if(particles[i]->type == particleType){
				int j = (log(particles[i]->momentum.norm()*gyroradius/plasma_period) - logMinMomentum)/deltaLogP;
				if( j >= 0 && j < pnumber){
					distribution[j] += particles[i]->weight;
					weight += particles[i]->weight;
				}
			}
		}
	}

	for(int i = 0; i < pnumber; ++i){
		distribution[i] /= (weight*(pgrid[i+1] - pgrid[i]));
	}

	for(int i = 0; i < pnumber; ++i){
		fprintf(outFile, "%20.15g %20.15g\n", (pgrid[i] + pgrid[i+1])/2, distribution[i]);
	}
}

void outputDistribution(FILE* outFile, std::vector<Particle*> particles, int particleType, double plasma_period, double gyroradius){
	double minMomentum = 0; //todo something else
	double maxMomentum = 0;
	for(int i = 0; i < particles.size(); ++i){
		if(particles[i]->type == particleType){
			if(minMomentum <= 0){
				minMomentum = particles[i]->momentum.norm()*gyroradius/plasma_period;
			}
			if(particles[i]->momentum.norm()*gyroradius/plasma_period < minMomentum){
				minMomentum = particles[i]->momentum.norm()*gyroradius/plasma_period;
			} else {
				if(particles[i]->momentum.norm()*gyroradius/plasma_period > maxMomentum){
					maxMomentum = particles[i]->momentum.norm()*gyroradius/plasma_period;
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
			int j = (log(particles[i]->momentum.norm()*gyroradius/plasma_period) - logMinMomentum)/deltaLogP;
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
		fprintf(outFile, "%20.15g %20.15g\n", (pgrid[i] + pgrid[i+1])/2, distribution[i]);
	}
}

void outputTraectory(FILE* outFile, Particle* particle, double time, double plasma_period, double gyroradius){
	fprintf(outFile, "%g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g\n", time*plasma_period, particle->x*gyroradius, particle->y*gyroradius, particle->z*gyroradius, particle->momentum.x*gyroradius/plasma_period, particle->momentum.y*gyroradius/plasma_period, particle->momentum.z*gyroradius/plasma_period, particle->momentum.norm()*gyroradius/plasma_period);
}

void outputGrid(FILE* outFile, double* grid, int number, double scale) {
	for(int i = 0; i < number+1; ++i) {
		fprintf(outFile, "%15.10g\n", grid[i]*scale);
	}
}

void outputFields(FILE* outEfile, FILE* outBfile, Vector3d* Efield, Vector3d* Bfield, int xnumber, double plasma_period, double gyroradius, double fieldScale) {
	double scale = fieldScale/(plasma_period*sqrt(gyroradius));
	for(int i = 0; i < xnumber; ++i) {
				fprintf(outBfile, "%15.10g %15.10g %15.10g\n", scale*Bfield[i].x, scale*Bfield[i].y, scale*Bfield[i].z);
	}

	for(int i = 0; i < xnumber + 1; ++i) {
				fprintf(outEfile, "%15.10g %15.10g %15.10g\n", scale*Efield[i].x, scale*Efield[i].y, scale*Efield[i].z);
	}
}

void outputConcentrations(FILE* outFile, double* electronConcentration, double* protonConcentration, double* chargeDensity, double* shiftChargeDensity, int xnumber, double plasma_period, double gyroradius, double fieldScale) {
	for(int i = 0; i < xnumber; ++i) {
		fprintf(outFile, "%15.10g %15.10g %15.10g %15.10g\n", electronConcentration[i]/cube(gyroradius), protonConcentration[i]/cube(gyroradius), chargeDensity[i]/(sqrt(cube(gyroradius))*plasma_period), shiftChargeDensity[i]/(sqrt(cube(gyroradius))*plasma_period));
	}
}

void outputVelocity(FILE* outFile, FILE* outElectronFile, Vector3d* velocity, Vector3d* electronVelocity, int xnumber, double plasma_period, double gyroradius) {
	for(int i = 0; i < xnumber; ++i) {
		fprintf(outFile, "%15.10g %15.10g %15.10g\n", velocity[i].x*gyroradius/plasma_period, velocity[i].y*gyroradius/plasma_period, velocity[i].z*gyroradius/plasma_period);
		fprintf(outElectronFile, "%15.10g %15.10g %15.10g\n", electronVelocity[i].x*gyroradius/plasma_period, electronVelocity[i].y*gyroradius/plasma_period, electronVelocity[i].z*gyroradius/plasma_period);
	}
}

void outputFlux(FILE* outFile, Vector3d* electricFlux, Vector3d* externalElectricFlux, int xnumber, double plasma_period, double gyroradius, double fieldScale) {
	double scale = plasma_period*plasma_period*sqrt(gyroradius);
	for(int i = 0; i < xnumber + 1; ++i) {
				fprintf(outFile, "%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n", electricFlux[i].x/scale, electricFlux[i].y/scale, electricFlux[i].z/scale, externalElectricFlux[i].x/scale, externalElectricFlux[i].y/scale, externalElectricFlux[i].z/scale);
	}
}

void outputArrayParameter(FILE* outFile, double* arrayParameter, int xnumber, double scale = 1.0) {
	for(int i = 0; i < xnumber; ++i){
				fprintf(outFile, "%g\n", arrayParameter[i]*scale);
	}
}

void outputGeneral(FILE* outFile, Simulation* simulation) {
	int particlesCount = simulation->particles.size();
	fprintf(outFile, "%d %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %d\n", simulation->currentIteration, simulation->time, simulation->time*simulation->plasma_period, simulation->particleEnergy,
		simulation->electricFieldEnergy, simulation->magneticFieldEnergy, simulation->energy, simulation->momentum.x, simulation->momentum.y, simulation->momentum.z, simulation->maxEfield.norm(), simulation->maxBfield.norm(), simulation->deltaT, particlesCount);
}

void outputDivergenceError(FILE* outFile, Simulation* simulation, double plasma_period, double gyroradius, double fieldScale) {
	for(int i = 0; i < simulation->xnumber; ++i) {
				double div = simulation->evaluateDivE(i);
				double div2 = simulation->evaluateDivTempE(i);
				fprintf(outFile, "%g %g %g\n", (4*pi*simulation->chargeDensity[i] - div)/(sqrt(cube(gyroradius))*plasma_period), div/(sqrt(cube(gyroradius))*plasma_period), 4*pi*simulation->chargeDensity[i]/(sqrt(cube(gyroradius))*plasma_period));
	}
}

void outputVectorArray(FILE* outFile, Vector3d* vector3d, int number, double scale) {
	for(int i = 0; i < number; ++i) {
		fprintf(outFile, "%g %g %g\n", vector3d[i].x*scale, vector3d[i].y*scale, vector3d[i].z*scale);
	}
}

void outputMatrixArray(FILE* outFile, Matrix3d* matrix3d, int number, double scale){
	for(int i = 0; i < number; ++i) {
		fprintf(outFile, "%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n", 
			matrix3d[i].matrix[0][0]*scale, matrix3d[i].matrix[0][1]*scale, matrix3d[i].matrix[0][2]*scale,
			matrix3d[i].matrix[1][0]*scale, matrix3d[i].matrix[1][1]*scale, matrix3d[i].matrix[1][2]*scale,
			matrix3d[i].matrix[2][0]*scale, matrix3d[i].matrix[2][1]*scale, matrix3d[i].matrix[2][2]*scale);
	}
}

void outputSimulationBackup(FILE* generalFile, FILE* Efile, FILE* Bfile, FILE* particlesFile, Simulation* simulation){

	fprintf(generalFile, "%d\n", simulation->xnumber);
	fprintf(generalFile, "%d\n", simulation->particlesNumber);
	fprintf(generalFile, "%d\n", simulation->particlesPerBin);

	fprintf(generalFile, "%15.10g\n", simulation->density);
	fprintf(generalFile, "%15.10g\n", simulation->temperature);
	fprintf(generalFile, "%15.10g\n", simulation->plasma_period);
	fprintf(generalFile, "%15.10g\n", simulation->gyroradius);
	fprintf(generalFile, "%15.10g\n", simulation->fieldScale);

	fprintf(generalFile, "%15.10g\n", simulation->time);
	fprintf(generalFile, "%15.10g\n", simulation->maxTime);

	fprintf(generalFile, "%d\n", simulation->currentIteration);
	fprintf(generalFile, "%d\n", simulation->maxIteration);

	fprintf(generalFile, "%15.10g\n", simulation->xsize);
	fprintf(generalFile, "%15.10g\n", simulation->theta);
	fprintf(generalFile, "%15.10g\n", simulation->eta);

	int debugMode = 0;
	if(simulation->debugMode){
		debugMode = 1;
	} else {
		debugMode = 0;
	}

	fprintf(generalFile, "%d\n", debugMode);

	int solverType = 0;
	if(simulation->solverType == IMPLICIT){
		solverType = 1;
	} else {
		solverType = 0;
	}

	fprintf(generalFile, "%d\n", solverType);

	int boundaryConditionType = 0;
	if(simulation->boundaryConditionType == PERIODIC){
		boundaryConditionType = 1;
	} else {
		boundaryConditionType = 0;
	}

	fprintf(generalFile, "%d\n", boundaryConditionType);

	fprintf(generalFile, "%d\n", simulation->maxwellEquationMatrixSize);

	fprintf(generalFile, "%15.10g\n", simulation->extJ);

	fprintf(generalFile, "%15.10g\n", simulation->V0.x);
	fprintf(generalFile, "%15.10g\n", simulation->V0.y);
	fprintf(generalFile, "%15.10g\n", simulation->V0.z);
	fprintf(generalFile, "%15.10g\n", simulation->E0.x);
	fprintf(generalFile, "%15.10g\n", simulation->E0.y);
	fprintf(generalFile, "%15.10g\n", simulation->E0.z);
	fprintf(generalFile, "%15.10g\n", simulation->B0.x);
	fprintf(generalFile, "%15.10g\n", simulation->B0.y);
	fprintf(generalFile, "%15.10g\n", simulation->B0.z);

	outputFields(Efile, Bfile, simulation->Efield, simulation->Bfield, simulation->xnumber, 1.0, 1.0, 1.0);
	outputBackupParticles(particlesFile, simulation);
}

void outputParticles(FILE* outProtonsFile, FILE* outElectronsFile, Simulation* simulation){
	for(int i = 0; i < simulation->particles.size(); ++i){
		Particle* particle = simulation->particles[i];
		double p = particle->momentum.norm()*simulation->gyroradius/simulation->plasma_period;
		if(particle->type == PROTON){
			fprintf(outProtonsFile, "%15.10g %15.10g %15.10g\n", particle->x*simulation->gyroradius, p, particle->momentum.x*simulation->gyroradius/simulation->plasma_period);
		} else if(particle->type == ELECTRON){
			fprintf(outElectronsFile, "%15.10g %15.10g %15.10g\n", particle->x*simulation->gyroradius, p,  particle->momentum.x*simulation->gyroradius/simulation->plasma_period);
		}
	}
}

void outputBackupParticles(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->particles.size(); ++i){
		Particle* particle = simulation->particles[i];
		outputBackupParticle(outFile, particle);
	}
}

void outputBackupParticle(FILE* outFile, Particle* particle){
	fprintf(outFile, "%d ", particle->number);
	fprintf(outFile, "%15.10g ", particle->mass);
	fprintf(outFile, "%15.10g ", particle->charge);
	fprintf(outFile, "%15.10g ", particle->weight);
	int type = 0;
	if(particle->type == PROTON){
		type = 1;
	} else if (particle->type == ELECTRON){
		type = 2;
	}
	fprintf(outFile, "%d ", type);
	fprintf(outFile, "%15.10g ", particle->x);
	fprintf(outFile, "%15.10g ", particle->y);
	fprintf(outFile, "%15.10g ", particle->z);

	fprintf(outFile, "%15.10g ", particle->momentum.x);
	fprintf(outFile, "%15.10g ", particle->momentum.y);
	fprintf(outFile, "%15.10g ", particle->momentum.z);

	fprintf(outFile, "%15.10g\n", particle->dx);
}