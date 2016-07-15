#include "stdio.h"
#include "math.h"
#include "vector"

#include "output.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "util.h"

void outputDistributionUpstream(FILE* outFile, std::vector<Particle*> particles, int particleType, double shockWavePoint, double plasma_period, double gyroradius) {
	double minMomentum = 0; //todo something else
	double maxMomentum = 0;
	for (int i = 0; i < particles.size(); ++i) {
		if (particles[i]->coordinates.x > shockWavePoint) {
			if (particles[i]->type == particleType) {
				if (minMomentum <= 0) {
					minMomentum = particles[i]->momentum.norm() * gyroradius / plasma_period;
				}
				if (particles[i]->momentum.norm() * gyroradius / plasma_period < minMomentum) {
					minMomentum = particles[i]->momentum.norm() * gyroradius / plasma_period;
				} else {
					if (particles[i]->momentum.norm() * gyroradius / plasma_period > maxMomentum) {
						maxMomentum = particles[i]->momentum.norm() * gyroradius / plasma_period;
					}
				}
			}
		}
	}

	double pgrid[pnumber + 1];
	double distribution[pnumber];
	double logMinMomentum = log(minMomentum);
	pgrid[0] = minMomentum;
	distribution[0] = 0;
	double deltaLogP = (log(maxMomentum) - log(minMomentum)) / (pnumber);
	for (int i = 1; i < pnumber; ++i) {
		distribution[i] = 0;
		pgrid[i] = exp(logMinMomentum + i * deltaLogP);
	}
	pgrid[pnumber] = maxMomentum;

	double weight = 0;

	for (int i = 0; i < particles.size(); ++i) {
		if (particles[i]->coordinates.x > shockWavePoint) {
			if (particles[i]->type == particleType) {
				int j = (log(particles[i]->momentum.norm() * gyroradius / plasma_period) - logMinMomentum) / deltaLogP;
				if (j >= 0 && j < pnumber) {
					distribution[j] += particles[i]->weight;
					weight += particles[i]->weight;
				}
			}
		}
	}

	for (int i = 0; i < pnumber; ++i) {
		distribution[i] /= (weight * (pgrid[i + 1] - pgrid[i]));
	}

	for (int i = 0; i < pnumber; ++i) {
		fprintf(outFile, "%20.15g %20.15g\n", (pgrid[i] + pgrid[i + 1]) / 2, distribution[i]);
	}
}

void outputDistribution(FILE* outFile, std::vector<Particle*> particles, int particleType, double gyroradius, double plasma_period) {
	double minMomentum = 0; //todo something else
	double maxMomentum = 0;
	for (int i = 0; i < particles.size(); ++i) {
		if (particles[i]->type == particleType) {
			if (minMomentum <= 0) {
				minMomentum = particles[i]->momentum.norm() * gyroradius / plasma_period;
			}
			if (particles[i]->momentum.norm() * gyroradius / plasma_period < minMomentum) {
				minMomentum = particles[i]->momentum.norm() * gyroradius / plasma_period;
			} else {
				if (particles[i]->momentum.norm() * gyroradius / plasma_period > maxMomentum) {
					maxMomentum = particles[i]->momentum.norm() * gyroradius / plasma_period;
				}
			}
		}
	}

	double pgrid[pnumber + 1];
	double distribution[pnumber];
	double logMinMomentum = log(minMomentum);
	pgrid[0] = minMomentum;
	distribution[0] = 0;
	double deltaLogP = (log(maxMomentum) - log(minMomentum)) / (pnumber);
	for (int i = 1; i < pnumber; ++i) {
		distribution[i] = 0;
		pgrid[i] = exp(logMinMomentum + i * deltaLogP);
	}
	pgrid[pnumber] = maxMomentum;

	double weight = 0;

	for (int i = 0; i < particles.size(); ++i) {
		if (particles[i]->type == particleType) {
			int j = (log(particles[i]->momentum.norm() * gyroradius / plasma_period) - logMinMomentum) / deltaLogP;
			if (j >= 0 && j < pnumber) {
				distribution[j] += particles[i]->weight;
				weight += particles[i]->weight;
			}
		}
	}

	for (int i = 0; i < pnumber; ++i) {
		distribution[i] /= (weight * (pgrid[i + 1] - pgrid[i]));
	}

	for (int i = 0; i < pnumber; ++i) {
		fprintf(outFile, "%20.15g %20.15g\n", (pgrid[i] + pgrid[i + 1]) / 2, distribution[i]);
	}
}

void outputAnisotropy(FILE* outFile, Simulation* simulation, int particleType, double gyroradius, double plasma_period) {
	for (int i = 0; i < simulation->xnumber; ++i) {
		for (int j = 0; j < simulation->ynumber; ++j) {
			for (int k = 0; k < simulation->znumber; ++k) {
				Vector3d meanV = Vector3d(0, 0, 0);
				double concentration = 0;
				int particleCount = 0;
				if(simulation->particlesInBbin[i][j][k].size() > 0) {
					for (int pcount = 0; pcount < simulation->particlesInBbin[i][j][k].size(); ++pcount) {
						Particle *particle = simulation->particlesInBbin[i][j][k][pcount];
						double correlation = simulation->correlationWithBbin(*particle, i, j,
																			 k) / simulation->volumeB(
							i, j, k);
						if (particle->type == particleType) {
							meanV = meanV + particle->velocity(
								simulation->speed_of_light_normalized) * particle->weight * correlation;
							concentration += particle->weight * correlation;
							particleCount++;
						}
					}
					if(particleCount > 0) {
						meanV = meanV / concentration;
						//printf("i = %d meanv = %g %g %g\n", i, meanV.x, meanV.y, meanV.z);

						double parallelV2 = 0;
						double normalV2 = 0;

						for (int pcount = 0; pcount < simulation->particlesInBbin[i][j][k].size(); ++pcount) {
							Particle *particle = simulation->particlesInBbin[i][j][k][pcount];
							double correlation = simulation->correlationWithBbin(*particle, i, j,
																				 k) / simulation->volumeB(
								i, j, k);
							if (particle->type == particleType) {
								parallelV2 = parallelV2 + sqr(particle->velocityX(
									simulation->speed_of_light_normalized) - meanV.x) * particle->weight * correlation;
								normalV2 = normalV2 + (sqr(
									particle->velocityY(simulation->speed_of_light_normalized) - meanV.y) + sqr(
									particle->velocityZ(
										simulation->speed_of_light_normalized)) - meanV.z) * particle->weight * correlation;
							}
						}

						fprintf(outFile, "%g\n", (0.5 * normalV2 / parallelV2) - (2 * parallelV2 / normalV2));
					} else {
						fprintf(outFile, "%g\n", 0.0);
					}
				} else {
					fprintf(outFile, "%g\n", 0.0);
				}
			}
		}
	}
}

void outputTrajectory(FILE* outFile, Particle* particle, double time, double plasma_period, double gyroradius) {
	fprintf(outFile, "%g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g\n", time, particle->coordinates.x, particle->coordinates.y, particle->coordinates.z, particle->momentum.x, particle->momentum.y, particle->momentum.z, particle->momentum.norm());
}

void outputGrid(FILE* outFile, double* grid, int number, double scale) {
	for (int i = 0; i <= number; ++i) {
		fprintf(outFile, "%15.10g\n", grid[i] * scale);
	}
}

void outputFields(FILE* outEfile, FILE* outBfile, Vector3d*** Efield, Vector3d*** Bfield, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius, double fieldScale) {
	double scale = 1.0 / (plasma_period * gyroradius);
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				fprintf(outBfile, "%15.10g %15.10g %15.10g\n", scale * Bfield[i][j][k].x, scale * Bfield[i][j][k].y, scale * Bfield[i][j][k].z);
			}
		}
	}

	for (int i = 0; i <= xnumber; ++i) {
		for (int j = 0; j <= ynumber; ++j) {
			for (int k = 0; k <= znumber; ++k) {
				fprintf(outEfile, "%15.10g %15.10g %15.10g\n", scale * Efield[i][j][k].x, scale * Efield[i][j][k].y, scale * Efield[i][j][k].z);
			}
		}
	}
}

void outputConcentrations(FILE* outFile, double*** electronConcentration, double*** protonConcentration, double*** chargeDensity, double*** shiftChargeDensity, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius, double fieldScale) {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%15.10g %15.10g %15.10g %15.10g\n", electronConcentration[i][j][k] / cube(gyroradius), protonConcentration[i][j][k] / cube(gyroradius), chargeDensity[i][j][k] / (sqrt(cube(gyroradius)) * plasma_period), shiftChargeDensity[i][j][k] / (sqrt(cube(gyroradius)) * plasma_period));
			}
		}
	}
}

void outputVelocity(FILE* outFile, FILE* outElectronFile, Vector3d*** velocity, Vector3d*** electronVelocity, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius) {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%15.10g %15.10g %15.10g\n", velocity[i][j][k].x * gyroradius / plasma_period, velocity[i][j][k].y * gyroradius / plasma_period, velocity[i][j][k].z * gyroradius / plasma_period);
				fprintf(outElectronFile, "%15.10g %15.10g %15.10g\n", electronVelocity[i][j][k].x * gyroradius / plasma_period, electronVelocity[i][j][k].y * gyroradius / plasma_period, electronVelocity[i][j][k].z * gyroradius / plasma_period);
			}
		}
	}
}

void outputFlux(FILE* outFile, Vector3d*** electricFlux, Vector3d*** externalElectricFlux, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius, double fieldScale) {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n", electricFlux[i][j][k].x, electricFlux[i][j][k].y, electricFlux[i][j][k].z, externalElectricFlux[i][j][k].x, externalElectricFlux[i][j][k].y, externalElectricFlux[i][j][k].z);
			}
		}
	}
}

void outputVectorArray(FILE* outFile, Vector3d*** vector3d, int xnumber, int ynumber, int znumber, double scale) {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%15.10g %15.10g %15.10g\n", vector3d[i][j][k].x * scale, vector3d[i][j][k].y * scale, vector3d[i][j][k].z * scale);
			}
		}
	}
}

void outputMatrixArray(FILE* outFile, Matrix3d*** matrix3d, int xnumber, int ynumber, int znumber, double scale) {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",
				        matrix3d[i][j][k].matrix[0][0] * scale, matrix3d[i][j][k].matrix[0][1] * scale, matrix3d[i][j][k].matrix[0][2] * scale,
				        matrix3d[i][j][k].matrix[1][0] * scale, matrix3d[i][j][k].matrix[1][1] * scale, matrix3d[i][j][k].matrix[1][2] * scale,
				        matrix3d[i][j][k].matrix[2][0] * scale, matrix3d[i][j][k].matrix[2][1] * scale, matrix3d[i][j][k].matrix[2][2] * scale);
			}
		}
	}
}

void outputGeneral(FILE* outFile, Simulation* simulation) {
	int particlesCount = simulation->particles.size();
	fprintf(outFile, "%d %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %d\n",
			simulation->currentIteration, simulation->time, simulation->time * simulation->plasma_period, simulation->particleEnergy,
	        simulation->electricFieldEnergy, simulation->magneticFieldEnergy, simulation->energy, simulation->momentum.x, simulation->momentum.y, simulation->momentum.z, simulation->theoreticalEnergy, simulation->theoreticalMomentum.x, simulation->theoreticalMomentum.y,
	        simulation->theoreticalMomentum.z, simulation->maxEfield.norm(), simulation->maxBfield.norm(), simulation->deltaT, particlesCount);
}

void outputDivergenceError(FILE* outFile, Simulation* simulation, double plasma_period, double gyroradius, double fieldScale) {
	for (int i = 0; i < simulation->xnumber; ++i) {
		for (int j = 0; j < simulation->ynumber; ++j) {
			for (int k = 0; k < simulation->znumber; ++k) {
				double div = simulation->evaluateDivE(i, j, k);
				double div2 = simulation->evaluateDivTempE(i, j, k);
				fprintf(outFile, "%g %g %g\n", (4 * pi * simulation->chargeDensity[i][j][k] - div) / (sqrt(cube(gyroradius)) * plasma_period), div / (sqrt(cube(gyroradius)) * plasma_period), 4 * pi * simulation->chargeDensity[i][j][k] / (sqrt(cube(gyroradius)) * plasma_period));
			}
		}
	}
}

void outputParticles(FILE* outProtonsFile, FILE* outElectronsFile, FILE* outPositronsFile, FILE* outAlphaFile, FILE* outDeuteriumFile, FILE* outHelium3File, Simulation* simulation) {
	for (int i = 0; i < simulation->particles.size(); ++i) {
		Particle* particle = simulation->particles[i];
		double p = particle->momentum.norm() * simulation->scaleFactor / simulation->plasma_period;
		if (particle->type == PROTON) {
			fprintf(outProtonsFile, "%15.10g %15.10g %15.10g\n", particle->coordinates.x * simulation->scaleFactor, p, particle->momentum.x * simulation->scaleFactor / simulation->plasma_period);
		} else if (particle->type == ELECTRON) {
			fprintf(outElectronsFile, "%15.10g %15.10g %15.10g\n", particle->coordinates.x * simulation->scaleFactor, p, particle->momentum.x * simulation->scaleFactor / simulation->plasma_period);
		} else if (particle->type == POSITRON) {
			fprintf(outPositronsFile, "%15.10g %15.10g %15.10g\n", particle->coordinates.x * simulation->scaleFactor, p, particle->momentum.x * simulation->scaleFactor / simulation->plasma_period);
		} else if (particle->type == ALPHA) {
			fprintf(outAlphaFile, "%15.10g %15.10g %15.10g\n", particle->coordinates.x * simulation->scaleFactor, p, particle->momentum.x * simulation->scaleFactor / simulation->plasma_period);
		} else if (particle->type == DEUTERIUM) {
			fprintf(outDeuteriumFile, "%15.10g %15.10g %15.10g\n", particle->coordinates.x * simulation->scaleFactor, p, particle->momentum.x * simulation->scaleFactor / simulation->plasma_period);
		} else if (particle->type == HELIUM3) {
			fprintf(outHelium3File, "%15.10g %15.10g %15.10g\n", particle->coordinates.x * simulation->scaleFactor, p, particle->momentum.x * simulation->scaleFactor / simulation->plasma_period);
		}
	}
}

void outputMaxwellEquationMatrixSimple(std::vector<MatrixElement>****& maxwellEquationMatrix, int xnumber, int ynumber, int znumber, int lnumber) {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					//printf("%d %d %d %d\n", i,j, k, l);
					for (int m = 0; m < maxwellEquationMatrix[i][j][k][l].size(); ++m) {
						printf("%15.7g", maxwellEquationMatrix[i][j][k][l][m].value);
					}
					printf("\n");
				}
			}
		}
	}
}

void outputMaxwellEquationMatrixFull(FILE* outFile, std::vector<MatrixElement>****& maxwellEquationMatrix, int xnumber, int ynumber, int znumber, int lnumber) {
	printf("outputingMatrix\n");
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					//printf("%d %d %d %d\n", i,j, k, l);
					for(int tempI = 0; tempI < xnumber; ++tempI) {
						for(int tempJ = 0; tempJ < ynumber; ++tempJ) {
							for(int tempK = 0; tempK < znumber; ++tempK) {
								for(int tempL = 0; tempL < lnumber; ++tempL){
									MatrixElement element = MatrixElement(0, tempI, tempJ, tempK, tempL);
									for (int m = 0; m < maxwellEquationMatrix[i][j][k][l].size(); ++m) {
										if(element.equalsIndex(maxwellEquationMatrix[i][j][k][l][m])) {
											element.value = maxwellEquationMatrix[i][j][k][l][m].value;
										}
									}
									fprintf(outFile, "%28.22g ", element.value)
;								}
							}
						}
					}
					fprintf(outFile, "\n");
				}
			}
		}
	}
}

void outputSimulationBackup(FILE* generalFile, FILE* Efile, FILE* Bfile, FILE* particlesFile, Simulation* simulation) {
	int inputType = 2;
	if(simulation->inputType == CGS){
		inputType = 0;
	} else if(simulation->inputType == Theoretical){
		inputType = 1;
	}
	fprintf(generalFile, "%d\n", inputType);
	fprintf(generalFile, "%d\n", simulation->xnumber);
	fprintf(generalFile, "%d\n", simulation->ynumber);
	fprintf(generalFile, "%d\n", simulation->znumber);
	fprintf(generalFile, "%d\n", simulation->particlesNumber);
	fprintf(generalFile, "%d\n", simulation->electronsPerBin);
	fprintf(generalFile, "%d\n", simulation->protonsPerBin);
	fprintf(generalFile, "%d\n", simulation->positronsPerBin);
	fprintf(generalFile, "%d\n", simulation->alphaPerBin);

	fprintf(generalFile, "%15.10g\n", simulation->density);
	fprintf(generalFile, "%15.10g\n", simulation->temperature);
	fprintf(generalFile, "%15.10g\n", simulation->plasma_period);
	fprintf(generalFile, "%15.10g\n", simulation->scaleFactor);
	fprintf(generalFile, "%15.10g\n", simulation->fieldScale);

	fprintf(generalFile, "%15.10g\n", simulation->time);
	fprintf(generalFile, "%15.10g\n", simulation->maxTime);

	fprintf(generalFile, "%d\n", simulation->currentIteration);
	fprintf(generalFile, "%d\n", simulation->maxIteration);

	fprintf(generalFile, "%15.10g\n", simulation->xsize);
	fprintf(generalFile, "%15.10g\n", simulation->ysize);
	fprintf(generalFile, "%15.10g\n", simulation->zsize);
	fprintf(generalFile, "%15.10g\n", simulation->theta);
	fprintf(generalFile, "%15.10g\n", simulation->eta);

	int debugMode = 0;
	if (simulation->debugMode) {
		debugMode = 1;
	} else {
		debugMode = 0;
	}

	fprintf(generalFile, "%d\n", debugMode);

	int preserveChargeLocal = 0;
	if (simulation->preserveChargeLocal) {
		preserveChargeLocal = 1;
	} else {
		preserveChargeLocal = 0;
	}

	fprintf(generalFile, "%d\n", preserveChargeLocal);

	int preserveChargeGlobal = 0;
	if (simulation->preserveChargeGlobal) {
		preserveChargeGlobal = 1;
	} else {
		preserveChargeGlobal = 0;
	}

	fprintf(generalFile, "%d\n", preserveChargeGlobal);

	int solverType = 0;
	if (simulation->solverType == IMPLICIT) {
		solverType = 1;
	} else {
		solverType = 0;
	}

	fprintf(generalFile, "%d\n", solverType);

	int boundaryConditionType = 0;
	if (simulation->boundaryConditionType == PERIODIC) {
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

	outputFields(Efile, Bfile, simulation->Efield, simulation->Bfield, simulation->xnumber, simulation->ynumber, simulation->znumber, 1.0, 1.0, 1.0);
	outputBackupParticles(particlesFile, simulation);
}

void outputBackupParticles(FILE* outFile, Simulation* simulation) {
	for (int i = 0; i < simulation->particles.size(); ++i) {
		Particle* particle = simulation->particles[i];
		outputBackupParticle(outFile, particle);
	}
}

void outputBackupParticle(FILE* outFile, Particle* particle) {
	fprintf(outFile, "%d ", particle->number);
	fprintf(outFile, "%15.10g ", particle->mass);
	fprintf(outFile, "%15.10g ", particle->charge);
	fprintf(outFile, "%15.10g ", particle->weight);
	int type = 0;
	if (particle->type == ELECTRON) {
		type = 1;
	} else if (particle->type == PROTON) {
		type = 2;
	} else if (particle->type == POSITRON) {
		type = 3;
	} else if (particle->type == ALPHA) {
		type = 4;
	}
	fprintf(outFile, "%d ", type);
	fprintf(outFile, "%15.10g ", particle->coordinates.x);
	fprintf(outFile, "%15.10g ", particle->coordinates.y);
	fprintf(outFile, "%15.10g ", particle->coordinates.z);

	fprintf(outFile, "%15.10g ", particle->momentum.x);
	fprintf(outFile, "%15.10g ", particle->momentum.y);
	fprintf(outFile, "%15.10g ", particle->momentum.z);

	fprintf(outFile, "%15.10g\n", particle->dx);
}
