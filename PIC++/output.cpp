#include <mpi.h>
#include "stdio.h"
#include "math.h"
#include "vector"
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "output.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "util.h"
#include "simulation.h"
#include "matrixElement.h"
#include "mpi_util.h"
#include "memory_util.h"
#include "paths.h"

void outputDistribution(const char* outFileName, std::vector < Particle * >& particles, int particleType,
                        double gyroradius,
                        double plasma_period, int verbosity) {
	int rank;
	int nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	double minMomentum = 0; //todo something else
	double maxMomentum = 0;
	for (int i = 0; i < particles.size(); ++i) {
		if (particles[i]->type == particleType) {
			Vector3d momentum = particles[i]->getMomentum();
			if (minMomentum <= 0) {
				minMomentum = momentum.norm() * gyroradius / plasma_period;
			}
			if ((momentum.norm() * gyroradius / plasma_period < minMomentum) && (momentum.norm() > 0)) {
				minMomentum = momentum.norm() * gyroradius / plasma_period;
			} else {
				if (momentum.norm() * gyroradius / plasma_period > maxMomentum) {
					maxMomentum = momentum.norm() * gyroradius / plasma_period;
				}
			}
		}
	}

	if ((rank == 0) && (verbosity > 2)) printf("send minmax p\n");
	double minMaxP[2];
	if (rank > 0) {
		minMaxP[0] = minMomentum;
		minMaxP[1] = maxMomentum;
		MPI_Send(minMaxP, 2, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD);
	} else {
		MPI_Status status;
		for (int i = 1; i < nprocs; ++i) {
			MPI_Recv(minMaxP, 2, MPI_DOUBLE, i, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD, &status);
			if (minMaxP[0] < minMomentum && minMaxP[0] > 0) {
				minMomentum = minMaxP[0];
			}
			if (minMaxP[1] > maxMomentum) {
				maxMomentum = minMaxP[1];
			}
		}

		if (maxMomentum - minMomentum < 1E-8 * maxMomentum) {
			maxMomentum = maxMomentum * 2;
			minMomentum = minMomentum / 2;
		}

		minMaxP[0] = minMomentum;
		minMaxP[1] = maxMomentum;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if ((rank == 0) && (verbosity > 2)) printf("bcast minmax p\n");

	MPI_Bcast(minMaxP, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	minMomentum = minMaxP[0];
	maxMomentum = minMaxP[1];

	if ((rank == 0) && (verbosity > 2)) printf("evaluate distribution\n");


	double pgrid[pnumber + 1];
	double distribution[pnumber];
	double tempDistribution[pnumber];
	double logMinMomentum = log(minMomentum);
	pgrid[0] = minMomentum;
	distribution[0] = 0;
	double deltaLogP = (log(maxMomentum) - log(minMomentum)) / (pnumber);
	for (int i = 1; i < pnumber; ++i) {
		distribution[i] = 0;
		pgrid[i] = exp(logMinMomentum + i * deltaLogP);
	}
	pgrid[pnumber] = maxMomentum;

	double weight[1];
	weight[0] = 0;

	for (int i = 0; i < particles.size(); ++i) {
		if (particles[i]->type == particleType) {
			int j = (log(particles[i]->getMomentum().norm() * gyroradius / plasma_period) - logMinMomentum) / deltaLogP;
			if (j >= 0 && j < pnumber) {
				distribution[j] += particles[i]->weight;
				weight[0] += particles[i]->weight;
			}
		}
	}

	if ((rank == 0) && (verbosity > 2)) printf("send weight\n");

	if (rank > 0) {
		MPI_Send(weight, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD);
	} else {
		for (int i = 1; i < nprocs; ++i) {
			MPI_Status status;
			double tempWeight[1];
			MPI_Recv(tempWeight, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD, &status);
			weight[0] += tempWeight[0];
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if ((rank == 0) && (verbosity > 2)) printf("send distribution\n");
	if (rank > 0) {
		MPI_Send(distribution, pnumber, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD);
	} else {
		for (int i = 1; i < nprocs; ++i) {
			MPI_Status status;
			MPI_Recv(tempDistribution, pnumber, MPI_DOUBLE, i, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD, &status);
			for (int p = 0; p < pnumber; ++p) {
				distribution[p] += tempDistribution[p];
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if ((rank == 0) && (verbosity > 2)) printf("write distribution\n");

	if (rank == 0) {
		for (int i = 0; i < pnumber; ++i) {
			distribution[i] /= (weight[0] * (pgrid[i + 1] - pgrid[i]));
			//todo check
			if (maxMomentum <= 0 || minMomentum <= 0) {
				distribution[i] = 0;
				pgrid[i] = 0;
			}
		}

		FILE* outFile = fopen(outFileName, "a");
		for (int i = 0; i < pnumber; ++i) {
			fprintf(outFile, "%22.15g %22.15g\n", (pgrid[i] + pgrid[i + 1]) / 2, distribution[i]);
		}
		fclose(outFile);
	}
}

void outputDistributionShiftedSystem(const char* outFileName, std::vector < Particle * >& particles, Vector3d& shiftV,
                                     double& speed_of_light_normalized, int particleType, double gyroradius,
                                     double plasma_period, int verbosity) {
	int rank;
	int nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	double minMomentum = 0; //todo something else
	double maxMomentum = 0;
	for (int i = 0; i < particles.size(); ++i) {
		if (particles[i]->type == particleType) {
			Particle tempParticle = *particles[i];
			tempParticle.addVelocity(shiftV, speed_of_light_normalized);
			Vector3d momentum = tempParticle.getMomentum();
			if (minMomentum <= 0) {
				minMomentum = momentum.norm() * gyroradius / plasma_period;
			}
			if ((momentum.norm() * gyroradius / plasma_period < minMomentum) && (momentum.norm() > 0)) {
				minMomentum = momentum.norm() * gyroradius / plasma_period;
			} else {
				if (momentum.norm() * gyroradius / plasma_period > maxMomentum) {
					maxMomentum = momentum.norm() * gyroradius / plasma_period;
				}
			}
		}
	}

	if ((rank == 0) && (verbosity > 2)) printf("send minmax p\n");
	double minMaxP[2];
	if (rank > 0) {
		minMaxP[0] = minMomentum;
		minMaxP[1] = maxMomentum;
		MPI_Send(minMaxP, 2, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD);
	} else {
		MPI_Status status;
		for (int i = 1; i < nprocs; ++i) {
			MPI_Recv(minMaxP, 2, MPI_DOUBLE, i, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD, &status);
			if (minMaxP[0] < minMomentum && minMaxP[0] > 0) {
				minMomentum = minMaxP[0];
			}
			if (minMaxP[1] > maxMomentum) {
				maxMomentum = minMaxP[1];
			}
		}

		if (maxMomentum - minMomentum < 1E-8 * maxMomentum) {
			maxMomentum = maxMomentum * 2;
			minMomentum = minMomentum / 2;
		}

		minMaxP[0] = minMomentum;
		minMaxP[1] = maxMomentum;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if ((rank == 0) && (verbosity > 2)) printf("bcast minmax p\n");

	MPI_Bcast(minMaxP, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	minMomentum = minMaxP[0];
	maxMomentum = minMaxP[1];

	if ((rank == 0) && (verbosity > 2)) printf("evaluate distribution\n");


	double pgrid[pnumber + 1];
	double distribution[pnumber];
	double tempDistribution[pnumber];
	double logMinMomentum = log(minMomentum);
	pgrid[0] = minMomentum;
	distribution[0] = 0;
	double deltaLogP = (log(maxMomentum) - log(minMomentum)) / (pnumber);
	for (int i = 1; i < pnumber; ++i) {
		distribution[i] = 0;
		pgrid[i] = exp(logMinMomentum + i * deltaLogP);
	}
	pgrid[pnumber] = maxMomentum;

	double weight[1];
	weight[0] = 0;

	for (int i = 0; i < particles.size(); ++i) {
		if (particles[i]->type == particleType) {
			Particle tempParticle = *particles[i];
			tempParticle.addVelocity(shiftV, speed_of_light_normalized);
			int j = (log(tempParticle.getMomentum().norm() * gyroradius / plasma_period) - logMinMomentum) / deltaLogP;
			if (j >= 0 && j < pnumber) {
				distribution[j] += particles[i]->weight;
				weight[0] += particles[i]->weight;
			}
		}
	}

	if ((rank == 0) && (verbosity > 2)) printf("send weight\n");

	if (rank > 0) {
		MPI_Send(weight, 1, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD);
	} else {
		for (int i = 1; i < nprocs; ++i) {
			MPI_Status status;
			double tempWeight[1];
			MPI_Recv(tempWeight, 1, MPI_DOUBLE, i, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD, &status);
			weight[0] += tempWeight[0];
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if ((rank == 0) && (verbosity > 2)) printf("send distribution\n");
	if (rank > 0) {
		MPI_Send(distribution, pnumber, MPI_DOUBLE, 0, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD);
	} else {
		for (int i = 1; i < nprocs; ++i) {
			MPI_Status status;
			MPI_Recv(tempDistribution, pnumber, MPI_DOUBLE, i, MPI_SEND_DOUBLE_ALL_TO_FIRST, MPI_COMM_WORLD, &status);
			for (int p = 0; p < pnumber; ++p) {
				distribution[p] += tempDistribution[p];
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if ((rank == 0) && (verbosity > 2)) printf("write distribution\n");

	if (rank == 0) {
		for (int i = 0; i < pnumber; ++i) {
			distribution[i] /= (weight[0] * (pgrid[i + 1] - pgrid[i]));
			//todo check
			if (maxMomentum <= 0 || minMomentum <= 0) {
				distribution[i] = 0;
				pgrid[i] = 0;
			}
		}

		FILE* outFile = fopen(outFileName, "a");
		for (int i = 0; i < pnumber; ++i) {
			fprintf(outFile, "%22.15g %22.15g\n", (pgrid[i] + pgrid[i + 1]) / 2, distribution[i]);
		}
		fclose(outFile);
	}
}

/*void outputAnisotropy(const char* outFileName, Simulation* simulation, int particleType, double gyroradius,
                      double plasma_period) {
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int minI = 1;
	if (rank == 0) {
		minI = 0;
	}
    int maxI = simulation->xnumber;
    if(rank == size - 1){
        maxI = simulation->xnumber + 1;
    }

	for (int procCount = 0; procCount < size; ++procCount) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (procCount == rank) {
			FILE* outFile = fopen(outFileName, "a");
			for (int i = minI; i < maxI; ++i) {
				for (int j = 0; j < simulation->ynumber; ++j) {
					for (int k = 0; k < simulation->znumber; ++k) {
						Vector3d meanV = Vector3d(0, 0, 0);
						double concentration = 0;
						int particleCount = 0;
						if (simulation->particlesInBbin[i][j][k].size() > 0) {
							for (int pcount = 0; pcount < simulation->particlesInBbin[i][j][k].size(); ++pcount) {
								Particle* particle = simulation->particlesInBbin[i][j][k][pcount];
								double correlation = simulation->correlationWithBbin(*particle, i, j,
								                                                     k) / simulation->volumeB(
									i, j, k);
								if (particle->type == particleType) {
									//meanV = meanV + particle->velocity(
									//	simulation->speed_of_light_normalized) * particle->weight * correlation;
									concentration += particle->weight * correlation;
									particleCount++;
								}
							}
							if (particleCount > 0) {
								meanV = meanV / concentration;
								//printf("i = %d meanv = %g %g %g\n", i, meanV.x, meanV.y, meanV.z);

								double parallelV2 = 0;
								double normalV2 = 0;

								for (int pcount = 0; pcount < simulation->particlesInBbin[i][j][k].size(); ++pcount) {
									Particle* particle = simulation->particlesInBbin[i][j][k][pcount];
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

								fprintf(outFile, "%g %d\n", (0.5 * normalV2 / parallelV2) - (2 * parallelV2 / normalV2), particleCount);
							} else {
								fprintf(outFile, "%g %d\n", 0.0, 0);
							}
						} else {
							fprintf(outFile, "%g\n", 0.0);
						}
					}
				}
			}
			fclose(outFile);
		}
	}
}*/

void outputTrajectoryByNumber(const char* outFileName, int number, const Simulation* simulation) {
	int rank;
	int nprocs;
	MPI_Comm_rank(simulation->cartComm, &rank);
	MPI_Comm_size(simulation->cartComm, &nprocs);
	const int parametersNumber = 9;
	double parameters[parametersNumber];
	for (int i = 0; i < parametersNumber; ++i) {
		parameters[i] = 0;
	}
	int tempRankWithParticle[1];
	int rankWithParticle = 0;
	tempRankWithParticle[0] = 0;
	parameters[0] = simulation->time;
	parameters[1] = simulation->time * simulation->plasma_period;

	for (int p = 0; p < simulation->particles.size(); ++p) {
		Particle* particle = simulation->particles[p];
		if (particle->number == number) {
			tempRankWithParticle[0] = rank;
			parameters[2] = particle->coordinates.x * simulation->scaleFactor;
			parameters[3] = particle->coordinates.y * simulation->scaleFactor;
			parameters[4] = particle->coordinates.z * simulation->scaleFactor;

			Vector3d momentum = particle->getMomentum();
			parameters[5] = momentum.x * simulation->scaleFactor / simulation->plasma_period;
			parameters[6] = momentum.y * simulation->scaleFactor / simulation->plasma_period;
			parameters[7] = momentum.z * simulation->scaleFactor / simulation->plasma_period;

			parameters[8] = momentum.norm() * simulation->scaleFactor / simulation->plasma_period;

			break;
		}
	}

	MPI_Barrier(simulation->cartComm);


	if (rank > 0) {
		MPI_Send(tempRankWithParticle, 1, MPI_INT, 0, MPI_SEND_INTEGER_ALL_TO_FIRST, simulation->cartComm);
	} else {
		for (int i = 1; i < nprocs; ++i) {
			MPI_Status status;
			MPI_Recv(tempRankWithParticle, 1, MPI_INT, i, MPI_SEND_INTEGER_ALL_TO_FIRST, simulation->cartComm, &status);
			if (tempRankWithParticle[0] != 0) {
				rankWithParticle = tempRankWithParticle[0];
			}
		}
	}

	tempRankWithParticle[0] = rankWithParticle;
	MPI_Bcast(tempRankWithParticle, 1, MPI_INT, 0, simulation->cartComm);
	rankWithParticle = tempRankWithParticle[0];

	MPI_Bcast(parameters, parametersNumber, MPI_DOUBLE, rankWithParticle, simulation->cartComm);

	if (rank == 0) {
		FILE* outFile = fopen(outFileName, "a");
		for (int i = 0; i < parametersNumber; ++i) {
			fprintf(outFile, "%20.15g  ", parameters[i]);
		}
		fprintf(outFile, "%d", rankWithParticle);
		fprintf(outFile, "\n");
		fclose(outFile);
	}
}

void outputTrajectory(const char* outFileName, Particle* particle, double time, double plasma_period,
                      double gyroradius) {
	FILE* outFile = fopen(outFileName, "a");
	Vector3d momentum = particle->getMomentum();
	fprintf(outFile, "%20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %d\n", time,
	        time * plasma_period, particle->coordinates.x * gyroradius, particle->coordinates.y * gyroradius,
	        particle->coordinates.z * gyroradius, momentum.x * gyroradius / plasma_period,
	        momentum.y * gyroradius / plasma_period, momentum.z * gyroradius / plasma_period,
	        momentum.norm() * gyroradius / plasma_period, particle->number);
	fclose(outFile);
}

void outputParticlesTrajectories(const char* outFileName, const char* electronOutFileName,
                                 std::vector < Particle* > particles, int** numbers, int size, double time,
                                 double plasma_period, double scaleFactor, Simulation* simulation) {
	int rank;
	int nprocs;
	MPI_Comm_rank(simulation->cartComm, &rank);
	MPI_Comm_size(simulation->cartComm, &nprocs);
	FILE* outFile = NULL;
	//FILE* electronOutFile;
	if (rank == 0) {
		outFile = fopen(outFileName, "a");
		//electronOutFile = fopen(electronOutFileName, "a");
		fprintf(outFile, "%15.10g %15.10g ", time, time * plasma_period);
		//fprintf(electronOutFile, "%15.10g %15.10g ", time, time * plasma_period);
	}
	Vector3d V = simulation->V0;
	double gamma = 1.0 / sqrt(1 - sqr(V.norm() / simulation->speed_of_light_normalized));
	for (int i = 0; i < size; ++i) {
		const int parametersNumber = 8;
		double parameters[parametersNumber];
		for (int j = 0; j < parametersNumber; ++j) {
			parameters[j] = 0;
		}
		parameters[0] = 2 * simulation->xsizeGeneral * scaleFactor;
		parameters[1] = 1.5 * simulation->ysizeGeneral * scaleFactor;
		parameters[2] = 1.5 * simulation->zsizeGeneral * scaleFactor;
		parameters[parametersNumber - 1] = numbers[i][1];
		Vector3d momentum = Vector3d(0, 0, 0);
		momentum.x = simulation->types[numbers[i][1]].mass * V.x * gamma * scaleFactor / plasma_period;
		momentum.y = simulation->types[numbers[i][1]].mass * V.y * gamma * scaleFactor / plasma_period;
		momentum.z = simulation->types[numbers[i][1]].mass * V.z * gamma * scaleFactor / plasma_period;
		parameters[3] = momentum.x;
		parameters[4] = momentum.y;
		parameters[5] = momentum.z;
		parameters[6] = momentum.norm();
		int tempRankWithParticle[1];
		int rankWithParticle = 0;
		tempRankWithParticle[0] = 0;

		for (int p = 0; p < particles.size(); ++p) {
			Particle* particle = particles[p];
			if (particle->number == numbers[i][0]) {
				tempRankWithParticle[0] = rank;
				parameters[0] = particle->coordinates.x * scaleFactor;
				parameters[1] = particle->coordinates.y * scaleFactor;
				parameters[2] = particle->coordinates.z * scaleFactor;

				momentum = particle->getMomentum();
				parameters[3] = momentum.x * scaleFactor / plasma_period;
				parameters[4] = momentum.y * scaleFactor / plasma_period;
				parameters[5] = momentum.z * scaleFactor / plasma_period;

				parameters[6] = momentum.norm() * scaleFactor / plasma_period;

				switch (particle->type) {
				case ELECTRON:
					parameters[7] = 0;
					break;
				case PROTON:
					parameters[7] = 1;
					break;
				case POSITRON:
					parameters[7] = 2;
					break;
				case ALPHA:
					parameters[7] = 3;
					break;
				case DEUTERIUM:
					parameters[7] = 4;
					break;
				case HELIUM3:
					parameters[7] = 5;
					break;
				case OXYGEN_PLUS3:
					parameters[7] = 6;
					break;
				case SILICON_PLUS1:
					parameters[7] = 7;
					break;
				}

				if (parameters[7] != numbers[i][1]) {
					printf("particle type does not match in read tracked particle\n");
				}

				break;
			}
		}

		MPI_Barrier(simulation->cartComm);


		if (rank > 0) {
			MPI_Send(tempRankWithParticle, 1, MPI_INT, 0, MPI_SEND_INTEGER_ALL_TO_FIRST, simulation->cartComm);
		} else {
			for (int j = 1; j < nprocs; ++j) {
				MPI_Status status;
				MPI_Recv(tempRankWithParticle, 1, MPI_INT, j, MPI_SEND_INTEGER_ALL_TO_FIRST, simulation->cartComm, &status);
				if (tempRankWithParticle[0] != 0) {
					rankWithParticle = tempRankWithParticle[0];
				}
			}
		}

		tempRankWithParticle[0] = rankWithParticle;
		MPI_Bcast(tempRankWithParticle, 1, MPI_INT, 0, simulation->cartComm);
		rankWithParticle = tempRankWithParticle[0];

		MPI_Bcast(parameters, parametersNumber, MPI_DOUBLE, rankWithParticle, simulation->cartComm);

		//if (parameters[7] > 0) {
		if (rank == 0) {
			//if(parameters[7] > 0){
			for (int j = 0; j < parametersNumber; ++j) {
				fprintf(outFile, "%15.10g  ", parameters[j]);
			}
			/*} else {
				for (int j = 0; j < parametersNumber; ++j) {
					fprintf(electronOutFile, "%15.10g  ", parameters[j]);
				}
			}*/
		}
		//}
	}
	if (rank == 0) {
		fprintf(outFile, "\n");
		fclose(outFile);
		//fprintf(electronOutFile, "\n");
		//fclose(electronOutFile);
	}
}

void outputGridSimple(const char* outFileName, double* grid, int number, double scale) {
	FILE* outFile = fopen(outFileName, "w");
	for (int i = 0; i <= number; ++i) {
		fprintf(outFile, "%15.10g\n", grid[i] * scale);
	}
	fclose(outFile);
}

void outputGridSimpleReduced(const char* outFileName, double* grid, int number, int step, double scale) {
	FILE* outFile = fopen(outFileName, "w");
	for (int i = 0; i <= number; ++i) {
		if (i % step == 0) {
			fprintf(outFile, "%15.10g\n", grid[i] * scale);
		}
	}
	fclose(outFile);
}

void outputGridX(const char* outFileName, double* grid, int number, int additionalBinNumber, MPI_Comm& cartComm,
                 int* cartCoord, int* cartDim, double scale) {
	int rank;
	int size;
	MPI_Comm_rank(cartComm, &rank);
	MPI_Comm_size(cartComm, &size);
	int minI = 2 + additionalBinNumber;
	if (cartCoord[0] == 0) {
		minI = 0;
	}
	int maxI = number - 1 - additionalBinNumber;
	if (cartCoord[0] == cartDim[0] - 1) {
		maxI = number;
	}

	FILE* outFile = NULL;

	if (rank == 0) {
		outFile = fopen(outFileName, "w");
		fclose(outFile);
	}

	for (int procCount = 0; procCount < cartDim[0]; ++procCount) {
		MPI_Barrier(cartComm);
		if (procCount == cartCoord[0] && cartCoord[1] == 0 && cartCoord[2] == 0) {
			outFile = fopen(outFileName, "a");
			for (int i = minI; i <= maxI; ++i) {
				fprintf(outFile, "%15.10g\n", grid[i] * scale);
			}
			fclose(outFile);
		}
	}
}

void outputGridY(const char* outFileName, double* grid, int number, int additionalBinNumber, MPI_Comm& cartComm,
                 int* cartCoord, int* cartDim, double scale) {
	int rank;
	int size;
	MPI_Comm_rank(cartComm, &rank);
	MPI_Comm_size(cartComm, &size);
	int minJ = 2 + additionalBinNumber;
	if (cartCoord[1] == 0) {
		minJ = 0;
	}
	int maxJ = number - 1 - additionalBinNumber;
	if (cartCoord[1] == cartDim[1] - 1) {
		maxJ = number;
	}

	FILE* outFile = NULL;

	if (rank == 0) {
		outFile = fopen(outFileName, "w");
		fclose(outFile);
	}

	for (int procCount = 0; procCount < cartDim[1]; ++procCount) {
		MPI_Barrier(cartComm);
		if (procCount == cartCoord[1] && cartCoord[0] == 0 && cartCoord[2] == 0) {
			outFile = fopen(outFileName, "a");
			for (int j = minJ; j <= maxJ; ++j) {
				fprintf(outFile, "%15.10g\n", grid[j] * scale);
			}
			fclose(outFile);
		}
	}
}

void outputGridZ(const char* outFileName, double* grid, int number, int additionalBinNumber, MPI_Comm& cartComm,
                 int* cartCoord, int* cartDim, double scale) {
	int rank;
	int size;
	MPI_Comm_rank(cartComm, &rank);
	MPI_Comm_size(cartComm, &size);
	int minK = 2 + additionalBinNumber;
	if (cartCoord[2] == 0) {
		minK = 0;
	}
	int maxK = number - 1 - additionalBinNumber;
	if (cartCoord[2] == cartDim[2] - 1) {
		maxK = number;
	}

	FILE* outFile = NULL;

	if (rank == 0) {
		outFile = fopen(outFileName, "w");
		fclose(outFile);
	}

	for (int procCount = 0; procCount < cartDim[2]; ++procCount) {
		MPI_Barrier(cartComm);
		if (procCount == cartCoord[2] && cartCoord[0] == 0 && cartCoord[1] == 0) {
			outFile = fopen(outFileName, "a");
			for (int k = minK; k <= maxK; ++k) {
				fprintf(outFile, "%15.10g\n", grid[k] * scale);
			}
			fclose(outFile);
		}
	}
}

void outputGridReducedX(const char* outFileName, double* grid, int number, int additionalBinNumber, int step, int rank,
                        int prevRank, int nextRank, MPI_Comm& cartComm, int* cartCoord, int* cartDim, double scale) {
	int shift[1];
	shift[0] = step - 1;
	if (cartCoord[0] > 0 && cartCoord[1] == 0 && cartCoord[2] == 0) {
		MPI_Status status;
		MPI_Recv(shift, 1, MPI_INT, prevRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
	}
	int minI = 2 + additionalBinNumber;
	if (cartCoord[0] == 0) {
		minI = 0;
	}
	int maxI = number - 1 - additionalBinNumber;
	if (cartCoord[0] == cartDim[0] - 1) {
		maxI = number;
	}

	FILE* outFile = NULL;

	if (rank == 0) {
		outFile = fopen(outFileName, "w");
		fclose(outFile);
	}

	if (cartCoord[1] == 0 && cartCoord[2] == 0) {
		outFile = fopen(outFileName, "a");
		for (int i = minI; i <= maxI; ++i) {
			if ((i + shift[0]) % step == 0) {
				fprintf(outFile, "%15.10g\n", grid[i] * scale);
			}
		}
		fclose(outFile);

		shift[0] = (number - 1 + shift[0]) % step;

		if (cartCoord[0] < cartDim[0] - 1) {
			MPI_Send(shift, 1, MPI_INT, nextRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm);
		}
	}
}

void outputGridReducedY(const char* outFileName, double* grid, int number, int additionalBinNumber, int step, int rank,
                        int prevRank, int nextRank, MPI_Comm& cartComm, int* cartCoord, int* cartDim, double scale) {
	int shift[1];
	shift[0] = step - 1;
	if (cartCoord[1] > 0 && cartCoord[0] == 0 && cartCoord[2] == 0) {
		MPI_Status status;
		MPI_Recv(shift, 1, MPI_INT, prevRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
	}
	int minI = 2 + additionalBinNumber;
	if (cartCoord[1] == 0) {
		minI = 0;
	}
	int maxI = number - 1 - additionalBinNumber;
	if (cartCoord[1] == cartDim[1] - 1) {
		maxI = number;
	}

	FILE* outFile = NULL;

	if (rank == 0) {
		outFile = fopen(outFileName, "w");
		fclose(outFile);
	}

	if (cartCoord[0] == 0 && cartCoord[2] == 0) {
		outFile = fopen(outFileName, "a");
		for (int i = minI; i <= maxI; ++i) {
			if ((i + shift[0]) % step == 0) {
				fprintf(outFile, "%15.10g\n", grid[i] * scale);
			}
		}
		fclose(outFile);

		shift[0] = (number - 1 + shift[0]) % step;

		if (cartCoord[1] < cartDim[1] - 1) {
			MPI_Send(shift, 1, MPI_INT, nextRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm);
		}
	}
}

void outputGridReducedZ(const char* outFileName, double* grid, int number, int additionalBinNumber, int step, int rank,
                        int prevRank, int nextRank, MPI_Comm& cartComm, int* cartCoord, int* cartDim, double scale) {
	int shift[1];
	shift[0] = step - 1;
	if (cartCoord[2] > 0 && cartCoord[1] == 0 && cartCoord[0] == 0) {
		MPI_Status status;
		MPI_Recv(shift, 1, MPI_INT, prevRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm, &status);
	}
	int minI = 2 + additionalBinNumber;
	if (cartCoord[2] == 0) {
		minI = 0;
	}
	int maxI = number - 1 - additionalBinNumber;
	if (cartCoord[2] == cartDim[2] - 1) {
		maxI = number;
	}

	FILE* outFile = NULL;

	if (rank == 0) {
		outFile = fopen(outFileName, "w");
		fclose(outFile);
	}

	if (cartCoord[1] == 0 && cartCoord[0] == 0) {
		outFile = fopen(outFileName, "a");
		for (int i = minI; i <= maxI; ++i) {
			if ((i + shift[0]) % step == 0) {
				fprintf(outFile, "%15.10g\n", grid[i] * scale);
			}
		}
		fclose(outFile);

		shift[0] = (number - 1 + shift[0]) % step;

		if (cartCoord[2] < cartDim[2] - 1) {
			MPI_Send(shift, 1, MPI_INT, nextRank, MPI_SEND_INTEGER_NUMBER_RIGHT, cartComm);
		}
	}
}

void outputFieldsCrossectionYZ(const char* outEfileName, const char* outBfileName, Vector3d*** Efield,
                               Vector3d*** Bfield, int xnumberAdded,
                               int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period,
                               double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim,
                               int xindex) {
	double scale = 1.0 / (plasma_period * sqrt(gyroradius));
	FILE* outFile = NULL;
	if (cartDim[1] == 1 && cartDim[2] == 1) {
		outFile = fopen(outEfileName, "a");
	}
	for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
		int minJ = 3 + 2 * additionalBinNumber;
		if (cartCoord[1] == 0) {
			minJ = 0;
		}
		int maxJ = ynumberAdded;
		if (cartCoord[1] == cartDim[1] - 1) {
			maxJ = ynumberAdded;
		}
		MPI_Barrier(cartComm);
		if (cartJ == cartCoord[1]) {
			if (cartDim[1] > 1 && cartDim[2] == 1) {
				outFile = fopen(outEfileName, "a");
			}
			for (int j = minJ; j <= maxJ; ++j) {
				for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
					int minK = 3 + 2 * additionalBinNumber;
					if (cartCoord[2] == 0) {
						minK = 0;
					}
					int maxK = znumberAdded;
					if (cartCoord[2] == cartDim[2] - 1) {
						maxK = znumberAdded;
					}
					MPI_Barrier(subCommY);
					if (cartK == cartCoord[2]) {
						if (cartDim[2] > 1) {
							outFile = fopen(outEfileName, "a");
						}
						for (int k = minK; k <= maxK; ++k) {
							fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Efield[xindex][j][k].x,
							        scale * Efield[xindex][j][k].y, scale * Efield[xindex][j][k].z);
						}
						if (cartDim[2] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] > 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
	if (cartDim[1] == 1 && cartDim[2] == 1) {
		fclose(outFile);
	}


	if (cartDim[1] == 1 && cartDim[2] == 1) {
		outFile = fopen(outBfileName, "a");
	}
	for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
		int minJ = 2 + 2 * additionalBinNumber;
		if (cartCoord[1] == 0) {
			minJ = 0;
		}
		int maxJ = ynumberAdded - 1;
		if (cartCoord[1] == cartDim[1] - 1) {
			maxJ = ynumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartJ == cartCoord[1]) {
			if (cartDim[1] > 1 && cartDim[2] == 1) {
				outFile = fopen(outBfileName, "a");
			}
			for (int j = minJ; j <= maxJ; ++j) {
				for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
					int minK = 2 + 2 * additionalBinNumber;
					if (cartCoord[2] == 0) {
						minK = 0;
					}
					int maxK = znumberAdded - 1;
					if (cartCoord[2] == cartDim[2] - 1) {
						maxK = znumberAdded - 1;
					}
					MPI_Barrier(subCommY);
					if (cartK == cartCoord[2]) {
						if (cartDim[2] > 1) {
							outFile = fopen(outBfileName, "a");
						}
						for (int k = minK; k <= maxK; ++k) {
							fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Bfield[xindex][j][k].x,
							        scale * Bfield[xindex][j][k].y, scale * Bfield[xindex][j][k].z);
						}
						if (cartDim[2] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] > 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
	if (cartDim[1] == 1 && cartDim[2] == 1) {
		fclose(outFile);
	}
}

void outputFieldsCrossectionXY(const char* outEfileName, const char* outBfileName, Vector3d*** Efield,
                               Vector3d*** Bfield, int xnumberAdded,
                               int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period,
                               double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim,
                               int zindex) {
	double scale = 1.0 / (plasma_period * sqrt(gyroradius));
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 3 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1) {
				outFile = fopen(outEfileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 3 + 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded;
					}
					MPI_Barrier(subCommY);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1) {
							outFile = fopen(outEfileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {

							fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Efield[i][j][zindex].x,
							        scale * Efield[i][j][zindex].y, scale * Efield[i][j][zindex].z);

						}
						if (cartDim[1] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1) {
				fclose(outFile);
			}
		}
	}

	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1) {
				outFile = fopen(outBfileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 2 + 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded - 1;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded - 1;
					}
					MPI_Barrier(subCommY);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1) {
							outFile = fopen(outBfileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Bfield[i][j][zindex].x,
							        scale * Bfield[i][j][zindex].y, scale * Bfield[i][j][zindex].z);
						}
						if (cartDim[1] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputFieldsCrossectionXZ(const char* outEfileName, const char* outBfileName, Vector3d*** Efield,
                               Vector3d*** Bfield, int xnumberAdded,
                               int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period,
                               double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY, int* cartCoord, int* cartDim,
                               int yindex) {
	double scale = 1.0 / (plasma_period * sqrt(gyroradius));
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 3 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[2] == 1) {
				outFile = fopen(outEfileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
					int minK = 3 + 2 * additionalBinNumber;
					if (cartCoord[2] == 0) {
						minK = 0;
					}
					int maxK = znumberAdded;
					if (cartCoord[2] == cartDim[2] - 1) {
						maxK = znumberAdded;
					}
					MPI_Barrier(subCommY);
					if (cartK == cartCoord[2]) {
						if (cartDim[2] > 1) {
							outFile = fopen(outEfileName, "a");
						}
						for (int k = minK; k <= maxK; ++k) {
							fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Efield[i][yindex][k].x,
							        scale * Efield[i][yindex][k].y, scale * Efield[i][yindex][k].z);
						}
						if (cartDim[2] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}

	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[2] == 1) {
				outFile = fopen(outBfileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
					int minK = 2 + 2 * additionalBinNumber;
					if (cartCoord[2] == 0) {
						minK = 0;
					}
					int maxK = znumberAdded - 1;
					if (cartCoord[2] == cartDim[2] - 1) {
						maxK = znumberAdded - 1;
					}
					MPI_Barrier(subCommY);
					if (cartK == cartCoord[2]) {
						if (cartDim[2] > 1) {
							outFile = fopen(outBfileName, "a");
						}
						for (int k = minK; k <= maxK; ++k) {
							fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Bfield[i][yindex][k].x,
							        scale * Bfield[i][yindex][k].y, scale * Bfield[i][yindex][k].z);
						}
						if (cartDim[2] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputFieldsLineX(const char* outEfileName, const char* outBfileName, Vector3d*** Efield, Vector3d*** Bfield,
                       int xnumberAdded,
                       int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period,
                       double gyroradius, MPI_Comm& subCommX, int* cartCoord, int* cartDim, int yindex, int zindex) {
	double scale = 1.0 / (plasma_period * sqrt(gyroradius));
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 3 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded;
		}
		MPI_Barrier(subCommX);
		if (cartCoord[0] == cartI) {
			outFile = fopen(outEfileName, "a");
			for (int i = minI; i <= maxI; ++i) {
				fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Efield[i][yindex][zindex].x,
				        scale * Efield[i][yindex][zindex].y, scale * Efield[i][yindex][zindex].z);
			}
			fclose(outFile);
		}
	}

	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(subCommX);
		if (cartCoord[0] == cartI) {
			outFile = fopen(outBfileName, "a");
			for (int i = minI; i <= maxI; ++i) {
				fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Bfield[i][yindex][zindex].x,
				        scale * Bfield[i][yindex][zindex].y, scale * Bfield[i][yindex][zindex].z);
			}
			fclose(outFile);
		}
	}
}

void outputFieldsLineY(const char* outEfileName, const char* outBfileName, Vector3d*** Efield, Vector3d*** Bfield,
                       int xnumberAdded,
                       int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period,
                       double gyroradius, MPI_Comm& subCommY, int* cartCoord, int* cartDim, int xindex, int zindex) {
	double scale = 1.0 / (plasma_period * sqrt(gyroradius));
	FILE* outFile = NULL;
	for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
		int minJ = 3 + 2 * additionalBinNumber;
		if (cartCoord[1] == 0) {
			minJ = 0;
		}
		int maxJ = ynumberAdded;
		if (cartCoord[1] == cartDim[1] - 1) {
			maxJ = ynumberAdded;
		}
		MPI_Barrier(subCommY);
		if (cartJ == cartCoord[1]) {
			outFile = fopen(outEfileName, "a");
			for (int j = minJ; j <= maxJ; ++j) {
				fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Efield[xindex][j][zindex].x,
				        scale * Efield[xindex][j][zindex].y, scale * Efield[xindex][j][zindex].z);
			}
			fclose(outFile);

		}
	}


	for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
		int minJ = 2 + 2 * additionalBinNumber;
		if (cartCoord[1] == 0) {
			minJ = 0;
		}
		int maxJ = ynumberAdded - 1;
		if (cartCoord[1] == cartDim[1] - 1) {
			maxJ = ynumberAdded - 1;
		}
		MPI_Barrier(subCommY);
		if (cartJ == cartCoord[1]) {
			outFile = fopen(outBfileName, "a");
			for (int j = minJ; j <= maxJ; ++j) {
				fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Bfield[xindex][j][zindex].x,
				        scale * Bfield[xindex][j][zindex].y, scale * Bfield[xindex][j][zindex].z);
			}
			fclose(outFile);
		}
	}
}

void outputFieldsLineZ(const char* outEfileName, const char* outBfileName, Vector3d*** Efield, Vector3d*** Bfield,
                       int xnumberAdded,
                       int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period,
                       double gyroradius, MPI_Comm& subCommZ, int* cartCoord, int* cartDim, int xindex, int yindex) {
	double scale = 1.0 / (plasma_period * sqrt(gyroradius));
	FILE* outFile = NULL;

	for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
		int minK = 3 + 2 * additionalBinNumber;
		if (cartCoord[2] == 0) {
			minK = 0;
		}
		int maxK = znumberAdded;
		if (cartCoord[2] == cartDim[2] - 1) {
			maxK = znumberAdded;
		}
		MPI_Barrier(subCommZ);
		if (cartK == cartCoord[2]) {
			outFile = fopen(outEfileName, "a");
			for (int k = minK; k <= maxK; ++k) {
				fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Efield[xindex][yindex][k].x,
				        scale * Efield[xindex][yindex][k].y, scale * Efield[xindex][yindex][k].z);
			}
			fclose(outFile);

		}
	}


	for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
		int minK = 2 + 2 * additionalBinNumber;
		if (cartCoord[2] == 0) {
			minK = 0;
		}
		int maxK = znumberAdded - 1;
		if (cartCoord[2] == cartDim[2] - 1) {
			maxK = znumberAdded - 1;
		}
		MPI_Barrier(subCommZ);
		if (cartK == cartCoord[2]) {
			outFile = fopen(outBfileName, "a");
			for (int k = minK; k <= maxK; ++k) {
				fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Bfield[xindex][yindex][k].x,
				        scale * Bfield[xindex][yindex][k].y, scale * Bfield[xindex][yindex][k].z);
			}
			fclose(outFile);
		}
	}

}

void outputFields(const char* outEfileName, const char* outBfileName, Vector3d*** Efield, Vector3d*** Bfield,
                  int xnumberAdded,
                  int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period, double gyroradius,
                  MPI_Comm& cartComm, int* cartCoord, int* cartDim) {
	double scale = 1.0 / (plasma_period * sqrt(gyroradius));
	int dims1[3];
	dims1[0] = 0;
	dims1[1] = 1;
	dims1[2] = 1;
	MPI_Comm subCommX;
	MPI_Cart_sub(cartComm, dims1, &subCommX);
	int dims2[3];
	dims2[0] = 0;
	dims2[1] = 0;
	dims2[2] = 1;
	MPI_Comm subCommY;
	MPI_Cart_sub(cartComm, dims2, &subCommY);
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 3 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				outFile = fopen(outEfileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 3 + 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded;
					}
					MPI_Barrier(subCommX);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							outFile = fopen(outEfileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
								int minK = 3 + 2 * additionalBinNumber;
								if (cartCoord[2] == 0) {
									minK = 0;
								}
								int maxK = znumberAdded;
								if (cartCoord[2] == cartDim[2] - 1) {
									maxK = znumberAdded;
								}
								MPI_Barrier(subCommY);
								if (cartK == cartCoord[2]) {
									if (cartDim[2] > 1) {
										outFile = fopen(outEfileName, "a");
									}
									for (int k = minK; k <= maxK; ++k) {
										fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Efield[i][j][k].x,
										        scale * Efield[i][j][k].y, scale * Efield[i][j][k].z);
									}
									if (cartDim[2] > 1) {
										fclose(outFile);
									}
								}
							}
						}
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}

	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				outFile = fopen(outBfileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 2 + 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded - 1;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded - 1;
					}
					MPI_Barrier(subCommX);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							outFile = fopen(outBfileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
								int minK = 2 + 2 * additionalBinNumber;
								if (cartCoord[2] == 0) {
									minK = 0;
								}
								int maxK = znumberAdded - 1;
								if (cartCoord[2] == cartDim[2] - 1) {
									maxK = znumberAdded - 1;
								}
								MPI_Barrier(subCommY);
								if (cartK == cartCoord[2]) {
									if (cartDim[2] > 1) {
										outFile = fopen(outBfileName, "a");
									}
									for (int k = minK; k <= maxK; ++k) {
										fprintf(outFile, "%30.20g %30.20g %30.20g\n", scale * Bfield[i][j][k].x,
										        scale * Bfield[i][j][k].y, scale * Bfield[i][j][k].z);
									}
									if (cartDim[2] > 1) {
										fclose(outFile);
									}
								}
							}
						}
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

/*void outputFieldsReduced(const char* outEfileName, const char* outBfileName, Vector3d*** Efield, Vector3d*** Bfield, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int stepX, int stepY, int stepZ, double plasma_period, double gyroradius) {
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	double scale = 1.0 / (plasma_period * sqrt(gyroradius));
	int shift[1];
	shift[0] = stepX - 1;
	if (rank > 0) {
		MPI_Status status;
		MPI_Recv(shift, 1, MPI_INT, rank - 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD, &status);
	}
	int minI = 1;
	if (rank == 0) {
		minI = 0;
	}
	int maxI = xnumberAdded;
	if (rank == size - 1) {
		maxI = xnumberAdded + 1;
	}

	FILE* outBfile = fopen(outBfileName, "a");
	for (int i = minI; i < maxI; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				if (((i + shift[0]) % stepX == 0) && (j % stepY == 0) && (k % stepZ == 0)) {
					fprintf(outBfile, "%15.10g %15.10g %15.10g\n", scale * Bfield[i][j][k].x,
					        scale * Bfield[i][j][k].y, scale * Bfield[i][j][k].z);
				}
			}
		}
	}
	fclose(outBfile);

	shift[0] = (xnumberAdded - 1 + shift[0]) % stepX;

	if (rank < size - 1) {
		MPI_Send(shift, 1, MPI_INT, rank + 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD);
	}


	MPI_Barrier(MPI_COMM_WORLD);

	shift[0] = stepX - 1;
	if (rank > 0) {
		MPI_Status status;
		MPI_Recv(shift, 1, MPI_INT, rank - 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD, &status);
	}

	minI = 2;
	if (rank == 0) {
		minI = 0;
	}

	FILE* outEfile = fopen(outEfileName, "a");
	for (int i = minI; i <= maxI; ++i) {
		for (int j = 0; j <= ynumberAdded; ++j) {
			for (int k = 0; k <= znumberAdded; ++k) {
				if (((i + shift[0]) % stepX == 0) && (j % stepY == 0) && (k % stepZ == 0)) {
					fprintf(outEfile, "%15.10g %15.10g %15.10g\n", scale * Efield[i][j][k].x,
					        scale * Efield[i][j][k].y, scale * Efield[i][j][k].z);
				}
			}
		}
	}
	fclose(outEfile);


	shift[0] = (xnumberAdded - 1 + shift[0]) % stepX;

	if (rank < size - 1) {
		MPI_Send(shift, 1, MPI_INT, rank + 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD);
	}
}*/

void outputConcentrationsCrossectionXY(const char* outFileName, double**** particleConcentrations,
                                       double*** chargeDensity,
                                       double*** shiftChargeDensity, int xnumberAdded, int ynumberAdded,
                                       int znumberAdded, int additionalBinNumber,
                                       int typesNumber, double plasma_period, double gyroradius, MPI_Comm& cartComm,
                                       MPI_Comm& subCommY, int* cartCoord, int* cartDim, int zindex) {
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		FILE* outFile = NULL;
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 2 + 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded - 1;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded - 1;
					}
					MPI_Barrier(subCommY);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							fprintf(outFile, "%15.10g %15.10g ",
							        chargeDensity[i][j][zindex] / (sqrt(cube(gyroradius)) * plasma_period),
							        shiftChargeDensity[i][j][zindex] / (sqrt(cube(gyroradius)) * plasma_period));
							for (int t = 0; t < typesNumber; ++t) {
								fprintf(outFile, "%15.10g ", particleConcentrations[t][i][j][zindex] / cube(gyroradius));
							}
							fprintf(outFile, "\n");
						}
						if (cartDim[1] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputConcentrationsCrossectionXZ(const char* outFileName, double**** particleConcentrations,
                                       double*** chargeDensity,
                                       double*** shiftChargeDensity, int xnumberAdded, int ynumberAdded,
                                       int znumberAdded, int additionalBinNumber,
                                       int typesNumber, double plasma_period, double gyroradius, MPI_Comm& cartComm,
                                       MPI_Comm& subCommY, int* cartCoord, int* cartDim, int yindex) {
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		FILE* outFile = NULL;
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
					int minK = 2 + 2 * additionalBinNumber;
					if (cartCoord[2] == 0) {
						minK = 0;
					}
					int maxK = znumberAdded - 1;
					if (cartCoord[2] == cartDim[2] - 1) {
						maxK = znumberAdded - 1;
					}
					MPI_Barrier(subCommY);
					if (cartK == cartCoord[2]) {
						if (cartDim[2] > 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int k = minK; k <= maxK; ++k) {
							fprintf(outFile, "%15.10g %15.10g ",
							        chargeDensity[i][yindex][k] / (sqrt(cube(gyroradius)) * plasma_period),
							        shiftChargeDensity[i][yindex][k] / (sqrt(cube(gyroradius)) * plasma_period));
							for (int t = 0; t < typesNumber; ++t) {
								fprintf(outFile, "%15.10g ", particleConcentrations[t][i][yindex][k] / cube(gyroradius));
							}
							fprintf(outFile, "\n");
						}
						if (cartDim[2] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputConcentrationsCrossectionYZ(const char* outFileName, double**** particleConcentrations,
                                       double*** chargeDensity,
                                       double*** shiftChargeDensity, int xnumberAdded, int ynumberAdded,
                                       int znumberAdded, int additionalBinNumber,
                                       int typesNumber, double plasma_period, double gyroradius, MPI_Comm& cartComm,
                                       MPI_Comm& subCommY, int* cartCoord, int* cartDim, int xindex) {

	FILE* outFile = NULL;
	if (cartDim[1] == 1 && cartDim[2] == 1) {
		outFile = fopen(outFileName, "a");
	}
	for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
		int minJ = 2 + 2 * additionalBinNumber;
		if (cartCoord[1] == 0) {
			minJ = 0;
		}
		int maxJ = ynumberAdded - 1;
		if (cartCoord[1] == cartDim[1] - 1) {
			maxJ = ynumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartJ == cartCoord[1]) {
			if (cartDim[1] > 1 && cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int j = minJ; j <= maxJ; ++j) {
				for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
					int minK = 2 + 2 * additionalBinNumber;
					if (cartCoord[2] == 0) {
						minK = 0;
					}
					int maxK = znumberAdded - 1;
					if (cartCoord[2] == cartDim[2] - 1) {
						maxK = znumberAdded - 1;
					}
					MPI_Barrier(subCommY);
					if (cartK == cartCoord[2]) {
						if (cartDim[2] > 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int k = minK; k <= maxK; ++k) {
							fprintf(outFile, "%15.10g %15.10g ",
							        chargeDensity[xindex][j][k] / (sqrt(cube(gyroradius)) * plasma_period),
							        shiftChargeDensity[xindex][j][k] / (sqrt(cube(gyroradius)) * plasma_period));
							for (int t = 0; t < typesNumber; ++t) {
								fprintf(outFile, "%15.10g ", particleConcentrations[t][xindex][j][k] / cube(gyroradius));
							}
							fprintf(outFile, "\n");
						}
						if (cartDim[2] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] > 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
	if (cartDim[1] == 1 && cartDim[2] == 1) {
		fclose(outFile);
	}
}

void outputConcentrationsLineX(const char* outFileName, double**** particleConcentrations, double*** chargeDensity,
                               double*** shiftChargeDensity, int xnumberAdded, int ynumberAdded, int znumberAdded,
                               int additionalBinNumber,
                               int typesNumber, double plasma_period, double gyroradius, MPI_Comm& subCommX,
                               int* cartCoord, int* cartDim, int yindex, int zindex) {
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(subCommX);
		if (cartCoord[0] == cartI) {
			outFile = fopen(outFileName, "a");
			for (int i = minI; i <= maxI; ++i) {
				fprintf(outFile, "%15.10g %15.10g ",
				        chargeDensity[i][yindex][zindex] / (sqrt(cube(gyroradius)) * plasma_period),
				        shiftChargeDensity[i][yindex][zindex] / (sqrt(cube(gyroradius)) * plasma_period));
				for (int t = 0; t < typesNumber; ++t) {
					fprintf(outFile, "%15.10g ", particleConcentrations[t][i][yindex][zindex] / cube(gyroradius));
				}
				fprintf(outFile, "\n");
			}
			fclose(outFile);
		}
	}
}

void outputConcentrationsLineY(const char* outFileName, double**** particleConcentrations, double*** chargeDensity,
                               double*** shiftChargeDensity, int xnumberAdded, int ynumberAdded, int znumberAdded,
                               int additionalBinNumber,
                               int typesNumber, double plasma_period, double gyroradius, MPI_Comm& subCommY,
                               int* cartCoord, int* cartDim, int xindex, int zindex) {
	FILE* outFile = NULL;
	for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
		int minJ = 2 + 2 * additionalBinNumber;
		if (cartCoord[1] == 0) {
			minJ = 0;
		}
		int maxJ = ynumberAdded - 1;
		if (cartCoord[1] == cartDim[1] - 1) {
			maxJ = ynumberAdded - 1;
		}
		MPI_Barrier(subCommY);
		if (cartJ == cartCoord[1]) {
			outFile = fopen(outFileName, "a");
			for (int j = minJ; j <= maxJ; ++j) {
				fprintf(outFile, "%15.10g %15.10g ", chargeDensity[xindex][j][zindex] / (sqrt(cube(gyroradius)) * plasma_period),
				        shiftChargeDensity[xindex][j][zindex] / (sqrt(cube(gyroradius)) * plasma_period));
				for (int t = 0; t < typesNumber; ++t) {
					fprintf(outFile, "%15.10g ", particleConcentrations[t][xindex][j][zindex] / cube(gyroradius));
				}
				fprintf(outFile, "\n");
			}
			fclose(outFile);
		}
	}
}

void outputConcentrationsLineZ(const char* outFileName, double**** particleConcentrations, double*** chargeDensity,
                               double*** shiftChargeDensity, int xnumberAdded, int ynumberAdded, int znumberAdded,
                               int additionalBinNumber,
                               int typesNumber, double plasma_period, double gyroradius, MPI_Comm& subCommZ,
                               int* cartCoord, int* cartDim, int xindex, int yindex) {
	FILE* outFile = NULL;
	for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
		int minK = 2 + 2 * additionalBinNumber;
		if (cartCoord[2] == 0) {
			minK = 0;
		}
		int maxK = znumberAdded - 1;
		if (cartCoord[2] == cartDim[2] - 1) {
			maxK = znumberAdded - 1;
		}
		MPI_Barrier(subCommZ);
		if (cartK == cartCoord[2]) {
			outFile = fopen(outFileName, "a");
			for (int k = minK; k <= maxK; ++k) {
				fprintf(outFile, "%15.10g %15.10g ", chargeDensity[xindex][yindex][k] / (sqrt(cube(gyroradius)) * plasma_period),
				        shiftChargeDensity[xindex][yindex][k] / (sqrt(cube(gyroradius)) * plasma_period));
				for (int t = 0; t < typesNumber; ++t) {
					fprintf(outFile, "%15.10g ", particleConcentrations[t][xindex][yindex][k] / cube(gyroradius));
				}
				fprintf(outFile, "\n");
			}
			fclose(outFile);
		}
	}
}

void outputConcentrations(const char* outFileName, double**** particleConcentrations, double*** chargeDensity,
                          double*** shiftChargeDensity, int xnumberAdded, int ynumberAdded, int znumberAdded,
                          int additionalBinNumber,
                          int typesNumber, double plasma_period, double gyroradius, MPI_Comm& cartComm, int* cartCoord,
                          int* cartDim) {
	int dims1[3];
	dims1[0] = 0;
	dims1[1] = 1;
	dims1[2] = 1;
	MPI_Comm subCommX;
	MPI_Cart_sub(cartComm, dims1, &subCommX);
	int dims2[3];
	dims2[0] = 0;
	dims2[1] = 0;
	dims2[2] = 1;
	MPI_Comm subCommY;
	MPI_Cart_sub(cartComm, dims2, &subCommY);
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 2 + 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded - 1;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded - 1;
					}
					MPI_Barrier(subCommX);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
								int minK = 2 + 2 * additionalBinNumber;
								if (cartCoord[2] == 0) {
									minK = 0;
								}
								int maxK = znumberAdded - 1;
								if (cartCoord[2] == cartDim[2] - 1) {
									maxK = znumberAdded - 1;
								}
								MPI_Barrier(subCommY);
								if (cartK == cartCoord[2]) {
									if (cartDim[2] > 1) {
										outFile = fopen(outFileName, "a");
									}
									for (int k = minK; k <= maxK; ++k) {
										fprintf(outFile, "%15.10g %15.10g ",
										        chargeDensity[i][j][k] / (sqrt(cube(gyroradius)) * plasma_period),
										        shiftChargeDensity[i][j][k] / (sqrt(cube(gyroradius)) * plasma_period));
										for (int t = 0; t < typesNumber; ++t) {
											fprintf(outFile, "%15.10g ", particleConcentrations[t][i][j][k] / cube(gyroradius));
										}
										fprintf(outFile, "\n");
									}
									if (cartDim[2] > 1) {
										fclose(outFile);
									}
								}
							}
						}
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputVelocityCrossectionXY(const char* outFileName, Vector3d**** velocity, ParticleTypeContainer* types,
                                 int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                                 int typesNumber,
                                 double plasma_period, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY,
                                 int* cartCoord, int* cartDim, int zindex) {
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		FILE* outFile = NULL;
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 2 + 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded - 1;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded - 1;
					}
					MPI_Barrier(subCommY);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							for (int t = 0; t < typesNumber; ++t) {
								fprintf(outFile, "%15.10g %15.10g %15.10g ",
								        velocity[t][i][j][zindex].x * gyroradius / plasma_period,
								        velocity[t][i][j][zindex].y * gyroradius / plasma_period,
								        velocity[t][i][j][zindex].z * gyroradius / plasma_period);
							}
							fprintf(outFile, "\n");
						}
						if (cartDim[1] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputVelocityCrossectionXZ(const char* outFileName, Vector3d**** velocity, ParticleTypeContainer* types,
                                 int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                                 int typesNumber,
                                 double plasma_period, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY,
                                 int* cartCoord, int* cartDim, int yindex) {
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		FILE* outFile = NULL;
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
					int minK = 2 + 2 * additionalBinNumber;
					if (cartCoord[2] == 0) {
						minK = 0;
					}
					int maxK = znumberAdded - 1;
					if (cartCoord[2] == cartDim[2] - 1) {
						maxK = znumberAdded - 1;
					}
					MPI_Barrier(subCommY);
					if (cartK == cartCoord[2]) {
						if (cartDim[2] > 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int k = minK; k <= maxK; ++k) {
							for (int t = 0; t < typesNumber; ++t) {
								fprintf(outFile, "%15.10g %15.10g %15.10g ",
								        velocity[t][i][yindex][k].x * gyroradius / plasma_period,
								        velocity[t][i][yindex][k].y * gyroradius / plasma_period,
								        velocity[t][i][yindex][k].z * gyroradius / plasma_period);
							}
							fprintf(outFile, "\n");
						}
						if (cartDim[2] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputVelocityCrossectionYZ(const char* outFileName, Vector3d**** velocity, ParticleTypeContainer* types,
                                 int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber,
                                 int typesNumber,
                                 double plasma_period, double gyroradius, MPI_Comm& cartComm, MPI_Comm& subCommY,
                                 int* cartCoord, int* cartDim, int xindex) {
	FILE* outFile = NULL;
	if (cartDim[1] == 1 && cartDim[2] == 1) {
		outFile = fopen(outFileName, "a");
	}
	for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
		int minJ = 2 + 2 * additionalBinNumber;
		if (cartCoord[1] == 0) {
			minJ = 0;
		}
		int maxJ = ynumberAdded - 1;
		if (cartCoord[1] == cartDim[1] - 1) {
			maxJ = ynumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartJ == cartCoord[1]) {
			if (cartDim[1] > 1 && cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int j = minJ; j <= maxJ; ++j) {
				for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
					int minK = 2 + 2 * additionalBinNumber;
					if (cartCoord[2] == 0) {
						minK = 0;
					}
					int maxK = znumberAdded - 1;
					if (cartCoord[2] == cartDim[2] - 1) {
						maxK = znumberAdded - 1;
					}
					MPI_Barrier(subCommY);
					if (cartK == cartCoord[2]) {
						if (cartDim[2] > 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int k = minK; k <= maxK; ++k) {
							for (int t = 0; t < typesNumber; ++t) {
								fprintf(outFile, "%15.10g %15.10g %15.10g ",
								        velocity[t][xindex][j][k].x * gyroradius / plasma_period,
								        velocity[t][xindex][j][k].y * gyroradius / plasma_period,
								        velocity[t][xindex][j][k].z * gyroradius / plasma_period);
							}
							fprintf(outFile, "\n");
						}
						if (cartDim[2] > 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] > 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
	if (cartDim[1] == 1 && cartDim[2] == 1) {
		fclose(outFile);
	}
}

void outputVelocityLineX(const char* outFileName, Vector3d**** velocity, ParticleTypeContainer* types,
                         int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int typesNumber,
                         double plasma_period, double gyroradius, MPI_Comm& subCommX, int* cartCoord, int* cartDim,
                         int yindex, int zindex) {
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(subCommX);
		if (cartCoord[0] == cartI) {
			outFile = fopen(outFileName, "a");
			for (int i = minI; i <= maxI; ++i) {
				for (int t = 0; t < typesNumber; ++t) {
					fprintf(outFile, "%15.10g %15.10g %15.10g ",
					        velocity[t][i][yindex][zindex].x * gyroradius / plasma_period,
					        velocity[t][i][yindex][zindex].y * gyroradius / plasma_period,
					        velocity[t][i][yindex][zindex].z * gyroradius / plasma_period);
				}
				fprintf(outFile, "\n");
			}
			fclose(outFile);
		}
	}
}

void outputVelocityLineY(const char* outFileName, Vector3d**** velocity, ParticleTypeContainer* types,
                         int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int typesNumber,
                         double plasma_period, double gyroradius, MPI_Comm& subCommY, int* cartCoord, int* cartDim,
                         int xindex, int zindex) {
	FILE* outFile = NULL;
	for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
		int minJ = 2 + 2 * additionalBinNumber;
		if (cartCoord[1] == 0) {
			minJ = 0;
		}
		int maxJ = ynumberAdded - 1;
		if (cartCoord[1] == cartDim[1] - 1) {
			maxJ = ynumberAdded - 1;
		}
		MPI_Barrier(subCommY);
		if (cartJ == cartCoord[1]) {
			outFile = fopen(outFileName, "a");
			for (int j = minJ; j <= maxJ; ++j) {
				for (int t = 0; t < typesNumber; ++t) {
					fprintf(outFile, "%15.10g %15.10g %15.10g ",
					        velocity[t][xindex][j][zindex].x * gyroradius / plasma_period,
					        velocity[t][xindex][j][zindex].y * gyroradius / plasma_period,
					        velocity[t][xindex][j][zindex].z * gyroradius / plasma_period);
				}
				fprintf(outFile, "\n");
			}
			fclose(outFile);
		}
	}
}

void outputVelocityLineZ(const char* outFileName, Vector3d**** velocity, ParticleTypeContainer* types,
                         int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int typesNumber,
                         double plasma_period, double gyroradius, MPI_Comm& subCommZ, int* cartCoord, int* cartDim,
                         int xindex, int yindex) {
	FILE* outFile = NULL;
	for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
		int minK = 2 + 2 * additionalBinNumber;
		if (cartCoord[2] == 0) {
			minK = 0;
		}
		int maxK = znumberAdded - 1;
		if (cartCoord[2] == cartDim[2] - 1) {
			maxK = znumberAdded - 1;
		}
		MPI_Barrier(subCommZ);
		if (cartK == cartCoord[2]) {
			outFile = fopen(outFileName, "a");
			for (int k = minK; k <= maxK; ++k) {
				for (int t = 0; t < typesNumber; ++t) {
					fprintf(outFile, "%15.10g %15.10g %15.10g ",
					        velocity[t][xindex][yindex][k].x * gyroradius / plasma_period,
					        velocity[t][xindex][yindex][k].y * gyroradius / plasma_period,
					        velocity[t][xindex][yindex][k].z * gyroradius / plasma_period);
				}
				fprintf(outFile, "\n");
			}
			fclose(outFile);
		}
	}
}

void outputVelocity(const char* outFileName, Vector3d**** velocity, ParticleTypeContainer* types,
                    int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int typesNumber,
                    double plasma_period, double gyroradius, MPI_Comm& cartComm, int* cartCoord, int* cartDim) {
	int dims1[3];
	dims1[0] = 0;
	dims1[1] = 1;
	dims1[2] = 1;
	MPI_Comm subCommX;
	MPI_Cart_sub(cartComm, dims1, &subCommX);
	int dims2[3];
	dims2[0] = 0;
	dims2[1] = 0;
	dims2[2] = 1;
	MPI_Comm subCommY;
	MPI_Cart_sub(cartComm, dims2, &subCommY);
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 2 + 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded - 1;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded - 1;
					}
					MPI_Barrier(subCommX);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
								int minK = 2 + 2 * additionalBinNumber;
								if (cartCoord[2] == 0) {
									minK = 0;
								}
								int maxK = znumberAdded - 1;
								if (cartCoord[2] == cartDim[2] - 1) {
									maxK = znumberAdded - 1;
								}
								MPI_Barrier(subCommY);
								if (cartK == cartCoord[2]) {
									if (cartDim[2] > 1) {
										outFile = fopen(outFileName, "a");
									}
									for (int k = minK; k <= maxK; ++k) {
										for (int t = 0; t < typesNumber; ++t) {
											fprintf(outFile, "%15.10g %15.10g %15.10g ",
											        velocity[t][i][j][k].x * gyroradius / plasma_period,
											        velocity[t][i][j][k].y * gyroradius / plasma_period,
											        velocity[t][i][j][k].z * gyroradius / plasma_period);
										}
										fprintf(outFile, "\n");
									}
									if (cartDim[2] > 1) {
										fclose(outFile);
									}
								}
							}
						}
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputFlux(const char* outFileName, Vector3d*** electricFlux, Vector3d*** externalElectricFlux, int xnumberAdded,
                int ynumberAdded, int znumberAdded, int additionalBinNumber, double plasma_period, double gyroradius,
                MPI_Comm& cartComm, int* cartCoord, int* cartDim) {

	int dims1[3];
	dims1[0] = 0;
	dims1[1] = 1;
	dims1[2] = 1;
	MPI_Comm subCommX;
	MPI_Cart_sub(cartComm, dims1, &subCommX);
	int dims2[3];
	dims2[0] = 0;
	dims2[1] = 0;
	dims2[2] = 1;
	MPI_Comm subCommY;
	MPI_Cart_sub(cartComm, dims2, &subCommY);
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 3 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 3 + 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded;
					}
					MPI_Barrier(subCommX);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
								int minK = 3 + 2 * additionalBinNumber;
								if (cartCoord[2] == 0) {
									minK = 0;
								}
								int maxK = znumberAdded;
								if (cartCoord[2] == cartDim[2] - 1) {
									maxK = znumberAdded;
								}
								MPI_Barrier(subCommY);
								if (cartK == cartCoord[2]) {
									if (cartDim[2] > 1) {
										outFile = fopen(outFileName, "a");
									}
									for (int k = minK; k <= maxK; ++k) {
										fprintf(outFile, "%15.10g %15.10g %15.10g\n", electricFlux[i][j][k].x,
										        electricFlux[i][j][k].y, electricFlux[i][j][k].z);
									}
									if (cartDim[2] > 1) {
										fclose(outFile);
									}
								}
							}
						}
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputVectorNodeArraySimple(const char* outFileName, Vector3d*** vector3d, int xnumberAdded, int ynumberAdded,
                                 int znumberAdded,
                                 int additionalBinNumber, double scale) {
	FILE* outFile = fopen(outFileName, "w");
	for (int i = 0; i < xnumberAdded + 2; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				fprintf(outFile, "%15.10g %15.10g %15.10g\n", vector3d[i][j][k].x * scale,
				        vector3d[i][j][k].y * scale, vector3d[i][j][k].z * scale);
			}
		}
	}
	fclose(outFile);
}

void outputVectorNodeArray(const char* outFileName, Vector3d*** vector3d, int xnumberAdded, int ynumberAdded,
                           int znumberAdded,
                           int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, double scale) {
	int dims1[3];
	dims1[0] = 0;
	dims1[1] = 1;
	dims1[2] = 1;
	MPI_Comm subCommX;
	MPI_Cart_sub(cartComm, dims1, &subCommX);
	int dims2[3];
	dims2[0] = 0;
	dims2[1] = 0;
	dims2[2] = 1;
	MPI_Comm subCommY;
	MPI_Cart_sub(cartComm, dims2, &subCommY);
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 3 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 3 + 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded;
					}
					MPI_Barrier(subCommX);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
								int minK = 3 + 2 * additionalBinNumber;
								if (cartCoord[2] == 0) {
									minK = 0;
								}
								int maxK = znumberAdded;
								if (cartCoord[2] == cartDim[2] - 1) {
									maxK = znumberAdded;
								}
								MPI_Barrier(subCommY);
								if (cartK == cartCoord[2]) {
									if (cartDim[2] > 1) {
										outFile = fopen(outFileName, "a");
									}
									for (int k = minK; k <= maxK; ++k) {
										fprintf(outFile, "%15.10g %15.10g %15.10g\n", scale * vector3d[i][j][k].x,
										        scale * vector3d[i][j][k].y, scale * vector3d[i][j][k].z);
									}
									if (cartDim[2] > 1) {
										fclose(outFile);
									}
								}
							}
						}
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

/*void outputVectorNodeArrayReduced(const char* outFileName, Vector3d*** vector3d, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int stepX, int stepY, int stepZ, double scale) {
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int shift[1];
	shift[0] = stepX - 1;
	if (rank > 0) {
		MPI_Status status;
		MPI_Recv(shift, 1, MPI_INT, rank - 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD, &status);
	}
	int minI = 2;
	if (rank == 0) {
		minI = 0;
	}
	int maxI = xnumberAdded;
	if (rank == size - 1) {
		maxI = xnumberAdded + 1;
	}

	FILE* outFile = fopen(outFileName, "a");
	for (int i = minI; i < maxI; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				if (((i + shift[0]) % stepX == 0) && (j % stepY == 0) && (k % stepZ == 0)) {
					fprintf(outFile, "%15.10g %15.10g %15.10g\n", vector3d[i][j][k].x * scale,
					        vector3d[i][j][k].y * scale, vector3d[i][j][k].z * scale);
				}
			}
		}
	}
	fclose(outFile);

	shift[0] = (xnumberAdded - 1 + shift[0]) % stepX;

	if (rank < size - 1) {
		MPI_Send(shift, 1, MPI_INT, rank + 1, MPI_SEND_INTEGER_NUMBER_RIGHT, MPI_COMM_WORLD);
	}
}*/

void outputVectorCellArray(const char* outFileName, Vector3d*** vector3d, int xnumberAdded, int ynumberAdded,
                           int znumberAdded,
                           int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, double scale) {
	int dims1[3];
	dims1[0] = 0;
	dims1[1] = 1;
	dims1[2] = 1;
	MPI_Comm subCommX;
	MPI_Cart_sub(cartComm, dims1, &subCommX);
	int dims2[3];
	dims2[0] = 0;
	dims2[1] = 0;
	dims2[2] = 1;
	MPI_Comm subCommY;
	MPI_Cart_sub(cartComm, dims2, &subCommY);
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded - 3;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded - 1;
					}
					MPI_Barrier(subCommX);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
								int minK = 2 + 2 * additionalBinNumber;
								if (cartCoord[2] == 0) {
									minK = 0;
								}
								int maxK = znumberAdded - 1;
								if (cartCoord[2] == cartDim[2] - 1) {
									maxK = znumberAdded - 1;
								}
								MPI_Barrier(subCommY);
								if (cartK == cartCoord[2]) {
									if (cartDim[2] > 1) {
										outFile = fopen(outFileName, "a");
									}
									for (int k = minK; k <= maxK; ++k) {
										fprintf(outFile, "%15.10g %15.10g %15.10g\n", scale * vector3d[i][j][k].x, scale * vector3d[i][j][k].y,
										        scale * vector3d[i][j][k].z);

									}
									if (cartDim[2] > 1) {
										fclose(outFile);
									}
								}
							}
						}
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputVectorCellArray(const char* outFileName, double**** vector3d, int xnumberAdded, int ynumberAdded,
                           int znumberAdded,
                           int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, double scale) {
	int dims1[3];
	dims1[0] = 0;
	dims1[1] = 1;
	dims1[2] = 1;
	MPI_Comm subCommX;
	MPI_Cart_sub(cartComm, dims1, &subCommX);
	int dims2[3];
	dims2[0] = 0;
	dims2[1] = 0;
	dims2[2] = 1;
	MPI_Comm subCommY;
	MPI_Cart_sub(cartComm, dims2, &subCommY);
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded - 3;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded - 1;
					}
					MPI_Barrier(subCommX);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
								int minK = 2 + 2 * additionalBinNumber;
								if (cartCoord[2] == 0) {
									minK = 0;
								}
								int maxK = znumberAdded - 1;
								if (cartCoord[2] == cartDim[2] - 1) {
									maxK = znumberAdded - 1;
								}
								MPI_Barrier(subCommY);
								if (cartK == cartCoord[2]) {
									if (cartDim[2] > 1) {
										outFile = fopen(outFileName, "a");
									}
									for (int k = minK; k <= maxK; ++k) {
										fprintf(outFile, "%15.10g %15.10g %15.10g\n", scale * vector3d[i][j][k][0], scale * vector3d[i][j][k][1],
										        scale * vector3d[i][j][k][2]);

									}
									if (cartDim[2] > 1) {
										fclose(outFile);
									}
								}
							}
						}
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputMatrixArray(const char* outFileName, Matrix3d*** matrix3d, int xnumberAdded, int ynumberAdded,
                       int znumberAdded,
                       int additionalBinNumber, MPI_Comm& cartComm, int* cartCoord, int* cartDim, double scale) {
	int dims1[3];
	dims1[0] = 0;
	dims1[1] = 1;
	dims1[2] = 1;
	MPI_Comm subCommX;
	MPI_Cart_sub(cartComm, dims1, &subCommX);
	int dims2[3];
	dims2[0] = 0;
	dims2[1] = 0;
	dims2[2] = 1;
	MPI_Comm subCommY;
	MPI_Cart_sub(cartComm, dims2, &subCommY);
	FILE* outFile = NULL;
	for (int cartI = 0; cartI < cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = xnumberAdded - 1;
		if (cartCoord[0] == cartDim[0] - 1) {
			maxI = xnumberAdded - 1;
		}
		MPI_Barrier(cartComm);
		if (cartCoord[0] == cartI) {
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < cartDim[1]; ++cartJ) {
					int minJ = 2 * additionalBinNumber;
					if (cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = ynumberAdded - 3;
					if (cartCoord[1] == cartDim[1] - 1) {
						maxJ = ynumberAdded - 1;
					}
					MPI_Barrier(subCommX);
					if (cartJ == cartCoord[1]) {
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							for (int cartK = 0; cartK < cartDim[2]; ++cartK) {
								int minK = 2 + 2 * additionalBinNumber;
								if (cartCoord[2] == 0) {
									minK = 0;
								}
								int maxK = znumberAdded - 1;
								if (cartCoord[2] == cartDim[2] - 1) {
									maxK = znumberAdded - 1;
								}
								MPI_Barrier(subCommY);
								if (cartK == cartCoord[2]) {
									if (cartDim[2] > 1) {
										outFile = fopen(outFileName, "a");
									}
									for (int k = minK; k <= maxK; ++k) {
										fprintf(outFile, "%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",
										        scale * matrix3d[i][j][k].matrix[0][0],
										        scale * matrix3d[i][j][k].matrix[0][1], scale * matrix3d[i][j][k].matrix[0][2],
										        scale * matrix3d[i][j][k].matrix[1][0], scale * matrix3d[i][j][k].matrix[1][1],
										        scale * matrix3d[i][j][k].matrix[1][2], scale * matrix3d[i][j][k].matrix[2][0],
										        scale * matrix3d[i][j][k].matrix[2][1], scale * matrix3d[i][j][k].matrix[2][2]);
									}
									if (cartDim[2] > 1) {
										fclose(outFile);
									}
								}
							}
						}
						if (cartDim[1] > 1 && cartDim[2] == 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (cartDim[1] == 1 && cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputGeneral(const char* outFileName, Simulation* simulation) {
	FILE* outFile = fopen(outFileName, "a");
	int particlesCount = simulation->particlesNumber;
	double fieldFactor = simulation->plasma_period * sqrt(simulation->scaleFactor);
	double omega2 = 0;
	for (int i = 0; i < simulation->typesNumber; ++i) {
		omega2 += 4 * pi * simulation->types[i].concentration * simulation->types[i].charge * simulation->types[i].charge /
			simulation->types[i].mass;
	}

	double gamma = 1.0 / sqrt(1 - simulation->V0.scalarMult(simulation->V0) / simulation->speed_of_light_normalized_sqr);
	double omega = sqrt(omega2 / (gamma * gamma * gamma));
	fprintf(outFile,
	        "%d %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %d %15.10g %15.10g %15.10g %15.10g ",
	        simulation->currentIteration, simulation->time, simulation->time * simulation->plasma_period,
	        simulation->particleEnergy,
	        simulation->electricFieldEnergy, simulation->magneticFieldEnergy, simulation->energy,
	        simulation->globalMomentum.x, simulation->globalMomentum.y, simulation->globalMomentum.z,
	        simulation->generalTheoreticalEnergy, simulation->generalTheoreticalMomentum.x,
	        simulation->generalTheoreticalMomentum.y,
	        simulation->generalTheoreticalMomentum.z, simulation->maxEfield.norm() / fieldFactor,
	        simulation->maxBfield.norm() / fieldFactor, simulation->deltaT, particlesCount, simulation->shockWaveX,
	        simulation->meanSquaredEfield[0] / fieldFactor, simulation->meanSquaredEfield[1] / fieldFactor,
	        simulation->meanSquaredEfield[2] / fieldFactor);
	fprintf(outFile,
	        "%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %d %15.10g %d %15.10g\n",
	        simulation->electricFieldEnergyX, simulation->electricFieldEnergyY, simulation->electricFieldEnergyZ,
	        simulation->magneticFieldEnergyX, simulation->magneticFieldEnergyY, simulation->magneticFieldEnergyZ,
	        simulation->electromagneticMomentum.x, simulation->electromagneticMomentum.y,
	        simulation->electromagneticMomentum.z, simulation->particleMomentum.x, simulation->particleMomentum.y,
	        simulation->particleMomentum.z, simulation->derExPoint,
	        simulation->derExPoint * simulation->deltaX * simulation->scaleFactor, simulation->constMeanElevelPoint,
	        simulation->derConcentrationPoint * simulation->deltaX * simulation->scaleFactor);
	fclose(outFile);
}

void outputGeneralAnisotropy(const char* outFileName, Simulation* simulation) {
	FILE* outFile = fopen(outFileName, "a");
	fprintf(outFile, "%d %15.10g %15.10g ",
	        simulation->currentIteration, simulation->time, simulation->time * simulation->plasma_period);
	for (int t = 0; t < simulation->typesNumber; ++ t) {
		fprintf(outFile, "%15.10g %15.10g %15.10g ", simulation->types[t].anisotropy,
		        simulation->types[t].parallelTemperatureEvaluated, simulation->types[t].normalTemperatureEvaluated);
	}
	fprintf(outFile, "\n");
	fclose(outFile);
}

void outputDivergenceError(const char* outFileName, Simulation* simulation, double plasma_period, double gyroradius) {
	int dims1[3];
	dims1[0] = 0;
	dims1[1] = 1;
	dims1[2] = 1;
	MPI_Comm subCommX;
	MPI_Cart_sub(simulation->cartComm, dims1, &subCommX);
	int dims2[3];
	dims2[0] = 0;
	dims2[1] = 0;
	dims2[2] = 1;
	MPI_Comm subCommY;
	MPI_Cart_sub(simulation->cartComm, dims2, &subCommY);
	FILE* outFile = NULL;
	double scale = (sqrt(cube(gyroradius))) * plasma_period;
	for (int cartI = 0; cartI < simulation->cartDim[0]; ++cartI) {
		int minI = 2 + 2 * additionalBinNumber;
		if (simulation->cartCoord[0] == 0) {
			minI = 0;
		}
		int maxI = simulation->xnumberAdded - 1;
		if (simulation->cartCoord[0] == simulation->cartDim[0] - 1) {
			maxI = simulation->xnumberAdded - 1;
		}
		MPI_Barrier(simulation->cartComm);
		if (simulation->cartCoord[0] == cartI) {
			if (simulation->cartDim[1] == 1 && simulation->cartDim[2] == 1) {
				outFile = fopen(outFileName, "a");
			}
			for (int i = minI; i <= maxI; ++i) {
				for (int cartJ = 0; cartJ < simulation->cartDim[1]; ++cartJ) {
					int minJ = 2 + 2 * additionalBinNumber;
					if (simulation->cartCoord[1] == 0) {
						minJ = 0;
					}
					int maxJ = simulation->ynumberAdded - 1;
					if (simulation->cartCoord[1] == simulation->cartDim[1] - 1) {
						maxJ = simulation->ynumberAdded - 1;
					}
					MPI_Barrier(subCommX);
					if (cartJ == simulation->cartCoord[1]) {
						if (simulation->cartDim[1] > 1 && simulation->cartDim[2] == 1) {
							outFile = fopen(outFileName, "a");
						}
						for (int j = minJ; j <= maxJ; ++j) {
							for (int cartK = 0; cartK < simulation->cartDim[2]; ++cartK) {
								int minK = 2 + 2 * additionalBinNumber;
								if (simulation->cartCoord[2] == 0) {
									minK = 0;
								}
								int maxK = simulation->znumberAdded - 1;
								if (simulation->cartCoord[2] == simulation->cartDim[2] - 1) {
									maxK = simulation->znumberAdded - 1;
								}
								MPI_Barrier(subCommY);
								if (cartK == simulation->cartCoord[2]) {
									if (simulation->cartDim[2] > 1) {
										outFile = fopen(outFileName, "a");
									}
									for (int k = minK; k <= maxK; ++k) {
										if (simulation->solverType == BUNEMAN) {
											double div = 0;
											if (i > 0 && j > 0 && k > 0) {
												div = simulation->evaluateDivBunemanE(i, j, k);
											}
											double divB = simulation->evaluateDivBunemanB(i, j, k);
											fprintf(outFile, "%g %g %g %g\n", (4 * pi * simulation->bunemanChargeDensity[i][j][k] - div) / scale,
											        div / scale,
											        4 * pi * simulation->bunemanChargeDensity[i][j][k] / scale, divB / scale);
										} else {
											double div = simulation->evaluateDivE(i, j, k);
											double divB = 0;
											if (i > 0 && j > 0 && k > 0) {
												divB = simulation->evaluateDivB(i, j, k);
											}
											fprintf(outFile, "%g %g %g %g\n", (4 * pi * simulation->chargeDensity[i][j][k] - div) / scale, div / scale,
											        4 * pi * simulation->chargeDensity[i][j][k] / scale, divB / scale);
										}
									}
									if (simulation->cartDim[2] > 1) {
										fclose(outFile);
									}
								}
							}
						}
						if (simulation->cartDim[1] > 1 && simulation->cartDim[2] == 1) {
							fclose(outFile);
						}
					}
				}
			}
			if (simulation->cartDim[1] == 1 && simulation->cartDim[2] == 1) {
				fclose(outFile);
			}
		}
	}
}

void outputParticles(const char* outFileName, Simulation* simulation, ParticleTypes type) {
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	FILE* outFile = NULL;

	if (rank == 0) {
		outFile = fopen(outFileName, "w");
		fclose(outFile);
	}

	for (int procCount = 0; procCount < size; ++procCount) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (procCount == rank) {
			outFile = fopen(outFileName, "a");

			int typeCount = 0;

			for (int i = 0; i < simulation->particles.size(); ++i) {
				Particle* particle = simulation->particles[i];
				Vector3d momentum = particle->getMomentum();
				double p = momentum.norm() * simulation->scaleFactor / simulation->plasma_period;
				if (particle->type == type) {
					if (typeCount % writeParticleNumber == 0) {
						fprintf(outFile, "%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %d\n",
						        particle->coordinates.x * simulation->scaleFactor,
						        momentum.x * simulation->scaleFactor / simulation->plasma_period,
						        particle->coordinates.y * simulation->scaleFactor,
						        momentum.y * simulation->scaleFactor / simulation->plasma_period,
						        particle->coordinates.z * simulation->scaleFactor,
						        momentum.z * simulation->scaleFactor / simulation->plasma_period,
						        p, particle->number);
					}
					typeCount++;
				}
			}
			fclose(outFile);
		}
	}
}

void outputAcceleratedParticlesNumbers(const char* outFileName, Simulation* simulation) {
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	FILE* outFile = NULL;

	if (rank == 0) {
		outFile = fopen(outFileName, "w");

		for (int t = 0; t < simulation->typesNumber; ++t) {
			if (simulation->types[t].particlesPerBin > 0) {
				std::list < std::pair < int, double > >::iterator it = simulation->mostAcceleratedParticlesNumbers[t].begin();
				for (int i = 0; i < mostFastParticlesNumber; ++i) {
					std::pair < int, double > pair = *it;
					++it;
					fprintf(outFile, "%d %15.10g %d\n", pair.first, pair.second, t);
				}
			}
		}

		fclose(outFile);

	}
}

void outputSimulationBackup(const char* generalFileName, const char* EfileName, const char* BfileName,
                            const char* particlesFileName, Simulation* simulation) {

	if (simulation->rank == 0) {
		FILE* generalFile = fopen(generalFileName, "w");
		FILE* Efile = fopen(EfileName, "w");
		FILE* Bfile = fopen(BfileName, "w");
		fclose(Efile);
		fclose(Bfile);

		int inputType = 2;
		if (simulation->inputType == CGS) {
			inputType = 0;
		} else if (simulation->inputType == Theoretical) {
			inputType = 1;
		}
		fprintf(generalFile, "%d\n", inputType);
		fprintf(generalFile, "%d\n", simulation->xnumberGeneral);
		fprintf(generalFile, "%d\n", simulation->ynumberGeneral);
		fprintf(generalFile, "%d\n", simulation->znumberGeneral);
		fprintf(generalFile, "%d\n", simulation->particlesNumber);
		fprintf(generalFile, "%d\n", simulation->types[0].particlesPerBin);
		fprintf(generalFile, "%d\n", simulation->types[1].particlesPerBin);
		fprintf(generalFile, "%d\n", simulation->types[2].particlesPerBin);
		fprintf(generalFile, "%d\n", simulation->types[3].particlesPerBin);
		fprintf(generalFile, "%d\n", simulation->types[4].particlesPerBin);
		fprintf(generalFile, "%d\n", simulation->types[5].particlesPerBin);

		fprintf(generalFile, "%lf\n", simulation->types[0].concentration);
		fprintf(generalFile, "%lf\n", simulation->types[1].concentration);
		fprintf(generalFile, "%lf\n", simulation->types[2].concentration);
		fprintf(generalFile, "%lf\n", simulation->types[3].concentration);
		fprintf(generalFile, "%lf\n", simulation->types[4].concentration);
		fprintf(generalFile, "%lf\n", simulation->types[5].concentration);

		fprintf(generalFile, "%15.10g\n", simulation->temperature);
		fprintf(generalFile, "%15.10g\n", simulation->plasma_period);
		fprintf(generalFile, "%15.10g\n", simulation->scaleFactor);

		fprintf(generalFile, "%15.10g\n", simulation->time);
		fprintf(generalFile, "%15.10g\n", simulation->maxTime);

		fprintf(generalFile, "%d\n", simulation->currentIteration);
		fprintf(generalFile, "%d\n", simulation->maxIteration);

		fprintf(generalFile, "%15.10g\n", simulation->xsizeGeneral);
		fprintf(generalFile, "%15.10g\n", simulation->ysizeGeneral);
		fprintf(generalFile, "%15.10g\n", simulation->zsizeGeneral);
		fprintf(generalFile, "%15.10g\n", simulation->theta);
		fprintf(generalFile, "%15.10g\n", simulation->eta);

		int debugMode = 0;
		if (simulation->debugMode) {
			debugMode = 1;
		} else {
			debugMode = 0;
		}

		fprintf(generalFile, "%d\n", debugMode);

		/*int preserveChargeLocal = 0;
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

		fprintf(generalFile, "%d\n", preserveChargeGlobal);*/

		int solverType = 0;
		if (simulation->solverType == IMPLICIT) {
			solverType = 1;
		} else {
			solverType = 0;
		}

		fprintf(generalFile, "%d\n", solverType);

		int boundaryConditionType = 0;
		if (simulation->boundaryConditionTypeX == PERIODIC) {
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

		fclose(generalFile);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	outputFields(EfileName, BfileName, simulation->Efield, simulation->Bfield, simulation->xnumberAdded,
	             simulation->ynumberAdded,
	             simulation->znumberAdded, additionalBinNumber, 1.0, 1.0, simulation->cartComm, simulation->cartCoord,
	             simulation->cartDim);
	outputBackupParticles(particlesFileName, simulation);
}

void outputBackupParticles(const char* outFileName, Simulation* simulation) {
	int size;
	int rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	for (int curRank = 0; curRank < size; ++curRank) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == curRank) {
			FILE* outFile = NULL;
			if (rank == 0) {
				outFile = fopen(outFileName, "w");
			} else {
				outFile = fopen(outFileName, "a");
			}
			for (int i = 0; i < simulation->particles.size(); ++i) {
				Particle* particle = simulation->particles[i];
				outputBackupParticle(outFile, particle);
			}
			fclose(outFile);
		}
	}
}

void outputBackupParticle(FILE* outFile, Particle* particle) {
	fprintf(outFile, "%d ", particle->number);
	fprintf(outFile, "%15.10g ", particle->mass);
	fprintf(outFile, "%15.10g ", particle->charge);
	fprintf(outFile, "%15.10g ", particle->weight);
	int type;
	if (particle->type == ELECTRON) {
		type = 0;
	} else if (particle->type == PROTON) {
		type = 1;
	} else if (particle->type == POSITRON) {
		type = 2;
	} else if (particle->type == ALPHA) {
		type = 3;
	}
	fprintf(outFile, "%d ", type);
	fprintf(outFile, "%15.10g ", particle->coordinates.x);
	fprintf(outFile, "%15.10g ", particle->coordinates.y);
	fprintf(outFile, "%15.10g ", particle->coordinates.z);
	Vector3d momentum = particle->getMomentum();
	fprintf(outFile, "%15.10g ", momentum.x);
	fprintf(outFile, "%15.10g ", momentum.y);
	fprintf(outFile, "%15.10g ", momentum.z);

	fprintf(outFile, "%15.10g\n", particle->dx);
}

void outputGeneralInitialParameters(const char* outFileName, const char* outFileNameWithText, Simulation* simulation) {
	if (simulation->rank == 0) {
		FILE* outFile = fopen(outFileName, "w");
		FILE* outFileWithText = fopen(outFileNameWithText, "w");
		fprintf(outFileWithText, "1 X size general = %g\n", simulation->xsizeGeneral * simulation->scaleFactor);
		fprintf(outFileWithText, "2 Y size general = %g\n", simulation->ysizeGeneral * simulation->scaleFactor);
		fprintf(outFileWithText, "3 Z size general = %g\n", simulation->zsizeGeneral * simulation->scaleFactor);
		fprintf(outFile, "%g\n", simulation->xsizeGeneral * simulation->scaleFactor);
		fprintf(outFile, "%g\n", simulation->ysizeGeneral * simulation->scaleFactor);
		fprintf(outFile, "%g\n", simulation->zsizeGeneral * simulation->scaleFactor);

		fprintf(outFileWithText, "4 X number general = %d\n", simulation->xnumberGeneral);
		fprintf(outFileWithText, "5 Y number general = %d\n", simulation->ynumberGeneral);
		fprintf(outFileWithText, "6 Z number general = %d\n", simulation->znumberGeneral);
		fprintf(outFile, "%d\n", simulation->xnumberGeneral);
		fprintf(outFile, "%d\n", simulation->ynumberGeneral);
		fprintf(outFile, "%d\n", simulation->znumberGeneral);

		fprintf(outFileWithText, "7 dx = %g\n", simulation->deltaX * simulation->scaleFactor);
		fprintf(outFileWithText, "8 dy = %g\n", simulation->deltaY * simulation->scaleFactor);
		fprintf(outFileWithText, "9 dz = %g\n", simulation->deltaZ * simulation->scaleFactor);
		fprintf(outFile, "%g\n", simulation->deltaX * simulation->scaleFactor);
		fprintf(outFile, "%g\n", simulation->deltaY * simulation->scaleFactor);
		fprintf(outFile, "%g\n", simulation->deltaZ * simulation->scaleFactor);

		fprintf(outFileWithText, "10 c = %g\n", speed_of_light);
		fprintf(outFile, "%g\n", speed_of_light);

		fprintf(outFileWithText, "11 Vx = %g\n", simulation->V0.x * simulation->scaleFactor / simulation->plasma_period);
		fprintf(outFileWithText, "12 Vy = %g\n", simulation->V0.y * simulation->scaleFactor / simulation->plasma_period);
		fprintf(outFileWithText, "13 Vz = %g\n", simulation->V0.z * simulation->scaleFactor / simulation->plasma_period);
		fprintf(outFile, "%g\n", simulation->V0.x * simulation->scaleFactor / simulation->plasma_period);
		fprintf(outFile, "%g\n", simulation->V0.y * simulation->scaleFactor / simulation->plasma_period);
		fprintf(outFile, "%g\n", simulation->V0.z * simulation->scaleFactor / simulation->plasma_period);

		double fieldScale = 1.0 / (simulation->plasma_period * sqrt(simulation->scaleFactor));
		fprintf(outFileWithText, "14 Ex = %g\n", simulation->E0.x * fieldScale);
		fprintf(outFileWithText, "15 Ey = %g\n", simulation->E0.y * fieldScale);
		fprintf(outFileWithText, "16 Ez = %g\n", simulation->E0.z * fieldScale);
		fprintf(outFile, "%g\n", simulation->E0.x * fieldScale);
		fprintf(outFile, "%g\n", simulation->E0.y * fieldScale);
		fprintf(outFile, "%g\n", simulation->E0.z * fieldScale);

		fprintf(outFileWithText, "17 Bx = %g\n", simulation->B0.x * fieldScale);
		fprintf(outFileWithText, "18 By = %g\n", simulation->B0.y * fieldScale);
		fprintf(outFileWithText, "19 Bz = %g\n", simulation->B0.z * fieldScale);
		fprintf(outFile, "%g\n", simulation->B0.x * fieldScale);
		fprintf(outFile, "%g\n", simulation->B0.y * fieldScale);
		fprintf(outFile, "%g\n", simulation->B0.z * fieldScale);

		Vector3d V = simulation->V0 * (simulation->scaleFactor / simulation->plasma_period);
		Vector3d B = simulation->B0;
		double gamma = 1 / sqrt(1 - V.scalarMult(V) / (speed_of_light * speed_of_light));

		double omegaPlasmaElectron = sqrt(
			4 * pi * simulation->electron_charge_normalized * simulation->electron_charge_normalized * simulation->types[0].
			concentration / (simulation->types[0].mass * (gamma))) / simulation->plasma_period;
		fprintf(outFileWithText, "20 plasma electron frequency relativistic = %g\n", omegaPlasmaElectron);
		fprintf(outFile, "%g\n", omegaPlasmaElectron);
		double omega2 = 0;
		for (int i = 0; i < simulation->typesNumber; ++i) {
			omega2 += 4 * pi * simulation->types[i].charge * simulation->types[i].charge * simulation->types[i].concentration /
				simulation->types[i].mass;
		}
		double omegaPlasma = sqrt(omega2 / (gamma)) / simulation->plasma_period;
		fprintf(outFileWithText, "21 plasma frequency relativistic = %g\n", omegaPlasma);
		fprintf(outFile, "%g\n", omegaPlasma);

		double d2 = 0;
		for (int i = 0; i < simulation->typesNumber; ++i) {
			d2 += 4 * pi * simulation->types[i].charge * simulation->types[i].charge * simulation->types[i].concentration / (
				simulation->kBoltzman_normalized * simulation->types[i].temperatureX);
		}
		double debyeLength = simulation->scaleFactor / sqrt(d2);
		fprintf(outFileWithText, "22 debye length = %g\n", debyeLength);
		fprintf(outFile, "%g\n", debyeLength);


		//todo relativistic temperature
		double omegaGyroElectron = (simulation->types[0].charge * simulation->B0.norm() / (simulation->types[0].mass *
			simulation->speed_of_light_normalized)) / simulation->plasma_period;
		double omegaGyroProton = (simulation->types[1].charge * simulation->B0.norm() / (simulation->types[1].mass *
			simulation->speed_of_light_normalized)) / simulation->plasma_period;
		//double gyroRadiusElectron = ((simulation->types[0].mass*simulation->V0.norm()*gamma + sqrt(simulation->types[0].mass*simulation->types[0].temperatureX*simulation->kBoltzman_normalized))*simulation->speed_of_light_normalized/(simulation->types[0].charge*simulation->B0.norm()))*simulation->scaleFactor;
		//double gyroRadiusProton = ((simulation->types[1].mass*simulation->V0.norm()*gamma + sqrt(simulation->types[1].mass*simulation->types[1].temperatureX*simulation->kBoltzman_normalized))*simulation->speed_of_light_normalized/(simulation->types[1].charge*simulation->B0.norm()))*simulation->scaleFactor;
		double gyroRadiusElectron = ((simulation->types[0].mass * simulation->V0.norm() * gamma) * simulation->
			speed_of_light_normalized / (simulation->types[0].charge * simulation->B0.norm())) * simulation->scaleFactor;
		double gyroRadiusProton = ((simulation->types[1].mass * simulation->V0.norm() * gamma) * simulation->
			speed_of_light_normalized / (simulation->types[1].charge * simulation->B0.norm())) * simulation->scaleFactor;

		fprintf(outFileWithText, "23 electron gyro frequency = %g\n", omegaGyroElectron);
		fprintf(outFileWithText, "24 proton gyro frequency = %g\n", omegaGyroProton);
		fprintf(outFileWithText, "25 electron gyro radius = %g\n", gyroRadiusElectron);
		fprintf(outFileWithText, "26 proton gyro radius = %g\n", gyroRadiusProton);
		fprintf(outFile, "%g\n", omegaGyroElectron);
		fprintf(outFile, "%g\n", omegaGyroProton);
		fprintf(outFile, "%g\n", gyroRadiusElectron);
		fprintf(outFile, "%g\n", gyroRadiusProton);

		double electronSkinDepth = speed_of_light / omegaPlasmaElectron;
		fprintf(outFileWithText, "27 electron skin depth = %g\n", electronSkinDepth);
		fprintf(outFile, "%g\n", electronSkinDepth);

		double density = 0;

		for (int i = 0; i < simulation->typesNumber; ++i) {
			density += simulation->types[i].concentration * simulation->types[i].mass;
		}

		double alfvenV = simulation->B0.norm() / sqrt(4 * pi * density);

		double alfvenMach = simulation->V0.norm() / alfvenV;

		fprintf(outFileWithText, "28 alfven Mach = %g\n", alfvenMach);
		fprintf(outFile, "%g\n", alfvenMach);

		fprintf(outFileWithText, "29 nprocs = %d\n", simulation->nprocs);
		fprintf(outFile, "%d\n", simulation->nprocs);

		double theta = atan2(sqrt(simulation->B0.y * simulation->B0.y + simulation->B0.z * simulation->B0.z),
		                     simulation->B0.x) * 180 / pi;
		fprintf(outFileWithText, "30 flow gamma = %g\n", gamma);
		fprintf(outFile, "%g\n", gamma);

		fprintf(outFileWithText, "31 magnetic field angle = %g\n", theta);
		fprintf(outFile, "%g\n", theta);

		Vector3d V0 = simulation->V0;
		//double gamma = 1.0/sqrt(1 - V0.scalarMult(V0)/simulation->speed_of_light_normalized_sqr);

		double sigma = simulation->B0.scalarMult(simulation->B0) / (simulation->density * simulation->
			speed_of_light_normalized_sqr * gamma * 4 * pi);

		fprintf(outFileWithText, "32 magnetization relativistic = %g\n", sigma);
		fprintf(outFile, "%g\n", sigma);

		double beta0 = fabs(V0.x) / simulation->speed_of_light_normalized;
		double betaShock = (sqrt(9 + 16 * beta0 * beta0) - 3) / (8 * beta0);
		if (beta0 == 0) {
			betaShock = 0;
		}
		double Vshock = betaShock * simulation->speed_of_light_normalized;

		fprintf(outFileWithText, "33 supposed shock velocity = %g\n", Vshock);
		fprintf(outFile, "%g\n", Vshock);

		double gammaShock = 1 / sqrt(1 - betaShock * betaShock);

		double thetaCritical = atan(1.0 / (gammaShock * (beta0 + betaShock))) * 180 / pi;
		double thetaCriticalUpstreamSystem = acos(betaShock) * 180 / pi;

		fprintf(outFileWithText, "34 critical theta = %g\n", thetaCritical);
		fprintf(outFile, "%g\n", thetaCritical);

		fprintf(outFileWithText, "35 critical theta in upstream system = %g\n", thetaCriticalUpstreamSystem);
		fprintf(outFile, "%g\n", thetaCriticalUpstreamSystem);

		fprintf(outFileWithText, "36 electron simulation mass = %g\n", simulation->massElectron);
		fprintf(outFile, "%g\n", simulation->massElectron);

		fclose(outFile);
		fclose(outFileWithText);
	}
}

void outputMemory(const char* outFileName, MPI_Comm& cartComm, int* cartCoord, int* cartDim) {
	if (cartCoord[0] == 0 && cartCoord[1] == 0 && cartCoord[2] == 0) {
		FILE* outFile = fopen(outFileName, "w");
		fclose(outFile);
	}
	for (int i = 0; i < cartDim[0]; ++i) {
		for (int j = 0; j < cartDim[1]; ++j) {
			for (int k = 0; k < cartDim[2]; ++k) {
				if ((i == cartCoord[0]) && (j == cartCoord[1]) && (k == cartCoord[2])) {
					FILE* outFile = fopen(outFileName, "a");
					double vm_usage;
					double resident_set;
					process_mem_usage(vm_usage, resident_set);
					fprintf(outFile, "process rank = %d %d %d\n", i, j, k);
					fprintf(outFile, "VM = %g RSS = %g\n", vm_usage, resident_set);
					fclose(outFile);
				}
				MPI_Barrier(cartComm);
			}
		}
	}
}
