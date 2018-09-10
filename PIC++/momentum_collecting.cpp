#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include <time.h>
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "mpi_util.h"
#include "util.h"
#include "output.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "particle.h"
#include "random.h"
#include "simulation.h"
#include "paths.h"

void Simulation::tristanUpdateFlux() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	//MPI_Barrier(cartComm);
	if ((rank == 0) && (verbosity > 1)) printf("reseting tristan flux\n");
	if ((rank == 0) && (verbosity > 1)) printLog("reseting tristan flux\n");

	int maxI = xnumberAdded;
	int maxJ = ynumberAdded;
	int maxK = znumberAdded;

	for (int i = 0; i <= maxI; ++i) {
		for (int j = 0; j <= maxJ; ++j) {
			for (int k = 0; k <= maxK; ++k) {
				if(i != maxI){
					bunemanJx[i][j][k] = 0;
				}
				if(j != maxJ) {
					bunemanJy[i][j][k] = 0;
				}
				if(k != maxK) {
					bunemanJz[i][j][k] = 0;
				}
			}
		}
	}
	/*for (int i = 0; i <= maxI; ++i) {
		for (int j = 0; j < maxJ; ++j) {
			for (int k = 0; k <= maxK; ++k) {
				bunemanJy[i][j][k] = 0;
			}
		}
	}
	for (int i = 0; i <= maxI; ++i) {
		for (int j = 0; j <= maxJ; ++j) {
			for (int k = 0; k < maxK; ++k) {
				bunemanJz[i][j][k] = 0;
			}
		}
	}*/

	//MPI_Barrier(cartComm);

	if ((rank == 0) && (verbosity > 0)) printf("updating tristan flux\n");
	if ((rank == 0) && (verbosity > 0)) printLog("updating tristan flux\n");

	for (int i = 0; i < particles.size(); ++i) {
		Particle* particle = particles[i];
	//for (auto particle : particles) {
		addParticleFluxZigzag(particle);
	}
	//MPI_Barrier(cartComm);
	if(verbosity > 2) printf("exchanging buneman flux rank = %d\n", rank);

	exchangeBunemanFlux();
	MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("evaluating electric flux time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::addParticleFluxZigzag(Particle* particle) {
	Vector3d velocity = particle->getVelocity();
	double fullChargeDensity = particle->charge*particle->weight/cellVolume;
	Vector3d localCoordinates = particle->coordinates - Vector3d(xgrid[0], ygrid[0], zgrid[0]);
	const int i = floor(localCoordinates.x/deltaX);
	const int j = floor(localCoordinates.y/deltaY);
	const int k = floor(localCoordinates.z/deltaZ);
	Vector3d prevLocalCoordinates = localCoordinates - velocity*deltaT;
	//really?
	if((cartCoord[0] == 0) && (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) && (i == reflectingWallPoint)) {
		if(velocity.x > 0) {
			if(prevLocalCoordinates.x < reflectingWallPoint*deltaX) {
				prevLocalCoordinates.x = 2*reflectingWallPoint*deltaX - prevLocalCoordinates.x;
			}
		}
	}
	const int previ = floor(prevLocalCoordinates.x/deltaX);
	const int prevj = floor(prevLocalCoordinates.y/deltaY);
	const int prevk = floor(prevLocalCoordinates.z/deltaZ);

	/*if(i < 0) {
		printf("particle i < 0\n");
		MPI_Finalize();
		exit(0);
	}
	if(i >= xnumberAdded) {
		printf("particle i >= xnumberAdded\n");
		MPI_Finalize();
		exit(0);
	}
	if(previ < 0) {
		printf("particle previ < 0\n");
		MPI_Finalize();
		exit(0);
	}
	if(previ >= xnumberAdded) {
		printf("particle previ >= xnumberAdded\n");
		MPI_Finalize();
		exit(0);
	}
	if(ynumberGeneral > 1) {
		if(j < 0) {
			printf("particle j < 0\n");
			MPI_Finalize();
			exit(0);
		}
		if(j >= ynumberAdded) {
			printf("particle j >= ynumberAdded\n");
			MPI_Finalize();
			exit(0);
		}
		if(prevj < 0) {
			printf("particle prevj < 0\n");
			MPI_Finalize();
			exit(0);
		}
		if(prevj >= ynumberAdded) {
			printf("particle prevj >= ynumberAdded\n");
			MPI_Finalize();
			exit(0);
		}
	}
	if(znumberGeneral > 1) {
		if(k < 0) {
			printf("particle k < 0\n");
			MPI_Finalize();
			exit(0);
		}
		if(k >= znumberAdded) {
			printf("particle k >= znumberAdded\n");
			MPI_Finalize();
			exit(0);
		}
		if(prevk < 0) {
			printf("particle prevk < 0\n");
			MPI_Finalize();
			exit(0);
		}
		if(prevk >= znumberAdded) {
			printf("particle prevk >= znumberAdded\n");
			MPI_Finalize();
			exit(0);
		}
	}*/

	double xr = (localCoordinates.x + prevLocalCoordinates.x)/2;
	double yr = (localCoordinates.y + prevLocalCoordinates.y)/2;
	double zr = (localCoordinates.z + prevLocalCoordinates.z)/2;

	switch(dimensionType){
	case THREE_D:{ 
		if(i != previ) {
			xr = max2(xgrid[previ], xgrid[i]) - xgrid[0];
		}
		if(j != prevj) {
			yr = max2(ygrid[prevj], ygrid[j]) - ygrid[0];
		}
		if(k != prevk) {
			zr = max2(zgrid[prevk], zgrid[k]) - zgrid[0];
		}

		//not optimize on vx because of reflection on the left wall
		double Fx1 = fullChargeDensity*(xr - prevLocalCoordinates.x)/deltaT;
		double Fy1 = fullChargeDensity*(yr - prevLocalCoordinates.y)/deltaT;
		double Fz1 = fullChargeDensity*(zr - prevLocalCoordinates.z)/deltaT;
		double Fx2 = fullChargeDensity*(localCoordinates.x - xr)/deltaT;
		double Fy2 = fullChargeDensity*velocity.y - Fy1;
		double Fz2 = fullChargeDensity*velocity.z - Fz1;

		double Wx1 = (prevLocalCoordinates.x + xr)/(2*deltaX) - previ;
		double Wy1 = (prevLocalCoordinates.y + yr)/(2*deltaY) - prevj;
		double Wz1 = (prevLocalCoordinates.z + zr)/(2*deltaZ) - prevk;
		double Wx2 = (localCoordinates.x + xr)/(2*deltaX) - i;
		double Wy2 = (localCoordinates.y + yr)/(2*deltaY) - j;
		double Wz2 = (localCoordinates.z + zr)/(2*deltaZ) - k;
		double onemWx1 = 1.0 - Wx1;
		double onemWy1 = 1.0 - Wy1;
		double onemWz1 = 1.0 - Wz1;
		double onemWx2 = 1.0 - Wx2;
		double onemWy2 = 1.0 - Wy2;
		double onemWz2 = 1.0 - Wz2;

		bunemanJx[previ][prevj][prevk] += Fx1*onemWy1*onemWz1;
		bunemanJx[previ][prevj+1][prevk] += Fx1*Wy1*onemWz1;
		bunemanJx[previ][prevj][prevk+1] += Fx1*onemWy1*Wz1;
		bunemanJx[previ][prevj+1][prevk+1] += Fx1*Wy1*Wz1;
		bunemanJy[previ][prevj][prevk] += Fy1*onemWx1*onemWz1;
		bunemanJy[previ+1][prevj][prevk] += Fy1*Wx1*onemWz1;
		bunemanJy[previ][prevj][prevk+1] += Fy1*onemWx1*Wz1;
		bunemanJy[previ+1][prevj][prevk+1] += Fy1*Wx1*Wz1;
		bunemanJz[previ][prevj][prevk] += Fz1*onemWx1*onemWy1;
		bunemanJz[previ][prevj+1][prevk] += Fz1*onemWx1*Wy1;
		bunemanJz[previ+1][prevj][prevk] += Fz1*Wx1*onemWy1;
		bunemanJz[previ+1][prevj+1][prevk] += Fz1*Wx1*Wy1;

		bunemanJx[i][j][k] += Fx2*onemWy2*onemWz2;
		bunemanJx[i][j+1][k] += Fx2*Wy2*onemWz2;
		bunemanJx[i][j][k+1] += Fx2*onemWy2*Wz2;
		bunemanJx[i][j+1][k+1] += Fx2*Wy2*Wz2;
		bunemanJy[i][j][k] += Fy2*onemWx2*onemWz2;
		bunemanJy[i+1][j][k] += Fy2*Wx2*onemWz2;
		bunemanJy[i][j][k+1] += Fy2*onemWx2*Wz2;
		bunemanJy[i+1][j][k+1] += Fy2*Wx2*Wz2;
		bunemanJz[i][j][k] += Fz2*onemWx2*onemWy2;
		bunemanJz[i][j+1][k] += Fz2*onemWx2*Wy2;
		bunemanJz[i+1][j][k] += Fz2*Wx2*onemWy2;
		bunemanJz[i+1][j+1][k] += Fz2*Wx2*Wy2;
				 }
		break;
	case TWO_D_XY:{
	
		if(i != previ) {
			xr = max2(xgrid[previ], xgrid[i]) - xgrid[0];
		}
		if(j != prevj) {
			yr = max2(ygrid[prevj], ygrid[j]) - ygrid[0];
		}

		double Fx1 = fullChargeDensity*(xr - prevLocalCoordinates.x)/deltaT;
		double Fy1 = fullChargeDensity*(yr - prevLocalCoordinates.y)/deltaT;
		double Fz1 = fullChargeDensity*(zr - prevLocalCoordinates.z)/deltaT;
		double Fx2 = fullChargeDensity*(localCoordinates.x - xr)/deltaT;
		double Fy2 = fullChargeDensity*velocity.y - Fy1;
		double Fz2 = fullChargeDensity*velocity.z - Fz1;

		double Wx1 = (prevLocalCoordinates.x + xr)/(2*deltaX) - previ;
		double Wy1 = (prevLocalCoordinates.y + yr)/(2*deltaY) - prevj;
		double Wx2 = (localCoordinates.x + xr)/(2*deltaX) - i;
		double Wy2 = (localCoordinates.y + yr)/(2*deltaY) - j;

		double onemWx1 = 1.0 - Wx1;
		double onemWy1 = 1.0 - Wy1;
		double onemWx2 = 1.0 - Wx2;
		double onemWy2 = 1.0 - Wy2;

		bunemanJx[previ][prevj][0] += Fx1*onemWy1;
		bunemanJx[previ][prevj+1][0] += Fx1*Wy1;
		bunemanJy[previ][prevj][0] += Fy1*onemWx1;
		bunemanJy[previ+1][prevj][0] += Fy1*Wx1;

		bunemanJz[previ][prevj][0] += Fz1*onemWx1*onemWy1;
		bunemanJz[previ][prevj+1][0] += Fz1*onemWx1*Wy1;
		bunemanJz[previ+1][prevj][0] += Fz1*Wx1*onemWy1;
		bunemanJz[previ+1][prevj+1][0] += Fz1*Wx1*Wy1;

		bunemanJx[i][j][0] += Fx2*onemWy2;
		bunemanJx[i][j+1][0] += Fx2*Wy2;
		bunemanJy[i][j][0] += Fy2*onemWx2;
		bunemanJy[i+1][j][0] += Fy2*Wx2;

		bunemanJz[i][j][0] += Fz2*onemWx2*onemWy2;
		bunemanJz[i][j+1][0] += Fz2*onemWx2*Wy2;
		bunemanJz[i+1][j][0] += Fz2*Wx2*onemWy2;
		bunemanJz[i+1][j+1][0] += Fz2*Wx2*Wy2;
				  }
		break;
		
	case TWO_D_XZ:{
		if(i != previ) {
			xr = max2(xgrid[previ], xgrid[i]) - xgrid[0];
		}
		if(k != prevk) {
			zr = max2(zgrid[prevk], zgrid[k]) - zgrid[0];
		}

		double Fx1 = fullChargeDensity*(xr - prevLocalCoordinates.x)/deltaT;
		double Fy1 = fullChargeDensity*(yr - prevLocalCoordinates.y)/deltaT;
		double Fz1 = fullChargeDensity*(zr - prevLocalCoordinates.z)/deltaT;
		double Fx2 = fullChargeDensity*(localCoordinates.x - xr)/deltaT;
		double Fy2 = fullChargeDensity*velocity.y - Fy1;
		double Fz2 = fullChargeDensity*velocity.z - Fz1;

		double Wx1 = (prevLocalCoordinates.x + xr)/(2*deltaX) - previ;
		double Wz1 = (prevLocalCoordinates.z + zr)/(2*deltaZ) - prevk;
		double Wx2 = (localCoordinates.x + xr)/(2*deltaX) - i;
		double Wz2 = (localCoordinates.z + zr)/(2*deltaZ) - k;
		double onemWx1 = 1.0 - Wx1;
		double onemWz1 = 1.0 - Wz1;
		double onemWx2 = 1.0 - Wx2;
		double onemWz2 = 1.0 - Wz2;

		bunemanJx[previ][0][prevk] += Fx1*onemWz1;
		bunemanJx[previ][0][prevk+1] += Fx1*Wz1;
		bunemanJz[previ][0][prevk] += Fz1*onemWx1;
		bunemanJz[previ+1][0][prevk] += Fz1*Wx1;

		bunemanJy[previ][0][prevk] += Fy1*onemWx1*onemWz1;
		bunemanJy[previ][0][prevk+1] += Fy1*onemWx1*Wz1;
		bunemanJy[previ+1][0][prevk] += Fy1*Wx1*onemWz1;
		bunemanJy[previ+1][0][prevk+1] += Fy1*Wx1*Wz1;

		bunemanJx[i][0][k] += Fx2*onemWz2;
		bunemanJx[i][0][k+1] += Fx2*Wz2;
		bunemanJz[i][0][k] += Fz2*onemWx2;
		bunemanJz[i+1][0][k] += Fz2*Wx2;

		bunemanJy[i][0][k] += Fy2*onemWx2*onemWz2;
		bunemanJy[i][0][k+1] += Fy2*onemWx2*Wz2;
		bunemanJy[i+1][0][k] += Fy2*Wx2*onemWz2;
		bunemanJy[i+1][0][k+1] += Fy2*Wx2*Wz2;
				  }
		break;
	case ONE_D:{
		if(i != previ) {
			xr = max2(xgrid[previ], xgrid[i]) - xgrid[0];
		}

		double Fx1 = fullChargeDensity*(xr - prevLocalCoordinates.x)/deltaT;
		double Fy1 = fullChargeDensity*(yr - prevLocalCoordinates.y)/deltaT;
		double Fz1 = fullChargeDensity*(zr - prevLocalCoordinates.z)/deltaT;
		double Fx2 = fullChargeDensity*(localCoordinates.x - xr)/deltaT;
		double Fy2 = fullChargeDensity*velocity.y - Fy1;
		double Fz2 = fullChargeDensity*velocity.z - Fz1;

		double Wx1 = (prevLocalCoordinates.x + xr)/(2*deltaX) - previ;
		double Wx2 = (localCoordinates.x + xr)/(2*deltaX) - i;

		bunemanJx[previ][0][0] += Fx1;

		bunemanJy[previ][0][0] += Fy1*(1.0 - Wx1);
		bunemanJy[previ+1][0][0] += Fy1*Wx1;

		bunemanJz[previ][0][0] += Fz1*(1 - Wx1);
		bunemanJz[previ+1][0][0] += Fz1*Wx1;

		bunemanJx[i][0][0] += Fx2;

		bunemanJy[i][0][0] += Fy2*(1.0 - Wx2);
		bunemanJy[i+1][0][0] += Fy2*Wx2;

		bunemanJz[i][0][0] += Fz2*(1 - Wx2);
		bunemanJz[i+1][0][0] += Fz2*Wx2;
			   }
		break;
	default:
		printf("wrong dimension type\n");
		MPI_Finalize();
		exit(0);
	}
}

void Simulation::updateElectroMagneticParameters() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
		//MPI_Barrier(cartComm);
		if ((rank == 0) && (verbosity > 0)) printf("updating flux, density snd dielectric tensor\n");
		if ((rank == 0) && (verbosity > 0)) printLog("updating flux, density and dielectric tensor\n");

		int particlePartsCount = 0;
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					electricFlux[i][j][k] = Vector3d(0, 0, 0);
					electricFluxMinus[i][j][k] = Vector3d(0, 0, 0);
					dielectricTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
					divPressureTensor[i][j][k] = Vector3d(0, 0, 0);
					divPressureTensorMinus[i][j][k] = Vector3d(0, 0, 0);
				}
			}
		}

		int crossBinNumberX = splineOrder + 2;
		int crossBinNumberY = splineOrder + 2;
		int crossBinNumberZ = splineOrder + 2;
		if(ynumberGeneral == 1) {
			crossBinNumberY = 1;
		}
		if(znumberGeneral == 1) {
			crossBinNumberZ = 1;
		}

		for (int pcount = 0; pcount < particles.size(); ++pcount) {
			Particle* particle = particles[pcount];

			Vector3d velocity = particle->getVelocity();
			double gamma = particle->gammaFactor();

			Vector3d rotatedVelocity = particle->rotationTensor * (velocity * gamma);
			
			double particleOmega = particle->weight * theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge /
				particle->mass;


			for (int i = 0; i < crossBinNumberX; ++i) {
				for (int j = 0; j < crossBinNumberY; ++j) {
					for (int k = 0; k < crossBinNumberZ; ++k) {
						int curI = particle->correlationMapNode.xindex[i];
						int curJ = particle->correlationMapNode.yindex[j];
						if(ynumberGeneral == 1) {
							curJ = 0;
						}
						int curK = particle->correlationMapNode.zindex[k];
						if(znumberGeneral == 1) {
							curK = 0;
						}

						double correlation;
						switch(Simulation::dimensionType){
							case DimensionType::THREE_D: correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.ycorrelation[j] * particle->correlationMapNode.zcorrelation[k] / volumeE();
								break;
							case DimensionType::TWO_D_XY: correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.ycorrelation[j]/ volumeE();
								break;
							case DimensionType::TWO_D_XZ: correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.zcorrelation[k] / volumeE();
								break;
							case DimensionType::ONE_D: correlation =  particle->correlationMapNode.xcorrelation[i]/ volumeE();
								break;
							default: correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.ycorrelation[j] * particle->correlationMapNode.zcorrelation[k] / volumeE();
								break;				
						}

						double particleCharge = particle->charge * particle->weight;
						if (curI <= additionalBinNumber && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
							Particle tempParticle = *particle;
							tempParticle.reflectMomentumX();
							Vector3d tempVelocity = velocity;
							tempVelocity.x = -tempVelocity.x;
							//tempParticle.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);
							Vector3d tempRotatedVelocity = tempParticle.rotationTensor * (tempVelocity * gamma);
							if (solverType == IMPLICIT) {
								if (particleCharge > 0) {
									electricFlux[2 + 2 * additionalBinNumber - curI][curJ][curK] += tempRotatedVelocity * (particleCharge * correlation);																
								} else {
									electricFluxMinus[2 + 2 * additionalBinNumber - curI][curJ][curK] += tempRotatedVelocity * (particleCharge *correlation);														
								}								
							} else {
								if (particleCharge > 0) {
									electricFlux[2 + 2 * additionalBinNumber - curI][curJ][curK] += tempVelocity * particleCharge * correlation;
								} else {
									electricFluxMinus[2 + 2 * additionalBinNumber - curI][curJ][curK] += tempVelocity * particleCharge * correlation;
								}
							}
						} else {

							if (solverType == IMPLICIT) {
								if (particleCharge > 0) {
									electricFlux[curI][curJ][curK] += rotatedVelocity * (particleCharge * correlation);
								} else {
									electricFluxMinus[curI][curJ][curK] += rotatedVelocity * (particleCharge * correlation);
								}
								dielectricTensor[curI][curJ][curK] = dielectricTensor[curI][curJ][curK] - particle->rotationTensor * (
									particleOmega * correlation);


							} else if (solverType == EXPLICIT) {
								if (particleCharge > 0) {
									electricFlux[curI][curJ][curK] += velocity * particleCharge * correlation;
								} else {
									electricFluxMinus[curI][curJ][curK] += velocity * particleCharge * correlation;
								}
							}
						}
					}
				}
			}
		}

		//for debug onle
		/*for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					electricFlux[i][j][k] = Vector3d(0, 0, 0);
					electricFluxMinus[i][j][k] = Vector3d(0, 0, 0);
					dielectricTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
					divPressureTensor[i][j][k] = Vector3d(0, 0, 0);
				}
			}
		}*/

		///


		if (cartCoord[0] == 0 && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					electricFlux[1 + additionalBinNumber][j][k] = Vector3d(0, 0, 0);
					electricFluxMinus[1 + additionalBinNumber][j][k] = Vector3d(0, 0, 0);
					divPressureTensor[1 + additionalBinNumber][j][k] = Vector3d(0, 0, 0);
					dielectricTensor[1 + additionalBinNumber][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				}
			}
		}

		if ((verbosity > 1)) printf("updating electricDensityHat and pressure tensor\n");

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					chargeDensityHat[i][j][k] = 0;
					chargeDensity[i][j][k] = 0;
					chargeDensityMinus[i][j][k] = 0;
					pressureTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				}
			}
		}

		for (int pcount = 0; pcount < particles.size(); ++pcount) {
			Particle* particle = particles[pcount];
			Vector3d velocity = particle->getVelocity();
			double gamma = particle->gammaFactor();
			Vector3d rotatedVelocity = particle->rotationTensor * (velocity * gamma);
			Matrix3d tensor = rotatedVelocity.selfTensorMult();

			for (int i = 0; i < crossBinNumberX; ++i) {
				for (int j = 0; j < crossBinNumberY; ++j) {
					for (int k = 0; k < crossBinNumberZ; ++k) {
						int curI = particle->correlationMapCell.xindex[i];
						int curJ = particle->correlationMapCell.yindex[j];
						if(ynumberGeneral == 1) {
							curJ = 0;
						}
						int curK = particle->correlationMapCell.zindex[k];
						if(znumberGeneral == 1) {
							curK = 0;
						}

						double correlation;
						switch(Simulation::dimensionType){
							case DimensionType::THREE_D: correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.ycorrelation[j] * particle->correlationMapCell.zcorrelation[k] / volumeB();
								break;
							case DimensionType::TWO_D_XY: correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.ycorrelation[j]/ volumeB();
								break;
							case DimensionType::TWO_D_XZ: correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.zcorrelation[k] / volumeB();
								break;
							case DimensionType::ONE_D: correlation =  particle->correlationMapCell.xcorrelation[i]/ volumeB();
								break;
							default: correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.ycorrelation[j] * particle->correlationMapCell.zcorrelation[k] / volumeB();
								break;				
						}
						double particleCharge = particle->charge * particle->weight;


						if (curI <= additionalBinNumber && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
							if (particleCharge > 0) {
								chargeDensityHat[1 + 2 * additionalBinNumber - curI][curJ][curK] += particleCharge * correlation;
							} else {
								chargeDensityMinus[1 + 2 * additionalBinNumber - curI][curJ][curK] += particleCharge * correlation;
							}
							for(int l = 0; l < 3; ++l) {
								for(int m = 0; m < 3; ++m) {
									pressureTensor[1 + 2 * additionalBinNumber - curI][j][k].matrix[l][m] += particleCharge*correlation*tensor.matrix[l][m];
								}
							}
						} else {

							if (particleCharge > 0) {
								chargeDensityHat[curI][curJ][curK] += particleCharge * correlation;
							} else {
								chargeDensityMinus[curI][curJ][curK] += particleCharge * correlation;
							}
							for(int l = 0; l < 3; ++l) {
								for(int m = 0; m < 3; ++m) {
									pressureTensor[curI][j][k].matrix[l][m] += particleCharge*correlation*tensor.matrix[l][m];
								}
							}
						}
					}
				}
			}
		}

		if (rank == 0 && (verbosity > 0)) printf("start exchange parameters\n");
		if ((verbosity > 1)) printf("sum charge density hat\n");
		////todo geometry!!!!
		//MPI_Barrier(cartComm);
		sumChargeDensityHatX();
		sumChargeDensityHatY();
		sumChargeDensityHatZ();
		if ((verbosity > 1)) printf("sumcell matrix parameters\n");
		//MPI_Barrier(cartComm);
		sumCellMatrixParametersX();
		sumCellMatrixParametersY();
		sumCellMatrixParametersZ();

		if ((verbosity > 1)) printf("sum node vector parameters\n");
		//MPI_Barrier(cartComm);
		sumNodeVectorParametersX();
		sumNodeVectorParametersY();
		sumNodeVectorParametersZ();

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					electricFlux[i][j][k] += electricFluxMinus[i][j][k];
					divPressureTensor[i][j][k] += divPressureTensorMinus[i][j][k];
				}
			}
		}


		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					chargeDensityHat[i][j][k] += chargeDensityMinus[i][j][k];
				}
			}
		}

		if (solverType == IMPLICIT) {
			for (int i = 0; i < xnumberAdded + 1; ++i) {
				for (int j = 0; j < ynumberAdded + 1; ++j) {
					for (int k = 0; k < znumberAdded + 1; ++k) {
						//todo
						divPressureTensor[i][j][k] = evaluateDivPressureTensor(i, j, k);
						electricFlux[i][j][k] = electricFlux[i][j][k] - divPressureTensor[i][j][k] * eta * deltaT;
						//alertNaNOrInfinity(electricFlux[i][j][k].x, "electricFlux[i][j][k].x = NaN");
						//alertNaNOrInfinity(electricFlux[i][j][k].y, "electricFlux[i][j][k].y = NaN");
						//alertNaNOrInfinity(electricFlux[i][j][k].z, "electricFlux[i][j][k].z = NaN");
					}
				}
			}
		}

		//todo realy?
		if (boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					for (int i = 0; i <= additionalBinNumber + 1; ++i) {
						electricFlux[i][j][k] = Vector3d(0, 0, 0);
					}
				}
			}
		}

		//here for evaluating div J

		if (solverType == IMPLICIT) {
			int minI = 0;
			/*if(boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0){
				minI = 2 + additionalBinNumber;
			}*/
			for (int i = minI; i < xnumberAdded; ++i) {
				for (int j = 0; j < ynumberAdded; ++j) {
					for (int k = 0; k < znumberAdded; ++k) {
						double divJ = evaluateDivFlux(i, j, k);

						chargeDensityHat[i][j][k] -= deltaT * theta * divJ;
						//chargeDensityHat[i][j][k] = 0;
					}
				}
			}
		}

		//zero densities!!!
		/*for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					chargeDensityHat[i][j][k] = 0;
					chargeDensity[i][j][k] = 0;
					chargeDensityMinus[i][j][k] = 0;
				}
			}
		}*/
		//////

		///zero flux
		/*for(int i = 0; i < xnumberAdded + 1; ++i){
			for(int j = 0; j < ynumberAdded + 1; ++j){
				for(int k = 0; k < znumberAdded + 1; ++k){
					electricFlux[i][j][k] = Vector3d(0, 0, 0);
				}
			}
		}*/
		////

		///zero tensor
		/*for(int i = 0; i < xnumberAdded + 1; ++i){
			for(int j = 0; j < ynumberAdded + 1; ++j){
				for(int k = 0; k < znumberAdded + 1; ++k){
					dielectricTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				}
			}
		}*/
		////
		//fprintf(outputEverythingFile, "density %d afterFlux = %28.22g\n", chargeDensityHat[debugPoint - 1][0][0]);
		//fprintf(outputEverythingFile, "density %d afterFlux = %28.22g\n", chargeDensityHat[debugPoint][0][0]);
		//fprintf(outputEverythingFile, "density %d afterFlux = %28.22g\n", chargeDensityHat[debugPoint + 1][0][0]);

		/*FILE* tempDensityFile = fopen((outputDir + "tempDensityFile.dat").c_str(), "w");
		for(int i = 0; i < xnumber; ++i){
		    fprintf(tempDensityFile, "%28.22g\n", chargeDensityHat[i][0][0]);
		}
		fclose(tempDensityFile);*/
		//MPI_Barrier(cartComm);
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock() - procTime;
			printf("updating electromagnetic parameters time = %g sec\n", procTime / CLOCKS_PER_SEC);
		}
		if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
			procTime = clock();
		}
		//MPI_Barrier(cartComm);

		/*for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].x, "electricFlux[i][j][k].x = NaN");
					if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].y, "electricFlux[i][j][k].y = NaN");
					if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].z, "electricFlux[i][j][k].z = NaN");
					alertLargeModlue(electricFlux[i][j][k].x, 1E100, "electric Flux x is too large\n");
					alertLargeModlue(electricFlux[i][j][k].y, 1E100, "electric Flux y is too large\n");
					alertLargeModlue(electricFlux[i][j][k].z, 1E100, "electric Flux z is too large\n");
				}
			}
		}*/

		if ((verbosity > 1)) printf("sum node matrix parameters\n");
		//MPI_Barrier(cartComm);
		sumNodeMatrixParametersX();
		sumNodeMatrixParametersY();
		sumNodeMatrixParametersZ();
		if ((verbosity > 1)) printf("update external flux\n");
		//updateExternalFlux();

		/*for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded + 1; ++j) {
				for (int k = 0; k < znumberAdded + 1; ++k) {
					electricFlux[i][j][k] = electricFlux[i][j][k] + externalElectricFlux[i][j][k];
					//electricFlux[i][j][k].x = 0;
					//electricFlux[i][j][k].y = 0;
					//electricFlux[i][j][k].z = 0;
					if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].x, "electricFlux[i][j][k].x = NaN");
					if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].y, "electricFlux[i][j][k].y = NaN");
					if (debugMode) alertNaNOrInfinity(electricFlux[i][j][k].z, "electricFlux[i][j][k].z = NaN");
				}
			}
		}*/
	if ((verbosity > 1)) printf("finish updating electromagnetic parameters\n");
	//MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("exchanging electromagnetic parameters time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::updateDensityParameters() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if ((verbosity > 0)) {
		printf("updating densityParameters\n");
	}
	if (solverType == BUNEMAN) {
		updateBunemanChargeDensity();
	}
	//FILE* debugFile = fopen("./output/particleCorrelations.dat","w");
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int t = 0; t < typesNumber; ++t) {
					particleConcentrations[t][i][j][k] = 0;
					particleEnergies[t][i][j][k] = 0;
					particleBulkVelocities[t][i][j][k] = Vector3d(0, 0, 0);
				}
				chargeDensity[i][j][k] = 0;
				chargeDensityMinus[i][j][k] = 0;
			}
		}
	}

	int crossBinNumberX = splineOrder + 2;
	int crossBinNumberY = splineOrder + 2;
	int crossBinNumberZ = splineOrder + 2;
	if(ynumberGeneral == 1) {
		crossBinNumberY = 1;
	}
	if(znumberGeneral == 1) {
		crossBinNumberZ = 1;
	}

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		for (int i = 0; i < crossBinNumberX; ++i) {
			//todo
			for (int j = 0; j < crossBinNumberY; ++j) {
				for (int k = 0; k < crossBinNumberZ; ++k) {
					int curI = particle->correlationMapCell.xindex[i];
					int curJ = particle->correlationMapCell.yindex[j];
					if(ynumberGeneral == 1) {
						curJ = 0;
					}
					int curK = particle->correlationMapCell.zindex[k];
					if(znumberGeneral == 1) {
						curK = 0;
					}
					if (curI < 0) {
						printf("curI < 0 %d particle number = %d\n", curI, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					if (curI > xnumberAdded - 1) {
						printf("curI > xnumberAdded - 1 %d particle number = %d\n", curI, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					if(ynumberGeneral > 1){
					if (curJ < 0) {
						printf("curJ < 0 %d particle number = %d\n", curJ, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					if (curJ > ynumberAdded - 1) {
						printf("curJ > ynumberAdded - 1 %d particle number = %d\n", curJ, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					}
					if(znumberGeneral > 1){
					if (curK < 0) {
						printf("curK < 0 %d particle number = %d\n", curK, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					if (curK > znumberAdded - 1) {
						printf("curK > znumberAdded - 1 %d particle number = %d\n", curK, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					}
					double correlation;
					switch(Simulation::dimensionType){
						case DimensionType::THREE_D: correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.ycorrelation[j] * particle->correlationMapCell.zcorrelation[k] / volumeB();
							break;
						case DimensionType::TWO_D_XY: correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.ycorrelation[j]/ volumeB();
							break;
						case DimensionType::TWO_D_XZ: correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.zcorrelation[k] / volumeB();
							break;
						case DimensionType::ONE_D: correlation =  particle->correlationMapCell.xcorrelation[i]/ volumeB();
							break;
						default: correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.ycorrelation[j] * particle->correlationMapCell.zcorrelation[k] / volumeB();
							break;				
					}
					double particleCharge = particle->charge * particle->weight;
					int typeN = getTypeNumber(particle);

					if (curI <= additionalBinNumber && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
						Vector3d reflectedMomentum = particle->getMomentum();
						reflectedMomentum.x = -reflectedMomentum.x;
						particleBulkVelocities[typeN][1 + 2 * additionalBinNumber - curI][curJ][curK] += reflectedMomentum * particle->
							weight * correlation;
						if (particleCharge > 0) {
							chargeDensity[1 + 2 * additionalBinNumber - curI][curJ][curK] += particleCharge * correlation;
						} else {
							chargeDensityMinus[1 + 2 * additionalBinNumber - curI][curJ][curK] += particleCharge * correlation;
						}
						particleConcentrations[typeN][1 + 2 * additionalBinNumber - curI][curJ][curK] += particle->weight * correlation;
						particleEnergies[typeN][1 + 2 * additionalBinNumber - curI][curJ][curK] += particle->weight * particle->
							fullEnergy() * correlation;
					} else {
						particleBulkVelocities[typeN][curI][curJ][curK] += particle->getMomentum() * particle->weight * correlation;
						if (particleCharge > 0) {
							chargeDensity[curI][curJ][curK] += particleCharge * correlation;
						} else {
							chargeDensityMinus[curI][curJ][curK] += particleCharge * correlation;
						}
						particleConcentrations[typeN][curI][curJ][curK] += particle->weight * correlation;
						particleEnergies[typeN][curI][curJ][curK] += particle->weight * particle->fullEnergy() *
							correlation;
					}
				}
			}
		}
	}

	//MPI_Barrier(cartComm);
	if ((verbosity > 0)) {
		printf("sum densityParameters\n");
	}

	if ((verbosity > 2)) {
		printf("sum cell parameters x\n");
	}
	sumCellParametersX();
	//MPI_Barrier(cartComm);
	if ((rank == 0) && (verbosity > 0)) {
		printf("finish sum cell parameters x\n");
	}
	if ((verbosity > 2)) {
		printf("sum cell parameters y\n");
	}
	sumCellParametersY();
	//MPI_Barrier(cartComm);
	if ((rank == 0) && (verbosity > 0)) {
		printf("finish sum cell parameters y\n");
	}
	if ((verbosity > 2)) {
		printf("sum cell parameters z\n");
	}
	sumCellParametersZ();
	//MPI_Barrier(cartComm);
	if ((rank == 0) && (verbosity > 0)) {
		printf("finish sum cell parameters z\n");
	}
	if ((verbosity > 2)) {
		printf("sum cell  xector parameters x\n");
	}
	sumCellVectorParametersX();
	//MPI_Barrier(cartComm);
	if ((rank == 0) && (verbosity > 1)) {
		printf("finish sum cell vector parameters x\n");
	}
	if ((verbosity > 2)) {
		printf("sum cell  xector parameters y\n");
	}
	sumCellVectorParametersY();
	if ((rank == 0) && (verbosity > 1)) {
		printf("finish sum cell vector parameters y\n");
	}
	if ((verbosity > 2)) {
		printf("sum cell  xector parameters z\n");
	}
	sumCellVectorParametersZ();
	if ((rank == 0) && (verbosity > 1)) {
		printf("finish sum cell vector parameters z\n");
	}


	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				chargeDensity[i][j][k] += chargeDensityMinus[i][j][k];
			}
		}
	}

	//zero density!!!
	/*for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				chargeDensity[i][j][k] = 0;
				chargeDensityMinus[i][j][k] = 0;
			}
		}
	}
	*/
	//////

	if ((verbosity > 0)) {
		printf("evaluating velocity\n");
	}
	//todo v = c^2 * p/E
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				//fprintf(debugFile, "charge %15.10g proton %15.10g electron %15.10g\n", chargeDensity[i][j][k], protonConcentration[i][j][k], electronConcentration[i][j][k]);
				for (int t = 0; t < typesNumber; ++t) {
					if (particleEnergies[t][i][j][k] > 0) {
						/*particleBulkVelocities[t][i][j][k] = particleBulkVelocities[t][i][j][k] / (particleConcentrations[t][i][j][k] * types[t].mass);
						double gamma = sqrt((particleBulkVelocities[t][i][j][k].scalarMult(
							particleBulkVelocities[t][i][j][k]) / speed_of_light_normalized_sqr) + 1);*/
						//particleBulkVelocities[t][i][j][k] = particleBulkVelocities[t][i][j][k] * (speed_of_light_normalized_sqr /particleEnergies[t][i][j][k]);
						particleBulkVelocities[t][i][j][k] = particleBulkVelocities[t][i][j][k] / particleEnergies[t][i][j][k];
					} else {
						particleBulkVelocities[t][i][j][k] = particleBulkVelocities[t][i][j][k] = Vector3d(0, 0, 0);
					}
				}
			}
		}
	}
	if ((rank == 0) && (verbosity > 0)) {
		printf("finish evaluating velocity\n");
	}
	if ((verbosity > 2)) {
		printf("finish update density parameters\n");
	}
	//MPI_Barrier(cartComm);
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("updating density parameters time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::updateBunemanChargeDensity() {
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanChargeDensity[i][j][k] = 0;
			}
		}
	}

	int crossBinNumberX = splineOrder + 2;
	int crossBinNumberY = splineOrder + 2;
	int crossBinNumberZ = splineOrder + 2;
	if(ynumberGeneral == 1) {
		crossBinNumberY = 1;
	}
	if(znumberGeneral == 1) {
		crossBinNumberZ = 1;
	}

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		for (int i = 0; i < crossBinNumberX; ++i) {
			//todo
			for (int j = 0; j < crossBinNumberY; ++j) {
				for (int k = 0; k < crossBinNumberZ; ++k) {
					int curI = particle->correlationMapNode.xindex[i];
					int curJ = particle->correlationMapNode.yindex[j];
					if(ynumberGeneral == 1) {
						curJ = 0;
					}
					int curK = particle->correlationMapNode.zindex[k];
					if(znumberGeneral == 1) {
						curK = 0;
					}
					if (curI < 0) {
						printf("curI < 0 %d particle number = %d\n", curI, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					if (curI > xnumberAdded) {
						printf("curI > xnumberAdded - 1 %d particle number = %d\n", curI, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					if (curJ < 0) {
						printf("curJ < 0 %d particle number = %d\n", curJ, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					if (curJ > ynumberAdded) {
						printf("curJ > ynumberAdded - 1 %d particle number = %d\n", curJ, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					if (curK < 0) {
						printf("curK < 0 %d particle number = %d\n", curK, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					if (curK > znumberAdded) {
						printf("curK > znumberAdded - 1 %d particle number = %d\n", curK, particle->number);
						printf("particle.x = %g particle.y = %g particle.z = %g, z[0] = %g z[znumberAdded] = %g", particle->coordinates.x,
						       particle->coordinates.y, particle->coordinates.z, zgrid[0], zgrid[znumberAdded]);
						Vector3d velocity = particle->getVelocity();
						printf("particle.vx = %g particle.vy = %g particle.vz = %g", velocity.x, velocity.y, velocity.z);
						MPI_Finalize();
						exit(0);
					}
					double correlation;
					switch(Simulation::dimensionType){
						case DimensionType::THREE_D: correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.ycorrelation[j] * particle->correlationMapNode.zcorrelation[k] / volumeE();
							break;
						case DimensionType::TWO_D_XY: correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.ycorrelation[j]/ volumeE();
							break;
						case DimensionType::TWO_D_XZ: correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.zcorrelation[k] / volumeE();
							break;
						case DimensionType::ONE_D: correlation =  particle->correlationMapNode.xcorrelation[i]/ volumeE();
							break;
						default: correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.ycorrelation[j] * particle->correlationMapNode.zcorrelation[k] / volumeE();
							break;				
					}
					double particleCharge = particle->charge * particle->weight;

					if (curI <= additionalBinNumber && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
						bunemanChargeDensity[2 + 2 * additionalBinNumber - curI][curJ][curK] += particleCharge * correlation;
					} else {
						bunemanChargeDensity[curI][curJ][curK] += particleCharge * correlation;
					}
				}
			}
		}
	}
	sumNodeParametersX();
	sumNodeParametersY();
	sumNodeParametersZ();
}

void Simulation::updateExternalFlux() {
	double alfvenV;
	if (density <= 1E-100) {
		//alfvenV = speed_of_light_normalized;
		alfvenV =1.0;
	} else {
		alfvenV = B0.norm() / sqrt(4 * pi * density);
	}
	double concentration = density / (massProton + massElectron);
	double phaseV = 2 * alfvenV;
	double kw = 2 * pi / xsize;
	double omega = kw * phaseV;

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				externalElectricFlux[i][j][k] = Vector3d(0, 0, 1.0) * extJ * cos(kw * xgrid[i] - omega * time);
				//alertNaNOrInfinity(externalElectricFlux[i][j][k].x, "externalFlux.x = NaN\n");
			}
		}
	}
}

Vector3d Simulation::getBunemanFlux(int i, int j, int k) {
	return Vector3d(bunemanJx[i][j][k], bunemanJy[i][j][k], bunemanJz[i][j][k]);
}
