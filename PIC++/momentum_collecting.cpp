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
	if ((rank == 0) && (verbosity > 0)) printf("updating tristan flux\n");
	if ((rank == 0) && (verbosity > 0)) printLog("updating tristan flux\n");
	int particlePartsCount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanJx[i][j][k] = 0;
			}
		}
	}
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanJy[i][j][k] = 0;
			}
		}
	}
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				bunemanJz[i][j][k] = 0;
			}
		}
	}

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		addParticleFluxZigzag(particle);
	}


	exchangeBunemanFlux();

	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock() - procTime;
		printf("evaluating electric flux time = %g sec\n", procTime / CLOCKS_PER_SEC);
	}
}

void Simulation::addParticleFluxZigzag(Particle* particle) {
	Vector3d velocity = particle->getVelocity();
	Vector3d localCoordinates = particle->coordinates - Vector3d(xgrid[0], ygrid[0], zgrid[0]);
	Vector3d prevLocalCoordinates = localCoordinates - velocity*deltaT;
	double fullCharge = particle->charge*particle->weight;
	const int i = floor(localCoordinates.x/deltaX);
	const int j = floor(localCoordinates.y/deltaY);
	const int k = floor(localCoordinates.z/deltaZ);
	const int previ = floor(prevLocalCoordinates.x/deltaX);
	const int prevj = floor(prevLocalCoordinates.y/deltaY);
	const int prevk = floor(prevLocalCoordinates.z/deltaZ);

	double xr = (localCoordinates.x + prevLocalCoordinates.x)/2;
	double yr = (localCoordinates.y + prevLocalCoordinates.y)/2;
	double zr = (localCoordinates.z + prevLocalCoordinates.z)/2;

	if(dimensionType == THREE_D) {
		if(i != previ) {
			xr = max2(xgrid[previ], xgrid[i]);
		}
		if(j != prevj) {
			yr = max2(ygrid[prevj], ygrid[j]);
		}
		if(k != prevk) {
			zr = max2(zgrid[prevk], zgrid[k]);
		}

		double Fx1 = fullCharge*(xr - prevLocalCoordinates.x)/deltaT;
		double Fy1 = fullCharge*(yr - prevLocalCoordinates.y)/deltaT;
		double Fz1 = fullCharge*(zr - prevLocalCoordinates.z)/deltaT;
		double Fx2 = fullCharge*velocity.x - Fx1;
		double Fy2 = fullCharge*velocity.y - Fy1;
		double Fz2 = fullCharge*velocity.z - Fz1;

		double Wx1 = (prevLocalCoordinates.x + xr)/(2*deltaX) - previ;
		double Wy1 = (prevLocalCoordinates.y + yr)/(2*deltaY) - prevj;
		double Wz1 = (prevLocalCoordinates.z + zr)/(2*deltaZ) - prevk;
		double Wx2 = (localCoordinates.x + xr)/(2*deltaX) - i;
		double Wy2 = (localCoordinates.y + yr)/(2*deltaY) - j;
		double Wz2 = (localCoordinates.z + zr)/(2*deltaZ) - k;

		bunemanJx[previ][prevj][prevk] += Fx1*(1.0 - Wy1)*(1.0 - Wz1)/cellVolume;
		bunemanJx[previ][prevj+1][prevk] += Fx1*Wy1*(1 - Wz1)/cellVolume;
		bunemanJx[previ][prevj][prevk+1] += Fx1*(1.0 - Wy1)*Wz1/cellVolume;
		bunemanJx[previ][prevj+1][prevk+1] += Fx1*Wy1*Wz1/cellVolume;
		bunemanJy[previ][prevj][prevk] += Fy1*(1.0 - Wx1)*(1.0 - Wz1)/cellVolume;
		bunemanJy[previ+1][prevj][prevk] += Fy1*Wx1*(1.0 - Wz1)/cellVolume;
		bunemanJy[previ][prevj][prevk+1] += Fy1*(1.0 - Wx1)*Wz1/cellVolume;
		bunemanJy[previ+1][prevj][prevk+1] += Fy1*Wx1*Wz1/cellVolume;
		bunemanJz[previ][prevj][prevk] = Fz1*(1.0 - Wx1)*(1.0 - Wy1)/cellVolume;
		bunemanJz[previ][prevj+1][prevk] = Fz1*(1.0 - Wx1)*Wy1/cellVolume;
		bunemanJz[previ+1][prevj][prevk] = Fz1*Wx1*(1.0 - Wy1)/cellVolume;
		bunemanJz[previ+1][prevj+1][prevk] = Fz1*Wx1*Wy1/cellVolume;

		bunemanJx[i][j][k] += Fx2*(1.0 - Wy2)*(1.0 - Wz2)/cellVolume;
		bunemanJx[i][j+1][k] += Fx2*Wy2*(1.0 - Wz2)/cellVolume;
		bunemanJx[i][j][k+1] += Fx2*(1.0 - Wy2)*Wz2/cellVolume;
		bunemanJx[i][j+1][k+1] += Fx2*Wy2*Wz2/cellVolume;
		bunemanJy[i][j][k] += Fy2*(1.0 - Wx2)*(1.0 - Wz2)/cellVolume;
		bunemanJy[i+1][j][k] += Fy2*Wx2*(1.0 - Wz2)/cellVolume;
		bunemanJy[i][j][k+1] += Fy2*(1.0 - Wx2)*Wz2/cellVolume;
		bunemanJy[i+1][j][k+1] += Fy2*Wx2*Wz2/cellVolume;
		bunemanJz[i][j][k] = Fz2*(1.0 - Wx2)*(1.0 - Wy2)/cellVolume;
		bunemanJz[i][j+1][k] = Fz2*(1.0 - Wx2)*Wy2/cellVolume;
		bunemanJz[i+1][j][k] = Fz2*Wx2*(1.0 - Wy2)/cellVolume;
		bunemanJz[i+1][j+1][k] = Fz2*Wx2*Wy2/cellVolume;
		
	} else if(dimensionType == TWO_D_XY) {
	
		if(i != previ) {
			xr = max2(xgrid[previ], xgrid[i]);
		}
		if(j != prevj) {
			yr = max2(ygrid[prevj], ygrid[j]);
		}

		double Fx1 = fullCharge*(xr - prevLocalCoordinates.x)/deltaT;
		double Fy1 = fullCharge*(yr - prevLocalCoordinates.y)/deltaT;
		double Fz1 = fullCharge*(zr - prevLocalCoordinates.z)/deltaT;
		double Fx2 = fullCharge*velocity.x - Fx1;
		double Fy2 = fullCharge*velocity.y - Fy1;
		double Fz2 = fullCharge*velocity.z - Fz1;

		double Wx1 = (prevLocalCoordinates.x + xr)/(2*deltaX) - previ;
		double Wy1 = (prevLocalCoordinates.y + yr)/(2*deltaY) - prevj;
		double Wx2 = (localCoordinates.x + xr)/(2*deltaX) - i;
		double Wy2 = (localCoordinates.y + yr)/(2*deltaY) - j;

		bunemanJx[previ][prevj][0] += Fx1*(1.0 - Wy1)/cellVolume;
		bunemanJx[previ][prevj+1][0] += Fx1*Wy1/cellVolume;
		bunemanJy[previ][prevj][0] += Fy1*(1.0 - Wx1)/cellVolume;
		bunemanJy[previ+1][prevj][0] += Fy1*Wx1/cellVolume;

		bunemanJz[previ][prevj][0] = Fz1*(1 - Wx1)*(1 - Wy1)/cellVolume;
		bunemanJz[previ][prevj+1][0] = Fz1*(1 - Wx1)*Wy1/cellVolume;
		bunemanJz[previ+1][prevj][0] = Fz1*Wx1*(1 - Wy1)/cellVolume;
		bunemanJz[previ+1][prevj+1][0] = Fz1*Wx1*Wy1/cellVolume;

		bunemanJx[i][j][0] += Fx2*(1.0 - Wy2)/cellVolume;
		bunemanJx[i][j+1][0] += Fx2*Wy2/cellVolume;
		bunemanJy[i][j][0] += Fy2*(1.0 - Wx2)/cellVolume;
		bunemanJy[i+1][j][0] += Fy2*Wx2/cellVolume;

		bunemanJz[i][j][0] = Fz2*(1.0 - Wx2)*(1.0 - Wy2)/cellVolume;
		bunemanJz[i][j+1][0] = Fz2*(1.0 - Wx2)*Wy2/cellVolume;
		bunemanJz[i+1][j][0] = Fz2*Wx2*(1.0 - Wy2)/cellVolume;
		bunemanJz[i+1][j+1][0] = Fz2*Wx2*Wy2/cellVolume;
		
	} else if(dimensionType == TWO_D_XZ) {
		if(i != previ) {
			xr = max2(xgrid[previ], xgrid[i]);
		}
		if(k != prevk) {
			zr = max2(zgrid[prevk], zgrid[k]);
		}

		double Fx1 = fullCharge*(xr - prevLocalCoordinates.x)/deltaT;
		double Fy1 = fullCharge*(yr - prevLocalCoordinates.y)/deltaT;
		double Fz1 = fullCharge*(zr - prevLocalCoordinates.z)/deltaT;
		double Fx2 = fullCharge*velocity.x - Fx1;
		double Fy2 = fullCharge*velocity.y - Fy1;
		double Fz2 = fullCharge*velocity.z - Fz1;

		double Wx1 = (prevLocalCoordinates.x + xr)/(2*deltaX) - previ;
		double Wz1 = (prevLocalCoordinates.z + zr)/(2*deltaZ) - prevk;
		double Wx2 = (localCoordinates.x + xr)/(2*deltaX) - i;
		double Wz2 = (localCoordinates.z + zr)/(2*deltaZ) - k;

		bunemanJx[previ][0][prevk] += Fx1*(1.0 - Wz1)/cellVolume;
		bunemanJx[previ][0][prevk+1] += Fx1*Wz1/cellVolume;
		bunemanJz[previ][0][prevk] += Fy1*(1.0 - Wx1)/cellVolume;
		bunemanJz[previ+1][0][prevk] += Fy1*Wx1/cellVolume;

		bunemanJy[previ][0][prevk] = Fz1*(1 - Wx1)*(1 - Wz1)/cellVolume;
		bunemanJy[previ][0][prevk+1] = Fz1*(1 - Wx1)*Wz1/cellVolume;
		bunemanJy[previ+1][0][prevk] = Fz1*Wx1*(1 - Wz1)/cellVolume;
		bunemanJy[previ+1][0][prevj+k] = Fz1*Wx1*Wz1/cellVolume;

		bunemanJx[i][0][k] += Fx2*(1.0 - Wz2)/cellVolume;
		bunemanJx[i][0][k+1] += Fx2*Wz2/cellVolume;
		bunemanJz[i][0][k] += Fy2*(1.0 - Wx2)/cellVolume;
		bunemanJz[i+1][0][k] += Fy2*Wx2/cellVolume;

		bunemanJy[i][0][k] = Fz2*(1.0 - Wx2)*(1.0 - Wz2)/cellVolume;
		bunemanJy[i][0][k+1] = Fz2*(1.0 - Wx2)*Wz2/cellVolume;
		bunemanJy[i+1][0][k] = Fz2*Wx2*(1.0 - Wz2)/cellVolume;
		bunemanJy[i+1][0][k+1] = Fz2*Wx2*Wz2/cellVolume;
	} else if(dimensionType == ONE_D) {
		if(i != previ) {
			xr = max2(xgrid[previ], xgrid[i]);
		}

		double Fx1 = fullCharge*(xr - prevLocalCoordinates.x)/deltaT;
		double Fy1 = fullCharge*(yr - prevLocalCoordinates.y)/deltaT;
		double Fz1 = fullCharge*(zr - prevLocalCoordinates.z)/deltaT;
		double Fx2 = fullCharge*velocity.x - Fx1;
		double Fy2 = fullCharge*velocity.y - Fy1;
		double Fz2 = fullCharge*velocity.z - Fz1;

		double Wx1 = (prevLocalCoordinates.x + xr)/(2*deltaX) - previ;
		double Wx2 = (localCoordinates.x + xr)/(2*deltaX) - i;

		bunemanJx[previ][0][0] += Fx1/cellVolume;

		bunemanJy[previ][0][0] += Fy1*(1.0 - Wx1)/cellVolume;
		bunemanJy[previ+1][0][0] += Fy1*Wx1/cellVolume;

		bunemanJz[previ][0][0] = Fz1*(1 - Wx1)/cellVolume;
		bunemanJz[previ+1][0][0] = Fz1*Wx1/cellVolume;

		bunemanJx[i][0][0] += Fx2/cellVolume;

		bunemanJy[i][0][0] += Fy2*(1.0 - Wx2)/cellVolume;
		bunemanJy[i+1][0][0] += Fy2*Wx2/cellVolume;

		bunemanJz[i][0][0] = Fz2*(1 - Wx2)/cellVolume;
		bunemanJz[i+1][0][0] = Fz2*Wx2/cellVolume;
	}
}

void Simulation::updateElectroMagneticParameters() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	if (solverType == BUNEMAN) {
		updateBunemanFlux();
	} else {
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

					if(solverType == IMPLICIT_EC){
					for (int tempI = 0; tempI < 2 * splineOrder + 3; ++tempI) {
						for (int tempJ = 0; tempJ < 2 * splineOrder + 3; ++tempJ) {
							for (int tempK = 0; tempK < 2 * splineOrder + 3; ++tempK) {
								for (int curI = 0; curI < 3; ++curI) {
									for (int curJ = 0; curJ < 3; ++curJ) {
										massMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[curI][curJ] = 0;
										tempMassMatrix[i][j][k].matrix[tempI][tempJ][tempK].matrix[curI][curJ] = 0;
									}
								}
							}
						}
					}
					}
				}
			}
		}

		int crossBinNumberX = splineOrder + 2;
		int crossBinNumberY = splineOrder + 2;
		int crossBinNumberZ = splineOrder + 2;

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
						int curK = particle->correlationMapNode.zindex[k];

						double correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.ycorrelation[j] *
							particle->correlationMapNode.zcorrelation[k] / volumeE();
						double particleCharge = particle->charge * particle->weight;
						if (curI <= additionalBinNumber && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
							Particle tempParticle = *particle;
							tempParticle.reflectMomentumX();
							Vector3d tempVelocity = velocity;
							tempVelocity.x = -tempVelocity.x;
							//tempParticle.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, gamma, oldE, oldB);
							Vector3d tempRotatedVelocity = tempParticle.rotationTensor * (tempVelocity * gamma);
							if (solverType == IMPLICIT || solverType == IMPLICIT_EC) {
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

							if (solverType == IMPLICIT || solverType == IMPLICIT_EC) {
								if (solverType == IMPLICIT_EC) {
									for (int tempI = 0; tempI < splineOrder + 2; ++tempI) {
										for (int tempJ = 0; tempJ < splineOrder + 2; ++tempJ) {
											for (int tempK = 0; tempK < splineOrder + 2; ++tempK) {
												int xindex = particle->correlationMapNode.xindex[tempI];
												int yindex = particle->correlationMapNode.yindex[tempJ];
												int zindex = particle->correlationMapNode.zindex[tempK];

												double tempCorrelation1 = particle->correlationMapNode.xcorrelation[tempI] * particle->correlationMapNode.
												                              ycorrelation[tempJ] * particle->correlationMapNode.zcorrelation[tempK];

												for (int tempI1 = 0; tempI1 < splineOrder + 2; ++tempI1) {
													for (int tempJ1 = 0; tempJ1 < splineOrder + 2; ++tempJ1) {
														for (int tempK1 = 0; tempK1 < splineOrder + 2; ++tempK1) {
															double tempCorrelation2 = particle->correlationMapNode.xcorrelation[tempI1] * particle->correlationMapNode.
															                              ycorrelation[tempJ1]*particle->correlationMapNode.zcorrelation[tempK1];
															double doubleCorrelation = tempCorrelation1 * tempCorrelation2 / volumeE();
															for (int l = 0; l < 3; ++l) {
																for (int m = 0; m < 3; ++ m) {
																	massMatrix[xindex][yindex][zindex].matrix[splineOrder + 1 - tempI + tempI1][splineOrder + 1 - tempJ + tempJ1][
																			splineOrder + 1 - tempK + tempK1].matrix[l][m] += particle->rotationTensor.matrix[l][m] *particleOmega * doubleCorrelation;
																}
															}
														}
													}
												}
											}
										}
									}
								}

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
						int curK = particle->correlationMapCell.zindex[k];

						double correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.ycorrelation[j] *
							particle->correlationMapCell.zcorrelation[k] / volumeB();
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

		if (solverType == IMPLICIT_EC) {
			sumNodeMassMatrixParametersX();
			sumNodeMassMatrixParametersY();
			sumNodeMassMatrixParametersZ();
		}

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

		if (solverType == IMPLICIT || solverType == IMPLICIT_EC) {
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

		if (solverType == IMPLICIT || solverType == IMPLICIT_EC) {
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
	}
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

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		for (int i = 0; i < splineOrder + 2; ++i) {
			//todo
			for (int j = 0; j < splineOrder + 2; ++j) {
				for (int k = 0; k < splineOrder + 2; ++k) {
					int curI = particle->correlationMapCell.xindex[i];
					int curJ = particle->correlationMapCell.yindex[j];
					int curK = particle->correlationMapCell.zindex[k];
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
					double correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapCell.ycorrelation[j] *
						particle->correlationMapCell.zcorrelation[k] / volumeB();
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
						particleBulkVelocities[t][i][j][k] = particleBulkVelocities[t][i][j][k] * (speed_of_light_normalized_sqr /
							particleEnergies[t][i][j][k]);
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

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		for (int i = 0; i < splineOrder + 2; ++i) {
			//todo
			for (int j = 0; j < splineOrder + 2; ++j) {
				for (int k = 0; k < splineOrder + 2; ++k) {
					int curI = particle->correlationMapNode.xindex[i];
					int curJ = particle->correlationMapNode.yindex[j];
					int curK = particle->correlationMapNode.zindex[k];
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
					double correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.ycorrelation[j] *
						particle->correlationMapNode.zcorrelation[k] / volumeB();
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

void Simulation::smoothDensity() {
	//todo
	for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
		for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
			double newLeftDensity = (chargeDensityHat[0][j][k] + chargeDensityHat[1][j][k]) / 2.0;
			double newRightDensity = (chargeDensityHat[xnumberAdded - 1][j][k] + chargeDensityHat[xnumberAdded - 2][j][k]) / 2.0;
			if (boundaryConditionTypeX == PERIODIC) {
				newLeftDensity = (chargeDensityHat[xnumberAdded - 1][j][k] + 2.0 * chargeDensityHat[0][j][k] + chargeDensityHat[1][j
				][k]) / 4.0;
				newRightDensity = (chargeDensityHat[xnumberAdded - 2][j][k] + 2.0 * chargeDensityHat[xnumberAdded - 1][j][k] +
					chargeDensityHat[0][j][k]) / 4.0;
			}
			double prevDensity = chargeDensityHat[0][j][k];
			for (int i = 1; i < xnumberAdded - 1; ++i) {
				double tempDensity = chargeDensityHat[i][j][k];

				chargeDensityHat[i][j][k] = (prevDensity + 2.0 * chargeDensityHat[i][j][k] + chargeDensityHat[i + 1][j][k]) / 4.0;

				prevDensity = tempDensity;
			}

			chargeDensityHat[0][j][k] = newLeftDensity;
			chargeDensityHat[xnumberAdded - 1][j][k] = newRightDensity;
		}
	}
}

void Simulation::updateExternalFlux() {
	double alfvenV;
	if (density <= 1E-100) {
		alfvenV = speed_of_light_normalized;
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

void Simulation::updateBunemanFlux() {
	double procTime = 0;
	if (timing && (rank == 0) && (currentIteration % writeParameter == 0)) {
		procTime = clock();
	}
	//MPI_Barrier(cartComm);
	if ((rank == 0) && (verbosity > 0)) printf("updating buneman flux\n");
	if ((rank == 0) && (verbosity > 0)) printLog("updating buneman flux\n");
	int particlePartsCount = 0;
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanJx[i][j][k] = 0;
			}
		}
	}
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				bunemanJy[i][j][k] = 0;
			}
		}
	}
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				bunemanJz[i][j][k] = 0;
			}
		}
	}

	int crossBinNumberX = splineOrder + 2;
	int crossBinNumberY = splineOrder + 2;
	int crossBinNumberZ = splineOrder + 2;

	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		Vector3d uVelocity = (particle->getMomentum() + particle->prevMomentum) / (2.0 * particle->mass);
		double tempGamma = sqrt(1.0 + uVelocity.scalarMult(uVelocity) / speed_of_light_normalized_sqr);
		Vector3d velocity = uVelocity / tempGamma;
		double particleCharge = particle->charge * particle->weight;
		for (int i = 0; i < crossBinNumberX; ++i) {
			for (int j = 0; j < crossBinNumberY; ++j) {
				for (int k = 0; k < crossBinNumberZ; ++k) {
					int curI = particle->correlationMapCell.xindex[i];
					int curJ = particle->correlationMapNode.yindex[j];
					int curK = particle->correlationMapNode.zindex[k];

					double correlation = particle->correlationMapCell.xcorrelation[i] * particle->correlationMapNode.ycorrelation[j] *
						particle->correlationMapNode.zcorrelation[k] / volumeE();

					if (curI <= additionalBinNumber && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
						velocity.x = -velocity.x;
						bunemanJx[2 + 2 * additionalBinNumber - curI][curJ][curK] += velocity.x * (particleCharge * correlation);
					} else {
						bunemanJx[curI][curJ][curK] += velocity.x * (particleCharge * correlation);
					}

					curI = particle->correlationMapNode.xindex[i];
					curJ = particle->correlationMapCell.yindex[j];
					curK = particle->correlationMapNode.zindex[k];

					correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapCell.ycorrelation[j] *
						particle->correlationMapNode.zcorrelation[k] / volumeE();
					if (curI <= additionalBinNumber && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
						bunemanJy[2 + 2 * additionalBinNumber - curI][curJ][curK] += velocity.y * (particleCharge * correlation);
					} else {
						bunemanJy[curI][curJ][curK] += velocity.y * (particleCharge * correlation);
					}

					curI = particle->correlationMapNode.xindex[i];
					curJ = particle->correlationMapNode.yindex[j];
					curK = particle->correlationMapCell.zindex[k];

					correlation = particle->correlationMapNode.xcorrelation[i] * particle->correlationMapNode.ycorrelation[j] *
						particle->correlationMapCell.zcorrelation[k] / volumeE();
					if (curI <= additionalBinNumber && boundaryConditionTypeX == SUPER_CONDUCTOR_LEFT && cartCoord[0] == 0) {
						bunemanJz[2 + 2 * additionalBinNumber - curI][curJ][curK] += velocity.z * (particleCharge * correlation);
					} else {
						bunemanJz[curI][curJ][curK] += velocity.z * (particleCharge * correlation);
					}
				}
			}
		}
	}

	exchangeBunemanFlux();

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				electricFlux[i][j][k] = getBunemanFlux(i, j, k);
			}
		}
	}
}

Vector3d Simulation::getBunemanFlux(int i, int j, int k) {
	return Vector3d(bunemanJx[i][j][k], bunemanJy[i][j][k], bunemanJz[i][j][k]);
}
