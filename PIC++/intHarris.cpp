#include <math.h>
#include <stdio.h>

#include "particle.h"
#include "simulation.h"
#include "vector3d.h"
#include "boundaryFieldEvaluator.h"
#include "util.h"

void Simulation::initializeHarris() {
	//boundaryConditionTypeX = SUPER_CONDUCTOR_LEFT;
	//boundaryConditionTypeX = FREE_BOTH;
	boundaryConditionTypeX = FREE_MIRROR_BOTH;
	double l = 4*deltaX;
    //E0 = Vector3d(0, 0, 0);
	//B0 = Vector3d(0, 0, 0);
	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				Bfield[i][j][k] = B0*tanh((xgrid[i] - 1.5*xsizeGeneral)/l);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
	createParticlesHarris(l);
	E0 = Vector3d(0, 0, 0);

	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded + 1; ++j) {
			for (int k = 0; k < znumberAdded + 1; ++k) {
				Efield[i][j][k] = E0;
				newEfield[i][j][k] = Efield[i][j][k];
			}
		}
	}

	rightBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(E0, B0);
	Vector3d B1 = B0*(-1.0);
	leftBoundaryFieldEvaluator = new ConstantBoundaryFieldEvaluator(E0, B1);

}

void Simulation::createParticlesHarris(double harrisWidth) {
	evaluateParticleTypesAlpha();
	if (rank == 0) printf("creating particles harris\n");
	fflush(stdout);
	if (rank == 0) printLog("creating particles harris\n");
	double n0 = types[0].concentration;
	int n = 0;
	double ionTemperatureFactor = 5;
	for(int t = 1; t < typesNumber; ++t) {
		types[t].temperatureX = types[0].temperatureX*ionTemperatureFactor;
		types[t].temperatureY = types[0].temperatureY*ionTemperatureFactor;
		types[t].temperatureZ = types[0].temperatureZ*ionTemperatureFactor;
	}
	//int density = ; %задаем значение плотности
	//for (int i = 0; i < xnumber; ++i) {
	for (int i = 1 + additionalBinNumber; i < xnumberAdded - additionalBinNumber - 1; ++i) {
		for (int j = 1 + additionalBinNumber; j < ynumberAdded - additionalBinNumber - 1; ++j) {
			for (int k = 1 + additionalBinNumber; k < znumberAdded - additionalBinNumber - 1; ++k) {
				//int maxParticlesPerBin = types[0].particlesPerBin;
				double x = xgrid[i] + 0.0001 * deltaX;
				double y = ygrid[j] + 0.0001 * deltaY;
				double z = zgrid[k] + 0.0001 * deltaZ;
				double ch = cosh((xgrid[i] - 1.5*xsizeGeneral)/harrisWidth);
				double concentr = n0/(ch*ch);
				//for (int l = 0; l < maxParticlesPerBin; ++l) {
				for (int typeCounter = 0; typeCounter < typesNumber; ++typeCounter) {
					double weight = (concentr* volumeB()) / types[typeCounter].particlesPerBin ;
					double deltaXParticles = types[typeCounter].particesDeltaX;
					double deltaYParticles = types[typeCounter].particesDeltaY;
					double deltaZParticles = types[typeCounter].particesDeltaZ;
					//if (l < types[typeCounter].particlesPerBin) {
					for (int l = 0; l < types[typeCounter].particlesPerBin; ++l) {
						ParticleTypes type = types[typeCounter].type;
						Particle* particle = createParticle(n, i, j, k, weight, type, types[typeCounter],
						                                    types[typeCounter].temperatureX,
						                                    types[typeCounter].temperatureY,
						                                    types[typeCounter].temperatureZ);
						n++;
						particle->coordinates.x = x + deltaXParticles * l;
						particle->coordinates.y = y + deltaYParticles * l;
						particle->coordinates.z = z + deltaZParticles * l;
						//particle->coordinates.x = middleXgrid[i];
						//particle->coordinates.y = middleYgrid[j];
						//particle->coordinates.z = middleZgrid[k];

						//particle->coordinates.x = xgrid[i] + uniformDistribution()*deltaX;
						//particle->coordinates.y = ygrid[j] + uniformDistribution()*deltaY;
						//particle->coordinates.z = zgrid[k] + uniformDistribution()*deltaZ;
						particle->initialCoordinates = particle->coordinates;

						double P=harrisWidth*B0.norm()*types[typeCounter].charge;
						Vector3d V = Vector3d(0, 1.0, 0) * (-2*kBoltzman_normalized*speed_of_light_normalized*types[typeCounter].temperatureX/P);

							particle->addVelocity(V, speed_of_light_normalized);

						Vector3d momentum = particle->getMomentum();
						particle->initialMomentum = momentum;
						particle->prevMomentum = momentum;
						particles.push_back(particle);
						particlesNumber++;
						if (particlesNumber % 1000 == 0) {
							if ((rank == 0) && (verbosity > 0))printf("create particle number %d\n", particlesNumber);
						}
						alertNaNOrInfinity(particle->coordinates.x, "particle.x = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.y, "particle.y = NaN in createParticles\n");
						alertNaNOrInfinity(particle->coordinates.z, "particle.z = NaN in createParticles\n");
					}
				}
			}
		}
	}

}