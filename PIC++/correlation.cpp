#include <cmath>

#include "util.h"
#include "simulation.h"
#include "constants.h"

void Simulation::collectParticlesIntoBins() {

	for (int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				particlesInBbin[i][j][k].clear();
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				particlesInEbin[i][j][k].clear();
			}
		}
	}
	double fullSum = 0;
	int pcount = 0;

	#pragma omp parallel for private(pcount) 
	for (pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		checkParticleInBox(*particle);

		int xcount = floor((particle->coordinates.x - xgrid[0]) / deltaX);
		int ycount = floor((particle->coordinates.y - ygrid[0]) / deltaY);
		int zcount = floor((particle->coordinates.z - zgrid[0]) / deltaZ);

		bool crossPrevX = particleCrossBbinX(*particle, xcount - 1);
		bool crossNextX = particleCrossBbinX(*particle, xcount + 1);

		bool crossPrevY = particleCrossBbinY(*particle, ycount - 1) && (ynumber > 1);
		bool crossNextY = particleCrossBbinY(*particle, ycount + 1) && (ynumber > 1);

		bool crossPrevZ = particleCrossBbinZ(*particle, zcount - 1) && (znumber > 1);
		bool crossNextZ = particleCrossBbinZ(*particle, zcount + 1) && (znumber > 1);

		pushParticleIntoBbin(particle, xcount, ycount, zcount);

		if(crossPrevY){
			pushParticleIntoBbin(particle, xcount, ycount - 1, zcount);
		}

		if(crossPrevZ){
			pushParticleIntoBbin(particle, xcount, ycount, zcount - 1);
		}

		if(crossNextY){
			pushParticleIntoBbin(particle, xcount, ycount + 1, zcount);
		}

		if(crossNextZ){
			pushParticleIntoBbin(particle, xcount, ycount, zcount + 1);
		}

		if(crossPrevY && crossPrevZ){
			pushParticleIntoBbin(particle, xcount, ycount - 1, zcount - 1);
		}

		if(crossNextY && crossPrevZ){
			pushParticleIntoBbin(particle, xcount, ycount + 1, zcount - 1);
		}

		if(crossPrevY && crossNextZ){
			pushParticleIntoBbin(particle, xcount, ycount - 1, zcount + 1);
		}

		if(crossNextY && crossNextZ){
			pushParticleIntoBbin(particle, xcount, ycount + 1, zcount + 1);
		}

		if(crossPrevX){
			pushParticleIntoBbin(particle, xcount - 1, ycount, zcount);

			if(crossPrevY){
				pushParticleIntoBbin(particle, xcount - 1, ycount - 1, zcount);
			}

			if(crossPrevZ){
				pushParticleIntoBbin(particle, xcount - 1, ycount, zcount - 1);
			}

			if(crossNextY){
				pushParticleIntoBbin(particle, xcount - 1, ycount + 1, zcount);
			}

			if(crossNextZ){
				pushParticleIntoBbin(particle, xcount - 1, ycount, zcount + 1);
			}

			if(crossPrevY && crossPrevZ){
				pushParticleIntoBbin(particle, xcount - 1, ycount - 1, zcount - 1);
			}

			if(crossNextY && crossPrevZ){
				pushParticleIntoBbin(particle, xcount - 1, ycount + 1, zcount - 1);
			}

			if(crossPrevY && crossNextZ){
				pushParticleIntoBbin(particle, xcount - 1, ycount - 1, zcount + 1);
			}

			if(crossNextY && crossNextZ){
				pushParticleIntoBbin(particle, xcount - 1, ycount + 1, zcount + 1);
			}
		}

		if(crossNextX){
			pushParticleIntoBbin(particle, xcount + 1, ycount, zcount);

			if(crossPrevY){
				pushParticleIntoBbin(particle, xcount + 1, ycount - 1, zcount);
			}

			if(crossPrevZ){
				pushParticleIntoBbin(particle, xcount + 1, ycount, zcount - 1);
			}

			if(crossNextY){
				pushParticleIntoBbin(particle, xcount + 1, ycount + 1, zcount);
			}

			if(crossNextZ){
				pushParticleIntoBbin(particle, xcount + 1, ycount, zcount + 1);
			}

			if(crossPrevY && crossPrevZ){
				pushParticleIntoBbin(particle, xcount + 1, ycount - 1, zcount - 1);
			}

			if(crossNextY && crossPrevZ){
				pushParticleIntoBbin(particle, xcount + 1, ycount + 1, zcount - 1);
			}

			if(crossPrevY && crossNextZ){
				pushParticleIntoBbin(particle, xcount + 1, ycount - 1, zcount + 1);
			}

			if(crossNextY && crossNextZ){
				pushParticleIntoBbin(particle, xcount + 1, ycount + 1, zcount + 1);
			}
		}

		xcount = floor(((particle->coordinates.x - xgrid[0]) / deltaX) + 0.5);
		ycount = floor(((particle->coordinates.y - ygrid[0]) / deltaY) + 0.5);
		zcount = floor(((particle->coordinates.z - zgrid[0]) / deltaZ) + 0.5);

		int prevX = xcount - 1;
		if(prevX == -1  && boundaryConditionType == PERIODIC){
			prevX = xnumber - 1;
		}
		int nextX = xcount + 1;

		if(nextX == xnumber + 1 && boundaryConditionType == PERIODIC){
			nextX = 1;
		}

		int prevY = ycount - 1;
		if(prevY == -1){
			prevY = ynumber - 1;
		}
		int nextY = ycount + 1;
		if(nextY == ynumber + 1){
			nextY = 1;
		}

		int prevZ = zcount - 1;
		if(prevZ == -1){
			prevZ = znumber - 1;
		}
		int nextZ = zcount + 1;
		if(nextZ == znumber + 1){
			nextZ = 1;
		}

		crossPrevX = particleCrossEbinX(*particle, prevX);
		crossNextX = particleCrossEbinX(*particle, nextX);

		crossPrevY = particleCrossEbinY(*particle, prevY) && (ynumber > 1);
		crossNextY = particleCrossEbinY(*particle, nextY) && (ynumber > 1);

		crossPrevZ = particleCrossEbinZ(*particle, prevZ) && (znumber > 1);
		crossNextZ = particleCrossEbinZ(*particle, nextZ) && (znumber > 1);

		pushParticleIntoEbin(particle, xcount, ycount, zcount);

		if(crossPrevY){
			pushParticleIntoEbin(particle, xcount, ycount - 1, zcount);
		}

		if(crossPrevZ){
			pushParticleIntoEbin(particle, xcount, ycount, zcount - 1);
		}

		if(crossNextY){
			pushParticleIntoEbin(particle, xcount, ycount + 1, zcount);
		}

		if(crossNextZ){
			pushParticleIntoEbin(particle, xcount, ycount, zcount + 1);
		}

		if(crossPrevY && crossPrevZ){
			pushParticleIntoEbin(particle, xcount, ycount - 1, zcount - 1);
		}

		if(crossNextY && crossPrevZ){
			pushParticleIntoEbin(particle, xcount, ycount + 1, zcount - 1);
		}

		if(crossPrevY && crossNextZ){
			pushParticleIntoEbin(particle, xcount, ycount - 1, zcount + 1);
		}

		if(crossNextY && crossNextZ){
			pushParticleIntoEbin(particle, xcount, ycount + 1, zcount + 1);
		}

		if(crossPrevX){
			pushParticleIntoEbin(particle, xcount - 1, ycount, zcount);

			if(crossPrevY){
				pushParticleIntoEbin(particle, xcount - 1, ycount - 1, zcount);
			}

			if(crossPrevZ){
				pushParticleIntoEbin(particle, xcount - 1, ycount, zcount - 1);
			}

			if(crossNextY){
				pushParticleIntoEbin(particle, xcount - 1, ycount + 1, zcount);
			}

			if(crossNextZ){
				pushParticleIntoEbin(particle, xcount - 1, ycount, zcount + 1);
			}

			if(crossPrevY && crossPrevZ){
				pushParticleIntoEbin(particle, xcount - 1, ycount - 1, zcount - 1);
			}

			if(crossNextY && crossPrevZ){
				pushParticleIntoEbin(particle, xcount - 1, ycount + 1, zcount - 1);
			}

			if(crossPrevY && crossNextZ){
				pushParticleIntoEbin(particle, xcount - 1, ycount - 1, zcount + 1);
			}

			if(crossNextY && crossNextZ){
				pushParticleIntoEbin(particle, xcount - 1, ycount + 1, zcount + 1);
			}
		}

		if(crossNextX){
			pushParticleIntoEbin(particle, xcount + 1, ycount, zcount);

			if(crossPrevY){
				pushParticleIntoEbin(particle, xcount + 1, ycount - 1, zcount);
			}

			if(crossPrevZ){
				pushParticleIntoEbin(particle, xcount + 1, ycount, zcount - 1);
			}

			if(crossNextY){
				pushParticleIntoEbin(particle, xcount + 1, ycount + 1, zcount);
			}

			if(crossNextZ){
				pushParticleIntoEbin(particle, xcount + 1, ycount, zcount + 1);
			}

			if(crossPrevY && crossPrevZ){
				pushParticleIntoEbin(particle, xcount + 1, ycount - 1, zcount - 1);
			}

			if(crossNextY && crossPrevZ){
				pushParticleIntoEbin(particle, xcount + 1, ycount + 1, zcount - 1);
			}

			if(crossPrevY && crossNextZ){
				pushParticleIntoEbin(particle, xcount + 1, ycount - 1, zcount + 1);
			}

			if(crossNextY && crossNextZ){
				pushParticleIntoEbin(particle, xcount + 1, ycount + 1, zcount + 1);
			}
		}
	}
}

void Simulation::pushParticleIntoEbin(Particle* particle, int i, int j, int k) {
	if (i < 0) return;
	if (i > xnumber) return;

	if (j < 0) return;
	if (j > ynumber) return;

	if (k < 0) return;
	if (k > znumber) return;

	particlesInEbin[i][j][k].push_back(particle);
}

void Simulation::pushParticleIntoBbin(Particle* particle, int i, int j, int k) {
	if (i < 0){
		if(boundaryConditionType == PERIODIC){
			i = xnumber - 1;
		} else {
			return;
		}
	}
	if (i >= xnumber){
		if(boundaryConditionType == PERIODIC){
			i = 0;
		} else {
			return;
		}
	}

	if(j < 0){
		j = ynumber - 1;
	}
	if(j >= ynumber){
		j = 0;
	}
	if(k < 0){
		k = znumber -1;
	}
	if(k >= znumber){
		k = 0;
	}

	particlesInBbin[i][j][k].push_back(particle);
}

bool Simulation::particleCrossBbinX(Particle& particle, int i) {
	if(boundaryConditionType == PERIODIC){
		if(i < 0) {
			i = xnumber - 1;
		} else if(i >= xnumber) {
			i = 0;
		}	

		if (i == 0) {
			if ((xgrid[i + 1] < particle.coordinates.x - particle.dx) && (xgrid[xnumber] > particle.coordinates.x + particle.dx))
				return false;
		} else if (i == xnumber - 1) {
			if ((xgrid[i] > particle.coordinates.x + particle.dx) && (xgrid[0] < particle.coordinates.x - particle.dx))
				return false;
		} else {
			if ((xgrid[i] > particle.coordinates.x + particle.dx) || (xgrid[i + 1] < particle.coordinates.x - particle.dx))
				return false;
		}

		return true;
	} else if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
		if(i < 0) {
			return false;
		} else if(i >= xnumber) {
			return false;
		}	

		if (i == 0) {
			if (xgrid[i + 1] < particle.coordinates.x - particle.dx)
				return false;
		} else if (i == xnumber - 1) {
			if (xgrid[i] > particle.coordinates.x + particle.dx)
				return false;
		} else {
			if ((xgrid[i] > particle.coordinates.x + particle.dx) || (xgrid[i + 1] < particle.coordinates.x - particle.dx))
				return false;
		}
		return true;

	}
	return false;
}

bool Simulation::particleCrossBbinY(Particle& particle, int j) {
	if(j < 0) {
		j = ynumber - 1;
	} else if(j >= ynumber) {
		j = 0;
	}	
	if (j == 0) {
		if ((ygrid[j + 1] < particle.coordinates.y - particle.dy) && (ygrid[ynumber] > particle.coordinates.y + particle.dy))
			return false;
	} else if (j == ynumber - 1) {
		if ((ygrid[j] > particle.coordinates.y + particle.dy) && (ygrid[0] < particle.coordinates.y - particle.dy))
			return false;
	} else {
		if ((ygrid[j] > particle.coordinates.y + particle.dy) || (ygrid[j + 1] < particle.coordinates.y - particle.dy))
			return false;
	}
	return true;
}

bool Simulation::particleCrossBbinZ(Particle& particle, int k) {
	if(k < 0) {
		k = znumber - 1;
	} else if(k >= znumber) {
		k = 0;
	}	
	if (k == 0) {
		if ((zgrid[k + 1] < particle.coordinates.z - particle.dz) && (zgrid[xnumber] > particle.coordinates.z + particle.dz))
			return false;
	} else if (k == znumber - 1) {
		if ((zgrid[k] > particle.coordinates.z + particle.dz) && (zgrid[0] < particle.coordinates.z - particle.dz))
			return false;
	} else {
		if ((zgrid[k] > particle.coordinates.z + particle.dz) || (zgrid[k + 1] < particle.coordinates.z - particle.dz))
			return false;
	}

	return true;
}

bool Simulation::particleCrossBbin(Particle& particle, int i, int j, int k) {
	bool crossX = particleCrossBbinX(particle, i);

	if(! crossX){
		return false;
	}

	bool crossY = particleCrossBbinY(particle, j);

	if(! crossY){
		return false;
	}

	return particleCrossBbinZ(particle, k);
}

bool Simulation::particleCrossEbinX(Particle& particle, int i){
	if(boundaryConditionType == PERIODIC){
		if(i == 0){
			if(xgrid[0] + (deltaX/2) < particle.coordinates.x - particle.dx && xgrid[xnumber] - (deltaX/2) > particle.coordinates.x + particle.dx)
				return false;
		} else if(i == xnumber){
			if(xgrid[0] + (deltaX/2) < particle.coordinates.x - particle.dx && xgrid[xnumber] - (deltaX/2) > particle.coordinates.x + particle.dx)
				return false;
		} else {
			if ((xgrid[i] - (deltaX / 2) > particle.coordinates.x + particle.dx) || (xgrid[i + 1] - (deltaX / 2) < particle.coordinates.x - particle.dx))
				return false;
		}
		return true;
	} else if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
		if(i == 0){
			if(xgrid[0] + (deltaX/2) < particle.coordinates.x - particle.dx)
				return false;
		} else if(i == xnumber){
			if(xgrid[xnumber] - (deltaX/2) > particle.coordinates.x + particle.dx)
				return false;
		} else {
			if ((xgrid[i] - (deltaX / 2) > particle.coordinates.x + particle.dx) || (xgrid[i + 1] - (deltaX / 2) < particle.coordinates.x - particle.dx))
				return false;
		}
		return true;
	}

	return false;
}

bool Simulation::particleCrossEbinY(Particle& particle, int j){
	if(j == 0){
		if(ygrid[0] + (deltaY/2) < particle.coordinates.y - particle.dy && ygrid[ynumber] - (deltaY/2) > particle.coordinates.y + particle.dy)
			return false;
	} else if(j == ynumber){
		if(ygrid[0] + (deltaY/2) < particle.coordinates.y - particle.dy && ygrid[ynumber] - (deltaY/2) > particle.coordinates.y + particle.dy)
			return false;
	} else {
		if ((ygrid[j] - (deltaY / 2) > particle.coordinates.y + particle.dy) || (ygrid[j + 1] - (deltaY / 2) < particle.coordinates.y - particle.dy))
			return false;
	}

	return true;
}

bool Simulation::particleCrossEbinZ(Particle& particle, int k){
	if(k == 0){
		if(zgrid[0] + (deltaZ/2) < particle.coordinates.z - particle.dz && zgrid[znumber] - (deltaZ/2) > particle.coordinates.z + particle.dz)
			return false;
	} else if(k == znumber){
		if(zgrid[0] + (deltaZ/2) < particle.coordinates.z - particle.dz && zgrid[znumber] - (deltaZ/2) > particle.coordinates.z + particle.dz)
			return false;
	} else {
		if ((zgrid[k] - (deltaZ / 2) > particle.coordinates.z + particle.dz) || (zgrid[k + 1] - (deltaZ / 2) < particle.coordinates.z - particle.dz))
			return false;
	}

	return true;
}

bool Simulation::particleCrossEbin(Particle& particle, int i, int j, int k) {
	bool crossX = particleCrossEbinX(particle, i);

	if(! crossX){
		return false;
	}

	bool crossY = particleCrossEbinY(particle, j);

	if(! crossY){
		return false;
	}

	return particleCrossEbinZ(particle, k);
}


Vector3d Simulation::correlationTempEfield(Particle* particle) {
	return correlationTempEfield(*particle);
}

Vector3d Simulation::correlationNewEfield(Particle* particle) {
	return correlationNewEfield(*particle);
}

Vector3d Simulation::correlationBfield(Particle* particle) {
	return correlationBfield(*particle);
}

Vector3d Simulation::correlationNewBfield(Particle* particle) {
	return correlationNewBfield(*particle);
}

Vector3d Simulation::correlationEfield(Particle* particle) {
	return correlationEfield(*particle);
}

Vector3d Simulation::correlationEfield(Particle& particle) {
	//checkParticleInBox(particle);

	/*int xcount = floor((particle.x / deltaX) + 0.5);

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = xcount - 1; i <= xcount + 1; ++i) {
				result = result + correlationFieldWithEbin(particle, i);
	}

	return result;*/

	int xcount = floor((particle.coordinates.x - xgrid[0])/deltaX);
	int ycount = floor((particle.coordinates.y - ygrid[0])/deltaY);
	int zcount = floor((particle.coordinates.z - zgrid[0])/deltaZ);

	double rightWeightX = (particle.coordinates.x - xgrid[xcount])/deltaX;
	double leftWeightX = (xgrid[xcount + 1] - particle.coordinates.x)/deltaX;

	double rightWeightY = (particle.coordinates.y - ygrid[ycount])/deltaY;
	double leftWeightY = (ygrid[ycount + 1] - particle.coordinates.y)/deltaY;

	double rightWeightZ = (particle.coordinates.z - zgrid[zcount])/deltaZ;
	double leftWeightZ = (zgrid[zcount + 1] - particle.coordinates.z)/deltaZ;

	return ((Efield[xcount][ycount][zcount]*leftWeightZ + Efield[xcount][ycount][zcount+1]*rightWeightZ)*leftWeightY + (Efield[xcount][ycount+1][zcount]*leftWeightZ + Efield[xcount][ycount+1][zcount+1]*rightWeightZ)*rightWeightY)*leftWeightX 
		+ ((Efield[xcount+1][ycount][zcount]*leftWeightZ + Efield[xcount+1][ycount][zcount+1]*rightWeightZ)*leftWeightY + (Efield[xcount+1][ycount+1][zcount]*leftWeightZ + Efield[xcount+1][ycount+1][zcount+1]*rightWeightZ)*rightWeightY)*rightWeightX;
}

Vector3d Simulation::correlationTempEfield(Particle& particle) {
	//checkParticleInBox(particle);

	/*int xcount = floor((particle.x / deltaX) + 0.5);

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = xcount - 1; i <= xcount + 1; ++i) {
				result = result + correlationFieldWithTempEbin(particle, i);
	}

	return result;*/
	int xcount = floor((particle.coordinates.x - xgrid[0])/deltaX);
	int ycount = floor((particle.coordinates.y - ygrid[0])/deltaY);
	int zcount = floor((particle.coordinates.z - zgrid[0])/deltaZ);

	double rightWeightX = (particle.coordinates.x - xgrid[xcount])/deltaX;
	double leftWeightX = (xgrid[xcount + 1] - particle.coordinates.x)/deltaX;

	double rightWeightY = (particle.coordinates.y - ygrid[ycount])/deltaY;
	double leftWeightY = (ygrid[ycount + 1] - particle.coordinates.y)/deltaY;

	double rightWeightZ = (particle.coordinates.z - zgrid[zcount])/deltaZ;
	double leftWeightZ = (zgrid[zcount + 1] - particle.coordinates.z)/deltaZ;

	return ((tempEfield[xcount][ycount][zcount]*leftWeightZ + tempEfield[xcount][ycount][zcount+1]*rightWeightZ)*leftWeightY + (tempEfield[xcount][ycount+1][zcount]*leftWeightZ + tempEfield[xcount][ycount+1][zcount+1]*rightWeightZ)*rightWeightY)*leftWeightX 
		+ ((tempEfield[xcount+1][ycount][zcount]*leftWeightZ + tempEfield[xcount+1][ycount][zcount+1]*rightWeightZ)*leftWeightY + (tempEfield[xcount+1][ycount+1][zcount]*leftWeightZ + tempEfield[xcount+1][ycount+1][zcount+1]*rightWeightZ)*rightWeightY)*rightWeightX;
}

Vector3d Simulation::correlationNewEfield(Particle& particle) {
	//checkParticleInBox(particle);

	/*int xcount = floor((particle.x / deltaX) + 0.5);

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = xcount - 1; i <= xcount + 1; ++i) {
				result = result + correlationFieldWithNewEbin(particle, i);
	}

	return result;*/
	int xcount = floor((particle.coordinates.x - xgrid[0])/deltaX);
	int ycount = floor((particle.coordinates.y - ygrid[0])/deltaY);
	int zcount = floor((particle.coordinates.z - zgrid[0])/deltaZ);

	double rightWeightX = (particle.coordinates.x - xgrid[xcount])/deltaX;
	double leftWeightX = (xgrid[xcount + 1] - particle.coordinates.x)/deltaX;

	double rightWeightY = (particle.coordinates.y - ygrid[ycount])/deltaY;
	double leftWeightY = (ygrid[ycount + 1] - particle.coordinates.y)/deltaY;

	double rightWeightZ = (particle.coordinates.z - zgrid[zcount])/deltaZ;
	double leftWeightZ = (zgrid[zcount + 1] - particle.coordinates.z)/deltaZ;

	return ((newEfield[xcount][ycount][zcount]*leftWeightZ + newEfield[xcount][ycount][zcount+1]*rightWeightZ)*leftWeightY + (newEfield[xcount][ycount+1][zcount]*leftWeightZ + newEfield[xcount][ycount+1][zcount+1]*rightWeightZ)*rightWeightY)*leftWeightX 
		+ ((newEfield[xcount+1][ycount][zcount]*leftWeightZ + newEfield[xcount+1][ycount][zcount+1]*rightWeightZ)*leftWeightY + (newEfield[xcount+1][ycount+1][zcount]*leftWeightZ + newEfield[xcount+1][ycount+1][zcount+1]*rightWeightZ)*rightWeightY)*rightWeightX;
}

Vector3d Simulation::correlationBfield(Particle& particle) {

	double x = particle.coordinates.x;
	if(x > xgrid[xnumber]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}
	if(x < xgrid[0]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}

	double y = particle.coordinates.y;
	if(y > ygrid[ynumber]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}
	if(y < ygrid[0]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}

	double z = particle.coordinates.z;
	if(z > zgrid[znumber]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}
	if(z < zgrid[0]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}

	int xcount = floor(((x - xgrid[0]) / deltaX) + 0.5);
	int ycount = floor(((y - ygrid[0]) / deltaY) + 0.5);
	int zcount = floor(((z - zgrid[0]) / deltaZ) + 0.5);

	Vector3d leftField = Vector3d(0, 0, 0);
	Vector3d rightField = Vector3d(0, 0, 0);

	double rightWeightX;
	double leftWeightX;
	double rightWeightY;
	double leftWeightY;
	double rightWeightZ;
	double leftWeightZ;

	int curI = xcount;
	int prevI = xcount - 1;
	int curJ = ycount;
	int prevJ = ycount - 1;
	int curK = zcount;
	int prevK = zcount - 1;


	if(xcount > 0 && xcount < xnumber){
		curI = xcount;
		prevI = xcount - 1;
		rightWeightX = (particle.coordinates.x - middleXgrid[xcount-1])/deltaX;
		leftWeightX = (middleXgrid[xcount] - particle.coordinates.x)/deltaX;

	} else if(xcount == 0){
		curI = xcount;

		leftWeightX = (middleXgrid[xcount] - particle.coordinates.x)/deltaX;
		rightWeightX = 1.0 - leftWeightX;

		if(boundaryConditionType == PERIODIC){
			prevI = xnumber - 1;
		} else {
			prevI = 0;
		}
	} else if(xcount == xnumber){
		rightWeightX = (particle.coordinates.x - middleXgrid[xcount-1])/deltaX;
		leftWeightX = 1 - rightWeightX;
		
		prevI = xcount - 1;
		if(boundaryConditionType == PERIODIC){
			curI = 0;
		} else {
			curI = xcount-1;
		}
	}

	if(ycount > 0 && ycount < ynumber){
		curJ = ycount;
		prevJ = ycount - 1;
		rightWeightY = (particle.coordinates.y - middleYgrid[ycount-1])/deltaY;
		leftWeightY = (middleYgrid[ycount] - particle.coordinates.y)/deltaY;

	} else if(ycount == 0){
		curJ = ycount;
		prevJ = ynumber - 1;

		leftWeightY = (middleYgrid[ycount] - particle.coordinates.y)/deltaY;
		rightWeightY = 1.0 - leftWeightY;

	} else if(ycount == ynumber){
		rightWeightY = (particle.coordinates.y - middleYgrid[ycount-1])/deltaY;
		leftWeightY = 1 - rightWeightY;
		
		prevJ = ycount - 1;
		curJ = 0;
	}

	if(zcount > 0 && zcount < znumber){
		curK = zcount;
		prevK = zcount - 1;
		rightWeightZ = (particle.coordinates.z - middleZgrid[zcount-1])/deltaZ;
		leftWeightZ = (middleZgrid[zcount] - particle.coordinates.z)/deltaZ;

	} else if(zcount == 0){
		curK = zcount;
		prevK = znumber - 1;

		leftWeightZ = (middleZgrid[zcount] - particle.coordinates.z)/deltaZ;
		rightWeightZ = 1.0 - leftWeightZ;

	} else if(zcount == znumber){
		rightWeightZ = (particle.coordinates.z - middleZgrid[zcount-1])/deltaZ;
		leftWeightZ = 1 - rightWeightZ;
		
		prevK = zcount - 1;
		curK = 0;
	}

	return ((Bfield[prevI][prevJ][prevK]*leftWeightZ + Bfield[prevI][prevJ][curK]*rightWeightZ)*leftWeightY + ((Bfield[prevI][curJ][prevK]*leftWeightZ + Bfield[prevI][curJ][curK]*rightWeightZ)*rightWeightY))*leftWeightX
		+ ((Bfield[curI][prevJ][prevK]*leftWeightZ + Bfield[curI][prevJ][curK]*rightWeightZ)*leftWeightY + ((Bfield[curI][curJ][prevK]*leftWeightZ + Bfield[curI][curJ][curK]*rightWeightZ)*rightWeightY))*rightWeightX;
}

Vector3d Simulation::correlationNewBfield(Particle& particle) {
	double x = particle.coordinates.x;
	if(x > xgrid[xnumber]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}
	if(x < xgrid[0]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}

	double y = particle.coordinates.y;
	if(y > ygrid[ynumber]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}
	if(y < ygrid[0]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}

	double z = particle.coordinates.z;
	if(z > zgrid[znumber]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}
	if(z < zgrid[0]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}

	int xcount = floor(((x - xgrid[0]) / deltaX) + 0.5);
	int ycount = floor(((y - ygrid[0]) / deltaY) + 0.5);
	int zcount = floor(((z - xgrid[0]) / deltaZ) + 0.5);

	Vector3d leftField = Vector3d(0, 0, 0);
	Vector3d rightField = Vector3d(0, 0, 0);

	double rightWeightX;
	double leftWeightX;
	double rightWeightY;
	double leftWeightY;
	double rightWeightZ;
	double leftWeightZ;

	int curI = xcount;
	int prevI = xcount - 1;
	int curJ = ycount;
	int prevJ = ycount - 1;
	int curK = zcount;
	int prevK = zcount - 1;


	if(xcount > 0 && xcount < xnumber){
		curI = xcount;
		prevI = xcount - 1;
		rightWeightX = (particle.coordinates.x - middleXgrid[xcount-1])/deltaX;
		leftWeightX = (middleXgrid[xcount] - particle.coordinates.x)/deltaX;

	} else if(xcount == 0){
		curI = xcount;

		leftWeightX = (middleXgrid[xcount] - particle.coordinates.x)/deltaX;
		rightWeightX = 1.0 - leftWeightX;

		if(boundaryConditionType == PERIODIC){
			prevI = xnumber - 1;
		} else {
			prevI = 0;
		}
	} else if(xcount == xnumber){
		rightWeightX = (particle.coordinates.x - middleXgrid[xcount-1])/deltaX;
		leftWeightX = 1 - rightWeightX;
		
		prevI = xcount - 1;
		if(boundaryConditionType == PERIODIC){
			curI = 0;
		} else {
			curI = xcount-1;
		}
	}

	if(ycount > 0 && ycount < ynumber){
		curJ = ycount;
		prevJ = ycount - 1;
		rightWeightY = (particle.coordinates.y - middleYgrid[ycount-1])/deltaY;
		leftWeightY = (middleYgrid[ycount] - particle.coordinates.y)/deltaY;

	} else if(ycount == 0){
		curJ = ycount;
		prevJ = ynumber - 1;

		leftWeightY = (middleYgrid[ycount] - particle.coordinates.y)/deltaY;
		rightWeightY = 1.0 - leftWeightY;

	} else if(ycount == ynumber){
		rightWeightY = (particle.coordinates.y - middleYgrid[ycount-1])/deltaY;
		leftWeightY = 1 - rightWeightY;
		
		prevJ = ycount - 1;
		curJ = 0;
	}

	if(zcount > 0 && zcount < znumber){
		curK = zcount;
		prevK = zcount - 1;
		rightWeightZ = (particle.coordinates.z - middleZgrid[zcount-1])/deltaZ;
		leftWeightZ = (middleZgrid[zcount] - particle.coordinates.z)/deltaZ;

	} else if(zcount == 0){
		curK = zcount;
		prevK = znumber - 1;

		leftWeightZ = (middleZgrid[zcount] - particle.coordinates.z)/deltaZ;
		rightWeightZ = 1.0 - leftWeightZ;

	} else if(zcount == znumber){
		rightWeightZ = (particle.coordinates.z - middleZgrid[zcount-1])/deltaZ;
		leftWeightZ = 1 - rightWeightZ;
		
		prevK = zcount - 1;
		curK = 0;
	}

	return ((newBfield[prevI][prevJ][prevK]*leftWeightZ + newBfield[prevI][prevJ][curK]*rightWeightZ)*leftWeightY + ((newBfield[prevI][curJ][prevK]*leftWeightZ + newBfield[prevI][curJ][curK]*rightWeightZ)*rightWeightY))*leftWeightX
		+ ((newBfield[curI][prevJ][prevK]*leftWeightZ + newBfield[curI][prevJ][curK]*rightWeightZ)*leftWeightY + ((newBfield[curI][curJ][prevK]*leftWeightZ + newBfield[curI][curJ][curK]*rightWeightZ)*rightWeightY))*rightWeightX;
}

double Simulation::correlationWithBbin(Particle& particle, int i, int j, int k) {
	//if (! particleCrossBbin(particle, i, j, k))
		//return 0.0;

	double x = particle.coordinates.x;
	double y = particle.coordinates.y;
	double z = particle.coordinates.z;

	double leftx;
	double rightx;

	double lefty;
	double righty;

	double leftz;
	double rightz;

	if(i == -1) {
		if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
			return 0.0;
		}
		if (particle.coordinates.x - particle.dx > xgrid[0]) {
			leftx = xgrid[xnumber-1];
			rightx = xgrid[xnumber];
		} else {
			leftx = xgrid[0] - deltaX;
			rightx = xgrid[0];
		}	
	} else if(i == 0){
		if(boundaryConditionType == PERIODIC){
			if (particle.coordinates.x - particle.dx > xgrid[1]) {
				leftx = xgrid[xnumber];
				rightx = xgrid[xnumber] + deltaX;
			} else {
				leftx = xgrid[0];
				rightx = xgrid[1];
			}
		} else {
			leftx = -deltaX;
			rightx = xgrid[1];
		}
	} else if(i == xnumber - 1){
		if(boundaryConditionType == PERIODIC){
			if (particle.coordinates.x + particle.dx < xgrid[xnumber - 1]) {
				leftx = xgrid[0] - deltaX;
				rightx = xgrid[0];
			} else {
				leftx = xgrid[xnumber - 1];
				rightx = xgrid[xnumber];
			}
		} else {
			leftx = xgrid[xnumber - 1];
			rightx = xgrid[xnumber];
		}
	} else if(i == xnumber) {
		if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
			return 0.0;
		}
		if (particle.coordinates.x + particle.dx < xgrid[xnumber]) {
			leftx = xgrid[0];
			rightx = xgrid[1];
		} else {
			leftx = xgrid[xnumber];
			rightx = xgrid[xnumber] + deltaX;
		}
	} else {
		leftx = xgrid[i];
		rightx = xgrid[i + 1];
	}

	if(j == -1) {
		if (particle.coordinates.y - particle.dy > ygrid[0]) {
			lefty = ygrid[ynumber-1];
			righty = ygrid[ynumber];
		} else {
			lefty = ygrid[0] - deltaY;
			righty = ygrid[0];
		}	
	} else if(j == 0){
		if (particle.coordinates.y - particle.dy > ygrid[1]) {
			lefty = ygrid[ynumber];
			righty = ygrid[ynumber] + deltaY;
		} else {
			lefty = ygrid[0];
			righty = ygrid[1];
		}
	} else if(j == ynumber - 1){
		if (particle.coordinates.y + particle.dy < ygrid[ynumber - 1]) {
			lefty = ygrid[0] - deltaY;
			righty = ygrid[0];
		} else {
			lefty = ygrid[ynumber - 1];
			righty = ygrid[ynumber];
		}
	} else if(j == ynumber) {
		if (particle.coordinates.y + particle.dy < ygrid[ynumber]) {
			lefty = ygrid[0];
			righty = ygrid[1];
		} else {
			lefty = ygrid[ynumber];
			righty = ygrid[ynumber] + deltaY;
		}
	} else {
		lefty = ygrid[j];
		righty = ygrid[j + 1];
	}

	if(k == -1) {
		if (particle.coordinates.z - particle.dz > zgrid[0]) {
			leftz = zgrid[znumber-1];
			rightz = zgrid[znumber];
		} else {
			leftz = zgrid[0] - deltaZ;
			rightz = zgrid[0];
		}	
	} else if(k == 0){
		if (particle.coordinates.z - particle.dz > zgrid[1]) {
			leftz = ygrid[znumber];
			rightz = ygrid[znumber] + deltaZ;
		} else {
			leftz = zgrid[0];
			rightz = zgrid[1];
		}
	} else if(k == znumber - 1){
		if (particle.coordinates.z + particle.dz < zgrid[znumber - 1]) {
			leftz = zgrid[0] - deltaZ;
			rightz = zgrid[0];
		} else {
			leftz = zgrid[znumber - 1];
			rightz = zgrid[znumber];
		}
	} else if(k == znumber) {
		if (particle.coordinates.z + particle.dz < zgrid[znumber]) {
			leftz = zgrid[0];
			rightz = zgrid[1];
		} else {
			leftz = zgrid[znumber];
			rightz = zgrid[znumber] + deltaZ;
		}
	} else {
		leftz = zgrid[k];
		rightz = zgrid[k + 1];
	}



	double correlationx = correlationBspline(x, particle.dx, leftx, rightx);
	double correlationy = correlationBspline(y, particle.dy, lefty, righty);
	double correlationz = correlationBspline(z, particle.dz, leftz, rightz);

	return correlationx*correlationy*correlationz;
}

double Simulation::correlationWithEbin(Particle& particle, int i, int j, int k) {
	double x = particle.coordinates.x;
	double y = particle.coordinates.y;
	double z = particle.coordinates.z;

	double leftx;
	double rightx;
	double lefty;
	double righty;
	double leftz;
	double rightz;

	if(i == -1) {
		if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
			return 0.0;
		}
		if (particle.coordinates.x - particle.dx > xgrid[0] - deltaX/2) {
			leftx = xgrid[xnumber - 1] - deltaX/2;
			rightx = xgrid[xnumber - 1] + deltaX/2;
		} else {
			leftx = xgrid[0] - 3*deltaX/2;
			rightx = xgrid[0] - deltaX/2;
		}		
	} else if(i == 0){
		if(boundaryConditionType == PERIODIC){
			if (particle.coordinates.x - particle.dx > xgrid[0] + deltaX/2) {
				leftx = xgrid[xnumber] - deltaX/2;
				rightx = xgrid[xnumber] + deltaX/2;
			} else {
				leftx = xgrid[0] - deltaX/2;
				rightx = xgrid[0] + deltaX/2;
			}
		} else {
			leftx = xgrid[0] - deltaX/2;
			rightx = xgrid[0] + deltaX/2;
		}
	} else if (i == xnumber) {
		if(boundaryConditionType == PERIODIC){
			if (particle.coordinates.x + particle.dx < xgrid[xnumber] - deltaX/2) {
				leftx = xgrid[0] - deltaX/2;
				rightx = xgrid[0] + deltaX/2;
			} else {
				leftx = xgrid[xnumber] - deltaX/2;
				rightx = xgrid[xnumber] + deltaX/2;
			}
		} else {
			leftx = xgrid[xnumber] - deltaX/2;
			rightx = xgrid[xnumber] + deltaX/2;
		}
	} else if(i == xnumber + 1) {
		if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
			return 0.0;
		}
		if (particle.coordinates.x + particle.dx < xgrid[xnumber] + deltaX/2) {
			leftx = xgrid[1] - deltaX/2;
			rightx = xgrid[1] + deltaX/2;
		} else {
			leftx = xgrid[xnumber] + deltaX/2;
			rightx = xgrid[xnumber] + 3*deltaX/2;
		}
	} else {
		leftx = xgrid[i] - (deltaX / 2);
		rightx = xgrid[i] + (deltaX / 2);
	}

	if(j == -1) {
		if (particle.coordinates.y - particle.dy > ygrid[0] - deltaY/2) {
			lefty = ygrid[ynumber - 1] - deltaY/2;
			righty = ygrid[ynumber - 1] + deltaY/2;
		} else {
			lefty = ygrid[0] - 3*deltaY/2;
			righty = ygrid[0] - deltaY/2;
		}		
	} else if(j == 0){
			if (particle.coordinates.y - particle.dy > ygrid[0] + deltaY/2) {
				lefty = ygrid[ynumber] - deltaY/2;
				righty = ygrid[ynumber] + deltaY/2;
			} else {
				lefty = ygrid[0] - deltaY/2;
				righty = ygrid[0] + deltaY/2;
			}
	} else if (j == ynumber) {
			if (particle.coordinates.y + particle.dy < ygrid[ynumber] - deltaY/2) {
				lefty = ygrid[0] - deltaY/2;
				righty = ygrid[0] + deltaY/2;
			} else {
				lefty = ygrid[ynumber] - deltaY/2;
				righty = ygrid[ynumber] + deltaY/2;
			}
	} else if(j == ynumber + 1) {
		if (particle.coordinates.y + particle.dy < ygrid[ynumber] + deltaY/2) {
			lefty = ygrid[1] - deltaY/2;
			righty = ygrid[1] + deltaY/2;
		} else {
			lefty = ygrid[ynumber] + deltaY/2;
			righty = ygrid[ynumber] + 3*deltaY/2;
		}
	} else {
		lefty = ygrid[j] - (deltaY / 2);
		righty = ygrid[j] + (deltaY / 2);
	}

	if(k == -1) {
		if (particle.coordinates.z - particle.dz > zgrid[0] - deltaZ/2) {
			leftz = zgrid[ynumber - 1] - deltaZ/2;
			rightz = zgrid[ynumber - 1] + deltaZ/2;
		} else {
			leftz = zgrid[0] - 3*deltaZ/2;
			rightz = zgrid[0] - deltaZ/2;
		}		
	} else if(k == 0){
			if (particle.coordinates.z - particle.dz > zgrid[0] + deltaZ/2) {
				leftz = zgrid[znumber] - deltaZ/2;
				rightz = zgrid[znumber] + deltaZ/2;
			} else {
				leftz = zgrid[0] - deltaZ/2;
				rightz = zgrid[0] + deltaZ/2;
			}
	} else if (k == znumber) {
			if (particle.coordinates.z + particle.dz < zgrid[znumber] - deltaZ/2) {
				leftz = zgrid[0] - deltaZ/2;
				rightz = zgrid[0] + deltaZ/2;
			} else {
				leftz = zgrid[znumber] - deltaZ/2;
				rightz = zgrid[znumber] + deltaZ/2;
			}
	} else if(k == znumber + 1) {
		if (particle.coordinates.z + particle.dz < zgrid[znumber] + deltaZ/2) {
			leftz = zgrid[1] - deltaZ/2;
			rightz = zgrid[1] + deltaZ/2;
		} else {
			leftz = zgrid[znumber] + deltaZ/2;
			rightz = zgrid[znumber] + 3*deltaZ/2;
		}
	} else {
		leftz = zgrid[k] - (deltaZ / 2);
		rightz = zgrid[k] + (deltaZ / 2);
	}

	double correlationx = correlationBspline(x, particle.dx, leftx, rightx);
	double correlationy = correlationBspline(y, particle.dy, lefty, righty);
	double correlationz = correlationBspline(z, particle.dz, leftz, rightz);

	return correlationx*correlationy*correlationz;
}

double Simulation::correlationBspline(const double& x, const double& dx, const double& leftx, const double& rightx) {

	if (rightx < leftx) {
		printf("rightx < leftx\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "rightx = %15.10g < leftx = %15.10g\n", rightx, leftx);
		fclose(errorLogFile);
		exit(0);
	}
	if (dx > rightx - leftx) {
		printf("dx > rightx - leftx\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "dx = %15.10g > rightx - leftx = %15.10g\n", dx, rightx - leftx);
		fclose(errorLogFile);
		exit(0);
	}

	double correlation = 0;

	if (x < leftx - dx)
		return 0;
	if (x > rightx + dx)
		return 0;

	switch (splineOrder){
		case 0:
			if ( x < leftx + deltaX){
				correlation = 0.5*(x + deltaX - leftx)/deltaX;
			} else if( x > rightx - deltaX){
				correlation = 0.5*(rightx - (x - deltaX))/deltaX;
			} else {
				correlation = 1;
			}
			break;
		case 1:
			if( x < leftx){
				correlation = 0.5*sqr(x + deltaX - leftx)/deltaX2;
			} else if (x < leftx + deltaX){
				correlation = 1 - 0.5*sqr(leftx - (x - deltaX))/deltaX2;
			} else if (x > rightx){
				correlation = 0.5*sqr(rightx - (x - deltaX))/deltaX2;
			} else if (x > rightx - deltaX){
				correlation = 1 - 0.5*sqr(x + deltaX - rightx)/deltaX2;
			} else {
				correlation = 1;
			}
			break;
		case 2:
			if (x < leftx - dx/2) {
				correlation = 2*cube(x + dx - leftx)/(3*cube(dx));
			} else if(x < leftx){
				correlation = (1.0/12.0) + ((x + dx/2 - leftx)/dx) - 2*(cube(dx/2) - cube(leftx - x))/(3*cube(dx));
			} else if (x > rightx + dx/2) {
				correlation = 2*cube(rightx - (x - dx))/(3*cube(dx));
			} else if(x > rightx){
				correlation = (1.0/12.0) + ((-(x - dx/2) + rightx)/dx) - 2*(cube(dx/2) - cube(x - rightx))/(3*cube(dx));
			} else if (x < leftx + dx/2) {
				correlation = 0.5 + ((x - leftx)/dx) - 2*(cube(x - leftx))/(3*cube(dx));
			} else if(x < leftx + dx){
				correlation = 11.0/12.0 + 2*(cube(dx/2) - cube(leftx - (x - dx)))/(3*cube(dx));
			} else if (x > rightx - dx/2) {
				correlation = 0.5 + ((rightx - x)/dx) - 2*(cube(rightx - x))/(3*cube(dx));
			} else if(x > rightx - dx) {
				correlation = 11.0/12.0 + 2*(cube(dx/2) - cube(x + dx - rightx))/(3*cube(dx));
			}else {
				correlation = 1;
			}
			break;
		default:
			errorLogFile = fopen("./output/errorLog.dat", "w");
			fprintf(errorLogFile, "spline order is wrong");
			fclose(errorLogFile);
			exit(0);
	}

	return correlation;
}