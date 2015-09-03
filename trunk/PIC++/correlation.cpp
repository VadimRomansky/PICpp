#include <cmath>

#include "util.h"
#include "simulation.h"

void Simulation::collectParticlesIntoBins() {
	double*** correlations = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		correlations[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j){
			correlations[i][j] = new double[znumber];
			for(int k = 0; k < znumber; ++k){
				correlations[i][j][k] = 0;
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				particlesInBbin[i][j][k].clear();
			}
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				particlesInEbin[i][j][k].clear();
			}
		}
	}

	//FILE* debugFile = fopen("./output/particleCorrelations1.dat","w");
	double fullSum = 0;
	#pragma omp parallel for
	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		checkParticleInBox(*particle);

		int xcount = floor(particle->coordinates.x / deltaX);
		int ycount = floor(particle->coordinates.y / deltaY);
		int zcount = floor(particle->coordinates.z / deltaZ);

		double correlationSum = 0;
		int counter = 0;
		for (int i = xcount - 1; i <= xcount + 1; ++i) {
			for (int j = ycount - 1; j <= ycount + 1; ++j) {
				for (int k = zcount - 1; k <= zcount + 1; ++k) {
					if (particleCrossBbin(*particle, i, j, k)) {
						double correlation = correlationWithBbin(*particle, i, j, k);
						pushParticleIntoBbin(particle, i, j, k);
						int tempi = i;
						if(i < 0){
							tempi = xnumber - 1;
						}
						if(i >= xnumber){
							tempi = 0;
						}
						int tempj = j;
						if(j < 0){
							tempj = ynumber - 1;
						}
						if(j >= ynumber){
							tempj = 0;
						}
						int tempk = k;
						if(k < 0){
							tempk = znumber - 1;
						}
						if(k >= znumber){
							tempk = 0;
						}
						//fprintf(debugFile, "particle %d i %d j %d k %d correlation %15.10g\n", particle->number, tempi, tempj, tempk, correlation);
						correlations[tempi][tempj][tempk] += correlation*particle->charge*particle->weight/(deltaX*deltaY*deltaZ);
						/*if((i == xcount - 1) || (j = ycount - 1) || (k == zcount - 1)) {
							printf("aaa\n");
						}
						bool f = particleCrossBbin(*particle, i, j, k);*/
						correlationSum += correlation;
						counter++;
						if(counter > 8){
							printf("aaaa\n");
						}
					}
				}
			}
		}

		if(abs(correlationSum - 1) > 0.000000001){
			printf("aaa\n");
		}

		fullSum += correlationSum*particle->weight*particle->charge/(xsize*ysize*zsize);

		xcount = floor((particle->coordinates.x / deltaX) + 0.5);
		ycount = floor((particle->coordinates.y / deltaY) + 0.5);
		zcount = floor((particle->coordinates.z / deltaZ) + 0.5);

		for (int i = xcount - 1; i <= xcount + 1; ++i) {
			for (int j = ycount - 1; j <= ycount + 1; ++j) {
				for (int k = zcount - 1; k <= zcount + 1; ++k) {
					if (particleCrossEbin(*particle, i, j, k)) {
						pushParticleIntoEbin(particle, i, j, k);
					}
				}
			}
		}
	}

	double fullSum2 = 0;
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				//fprintf(debugFile, "charge %15.10g\n", correlations[i][j][k]);
				fullSum2 += correlations[i][j][k]*deltaX*deltaY*deltaZ;
			}
		}
	}
	fullSum2 /= (xsize*ysize*zsize);

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			delete[] correlations[i][j];
		}
		delete[] correlations[i];
	}
	delete[] correlations;
	//fclose(debugFile);
}

void Simulation::pushParticleIntoEbin(Particle* particle, int i, int j, int k) {
	if (i < 0) return;
	if (i > xnumber) return;
	if (j < 0) {
		j = ynumber - 1;
	} else if (j > ynumber) {
		j = 1;
	}
	if (k < 0) {
		k = znumber - 1;
	} else if (k > znumber) {
		k = 1;
	}
	particlesInEbin[i][j][k].push_back(particle);
}

void Simulation::pushParticleIntoBbin(Particle* particle, int i, int j, int k) {
	if (i < 0){
		if(boundaryConditionType == SUPERCONDUCTERLEFT){
			return;
		}
		if(boundaryConditionType == PERIODIC){
			i = xnumber - 1;
		}
	}
	if (i >= xnumber){
		if(boundaryConditionType == SUPERCONDUCTERLEFT){
			return;
		}
		if(boundaryConditionType == PERIODIC){
			i = 0;
		}
	}
	if (j < 0) {
		j = ynumber - 1;
	} else if (j >= ynumber) {
		j = 0;
	}
	if (k < 0) {
		k = znumber - 1;
	} else if (k >= znumber) {
		k = 0;
	}
	particlesInBbin[i][j][k].push_back(particle);
}

bool Simulation::particleCrossBbin(Particle& particle, int i, int j, int k) {
	if(j < 0) {
		j = ynumber - 1;
	} else if(j >= ynumber) {
		j = 0;
	}

	if(k < 0) {
		k = znumber - 1;
	} else if(k >= znumber) {
		k = 0;
	}

	if(boundaryConditionType == PERIODIC){
		if(i < 0) {
				i = xnumber - 1;
		} else if(i >= xnumber) {
			i = 0;
		}
	}
	
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		if(i < 0) {
			if(particle.coordinates.x - particle.dx > 0)
				return false;
		} else {
			if ((xgrid[i] > particle.coordinates.x + particle.dx) || (xgrid[i + 1] < particle.coordinates.x - particle.dx))
				return false;
		}
	}
	if(boundaryConditionType == PERIODIC){
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

	if (k == 0) {
		if ((zgrid[k + 1] < particle.coordinates.z - particle.dz) && (zgrid[znumber] > particle.coordinates.z + particle.dz))
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

bool Simulation::particleCrossEbin(Particle& particle, int i, int j, int k) {
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		if ((xgrid[i] - (deltaX / 2) > particle.coordinates.x + particle.dx) || (xgrid[i] + (deltaX / 2) < particle.coordinates.x - particle.dx))
			return false;
	}
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
	}
	if ((ygrid[j] - (deltaY / 2) > particle.coordinates.y + particle.dy) || (ygrid[j + 1] - (deltaY / 2) < particle.coordinates.y - particle.dy))
		return false;
	if ((zgrid[k] - (deltaZ / 2) > particle.coordinates.z + particle.dz) || (zgrid[k + 1] - (deltaZ / 2) < particle.coordinates.z - particle.dz))
		return false;
	return true;
}


Vector3d Simulation::correlationTempEfield(Particle* particle) {
	return correlationTempEfield(*particle);
}

Vector3d Simulation::correlationBfield(Particle* particle) {
	return correlationBfield(*particle);
}

Vector3d Simulation::correlationEfield(Particle* particle) {
	return correlationEfield(*particle);
}

Vector3d Simulation::correlationTempEfield(Particle& particle) {
	//checkParticleInBox(particle);

	int xcount = floor((particle.coordinates.x / deltaX) + 0.5);
	int ycount = floor((particle.coordinates.y / deltaY) + 0.5);
	int zcount = floor((particle.coordinates.z / deltaZ) + 0.5);

	/*if (xcount < 0) {
		printf("xcount < 0\n");
		//exit(0);
	}

	if (xcount > xnumber) {
		printf("xcount > xnumber\n");
		//exit(0);
	}

	if (ycount < 0) {
		printf("ycount < 0\n");
		//exit(0);
	}

	if (ycount > ynumber) {
		printf("ycount > ynumber\n");
		//exit(0);
	}

	if (zcount < 0) {
		printf("zcount < 0\n");
		//exit(0);
	}

	if (zcount > znumber) {
		printf("zcount > znumber\n");
		//exit(0);
	}*/

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = xcount - 1; i <= xcount + 1; ++i) {
		for (int j = ycount - 1; j <= ycount + 1; ++j) {
			for (int k = zcount - 1; k <= zcount + 1; ++k) {
				result = result + correlationFieldWithTempEbin(particle, i, j, k);
			}
		}
	}

	return result;
}

Vector3d Simulation::correlationBfield(Particle& particle) {
	//checkParticleInBox(particle);

	int xcount = floor(particle.coordinates.x / deltaX);
	int ycount = floor(particle.coordinates.y / deltaY);
	int zcount = floor(particle.coordinates.z / deltaZ);

	/*if (xcount < 0) {
		printf("xcount < 0\n");
		//exit(0);
	}

	if (xcount >= xnumber) {
		printf("xcount >= xnumber\n");
		//exit(0);
	}

	if (ycount < 0) {
		printf("ycount < 0\n");
		//exit(0);
	}

	if (ycount >= ynumber) {
		printf("ycount >= ynumber\n");
		//exit(0);
	}

	if (zcount < 0) {
		printf("zcount < 0\n");
		//exit(0);
	}

	if (zcount >= znumber) {
		printf("zcount >= znumber\n");
		//exit(0);
	}*/

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = xcount - 1; i <= xcount + 1; ++i) {
		for (int j = ycount - 1; j <= ycount + 1; ++j) {
			for (int k = zcount - 1; k <= zcount + 1; ++k) {
				result = result + correlationFieldWithBbin(particle, i, j, k);
			}
		}
	}

	return result;
}

Vector3d Simulation::correlationEfield(Particle& particle) {
	//checkParticleInBox(particle);

	int xcount = floor((particle.coordinates.x / deltaX) + 0.5);
	int ycount = floor((particle.coordinates.y / deltaY) + 0.5);
	int zcount = floor((particle.coordinates.z / deltaZ) + 0.5);

	/*if (xcount < 0) {
		printf("xcount < 0\n");
		//exit(0);
	}

	if (xcount > xnumber) {
		printf("xcount > xnumber\n");
		//exit(0);
	}

	if (ycount < 0) {
		printf("ycount < 0\n");
		//exit(0);
	}

	if (ycount > ynumber) {
		printf("ycount > ynumber\n");
		//exit(0);
	}

	if (zcount < 0) {
		printf("zcount < 0\n");
		//exit(0);
	}

	if (zcount > znumber) {
		printf("zcount > znumber\n");
		//exit(0);
	}*/

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = xcount - 1; i <= xcount + 1; ++i) {
		for (int j = ycount - 1; j <= ycount + 1; ++j) {
			for (int k = zcount - 1; k <= zcount + 1; ++k) {
				result = result + correlationFieldWithEbin(particle, i, j, k);
			}
		}
	}

	return result;
}


Vector3d Simulation::correlationFieldWithBbin(Particle& particle, int i, int j, int k) {

	Vector3d _field;

	/*if (i < 0) {
		//_field = Vector3d(0,0,0);
		//note: not zero because reflection
		if (j < 0) {
			if (k < 0) {
				_field = Bfield[0][ynumber - 1][znumber - 1];
			} else if (k >= znumber) {
				_field = Bfield[0][ynumber - 1][0];
			} else {
				_field = Bfield[0][ynumber - 1][k];
			}
		} else if (j >= ynumber) {
			if (k < 0) {
				_field = Bfield[0][0][znumber - 1];
			} else if (k >= znumber) {
				_field = Bfield[0][0][0];
			} else {
				_field = Bfield[0][0][k];
			}
		} else {
			if (k < 0) {
				_field = Bfield[0][j][znumber - 1];
			} else if (k >= znumber) {
				_field = Bfield[0][j][0];
			} else {
				_field = Bfield[0][j][k];
			}
		}
	} else*/

	double correlation = correlationWithBbin(particle, i, j, k);

	if(i < 0) {
		i = 0;
	}

	if(j < 0) {
		j = j + ynumber;
	}

	if(j >= ynumber) {
		j = j - ynumber;
	}

	if(k < 0) {
		k  = k + znumber;
	}

	if(k >= znumber) {
		k = k - znumber;
	}

	if (i >= xnumber) {
		_field = B0;
	} else {
		_field = Bfield[i][j][k];
	}

	return _field * correlation;
}

Vector3d Simulation::correlationFieldWithEbin(Particle& particle, int i, int j, int k) {

	Vector3d _field = getEfield(i, j, k);

	double correlation = correlationWithEbin(particle, i, j, k);

	return _field * correlation;
}

Vector3d Simulation::correlationFieldWithTempEbin(Particle& particle, int i, int j, int k) {

	Vector3d _field = getTempEfield(i, j, k);

	double correlation = correlationWithEbin(particle, i, j, k);

	return _field * correlation;
}

double Simulation::correlationWithBbin(Particle& particle, int i, int j, int k) {
	if (! particleCrossBbin(particle, i, j, k))
		return 0.0;

	double x = particle.coordinates.x;
	double y = particle.coordinates.y;
	double z = particle.coordinates.z;

	double leftx;
	double rightx;
	double lefty;
	double righty;
	double leftz;
	double rightz;

	if (i < 0) {
		leftx = particle.coordinates.x - deltaX;
		rightx = xgrid[0];
	} else if (i >= xnumber) {
		leftx = xgrid[xnumber];
		rightx = particle.coordinates.x + deltaX;
	} else if(i == 0 && boundaryConditionType == PERIODIC){
		if (particle.coordinates.x - particle.dx > xgrid[1]) {
			leftx = xgrid[xnumber];
			rightx = xgrid[xnumber] + deltaX;
		} else {
			leftx = xgrid[0];
			rightx = xgrid[1];
		}
	} else if(i == xnumber - 1 && boundaryConditionType == PERIODIC){
		if (particle.coordinates.x + particle.dx < xgrid[xnumber - 1]) {
			leftx = xgrid[0] - deltaX;
			rightx = xgrid[0];
		} else {
			leftx = xgrid[xnumber - 1];
			rightx = xgrid[xnumber];
		}
	} else {
		leftx = xgrid[i];
		rightx = xgrid[i + 1];

	}


	if (j < 0) {
		lefty = ygrid[0] - deltaY;
		righty = ygrid[0];
	} else if (j >= ynumber) {
		lefty = ygrid[ynumber];
		righty = ygrid[ynumber] + deltaY;
	} else if (j == 0) {
		if (particle.coordinates.y - particle.dy > ygrid[1]) {
			lefty = ygrid[ynumber];
			righty = ygrid[ynumber] + deltaY;
		} else {
			lefty = ygrid[0];
			righty = ygrid[1];
		}
	} else if (j == ynumber - 1) {
		if (particle.coordinates.y + particle.dy < ygrid[ynumber - 1]) {
			lefty = ygrid[0] - deltaY;
			righty = ygrid[0];
		} else {
			lefty = ygrid[ynumber - 1];
			righty = ygrid[ynumber];
		}
	} else {
		lefty = ygrid[j];
		righty = ygrid[j + 1];
	}

	if (k < 0) {
		leftz = zgrid[0] - deltaZ;
		rightz = zgrid[0];
	} else if (k >= znumber) {
		leftz = zgrid[znumber];
		rightz = zgrid[znumber] + deltaZ;
	} else if (k == 0) {
		if (particle.coordinates.z - particle.dz > zgrid[1]) {
			leftz = zgrid[ynumber];
			rightz = zgrid[ynumber] + deltaZ;
		} else {
			leftz = zgrid[0];
			rightz = zgrid[1];
		}
	} else if (k == znumber - 1) {
		if (particle.coordinates.z + particle.dz < zgrid[znumber - 1]) {
			leftz = zgrid[0] - deltaZ;
			rightz = zgrid[0];
		} else {
			leftz = zgrid[znumber - 1];
			rightz = zgrid[znumber];
		}
	} else {
		leftz = zgrid[k];
		rightz = zgrid[k + 1];
	}


	double correlation = 1;

	double correlationx = correlationBspline(x, particle.dx, leftx, rightx);
	double correlationy = correlationBspline(y, particle.dy, lefty, righty);
	double correlationz = correlationBspline(z, particle.dz, leftz, rightz);

	correlation = correlationx * correlationy * correlationz;

	return correlation;
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

	if (i < 0) {
		leftx = particle.coordinates.x - 2 * deltaX;
		rightx = xgrid[0];
	} else if(i == 0){
		leftx = xgrid[0];
		rightx = xgrid[0] + deltaX/2;
	} else if (i > xnumber) {
		leftx = xgrid[xnumber] + (deltaX / 2);
		rightx = particle.coordinates.x + 2 * deltaX;
		/*} else if(i == 0){
		////note: needs because in the middle of 0 Ebin is wall
		leftx = 0;
		rightx = deltaX/2;*/
	} else {
		leftx = xgrid[i] - (deltaX / 2);
		rightx = xgrid[i] + (deltaX / 2);
	}


	if (j < 0) {
		lefty = particle.coordinates.y - (2 * deltaY);
		righty = ygrid[0] - (deltaY / 2);
	} else if (j > ynumber) {
		lefty = ygrid[ynumber] + (deltaY / 2);
		righty = particle.coordinates.y + (2 * deltaY);
	} else {
		lefty = ygrid[j] - (deltaY / 2);
		righty = ygrid[j] + (deltaY / 2);
	}

	if (k < 0) {
		leftz = particle.coordinates.z - 2 * deltaZ;
		rightz = zgrid[0] - (deltaZ/2);
	} else if (k > znumber) {
		leftz = zgrid[znumber] + (deltaZ / 2);
		rightz = particle.coordinates.z + 2 * deltaZ;
	} else {
		leftz = zgrid[k] - (deltaZ / 2);
		rightz = zgrid[k] + (deltaZ / 2);
	}


	double correlation = 1;

	double correlationx = correlationBspline(x, particle.dx, leftx, rightx);
	double correlationy = correlationBspline(y, particle.dy, lefty, righty);
	double correlationz = correlationBspline(z, particle.dz, leftz, rightz);

	correlation = correlationx * correlationy * correlationz;

	return correlation;
}

double Simulation::correlationBspline(const double& x, const double& dx, const double& leftx, const double& rightx) {

	if (rightx < leftx) {
		printf("rightx < leftx\n");
		exit(0);
	}
	if (dx > rightx - leftx) {
		printf("dx < rightx - leftx\n");
		exit(0);
	}

	double correlation = 0;

	if (x < leftx - dx)
		return 0;
	if (x > rightx + dx)
		return 0;

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

	//correlation /= dx * dx;

	return correlation;
}