#include <cmath>

#include "util.h"
#include "simulation.h"
#include "constants.h"

void Simulation::collectParticlesIntoBins() {

	for (int i = 0; i < xnumber; ++i) {
				particlesInBbin[i].clear();
	}

	for (int i = 0; i < xnumber + 1; ++i) {
				particlesInEbin[i].clear();
	}
	double fullSum = 0;
	int pcount = 0;

	#pragma omp parallel for private(pcount) 
	for (pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		checkParticleInBox(*particle);

		int xcount = floor((particle->x - xgrid[0]) / deltaX);

		double correlationSum = 0;
		for (int i = xcount - 1; i <= xcount + 1; ++i) {
			if (particleCrossBbin(*particle, i)) {
				pushParticleIntoBbin(particle, i);
			}
		}

		xcount = floor(((particle->x - xgrid[0]) / deltaX) + 0.5);

		for (int i = xcount - 1; i <= xcount + 1; ++i) {
			int tempi = i;
			if( tempi == -1) {
				tempi = xnumber - 1;
			}
			if(tempi == xnumber + 1) {
				tempi = 1;
			}
			if (particleCrossEbin(*particle, tempi)) {
				pushParticleIntoEbin(particle, tempi);
			}
		}
	}
}

void Simulation::pushParticleIntoEbin(Particle* particle, int i) {
	if (i < 0) return;
	if (i > xnumber) return;

	particlesInEbin[i].push_back(particle);
}

void Simulation::pushParticleIntoBbin(Particle* particle, int i) {
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

	particlesInBbin[i].push_back(particle);
}

bool Simulation::particleCrossBbin(Particle& particle, int i) {
	if(boundaryConditionType == PERIODIC){
		if(i < 0) {
			i = xnumber - 1;
		} else if(i >= xnumber) {
			i = 0;
		}	

		if (i == 0) {
			if ((xgrid[i + 1] < particle.x - particle.dx) && (xgrid[xnumber] > particle.x + particle.dx))
				return false;
		} else if (i == xnumber - 1) {
			if ((xgrid[i] > particle.x + particle.dx) && (xgrid[0] < particle.x - particle.dx))
				return false;
		} else {
			if ((xgrid[i] > particle.x + particle.dx) || (xgrid[i + 1] < particle.x - particle.dx))
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
			if (xgrid[i + 1] < particle.x - particle.dx)
				return false;
		} else if (i == xnumber - 1) {
			if (xgrid[i] > particle.x + particle.dx)
				return false;
		} else {
			if ((xgrid[i] > particle.x + particle.dx) || (xgrid[i + 1] < particle.x - particle.dx))
				return false;
		}
		return true;

	}
	return false;
}

bool Simulation::particleCrossEbin(Particle& particle, int i) {
	if(boundaryConditionType == PERIODIC){
		if(i == 0){
			if(xgrid[0] + (deltaX/2) < particle.x - particle.dx && xgrid[xnumber] - (deltaX/2) > particle.x + particle.dx)
				return false;
		} else if(i == xnumber){
			if(xgrid[0] + (deltaX/2) < particle.x - particle.dx && xgrid[xnumber] - (deltaX/2) > particle.x + particle.dx)
				return false;
		} else {
			if ((xgrid[i] - (deltaX / 2) > particle.x + particle.dx) || (xgrid[i + 1] - (deltaX / 2) < particle.x - particle.dx))
				return false;
		}
		return true;
	} else if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
		if(i == 0){
			if(xgrid[0] + (deltaX/2) < particle.x - particle.dx)
				return false;
		} else if(i == xnumber){
			if(xgrid[xnumber] - (deltaX/2) > particle.x + particle.dx)
				return false;
		} else {
			if ((xgrid[i] - (deltaX / 2) > particle.x + particle.dx) || (xgrid[i + 1] - (deltaX / 2) < particle.x - particle.dx))
				return false;
		}
		return true;
	}

	return false;
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

	int xcount = floor((particle.x - xgrid[0])/deltaX);
	double rightWeight = (particle.x - xgrid[xcount])/deltaX;
	double leftWeight = (xgrid[xcount + 1] - particle.x)/deltaX;

	return Efield[xcount]*leftWeight + Efield[xcount + 1]*rightWeight;
}

Vector3d Simulation::correlationTempEfield(Particle& particle) {
	//checkParticleInBox(particle);

	/*int xcount = floor((particle.x / deltaX) + 0.5);

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = xcount - 1; i <= xcount + 1; ++i) {
				result = result + correlationFieldWithTempEbin(particle, i);
	}

	return result;*/
	int xcount = floor((particle.x - xgrid[0])/deltaX);
	double rightWeight = (particle.x - xgrid[xcount])/deltaX;
	double leftWeight = (xgrid[xcount + 1] - particle.x)/deltaX;

	return tempEfield[xcount]*leftWeight + tempEfield[xcount + 1]*rightWeight;
}

Vector3d Simulation::correlationNewEfield(Particle& particle) {
	//checkParticleInBox(particle);

	/*int xcount = floor((particle.x / deltaX) + 0.5);

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = xcount - 1; i <= xcount + 1; ++i) {
				result = result + correlationFieldWithNewEbin(particle, i);
	}

	return result;*/
	int xcount = floor((particle.x - xgrid[0])/deltaX);
	double rightWeight = (particle.x - xgrid[xcount])/deltaX;
	double leftWeight = (xgrid[xcount + 1] - particle.x)/deltaX;

	return newEfield[xcount]*leftWeight + newEfield[xcount + 1]*rightWeight;
}

Vector3d Simulation::correlationBfield(Particle& particle) {

	double x = particle.x;
	if(x > xgrid[xnumber]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}
	if(x < xgrid[0]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}

	int xcount = floor(((x - xgrid[0]) / deltaX) + 0.5);
	Vector3d leftField = Vector3d(0, 0, 0);
	Vector3d rightField = Vector3d(0, 0, 0);

	double rightWeight;
	double leftWeight;

	if(xcount > 0 && xcount < xnumber){
		rightWeight = (particle.x - middleXgrid[xcount-1])/deltaX;
		leftWeight = (middleXgrid[xcount] - particle.x)/deltaX;
		leftField = Bfield[xcount-1];
		rightField = Bfield[xcount];
	} else if(xcount == 0){
		leftWeight = (middleXgrid[xcount] - particle.x)/deltaX;
		rightWeight = 1.0 - leftWeight;
		rightField = Bfield[xcount];
		if(boundaryConditionType == PERIODIC){
			leftField = Bfield[xnumber - 1];
		} else {
			leftField = Bfield[0];
		}
	} else if(xcount == xnumber){
		rightWeight = (particle.x - middleXgrid[xcount-1])/deltaX;
		leftWeight = 1 - rightWeight;
		leftField = Bfield[xcount - 1];
		if(boundaryConditionType == PERIODIC){
			rightField = Bfield[0];
		} else {
			rightField = B0;
		}
	}

	return leftField*leftWeight + rightField*rightWeight;

	/*double leftX;
	double rightX;
	Vector3d result;
	Vector3d leftField;
	Vector3d rightField;
	if( xcount == 0) {
		leftX = middleXgrid[0] - deltaX;
		rightX = middleXgrid[0];
		leftField = Bfield[xnumber - 1];
		rightField = Bfield[0];
	} else if(xcount == xnumber) {
		leftX = middleXgrid[xnumber - 1];
		rightX = middleXgrid[xnumber - 1] + deltaX;
		leftField = Bfield[xnumber - 1];
		rightField = Bfield[0];
	} else {
		leftX = middleXgrid[xcount - 1];
		rightX = middleXgrid[xcount];
		leftField = Bfield[xcount-1];
		rightField = Bfield[xcount];
	}

	double rightWeight = (x - leftX)/deltaX;
	double leftWeight = (rightX - x)/deltaX;

	if(rightWeight > 1 || rightWeight < 0){
		printf("rightWeight is invalid\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "rightWeight = %15.10g is invalid\n", rightWeight);
		fclose(errorLogFile);
		exit(0);
	}
	if(leftWeight > 1 || leftWeight < 0){
		printf("leftWeight is invalid\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "leftWeight = %15.10g is invalid\n", leftWeight);
		fclose(errorLogFile);
		exit(0);
	}

	result = rightField*rightWeight + leftField*leftWeight;

	return result;*/
}

Vector3d Simulation::correlationNewBfield(Particle& particle) {

	double x = particle.x;

	int xcount = floor(((x - xgrid[0]) / deltaX) + 0.5);

	if(x > xgrid[xnumber]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}
	if(x < xgrid[0]) {
		printf("particle out of box in correlationBfield\n");
		exit(0);
	}

	Vector3d leftField;
	Vector3d rightField;

	double rightWeight;
	double leftWeight;

	if(xcount > 0 && xcount < xnumber){
		rightWeight = (particle.x - middleXgrid[xcount-1])/deltaX;
		leftWeight = (middleXgrid[xcount] - particle.x)/deltaX;
		leftField = newBfield[xcount-1];
		rightField = newBfield[xcount];
	} else if(xcount == 0){
		leftWeight = (middleXgrid[xcount] - particle.x)/deltaX;
		rightWeight = 1.0 - leftWeight;
		rightField = newBfield[xcount];
		if(boundaryConditionType == PERIODIC){
			leftField = newBfield[xnumber - 1];
		} else {
			leftField = newBfield[0];
		}
	} else if(xcount == xnumber){
		rightWeight = (particle.x - middleXgrid[xcount-1])/deltaX;
		leftWeight = 1 - rightWeight;
		leftField = newBfield[xcount - 1];
		if(boundaryConditionType == PERIODIC){
			rightField = newBfield[0];
		} else {
			rightField = B0;
		}
	}

	return leftField*leftWeight + rightField*rightWeight;
	/*Vector3d result = Vector3d(0, 0, 0);

	double leftX;
	double rightX;
	Vector3d leftField;
	Vector3d rightField;
	if( xcount == 0) {
		leftX = middleXgrid[0] - deltaX;
		rightX = middleXgrid[0];
		leftField = newBfield[xnumber - 1];
		rightField = newBfield[0];
	} else if(xcount == xnumber) {
		leftX = middleXgrid[xnumber - 1];
		rightX = middleXgrid[xnumber - 1] + deltaX;
		leftField = newBfield[xnumber - 1];
		rightField = newBfield[0];
	} else {
		leftX = middleXgrid[xcount - 1];
		rightX = middleXgrid[xcount];
		leftField = newBfield[xcount-1];
		rightField = newBfield[xcount];
	}

	double rightWeight = (x - leftX)/deltaX;
	double leftWeight = (rightX - x)/deltaX;

	if(rightWeight > 1 || rightWeight < 0){
		printf("rightWeight is invalid\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "rightWeight = %15.10g is invalid\n", rightWeight);
		fclose(errorLogFile);
		exit(0);
	}
	if(leftWeight > 1 || leftWeight < 0){
		printf("leftWeight is invalid\n");
		errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "rleftWeight = %15.10g is invalid\n", leftWeight);
		fclose(errorLogFile);
		exit(0);
	}

	result = rightField*rightWeight + leftField*leftWeight;

	return result;*/
}

double Simulation::correlationWithBbin(Particle& particle, int i) {
	if (! particleCrossBbin(particle, i))
		return 0.0;

	double x = particle.x;

	double leftx;
	double rightx;

	if(i == -1) {
		if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
			return 0.0;
		}
		if (particle.x - particle.dx > xgrid[0]) {
			leftx = xgrid[xnumber-1];
			rightx = xgrid[xnumber];
		} else {
			leftx = xgrid[0] - deltaX;
			rightx = xgrid[0];
		}	
	} else if(i == 0){
		if(boundaryConditionType == PERIODIC){
			if (particle.x - particle.dx > xgrid[1]) {
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
			if (particle.x + particle.dx < xgrid[xnumber - 1]) {
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
		if (particle.x + particle.dx < xgrid[xnumber]) {
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



	double correlationx = correlationBspline(x, particle.dx, leftx, rightx);

	return correlationx;
}

double Simulation::correlationWithEbin(Particle& particle, int i) {
	double x = particle.x;

	double leftx;
	double rightx;

	if(i == -1) {
		if(boundaryConditionType == SUPER_CONDUCTOR_LEFT){
			leftx = xgrid[0] - deltaX;
			rightx = xgrid[0];
		}
		if (particle.x - particle.dx > xgrid[0] - deltaX/2) {
			leftx = xgrid[xnumber - 1] - deltaX/2;
			rightx = xgrid[xnumber - 1] + deltaX/2;
		} else {
			leftx = xgrid[0] - 3*deltaX/2;
			rightx = xgrid[0] - deltaX/2;
		}		
	} else if(i == 0){
		if(boundaryConditionType == PERIODIC){
			if (particle.x - particle.dx > xgrid[0] + deltaX/2) {
				leftx = xgrid[xnumber] - deltaX/2;
				rightx = xgrid[xnumber] + deltaX/2;
			} else {
				leftx = xgrid[0] - deltaX/2;
				rightx = xgrid[0] + deltaX/2;
			}
		} else if(boundaryConditionType == FREE_BOTH){
			leftx = xgrid[0] - deltaX/2;
			rightx = xgrid[0] + deltaX/2;
		} else {
			leftx = xgrid[0];
			rightx = xgrid[0] + deltaX/2;
		}
	} else if (i == xnumber) {
		if(boundaryConditionType == PERIODIC){
			if (particle.x + particle.dx < xgrid[xnumber] - deltaX/2) {
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
		if (particle.x + particle.dx < xgrid[xnumber] + deltaX/2) {
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

	double correlationx = correlationBspline(x, particle.dx, leftx, rightx);

	return correlationx;
}

double Simulation::correlationBspline(const double& x, const double& dx, const double& leftx, const double& rightx) {
	double dx2 = dx*dx;
	double dx3 = dx*dx*dx;
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
			if ( x < leftx + dx){
				correlation = 0.5*(x + dx - leftx)/dx;
			} else if( x > rightx - dx){
				correlation = 0.5*(rightx - (x - dx))/dx;
			} else {
				correlation = 1;
			}
			break;
		case 1:
			if( x < leftx){
				correlation = 0.5*sqr(x + dx - leftx)/dx2;
			} else if (x < leftx + deltaX){
				correlation = 1 - 0.5*sqr(leftx - (x - dx))/dx2;
			} else if (x > rightx){
				correlation = 0.5*sqr(rightx - (x - dx))/dx2;
			} else if (x > rightx - deltaX){
				correlation = 1 - 0.5*sqr(x + dx - rightx)/dx2;
			} else {
				correlation = 1;
			}
			break;
		case 2:
			if (x < leftx - dx/2) {
				correlation = 2*cube(x + dx - leftx)/(3*dx3);
			} else if(x < leftx){
				correlation = (1.0/12.0) + ((x + dx/2 - leftx)/dx) - 2*(cube(dx/2) - cube(leftx - x))/(3*dx3);
			} else if (x > rightx + dx/2) {
				correlation = 2*cube(rightx - (x - dx))/(3*dx3);
			} else if(x > rightx){
				correlation = (1.0/12.0) + ((-(x - dx/2) + rightx)/dx) - 2*(cube(dx/2) - cube(x - rightx))/(3*dx3);
			} else if (x < leftx + dx/2) {
				correlation = 0.5 + ((x - leftx)/dx) - 2*(cube(x - leftx))/(3*dx3);
			} else if(x < leftx + dx){
				correlation = 11.0/12.0 + 2*(cube(dx/2) - cube(leftx - (x - dx)))/(3*dx3);
			} else if (x > rightx - dx/2) {
				correlation = 0.5 + ((rightx - x)/dx) - 2*(cube(rightx - x))/(3*dx3);
			} else if(x > rightx - dx) {
				correlation = 11.0/12.0 + 2*(cube(dx/2) - cube(x + dx - rightx))/(3*dx3);
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