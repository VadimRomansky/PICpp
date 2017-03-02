#include <cmath>
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "util.h"
#include "constants.h"
#include "particle.h"
#include "vector3d.h"
#include "simulation.h"


Vector3d Simulation::correlationTempEfield(Particle* particle) {
	return correlationTempEfield(*particle);
}

Vector3d Simulation::correlationNewEfield(Particle* particle) {
	return correlationNewEfield(*particle);
}

Vector3d Simulation::correlationBfield(Particle* particle) const {
	return correlationBfield(*particle);
}

Vector3d Simulation::correlationEfield(Particle* particle) {
	return correlationEfield(*particle);
}

Vector3d Simulation::correlationGeneralEfield(Particle& particle, Vector3d*** field, Vector3d*** additionalFieldLeft, Vector3d*** additionalFieldRight) {
	//alertNaNOrInfinity(particle.coordinates.x, "particle.x = NaN in correlationGeneralEfield\n");
	int xcount = floor(((particle.coordinates.x - xgrid[0]) / deltaX) + 0.5);
	int ycount = floor(((particle.coordinates.y - ygrid[0]) / deltaY) + 0.5);
	int zcount = floor(((particle.coordinates.z - zgrid[0]) / deltaZ) + 0.5);
	if (xcount < 0) {
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		printf("xcount < 0 in correlationGeneralEfield\n");
		fprintf(errorLogFile, "xcount < 0 in correlationGeneralEfield\n");
		printf("xgrid[0] = %g particle.x = %g number = %d\n", xgrid[0], particle.coordinates.x, particle.number);
		fprintf(errorLogFile, "xgrid[0] = %g particle.x = %g number = %d\n", xgrid[0], particle.coordinates.x, particle.number);
		exit(0);
	}
	if (xcount > xnumber + 1) {
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		printf("xcount > xnumber + 1 in correlationGeneralEfield\n");
		fprintf(errorLogFile, "xcount > xnumber + 1 in correlationGeneralEfield\n");
		printf("xgrid[xnumber + 1] = %g particle.x = %g\n", xgrid[xnumber + 1], particle.coordinates.x);
		fprintf(errorLogFile, "xgrid[xnumber + 1] = %g particle.x = %g\n", xgrid[xnumber + 1], particle.coordinates.x);
		exit(0);
	}

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = 0; i < splineOrder + 2; ++i) {
		for (int j = 0; j < min2(splineOrder + 2, ynumber); ++j) {
			for (int k = 0; k < min2(splineOrder + 2, znumber); ++k) {
				int curI = particle.correlationMapNode.xindex[i];
				int curJ = particle.correlationMapNode.yindex[j];
				int curK = particle.correlationMapNode.zindex[k];

				while (curJ > ynumber) {
					curJ = curJ - ynumber;
				}
				while (curJ < 0) {
					curJ = curJ + ynumber;
				}
				while (curK > znumber) {
					curK = curK - znumber;
				}
				while (curK < 0) {
					curK = curK + znumber;
				}
				Vector3d E = Vector3d(0, 0, 0);
				if (curI < 0) {
					E = additionalFieldLeft[-curI - 1][curJ][curK];
				} else if (curI > xnumber + 1) {
					E = additionalFieldRight[curI - xnumber - 2][curJ][curK];
				} else {
					E = field[curI][curJ][curK];
				}
				double correlation = particle.correlationMapNode.xcorrelation[i] * particle.correlationMapNode.ycorrelation[j] * particle.correlationMapNode.zcorrelation[k];
				result = result + E * correlation;
			}
		}
	}

	return result;
}

Vector3d Simulation::correlationEfield(Particle& particle) {
	return correlationGeneralEfield(particle, Efield, additionalEfieldLeft, additionalEfieldRight);
}

Vector3d Simulation::correlationTempEfield(Particle& particle) {
	return correlationGeneralEfield(particle, tempEfield, additionalTempEfieldLeft, additionalTempEfieldRight);
}

Vector3d Simulation::correlationNewEfield(Particle& particle) {
	return correlationGeneralEfield(particle, newEfield, additionalNewEfieldLeft, additionalNewEfieldRight);
}

Vector3d Simulation::correlationBfield(Particle& particle) const {
	Vector3d result = Vector3d(0, 0, 0);
	for (int i = 0; i < splineOrder + 2; ++i) {
		for (int j = 0; j < min2(splineOrder + 2, ynumber); ++j) {
			for (int k = 0; k < min2(splineOrder + 2, znumber); ++k) {
				int curI = particle.correlationMapCell.xindex[i];
				int curJ = particle.correlationMapCell.yindex[j];
				int curK = particle.correlationMapCell.zindex[k];

				while (curJ >= ynumber) {
					curJ = curJ - ynumber;
				}
				while (curJ < 0) {
					curJ = curJ + ynumber;
				}
				while (curK >= znumber) {
					curK = curK - znumber;
				}
				while (curK < 0) {
					curK = curK + znumber;
				}
				Vector3d B = Vector3d(0, 0, 0);
				if (curI < 0) {
					B = additionalBfieldLeft[-curI - 1][curJ][curK];
				} else if (curI > xnumber) {
					B = additionalBfieldRight[curI - xnumber - 1][curJ][curK];
				} else {
					B = Bfield[curI][curJ][curK];
				}
				double correlation = particle.correlationMapCell.xcorrelation[i] * particle.correlationMapCell.ycorrelation[j] * particle.correlationMapCell.zcorrelation[k];
				result = result + B * correlation;
			}
		}
	}

	return result;
}

double Simulation::correlationWithBbin(Particle& particle, int i, int j, int k) {
	int tempI = -1;
	for (int index = 0; index < splineOrder + 2; ++index) {
		if (particle.correlationMapCell.xindex[index] == i) {
			tempI = index;
			break;
		}
	}
	if (tempI == -1) {
		return 0.0;
	}
	double correlationX = particle.correlationMapCell.xcorrelation[tempI];
	double correlationY = 0;
	if (ynumber == 1) {
		if (j == 0) {
			correlationY = 1.0;
		} else {
			return 0.0;
		}
	} else {
		tempI = -1;
		for (int index = 0; index < min2(splineOrder + 2, ynumber); ++index) {
			if (particle.correlationMapCell.yindex[index] == j) {
				tempI = index;
				break;
			}
		}
		if (tempI == -1) {
			return 0.0;
		}
		correlationY = particle.correlationMapCell.ycorrelation[tempI];
	}

	double correlationZ = 0;
	if (znumber == 1) {
		if (k == 0) {
			correlationZ = 1.0;
		} else {
			return 0.0;
		}
	} else {
		tempI = -1;
		for (int index = 0; index < min2(splineOrder + 2, znumber); ++index) {
			if (particle.correlationMapCell.zindex[index] == k) {
				tempI = index;
				break;
			}
		}
		if (tempI == -1) {
			return 0.0;
		}
		correlationZ = particle.correlationMapCell.zcorrelation[tempI];
	}

	return correlationX * correlationY * correlationZ;
}

double Simulation::correlationWithEbin(Particle& particle, int i, int j, int k) {
	int tempI = -1;
	for (int index = 0; index < splineOrder + 2; ++index) {
		if (particle.correlationMapNode.xindex[index] == i) {
			tempI = index;
			break;
		}
	}
	if (tempI == -1) {
		return 0.0;
	}
	double correlationX = particle.correlationMapNode.xcorrelation[tempI];
	double correlationY = 0;
	if (ynumber == 1) {
		if (j == 0) {
			correlationY = 1.0;
		} else {
			return 0.0;
		}
	} else {
		tempI = -1;
		for (int index = 0; index < min2(splineOrder + 2, ynumber); ++index) {
			if (particle.correlationMapNode.yindex[index] == j) {
				tempI = index;
				break;
			}
		}
		if (tempI == -1) {
			return 0.0;
		}
		correlationY = particle.correlationMapNode.ycorrelation[tempI];
	}

	double correlationZ = 0;
	if (znumber == 1) {
		if (k == 0) {
			correlationZ = 1.0;
		} else {
			return 0.0;
		}
	} else {
		tempI = -1;
		for (int index = 0; index < min2(splineOrder + 2, znumber); ++index) {
			if (particle.correlationMapNode.zindex[index] == k) {
				tempI = index;
				break;
			}
		}
		if (tempI == -1) {
			return 0.0;
		}
		correlationZ = particle.correlationMapNode.zcorrelation[tempI];
	}

	return correlationX * correlationY * correlationZ;
}

double Simulation::correlationBspline(const double& x, const double& dx, const double& leftx, const double& rightx) {
	double dx2 = dx * dx;
	double dx3 = dx * dx * dx;
	if (rightx < leftx) {
		printf("rightx < leftx\n");
		fflush(stdout);
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "rightx = %15.10g < leftx = %15.10g\n", rightx, leftx);
		fclose(errorLogFile);
		exit(0);
	}

	double correlation = 0;

	if (x < leftx - dx)
		return 0;
	if (x > rightx + dx)
		return 0;

	switch (splineOrder) {
	case 0:
		if (x < leftx + dx) {
			correlation = 0.5 * (x + dx - leftx) / dx;
		} else if (x > rightx - dx) {
			correlation = 0.5 * (rightx - (x - dx)) / dx;
		} else {
			correlation = 1;
		}
		if (correlation < 0 && correlation > -1E-14) {
			correlation = 0;
		}
		break;
	case 1:
		if (x < leftx) {
			correlation = 0.5 * sqr(x + dx - leftx) / dx2;
		} else if (x > rightx) {
			correlation = 0.5 * sqr(rightx - (x - dx)) / dx2;
		} else {
			correlation = 1 - 0.5 * (sqr(x + dx - rightx) + sqr(leftx - (x - dx))) / dx2;
		}
		if (correlation < 0 && correlation > -1E-14) {
			correlation = 0;
		}
		break;
	case 2:

		if (x + dx < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			correlation = (-9.0 * t3 / 16.0) + (27.0 * t2 / 16.0) + (-27.0 * t / 16.0) + 9.0 / 16.0;
		} else if (x + dx / 3 < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			correlation = (27.0 * t3 / 16.0) + (-9.0 * t2 / 16.0) + (-15.0 * t / 16.0) + 23.0 / 48.0;
		} else if (x - dx / 3 < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			correlation = (-27.0 * t3 / 16.0) + (-63.0 * t2 / 16.0) + (-33.0 * t / 16.0) + 17.0 / 48.0;
		} else if (x - dx < rightx) {
			double t = (rightx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			correlation = (9.0 * t3 / 16.0) + (27.0 * t2 / 16.0) + (27.0 * t / 16.0) + 9.0 / 16.0;
		} else {
			return 0;
		}
		if (correlation < 0 && correlation > -1E-14) {
			correlation = 0;
		}
		break;
	case 3:
		if (x + dx < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			double t4 = t * t3;
			correlation = (2.0 * t4 / 3.0) + (-8.0 * t3 / 3.0) + (4 * t2) + (-8.0 * t / 3.0) + 2.0 / 3.0;
		} else if (x + dx / 2 < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			double t4 = t * t3;
			correlation = (-8.0 * t4 / 3.0) + (4 * t3) + (-t2) + (-t) + 11.0 / 24.0;
		} else if (x < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			double t4 = t * t3;
			correlation = (4 * t4) + (4 * t3) + (-t2) + (-t) + 11.0 / 24.0;
		} else if (x - dx / 2 < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			double t4 = t * t3;
			correlation = (-8.0 * t4 / 3.0) + (-28.0 * t3 / 3.0) + (-11 * t2) + (-13.0 * t / 3) + 1.0 / 24;
		} else if (x - dx < rightx) {
			double t = (rightx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			double t4 = t * t3;
			correlation = (2.0 * t4 / 3.0) + (8.0 * t3 / 3.0) + (4 * t2) + (8.0 * t / 3.0) + 2.0 / 3.0;
		} else {
			return 0;
		}
		if (correlation < 0 && correlation > -1E-14) {
			correlation = 0;
		}
		break;
	case 4:
		if (x + dx < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			double t4 = t * t3;
			double t5 = t * t4;
			correlation = (-625.0 * t5 / 768.0) + (3125.0 * t4 / 768.0) + (-3125.0 * t3 / 384.0) + (3125.0 * t2 / 384.0) + (-3125.0 * t / 768.0) + 625.0 / 768.0;
			/*if(correlation < 0 && correlation > -1E-14){
			    correlation = 0;
			}
			if(correlation < 0) {

			    errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			    fprintf(errorLogFile, "correlation < 0");
			    printf("correlation < 0");
			    fclose(errorLogFile);
			    MPI_Finalize();
			    exit(0);

			}*/
		} else if (x + 3 * dx / 5 < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			double t4 = t * t3;
			double t5 = t * t4;
			correlation = (3125.0 * t5 / 768.0) + (-8125.0 * t4 / 768.0) + (3625.0 * t3 / 384.0) + (-925.0 * t2 / 384.0) + (-695.0 * t / 768.0) + 1667.0 / 3840;
			/*if(correlation < 0 && correlation > -1E-14){
			    correlation = 0;
			}
			if(correlation < 0) {

			    errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			    fprintf(errorLogFile, "correlation < 0");
			    printf("correlation < 0");
			    fclose(errorLogFile);
			    MPI_Finalize();
			    exit(0);

			}*/
		} else if (x + dx / 5 < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			double t4 = t * t3;
			double t5 = t * t4;
			correlation = (-3125.0 * t5 / 384.0) + (625.0 * t4 / 384.0) + (875.0 * t3 / 192.0) + (-275.0 * t2 / 192.0) + (-385.0 * t / 384.0) + 841.0 / 1920.0;
			/*if(correlation < 0 && correlation > -1E-14){
			    correlation = 0;
			}
			if(correlation < 0) {

			    errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			    fprintf(errorLogFile, "correlation < 0");
			    printf("correlation < 0");
			    fclose(errorLogFile);
			    MPI_Finalize();
			    exit(0);

			}*/
		} else if (x - dx / 5 < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			double t4 = t * t3;
			double t5 = t * t4;
			correlation = (3125.0 * t5 / 384.0) + (6875.0 * t4 / 384.0) + (2125.0 * t3 / 192.0) + (-25.0 * t2 / 192.0) - (335.0 * t / 384.0) + 851.0 / 1920.0;
			/*if(correlation < 0 && correlation > -1E-14){
			    correlation = 0;
			}
			if(correlation < 0) {

			    errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			    fprintf(errorLogFile, "correlation < 0");
			    printf("correlation < 0");
			    fclose(errorLogFile);
			    MPI_Finalize();
			    exit(0);

			}*/
		} else if (x - 3 * dx / 5 < rightx) {
			double t = (leftx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			double t4 = t * t3;
			double t5 = t * t4;
			correlation = (-3125.0 * t5 / 768.0) + (-14375.0 * t4 / 768.0) + (-12625.0 * t3 / 384.0) + (-10175.0 * t2 / 384.0) + (-6745.0 * t / 768.0) - 1943.0 / 3840.0;
			/*if(correlation < 0 && correlation > -1E-14){
			    correlation = 0;
			}
			if(correlation < 0) {

			    errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			    fprintf(errorLogFile, "correlation < 0");
			    printf("correlation < 0");
			    fclose(errorLogFile);
			    MPI_Finalize();
			    exit(0);

			}*/
		} else if (x - dx < rightx) {
			double t = (rightx - x) / dx;
			double t2 = t * t;
			double t3 = t * t2;
			double t4 = t * t3;
			double t5 = t * t4;
			correlation = (625.0 * t5 / 768.0) + (3125.0 * t4 / 768.0) + (3125.0 * t3 / 384.0) + (3125.0 * t2 / 384.0) + (3125.0 * t / 768.0) + 625.0 / 768.0;
			/*if(correlation < 0 && correlation > -1E-14){
			    correlation = 0;
			}
			if(correlation < 0) {

			    errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
			    fprintf(errorLogFile, "correlation < 0");
			    printf("correlation < 0");
			    fclose(errorLogFile);
			    MPI_Finalize();
			    exit(0);

			}*/
		} else {
			return 0;
		}
		if (correlation < 0 && correlation > -1E-14) {
			correlation = 0;
		}
		break;
	default:
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "spline order is wrong");
		fclose(errorLogFile);
		exit(0);
	}


	if (correlation > 1) {
		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "correlation > 1");
		printf("correlation > 1");
		fclose(errorLogFile);
		exit(0);
	}
	if (correlation < 0) {

		errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "correlation < 0");
		printf("correlation < 0");
		fclose(errorLogFile);
		exit(0);

	}

	return correlation;
}

void Simulation::updateParticleCorrelationMaps() {
	for (int pcount = 0; pcount < particles.size(); pcount++) {
		Particle* particle = particles[pcount];
		updateCorrelationMaps(particle);
	}
}

void Simulation::updateCorrelationMaps(Particle* particle) {
	updateCorrelationMapCell(particle);
	updateCorrelationMapNode(particle);
}

void Simulation::updateCorrelationMaps(Particle& particle) {
	updateCorrelationMapCell(particle);
	updateCorrelationMapNode(particle);
}

void Simulation::updateCorrelationMapsX(Particle& particle) {
	updateCorrelationMapCellX(particle);
	updateCorrelationMapNodeX(particle);
}

void Simulation::updateCorrelationMapsY(Particle& particle) {
	updateCorrelationMapCellY(particle);
	updateCorrelationMapNodeY(particle);
}

void Simulation::updateCorrelationMapsZ(Particle& particle) {
	updateCorrelationMapCellZ(particle);
	updateCorrelationMapNodeZ(particle);
}

void Simulation::updateCorrelationMapCell(Particle* particle) {
	updateCorrelationMapCell(*particle);
}

void Simulation::updateCorrelationMapCellX(Particle& particle) {
	int xcount = floor((particle.coordinates.x - xgrid[0]) / deltaX);
	if ((splineOrder % 2) == 0) {
		bool leftSideX = particle.coordinates.x < middleXgrid[xcount];
		int leftShiftX = leftSideX ? -1 : 0;
		int tempIndex = 0;
		for (int i = xcount + leftShiftX - (splineOrder / 2); i <= xcount + leftShiftX + (splineOrder / 2) + 1; ++i) {
			particle.correlationMapCell.xindex[tempIndex] = i;
			particle.correlationMapCell.xcorrelation[tempIndex] = correlationBspline(particle.coordinates.x, particle.dx, xgrid[0] + (i) * deltaX, xgrid[0] + (i + 1) * deltaX);
			tempIndex++;
		}
	} else {
		int tempIndex = 0;
		for (int i = xcount - (splineOrder / 2) - 1; i <= xcount + (splineOrder / 2) + 1; ++i) {
			particle.correlationMapCell.xindex[tempIndex] = i;
			particle.correlationMapCell.xcorrelation[tempIndex] = correlationBspline(particle.coordinates.x, particle.dx, xgrid[0] + (i) * deltaX, xgrid[0] + (i + 1) * deltaX);
			tempIndex++;
		}
	}
}

void Simulation::updateCorrelationMapCellY(Particle& particle) {
	int ycount = floor((particle.coordinates.y - ygrid[0]) / deltaY);
	if ((splineOrder % 2) == 0) {
		bool leftSideY = particle.coordinates.y < middleYgrid[ycount];
		int leftShiftY = leftSideY ? -1 : 0;
		int tempIndex = 0;
		if (ynumber > 1) {
			for (int j = ycount + leftShiftY - (splineOrder / 2); j <= ycount + leftShiftY + (splineOrder / 2) + 1; ++j) {
				particle.correlationMapCell.yindex[tempIndex] = j;
				particle.correlationMapCell.ycorrelation[tempIndex] = correlationBspline(particle.coordinates.y, particle.dy, ygrid[0] + (j) * deltaY, ygrid[0] + (j + 1) * deltaY);
				tempIndex++;
			}
		} else {
			for (int j = 0; j < splineOrder + 2; ++j) {
				particle.correlationMapCell.yindex[j] = j;
				particle.correlationMapCell.ycorrelation[j] = 0.0;
			}
			particle.correlationMapCell.ycorrelation[0] = 1.0;
		}
	} else {
		int tempIndex = 0;
		if (ynumber > 1) {
			for (int j = ycount - (splineOrder / 2) - 1; j <= ycount + (splineOrder / 2) + 1; ++j) {
				particle.correlationMapCell.yindex[tempIndex] = j;
				particle.correlationMapCell.ycorrelation[tempIndex] = correlationBspline(particle.coordinates.y, particle.dy, ygrid[0] + (j) * deltaY, ygrid[0] + (j + 1) * deltaY);
				tempIndex++;
			}
		} else {
			for (int j = 0; j < splineOrder + 2; ++j) {
				particle.correlationMapCell.yindex[j] = j;
				particle.correlationMapCell.ycorrelation[j] = 0.0;
			}
			particle.correlationMapCell.ycorrelation[0] = 1.0;
		}
	}
}

void Simulation::updateCorrelationMapCellZ(Particle& particle) {
	int zcount = floor((particle.coordinates.z - zgrid[0]) / deltaZ);
	if ((splineOrder % 2) == 0) {
		int tempIndex = 0;
		bool leftSideZ = particle.coordinates.z < middleZgrid[zcount];
		int leftShiftZ = leftSideZ ? -1 : 0;
		if (znumber > 1) {
			for (int k = zcount + leftShiftZ - (splineOrder / 2); k <= zcount + leftShiftZ + (splineOrder / 2) + 1; ++k) {
				particle.correlationMapCell.zindex[tempIndex] = k;
				particle.correlationMapCell.zcorrelation[tempIndex] = correlationBspline(particle.coordinates.z, particle.dz, zgrid[0] + (k) * deltaZ, zgrid[0] + (k + 1) * deltaZ);
				tempIndex++;
			}
		} else {
			for (int k = 0; k < splineOrder + 2; ++k) {
				particle.correlationMapCell.zindex[k] = k;
				particle.correlationMapCell.zcorrelation[k] = 0.0;
			}
			particle.correlationMapCell.zcorrelation[0] = 1.0;
		}
	} else {
		int tempIndex = 0;
		if (znumber > 1) {
			for (int k = zcount - (splineOrder / 2) - 1; k <= zcount + (splineOrder / 2) + 1; ++k) {
				particle.correlationMapCell.zindex[tempIndex] = k;
				particle.correlationMapCell.zcorrelation[tempIndex] = correlationBspline(particle.coordinates.z, particle.dz, zgrid[0] + (k) * deltaZ, zgrid[0] + (k + 1) * deltaZ);
				tempIndex++;
			}
		} else {
			for (int k = 0; k < splineOrder + 2; ++k) {
				particle.correlationMapCell.zindex[k] = k;
				particle.correlationMapCell.zcorrelation[k] = 0.0;
			}
			particle.correlationMapCell.zcorrelation[0] = 1.0;
		}
	}
}

void Simulation::updateCorrelationMapCell(Particle& particle) {
	int xcount = floor((particle.coordinates.x - xgrid[0]) / deltaX);
	int ycount = floor((particle.coordinates.y - ygrid[0]) / deltaY);
	int zcount = floor((particle.coordinates.z - zgrid[0]) / deltaZ);
	if ((splineOrder % 2) == 0) {
		bool leftSideX = particle.coordinates.x < middleXgrid[xcount];
		bool leftSideY = particle.coordinates.y < middleYgrid[ycount];
		bool leftSideZ = particle.coordinates.z < middleZgrid[zcount];
		int leftShiftX = leftSideX ? -1 : 0;
		int leftShiftY = leftSideY ? -1 : 0;
		int leftShiftZ = leftSideZ ? -1 : 0;
		int tempIndex = 0;
		for (int i = xcount + leftShiftX - (splineOrder / 2); i <= xcount + leftShiftX + (splineOrder / 2) + 1; ++i) {
			particle.correlationMapCell.xindex[tempIndex] = i;
			particle.correlationMapCell.xcorrelation[tempIndex] = correlationBspline(particle.coordinates.x, particle.dx, xgrid[0] + (i) * deltaX, xgrid[0] + (i + 1) * deltaX);
			tempIndex++;
		}
		tempIndex = 0;
		if (ynumber > 1) {
			for (int j = ycount + leftShiftY - (splineOrder / 2); j <= ycount + leftShiftY + (splineOrder / 2) + 1; ++j) {
				particle.correlationMapCell.yindex[tempIndex] = j;
				particle.correlationMapCell.ycorrelation[tempIndex] = correlationBspline(particle.coordinates.y, particle.dy, ygrid[0] + (j) * deltaY, ygrid[0] + (j + 1) * deltaY);
				tempIndex++;
			}
		} else {
			for (int j = 0; j < splineOrder + 2; ++j) {
				particle.correlationMapCell.yindex[j] = j;
				particle.correlationMapCell.ycorrelation[j] = 0.0;
			}
			particle.correlationMapCell.ycorrelation[0] = 1.0;
		}
		tempIndex = 0;
		if (znumber > 1) {
			for (int k = zcount + leftShiftZ - (splineOrder / 2); k <= zcount + leftShiftZ + (splineOrder / 2) + 1; ++k) {
				particle.correlationMapCell.zindex[tempIndex] = k;
				particle.correlationMapCell.zcorrelation[tempIndex] = correlationBspline(particle.coordinates.z, particle.dz, zgrid[0] + (k) * deltaZ, zgrid[0] + (k + 1) * deltaZ);
				tempIndex++;
			}
		} else {
			for (int k = 0; k < splineOrder + 2; ++k) {
				particle.correlationMapCell.zindex[k] = k;
				particle.correlationMapCell.zcorrelation[k] = 0.0;
			}
			particle.correlationMapCell.zcorrelation[0] = 1.0;
		}
	} else {
		int tempIndex = 0;
		for (int i = xcount - (splineOrder / 2) - 1; i <= xcount + (splineOrder / 2) + 1; ++i) {
			particle.correlationMapCell.xindex[tempIndex] = i;
			particle.correlationMapCell.xcorrelation[tempIndex] = correlationBspline(particle.coordinates.x, particle.dx, xgrid[0] + (i) * deltaX, xgrid[0] + (i + 1) * deltaX);
			tempIndex++;
		}
		tempIndex = 0;
		if (ynumber > 1) {
			for (int j = ycount - (splineOrder / 2) - 1; j <= ycount + (splineOrder / 2) + 1; ++j) {
				particle.correlationMapCell.yindex[tempIndex] = j;
				particle.correlationMapCell.ycorrelation[tempIndex] = correlationBspline(particle.coordinates.y, particle.dy, ygrid[0] + (j) * deltaY, ygrid[0] + (j + 1) * deltaY);
				tempIndex++;
			}
		} else {
			for (int j = 0; j < splineOrder + 2; ++j) {
				particle.correlationMapCell.yindex[j] = j;
				particle.correlationMapCell.ycorrelation[j] = 0.0;
			}
			particle.correlationMapCell.ycorrelation[0] = 1.0;
		}
		tempIndex = 0;
		if (znumber > 1) {
			for (int k = zcount - (splineOrder / 2) - 1; k <= zcount + (splineOrder / 2) + 1; ++k) {
				particle.correlationMapCell.zindex[tempIndex] = k;
				particle.correlationMapCell.zcorrelation[tempIndex] = correlationBspline(particle.coordinates.z, particle.dz, zgrid[0] + (k) * deltaZ, zgrid[0] + (k + 1) * deltaZ);
				tempIndex++;
			}
		} else {
			for (int k = 0; k < splineOrder + 2; ++k) {
				particle.correlationMapCell.zindex[k] = k;
				particle.correlationMapCell.zcorrelation[k] = 0.0;
			}
			particle.correlationMapCell.zcorrelation[0] = 1.0;
		}
	}

	double fullCorrelation = 0;
	for(int i = 0; i < splineOrder+2; ++i) {
		for(int j = 0; j < splineOrder+2; ++j) {
			for(int k = 0; k <splineOrder+2; ++k) {
				fullCorrelation += particle.correlationMapCell.xcorrelation[i]*particle.correlationMapCell.ycorrelation[j]*particle.correlationMapCell.zcorrelation[k];
			}
		}
	}
	if((fullCorrelation < 1.0 - 1E-10) || (fullCorrelation > 1.0 + 1E-10)){
		printf("full correlation cell = %20.15g\n", fullCorrelation);
	}
}

void Simulation::updateCorrelationMapNode(Particle* particle) {
	updateCorrelationMapNode(*particle);
}

void Simulation::updateCorrelationMapNodeX(Particle& particle) {
	int xcount = floor(((particle.coordinates.x - xgrid[0]) / deltaX) + 0.5);
	if ((splineOrder % 2) == 0) {
		bool leftSideX = particle.coordinates.x < xgrid[xcount];
		int leftShiftX = leftSideX ? -1 : 0;
		int tempIndex = 0;
		for (int i = xcount + leftShiftX - (splineOrder / 2); i <= xcount + leftShiftX + (splineOrder / 2) + 1; ++i) {
			particle.correlationMapNode.xindex[tempIndex] = i;
			particle.correlationMapNode.xcorrelation[tempIndex] = correlationBspline(particle.coordinates.x, particle.dx, xgrid[0] + (i - 0.5) * deltaX, xgrid[0] + (i + 0.5) * deltaX);
			tempIndex++;
		}
	} else {
		int tempIndex = 0;
		for (int i = xcount - (splineOrder / 2) - 1; i <= xcount + (splineOrder / 2) + 1; ++i) {
			particle.correlationMapNode.xindex[tempIndex] = i;
			particle.correlationMapNode.xcorrelation[tempIndex] = correlationBspline(particle.coordinates.x, particle.dx, xgrid[0] + (i - 0.5) * deltaX, xgrid[0] + (i + 0.5) * deltaX);
			tempIndex++;
		}
	}
}

void Simulation::updateCorrelationMapNodeY(Particle& particle) {
	int ycount = floor(((particle.coordinates.y - ygrid[0]) / deltaY) + 0.5);
	if ((splineOrder % 2) == 0) {
		bool leftSideY = particle.coordinates.y < ygrid[ycount];
		int leftShiftY = leftSideY ? -1 : 0;
		int tempIndex = 0;
		if (ynumber > 1) {
			for (int j = ycount + leftShiftY - (splineOrder / 2); j <= ycount + leftShiftY + (splineOrder / 2) + 1; ++j) {
				particle.correlationMapNode.yindex[tempIndex] = j;
				particle.correlationMapNode.ycorrelation[tempIndex] = correlationBspline(particle.coordinates.y, particle.dy, ygrid[0] + (j - 0.5) * deltaY, ygrid[0] + (j + 0.5) * deltaY);
				tempIndex++;
			}
		} else {
			for (int j = 0; j < splineOrder + 2; ++j) {
				particle.correlationMapNode.yindex[j] = j;
				particle.correlationMapNode.ycorrelation[j] = 0.0;
			}
			particle.correlationMapNode.ycorrelation[0] = 1.0;
		}
	} else {
		int tempIndex = 0;
		if (ynumber > 1) {
			for (int j = ycount - (splineOrder / 2) - 1; j <= ycount + (splineOrder / 2) + 1; ++j) {
				particle.correlationMapNode.yindex[tempIndex] = j;
				particle.correlationMapNode.ycorrelation[tempIndex] = correlationBspline(particle.coordinates.y, particle.dy, ygrid[0] + (j - 0.5) * deltaY, ygrid[0] + (j + 0.5) * deltaY);
				tempIndex++;
			}
		} else {
			for (int j = 0; j < splineOrder + 2; ++j) {
				particle.correlationMapNode.yindex[j] = j;
				particle.correlationMapNode.ycorrelation[j] = 0.0;
			}
			particle.correlationMapNode.ycorrelation[0] = 1.0;
		}
	}
}

void Simulation::updateCorrelationMapNodeZ(Particle& particle) {
	int zcount = floor(((particle.coordinates.z - zgrid[0]) / deltaZ) + 0.5);
	if ((splineOrder % 2) == 0) {
		bool leftSideZ = particle.coordinates.z < zgrid[zcount];
		int leftShiftZ = leftSideZ ? -1 : 0;
		int tempIndex = 0;
		if (znumber > 1) {
			for (int k = zcount + leftShiftZ - (splineOrder / 2); k <= zcount + leftShiftZ + (splineOrder / 2) + 1; ++k) {
				particle.correlationMapNode.zindex[tempIndex] = k;
				particle.correlationMapNode.zcorrelation[tempIndex] = correlationBspline(particle.coordinates.z, particle.dz, zgrid[0] + (k - 0.5) * deltaZ, zgrid[0] + (k + 0.5) * deltaZ);
				tempIndex++;
			}
		} else {
			for (int k = 0; k < splineOrder + 2; ++k) {
				particle.correlationMapNode.zindex[k] = k;
				particle.correlationMapNode.zcorrelation[k] = 0.0;
			}
			particle.correlationMapNode.zcorrelation[0] = 1.0;
		}
	} else {
		int tempIndex = 0;
		if (znumber > 1) {
			for (int k = zcount - (splineOrder / 2) - 1; k <= zcount + (splineOrder / 2) + 1; ++k) {
				particle.correlationMapNode.zindex[tempIndex] = k;
				particle.correlationMapNode.zcorrelation[tempIndex] = correlationBspline(particle.coordinates.z, particle.dz, zgrid[0] + (k - 0.5) * deltaZ, zgrid[0] + (k + 0.5) * deltaZ);
				tempIndex++;
			}
		} else {
			for (int k = 0; k < splineOrder + 2; ++k) {
				particle.correlationMapNode.zindex[k] = k;
				particle.correlationMapNode.zcorrelation[k] = 0.0;
			}
			particle.correlationMapNode.zcorrelation[0] = 1.0;
		}
	}
}

void Simulation::updateCorrelationMapNode(Particle& particle) {
	int xcount = floor(((particle.coordinates.x - xgrid[0]) / deltaX) + 0.5);
	int ycount = floor(((particle.coordinates.y - ygrid[0]) / deltaY) + 0.5);
	int zcount = floor(((particle.coordinates.z - zgrid[0]) / deltaZ) + 0.5);

	if ((splineOrder % 2) == 0) {
		bool leftSideX = particle.coordinates.x < xgrid[xcount];
		bool leftSideY = particle.coordinates.y < ygrid[ycount];
		bool leftSideZ = particle.coordinates.z < zgrid[zcount];
		int leftShiftX = leftSideX ? -1 : 0;
		int leftShiftY = leftSideY ? -1 : 0;
		int leftShiftZ = leftSideZ ? -1 : 0;
		int tempIndex = 0;
		for (int i = xcount + leftShiftX - (splineOrder / 2); i <= xcount + leftShiftX + (splineOrder / 2) + 1; ++i) {
			particle.correlationMapNode.xindex[tempIndex] = i;
			particle.correlationMapNode.xcorrelation[tempIndex] = correlationBspline(particle.coordinates.x, particle.dx, xgrid[0] + (i - 0.5) * deltaX, xgrid[0] + (i + 0.5) * deltaX);
			tempIndex++;
		}
		tempIndex = 0;
		if (ynumber > 1) {
			for (int j = ycount + leftShiftY - (splineOrder / 2); j <= ycount + leftShiftY + (splineOrder / 2) + 1; ++j) {
				particle.correlationMapNode.yindex[tempIndex] = j;
				particle.correlationMapNode.ycorrelation[tempIndex] = correlationBspline(particle.coordinates.y, particle.dy, ygrid[0] + (j - 0.5) * deltaY, ygrid[0] + (j + 0.5) * deltaY);
				tempIndex++;
			}
		} else {
			for (int j = 0; j < splineOrder + 2; ++j) {
				particle.correlationMapNode.yindex[j] = j;
				particle.correlationMapNode.ycorrelation[j] = 0.0;
			}
			particle.correlationMapNode.ycorrelation[0] = 1.0;
		}
		tempIndex = 0;
		if (znumber > 1) {
			for (int k = zcount + leftShiftZ - (splineOrder / 2); k <= zcount + leftShiftZ + (splineOrder / 2) + 1; ++k) {
				particle.correlationMapNode.zindex[tempIndex] = k;
				particle.correlationMapNode.zcorrelation[tempIndex] = correlationBspline(particle.coordinates.z, particle.dz, zgrid[0] + (k - 0.5) * deltaZ, zgrid[0] + (k + 0.5) * deltaZ);
				tempIndex++;
			}
		} else {
			for (int k = 0; k < splineOrder + 2; ++k) {
				particle.correlationMapNode.zindex[k] = k;
				particle.correlationMapNode.zcorrelation[k] = 0.0;
			}
			particle.correlationMapNode.zcorrelation[0] = 1.0;
		}
	} else {
		int tempIndex = 0;
		for (int i = xcount - (splineOrder / 2) - 1; i <= xcount + (splineOrder / 2) + 1; ++i) {
			particle.correlationMapNode.xindex[tempIndex] = i;
			particle.correlationMapNode.xcorrelation[tempIndex] = correlationBspline(particle.coordinates.x, particle.dx, xgrid[0] + (i - 0.5) * deltaX, xgrid[0] + (i + 0.5) * deltaX);
			tempIndex++;
		}
		tempIndex = 0;
		if (ynumber > 1) {
			for (int j = ycount - (splineOrder / 2) - 1; j <= ycount + (splineOrder / 2) + 1; ++j) {
				particle.correlationMapNode.yindex[tempIndex] = j;
				particle.correlationMapNode.ycorrelation[tempIndex] = correlationBspline(particle.coordinates.y, particle.dy, ygrid[0] + (j - 0.5) * deltaY, ygrid[0] + (j + 0.5) * deltaY);
				tempIndex++;
			}
		} else {
			for (int j = 0; j < splineOrder + 2; ++j) {
				particle.correlationMapNode.yindex[j] = j;
				particle.correlationMapNode.ycorrelation[j] = 0.0;
			}
			particle.correlationMapNode.ycorrelation[0] = 1.0;
		}
		tempIndex = 0;
		if (znumber > 1) {
			for (int k = zcount - (splineOrder / 2) - 1; k <= zcount + (splineOrder / 2) + 1; ++k) {
				particle.correlationMapNode.zindex[tempIndex] = k;
				particle.correlationMapNode.zcorrelation[tempIndex] = correlationBspline(particle.coordinates.z, particle.dz, zgrid[0] + (k - 0.5) * deltaZ, zgrid[0] + (k + 0.5) * deltaZ);
				tempIndex++;
			}
		} else {
			for (int k = 0; k < splineOrder + 2; ++k) {
				particle.correlationMapNode.zindex[k] = k;
				particle.correlationMapNode.zcorrelation[k] = 0.0;
			}
			particle.correlationMapNode.zcorrelation[0] = 1.0;
		}
	}

	double fullCorrelation = 0;
	for(int i = 0; i < splineOrder+2; ++i) {
		for(int j = 0; j < splineOrder+2; ++j) {
			for(int k = 0; k <splineOrder+2; ++k) {
				fullCorrelation += particle.correlationMapNode.xcorrelation[i]*particle.correlationMapNode.ycorrelation[j]*particle.correlationMapNode.zcorrelation[k];
			}
		}
	}
	if((fullCorrelation < 1.0 - 1E-10) || (fullCorrelation > 1.0 + 1E-10)){
		printf("full correlation node = %20.15g\n", fullCorrelation);
	}
}
