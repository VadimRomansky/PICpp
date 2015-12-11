#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>

#include "simulation.h"
#include "util.h"
#include "output.h"
#include "constants.h"
#include "matrix3d.h"

Simulation::Simulation() {
	debugMode = true;

	currentIteration = 0;
	time = 0;
	particlesNumber = 0;

	theta = 0.6;

	particleEnergy = 0;
	electricFieldEnergy = 0;
	magneticFieldEnergy = 0;

	momentum = Vector3d(0, 0, 0);

	//read input!
	xnumber = 10;

	xsize = 1E15;

	temperature = 0.1E14;
	density = 1.6E-24;

	maxIteration = 1E7;
	maxTime = 1E7;

	particlesPerBin = 10;

	B0 = Vector3d(1E-5, 0, 0);
	E0 = Vector3d(1E-15, 0, 0);

	double concentration = density / massProton;

	plasma_period = sqrt(massElectron / (4 * pi * concentration * sqr(electron_charge))) / (2 * pi);

	double thermal_momentum;
	if(kBoltzman*temperature > massElectron*speed_of_light*speed_of_light) {
		thermal_momentum = kBoltzman*temperature/speed_of_light;
	} else {
		thermal_momentum = sqrt(2*massElectron*kBoltzman*temperature);
	}
	gyroradius = thermal_momentum*speed_of_light / (electron_charge * B0.norm());

	kBoltzman_normalized = kBoltzman * gyroradius * gyroradius / (plasma_period * plasma_period);
	speed_of_light_normalized = speed_of_light * plasma_period / gyroradius;
	speed_of_light_normalized_sqr = speed_of_light_normalized * speed_of_light_normalized;
	electron_charge_normalized = electron_charge * plasma_period / sqrt(cube(gyroradius));

	E0 = E0 / (plasma_period * gyroradius);
	B0 = B0 / (plasma_period * gyroradius);

	density = density * cube(gyroradius);

	Kronecker = Matrix3d(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				LeviCivita[i][j][k] = 0;
			}
		}
	}
	LeviCivita[0][1][2] = 1.0;
	LeviCivita[0][2][1] = -1.0;
	LeviCivita[1][0][2] = -1.0;
	LeviCivita[1][2][0] = 1.0;
	LeviCivita[2][0][1] = 1.0;
	LeviCivita[2][1][0] = -1.0;
}

Simulation::Simulation(double xn, double xsizev, double temp, double rho, double Ex, double Ey, double Ez, double Bx, double By, double Bz, int maxIterations, double maxTimeV, int particlesPerBinV) {
	debugMode = true;

	currentIteration = 0;
	time = 0;
	particlesNumber = 0;

	particleEnergy = 0;
	electricFieldEnergy = 0;
	magneticFieldEnergy = 0;

	momentum = Vector3d(0, 0, 0);


	theta = 0.6;

	xnumber = xn;

	xsize = xsizev;

	temperature = temp;
	density = rho;

	maxIteration = maxIterations;
	maxTime = maxTimeV;

	particlesPerBin = particlesPerBinV;

	B0 = Vector3d(Bx, By, Bz);
	E0 = Vector3d(Ex, Ey, Ez);

	double concentration = density / (massProton + massElectron);

	plasma_period = sqrt(massElectron / (4 * pi * concentration * sqr(electron_charge))) / (2 * pi);
	plasma_period2 = plasma_period;
	double thermal_momentum;
	if(kBoltzman*temperature > massElectron*speed_of_light*speed_of_light) {
		thermal_momentum = kBoltzman*temperature/speed_of_light;
	} else {
		thermal_momentum = sqrt(2*massElectron*kBoltzman*temperature);
	}
	gyroradius = thermal_momentum*speed_of_light / (electron_charge * B0.norm());

	plasma_period = 1.0;
	gyroradius = 1.0;

	kBoltzman_normalized = kBoltzman * plasma_period * plasma_period / (gyroradius * gyroradius);
	speed_of_light_normalized = speed_of_light * plasma_period / gyroradius;
	speed_of_light_normalized_sqr = speed_of_light_normalized * speed_of_light_normalized;
	electron_charge_normalized = electron_charge * plasma_period / sqrt(cube(gyroradius));

	E0 = E0 * (plasma_period * gyroradius);
	B0 = B0 * (plasma_period * gyroradius);

	density = density * cube(gyroradius);

	xsize /= gyroradius;
	printf("xsize/gyroradius = %lf\n", xsize);

	Kronecker = Matrix3d(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				LeviCivita[i][j][k] = 0;
			}
		}
	}
	LeviCivita[0][1][2] = 1.0;
	LeviCivita[0][2][1] = -1.0;
	LeviCivita[1][0][2] = -1.0;
	LeviCivita[1][2][0] = 1.0;
	LeviCivita[2][0][1] = 1.0;
	LeviCivita[2][1][0] = -1.0;

}

Simulation::~Simulation() {

	for (int i = 0; i < xnumber; ++i) {
		delete[] maxwellEquationMatrix[i];
		delete[] maxwellEquationRightPart[i];
	}
	delete[] maxwellEquationMatrix;
	delete[] maxwellEquationRightPart;

	for(int i = 0; i < xnumber; ++i) {
		delete[] divergenceCleaningField[i];
		delete[] divergenceCleaningPotential[i];
		delete[] divergenceCleanUpMatrix[i];
		delete[] divergenceCleanUpRightPart[i];
	}
	delete[] divergenceCleaningField;
	delete[] divergenceCleaningPotential;
	delete[] divergenceCleanUpMatrix;
	delete[] divergenceCleanUpRightPart;

	delete[] Efield;
	delete[] newEfield;
	delete[] tempEfield;
	delete[] implicitField;
	delete[] Bfield;
	delete[] newBfield;

	delete[] xgrid;
	delete[] middleXgrid;
}

void Simulation::initialize() {
	printf("initialization\n");

	deltaX = xsize / (xnumber);

	deltaX2 = deltaX*deltaX;

	for (int i = 0; i <= xnumber; ++i) {
		xgrid[i] = i * deltaX;
	}

	for (int i = 0; i < xnumber; ++i) {
		middleXgrid[i] = (xgrid[i] + xgrid[i + 1]) / 2;
	}

	for (int i = 0; i < xnumber + 1; ++i) {
				Efield[i] = E0;
				newEfield[i] = Efield[i];
				tempEfield[i] = Efield[i];
				implicitField[i] = Efield[i];
	}

	for (int i = 0; i < xnumber; ++i) {
				Bfield[i] = B0;
				Bfield[i].x = B0.x;
				newBfield[i] = Bfield[i];
	}

	checkDebyeParameter();

	double concentration = density / (massProton + massElectron);

	omegaPlasmaProton = sqrt(4*pi*concentration*electron_charge*electron_charge/massProton);
	omegaPlasmaElectron = sqrt(4*pi*concentration*electron_charge*electron_charge/massElectron);

	omegaGyroProton = electron_charge*B0.norm()/(massProton*speed_of_light);
	omegaGyroProton = electron_charge*B0.norm()/(massElectron*speed_of_light);

}

void Simulation::initializeSimpleElectroMagneticWave() {
	E0 = Vector3d(0, 0, 0);
	B0 = Vector3d(0, 0, 0);
	for (int i = 0; i < xnumber; ++i) {
				Bfield[i] = Vector3d(0, 0, 0);
				newBfield[i] = Bfield[i];
	}
	double kw = 2 * pi / xsize;
	double E = 1E-5;

	for (int i = 0; i < xnumber; ++i) {
				Efield[i].x = 0;
				Efield[i].y = E * sin(kw * xgrid[i]);
				Efield[i].z = 0;
				tempEfield[i] = Efield[i];
				newEfield[i] = tempEfield[i];
				implicitField[i] =Efield[i];
	}
	Efield[xnumber] = Efield[0];
	tempEfield[xnumber] = Efield[0];
	newEfield[xnumber] = Efield[0];
	implicitField[xnumber] = Efield[0];

	for(int i = 0; i < xnumber; ++i) {
		Bfield[i].z = E*sin(kw*middleXgrid[i]);
	}

	//double t = 2 * pi / (kw * speed_of_light_normalized);
}

void Simulation::createArrays() {
	printf("creating arrays\n");
	xgrid = new double[xnumber + 1];

	middleXgrid = new double[xnumber];

	Efield = new Vector3d[xnumber + 1];
	newEfield = new Vector3d[xnumber + 1];
	tempEfield = new Vector3d[xnumber + 1];
	implicitField = new Vector3d[xnumber + 1];
	Bfield = new Vector3d[xnumber];
	newBfield = new Vector3d[xnumber];
	//tempBfield = new Vector3d**[xnumber];



	divergenceCleaningField = new double*[xnumber];
	divergenceCleaningPotential = new double*[xnumber];

	for (int i = 0; i < xnumber; ++i) {
		Bfield[i] = Vector3d(0, 0, 0);
		newBfield[i] = Vector3d(0, 0, 0);
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		Efield[i] = Vector3d(0, 0, 0);
		newEfield[i] = Vector3d(0, 0, 0);
		tempEfield[i] = Vector3d(0, 0, 0);
	}

	maxwellEquationMatrix = new std::vector<MatrixElement>*[xnumber];
	maxwellEquationRightPart = new double*[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		maxwellEquationMatrix[i] = new std::vector<MatrixElement>[3];
		maxwellEquationRightPart[i] = new double[3];
	}

	divergenceCleanUpMatrix = new std::vector<MatrixElement>*[xnumber];
	divergenceCleanUpRightPart = new double*[xnumber];

	for(int i = 0; i < xnumber; ++i) {
				divergenceCleaningField[i] = new double[3];
				divergenceCleaningPotential[i] = new double[1];
				divergenceCleanUpMatrix[i] = new std::vector<MatrixElement>[3];
				divergenceCleanUpRightPart[i] = new double[3];
			}
}

void Simulation::createFiles() {
	printf("creating files\n");
	EfieldFile = fopen("./output/Efield.dat", "w");
	fclose(EfieldFile);
	BfieldFile = fopen("./output/Bfield.dat", "w");
	fclose(BfieldFile);
	Xfile = fopen("./output/Xfile.dat", "w");
	fclose(Xfile);
	generalFile = fopen("./output/general.dat", "w");
	fclose(generalFile);
	divergenceErrorFile = fopen("./output/divergence_error.dat", "w");
	fclose(divergenceErrorFile);
	informationFile = fopen("./output/information.dat", "w");
	fclose(informationFile);
}

void Simulation::simulate() {
	createArrays();
	initialize();
	createFiles();
	//createParticles();
	//initializeAlfvenWave();
	initializeSimpleElectroMagneticWave();
	//updateElectroMagneticParameters();
	//cleanupDivergence();
	updateEnergy();

	//updateDeltaT();
	//deltaT = 0;

	while (time*plasma_period < maxTime && currentIteration < maxIteration) {
		printf("iteration number %d time = %15.10g\n", currentIteration, time * plasma_period);

		if (currentIteration % writeParameter == 0) {
			output();
		}
		updateDeltaT();
		evaluateFields();
		//cleanupDivergence();
		evaluateMagneticField();
		updateFields();
		updateEnergy();

		time += deltaT;
		currentIteration++;

	}
}

void Simulation::output() {
	printf("outputing\n");
	EfieldFile = fopen("./output/Efield.dat", "a");
	BfieldFile = fopen("./output/Bfield.dat", "a");
	outputFields(EfieldFile, BfieldFile, Efield, Bfield, xnumber, plasma_period, gyroradius);
	fclose(EfieldFile);
	fclose(BfieldFile);

	Xfile = fopen("./output/Xfile.dat", "w");
	outputGrid(Xfile, xgrid, xnumber);
	fclose(Xfile);

	divergenceErrorFile = fopen("./output/divergence_error.dat", "a");
	outputDivergenceError(divergenceErrorFile, this);
	fclose(divergenceErrorFile);
}

void Simulation::checkDebyeParameter(){
	informationFile = fopen("./output/information.dat", "a");
	double concentration = density/(massProton + massElectron);
	double weight = concentration*volume(0)/particlesPerBin;
	double superParticleCharge = electron_charge_normalized*weight;
	double superParticleConcentration = concentration/weight;
	double superParticleTemperature = temperature*weight;

	double debyeLength = 1/sqrt(4*pi*electron_charge*electron_charge_normalized*concentration/(kBoltzman_normalized*temperature));
	double debyeNumber = 4*pi*cube(debyeLength)*concentration/3;

	if(debyeNumber < 1.0){\
		printf("debye number < 1\n");
		fprintf(informationFile, "debye number < 1\n");
	} else if(debyeNumber < 100.0){
		printf("debye number < 100\n");
		fprintf(informationFile, "debye number < 100\n");
	}
	printf("debye number = %g\n", debyeNumber);
	fprintf(informationFile, "debye number = %g\n", debyeNumber);

	double superParticleDebyeLength = 1/sqrt(4*pi*superParticleCharge*superParticleCharge*superParticleConcentration/(kBoltzman_normalized*superParticleTemperature));
	double superParticleDebyeNumber = 4*pi*cube(superParticleDebyeLength)*superParticleConcentration/3;

	if(superParticleDebyeNumber < 1.0){\
		printf("superparticle debye number < 1\n");
		fprintf(informationFile, "superparticle debye number < 1\n");
	} else if(superParticleDebyeNumber < 100.0){
		printf("superparticle debye number < 100\n");
		fprintf(informationFile, "superparticle debye number < 100\n");
	}
	printf("superparticle debye number = %g\n", superParticleDebyeNumber);
	fprintf(informationFile, "superparticle debye number = %g\n", superParticleDebyeNumber);
	fclose(informationFile);
}

void Simulation::checkCollisionTime(double omega){
	informationFile = fopen("./output/information.dat", "a");
	double concentration = density/(massProton + massElectron);
	double weight = concentration*volume(0)/particlesPerBin;
	double superParticleCharge = electron_charge_normalized*weight;
	double superParticleConcentration = concentration/weight;
	double superParticleTemperature = temperature*weight;

	double qLog = 15;
	double nuElectronIon = 4*sqrt(2*pi)*qLog*power(electron_charge_normalized, 4)*(concentration)/(3*sqrt(massElectron)*sqrt(cube(kBoltzman_normalized*temperature)));
	double collisionlessParameter = omega/nuElectronIon;

	if(collisionlessParameter < 1.0){
		printf("collisionlessParameter < 1\n");
		fprintf(informationFile, "collisionlessParameter < 1\n");
	} else if(collisionlessParameter < 100.0){
		printf("collisionlessParameter < 100\n");
		fprintf(informationFile, "collisionlessParameter < 100\n");
	}
	printf("collisionlessParameter = %g\n", collisionlessParameter);
	fprintf(informationFile, "collisionlessParameter = %g\n", collisionlessParameter);

	double superParticleNuElectronIon = 4*sqrt(2*pi)*qLog*power(electron_charge_normalized*weight, 4)*(superParticleConcentration)/(3*sqrt(massElectron*weight)*sqrt(cube(kBoltzman_normalized*superParticleTemperature)));
	double superParticleCollisionlessParameter = omega/superParticleNuElectronIon;

	if(superParticleCollisionlessParameter < 1.0){
		printf("superParticleCollisionlessParameter < 1\n");
		fprintf(informationFile, "superParticleCollisionlessParameter < 1\n");
	} else if(superParticleCollisionlessParameter < 100.0){
		printf("superParticleCollisionlessParameter < 100\n");
		fprintf(informationFile, "superParticleCollisionlessParameter < 100\n");
	}
	printf("superParticleCollisionlessParameter = %g\n", superParticleCollisionlessParameter);
	fprintf(informationFile, "superParticleCollisionlessParameter = %g\n", superParticleCollisionlessParameter);
	fclose(informationFile);
}

void Simulation::checkMagneticReynolds(double v){
	informationFile = fopen("./output/information.dat", "a");
	double concentration = density/(massProton + massElectron);
	double weight = concentration*volume(0)/particlesPerBin;
	double superParticleCharge = electron_charge_normalized*weight;
	double superParticleConcentration = concentration/weight;
	double superParticleTemperature = temperature*weight;

	double qLog = 15;
	double conductivity = (3*massProton*sqrt(massElectron)*sqrt(cube(kBoltzman_normalized*temperature)))/(sqr(massElectron)*4*sqrt(2*pi)*qLog*sqr(electron_charge_normalized));

	double magneticReynolds = conductivity*v*0.1*xsize;

	if(magneticReynolds < 1.0){
		printf("magneticReynolds < 1\n");
		fprintf(informationFile, "magneticReynolds < 1\n");
	} else if(magneticReynolds <100.0){
		printf("magneticReynolds < 100\n");
		fprintf(informationFile, "magneticReynolds < 100\n");
	}
	printf("magnetic Reynolds = %g\n", magneticReynolds);
	fprintf(informationFile, "magnetic Reynolds = %g\n", magneticReynolds);

	double superParticleConductivity = (3*massProton*weight*sqrt(massElectron*weight)*sqrt(cube(kBoltzman_normalized*superParticleTemperature)))/(sqr(massElectron*weight)*4*sqrt(2*pi)*qLog*sqr(electron_charge_normalized*weight));

	double superParticleMagneticReynolds = superParticleConductivity*v*0.1*xsize;

	if(superParticleMagneticReynolds < 1.0){
		printf("superParticleMagneticReynolds < 1\n");
		fprintf(informationFile, "superParticleMagneticReynolds < 1\n");
	} else if(superParticleMagneticReynolds <100.0){
		printf("superParticleMagneticReynolds < 100\n");
		fprintf(informationFile, "superParticleMagneticReynolds < 100\n");
	}
	printf("superparticle magnetic Reynolds = %g\n", superParticleMagneticReynolds);
	fprintf(informationFile, "superparticle magnetic Reynolds = %g\n", superParticleMagneticReynolds);
	fclose(informationFile);
}

void Simulation::checkDissipation(double k, double alfvenV){
	informationFile = fopen("./output/information.dat", "a");
	double omega = k*alfvenV;

	double concentration = density/(massProton + massElectron);
	double weight = concentration*volume(0)/particlesPerBin;
	double superParticleCharge = electron_charge_normalized*weight;
	double superParticleConcentration = concentration/weight;
	double superParticleTemperature = temperature*weight;

	double qLog = 15;
	double conductivity = (3*massProton*sqrt(massElectron)*sqrt(cube(kBoltzman_normalized*temperature)))/(sqr(massElectron)*4*sqrt(2*pi)*qLog*sqr(electron_charge_normalized));

	double nuMagnetic = speed_of_light_normalized_sqr/(4*pi*conductivity);

	double kdissipation = omega*omega*nuMagnetic/(2*cube(alfvenV));

	if(kdissipation > k){
		printf("kdissipation > k\n");
		fprintf(informationFile, "kdissipation > k\n");
	} else if(kdissipation > 0.1*k){
		printf("kdissipation > 0.1*k\n");
		fprintf(informationFile, "kdissipation > 0.1*k\n");
	}
	printf("kdissipation/k = %g\n", kdissipation/k);
	fprintf(informationFile, "kdissipation/k = %g\n", kdissipation/k);

	double superParticleConductivity = (3*massProton*weight*sqrt(massElectron*weight)*sqrt(cube(kBoltzman_normalized*superParticleTemperature)))/(sqr(massElectron*weight)*4*sqrt(2*pi)*qLog*sqr(electron_charge_normalized*weight));

	double superParticleNuMagnetic = speed_of_light_normalized_sqr/(4*pi*superParticleConductivity);

	double superParticleKdissipation = omega*omega*superParticleNuMagnetic/(2*cube(alfvenV));

	if(superParticleKdissipation > k){
		printf("super particle kdissipation > k\n");
		fprintf(informationFile, "super particle kdissipation > k\n");
	} else if(superParticleKdissipation > 0.1*k){
		printf("super particle kdissipation > 0.1*k\n");
		fprintf(informationFile, "super particle kdissipation > 0.1*k\n");
	}
	printf("super particle kdissipation/k = %g\n", superParticleKdissipation/k);
	fprintf(informationFile, "super particle kdissipation/k = %g\n", superParticleKdissipation/k);
	fclose(informationFile);
}

void Simulation::updateDeltaT() {
	printf("updating time step\n");
	double delta = deltaX;
	deltaT = 0.1 * delta / speed_of_light_normalized;
	double B = B0.norm();
	double E = E0.norm();
	for(int i = 0; i < xnumber; ++i) {
				if(Bfield[i].norm() > B) {
					B = Bfield[i].norm();
				}
	}
	for(int i = 0; i < xnumber; ++i) {
				if(Efield[i].norm() > E) {
					E = Efield[i].norm();
				}
	}

	double concentration = density/(massProton + massElectron);

	omegaPlasmaProton = sqrt(4*pi*concentration*electron_charge*electron_charge/massProton);
	omegaPlasmaElectron = sqrt(4*pi*concentration*electron_charge*electron_charge/massElectron);

	omegaGyroProton = electron_charge*B/(massProton*speed_of_light);
	omegaGyroElectron = electron_charge*B/(massElectron*speed_of_light);

	if(B > 0){
		deltaT = min2(deltaT, 0.1 * massElectron * speed_of_light_normalized / (electron_charge_normalized * B));
	}
	if(E > 0) {
		deltaT = min2(deltaT, 0.1*massElectron * speed_of_light_normalized/(electron_charge_normalized*E));
	}
	//deltaT = 0.005 * massElectron * speed_of_light_normalized / (electron_charge_normalized * B0.norm());
	//deltaT = min2(deltaT, 0.02);
	deltaT = min2(deltaT, 0.05/(omegaPlasmaElectron));
	//deltaT = min2(deltaT, 1E-1);
}

double Simulation::volume(int i) {
	return deltaX;
}

void Simulation::updateEnergy() {
	particleEnergy = 0;
	electricFieldEnergy = 0;
	magneticFieldEnergy = 0;
	energy = 0;

	momentum = Vector3d(0, 0, 0);
	for(int i = 0; i < xnumber+1; ++i) {
				double factor = 1;
				if(i == 0 || i == xnumber) {
					factor = factor/2;
				}
				Vector3d E = Efield[i];
				electricFieldEnergy += E.scalarMult(E)*volume(i)*factor/(8*pi);
	}

	for(int i = 0; i < xnumber; ++i) {
				Vector3d B = Bfield[i];
				magneticFieldEnergy += B.scalarMult(B)*volume(i)/(8*pi);
	}

	for(int i = 0; i < xnumber; ++i) {
				Vector3d E = (Efield[i]  + Efield[i] + Efield[i] + Efield[i] + Efield[i+1] + Efield[i+1] + Efield[i+1] + Efield[i+1])/8;
				Vector3d B = Bfield[i];

				momentum += (E.vectorMult(B)/(4*pi*speed_of_light_normalized))*volume(i);
	}


	particleEnergy *= sqr(gyroradius/plasma_period);
	electricFieldEnergy *= sqr(gyroradius/plasma_period);
	magneticFieldEnergy *= sqr(gyroradius/plasma_period);
	momentum = momentum * gyroradius/plasma_period;


	energy = particleEnergy + electricFieldEnergy + magneticFieldEnergy;
}

Vector3d Simulation::getBfield(int i) {
		if(i == -1){
			i = xnumber - 1;
		} else if (i == xnumber){
			i = 0;
		} else if(i < -1 || i > xnumber){
			printf("i < -1 || i > xnumber in getBfied i = %d\n", i);
			exit(0);
		}

	return Bfield[i];
}

Vector3d Simulation::getTempEfield(int i) {
	if (i < 0) {
		//return Vector3d(0.0, 0.0, 0.0);
		i = 0;
	} else if (i > xnumber) {
		return E0;
	}

	return tempEfield[i];
}

Vector3d Simulation::getEfield(int i) {
	if (i < 0) {
		return Vector3d(0.0, 0.0, 0.0);
	} else if (i > xnumber) {
		return E0;
	}

	return Efield[i];
}