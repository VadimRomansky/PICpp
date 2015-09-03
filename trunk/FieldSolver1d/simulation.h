#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "matrix3d.h"
#include "matrixElement.h"


class Simulation{
public:
	int xnumber;
	int particlesNumber;
	int particlesPerBin;

	double density;
	double temperature;

	double plasma_period;
	double plasma_period2;
	double gyroradius;

	double speed_of_light_normalized;
	double speed_of_light_normalized_sqr;
	double kBoltzman_normalized;
	double electron_charge_normalized;

	double time;
	double maxTime;
	int currentIteration;
	int maxIteration;
	double xsize;

	double deltaT;

	double deltaX;

	double deltaX2;

	double theta;

	bool debugMode;

	double particleEnergy;
	double electricFieldEnergy;
	double magneticFieldEnergy;
	double energy;

	Vector3d momentum;


	Vector3d V0;

	Vector3d B0;
	Vector3d E0;
	Matrix3d pressureTensor0;

	double omegaPlasmaProton;
	double omegaPlasmaElectron;
	double omegaGyroProton;
	double omegaGyroElectron;

	double* xgrid;

	double* middleXgrid;


	std::vector<MatrixElement>** maxwellEquationMatrix;
	double** maxwellEquationRightPart;

	std::vector<MatrixElement>** divergenceCleanUpMatrix;
	double** divergenceCleanUpRightPart;

	Vector3d* Efield;
	Vector3d* Bfield;
	Vector3d* implicitField;

	Vector3d* newEfield;
	Vector3d* newBfield;

	Vector3d* tempEfield;
	//Vector3d*** tempBfield;

	double** divergenceCleaningField;
	double** divergenceCleaningPotential;

	Matrix3d Kronecker;
	double LeviCivita[3][3][3];

	FILE* EfieldFile;
	FILE* BfieldFile;
	FILE* Xfile;
	FILE* generalFile;
	FILE* divergenceErrorFile;
	FILE* informationFile;

	Simulation();
	Simulation(double xn, double xsizev, double temp, double rho, double Ex, double Ey, double Ez, double Bx, double By, double Bz, int maxIterations, double maxTimeV, int particlesPerBinV);
	~Simulation();

	void initialize();
	void initializeSimpleElectroMagneticWave();
	void createArrays();
	void createFiles();
	void simulate();
	void output();

	void checkDebyeParameter();
	void checkCollisionTime(double omega);
	void checkMagneticReynolds(double v);
	void checkDissipation(double k, double alfvenV);

	void updateDeltaT();

	void evaluateFields();
	void updateEfield();
	void updateBfield();
	void checkEquationMatrix(std::vector<MatrixElement>** matrix, int lnumber);
	void createPerfectConductaryBoundaryCondition();
	void createInternalEquation(int i);
	void evaluateMaxwellEquationMatrix();
	void evaluateMagneticField();
	void updateBoundaries();
	void updateBoundariesOldField();
	void updateBoundariesNewField();
	void cleanupDivergence();
	void updateFieldByCleaning();
	void evaluateDivergenceCleaningField();
	void createDivergenceCleanupInternalEquation(int i);
	void createDivergenceCleanupLeftEquation();
	void createDivergenceCleanupRightEquation();
	double cleanUpRightPart(int i);

	double volume(int i);

	void updateEnergy();
	void updateFields();
	Vector3d evaluateRotB(int i);
	Vector3d evaluateRotTempE(int i);
	Vector3d evaluateRotE(int i);
	Vector3d evaluateRotNewE(int i);
	double evaluateDivE(int i);
	double evaluateDivCleaningE(int i);
	double evaluateDivTempE(int i);
	double evaluateDivNewE(int i);

	Vector3d getBfield(int i);
	Vector3d getTempEfield(int i);
	Vector3d getEfield(int i);
};

#endif
