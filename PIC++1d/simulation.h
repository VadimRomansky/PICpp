#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "matrix3d.h"
#include "matrixElement.h"
#include "particle.h"

enum SolverType {EXPLICIT, IMPLICIT};

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
	double fieldScale;

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
	double eta;

	bool debugMode;

	bool newlyStarted;

	SolverType solverType;
	int maxwellEquationMatrixSize;

	double particleEnergy;
	double electricFieldEnergy;
	double magneticFieldEnergy;
	double energy;

	double extJ;

	Vector3d momentum;

	double* electronConcentration;
	double* protonConcentration;
	double* chargeDensity;
	Vector3d* velocityBulk;
	Vector3d* velocityBulkElectron;

	Vector3d maxEfield;
	Vector3d maxBfield;


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

	Vector3d* electricFlux;
	Vector3d* externalElectricFlux;
	double* electricDensity;
	Matrix3d* dielectricTensor;
	Matrix3d* pressureTensor;

	Vector3d* Efield;
	Vector3d* Bfield;

	Vector3d* newEfield;
	Vector3d* newBfield;

	Vector3d* tempEfield;
	//Vector3d*** tempBfield;

	Vector3d* explicitEfield;
	Vector3d* rotB;
	Vector3d* Ederivative;

	double** divergenceCleaningField;
	double** divergenceCleaningPotential;

	std::vector<Particle*> particles;

	std::vector<Particle*>* particlesInEbin;
	std::vector<Particle*>* particlesInBbin;

	Matrix3d Kronecker;
	double LeviCivita[3][3][3];

	FILE* protonTraectoryFile;
	FILE* electronTraectoryFile;
	FILE* distributionFile;
	FILE* EfieldFile;
	FILE* BfieldFile;
	FILE* Xfile;
	FILE* Yfile;
	FILE* Zfile;
	FILE* generalFile;
	FILE* densityFile;
	FILE* divergenceErrorFile;
	FILE* velocityFile;
	FILE* velocityElectronFile;
	FILE* fluxFile;
	FILE* dielectricTensorFile;
	FILE* informationFile;

	FILE* rotBFile;
	FILE* EderivativeFile;

	FILE* backupGeneralFile;
	FILE* backupParticlesFile;
	FILE* backupEfieldFile;
	FILE* backupBfieldFile;

	FILE* errorLogFile;

	//Simulation();
	Simulation();
	Simulation(double xn, double xsizev, double temp, double rho, double Ex, double Ey, double Ez, double Bx, double By, double Bz, int maxIterations, double maxTimeV, int particlesPerBinV);
	~Simulation();

	void initialize();
	void initializeSimpleElectroMagneticWave();
	void checkFrequency(double omega);
	void initializeAlfvenWave();
	void initializeLangmuirWave();
	void initializeTwoStream();
	void initializeExternalFluxInstability();
	void createArrays();
	void createFiles();
	void simulate();
	void output();
	void outputBackup();
	void rescaleConstants();

	void checkDebyeParameter();
	void checkCollisionTime(double omega);
	void checkMagneticReynolds(double v);
	void checkDissipation(double k, double alfvenV);
	Matrix3d evaluateAlphaRotationTensor(double beta, Vector3d velocity, Vector3d EField, Vector3d BField);
	void updateDeltaT();

	void evaluateFields();
	void smoothEfield();
	void updateEfield();
	void updateBfield();
	void evaluateExplicitDerivative();
	void checkEquationMatrix(std::vector<MatrixElement>** matrix, int lnumber);
	void createPerfectConductaryBoundaryCondition();
	void createInternalEquationX(int i);
	void createInternalEquationY(int i);
	void createInternalEquationZ(int i);
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

	void resetNewTempFields();

	double volume(int i);
	void checkParticleInBox(Particle& particle);
	void checkParticlesInBin();
	void updateElectroMagneticParameters();
	void updateDensityParameters();
	void updateEnergy();
	void updateFields();
	void updateParameters();
	void updateExternalFlux();
	Vector3d evaluateRotB(int i);
	Vector3d evaluateRotTempE(int i);
	Vector3d evaluateRotE(int i);
	Vector3d evaluateRotNewE(int i);
	double evaluateDivE(int i);
	double evaluateDivCleaningE(int i);
	double evaluateDivTempE(int i);
	double evaluateDivNewE(int i);
	double evaluateDivFlux(int i);
	Vector3d evaluateDivPressureTensor(int i);
	Vector3d evaluateGradDensity(int i);

	Vector3d getBfield(int i);
	Vector3d getTempEfield(int i);
	Vector3d getNewEfield(int i);
	Vector3d getEfield(int i);
	Matrix3d getPressureTensor(int i);
	double getDensity(int i);
	void smoothFlux();
	void smoothEderivative();

	void createParticles();
	Particle* getFirstProton();
	Particle* getFirstElectron();
	Particle* createParticle(int n, int i, double weight, ParticleTypes type);

	void moveParticles();
	void moveParticle(Particle* particle);
	void correctParticlePosition(Particle* particle);
	void correctParticlePosition(Particle& particle);
	void moveParticleNewtonIteration(Particle* particle, double* const oldCoordinates, double* const tempCoordinates, double* const newCoordinates);
	void evaluateParticlesRotationTensor();

	void collectParticlesIntoBins();
	void pushParticleIntoEbin(Particle* particle, int i);
	void pushParticleIntoBbin(Particle* particle, int i);
	bool particleCrossBbin(Particle& particle, int i);
	bool particleCrossEbin(Particle& particle, int i);
	Vector3d correlationTempEfield(Particle* particle);
	Vector3d correlationNewEfield(Particle* particle);
	Vector3d correlationBfield(Particle* particle);
	Vector3d correlationNewBfield(Particle* particle);
	Vector3d correlationEfield(Particle* particle);
	Vector3d correlationTempEfield(Particle& particle);
	Vector3d correlationNewEfield(Particle& particle);
	Vector3d correlationBfield(Particle& particle);
	Vector3d correlationNewBfield(Particle& particle);
	Vector3d correlationEfield(Particle& particle);
	Vector3d correlationFieldWithBbin(Particle& particle, int i);
	Vector3d correlationFieldWithEbin(Particle& particle, int i);
	Vector3d correlationFieldWithTempEbin(Particle& particle, int i);
	Vector3d correlationFieldWithNewEbin(Particle& particle, int i);
	double correlationWithBbin(Particle& particle, int i);
	double correlationWithEbin(Particle& particle, int i);
	double correlationBspline(const double& x, const double& dx, const double& leftx, const double& rightx);
};

#endif
