#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "matrix3d.h"
#include "matrixElement.h"
#include "particle.h"

enum ParticleTypes;

enum BoundaryConditionTypes {SUPERCONDUCTERLEFT, PERIODIC};


class Simulation{
public:
	int xnumber;
	int ynumber;
	int znumber;
	int particlesNumber;
	int particlesPerBin;

	double density;
	double temperature;

	BoundaryConditionTypes boundaryConditionType;

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
	double ysize;
	double zsize;

	double deltaT;

	double deltaX;
	double deltaY;
	double deltaZ;

	double deltaX2;
	double deltaY2;
	double deltaZ2;

	double theta;

	bool debugMode;

	double particleEnergy;
	double electricFieldEnergy;
	double magneticFieldEnergy;
	double energy;

	Vector3d momentum;

	double*** electronConcentration;
	double*** protonConcentration;
	double*** chargeDensity;
	Vector3d*** velocityBulk;
	Vector3d*** velocityBulkElectron;


	Vector3d V0;

	Vector3d B0;
	Vector3d E0;
	Matrix3d pressureTensor0;

	double omegaPlasmaProton;
	double omegaPlasmaElectron;
	double omegaGyroProton;
	double omegaGyroElectron;

	double* xgrid;
	double* ygrid;
	double* zgrid;
	double* middleXgrid;
	double* middleYgrid;
	double* middleZgrid;

	std::vector<MatrixElement>**** maxwellEquationMatrix;
	double**** maxwellEquationRightPart;

	std::vector<MatrixElement>**** divergenceCleanUpMatrix;
	double**** divergenceCleanUpRightPart;

	Vector3d*** electricFlux;
	double*** electricDensity;
	Matrix3d*** dielectricTensor;
	Matrix3d*** pressureTensor;

	Vector3d*** Efield;
	Vector3d*** Bfield;

	Vector3d*** newEfield;
	Vector3d*** newBfield;

	Vector3d*** tempEfield;
	//Vector3d*** tempBfield;

	double**** divergenceCleaningField;
	double**** divergenceCleaningPotential;

	std::vector<Particle*> particles;

	std::vector<Particle*>*** particlesInEbin;
	std::vector<Particle*>*** particlesInBbin;

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
	FILE* informationFile;

	Simulation();
	Simulation(double xn, double yn, double zn, double xsizev, double ysizev, double zsizev, double temp, double rho, double Ex, double Ey, double Ez, double Bx, double By, double Bz, int maxIterations, double maxTimeV, int particlesPerBinV);
	~Simulation();

	void initialize();
	void initializeSimpleElectroMagneticWave();
	void initializeAlfvenWave();
	void createArrays();
	void createFiles();
	void simulate();
	void output();

	void checkDebyeParameter();
	void checkCollisionTime(double omega);
	void checkMagneticReynolds(double v);
	void checkDissipation(double k, double alfvenV);

	void updateDeltaT();
	void createParticles();
	Particle* createParticle(int n, int i, int j, int k, double weight, ParticleTypes type);
	Particle* getFirstProton();
	Particle* getFirstElectron();

	Vector3d correlationTempEfield(Particle* particle);
	Vector3d correlationBfield(Particle* particle);
	Vector3d correlationTempEfield(Particle& particle);
	Vector3d correlationBfield(Particle& particle);
	Vector3d correlationEfield(Particle* particle);
	Vector3d correlationEfield(Particle& particle);
	
	Vector3d correlationFieldWithBbin(Particle& particle, int i, int j, int k);
	Vector3d correlationFieldWithEbin(Particle& particle, int i, int j, int k);
	Vector3d correlationFieldWithTempEbin(Particle& particle, int i, int j, int k);
	double correlationWithBbin(Particle& particle, int i, int j, int k);
	double correlationWithEbin(Particle& particle, int i, int , int k);
	double correlationBspline(const double& x, const double&  dx, const double& leftx, const double& rightx);

    Matrix3d evaluateAlphaRotationTensor(double beta, Vector3d velocity, Vector3d EField, Vector3d BField); //see Noguchi

	void moveParticles();
	void correctParticlePosition(Particle* particle);
	void moveParticle(Particle* particle);
	void moveParticleNewtonIteration(Particle* particle, double* const oldCoordinates, double* const tempCoordinates, double* const newCoordinates);
	void evaluateParticlesRotationTensor();

	void evaluateFields();
	void updateEfield();
	void updateBfield();
	void checkEquationMatrix(std::vector<MatrixElement>**** matrix, int lnumber);
	void createPerfectConductaryBoundaryCondition(int j, int k);
	void createInternalEquationX(int i, int j, int k, Vector3d& rightPart);
	void createInternalEquationY(int i, int j, int k, Vector3d& rightPart);
	void createInternalEquationZ(int i, int j, int k, Vector3d& rightPart);
	void createInternalEquation(int i, int j, int k);
	void evaluateMaxwellEquationMatrix();
	void evaluateMagneticField();
	void updateBoundaries();
	void updateBoundariesOldField();
	void updateBoundariesNewField();
	void cleanupDivergence();
	void updateFieldByCleaning();
	void evaluateDivergenceCleaningField();
	void createDivergenceCleanupInternalEquation(int i, int j, int k);
	void createDivergenceCleanupLeftEquation(int j, int k);
	void createDivergenceCleanupRightEquation(int j, int k);
	double cleanUpRightPart(int i, int j, int k);

	double volume(int i, int j, int k);

	void collectParticlesIntoBins();
	void pushParticleIntoEbin(Particle* particle, int i, int j, int k);
	void pushParticleIntoBbin(Particle* particle, int i, int j, int k);
	bool particleCrossBbin(Particle& particle, int i, int j, int k);
	bool particleCrossEbin(Particle& particle, int i, int j, int k);
	void checkParticleInBox(Particle& particle);

	void updateElectroMagneticParameters();
	void updateDensityParameters();
	void updateEnergy();
	void updateFields();
	double evaluateDivFlux(int i, int j, int k);
	Vector3d evaluateRotB(int i, int j, int k);
	Vector3d evaluateRotTempE(int i, int j, int k);
	Vector3d evaluateRotE(int i, int j, int k);
	Vector3d evaluateRotNewE(int i, int j, int k);
	double evaluateDivE(int i, int j, int k);
	double evaluateDivCleaningE(int i, int j, int k);
	double evaluateDivTempE(int i, int j, int k);
	double evaluateDivNewE(int i, int j, int k);
	Vector3d evaluateDivPressureTensor(int i, int j, int k);
	Vector3d evaluateGradDensity(int i, int j, int k);
	//Vector3d evaluateGradPotential(int i, int j, int k);

	Vector3d getBfield(int i, int j, int k);
	Vector3d getTempEfield(int i, int j, int k);
	Vector3d getEfield(int i, int j, int k);
	Matrix3d getPressureTensor(int i, int j, int k);
	double getDensity(int i, int j, int k);
};

#endif
