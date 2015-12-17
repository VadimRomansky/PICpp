#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "list"
#include "matrix3d.h"
#include "matrixElement.h"
#include "particle.h"

enum SolverType {EXPLICIT, IMPLICIT};

enum BoundaryConditionType {PERIODIC, SUPER_CONDUCTOR_LEFT, FREE_BOTH};

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

	int shockWavePoint;

	bool debugMode;

	bool newlyStarted;

	SolverType solverType;
	BoundaryConditionType boundaryConditionType;
	int maxwellEquationMatrixSize;

	double particleEnergy;
	double electricFieldEnergy;
	double magneticFieldEnergy;
	double energy;

	double theoreticalEnergy;
	Vector3d theoreticalMomentum;

	double extJ;

	Vector3d momentum;

	double* electronConcentration;
	double* protonConcentration;
	double* chargeDensity;
	Vector3d* velocityBulkProton;
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

	//for debug alven wave only/////
	double omega;
	double VyamplitudeProton;
	double VyamplitudeElectron;
	double VzamplitudeProton;
	double VzamplitudeElectron;
	double Eyamplitude;
	double Ezamplitude;
	double Byamplitude;
	double Bzamplitude;
	/////////////////////////

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
	Vector3d* divPressureTensor;

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
	std::vector<Particle*> escapedParticles;

	std::vector<Particle*>* particlesInEbin;
	std::vector<Particle*>* particlesInBbin;

	Matrix3d Kronecker;
	double LeviCivita[3][3][3];

	FILE* protonTraectoryFile;
	FILE* electronTraectoryFile;
	FILE* distributionFileProton;
	FILE* distributionFileElectron;
	FILE* distributionFileProtonUpstream;
	FILE* distributionFileElectronUpstream;
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
	FILE* particleProtonsFile;
	FILE* particleElectronsFile;

	FILE* rotBFile;
	FILE* EderivativeFile;

	FILE* backupGeneralFile;
	FILE* backupParticlesFile;
	FILE* backupEfieldFile;
	FILE* backupBfieldFile;

	FILE* errorLogFile;

	//Simulation();
	Simulation();
	Simulation(int xn, double xsizev, double temp, double rho, double Vx, double Vy, double Vz, double Ex, double Ey, double Ez, double Bx, double By, double Bz, int maxIterations, double maxTimeV, int particlesPerBinV);
	~Simulation();

	void initialize();
	void initializeSimpleElectroMagneticWave();
	void checkFrequency(double omega);
	void initializeAlfvenWave(int wavesCount, double amplitudeRelation);
	void initializeLangmuirWave();
	void initializeTwoStream();
	void initializeExternalFluxInstability();
	void initializeFluxFromRight();
	void fieldsLorentzTransitionX(const double& v);
	void initializeShockWave();
	void initializeKolmogorovSpectrum(int start, int end);
	void createArrays();
	void createFiles();
	void simulate();
	void output();
	void outputBackup();
	void rescaleConstants();

	void checkDebyeParameter();
	void checkGyroRadius();
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
	void createSuperConductorLeftEquation();
	void createFreeRightEquation();
	void createFreeLeftEquation();
	void createFreeRightEquationX(Vector3d& rightPart);
	void createFreeRightEquationY(Vector3d& rightPart);
	void createFreeRightEquationZ(Vector3d& rightPart);
	void createInternalEquationX(int i);
	void createInternalEquationY(int i);
	void createInternalEquationZ(int i);
	void createInternalEquation(int i);
	void evaluateMaxwellEquationMatrix();
	void evaluateMagneticField();


	void updateBoundaries();
	void updateBoundariesOldField();
	void updateBoundariesNewField();
	void resetNewTempFields();


	void cleanupDivergence();
	void updateFieldByCleaning();
	void evaluateDivergenceCleaningField();
	void createDivergenceCleanupInternalEquation(int i);
	void createDivergenceCleanupLeftEquation();
	void createDivergenceCleanupRightEquation();
	double cleanUpRightPart(int i);


	void fourierFilter();
	void fourierFilter(int startPoint, int endPoint);



	double volumeE(int i);
	double volumeB(int i);
	void checkParticleInBox(Particle& particle);
	void checkParticlesInBin();
	void updateElectroMagneticParameters();
	void smoothDensity();
	void addReflectedParticleToElectroMagneticParameters(const Particle* particle);
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
	Particle* getLastProton();
	Particle* getLastElectron();
	Particle* createParticle(int n, int i, double weight, ParticleTypes type, double localTemperature);

	void moveParticles();
	void removeEscapedParticles();
	void moveParticle(Particle* particle);
	void correctParticlePosition(Particle* particle);
	void correctParticlePosition(Particle& particle);
	void moveParticleNewtonIteration(Particle* particle, double* const oldCoordinates, double* const tempCoordinates, double* const newCoordinates);
	void evaluateParticlesRotationTensor();
	void injectNewParticles(int count, double length);
	void scatterParticle(Particle* particle);

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
	double correlationWithBbin(Particle& particle, int i);
	double correlationWithEbin(Particle& particle, int i);
	double correlationBspline(const double& x, const double& dx, const double& leftx, const double& rightx);
};

#endif
