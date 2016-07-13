#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <string>
#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "list"
#include "matrix3d.h"
#include "matrixElement.h"
#include "particle.h"
#include "complex.h"

enum SolverType {EXPLICIT, IMPLICIT};
enum InputType {CGS, Theoretical};

enum BoundaryConditionType {PERIODIC, SUPER_CONDUCTOR_LEFT, FREE_BOTH};

class Simulation{
public:
	std::string outputDir;
	int xnumber;
	int ynumber;
	int znumber;
	int particlesNumber;
	int electronsPerBin;
	int protonsPerBin;
	int positronsPerBin;
	int alphaPerBin;

	double massProton;
	double massElectron;
	double massAlpha;
	double massDeuterium;
	double massHelium3;

	double* concentrations;
	int* particlesPerBin;

	double density;
	double temperature;

	double plasma_period;
	double plasma_period2;
	double scaleFactor;
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
	double eta;

	int shockWavePoint;

	bool debugMode;

	bool newlyStarted;
	bool preserveChargeLocal;
	bool preserveChargeGlobal;

	InputType inputType;
	SolverType solverType;
	BoundaryConditionType boundaryConditionType;
	int maxwellEquationMatrixSize;

	double particleEnergy;
	double electricFieldEnergy;
	double magneticFieldEnergy;
	double energy;

	double theoreticalEnergy;
	Vector3d theoreticalMomentum;

	int chargeBalance;

	double extJ;

	Vector3d momentum;

	ParticleTypeContainer* types;
	int typesNumber;

	double*** electronConcentration;
	double*** protonConcentration;
	double*** chargeDensity;
	Vector3d*** velocityBulkProton;
	Vector3d*** velocityBulkElectron;

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
	Vector3d*** externalElectricFlux;
	double*** electricDensity;
	Matrix3d*** dielectricTensor;
	Matrix3d*** pressureTensor;
	Vector3d*** divPressureTensor;

	Vector3d*** Efield;
	Vector3d*** Bfield;

	Vector3d*** newEfield;
	Vector3d*** newBfield;

	Vector3d*** tempEfield;
	//Vector3d*** tempBfield;

	Vector3d*** explicitEfield;
	Vector3d*** rotB;
	Vector3d*** Ederivative;

	double**** divergenceCleaningField;
	double**** divergenceCleaningPotential;

	double*** divergenceCleaningPotentialFourier;

	std::vector<Particle*> particles;
	std::vector<Particle*> escapedParticles;

	std::vector<Particle*>*** particlesInEbin;
	std::vector<Particle*>*** particlesInBbin;

	Matrix3d Kronecker;
	double LeviCivita[3][3][3];

	FILE* protonTraectoryFile;
	FILE* electronTraectoryFile;
	FILE* distributionFileProton;
	FILE* distributionFileElectron;
	FILE* distributionFileProtonUpstream;
	FILE* distributionFileElectronUpstream;
	FILE* anisotropyFileProton;
	FILE* anisotropyFileElectron;
	FILE* anisotropyFileAlpha;
	FILE* anisotropyFilePositron;
	FILE* anisotropyFileDeuterium;
	FILE* anisotropyFileHelium3;
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
	FILE* particlePositronsFile;
	FILE* particleAlphaFile;
	FILE* particleDeuteriumFile;
	FILE* particleHelium3File;

	FILE* rotBFile;
	FILE* EderivativeFile;

	FILE* maxwellMatrixFile;

	FILE* backupGeneralFile;
	FILE* backupParticlesFile;
	FILE* backupEfieldFile;
	FILE* backupBfieldFile;

	//FILE* outputEverythingFile;

	FILE* errorLogFile;

	//Simulation();
	Simulation();
	Simulation(int xn, int yn, int zn, double xsizev, double ysizev, double zsizev, double temp,
			   double Vx, double Vy, double Vz, double Ex, double Ey, double Ez, double Bx, double By,
			   double Bz, int maxIterations, double maxTimeV, int typesNumberV, int* particlesperBin, double* concentrations,
			   int inputType);
	~Simulation();

	void initialize();
	void initializeSimpleElectroMagneticWave();
	void initializeRotatedSimpleElectroMagneticWave(int wavesCount);
	void checkFrequency(double omega);
	void initializeAlfvenWaveX(int wavesCount, double amplitudeRelation);
	void initializeAlfvenWaveY(int wavesCount, double amplitudeRelation);
	void initializeRotatedAlfvenWave(int wavesCount, double amplitudeRelation);
	void initializeLangmuirWave();
	void initializeTwoStream();
	void initializeExternalFluxInstability();
	void initializeFluxFromRight();
	void fieldsLorentzTransitionX(const double& v);
	void initializeShockWave();
	void initializeAnisotropic();
	void initializeKolmogorovSpectrum(int start, int end);
	void createArrays();
	void createParticleTypes(double *concentrations, int *particlesPerBin);
	void createFiles();
	void simulate();
	void output();
	void outputBackup();
	void rescaleConstants();

	void rescaleConstantsToTheoretical();
	void checkDebyeParameter();
	void checkGyroRadius();
	void checkCollisionTime(double omega);
	void checkMagneticReynolds(double v);
	void checkDissipation(double k, double alfvenV);

	Matrix3d evaluateAlphaRotationTensor(double beta, Vector3d velocity, Vector3d EField, Vector3d BField);
	void updateDeltaT();
	void evaluateFields();
	void updateEfield();
	void updateBfield();
	void evaluateExplicitDerivative();
	void checkEquationMatrix(std::vector<MatrixElement>**** matrix, int lnumber);
	void createSuperConductorLeftEquation(int j, int k);
	void createFreeRightEquation(int j, int k);
	void createFreeRightEquationX(int j, int k, Vector3d& rightPart);
	void createFreeRightEquationY(int j, int k, Vector3d& rightPart);
	void createFreeRightEquationZ(int j, int k, Vector3d& rightPart);
	void createInternalEquationX(int i, int j, int k);
	void createInternalEquationY(int i, int j, int k);
	void createInternalEquationZ(int i, int j, int k);
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

	Complex*** evaluateFourierTranslation(double*** a);

	double*** evaluateReverceFourierTranslation(Complex*** a);
	void resetNewTempFields();
	double volumeE(int i, int j, int k);
	double volumeB(int i, int j, int k);
	void checkParticleInBox(Particle& particle);
	void checkParticlesInBin();
	void updateElectroMagneticParameters();
	void smoothDensity();
	void addReflectedParticleToElectroMagneticParameters(const Particle* particle, int j, int k);
	void updateDensityParameters();
	void updateEnergy();
	void updateFields();
	void updateParameters();
	void updateExternalFlux();
	Vector3d evaluateRotB(int i,int j, int k);
	Vector3d evaluateRotTempE(int i, int j, int k);
	Vector3d evaluateRotE(int i, int j, int k);
	Vector3d evaluateRotNewE(int i, int j, int k);
	double evaluateDivE(int i, int j, int k);
	double evaluateDivCleaningE(int i, int j, int k);
	double evaluateDivTempE(int i, int j, int k);
	double evaluateDivNewE(int i, int j, int k);
	double evaluateDivFlux(int i, int j, int k);

	Vector3d evaluateDivPressureTensor(int i, int j, int k);
	Vector3d evaluateGradDensity(int i, int j, int k);
	Vector3d getBfield(int i, int j, int k);
	Vector3d getTempEfield(int i, int j, int k);
	Vector3d getNewEfield(int i, int j, int k);
	Vector3d getEfield(int i, int j, int k);

	Matrix3d getPressureTensor(int i, int j, int k);
	double getDensity(int i, int j, int k);
	void createParticles();
	void moveToPreserveChargeLocal();
	void addToPreserveChargeGlobal();
	Particle* getFirstProton();
	Particle* getFirstElectron();
	Particle* getLastProton();
	Particle* getLastElectron();
	Particle* getProton(int n);

	Particle* getElectron(int n);
	Particle* createParticle(int n, int i, int j, int k, double weight, ParticleTypes type, ParticleTypeContainer typeContainer, double temperatureX, double temperatureY, double temperatureZ);
	void moveParticles();
	void removeEscapedParticles();
	void moveParticle(Particle* particle);
	void correctParticlePosition(Particle* particle);
	void correctParticlePosition(Particle& particle);


	void evaluateParticlesRotationTensor();
	void injectNewParticles(int count, ParticleTypeContainer typeContainer, double length);
	void collectParticlesIntoBins();
	void pushParticleIntoEbin(Particle* particle, int i, int j, int k);
	void pushParticleIntoBbin(Particle* particle, int i, int j, int k);
	bool particleCrossBbinX(Particle& particle, int i);
	bool particleCrossBbinY(Particle& particle, int j);
	bool particleCrossBbinZ(Particle& particle, int k);
	bool particleCrossBbin(Particle& particle, int i, int j, int k);
	bool particleCrossEbinX(Particle& particle, int i);
	bool particleCrossEbinY(Particle& particle, int j);
	bool particleCrossEbinZ(Particle& particle, int k);
	bool particleCrossEbin(Particle& particle, int i, int j, int k);
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
	double correlationWithBbin(Particle& particle, int i, int j, int k);

	double correlationWithEbin(Particle& particle, int i, int j, int k);

	double correlationBspline(const double& x, const double& dx, const double& leftx, const double& rightx);
};

#endif
