#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <string>
#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "list"
#include "largeVectorBasis.h"

class MatrixElement;
class Complex;

enum SolverType {EXPLICIT, IMPLICIT};

enum InputType {CGS, Theoretical};

enum BoundaryConditionType {PERIODIC, SUPER_CONDUCTOR_LEFT, FREE_BOTH};

class Matrix3d;
class Vector3d;
class Particle;
class ParticleTypeContainer;
enum ParticleTypes;

class Simulation {
private:
	bool arrayCreated;
public:

	std::string outputDir;
	int xnumberGeneral;
	int ynumberGeneral;
	int znumberGeneral;
	int firstAbsoluteXindex;
	double leftX;
	double rightX;

	int xnumber;
	int ynumber;
	int znumber;
	int particlesNumber;

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

	double speed_of_light_normalized;
	double speed_of_light_normalized_sqr;
	double kBoltzman_normalized;
	double electron_charge_normalized;

	double time;
	double maxTime;
	int currentIteration;
	int maxIteration;
	double xsizeGeneral;
	double ysizeGeneral;
	double zsizeGeneral;
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
	double shockWaveX;

	bool debugMode;

	bool newlyStarted;
	//bool preserveChargeLocal;
	bool preserveChargeGlobal;

	bool timing;

	int verbosity;

	InputType inputType;
	SolverType solverType;
	BoundaryConditionType boundaryConditionType;
	int maxwellEquationMatrixSize;

	double particleEnergy;
	double electricFieldEnergy;
	double magneticFieldEnergy;
	double energy;

	double theoreticalEnergy;
	double generalTheoreticalEnergy;
	Vector3d theoreticalMomentum;
	Vector3d generalTheoreticalMomentum;

	int chargeBalance;

	double extJ;

	Vector3d globalMomentum;

	ParticleTypeContainer* types;
	int typesNumber;

	double**** particleConcentrations;
	double**** additionalParticleConcentrationsLeft;
	double**** additionalParticleConcentrationsRight;
	double*** chargeDensity;
	double*** chargeDensityMinus;
	double*** tempCellParameter;
	double*** additionalChargeDensityLeft;
	double*** additionalChargeDensityMinusLeft;
	double*** additionalChargeDensityRight;
	double*** additionalChargeDensityMinusRight;
	Vector3d**** particleBulkVelocities;
	Vector3d**** additionalParticleBulkVelocitiesLeft;
	Vector3d**** additionalParticleBulkVelocitiesRight;

	Vector3d maxEfield;
	Vector3d maxBfield;
	double meanSquaredEfield[3];


	Vector3d V0;

	Vector3d B0;
	Vector3d E0;

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

	double**** gmresOutput;
	LargeVectorBasis* gmresMaxwellBasis;
	LargeVectorBasis* gmresCleanupBasis;


	std::vector<MatrixElement>**** maxwellEquationMatrix;
	double**** maxwellEquationRightPart;

	std::vector<MatrixElement>**** divergenceCleanUpMatrix;
	double**** divergenceCleanUpRightPart;

	Vector3d*** electricFlux;
	Vector3d*** electricFluxMinus;
	Vector3d*** additionalElectricFluxLeft;
	Vector3d*** additionalElectricFluxMinusLeft;
	Vector3d*** additionalElectricFluxRight;
	Vector3d*** additionalElectricFluxMinusRight;
	Vector3d*** externalElectricFlux;
	double*** chargeDensityHat;
	double*** additionalChargeDensityHatLeft;
	double*** additionalChargeDensityHatRight;
	Matrix3d*** dielectricTensor;
	Matrix3d*** additionalDielectricTensorLeft;
	Matrix3d*** additionalDielectricTensorRight;
	Matrix3d*** pressureTensor;
	Matrix3d*** additionalPressureTensorLeft;
	Matrix3d*** additionalPressureTensorRight;
	Vector3d*** divPressureTensor;
	Vector3d*** additionalDivPressureTensorLeft;
	Vector3d*** additionalDivPressureTensorRight;

	int additionalBinNumber;

	Vector3d*** Efield;
	Vector3d*** additionalEfieldLeft;
	Vector3d*** additionalEfieldRight;
	Vector3d*** Bfield;
	Vector3d*** additionalBfieldLeft;
	Vector3d*** additionalBfieldRight;

	Vector3d*** newEfield;
	Vector3d*** additionalNewEfieldLeft;
	Vector3d*** additionalNewEfieldRight;
	Vector3d*** newBfield;
	Vector3d*** additionalNewBfieldLeft;
	Vector3d*** additionalNewBfieldRight;

	Vector3d*** tempEfield;
	Vector3d*** additionalTempEfieldLeft;
	Vector3d*** additionalTempEfieldRight;

	Vector3d*** smoothingEfield;
	Vector3d*** smoothingBfield;

	//Vector3d*** tempBfield;

	Vector3d*** explicitEfield;
	Vector3d*** rotB;
	Vector3d*** Ederivative;

	Vector3d*** rotE;
	Vector3d*** Bderivative;

	double**** divergenceCleaningField;
	double**** divergenceCleaningPotential;
	double**** tempDivergenceCleaningPotential;

	double*** divergenceCleaningPotentialFourier;

	double* leftOutVectorNodeBuffer;
	double* rightOutVectorNodeBuffer;
	double* leftInVectorNodeBuffer;
	double* rightInVectorNodeBuffer;

	double* leftOutVectorCellBuffer;
	double* rightOutVectorCellBuffer;
	double* leftInVectorCellBuffer;
	double* rightInVectorCellBuffer;

	double* leftOutGmresBuffer;
	double* rightOutGmresBuffer;
	double* leftInGmresBuffer;
	double* rightInGmresBuffer;

	std::vector<Particle*> particles;
	std::vector<Particle*> escapedParticlesLeft;
	std::vector<Particle*> escapedParticlesRight;

	double*** tempCellParameterLeft;
	double*** tempCellParameterRight;
	double*** tempNodeParameterLeft;
	double*** tempNodeParameterRight;

	Vector3d*** tempCellVectorParameterLeft;
	Vector3d*** tempCellVectorParameterRight;
	Vector3d*** tempNodeVectorParameterLeft;
	Vector3d*** tempNodeVectorParameterRight;

	Matrix3d*** tempCellMatrixParameterLeft;
	Matrix3d*** tempCellMatrixParameterRight;
	Matrix3d*** tempNodeMatrixParameterLeft;
	Matrix3d*** tempNodeMatrixParameterRight;

	Matrix3d Kronecker;
	int LeviCivita[3][3][3];

	FILE* particleTypesFile;
	FILE* incrementFile;
	FILE* protonTraectoryFile;
	FILE* electronTraectoryFile;
	FILE* distributionFileProton;
	FILE* distributionFileElectron;
	FILE* distributionFileAlpha;
	FILE* distributionFilePositron;
	FILE* distributionFileHelium3;
	FILE* distributionFileDeuterium;
	FILE* distributionFileOxygen;
	FILE* distributionFileSilicon;
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
	FILE* generalAnisotropyFile;
	FILE* densityFile;
	FILE* divergenceErrorFile;
	FILE* velocityFile;
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
	FILE* rotEFile;
	FILE* EderivativeFile;

	FILE* maxwellMatrixFile;

	FILE* backupGeneralFile;
	FILE* backupParticlesFile;
	FILE* backupEfieldFile;
	FILE* backupBfieldFile;

    int protonNumber;
    int protonNumber1;
    int protonNumber2;
    int protonNumber3;
    int protonNumber4;
    int protonNumber5;
    int protonNumber6;
    int protonNumber7;
    int protonNumber8;
    int protonNumber9;
    int electronNumber;
    int electronNumber1;
    int electronNumber2;
    int electronNumber3;
    int electronNumber4;
    int electronNumber5;
    int electronNumber6;
    int electronNumber7;
    int electronNumber8;
    int electronNumber9;

	//FILE* outputEverythingFile;

	FILE* errorLogFile;

	//Simulation();
	Simulation();
	void setSpaceForProc();
	Simulation(int xn, int yn, int zn, double xsizev, double ysizev, double zsizev, double temp, double Vx,
                   double Vy, double Vz, double Ex, double Ey, double Ez, double Bx, double By, double Bz,
                   int maxIterations, double maxTimeV, int typesNumberV, int *particlesperBin,
                   double *concentrations, int inputType, int verbosityV);
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
	void initializeAnisotropicSilicon();
	void initializeWeibel();
	void initializeRingWeibel();
	void initializeHomogenouseFlow();
	void initializeKolmogorovSpectrum(int start, int end, double turbulenceFraction);
	void synchronizeParticleNumber();
	void createArrays();
	void createParticleTypes(double* concentrations, int* particlesPerBin);
    void evaluateParticleTypesAlpha();
	void createFiles();
	void simulate();
	void output();
	void outputTrajectories();
	void outputBackup();
	void rescaleConstants();

	void rescaleConstantsToTheoretical();
	void checkDebyeParameter();
	void checkGyroRadius();

	Matrix3d evaluateAlphaRotationTensor(double beta, Vector3d& velocity, double& gamma, Vector3d& EField, Vector3d& BField);
	void updateDeltaT();
	void evaluateElectricField();
	void updateEfield();
	void updateBfield();
	void updateShockWaveX();
	void evaluateExplicitDerivative();
	void checkEquationMatrix(std::vector<MatrixElement>**** matrix, int lnumber);
	void createSuperConductorLeftEquation(int i, int j, int k);
	void createFreeRightEquation(int i, int j, int k);
	void createFreeRightEquationX(int j, int k, Vector3d& rightPart);
	void createFreeRightEquationY(int j, int k, Vector3d& rightPart);
	void createFreeRightEquationZ(int j, int k, Vector3d& rightPart);
	void createLeftFakeEquation(int i, int j, int k);
	void createRightFakeEquation(int i, int j, int k);
	void createInternalEquationX(int i, int j, int k);
	void createInternalEquationY(int i, int j, int k);
	void createInternalEquationZ(int i, int j, int k);
	void createInternalEquation(int i, int j, int k);

	void smoothChargeDensityHat();
	void smoothChargeDensity();
	void smoothTempEfield();
	void smoothNewEfield();
	void smoothBfield();
	void smoothNewBfield();


	void evaluateMaxwellEquationMatrix();
	void evaluateMagneticField();
	void updateBoundaries();
	void updateBoundariesOldField();
	void updateBoundariesNewField();

	void substractMeanDensity();
	void cleanupDivergence();
	void cleanupDivergence1d();
	void substractMeanEfield();
	void updateFieldByCleaning();
	void evaluateDivergenceCleaningField();
	void createDivergenceCleanupInternalEquation(int i, int j, int k);

	void createDivergenceCleanupLeftFakeEquation(int i, int j, int k);
	void createDivergenceCleanupRightFakeEquation(int i, int j, int k);
	void createDivergenceCleanupSuperConductorEquation(int i, int j, int k);
	double cleanUpRightPart(int i, int j, int k);

	void exchangeEfield();

	void exchangeGeneralEfield(Vector3d*** field, Vector3d*** additionalFieldLeft, Vector3d*** additionalFieldRight);
	void exchangeGeneralBfield(Vector3d*** field, Vector3d*** additionalFieldLeft, Vector3d*** additionalFieldRight);

	Complex*** evaluateFourierTranslation(double*** a);

	double*** evaluateReverceFourierTranslation(Complex*** a);
	void resetNewTempFields();
	double volumeE(int i, int j, int k);
	double volumeB(int i, int j, int k);
	void checkParticleInBox(Particle& particle);

	void updateElectroMagneticParameters();
	void smoothDensity();
	void updateDensityParameters();

	void updateEnergy();
	void updateFields();
	void updateParameters();
	void updateAnisotropy();
	void updateExternalFlux();
	Vector3d evaluateRotB(int i, int j, int k);
	Vector3d evaluateRotTempE(int i, int j, int k);
	Vector3d evaluateRotE(int i, int j, int k);
	Vector3d evaluateRotNewE(int i, int j, int k);
	double evaluateDivE(int i, int j, int k);
	double evaluateDivCleaningE(int i, int j, int k);
	double evaluateDivTempE(int i, int j, int k);
	double evaluateDivNewE(int i, int j, int k);
	double evaluateDivFlux(int i, int j, int k);

	Vector3d getElectricFlux(int i, int j, int k);

	//Vector3d evaluateDivPressureTensor(int i, int j, int k);
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
    int getParticleNumber(int n, const ParticleTypes& type);
	int getTypeNumber(Particle* particle);

	Particle* getElectron(int n);
	Particle* createParticle(int n, int i, int j, int k, const double& weight, const ParticleTypes& type, const ParticleTypeContainer& typeContainer, const double& temperatureX, const double& temperatureY, const double& temperatureZ);
	void moveParticles();
	void removeEscapedParticles();
	void moveParticle(Particle* particle);
	void correctParticlePosition(Particle* particle);
	void correctParticlePosition(Particle& particle);
	void exchangeParticles();


	void evaluateParticlesRotationTensor();
	void injectNewParticles(int count, ParticleTypeContainer typeContainer, double length);

	void updateParticleCorrelationMaps();
	void updateCorrelationMaps(Particle* particle);
	void updateCorrelationMapCell(Particle* particle);
	void updateCorrelationMapNode(Particle* particle);
	void updateCorrelationMaps(Particle& particle);
	void updateCorrelationMapsX(Particle& particle);
	void updateCorrelationMapsY(Particle& particle);
	void updateCorrelationMapsZ(Particle& particle);
	void updateCorrelationMapCell(Particle& particle);
	void updateCorrelationMapCellX(Particle& particle);
	void updateCorrelationMapCellY(Particle& particle);
	void updateCorrelationMapCellZ(Particle& particle);
	void updateCorrelationMapNode(Particle& particle);
	void updateCorrelationMapNodeX(Particle& particle);
	void updateCorrelationMapNodeY(Particle& particle);
	void updateCorrelationMapNodeZ(Particle& particle);

	Vector3d correlationTempEfield(Particle* particle);
	Vector3d correlationNewEfield(Particle* particle);
	Vector3d correlationBfield(Particle* particle) const;
	Vector3d correlationEfield(Particle* particle);
	Vector3d correlationTempEfield(Particle& particle);
	Vector3d correlationEfield(Particle& particle);
	Vector3d correlationGeneralEfield(Particle& particle, Vector3d*** field, Vector3d*** additionalFieldLeft, Vector3d*** additionalFieldRight);
	Vector3d correlationNewEfield(Particle& particle);

	Vector3d correlationBfield(Particle& particle) const;

	double correlationWithBbin(Particle& particle, int i, int j, int k);

	double correlationWithEbin(Particle& particle, int i, int j, int k);

	double correlationBspline(const double& x, const double& dx, const double& leftx, const double& rightx);

	void eraseEscapedPaticles();

	void splitParticles();
	void splitParticle(Particle* particle);


	void sumCellParameters();
	void sumCellMatrixParameters();
	void sumNodeVectorParameters();
	void sumNodeMatrixParameters();
    void sumChargeDensityHat();

    void sumCellVectorParameters();

	void sumCellTempVectorParameters(Vector3d*** array);
	void sumCellTempParameters(double*** array);
	void sumCellTempMatrixParameters(Matrix3d*** array);

	void sumTempNodeParameters(double*** array);
	void sumTempNodeVectorParameters(Vector3d*** array);
	void sumTempNodeMatrixParameters(Matrix3d*** array);
};

#endif

