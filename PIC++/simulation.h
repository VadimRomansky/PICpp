#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <mpi.h>
#include <string>
#include "stdlib.h"
#include "constants.h"
#include "particleTypes.h"
#include "stdio.h"
#include "vector"
#include "list"
#include "largeVectorBasis.h"
#include "massMatrix.h"
#include "boundaryFieldEvaluator.h"

class MatrixElement;
class Complex;

enum SolverType {EXPLICIT, IMPLICIT, BUNEMAN};

enum InputType {CGS, Theoretical};

enum BoundaryConditionType {PERIODIC, SUPER_CONDUCTOR_LEFT, FREE_BOTH, FREE_MIRROR_BOTH};

enum DimensionType {ONE_D, TWO_D_XY, TWO_D_XZ, THREE_D};

class Matrix3d;
class Vector3d;
class Particle;
class ParticleTypeContainer;

class Simulation {
private:
	bool arrayCreated;
public:
	static double massProton;
	static double massElectron;
	static double massAlpha;
	static double massDeuterium;
	static double massHelium3;
	static double massOxygen;
	static double massSilicon;

	//static double speed_of_light_normalized;
	//static double speed_of_light_normalized_sqr;
	static double kBoltzman_normalized;
	static double electron_charge_normalized;
	static DimensionType dimensionType;

	static int initialRandom;

	int rank;
	int cartCoord[MPI_dim];
	int cartDim[MPI_dim];
	MPI_Comm cartComm;
	MPI_Comm cartCommX;
	MPI_Comm cartCommY;
	MPI_Comm cartCommZ;
	MPI_Comm cartCommXY;
	MPI_Comm cartCommYZ;
	MPI_Comm cartCommXZ;
	int nprocs;
	int leftRank;
	int rightRank;
	int frontRank;
	int backRank;
	int bottomRank;
	int topRank;

	std::string outputDir;
	std::string inputDir;
	std::string reducedOutputDir;
	int xnumberGeneral;
	int ynumberGeneral;
	int znumberGeneral;
	int firstAbsoluteXindex;
	int firstAbsoluteYindex;
	int firstAbsoluteZindex;
	double leftX;
	double rightX;
	double leftY;
	double rightY;
	double leftZ;
	double rightZ;

	int xnumber;
	int ynumber;
	int znumber;
	int xnumberAdded;
	int ynumberAdded;
	int znumberAdded;
	int particlesNumber;

	int resistiveLayerWidth;
	double fakeCondactivity;

	double* concentrations;
	int* particlesPerBin;

	int innputType;
	double initialMagnetization;
	double initialElectronConcentration;

	double Btheta;
	double Bphi;

	double density;
	double temperature;

	double plasma_period;
	double plasma_period2;
	double scaleFactor;

	double time;
	double maxTime;
	int currentIteration;
	int maxIteration;
	int writeParameter;
	int writeGeneralParameter;
	int writeTrajectoryNumber;
	int writeParticleNumber;
	int currentWriteNumber;
	double smoothingParameter;
	int smoothingCount;
	bool multiplyFileOutput;
	double xsizeGeneral;
	double ysizeGeneral;
	double zsizeGeneral;
	double xsize;
	double ysize;
	double zsize;

	double deltaT;
	double preferedDeltaT;
	double electronMassInput;

	double deltaX;
	double deltaY;
	double deltaZ;

	double deltaX2;
	double deltaY2;
	double deltaZ2;

	double cellVolume;

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
	BoundaryConditionType boundaryConditionTypeX;
	BoundaryConditionType boundaryConditionTypeY;
	BoundaryConditionType boundaryConditionTypeZ;
	int maxwellEquationMatrixSize;

	double particleEnergy;
	double electricFieldEnergy;
	double electricFieldEnergyX;
	double electricFieldEnergyY;
	double electricFieldEnergyZ;
	double magneticFieldEnergy;
	double magneticFieldEnergyX;
	double magneticFieldEnergyY;
	double magneticFieldEnergyZ;
	double energy;

	double theoreticalEnergy;
	double generalTheoreticalEnergy;
	Vector3d theoreticalMomentum;
	Vector3d electromagneticMomentum;
	Vector3d particleMomentum;
	Vector3d generalTheoreticalMomentum;

	int chargeBalance;

	double extJ;

	Vector3d globalMomentum;

	ParticleTypeContainer* types;
	int typesNumber;

	double**** particleConcentrations;
	double**** particleEnergies;

	double*** chargeDensity;
	double*** chargeDensityMinus;

	Vector3d**** particleBulkVelocities;

	MassMatrix*** massMatrix;
	MassMatrix*** tempMassMatrix;

	Vector3d maxEfield;
	Vector3d maxBfield;
	double meanSquaredEfield[3];

	Vector3d V0;

	Vector3d B0;
	Vector3d E0;

	Vector3d rightEfield;
	Vector3d leftEfield;
	Vector3d rightBfield;
	Vector3d leftBfield;

	BoundaryFieldEvaluator* leftBoundaryFieldEvaluator;
	BoundaryFieldEvaluator* rightBoundaryFieldEvaluator;

	double omegaPlasmaProton;
	double omegaPlasmaElectron;
	double omegaPlasmaTotal;
	double omegaGyroProton;
	double omegaGyroElectron;

	//for turbulence
	double turbulenceAmplitude;
	double turbulenceFraction;
	int turbulenceRandomSeed;

	int minTurbulenceLengthX, maxTurbulenceLengthX;
	int minTurbulenceLengthY, maxTurbulenceLengthY;
	int minTurbulenceLengthZ, maxTurbulenceLengthZ;

	//for debug alven wave only/////
	double omegaAlfven;
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

	//std::vector<MatrixElement>****bunemanDivergenceCleanUpMatrix;
	//double**** bunemanDivergenceCleanUpRightPart;

	Vector3d*** electricFlux;
	Vector3d*** electricFluxMinus;

	Vector3d*** externalElectricFlux;

	double*** chargeDensityHat;
	Matrix3d*** dielectricTensor;
	Matrix3d*** pressureTensor;
	Vector3d*** divPressureTensor;
	Vector3d*** divPressureTensorMinus;

	Vector3d*** Efield;
	Vector3d*** Bfield;

	Vector3d*** newEfield;
	Vector3d*** newBfield;

	Vector3d*** tempEfield;

	//Vector3d*** tempBfield;

	Vector3d*** explicitEfield;
	Vector3d*** rotB;
	Vector3d*** Ederivative;

	Vector3d*** rotE;
	Vector3d*** Bderivative;

	double*** bunemanJx;
	double*** bunemanJy;
	double*** bunemanJz;

	double*** bunemanEx;
	double*** bunemanEy;
	double*** bunemanEz;

	double*** bunemanDivCleaningEx;
	double*** bunemanDivCleaningEy;
	double*** bunemanDivCleaningEz;

	double*** bunemanNewEx;
	double*** bunemanNewEy;
	double*** bunemanNewEz;

	double*** tempBunemanExParameter;
	double*** tempBunemanEyParameter;
	double*** tempBunemanEzParameter;

	double*** bunemanBx;
	double*** bunemanBy;
	double*** bunemanBz;

	double*** bunemanNewBx;
	double*** bunemanNewBy;
	double*** bunemanNewBz;

	double*** tempBunemanBxParameter;
	double*** tempBunemanByParameter;
	double*** tempBunemanBzParameter;

	double*** bunemanDivCleaningBx;
	double*** bunemanDivCleaningBy;
	double*** bunemanDivCleaningBz;

	double**** divergenceCleaningField;
	double**** divergenceCleaningPotential;
	double**** tempDivergenceCleaningPotential;

	double*** bunemanChargeDensity;
	double**** bunemanDivergenceCleaningPotential;

	double*** divergenceCleaningPotentialFourier;

	Complex*** fourierInput;
	Complex*** fourierImage;
	Complex*** fourierOutput;

	double* leftOutMaximumBuffer;
	double* rightOutMaximumBuffer;
	double* leftInMaximumBuffer;
	double* rightInMaximumBuffer;

	double* frontOutMaximumBuffer;
	double* backOutMaximumBuffer;
	double* frontInMaximumBuffer;
	double* backInMaximumBuffer;

	double* topOutMaximumBuffer;
	double* bottomOutMaximumBuffer;
	double* topInMaximumBuffer;
	double* bottomInMaximumBuffer;

	double* leftOutGmresBuffer;
	double* rightOutGmresBuffer;
	double* leftInGmresBuffer;
	double* rightInGmresBuffer;

	double* frontOutGmresBuffer;
	double* backOutGmresBuffer;
	double* frontInGmresBuffer;
	double* backInGmresBuffer;

	double* topOutGmresBuffer;
	double* bottomOutGmresBuffer;
	double* topInGmresBuffer;
	double* bottomInGmresBuffer;

	double* leftOutDivergenceBuffer;
	double* rightOutDivergenceBuffer;
	double* leftInDivergenceBuffer;
	double* rightInDivergenceBuffer;

	double* frontOutDivergenceBuffer;
	double* backOutDivergenceBuffer;
	double* frontInDivergenceBuffer;
	double* backInDivergenceBuffer;

	double* topOutDivergenceBuffer;
	double* bottomOutDivergenceBuffer;
	double* topInDivergenceBuffer;
	double* bottomInDivergenceBuffer;


	////buneman E
	double* leftOutBunemanExBuffer;
	double* rightOutBunemanExBuffer;
	double* leftInBunemanExBuffer;
	double* rightInBunemanExBuffer;

	double* frontOutBunemanExBuffer;
	double* backOutBunemanExBuffer;
	double* frontInBunemanExBuffer;
	double* backInBunemanExBuffer;

	double* bottomOutBunemanExBuffer;
	double* topOutBunemanExBuffer;
	double* bottomInBunemanExBuffer;
	double* topInBunemanExBuffer;

	double* leftOutBunemanEyBuffer;
	double* rightOutBunemanEyBuffer;
	double* leftInBunemanEyBuffer;
	double* rightInBunemanEyBuffer;

	double* frontOutBunemanEyBuffer;
	double* backOutBunemanEyBuffer;
	double* frontInBunemanEyBuffer;
	double* backInBunemanEyBuffer;

	double* bottomOutBunemanEyBuffer;
	double* topOutBunemanEyBuffer;
	double* bottomInBunemanEyBuffer;
	double* topInBunemanEyBuffer;

	double* leftOutBunemanEzBuffer;
	double* rightOutBunemanEzBuffer;
	double* leftInBunemanEzBuffer;
	double* rightInBunemanEzBuffer;

	double* frontOutBunemanEzBuffer;
	double* backOutBunemanEzBuffer;
	double* frontInBunemanEzBuffer;
	double* backInBunemanEzBuffer;

	double* bottomOutBunemanEzBuffer;
	double* topOutBunemanEzBuffer;
	double* bottomInBunemanEzBuffer;
	double* topInBunemanEzBuffer;

	///buneman B
	double* leftOutBunemanBxBuffer;
	double* rightOutBunemanBxBuffer;
	double* leftInBunemanBxBuffer;
	double* rightInBunemanBxBuffer;

	double* frontOutBunemanBxBuffer;
	double* backOutBunemanBxBuffer;
	double* frontInBunemanBxBuffer;
	double* backInBunemanBxBuffer;

	double* bottomOutBunemanBxBuffer;
	double* topOutBunemanBxBuffer;
	double* bottomInBunemanBxBuffer;
	double* topInBunemanBxBuffer;

	double* leftOutBunemanByBuffer;
	double* rightOutBunemanByBuffer;
	double* leftInBunemanByBuffer;
	double* rightInBunemanByBuffer;

	double* frontOutBunemanByBuffer;
	double* backOutBunemanByBuffer;
	double* frontInBunemanByBuffer;
	double* backInBunemanByBuffer;

	double* bottomOutBunemanByBuffer;
	double* topOutBunemanByBuffer;
	double* bottomInBunemanByBuffer;
	double* topInBunemanByBuffer;

	double* leftOutBunemanBzBuffer;
	double* rightOutBunemanBzBuffer;
	double* leftInBunemanBzBuffer;
	double* rightInBunemanBzBuffer;

	double* frontOutBunemanBzBuffer;
	double* backOutBunemanBzBuffer;
	double* frontInBunemanBzBuffer;
	double* backInBunemanBzBuffer;

	double* bottomOutBunemanBzBuffer;
	double* topOutBunemanBzBuffer;
	double* bottomInBunemanBzBuffer;
	double* topInBunemanBzBuffer;

	std::vector<Particle> fakeParticles; //only to know possible size;
	std::vector<Particle*> particles;
	std::vector<Particle*> tempParticles;

	std::vector<Particle*> escapedParticlesLeft;
	std::vector<Particle*> escapedParticlesRight;
	std::vector<Particle*> escapedParticlesFront;
	std::vector<Particle*> escapedParticlesBack;
	std::vector<Particle*> escapedParticlesTop;
	std::vector<Particle*> escapedParticlesBottom;

	std::vector<Particle*> reservedParticles;

	std::list<std::pair<int, double> >* mostAcceleratedParticlesNumbers;

	int trackedParticlesNumber;
	int** trackedParticlesNumbers;

	double*** tempCellParameter;
	double*** tempCellParameterLeft;
	double*** tempCellParameterRight;
	double*** tempCellParameterFront;
	double*** tempCellParameterBack;
	double*** tempCellParameterBottom;
	double*** tempCellParameterTop;

	double*** tempNodeParameter;
	double*** tempNodeParameterLeft;
	double*** tempNodeParameterRight;
	double*** tempNodeParameterFront;
	double*** tempNodeParameterBack;
	double*** tempNodeParameterBottom;
	double*** tempNodeParameterTop;

	double*** tempBunemanJxLeft;
	double*** tempBunemanJxRight;
	double*** tempBunemanJxFront;
	double*** tempBunemanJxBack;
	double*** tempBunemanJxBottom;
	double*** tempBunemanJxTop;

	double*** tempBunemanJyLeft;
	double*** tempBunemanJyRight;
	double*** tempBunemanJyFront;
	double*** tempBunemanJyBack;
	double*** tempBunemanJyBottom;
	double*** tempBunemanJyTop;

	double*** tempBunemanJzLeft;
	double*** tempBunemanJzRight;
	double*** tempBunemanJzFront;
	double*** tempBunemanJzBack;
	double*** tempBunemanJzBottom;
	double*** tempBunemanJzTop;

	Vector3d*** tempCellVectorParameter;
	Vector3d*** tempCellVectorParameterLeft;
	Vector3d*** tempCellVectorParameterRight;
	Vector3d*** tempCellVectorParameterFront;
	Vector3d*** tempCellVectorParameterBack;
	Vector3d*** tempCellVectorParameterBottom;
	Vector3d*** tempCellVectorParameterTop;

	Vector3d*** tempNodeVectorParameter;
	Vector3d*** tempNodeVectorParameterLeft;
	Vector3d*** tempNodeVectorParameterRight;
	Vector3d*** tempNodeVectorParameterFront;
	Vector3d*** tempNodeVectorParameterBack;
	Vector3d*** tempNodeVectorParameterBottom;
	Vector3d*** tempNodeVectorParameterTop;

	Matrix3d*** tempCellMatrixParameter;
	Matrix3d*** tempCellMatrixParameterLeft;
	Matrix3d*** tempCellMatrixParameterRight;
	Matrix3d*** tempCellMatrixParameterFront;
	Matrix3d*** tempCellMatrixParameterBack;
	Matrix3d*** tempCellMatrixParameterBottom;
	Matrix3d*** tempCellMatrixParameterTop;

	Matrix3d*** tempNodeMatrixParameter;
	Matrix3d*** tempNodeMatrixParameterLeft;
	Matrix3d*** tempNodeMatrixParameterRight;
	Matrix3d*** tempNodeMatrixParameterFront;
	Matrix3d*** tempNodeMatrixParameterBack;
	Matrix3d*** tempNodeMatrixParameterBottom;
	Matrix3d*** tempNodeMatrixParameterTop;

	MassMatrix*** tempNodeMassMatrixParameterLeft;
	MassMatrix*** tempNodeMassMatrixParameterRight;
	MassMatrix*** tempNodeMassMatrixParameterFront;
	MassMatrix*** tempNodeMassMatrixParameterBack;
	MassMatrix*** tempNodeMassMatrixParameterBottom;
	MassMatrix*** tempNodeMassMatrixParameterTop;

	double**** residualBiconjugateDivE;
	double**** firstResidualBiconjugateDivE;
	double**** vBiconjugateDivE;
	double**** pBiconjugateDivE;
	double**** sBiconjugateDivE;
	double**** tBiconjugateDivE;

	double**** residualBiconjugateMaxwell;
	double**** firstResidualBiconjugateMaxwell;
	double**** vBiconjugateMaxwell;
	double**** pBiconjugateMaxwell;
	double**** sBiconjugateMaxwell;
	double**** tBiconjugateMaxwell;

	Complex*** fourierScalarInput;
	Complex*** fourierScalarOutput;
	Complex*** fourierScalarTempOutput;
	Complex*** fourierScalarTempOutput1;

	Complex*** fourierScalarMirrorInput;
	Complex*** fourierScalarMirrorOutput;
	Complex*** fourierScalarMirrorTempOutput;
	Complex*** fourierScalarMirrorTempOutput1;

	Complex* localFactorX;
	Complex* localFactorY;
	Complex* localFactorZ;

	int derExPoint;
	int derConcentrationPoint;
	Vector3d** leftElevel;
	Vector3d** rightElevel;

	Vector3d* rightMeanElevel;

	int constMeanElevelPoint;

	Matrix3d Kronecker;
	int LeviCivita[3][3][3];

	FILE* particleTypesFile;
	FILE* incrementFile;
	FILE* particlesTrajectoryFile;
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

	//FILE* outputEverythingFile;

	FILE* errorLogFile;

	//Simulation();
	Simulation();
	void setSpaceForProc();
	Simulation(int xn, int yn, int zn, double dxv, double temp, double Vx,
                   double Vy, double Vz, double sigmav, double Bthetav, double Bphiv, double E0x, double E0y, double E0z, double electronInitConc, 
                   int maxIterations, double maxTimeV, int writeIterationV, int writeGeneralV, int writeTrajectoryV, int writeParticleV, int smoothingCountV, double smoothingParameterV, bool multiplyFileOutputV, int typesNumberV, int *particlesperBin,
                   double *concentrations, int inputType, int nprocsV, int verbosityV, double preferedTimeStepV, double massElectronInputV, MPI_Comm& comm);
	Simulation(int xn, int yn, int zn, double dxv, double temp, double Vx,
                   double Vy, double Vz, double sigmav, double Bthetav, double Bphiv, double E0x, double E0y, double E0z, double electronInitConc, 
                   int maxIterations, double maxTimeV, int writeIterationV, int writeGeneralV, int writeTrajectoryV, int writeParticleV, int smoothingCountV, double smoothingParameterV, bool multiplyFileOutputV, int typesNumberV, int *particlesperBin,
                   double *concentrations, int inputType, int nprocsV, int verbosityV, double preferedTimeStepV, double massElectronInputV, double plasmaPeriodV,
                   double scaleFactorV, SolverType solverTypev, MPI_Comm& comm);
	/*Simulation(int xn, int yn, int zn, double xsizev, double ysizev, double zsizev, double temp, double Vx,
                   double Vy, double Vz, double Ex, double Ey, double Ez, double Bx, double By, double Bz,
                   int maxIterations, double maxTimeV, int typesNumberV, int *particlesperBin,
                   double *concentrations, int inputType, int nprocsV, int verbosityV, double preferedTimeStepV, double massElectronInputV, MPI_Comm& comm);
	Simulation(int xn, int yn, int zn, double xsizev, double ysizev, double zsizev, double temp, double Vx,
                   double Vy, double Vz, double Ex, double Ey, double Ez, double Bx, double By, double Bz,
                   int maxIterations, double maxTimeV, int typesNumberV, int *particlesperBin,
                   double *concentrations, int inputType, int nprocsV, int verbosityV, double preferedTimeStepV, double massElectronInputV, double plasmaPeriodV,
                   double scaleFactorV, SolverType solverTypev,  MPI_Comm& comm);*/
	~Simulation();

	void initialize();
	void initializeSimpleElectroMagneticWave();
	void initializeSimpleElectroMagneticWaveY();
	void initializeSimpleElectroMagneticWaveZ();
	void initializeRotatedSimpleElectroMagneticWave(int waveCountX, int waveCountY, int waveCountZ);
	void checkFrequency(double omega);
	void initializeAlfvenWaveX(int wavesCount, double amplitudeRelation);
	void initializeAlfvenWaveY(int wavesCount, double amplitudeRelation);
	void initializeAlfvenWaveZ(int wavesCount, double amplitudeRelation);
	void initializeRotatedAlfvenWave(int waveCountX, int waveCountY, int waveCountZ, double amplitudeRelation);
	void initializeLangmuirWave();
	void initializeTwoStream();
	void initializeExternalFluxInstability();
	void initializeFluxFromRight();
	void fieldsLorentzTransitionX(const double& v);
	void initializeShockWave();
	void initializeBell();
	void solveRankineHugoniot(double upstreamDensity, Vector3d upstreamVelocity, double upstreamPressure, Vector3d upstreamB, Vector3d upstreamE, double& downstreamDensity, Vector3d& downstreamVelocity, double& downstreamPressure, Vector3d& downstreamB, Vector3d& downstreamE, double adiabaticParameter, double compressionRatio);
	double evaluateTemperatureByPressure(double& pressure);
	double evaluatePressureByTemperature(double& temperature);
	void initializeAnisotropic();
	void initializeAnisotropicSilicon();
	void initializeWeibel();
	void initializeRingWeibel();
	void initializeHomogenouseFlow();
	void initializeKolmogorovSpectrum();
	double evaluateTurbulentB(int ki, int kj, int kk);
	void initializeRandomModes(int number, int minNumber, double energyFraction);
	void initializeFake();
	void initializeHarris();
	void initializeTestOneParticle();
	void synchronizeParticleNumber();
	void createArrays();
	void createParticleTypes(double* concentrations, int* particlesPerBin);
	void evaluateJuttnerFunctions();
	void evaluateJuttnerFunction(ParticleTypeContainer& typeContainer);
    void evaluateParticleTypesAlpha();
	void createFiles();
	void tristanEvaluateBhalfStep();
	void tristanEvaluateE();
	void addParticleFluxZigzag(Particle* particle);
	void tristanUpdateFlux();
	void injectNewParticles();
	void simulate();
	void output();
	void resetBunemanFieldToCellVectorParameter(double*** bunemanEx, double*** bunemanEy, double*** bunemanEz);
	void outputTrajectories();
	void outputBackup();
	void rescaleConstants();

	void rescaleConstantsToTheoretical();
	void checkDebyeParameter();
	void checkGyroRadius();

	int getCartCoordWithAbsoluteIndexX(int i);
	int getCartCoordWithAbsoluteIndexY(int j);
	int getCartCoordWithAbsoluteIndexZ(int k);

	int getLocalIndexByAbsoluteX(int i);
	int getLocalIndexByAbsoluteY(int j);
	int getLocalIndexByAbsoluteZ(int k);

	Matrix3d evaluateAlphaRotationTensor(const double& beta, Vector3d& velocity, double& gamma, Vector3d& EField, Vector3d& BField);
	void updateDeltaT();
	void evaluateElectricField();
	void updateEfield();
	void updateBfield();
	void updateShockWaveX();
	void evaluateExplicitDerivative();
	void checkEquationMatrix(std::vector<MatrixElement>**** matrix, int lnumber);
	void createSuperConductorLeftEquation(int i, int j, int k);
	void createFreeRightEquation(int i, int j, int k);
	void createFreeLeftEquation(int i, int j, int k);
	void createFakeEquation(int i, int j, int k);
	void createInternalEquationX(int i, int j, int k);
	void createInternalEquationY(int i, int j, int k);
	void createInternalEquationZ(int i, int j, int k);
	void createInternalECEquationX(int i, int j, int k);
	void createInternalECEquationY(int i, int j, int k);
	void createInternalECEquationZ(int i, int j, int k);
	void createInternalEquation(int i, int j, int k);
	bool isInResistiveLayer(int i, int j, int k);

	void createSmoothingMatrix(double smoothingmatrix[3][3][3]);

	void smoothChargeDensityHat();
	void smoothVectorNodeParameter(Vector3d*** E);
	void smoothMatrixNodeParameter(Matrix3d*** E);
	void smoothCellParameter(double*** array);
	void smoothChargeDensity();
	void smoothTempEfield();
	void smoothNewEfield();
	void smoothFlux();
	void smoothVectorCellParameter(Vector3d*** B);
	void smoothBfield();
	void smoothNewBfield();

	void smoothBunemanEfieldGeneral(double*** fieldX, double*** fieldY, double*** fieldZ);
	void smoothBunemanBfieldGeneral(double*** fieldX, double*** fieldY, double*** fieldZ);

	void fuckingStrangeSmoothingBunemanFields(double*** oldEx, double*** oldEy, double*** oldEz, double*** oldBx, double*** oldBy, double*** oldBz, double*** newEx,
	                                         double*** newEy, double*** newEz, double*** newBx, double*** newBy, double*** newBz);


	void evaluateMaxwellEquationMatrix();
	void evaluateMagneticField();

	void cleanupDivergence(Vector3d*** field, double*** density);
	void cleanupDivergence1d(Vector3d*** field, double*** density);
	void substractMeanEfield(double**** field);
	void updateFieldByCleaning(Vector3d*** field);
	void evaluateDivergenceCleaningField();
	void createDivergenceCleanupInternalEquation(int i, int j, int k, Vector3d*** field, double*** density);

	void createDivergenceFakeEquation(int i, int j, int k);
	void createDivergenceZeroEquation(int i, int j, int k);
	void createDivergenceFixEquation(int i, int j, int k, double*** density);
	double cleanUpRightPart(int i, int j, int k, Vector3d*** field, double*** density);
	double evaluateMeanChargeDensity();
	void substractMeanChargeDensity();

	void cleanupDivergenceMagnetic();
	void updateFieldByCleaningMagnetic();
	void evaluateDivergenceCleaningFieldMagnetic();
	void createDivergenceCleanupInternalEquationMagnetic(int i, int j, int k);
	double cleanUpRightPartMagnetic(int i, int j, int k);

	void cleanupDivergenceBuneman();
	void updateFieldByCleaningBuneman();
	void evaluateDivergenceCleaningFieldBuneman();
	void createDivergenceCleanupInternalEquationBuneman(int i, int j, int k);
	double cleanUpRightPartBuneman(int i, int j, int k);

	void cleanupDivergenceBunemanMagnetic();
	void updateFieldByCleaningBunemanMagnetic();
	void evaluateDivergenceCleaningFieldBunemanMagnetic();
	void createDivergenceCleanupInternalEquationBunemanMagnetic(int i, int j, int k);
	double cleanUpRightPartBunemanMagnetic(int i, int j, int k);

	void exchangeEfield();

	void exchangeGeneralEfield(Vector3d*** field);
	void exchangeGeneralEfieldX(Vector3d*** field);
	void exchangeGeneralEfieldY(Vector3d*** field);
	void exchangeGeneralEfieldZ(Vector3d*** field);

	void exchangeGeneralBfield(Vector3d*** field);
	void exchangeGeneralBfieldX(Vector3d*** field);
	void exchangeGeneralBfieldY(Vector3d*** field);
	void exchangeGeneralBfieldZ(Vector3d*** field);

	void exchangeGeneralScalarCellField(double**** field);
	void exchangeGeneralScalarCellFieldX(double**** field);
	void exchangeGeneralScalarCellFieldY(double**** field);
	void exchangeGeneralScalarCellFieldZ(double**** field);

	void exchangeGeneralScalarCellField(double*** field);
	void exchangeGeneralScalarCellFieldX(double*** field);
	void exchangeGeneralScalarCellFieldY(double*** field);
	void exchangeGeneralScalarCellFieldZ(double*** field);

	void exchangeGeneralScalarNodeField(double**** field);
	void exchangeGeneralScalarNodeFieldX(double**** field);
	void exchangeGeneralScalarNodeFieldY(double**** field);
	void exchangeGeneralScalarNodeFieldZ(double**** field);

	void exchangeGeneralMatrixNodeField(Matrix3d*** field);
	void exchangeGeneralMatrixNodeFieldX(Matrix3d*** field);
	void exchangeGeneralMatrixNodeFieldY(Matrix3d*** field);
	void exchangeGeneralMatrixNodeFieldZ(Matrix3d*** field);

	void exchangeBunemanEfield(double*** fieldX, double*** fieldY, double*** fieldZ);

	void exchangeBunemanExAlongX(double*** fieldX);
	void exchangeBunemanExAlongY(double*** fieldX);
	void exchangeBunemanExAlongZ(double*** fieldX);

	void exchangeBunemanEyAlongX(double*** fieldY);
	void exchangeBunemanEyAlongY(double*** fieldY);
	void exchangeBunemanEyAlongZ(double*** fieldY);

	void exchangeBunemanEzAlongX(double*** fieldZ);
	void exchangeBunemanEzAlongY(double*** fieldZ);
	void exchangeBunemanEzAlongZ(double*** fieldZ);

	void exchangeBunemanBfield(double*** fieldX, double*** fieldY, double*** fieldZ);

	void exchangeBunemanBxAlongX(double*** fieldX);
	void exchangeBunemanBxAlongY(double*** fieldX);
	void exchangeBunemanBxAlongZ(double*** fieldX);

	void exchangeBunemanByAlongX(double*** fieldY);
	void exchangeBunemanByAlongY(double*** fieldY);
	void exchangeBunemanByAlongZ(double*** fieldY);

	void exchangeBunemanBzAlongX(double*** fieldZ);
	void exchangeBunemanBzAlongY(double*** fieldZ);
	void exchangeBunemanBzAlongZ(double*** fieldZ);

	Complex*** evaluateFourierTranslation(double*** a);

	double*** evaluateReverceFourierTranslation(Complex*** a);
	void resetNewTempFields();
	double volumeE();
	double volumeB();
	void checkParticleInBox(Particle& particle);

	void updateElectroMagneticParameters();
	void updateDensityParameters();

	void updateEnergy();
	void updateTheoreticalEnergy();
	void updateFields();
	void updateParameters();
	void updateAnisotropy();
	void updateExternalFlux();
	Vector3d evaluateRotB(int i, int j, int k);
	Vector3d evaluateRotEgeneral(Vector3d*** E, int i, int j, int k);
	Vector3d evaluateRotTempE(int i, int j, int k);
	Vector3d evaluateRotE(int i, int j, int k);
	Vector3d evaluateRotNewE(int i, int j, int k);
	double evaluateDivEgeneral(Vector3d*** E, int i, int j, int k);
	double evaluateDivE(int i, int j, int k);
	double evaluateDivCleaningE(int i, int j, int k);
	double evaluateDivTempE(int i, int j, int k);
	double evaluateDivNewE(int i, int j, int k);
	double evaluateDivFlux(int i, int j, int k);
	Vector3d evaluateDivPressureTensor(int i, int j, int k);
	double evaluateDivB(int i, int j, int k);
	double evaluateDivNewB(int i, int j, int k);
	double evaluateDivBgeneral(Vector3d*** B, int i, int j, int k);

	double evaluateDivBunemanEgeneral(double*** fieldX, double*** fieldY, double*** fieldZ, int i, int j, int k);
	double evaluateDivBunemanNewE(int i, int j, int k);
	double evaluateDivBunemanE(int i, int j, int k);

	double evaluateDivBunemanBgeneral(double*** fieldX, double*** fieldY, double*** fieldZ, int i, int j, int k);
	double evaluateDivBunemanNewB(int i, int j, int k);
	double evaluateDivBunemanB(int i, int j, int k);

	double evaluateBunemanRotEx(int i, int j, int k);
	double evaluateBunemanRotEy(int i, int j, int k);
	double evaluateBunemanRotEz(int i, int j, int k);

	double evaluateBunemanRotBx(int i, int j, int k);
	double evaluateBunemanRotBy(int i, int j, int k);
	double evaluateBunemanRotBz(int i, int j, int k);

	Vector3d getBunemanElectricField(int i, int j, int k);
	Vector3d getBunemanMagneticField(int i, int j, int k);

	Vector3d getBunemanFlux(int i, int j, int k);


	void interpolateBunemanToLapentaEfield(double*** Ex, double*** Ey, double*** Ez, Vector3d*** E);
	void interpolateLapentaToBunemanEfield(double*** Ex, double*** Ey, double*** Ez, Vector3d*** E);

	void interpolateBunemanToLapentaBfield(double*** Bx, double*** By, double*** Bz, Vector3d*** B);
	void interpolateLapentaToBunemanBfield(double*** Bx, double*** By, double*** Bz, Vector3d*** B);

	double evaluateTurbulenceFieldAmplitude(const double& kx, const double& ky, const double& kz);

	void updateBunemanFields();
	void updateBunemanElectricField();
	void updateBunemanMagneticField();

	void updateBunemanChargeDensity();

	void filterFields(int cutWaveNumber);
	void filterFieldGeneral(Vector3d*** field, int cutWaveNumber);
	void filterFieldGeneralRight(Vector3d*** field, int cutWaveNumber, int startIndex);
	void filterFieldGeneralLeft(Vector3d*** field, int cutWaveNumber, int startIndex);
	void filterFieldGeneralRightMirror(Vector3d*** field, int cutWaveNumber, int startIndex);

	void filterFieldsLocal(int cutWaveNumber);
	void filterFieldGeneralLocal(Vector3d*** field, int cutWaveNumber);

	void updateMaxEderivativePoint();
	void updateMaxConcentrationDerivativePoint();
	void updateBoundaryLevelX();
	void substractStep(Vector3d*** field, Vector3d** left, Vector3d** right, int sign);
	void updateMeanLevel(Vector3d*** field);
	void updateLastMeanLevelPoint(double relativeError, int minValue);

	Vector3d averageFieldXY(Vector3d*** field, int i);
	Vector3d averageFieldXZ(Vector3d*** field, int i);
	Vector3d averageFieldYZ(Vector3d*** field, int i);

	double averageConcentrationYZ(double*** concentration, int i);


	Vector3d getElectricFlux(int i, int j, int k);

	//Vector3d evaluateDivPressureTensor(int i, int j, int k);
	Vector3d evaluateGradDensity(int i, int j, int k);
	void createParticles();
	void createParticlesHarris(double harrisWidth);
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
	void sortParticleToEscaped(Particle* particle);
	void removeEscapedParticles();
	void moveParticle(Particle* particle);
	void moveParticle(Particle* particle, int cur, int N);
	
	void moveParticleTristan(Particle* particle);
	void exchangeParticles();
	void collectMostAcceleratedParticles();


	void evaluateParticlesRotationTensor();
	void updateParticlesBeta();
	void injectNewParticles(int count, ParticleTypeContainer& typeContainer, const double& length);

	void updateParticleCorrelationMaps();
	void updateCorrelationMaps(Particle* particle);
	void updateCorrelationMapCell(Particle* particle);
	void updateCorrelationMapNode(Particle* particle);
	void updateCorrelationMapsX(Particle* particle);
	void updateCorrelationMapsY(Particle* particle);
	void updateCorrelationMapsZ(Particle* particle);
	void updateCorrelationMapCellX(Particle* particle);
	void updateCorrelationMapCellY(Particle* particle);
	void updateCorrelationMapCellZ(Particle* particle);
	void updateCorrelationMapNodeX(Particle* particle);
	void updateCorrelationMapNodeY(Particle* particle);
	void updateCorrelationMapNodeZ(Particle* particle);

	void updateBunemanCorrelationMaps(Particle* particle);
	void updateBunemanCorrelationMapCell(Particle* particle);
	void updateBunemanCorrelationMapNode(Particle* particle);
	void updateBunemanCorrelationMapsX(Particle* particle);
	void updateBunemanCorrelationMapsY(Particle* particle);
	void updateBunemanCorrelationMapsZ(Particle* particle);
	void updateBunemanCorrelationMapCellX(Particle* particle);
	void updateBunemanCorrelationMapCellY(Particle* particle);
	void updateBunemanCorrelationMapCellZ(Particle* particle);
	void updateBunemanCorrelationMapNodeX(Particle* particle);
	void updateBunemanCorrelationMapNodeY(Particle* particle);
	void updateBunemanCorrelationMapNodeZ(Particle* particle);

	Vector3d correlationTempEfield(Particle* particle);
	Vector3d correlationNewEfield(Particle* particle);
	Vector3d correlationBfield(Particle* particle) const;
	Vector3d correlationNewBfield(Particle* particle) const;
	Vector3d correlationEfield(Particle* particle);
	Vector3d correlationGeneralEfield(Particle* particle, Vector3d*** field);



	Vector3d correlationBunemanBfield(Particle* particle);
	Vector3d correlationBunemanNewBfield(Particle* particle);
	Vector3d correlationBunemanGeneralBfield(Particle* particle, double*** fieldX, double*** fieldY, double*** fieldZ);

	Vector3d correlationBunemanEfield(Particle* particle);
	Vector3d correlationBunemanNewEfield(Particle* particle);
	Vector3d correlationBunemanGeneralEfield(Particle* particle, double*** fieldX, double*** fieldY, double*** fieldZ);

	void correlationBunemanEBfields(Particle* particle, double*** Ex, double*** Ey, double*** Ez, double*** Bx, double*** By, double*** Bz, Vector3d& E, Vector3d& B);
	void correlationBunemanEBfieldsWithoutMaps(Particle* particle, double*** Ex, double*** Ey, double*** Ez, double*** Bx,
	                                           double*** By, double*** Bz, Vector3d& E, Vector3d& B);

	double correlationWithBbin(Particle& particle, int i, int j, int k);

	double correlationWithEbin(Particle& particle, int i, int j, int k);

	double correlationBspline(const double& x, const double& dx, const double& leftx, const double& rightx);

	void eraseEscapedPaticles();

	void splitParticles();
	void splitParticle(Particle* particle);

	void exchangeBunemanFlux();

	void sumBunemanJxAlongX();
	void sumBunemanJxAlongY();
	void sumBunemanJxAlongZ();

	void sumBunemanJyAlongX();
	void sumBunemanJyAlongY();
	void sumBunemanJyAlongZ();

	void sumBunemanJzAlongX();
	void sumBunemanJzAlongY();
	void sumBunemanJzAlongZ();

	void sumTempBunemanJxAlongX();
	void sumTempBunemanJxAlongY();
	void sumTempBunemanJxAlongZ();

	void sumTempBunemanJyAlongX();
	void sumTempBunemanJyAlongY();
	void sumTempBunemanJyAlongZ();

	void sumTempBunemanJzAlongX();
	void sumTempBunemanJzAlongY();
	void sumTempBunemanJzAlongZ();


	void sumCellParametersX();
	void sumCellParametersY();
	void sumCellParametersZ();
	void sumCellMatrixParametersX();
	void sumCellMatrixParametersY();
	void sumCellMatrixParametersZ();
	void sumNodeParametersX();
	void sumNodeParametersY();
	void sumNodeParametersZ();
	void sumNodeVectorParametersX();
	void sumNodeVectorParametersY();
	void sumNodeVectorParametersZ();
	void sumNodeMatrixParametersX();
	void sumNodeMatrixParametersY();
	void sumNodeMatrixParametersZ();
    void sumChargeDensityHatX();
    void sumChargeDensityHatY();
    void sumChargeDensityHatZ();

    void sumCellVectorParametersX();
    void sumCellVectorParametersY();
    void sumCellVectorParametersZ();

	void sumCellVectorTempParametersX(Vector3d*** array);
	void sumCellVectorTempParametersY(Vector3d*** array);
	void sumCellVectorTempParametersZ(Vector3d*** array);
	void sumCellMatrixParameterGeneralX(Matrix3d*** array, double* inBufferRight, double* outBufferRight, double* inBufferLeft, double* outBufferLeft);
	void sumCellMatrixParameterGeneralY(Matrix3d*** array, double* inBufferBack, double* outBufferBack, double* inBufferFront, double* outBufferFront);
	void sumCellMatrixParameterGeneralZ(Matrix3d*** array, double* inBufferTop, double* outBufferTop, double* inBufferBottom, double* outBufferBottom);
	void sumCellTempParametersX(double*** array);
	void sumCellTempParametersY(double*** array);
	void sumCellTempParametersZ(double*** array);
	void sumCellVectorParametersGeneralX(Vector3d*** array, double* inBufferRight, double* outBufferRight, double* inBufferLeft, double* outBufferLeft);
	void sumCellVectorParametersGeneralY(Vector3d*** array, double* inBufferBack, double* outBufferBack, double* inBufferFront, double* outBufferFront);
	void sumCellVectorParametersGeneralZ(Vector3d*** array, double* inBufferTop, double* outBufferTop, double* inBufferBottom, double* outBufferBottom);
	void sumCellTempMatrixParametersX(Matrix3d*** array);
	void sumCellTempMatrixParametersY(Matrix3d*** array);
	void sumCellTempMatrixParametersZ(Matrix3d*** array);

	void sumTempNodeParametersX(double*** array);
	void sumTempNodeParametersY(double*** array);
	void sumTempNodeParametersZ(double*** array);
	void sumCellParametersGeneralX(double*** array, double* inBufferRight, double* outBufferRight, double* inBufferLeft, double* outBufferLeft);
	void sumCellParametersGeneralY(double*** array, double* inBufferBack, double* outBufferBack, double* inBufferFront, double* outBufferFront);
	void sumCellParametersGeneralZ(double*** array, double* inBufferTop, double* outBufferTop, double* inBufferBottom, double* outBufferBottom);
	void sumTempNodeVectorParametersX(Vector3d*** array);
	void sumTempNodeVectorParametersY(Vector3d*** array);
	void sumTempNodeVectorParametersZ(Vector3d*** array);
	void sumTempNodeMatrixParametersX(Matrix3d*** array);
	void sumTempNodeMatrixParametersY(Matrix3d*** array);
	void sumTempNodeMatrixParametersZ(Matrix3d*** array);

	void sumNodeVectorParametersGeneralX(Vector3d*** vector, double* inBufferRight, double* outBufferRight,
	                                     double* inBufferLeft, double* outBufferLeft);
	void sumNodeVectorParametersGeneralY(Vector3d*** vector, double* inBufferBack, double* outBufferBack,
	                                     double* inBufferFront, double* outBufferFront);
	void sumNodeVectorParametersGeneralZ(Vector3d*** vector, double* inBufferTop, double* outBufferTop,
	                                     double* inBufferBottom, double* outBufferBottom);
};

#endif

