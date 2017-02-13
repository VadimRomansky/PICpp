#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <string>

#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"

class Matrix3d;
class Vector3d;

enum ParticleTypes {ELECTRON, PROTON, POSITRON, ALPHA, DEUTERIUM, HELIUM3, OXYGEN_PLUS3, SILICON_PLUS1};

class ParticleTypeContainer {
public:
	ParticleTypes type;
	int number;
	double mass;
	double charge;
	int chargeCount;
	int particlesPerBin;
	double particesDeltaX;
	double particesDeltaY;
	double particesDeltaZ;
	double concentration;
	double injectionLength;
	double temperatureX;
	double temperatureY;
	double temperatureZ;
	double alphaParallel;
	double alphaNormal;
	double anisotropy;
	double parallelTemperatureEvaluated;
	double normalTemperatureEvaluated;
	double generalWeight;

	std::string typeName;
};

class CorrelationMap {
public:
	int xindex[splineOrder+2];
	int yindex[splineOrder+2];
	int zindex[splineOrder+2];
	double xcorrelation[splineOrder+2];
	double ycorrelation[splineOrder+2];
	double zcorrelation[splineOrder+2];

	CorrelationMap();
};

class Particle{
private:
	Vector3d momentum;
	Vector3d velocity;
	double gamma;
	bool velocityCashed;
	bool gammaCashed;
public:
	int number;

	double mass;
	double charge;
	int chargeCount;
	double weight;
	ParticleTypes type;

	bool escaped;
	int crossBoundaryCount;

	Vector3d coordinates;
	Vector3d initialCoordinates;

	Vector3d initialMomentum;
	//Vector3d prevMomentum;

	Matrix3d rotationTensor;

	CorrelationMap correlationMapCell;
	CorrelationMap correlationMapNode;


	double dx;
	double dy;
	double dz;

	Particle(int n, double m, int qcount, double q, double w, ParticleTypes type, double x0, double y0,
                 double z0, double px0, double py0, double pz0, double dx0, double dy0, double dz0);
	Particle(const Particle& particle);

	double shapeFunctionX(const double& xvalue);
	double shapeFunctionY(const double& yvalue);
	double shapeFunctionZ(const double& zvalue);

	double shapeFunction(const double& xvalue, const double& yvalue, const double& zvalue);

	void setMomentum(Vector3d& p);
	void setMomentum(double& px, double& py, double& pz);
	void setMomentumX(double& px);
	void setMomentumY(double& py);
	void setMomentumZ(double& pz);
	void addMomentum(Vector3d dp);
	void reflectMomentumX();
	void reflectMomentumY();
	void reflectMomentumZ();
	void copyMomentum(const Particle& particle);
	void copyMomentum(const Particle* particle);
	Vector3d getMomentum() const;
	double momentumAbs() const;
	Vector3d getVelocity(const double& c);
	double velocityX(const double& c);
	double velocityY(const double& c);
	double velocityZ(const double& c);
	void addVelocity(const Vector3d &v, const double& c);
	double gammaFactor(const double& c);
	double energy(const double& c);

	void setMomentumByV(const Vector3d& v, const  double& c);
};

#endif