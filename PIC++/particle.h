#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <string>

#include "particleTypes.h"
#include "constants.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"

class Matrix3d;
class Vector3d;
class Simulation;

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
	double juttnerFunction[randomParameter];
	double juttnerValue[randomParameter];
	int jutnerNumber;

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
	double mc2;
public:
	int number;

	double mass;
	double charge;
	int chargeCount;
	double weight;
	ParticleTypes type;

	double beta;

	bool escaped;
	int crossBoundaryCount;

	Vector3d coordinates;
	Vector3d initialCoordinates;

	Vector3d initialMomentum;
	Vector3d prevMomentum;

	Matrix3d rotationTensor;

	CorrelationMap correlationMapCell;
	CorrelationMap correlationMapNode;


	double dx;
	double dy;
	double dz;

	Particle(int n, const double& m, int qcount, const double& q, const double& w, ParticleTypes type, const double& x0, const double& y0,
                 const double& z0, const double& px0, const double& py0, const double& pz0, const double& dx0, const double& dy0, const double& dz0, const double& deltaT);
	Particle(const Particle& particle);
	Particle(const Particle* const particle);

	double shapeFunctionX(const double& xvalue);
	double shapeFunctionY(const double& yvalue);
	double shapeFunctionZ(const double& zvalue);

	double shapeFunction(const double& xvalue, const double& yvalue, const double& zvalue);

	void setMomentum(const Vector3d& p);
	void setMomentum(const double& px, const double& py, const double& pz);
	void setMomentumX(const double& px);
	void setMomentumY(const double& py);
	void setMomentumZ(const double& pz);
	void addMomentum(const Vector3d& dp);
	void reflectMomentumX();
	void reflectMomentumY();
	void reflectMomentumZ();
	void copyMomentum(const Particle& particle);
	void copyMomentum(const Particle* particle);
	Vector3d getMomentum() const;
	double momentumAbs() const;
	Vector3d getVelocity();
	double velocityX();
	double velocityY();
	double velocityZ();
	void addVelocity(const Vector3d &v);
	double gammaFactor();
	double energy();
	double fullEnergy();

	void setMomentumByV(const Vector3d& v);

	static Vector3d evaluateVelocity(const Vector3d& p, const double& m);
};

#endif