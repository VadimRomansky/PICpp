#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "vector3d.h"

enum ParticleTypes {PROTON, ELECTRON, POSITRON, ALPHA};

class ParticleTypeContainer {
public:
	ParticleTypes type;
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
};

class Particle{
public:
	int number;

	double mass;
	double charge;
	int chargeCount;
	double weight;
	ParticleTypes type;
	ParticleTypeContainer typeContainer; //maybe unnecessary

	bool escaped;

	Vector3d coordinates;

	Vector3d momentum;
	Vector3d initialMomentum;

	Matrix3d rotationTensor;


	double dx;
	double dy;
	double dz;

	Particle(int n, double m, int qcount, double q, double w, ParticleTypes type, ParticleTypeContainer typeContainer, double x0, double y0, double z0, double px0, double py0, double pz0, double dx0, double dy0, double dz0);
	Particle(const Particle& particle);

	double shapeFunctionX(const double& xvalue);
	double shapeFunctionY(const double& yvalue);
	double shapeFunctionZ(const double& zvalue);

	double shapeFunction(const double& xvalue, const double& yvalue, const double& zvalue);

	double momentumAbs();
	Vector3d velocity(double c);
	double velocityX(double c);
	double velocityY(double c);
	double velocityZ(double c);
	void addVelocity(Vector3d& v, double c);
	double gammaFactor(double c);
	double energy(double c);

	void setMomentumByV(Vector3d v, double c);
};

#endif