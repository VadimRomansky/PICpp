#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "vector3d.h"

enum ParticleTypes {PROTON, ELECTRON, POSITRON, ALPHA};

class Particle{
public:
	int number;

	double mass;
	double charge;
	double weight;
	ParticleTypes type;

	bool escaped;

	double x;
	double y;
	double z;

	Vector3d momentum;

	double dx;

	Matrix3d rotationTensor;


	Particle(int n, double m, double q, double w, ParticleTypes type, double x0, double px0, double py0, double pz0, double dx0);
	Particle(const Particle& particle);

	double shapeFunctionX(const double& xvalue);

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