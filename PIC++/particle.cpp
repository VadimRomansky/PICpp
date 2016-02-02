#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "particle.h"
#include "util.h"

Particle::Particle(int n, double m, double q, double w, ParticleTypes t, ParticleTypeContainer type_container, double x0, double y0, double z0, double px0, double py0, double pz0, double dx0, double dy0, double dz0) {
	number = n;

	mass = m;
	charge = q;
	weight = w;
	type = t;
	typeContainer = type_container;

	coordinates.x = x0;
	coordinates.y = y0;
	coordinates.z = z0;

	momentum.x = px0;
	momentum.y = py0;
	momentum.z = pz0;
	initialMomentum = momentum;

	dx = dx0;
	dy = dy0;
	dz = dz0;

	escaped = false;
}

Particle::Particle(const Particle& particle) {
	number = -1; //todo

	mass = particle.mass;
	charge = particle.charge;
	weight = particle.weight;
	type = particle.type;

	coordinates = particle.coordinates;
	momentum = particle.momentum;
	initialMomentum = particle.momentum;
	rotationTensor = particle.rotationTensor;

	dx = particle.dx;
	dy = particle.dy;
	dz = particle.dz;

	escaped = particle.escaped;
}

double Particle::shapeFunctionX(const double& xvalue) {
	return Bspline(coordinates.x, dx, xvalue);
}

double Particle::shapeFunctionY(const double& yvalue) {
	return Bspline(coordinates.y, dy, yvalue);
}

double Particle::shapeFunctionZ(const double& zvalue) {
	return Bspline(coordinates.z, dz, zvalue);
}

double Particle::shapeFunction(const double& xvalue, const double& yvalue, const double& zvalue) {
	return shapeFunctionX(xvalue) * shapeFunctionY(yvalue) * shapeFunctionZ(zvalue);
}

double Particle::momentumAbs() {
	return momentum.norm();
}

Vector3d Particle::velocity(double c) {
	double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	double mc2 = mass * c * c;
	double gamma_factor = sqrt(p2 * c * c + mc2 * mc2) / mc2;

	return momentum / (mass * gamma_factor);
}

double Particle::velocityX(double c) {
	double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	double mc2 = mass * c * c;
	double gamma_factor = sqrt(p2 * c * c + mc2 * mc2) / mc2;
	return momentum.x / (mass * gamma_factor);
}

double Particle::velocityY(double c) {
	double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	double mc2 = mass * c * c;
	double gamma_factor = sqrt(p2 * c * c + mc2 * mc2) / mc2;
	return momentum.y / (mass * gamma_factor);
}

double Particle::velocityZ(double c) {
	double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	double mc2 = mass * c * c;
	double gamma_factor = sqrt(p2 * c * c + mc2 * mc2) / mc2;
	return momentum.z / (mass * gamma_factor);
}

void Particle::addVelocity(Vector3d& v, double c) {
	if (v.norm() > c) {
		printf("ERROR v > c\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "v/c > 1 in addVelocity\n");
		fclose(errorLogFile);
		exit(0);
	}

	double c2 = c * c;

	Vector3d vel = velocity(c);
	vel += v;

	Matrix3d* rotation = Matrix3d::createBasisByOneVector(v);
	Matrix3d* inverse = rotation->Inverse();

	Vector3d rotatedV = (*inverse) * velocity(c);

	double gamma = 1 / (sqrt(1 - v.scalarMult(v) / c2));
	double vnorm = v.norm();
	double denominator = 1 + vnorm * rotatedV.z / c2;

	Vector3d shiftedV;

	shiftedV.z = (vnorm + rotatedV.z) / (denominator);
	shiftedV.y = rotatedV.y / (gamma * denominator);
	shiftedV.x = rotatedV.x / (gamma * denominator);

	Vector3d vel1 = (*rotation) * shiftedV;

	//todo relativistic!
	//setMomentumByV(vel, c);
	setMomentumByV(vel1, c);

	delete rotation;
	delete inverse;
}

void Particle::setMomentumByV(Vector3d v, double c) {
	if (v.norm() > c) {
		printf("ERROR v > c\n");
		printf("v = %g\n", v.norm());
		printf("c = %g\n", c);
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "v/c > 1 in setMomentumByV\n");
		fclose(errorLogFile);
		exit(0);
	}
	double gamma_factor = 1 / sqrt(1 - v.scalarMult(v) / (c * c));
	momentum = v * (mass * gamma_factor);
}

double Particle::gammaFactor(double c) {
	double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	double mc2 = mass * c * c;
	double result = sqrt(p2 * c * c + mc2 * mc2) / mc2;
	alertNaNOrInfinity(result, "gamma = NaN");
	return result;
}

double Particle::energy(double c) {
	double gamma_factor = gammaFactor(c);
	//return gamma_factor*mass*c*c;
	return (gamma_factor - 1) * mass * c * c;
}
