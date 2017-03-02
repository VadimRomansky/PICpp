#include <string>
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "constants.h"
#include "particle.h"
#include "util.h"

Particle::Particle(int n, double m, int qcount, double q, double w, ParticleTypes t, double x0, double y0,
				   double z0, double px0, double py0, double pz0, double dx0, double dy0, double dz0) {
	number = n;

	mass = m;
	chargeCount = qcount;
	charge = q;
	weight = w;
	type = t;

	coordinates.x = x0;
	coordinates.y = y0;
	coordinates.z = z0;

	momentum.x = px0;
	momentum.y = py0;
	momentum.z = pz0;
	initialMomentum = momentum;
	//prevMomentum = momentum;

	dx = dx0;
	dy = dy0;
	dz = dz0;

	escaped = false;
	velocityCashed = false;
	gammaCashed = false;
    crossBoundaryCount = 0;
}

Particle::Particle(const Particle& particle) {
	number = particle.number; //todo

	mass = particle.mass;
	charge = particle.charge;
	chargeCount = particle.chargeCount;
	weight = particle.weight;
	type = particle.type;

	coordinates = particle.coordinates;
	momentum = particle.momentum;
	//prevMomentum = particle.prevMomentum;
	rotationTensor = particle.rotationTensor;
	initialMomentum = particle.initialMomentum;

	correlationMapNode = particle.correlationMapNode;
	correlationMapCell = particle.correlationMapCell;

	dx = particle.dx;
	dy = particle.dy;
	dz = particle.dz;

	escaped = particle.escaped;
	velocityCashed = false;
	gammaCashed = false;
    crossBoundaryCount = particle.crossBoundaryCount;
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

void Particle::setMomentum(Vector3d& p) {
	velocityCashed = false;
	gammaCashed = false;
	momentum = p;
}

void Particle::setMomentum(double& px, double& py, double& pz) {
	velocityCashed = false;
	gammaCashed = false;
	momentum.x = px;
	momentum.y = py;
	momentum.z = pz;
}

void Particle::setMomentumX(double& px) {
	velocityCashed = false;
	gammaCashed = false;
	momentum.x = px;
}

void Particle::setMomentumY(double& py) {
	velocityCashed = false;
	gammaCashed = false;
	momentum.y = py;
}

void Particle::setMomentumZ(double& pz) {
	velocityCashed = false;
	gammaCashed = false;
	momentum.z = pz;
}

void Particle::addMomentum(Vector3d dp) {
	velocityCashed = false;
	gammaCashed = false;
	momentum = momentum + dp;
}

void Particle::reflectMomentumX() {
	velocityCashed = false;
	momentum.x = - momentum.x;
}

void Particle::reflectMomentumY() {
	velocityCashed = false;
	momentum.y = - momentum.y;
}

void Particle::reflectMomentumZ() {
	velocityCashed = false;
	momentum.z = - momentum.z;
}

void Particle::copyMomentum(const Particle& particle) {
	velocityCashed = false;
	gammaCashed = false;
	momentum = particle.momentum;
}

void Particle::copyMomentum(const Particle* particle) {
	velocityCashed = false;
	gammaCashed = false;
	momentum = particle->momentum;
}

Vector3d Particle::getMomentum() const {
	return momentum;
}

double Particle::momentumAbs() const {
	return momentum.norm();
}

Vector3d Particle::getVelocity(const double& c) {
	if(!velocityCashed) {
		double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
		alertNaNOrInfinity(p2, "p2 = NaN in particle::getVelocity\n");
		double c2 = c*c;
		double mc2 = mass * c2;
		if (p2 < relativisticPrecision*relativisticPrecision * mass * mass * c2) {
			velocity = momentum / mass;
		} else {
			if(!gammaCashed){
				gamma = sqrt(p2 * c2 + mc2 * mc2) / mc2;
				gammaCashed = true;
			}
			velocity = momentum / (mass * gamma);
		}
		velocityCashed = true;
	}
	return velocity;
}

double Particle::velocityX(const double& c) {
	if(velocityCashed) {
		return velocity.x;
	}
	double c2 = c*c;
	double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	double mc2 = mass * c2;
	if (p2 < relativisticPrecision * mass * mass * c2) {
		return momentum.x / mass;
	}
	if(!gammaCashed){
		gamma = sqrt(p2 * c2 + mc2 * mc2) / mc2;
		gammaCashed = true;
	}
	return momentum.x / (mass * gamma);
}

double Particle::velocityY(const double& c) {
	if(velocityCashed) {
		return velocity.y;
	}
	double c2 = c*c;
	double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	double mc2 = mass * c2;
	if (p2 < relativisticPrecision * mass * mass * c2) {
		return momentum.y / mass;
	}
	if(!gammaCashed){
		gamma = sqrt(p2 * c2 + mc2 * mc2) / mc2;
		gammaCashed = true;
	}
	return momentum.y / (mass * gamma);
}

double Particle::velocityZ(const double& c) {
	if(velocityCashed) {
		return velocity.z;
	}
	double c2 = c*c;
	double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	double mc2 = mass * c2;
	if (p2 < relativisticPrecision * mass * mass * c2) {
		return momentum.z / mass;
	}
	if(!gammaCashed){
		gamma = sqrt(p2 * c2 + mc2 * mc2) / mc2;
		gammaCashed = true;
	}
	return momentum.z / (mass * gamma);
}

void Particle::addVelocity(const Vector3d &v, const double& c) {
	velocityCashed = false;
	gammaCashed = false;
	if (v.norm() > c) {
		printf("ERROR v > c in particle::addVelocity\n");
		std::string outputDir = outputDirectory;
		FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		//FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "v/c > 1 in addVelocity\n");
		fclose(errorLogFile);
		exit(0);
	}

	double c2 = c * c;

	//Vector3d vel = getVelocity(c);
	//vel += v;

	Matrix3d* rotation = Matrix3d::createBasisByOneVector(v);
	Matrix3d* inverse = rotation->Inverse();

	Vector3d rotatedV = (*inverse) * getVelocity(c);

	gamma = 1 / (sqrt(1 - v.scalarMult(v) / c2));
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

void Particle::setMomentumByV(const Vector3d& v, const double& c) {
	if (v.norm() > c) {
		printf("ERROR v > c in particle::setMomentumByV\n");
		printf("v = %g\n", v.norm());
		printf("c = %g\n", c);
		std::string outputDir = outputDirectory;
		//FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "v/c > 1 in setMomentumByV\n");
		fclose(errorLogFile);
		exit(0);
	}
	velocity = v;
	velocityCashed = true;
	double v2 = v.scalarMult(v);
	double c2 = c*c;
	if(v2 < relativisticPrecision*relativisticPrecision*c2){
		momentum = v * mass;
		return;
	}
	gamma = 1 / sqrt(1 - v2 / c2);
	gammaCashed = true;
	momentum = v * (mass * gamma);
}

double Particle::gammaFactor(const double& c) {
	if(!gammaCashed){
		double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
		double mc2 = mass * c * c;
		gamma = sqrt(p2 * c * c + mc2 * mc2) / mc2;
		alertNaNOrInfinity(gamma, "gamma = NaN");
		gammaCashed = true;
	}
	return gamma;
}

double Particle::energy(const double& c) {
	double gamma_factor = gammaFactor(c);
	//return gamma_factor*mass*c*c;
	return (gamma_factor - 1) * mass * c * c;
}
