#include <string>
#include <mpi.h>
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "simulation.h"
#include "constants.h"
#include "particle.h"
#include "util.h"
#include "paths.h"

Particle::Particle(int n, const double& m, int qcount, const double& q, const double& w, ParticleTypes t, const double& x0, const double& y0,
                   const double& z0, const double& px0, const double& py0, const double& pz0, const double& dx0, const double& dy0, const double& dz0, const double& deltaT) {
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
	prevMomentum = momentum;

	dx = dx0;
	dy = dy0;
	dz = dz0;

	escaped = false;
	velocityCashed = false;
	gammaCashed = false;
	crossBoundaryCount = 0;
	mc2 = mass;

	beta = charge * deltaT / (2.0 * mass * Simulation::speed_of_light_normalized);
}

Particle::Particle(const Particle& particle) {
	//private fields
	momentum = particle.momentum;
	velocity = particle.velocity;
	velocityCashed = particle.velocityCashed;
	gammaCashed = particle.gammaCashed;
	gamma = particle.gamma;
	mc2 = particle.mc2;
	//


	number = particle.number; //todo

	mass = particle.mass;
	charge = particle.charge;
	chargeCount = particle.chargeCount;
	weight = particle.weight;
	type = particle.type;

	escaped = particle.escaped;
	crossBoundaryCount = particle.crossBoundaryCount;

	coordinates = particle.coordinates;
	initialCoordinates = particle.initialCoordinates;

	prevMomentum = particle.prevMomentum;
	rotationTensor = particle.rotationTensor;
	initialMomentum = particle.initialMomentum;

	correlationMapNode = particle.correlationMapNode;
	correlationMapCell = particle.correlationMapCell;

	dx = particle.dx;
	dy = particle.dy;
	dz = particle.dz;

	beta = particle.beta;
}

Particle::Particle(const Particle* const particle) {
	//private fields
	momentum = particle->momentum;
	velocity = particle->velocity;
	velocityCashed = particle->velocityCashed;
	gammaCashed = particle->gammaCashed;
	gamma = particle->gamma;
	mc2 = particle->mc2;
	//


	number = particle->number; //todo

	mass = particle->mass;
	charge = particle->charge;
	chargeCount = particle->chargeCount;
	weight = particle->weight;
	type = particle->type;

	escaped = particle->escaped;
	crossBoundaryCount = particle->crossBoundaryCount;

	coordinates = particle->coordinates;
	initialCoordinates = particle->initialCoordinates;

	rotationTensor = particle->rotationTensor;
	initialMomentum = particle->initialMomentum;

	correlationMapNode = particle->correlationMapNode;
	correlationMapCell = particle->correlationMapCell;

	dx = particle->dx;
	dy = particle->dy;
	dz = particle->dz;

	beta = particle->beta;
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

void Particle::setMomentum(const Vector3d& p) {
	velocityCashed = false;
	gammaCashed = false;
	momentum = p;
}

void Particle::setMomentum(const double& px, const double& py, const double& pz) {
	velocityCashed = false;
	gammaCashed = false;
	momentum.x = px;
	momentum.y = py;
	momentum.z = pz;
}

void Particle::setMomentumX(const double& px) {
	velocityCashed = false;
	gammaCashed = false;
	momentum.x = px;
}

void Particle::setMomentumY(const double& py) {
	velocityCashed = false;
	gammaCashed = false;
	momentum.y = py;
}

void Particle::setMomentumZ(const double& pz) {
	velocityCashed = false;
	gammaCashed = false;
	momentum.z = pz;
}

void Particle::addMomentum(const Vector3d& dp) {
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

Vector3d Particle::getVelocity() {
	if (!velocityCashed) {
		double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
		//alertNaNOrInfinity(p2, "p2 = NaN in particle::getVelocity\n");
		if (p2 < relativisticPrecision * mass * mass * Simulation::speed_of_light_normalized_sqr) {
			velocity = momentum / mass;
		} else {
			if (!gammaCashed) {
				gamma = sqrt(p2 * Simulation::speed_of_light_normalized_sqr/(mc2*mc2) + 1.0);
				gammaCashed = true;
			}
			velocity = momentum / (mass * gamma);
		}
		velocityCashed = true;
	}
	return velocity;
}

double Particle::velocityX() {
	if (velocityCashed) {
		return velocity.x;
	}
	double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	if (p2 < relativisticPrecision * mass * mass * Simulation::speed_of_light_normalized_sqr) {
		return momentum.x / mass;
	}
	if (!gammaCashed) {
		gamma = sqrt(p2 * Simulation::speed_of_light_normalized_sqr/(mc2*mc2) + 1.0);
		gammaCashed = true;
	}
	return momentum.x / (mass * gamma);
}

double Particle::velocityY() {
	if (velocityCashed) {
		return velocity.y;
	}
	double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	if (p2 < relativisticPrecision * mass * mass * Simulation::speed_of_light_normalized_sqr) {
		return momentum.y / mass;
	}
	if (!gammaCashed) {
		gamma = sqrt(p2 * Simulation::speed_of_light_normalized_sqr/(mc2*mc2) + 1.0);
		gammaCashed = true;
	}
	return momentum.y / (mass * gamma);
}

double Particle::velocityZ() {
	if (velocityCashed) {
		return velocity.z;
	}
	double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
	if (p2 < relativisticPrecision * mass * mass * Simulation::speed_of_light_normalized_sqr) {
		return momentum.z / mass;
	}
	if (!gammaCashed) {
		gamma = sqrt(p2 * Simulation::speed_of_light_normalized_sqr/(mc2*mc2) + 1.0);
		gammaCashed = true;
	}
	return momentum.z / (mass * gamma);
}

void Particle::addVelocity(const Vector3d& v) {
	velocityCashed = false;
	gammaCashed = false;
	if (v.norm() > Simulation::speed_of_light_normalized) {
		printf("ERROR v > c in particle::addVelocity\n");
		std::string outputDir = outputDirectory;
		FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		//FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "v/c > 1 in addVelocity\n");
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}


	//Vector3d vel = getVelocity(c);
	//vel += v;

	Matrix3d* rotation = Matrix3d::createBasisByOneVector(v);
	Matrix3d* inverse = rotation->Inverse();

	Vector3d rotatedV = (*inverse) * getVelocity();

	gamma = 1 / (sqrt(1 - v.scalarMult(v) / Simulation::speed_of_light_normalized_sqr));
	double vnorm = v.norm();
	double denominator = 1 + vnorm * rotatedV.z / Simulation::speed_of_light_normalized_sqr;

	Vector3d shiftedV;

	shiftedV.z = (vnorm + rotatedV.z) / (denominator);
	shiftedV.y = rotatedV.y / (gamma * denominator);
	shiftedV.x = rotatedV.x / (gamma * denominator);

	Vector3d vel1 = (*rotation) * shiftedV;

	//todo relativistic!
	//setMomentumByV(vel, c);
	setMomentumByV(vel1);

	delete rotation;
	delete inverse;
}

void Particle::setMomentumByV(const Vector3d& v) {
	if (v.norm() > Simulation::speed_of_light_normalized) {
		printf("ERROR v > c in particle::setMomentumByV\n");
		printf("v = %g\n", v.norm());
		printf("c = %g\n", Simulation::speed_of_light_normalized);
		std::string outputDir = outputDirectory;
		//FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
		fprintf(errorLogFile, "v/c > 1 in setMomentumByV\n");
		fclose(errorLogFile);
		MPI_Finalize();
		exit(0);
	}
	velocity = v;
	velocityCashed = true;
	double v2 = v.scalarMult(v);
	if (v2 < relativisticPrecision * Simulation::speed_of_light_normalized_sqr) {
		momentum = v * mass;
		return;
	}
	gamma = 1 / sqrt(1 - v2 / Simulation::speed_of_light_normalized_sqr);
	gammaCashed = true;
	momentum = v * (mass * gamma);
}

Vector3d Particle::evaluateVelocity(const Vector3d& p, const double& m) {
	double p2 = p.x * p.x + p.y * p.y + p.z * p.z;
	double mc2 = m * Simulation::speed_of_light_normalized_sqr;
	if (p2 < relativisticPrecision * m * m * Simulation::speed_of_light_normalized_sqr) {
		return p / m;
	} else {
		double g = sqrt(p2 * Simulation::speed_of_light_normalized_sqr + mc2 * mc2) / mc2;;

		return p / (m * g);
	}
}

double Particle::gammaFactor() {
	if (!gammaCashed) {
		double p2 = momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z;
		gamma = sqrt(p2 * Simulation::speed_of_light_normalized_sqr/(mc2*mc2) + 1.0);
		//alertNaNOrInfinity(gamma, "gamma = NaN");
		gammaCashed = true;
	}
	return gamma;
}

double Particle::energy() {
	double gamma_factor = gammaFactor();
	//return gamma_factor*mass*c*c;
	return (gamma_factor - 1) * mc2;
}

double Particle::fullEnergy() {
	double gamma_factor = gammaFactor();
	return gamma_factor * mc2;
}
 