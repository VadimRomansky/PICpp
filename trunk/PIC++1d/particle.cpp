#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "particle.h"
#include "util.h"
#include "vector3d.h"
#include "matrix3d.h"

Particle::Particle(int n, double m, double q, double w, ParticleTypes t, double x0, double px0, double py0, double pz0, double dx0){
	number = n;

	mass = m;
	charge = q;
	weight = w;
	type = t;

	x = x0;

	momentum.x = px0;
	momentum.y = py0;
	momentum.z = pz0;

	dx = dx0;

	y = 0;
	z = 0;
}

Particle::Particle(const Particle& particle){
	number = -1; //todo

	mass = particle.mass;
	charge = particle.charge;
	weight = particle.weight;
	type = particle.type;

	x = particle.x;
	momentum = particle.momentum;
	rotationTensor = particle.rotationTensor;

	dx = particle.dx;

	y = particle.y;
	z = particle.z;
}

double Particle::shapeFunctionX(const double& xvalue){
	return Bspline(x, dx, xvalue);
}

double Particle::momentumAbs(){
	return momentum.norm();
}

Vector3d Particle::velocity(double c){
	double p2 = momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z;
	double mc2 = mass*c*c;
	double gamma_factor = sqrt(p2*c*c + mc2*mc2)/mc2;

	return momentum/(mass*gamma_factor);
}

double Particle::velocityX(double c){
	double p2 = momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z;
	double mc2 = mass*c*c;
	double gamma_factor = sqrt(p2*c*c + mc2*mc2)/mc2;
	return momentum.x/(mass*gamma_factor);
}

double Particle::velocityY(double c){
	double p2 = momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z;
	double mc2 = mass*c*c;
	double gamma_factor = sqrt(p2*c*c + mc2*mc2)/mc2;
	return momentum.y/(mass*gamma_factor);
}
double Particle::velocityZ(double c){
	double p2 = momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z;
	double mc2 = mass*c*c;
	double gamma_factor = sqrt(p2*c*c + mc2*mc2)/mc2;
	return momentum.z/(mass*gamma_factor);
}

void Particle::addVelocity(Vector3d& v, double c) {
	if(v.norm() > c) {
		printf("ERROR v > c\n");
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "v/c > 1 in addVelocity\n");
		fclose(errorLogFile);
		exit(0);		
	}
	double c2 = c*c;

	Vector3d vel = velocity(c);
	vel += v;

	Matrix3d* rotation = Matrix3d::createBasisByOneVector(v);
	Matrix3d* inverse = rotation->Inverse();

	Vector3d rotatedV = (*inverse)*velocity(c);

	double gamma = 1/(sqrt(1 - v.scalarMult(v)/c2));
	double vnorm = v.norm();
	double denominator = 1 + vnorm*rotatedV.x/c2;

	Vector3d shiftedV;

	shiftedV.z = (vnorm + rotatedV.z)/(denominator);
	shiftedV.y = rotatedV.y/(gamma*denominator);
	shiftedV.x = rotatedV.x/(gamma*denominator);

	Vector3d vel1 = (*rotation)*shiftedV;

	//todo relativistic!
	//setMomentumByV(vel, c);
	setMomentumByV(vel1, c);

	delete rotation;
	delete inverse;
}

void Particle::setMomentumByV(Vector3d v, double c){
	if(v.norm() > c){
		printf("ERROR v > c\n");
		printf("v = %g\n",v.norm());
		printf("c = %g\n", c);
		FILE* errorLogFile = fopen("./output/errorLog.dat", "w");
		fprintf(errorLogFile, "v/c > 1 in setMomentumByV\n");
		fclose(errorLogFile);
		exit(0);
	}
	double gamma_factor = 1/sqrt(1 - v.scalarMult(v)/(c*c));
	momentum = v*(mass*gamma_factor);
}

double Particle::gammaFactor(double c){
	double p2 = momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z;
	double mc2 = mass*c*c;
	double result =  sqrt(p2*c*c + mc2*mc2)/mc2;
	alertNaNOrInfinity(result, "gamma = NaN");
	return result;
}

double Particle::energy(double c) {
	double gamma_factor = gammaFactor(c);
	//return gamma_factor*mass*c*c;
	return (gamma_factor-1)*mass*c*c;
}