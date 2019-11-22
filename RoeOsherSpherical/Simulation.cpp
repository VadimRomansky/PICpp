#include <list>
#include <time.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

//главная функция
void Simulation::simulate(){
	printf("initialization\n");
	const char* suffix;
	if(!serialized) {
		suffix = "w";
		initializeProfile();
	} else {
		suffix = "a";
	}
	//updateShockWavePoint();
	//shockWavePoint = rgridNumber/100;
	//updateGrid();
	updateMaxSoundSpeed();
	updateParameters();
	updateDiffusionCoef();
	updateTimeStep();

	printf("creating files\n");
	std::string fileNumber = "";
	int currentWriteNumber = 0;
	fileNumber = std::string("_") + convertIntToString(currentWriteNumber);
	std::string outputDir = "./output/";
	FILE* outIteration = fopen("./output/iterations.dat",suffix);
	fclose(outIteration);
	FILE* outExtraIteration = fopen("./output/extra_iterations.dat",suffix);
	fclose(outExtraIteration);
	FILE* outTempGrid = fopen((outputDir + "temp_grid" + fileNumber + ".dat").c_str(),suffix);
	fclose(outTempGrid);
	FILE* outShockWave = fopen("./output/shock_wave.dat",suffix);
	fclose(outShockWave);
	FILE* outFile = fopen((outputDir + "temp_radial_profile" + fileNumber + ".dat").c_str(),suffix);
	output(outFile,this);
	fclose(outFile);
	FILE* outDistribution;
	FILE* outFullDistribution;
	FILE* outCoordinateDistribution;
	outDistribution  = fopen((outputDir + "distribution" + fileNumber + ".dat").c_str(),suffix);
	outFullDistribution  = fopen((outputDir + "fullDistribution" + fileNumber + ".dat").c_str(),suffix);
	outCoordinateDistribution  = fopen((outputDir + "coordinateDistribution" + fileNumber + ".dat").c_str(),suffix);
	outputDistributionP3(outDistribution, outFullDistribution, outCoordinateDistribution, this);
	//outputDistributionP3(outDistribution, outFullDistribution, outCoordinateDistribution, this);
	fclose(outCoordinateDistribution);
	fclose(outFullDistribution);
	fclose(outDistribution);
	FILE* outField = fopen((outputDir + "field" + fileNumber + ".dat").c_str(),suffix);
	fclose(outField);
	FILE* coordinateField = fopen((outputDir + "coordinate_field" + fileNumber + ".dat").c_str(),suffix);
	fclose(coordinateField);
	FILE* outFullField = fopen((outputDir + "full_field" + fileNumber + ".dat").c_str(),suffix);
	fclose(outFullField);
	FILE* outCoef = fopen((outputDir + "diff_coef" + fileNumber + ".dat").c_str(),suffix);
	fclose(outCoef);
	FILE* xFile = fopen("./output/xfile.dat",suffix);
	fclose(xFile);
	FILE* pFile = fopen("./output/pfile.dat",suffix);
	fclose(pFile);
	FILE* kFile = fopen("./output/kfile.dat",suffix);
	fclose(kFile);
	updateShockWavePoint();
	outShockWave = fopen("./output/shock_wave.dat",suffix);
	double shockWaveR = 0;
	double gasSpeed = 0;
	if(shockWavePoint >= 0 && shockWavePoint <= rgridNumber){
		shockWaveR = grid[shockWavePoint];
		gasSpeed = middleVelocity[shockWavePoint - 1];
	}

	fprintf(outShockWave, "%d %lf %d %lf %lf %lf\n", 0, myTime, shockWavePoint, shockWaveR, shockWaveSpeed, gasSpeed);
	fclose(outShockWave);
	deltaT = min2(minT, deltaT);

	updateDiffusionCoef();
	//основной цикл
	//fakeMoveShockWave();


	while(myTime < maxTime && currentIteration < iterationNumber){
		++currentIteration;
		printf("iteration %d\n", currentIteration);
		printf("time = %lf\n", myTime);

		//if(currentIteration < startCRevaluation){
			evaluateHydrodynamic();
		//}
		//fakeMoveShockWave();
		
		if(currentIteration > startCRevaluation){
			//deltaT /= 10;
			//for(int j = 0; j < 10; ++j){
				evaluateCR();
			//}
			//deltaT *= 10;
		}

		if(currentIteration > startCRevaluation){
			evaluateField();
		}

		myTime = myTime + deltaT;

		updateAll();

		updateMaxSoundSpeed();
		updateShockWavePoint();
		updateParameters();

		updateTimeStep();
		deltaT = min2(minT, deltaT);
		if(currentIteration % writeParameter == 0){
			fileNumber = std::string("_") + convertIntToString(currentWriteNumber);
			//вывод на некоторых итерациях
			printf("outputing\n");
			printf("iteration %d\n", currentIteration);
			outFile = fopen((outputDir + "tamc_radial_profile" + fileNumber + ".dat").c_str(),"w");
			output(outFile, this);
			fclose(outFile);

			outShockWave = fopen("./output/shock_wave.dat","a");
			shockWaveR = 0;
			gasSpeed = 0;
			if(shockWavePoint >= 0 && shockWavePoint <= rgridNumber){
				shockWaveR = grid[shockWavePoint];
				gasSpeed = middleVelocity[shockWavePoint - 1];
			}

			fprintf(outShockWave, "%d %lf %d %lf %lf %lf\n", currentIteration, myTime, shockWavePoint, shockWaveR, shockWaveSpeed, gasSpeed);
			fclose(outShockWave);

			outDistribution = fopen((outputDir + "distribution" + fileNumber + ".dat").c_str(),"w");
			outFullDistribution = fopen((outputDir + "fullDistribution" + fileNumber + ".dat").c_str(),"w");
			outCoordinateDistribution = fopen((outputDir + "coordinateDistribution" + fileNumber + ".dat").c_str(),"w");

			outputDistributionP3(outDistribution, outFullDistribution, outCoordinateDistribution, this);
			//outputDistribution(outDistribution, outFullDistribution, outCoordinateDistribution, this);

			fclose(outCoordinateDistribution);
			fclose(outFullDistribution);
			fclose(outDistribution);

			outIteration = fopen("./output/iterations.dat","a");
			fprintf(outIteration, "%d %g %g %g %g %g %g\n", currentIteration, myTime, mass, totalMomentum, totalEnergy, injectedParticles, totalParticles);
			fclose(outIteration);

			outExtraIteration = fopen("./output/extra_iterations.dat","a");
			fprintf(outExtraIteration, "%d %g %g %g %g %g %g %g %g %g %g %g %g\n", currentIteration, myTime, mass, totalMomentum, totalEnergy, totalKineticEnergy, totalTermalEnergy, totalParticleEnergy, totalMagneticEnergy, injectedParticles, totalParticles, injectedEnergy, uGradPEnergy);
			fclose(outExtraIteration);

			outTempGrid = fopen((outputDir + "temp_grid" + fileNumber + ".dat").c_str(),"w");
			outputNewGrid(outTempGrid, this);
			fclose(outTempGrid);

			outFullField = fopen((outputDir + "full_field" + fileNumber + ".dat").c_str(),"w");
			coordinateField = fopen((outputDir + "coordinate_field" + fileNumber + ".dat").c_str(),"w");
			outField = fopen((outputDir + "field" + fileNumber + ".dat").c_str(),"w");
			outCoef = fopen((outputDir + "diff_coef" + fileNumber + ".dat").c_str(),"w");
			xFile = fopen("./output/xfile.dat","w");
			kFile = fopen("./output/kfile.dat","w");
			pFile = fopen("./output/pfile.dat","w");
			outputField(outField, coordinateField, outFullField, outCoef, xFile, kFile, pFile, this);
			fclose(outField);
			fclose(coordinateField);
			fclose(outFullField);
			fclose(outCoef);
			fclose(xFile);
			fclose(kFile);
			fclose(pFile);
			currentWriteNumber++;
		}

		if(currentIteration % serializeParameter == 0){
			FILE* hydroFile = fopen("./save/hydro.dat","w");
			FILE* distributionFile = fopen("./save/distribution.dat", "w");
			FILE* fieldFile = fopen("./save/field.dat", "w");
			FILE* gridFile = fopen("./save/grid.dat", "w");
			FILE* pgridFile = fopen("./save/pgrid.dat", "w");
			FILE* kgridFile = fopen("./save/kgrid.dat", "w");
			FILE* infoFile = fopen("./save/info.dat", "w");

			serialize(hydroFile, distributionFile, fieldFile, gridFile, pgridFile, kgridFile, infoFile, this);

			fclose(hydroFile);
			fclose(distributionFile);
			fclose(fieldFile);
			fclose(gridFile);
			fclose(pgridFile);
			fclose(kgridFile);
			fclose(infoFile);
		}
	}

	printf("end\n");
}

//расчет гидродинамики
void Simulation::evaluateHydrodynamic() {
	//printf("evaluating hydrodynamic\n");

	solveDiscontinious();
	CheckNegativeDensity();

	#pragma omp parallel for
	for(int i = 0; i < rgridNumber; ++i){
		tempDensity[i] = middleDensity[i];
		tempMomentum[i] = momentum(i);
		tempEnergy[i] = energy(i);
	}

	evaluateFluxes();

	//for Einfeldt
	updateFluxes();

	TracPenRadial(tempDensity, dFlux);

	TracPenRadial(tempMomentum, mFlux);
	for(int i = 0; i < rgridNumber - 1; ++i){
		tempMomentum[i] += deltaT*2*middlePressure[i]/middleGrid[i];

		tempMomentum[i] -= deltaT*(cosmicRayPressure[i+1] - cosmicRayPressure[i])/(deltaR[i]);
		
		if(i > 0 && i < rgridNumber-1){
			for(int k = 0; k < kgridNumber; ++k){
				tempMomentum[i] -= deltaT*0.5*(magneticField[i][k] - magneticField[i-1][k])*kgrid[k]*deltaLogK/deltaR[i];
			}
		}
		
	}
	TracPenRadial(tempEnergy, eFlux);


    //#pragma omp for
	for(int i = 1; i < rgridNumber-2; ++i){
		double deltaE = 0;
		tempEnergy[i] -= deltaT*middleVelocity[i]*(cosmicRayPressure[i+1] - cosmicRayPressure[i])/deltaR[i];
		uGradPEnergy += deltaT*middleVelocity[i]*((cosmicRayPressure[i+1] - cosmicRayPressure[i])/deltaR[i])*volume(i);
		for(int k = 0; k < kgridNumber; ++k){
			deltaE += deltaT*growth_rate[i][k]*magneticField[i][k]*kgrid[k]*deltaLogK;
			tempEnergy[i] -= deltaT*0.5*middleVelocity[i]*(magneticField[i][k] - magneticField[i-1][k])*kgrid[k]*deltaLogK/deltaR[i];
			//tempEnergy[i] += deltaT*0.5*magneticField[i][k]*(middleVelocity[i] - middleVelocity[i-1])*kgrid[k]*deltaLogK/deltaR[i];
			//tempEnergy[i] -= deltaT*growth_rate[i][k]*magneticField[i][k]*kgrid[k]*deltaLogK;
			//alertNaNOrInfinity(tempEnergy[i], "energy = NaN");
			alertNegative(tempEnergy[i], "energy < 0");
		}
        /*if(abs2(cosmicRayPressure[i+1] - cosmicRayPressure[i]) > 1E-100){
            vscattering[i] = -deltaE*deltaR[i]/(cosmicRayPressure[i+1] - cosmicRayPressure[i]);
		} else {
			vscattering[i] = 0;
		}*/
		
	}

	/*if(shockWavePoint > 0){
		for(int i = min2(shockWavePoint + 20, rgridNumber-1); i >= shockWavePoint; --i){
			vscattering[i] = vscattering[i+1];
		}
	}*/

}


//расчет разрывов, задача Римана

void Simulation::solveDiscontinious(){
    //double t1 = omp_get_wtime();
    //#pragma omp parallel for schedule(static, 1)
	for(int ompi = 0; ompi < numThreads; ++ ompi){
		for(int i = ompi+1; i < rgridNumber; i = i + numThreads){
            //int n = omp_get_thread_num();
			//printf("thread %d\n",n);
			double rho1, rho2, u1, u2, c1, c2;
			rho1 = middleDensity[i-1];
			rho2 = middleDensity[i];
			u1 = middleVelocity[i-1];
			u2 = middleVelocity[i];
			c1 = sqrt(_gamma*middlePressure[i-1]/rho1);
			c2 = sqrt(_gamma*middlePressure[i]/rho2);
			double h1 = (energy(i-1) + middlePressure[i-1])/rho1;
			double h2 = (energy(i) + middlePressure[i])/rho2;
			pointDensity[i] = sqrt(rho1*rho2);
			pointVelocity[i] = (sqrt(rho1)*middleVelocity[i-1] + sqrt(middleDensity[i])*middleVelocity[i])/(sqrt(rho1) + sqrt(rho2));
			pointEnthalpy[i] = (sqrt(rho1)*h1 + sqrt(rho2)*h2)/(sqrt(rho1) + sqrt(rho2));
			pointSoundSpeed[i] = sqrt((sqrt(rho1)*c1*c1 + sqrt(rho2)*c2*c2)/(sqrt(rho1) + sqrt(rho2)) + 0.5*(_gamma - 1)*(sqrt(rho1*rho2)/sqr(sqrt(rho1)+ sqrt(rho2)))*sqr(u1 - u2));
		}
	}
    //double t2 = omp_get_wtime();
    //printf("parllel %lf\n",t2-t1);
	pointEnthalpy[0] = (middlePressure[0] + energy(0))/middleDensity[0];
	pointDensity[0] = middleDensity[0];
	pointVelocity[0] = 0;
	pointSoundSpeed[0] = sqrt(_gamma*middlePressure[0]/middleDensity[0]);
}


//проверка шага по времени, чтобы плотность не уменьшалась слишком сильно и не становилась отрицательной

void Simulation::CheckNegativeDensity(){
	double dt = deltaT;
	for(int i = 0; i < rgridNumber; ++i){
		if(middleDensity[i]*volume(i) - dt*4*pi*(gridsquare[i+1]*densityFlux(i+1) - gridsquare[i]*densityFlux(i))< 0){
			dt = 0.5*(middleDensity[i]*volume(i)/(4*pi*(gridsquare[i+1]*densityFlux(i+1) - gridsquare[i]*densityFlux(i))));
			alertNaNOrInfinity(dt, "dt = NaN");
		}
	}
	deltaT = dt;
}


//Трак и Пен, немного модифицированные

//с учетом сферичности
void Simulation::TracPenRadial(double* u, double* flux){

	tempU[0] = u[0] - deltaT*(gridsquare[1]*flux[1] - gridsquare[0]*flux[0])/(middleGrid[0]*middleGrid[0]*deltaR[0]);
	int i;
    //#pragma omp parallel for private(i)
		for(i = 1; i < rgridNumber - 1; ++i){
			tempU[i] = u[i] - deltaT*((gridsquare[i+1]*flux[i+1] - gridsquare[i]*flux[i])/(middleGrid[i]*middleGrid[i]*deltaR[i]));
			/*if(tempU[i] < 0){
				printf("temp U < 0\n");
			}*/
		}
	tempU[rgridNumber - 1] = u[rgridNumber - 1];

	for(i = 0; i < rgridNumber; ++i){
		u[i] = tempU[i];
	}
}

double Simulation::densityFlux(int i){
	if(i < 0){
		printf("i < 0");
		exit(0);
	} else if(i >= 0 && i <= rgridNumber) {
		if(i == 0) return middleVelocity[0]*middleDensity[0];
		return pointDensity[i]*pointVelocity[i];
		//return middleDensity[i-1]*middleVelocity[i-1];
	} else {
		printf("i > rgridNumber");
		exit(0);
	}
	return 0;
}

double Simulation::momentum(int i){
	if(i < 0){
		printf("i < 0");
		exit(0);
	} else if(i >= 0 && i < rgridNumber) {
		return middleDensity[i]*middleVelocity[i];
	} else {
		printf("i >= rgridNumber");
		exit(0);
	}
	return 0;
}

double Simulation::energy(int i){
	if(i < 0){
		printf("i < 0");
		exit(0);
	} else if(i >= 0 && i < rgridNumber) {
		return middlePressure[i]/(_gamma - 1) + middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2;
	} else {
		printf("i >= rgridNumber");
		exit(0);
	}
	return 0;
}

double Simulation::kineticEnergy(int i){
	if(i < 0){
		printf("i < 0");
		exit(0);
	} else if(i >= 0 && i < rgridNumber) {
		return middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2;
	} else {
		printf("i >= rgridNumber");
		exit(0);
	}
	return 0;
}

double Simulation::termalEnergy(int i){
	if(i < 0){
		printf("i < 0");
		exit(0);
	} else if(i >= 0 && i < rgridNumber) {
		return middlePressure[i]/(_gamma - 1);
	} else {
		printf("i >= rgridNumber");
		exit(0);
	}
	return 0;
}

double Simulation::temperatureIn(int i){
	if(i < 0){
		printf("i < 0");
		exit(0);
	} else if(i >= 0 && i < rgridNumber) {
		return middlePressure[i]*massProton/(kBoltzman*middleDensity[i]);
	} else {
		printf("i >= rgridNumber");
		exit(0);
	}
	return 0;
}

double Simulation::soundSpeed(int i){
	if(i < 0){
		printf("i < 0");
		exit(0);
	} else if(i >= 0 && i < rgridNumber) {
		return sqrt(_gamma*middlePressure[i]/middleDensity[i]);
	} else {
		printf("i >= rgridNumber");
		exit(0);
	}
	return 0;
}

void Simulation::evaluateFluxes(){
	//dFlux[0] = middleDensity[0]*middleVelocity[0];
	//mFlux[0] = dFlux[0]*middleVelocity[0] + middlePressure[0];
	//eFlux[0] = middleVelocity[0]*(energy(0) + middlePressure[0]);
	dFlux[0] = 0;
	mFlux[0] = middlePressure[0];
	eFlux[0] = 0;

	int ompi = 0;
    //#pragma omp parallel for private(ompi)
	for(ompi = 0; ompi < numThreads; ++ ompi){
		for(int i = ompi+1; i < rgridNumber; i = i + numThreads){
			double** vectors = new double*[3];
			double* deltaS = new double[3];
			double* lambdaPlus = new double[3];
			double* lambdaMinus = new double[3];
			double* lambdaMod = new double[3];

			for(int k = 0; k < 3; ++k){
				vectors[k] = new double[3];
			}
			double leftDflux = middleVelocity[i-1]*middleDensity[i-1];
			double rightDflux = middleVelocity[i]*middleDensity[i];
			double leftMflux = leftDflux*middleVelocity[i-1] + middlePressure[i-1];
			double rightMflux = rightDflux*middleVelocity[i] + middlePressure[i];
			double leftEflux = leftDflux*middleVelocity[i-1]*middleVelocity[i-1]/2 + middleVelocity[i-1]*middlePressure[i-1]*_gamma/(_gamma-1);
			double rightEflux = rightDflux*middleVelocity[i]*middleVelocity[i]/2 + middleVelocity[i]*middlePressure[i]*_gamma/(_gamma-1);

			lambdaMod[0] = min2(middleVelocity[i-1] - sqrt(_gamma*middlePressure[i-1]/middleDensity[i-1]), pointVelocity[i] - pointSoundSpeed[i]);
			lambdaMod[1] = pointVelocity[i];
			lambdaMod[2] = max2(middleVelocity[i] + sqrt(_gamma*middlePressure[i]/middleDensity[i]), pointVelocity[i] + pointSoundSpeed[i]);

			double lambda1 = pointVelocity[i] - pointSoundSpeed[i];
			double lambda2 = pointVelocity[i];
			double lambda3 = pointVelocity[i] + pointSoundSpeed[i];

			lambdaPlus[0] = lambda1*(lambda1 > 0);
			lambdaPlus[1] = lambda2*(lambda2 > 0);
			lambdaPlus[2] = lambda3*(lambda3 > 0);

			lambdaMinus[0] = lambda1*(lambda1 < 0);
			lambdaMinus[1] = lambda2*(lambda2 < 0);
			lambdaMinus[2] = lambda3*(lambda3 < 0);

			deltaS[0] = ((middlePressure[i] - middlePressure[i-1]) - pointDensity[i]*pointSoundSpeed[i]*(middleVelocity[i] - middleVelocity[i-1]))/(2*sqr(pointSoundSpeed[i]));
			deltaS[1] = (sqr(pointSoundSpeed[i])*(middleDensity[i] - middleDensity[i-1]) - (middlePressure[i] - middlePressure[i-1]))/(2*sqr(pointSoundSpeed[i]));
			deltaS[2] = ((middlePressure[i] - middlePressure[i-1]) + pointDensity[i]*pointSoundSpeed[i]*(middleVelocity[i] - middleVelocity[i-1]))/(2*sqr(pointSoundSpeed[i]));

			vectors[0][0] = 1;
			vectors[0][1] = 2;
			vectors[0][2] = 1;
			vectors[1][0] = pointVelocity[i] - pointSoundSpeed[i];
			vectors[1][1] = 2*pointVelocity[i];
			vectors[1][2] = pointVelocity[i] + pointSoundSpeed[i];
			vectors[2][0] = pointEnthalpy[i] - pointVelocity[i]*pointSoundSpeed[i];
			vectors[2][1] = sqr(pointVelocity[i]);
			vectors[2][2] = pointEnthalpy[i] + pointVelocity[i]*pointSoundSpeed[i];

			dFlux[i] = (leftDflux + rightDflux)/2;
			mFlux[i] = (leftMflux + rightMflux)/2;
			eFlux[i] = (leftEflux + rightEflux)/2;

			//for density
			for(int j = 0; j < 3; ++j){
                dFlux[i] -= 0.5*abs2(lambdaMod[j])*deltaS[j]*vectors[0][j];
				dFluxPlus[i][j] = lambdaPlus[j]*deltaS[j]*vectors[0][j];
				dFluxMinus[i][j] = lambdaMinus[j]*deltaS[j]*vectors[0][j];
			}
			//for momentum
			for(int j = 0; j < 3; ++j){
                mFlux[i] -= 0.5*abs2(lambdaMod[j])*deltaS[j]*vectors[1][j];
				mFluxPlus[i][j] = lambdaPlus[j]*deltaS[j]*vectors[1][j];
				mFluxMinus[i][j] = lambdaMinus[j]*deltaS[j]*vectors[1][j];
			}
			//for energy
			for(int j = 0; j < 3; ++j){
                eFlux[i] -= 0.5*abs2(lambdaMod[j])*deltaS[j]*vectors[2][j];
				eFluxPlus[i][j] = lambdaPlus[j]*deltaS[j]*vectors[2][j];
				eFluxMinus[i][j] = lambdaMinus[j]*deltaS[j]*vectors[2][j];
			}

			for(int k = 0; k < 3; ++k){
				delete[] vectors[k];
			}
			delete[] vectors;
			delete[] deltaS;
			delete[] lambdaPlus;
			delete[] lambdaMinus;
			delete[] lambdaMod;
		}
	
	}
}

void Simulation::updateFluxes(){
	updateFluxes(dFlux, dFluxPlus, dFluxMinus);
	updateFluxes(mFlux, mFluxPlus, mFluxMinus);
	updateFluxes(eFlux, eFluxPlus, eFluxMinus);
}

void Simulation::updateFluxes(double* flux, double** fluxPlus, double** fluxMinus){
	double beta = 0.5;
	double phi = 1.0/3;

	for(int i = 1; i < rgridNumber-1; ++i){
		for(int j = 0; j < 3; ++j){
			flux[i] += 0.25*(1+phi)*minmod(fluxPlus[i][j], beta*fluxPlus[i-1][j])
					   +0.25*(1-phi)*minmod(beta*fluxPlus[i-1][j], fluxPlus[i][j])
					   -0.25*(1+phi)*minmod(fluxMinus[i][j], beta*fluxMinus[i+1][j])
					   -0.25*(1-phi)*minmod(beta*fluxMinus[i][j], fluxMinus[i+1][j]);

		}
	}
}

double Simulation::volume(int i){
	if(i < 0){
		printf("i < 0");
		exit(0);
	} else if(i >= 0 && i <= rgridNumber) {
		return 4*pi*(cube(grid[i+1]) - cube(grid[i]))/3;
	} else {
		printf("i > rgridNumber");
		exit(0);
	}
	return 0;
}

double Simulation::minmod(double a, double b){
	if(a*b > 0){
        if(abs2(a) < abs2(b)){
			return a;
		} else {
			return b;
		}
	} else {
		return 0;
	}
}

double Simulation::superbee(double a, double b){
    if(abs2(a) >= abs2(b)){
		return minmod(a, 2*b);
	} else {
		return minmod(2*a, b);
	}
}

//пересчет шага по времени и максимальной скорости звука

void Simulation::updateMaxSoundSpeed(){
	maxSoundSpeed = (sqrt(_gamma*(middlePressure[0]+cosmicRayPressure[0] + 0.5*magneticEnergy[0])/middleDensity[0]) + abs2(middleVelocity[0]));
	double cs = maxSoundSpeed;
	for(int i = 1; i < rgridNumber - 1; ++i){
        cs = (sqrt(_gamma*(middlePressure[i]+max2(cosmicRayPressure[i],cosmicRayPressure[i+1]) + 0.5*magneticEnergy[i])/middleDensity[i]) + abs2(middleVelocity[i]));
		if(cs > maxSoundSpeed){
			maxSoundSpeed = cs;
		}
	}
    cs = (sqrt(_gamma*(middlePressure[rgridNumber-1]+cosmicRayPressure[rgridNumber-1] + 0.5*magneticEnergy[rgridNumber-1])/middleDensity[rgridNumber - 1]) + abs2(middleVelocity[rgridNumber - 1]));
	if(cs > maxSoundSpeed){
		maxSoundSpeed = cs;
	}
}

void Simulation::updateTimeStep(){
	//by dx
	double tempdt = min2(deltaR[0]/maxSoundSpeed, deltaR[1]/maxSoundSpeed);
	for(int i = 1; i < rgridNumber - 1; ++i){
		if(deltaR[i]/maxSoundSpeed < tempdt){
			tempdt = deltaR[i]/maxSoundSpeed;
		}
		if(deltaR[i-1]/maxSoundSpeed < tempdt){
			tempdt = deltaR[i-1]/maxSoundSpeed;
		}
		if(deltaR[i+1]/maxSoundSpeed < tempdt){
			tempdt = deltaR[i+1]/maxSoundSpeed;
		}
	}
	if(deltaR[rgridNumber - 1]/maxSoundSpeed < tempdt){
		tempdt = deltaR[rgridNumber - 1]/maxSoundSpeed;
	}
	if(deltaR[rgridNumber - 2]/maxSoundSpeed < tempdt){
		tempdt = deltaR[rgridNumber - 2]/maxSoundSpeed;
	}
	//by dp
	for(int i = 1; i < rgridNumber; ++i){
        if(abs2(middleGrid[i]*middleGrid[i]*middleVelocity[i] - middleGrid[i-1]*middleGrid[i-1]*middleVelocity[i-1])*tempdt > abs2(3.0*deltaLogP*middleDeltaR[i]*gridsquare[i])){
            tempdt = 0.5*abs2(3.0*deltaLogP*gridsquare[i]*middleDeltaR[i]/(middleGrid[i]*middleGrid[i]*middleVelocity[i] - middleGrid[i-1]*middleGrid[i-1]*middleVelocity[i-1]));
		}
	}

	double maxDiffusion = 0;
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			if(diffusionCoef[i][j] > maxDiffusion){
				maxDiffusion  = diffusionCoef[i][j];
			}
		}
	}

	for(int i = 1; i < rgridNumber; ++i){
		double dx = (grid[i+1] - grid[i-1])/2;
		double dxp=grid[i+1]-grid[i];
		double dxm=grid[i]-grid[i-1];
		double a = ((sqr(middleGrid[i-1])*diffusionCoef[i-1][kgridNumber-1]/dxm + sqr(middleGrid[i])*diffusionCoef[i][kgridNumber-1]/dxp) - (sqr(middleGrid[i-1])*diffusionCoef[i-1][kgridNumber-1]/dxm) -sqr(middleGrid[i])*(diffusionCoef[i][kgridNumber-1]/dxp))/(2*dx);
		if(gridsquare[i] + tempdt*a < 0){
            tempdt = 0.5*gridsquare[i]/abs2(a);
		}
	}

	for(int i = 1; i < rgridNumber-1; ++i){
		double dx = (grid[i+1] - grid[i-1])/2;
		double dxp=grid[i+1]-grid[i];
		double dxm=grid[i]-grid[i-1];
		double xp = (grid[i+1]+grid[i])/2;
		double xm = (grid[i]+grid[i-1])/2;
		double dV = (xp*xp*xp - xm*xm*xm)/3;
		for(int j = 1; j < pgridNumber; ++j){
			double gkp = distributionFunction[i][j];
			double gkm = distributionFunction[i][j-1];
			double der = (1/(2*dV))*(xp*xp*diffusionCoef[i][j]*(distributionFunction[i+1][j] - distributionFunction[i][j])/dxp
							- xm*xm*diffusionCoef[i-1][j]*(distributionFunction[i][j] - distributionFunction[i-1][j])/dxm)
							- (1/dV)*(xp*xp*middleVelocity[i]*distributionFunction[i][j] - xm*xm*middleVelocity[i-1]*distributionFunction[i-1][j]);
			if((distributionFunction[i][j] > 0) && (distributionFunction[i][j] + der*tempdt < 0)){
                //tempdt = 0.5*abs2(distributionFunction[i][j]/der);
				if(tempdt < 0.1){
					//printf("ooo\n");
				}
			}
		}
	}

	for(int i = 1; i < rgridNumber; ++i){
		for(int k = 0; k < kgridNumber; ++k){
			if(tempdt*growth_rate[i][k] > 2){
                tempdt = 0.5*abs2(1/growth_rate[i][k]);
			}
		}
	}
	deltaT = 0.5*tempdt;
}

//определение точки ударной волны

void Simulation::updateShockWavePoint(){
	int tempShockWavePoint = -1;
	shockWaveT += deltaT;
	//double maxGrad = density0;
	//double maxGrad = 0.001*U0/upstreamR;
	double maxGrad = 0;
	for(int i = 10; i < rgridNumber - 10; ++i){
        //double grad = abs2((middleDensity[i] - middleDensity[i + 1])/middleDeltaR[i+1]);
		double grad = (middleVelocity[i] - middleVelocity[i + 1])/middleDeltaR[i];

		//double grad = (middleDensity[i]);
		if(grad > maxGrad){
			maxGrad = grad;
			tempShockWavePoint = i+1;
		}
	}
	shockWaveMoved = (tempShockWavePoint != shockWavePoint);
	if(shockWaveMoved && (tempShockWavePoint > -1) && (shockWavePoint > -1)){
		prevShockWavePoint = shockWavePoint;
		shockWaveSpeed = (grid[tempShockWavePoint] - grid[prevShockWavePoint])/(shockWaveT);
		shockWaveT = 0;
		for(int i = rgridNumber-1; i > 1; --i){
			for(int j = 0; j < pgridNumber; ++j){
				double sigma = 0.9;
				distributionFunction[i][j] = ((1-sigma)*distributionFunction[i-1][j]*volume(i-1) + sigma*distributionFunction[i][j]*volume(i))/(volume(i));
				//distributionFunction[i][j] = distributionFunction[i+1][j]*middleDeltaR[i+1]/middleDeltaR[i];
			}
		}
	}
	shockWavePoint = tempShockWavePoint;
}


//подсчет полной массы энергии и импульса

void Simulation::updateParameters(){
	mass = 0;
	totalMomentum = 0;
	totalEnergy = 0;
	totalKineticEnergy = 0;
	totalTermalEnergy = 0;
	totalParticles = 0;
	totalMagneticEnergy = 0;
	totalParticleEnergy = 0;
	for(int i = 0; i < rgridNumber; ++i){
		mass += middleDensity[i]*volume(i);
		totalMomentum += momentum(i)*volume(i);
		totalKineticEnergy += kineticEnergy(i)*volume(i);
		totalTermalEnergy += termalEnergy(i)*volume(i);
		for(int k = 0; k < kgridNumber; ++k){
			totalMagneticEnergy += magneticField[i][k]*kgrid[k]*deltaLogK*volume(i);
		}
		for(int j = 0; j < pgridNumber; ++j){
			double dp;
			if(j == 0){
				dp = (pgrid[j + 1] - pgrid[j]);
			} else if(j == pgridNumber -1){
				dp = (pgrid[j] - pgrid[j - 1]);
			} else {
				dp = (pgrid[j + 1] - pgrid[j - 1])/2;
			}
			double dr = 0;
			if(i == 0){
				dr = deltaR[0];
			} else if(i ==rgridNumber-1){
				dr = deltaR[rgridNumber-1];
			} else {
				dr = middleGrid[i] - middleGrid[i-1];
			}
			totalParticles += distributionFunction[i][j]*volume(i)*deltaLogP;
			//totalParticles += distributionFunction[i][j]*volume(i)*dp;
			if(j > goodMomentum){
				totalParticleEnergy += speed_of_light*distributionFunction[i][j]*volume(i)*dp;
			}

			//totalParticles += distributionFunction[i][j]*volume(i)*deltaLogP*cube(pgrid[j]);
			//totalParticleEnergy += speed_of_light*distributionFunction[i][j]*volume(i)*dp*cube(pgrid[j]);
		}
	}
	mass -= myTime*(0 - middleDensity[rgridNumber-1]*middleVelocity[rgridNumber-1]);
	totalParticles +=mass/massProton;
	totalMomentum -=  myTime*(middlePressure[0] - middleDensity[rgridNumber-1]*sqr(middleVelocity[rgridNumber-1]) - middlePressure[rgridNumber-1]);
	for(int k = 0; k < kgridNumber; ++k){
		totalMagneticEnergy -=  myTime*(0 - magneticField[rgridNumber-1][k]*middleVelocity[rgridNumber-1])*kgrid[k]*deltaLogK;
	}
	totalKineticEnergy -=  myTime*(0 - middleDensity[rgridNumber-1]*cube(middleVelocity[rgridNumber-1]))/2;
	totalTermalEnergy -=  myTime*(0 - middlePressure[rgridNumber-1]*middleVelocity[rgridNumber-1])*_gamma/(_gamma-1);
	totalEnergy = totalTermalEnergy + totalKineticEnergy + totalParticleEnergy + totalMagneticEnergy;
}

void Simulation::updateAll(){
	//hydrodinamic

	int ompi = 0;
    //#pragma omp parallel for private(ompi)
	for(ompi = 0; ompi < numThreads; ++ ompi){
		for(int i = ompi; i < rgridNumber; i = i + numThreads){
			if(tempDensity[i] < 0){
				//printf("density < 0\n");
				middleDensity[i] *= epsilon;
			} else {
				middleDensity[i] = tempDensity[i];
			}
			alertNaNOrInfinity(middleDensity[i], "density = NaN");
			double middleMomentum = tempMomentum[i];
			alertNaNOrInfinity(middleMomentum, "momentum = NaN");
			double middleEnergy = tempEnergy[i];
			if(tempEnergy[i] < 0){
				printf("energy < 0\n");
				exit(0);
			}
			alertNaNOrInfinity(middleEnergy, "energy = NaN");
			middleVelocity[i] = middleMomentum/middleDensity[i];

			middlePressure[i] = (middleEnergy - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(_gamma - 1);
			if(middleDensity[i] <= epsilon*density0){
				middleDensity[i] = epsilon*density0;
				middleVelocity[i] = 0;
				middlePressure[i] = middleDensity[i]*kBoltzman*temperature/massProton;
			}
			if(middlePressure[i] <= 0){
				middlePressure[i] = epsilon*middleDensity[i]*kBoltzman*temperature/massProton;
			}
		}
	}

	//cosmic rays
	if(currentIteration > startCRevaluation){
    //	#pragma omp parallel for private(ompi)
		for(ompi = 0; ompi < numThreads; ++ ompi){
			for(int i = ompi; i < rgridNumber; i = i + numThreads){
				for(int j = 0; j < pgridNumber; ++j){
					distributionFunction[i][j] = tempDistributionFunction[i][j];
				}
			}
		}
		
		evaluateCosmicRayPressure();
		evaluateCRFlux();
	}


	//field

	if(currentIteration > startCRevaluation){
		for(int i = 0; i < rgridNumber; i = i + 1){
			for(int k = 0; k < kgridNumber; ++k){
				magneticField[i][k] = tempMagneticField[i][k];
			}
		}
		

		for(int i = 0; i < rgridNumber; ++i){
			magneticEnergy[i] = 0;
			for(int k = 0; k < kgridNumber; ++k){
				magneticEnergy[i] += magneticField[i][k]*kgrid[k]*deltaLogK;
				largeScaleField[i][k] = sqrt(4*pi*magneticEnergy[i] + B0*B0);
			}
			magneticInductionSum[i] = sqrt(4*pi*magneticEnergy[i] + B0*B0);
			/*if(!stopAmplification){
				if(magneticInductionSum[i] > 5E-4){
					stopAmplification = true;
					setGrowthRateToZero();
				}
			}*/
		}
		
		if(currentIteration > startFieldEvaluation){
			updateDiffusionCoef();
			if(!stopAmplification){
				growthRate();
			}
		}
	}
}

void Simulation::fakeMoveShockWave(){
	double a = 1E15;
	double t0 = 1E7;
	double r = a*power(myTime + t0, 2.0/5);
	double u = 0.75*0.4*a*power(myTime + t0, -3.0/5);
	for(int i = 0; i < rgridNumber; ++i){
		if(middleGrid[i] < r){
			middleVelocity[i] = u;
		} else {
			middleVelocity[i] = 0;
		}
	}
}

void Simulation::setGrowthRateToZero(){
	for(int i = 0; i < rgridNumber; ++i){
		for(int k = 0; k < kgridNumber; ++k){
			growth_rate[i][k] = 0;
		}
	}
}
