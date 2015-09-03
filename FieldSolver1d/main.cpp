#include "stdio.h"
#include <stdlib.h>
#include <time.h>
#include <crtdbg.h>
#include <omp.h>

#include "constants.h"
#include "matrix3d.h"
#include "util.h"
#include "input.h"
#include "simulation.h"

//test solving matrix

/*int main(){
	double** leftHalf;
	double* rightPart;
	double* output;

	leftHalf = new double*[6];
	rightPart = new double[6];
	output = new double[6];

	for(int i = 0; i < 6; ++i){
		leftHalf[i] = new double[3];
		for(int j = 0; j < 3; ++j){
			leftHalf[i][j] = 0;
		}
		rightPart[i] = 1;
		output[i] = 0;
	}
	
	
	//	1 2 3 0 0 0      1
	//	2 1 3 0 0 0      1
	//	4 1 1 0 0 0      1
	//  0 3 0 1 0 0      1
	//	0 0 2 0 1 0      1
	//	2 0 0 0 0 0      1

	//	x = 
	
	//	0.16667
	//	0.16667
	//	0.16667
	//	0.5
	//	0.66667
	//	0.66667

	leftHalf[0][0] = 1;
	leftHalf[1][1] = 1;
	leftHalf[2][2] = 1;
	leftHalf[5][0] = 2;
	leftHalf[3][1] = 3;
	leftHalf[4][2] = 2;
	leftHalf[0][1] = 2;
	leftHalf[0][2] = 3;
	leftHalf[1][0] = 2;
	leftHalf[1][2] = 3;
	leftHalf[2][0] = 4;
	leftHalf[2][1] = 1;


	solveSpecialMatrix(leftHalf, rightPart, output);

	for(int i = 0; i < 6; ++i){
		printf("%lf\n", output[i]);
	}

	for(int i = 0; i < 6; ++i){
		delete[] leftHalf[i];
	}
	delete[] leftHalf;
	delete[] rightPart;
	delete[] output;

	return 0;
}*/


//test McDonald
/*int main()
{
	srand (time(NULL));
	//double x = uniformDistribution()*10;
	double x = 1;
	double index = 2;

	double result = McDonaldFunction(x, index);

	printf("K(%lf, %lf) = %g\n", index, x, result);
}*/

//test distribution
/*int main(){
	srand (time(NULL));

	const int binCount = 100;
	double temperature = 1E14;
	double distribution[binCount];
	double xgrid[binCount];
	xgrid[0] = 0;
	double dx = 20*kBoltzman*temperature/binCount;
	distribution[0] = 0;
	int partCount = 100000;
	for(int i = 1; i < binCount; ++i){
		xgrid[i] = xgrid[i-1] + dx;
		distribution[i] = 0;
	}

	double weight = 0;

	for(int i = 0; i < partCount; ++i){
		printf("%d\n", i);
		//double x = normalDistribution();
		//double x = maxwellDistribution(temperature);
		double x = maxwellJuttnerDistribution(temperature, massProton);
		int j = (x - xgrid[0])/dx;
		if(j > 0 && j < binCount){
			distribution[j] += 1;
		}
	}

	for(int i = 0; i < binCount; ++i){
		distribution[i] /= partCount*dx;
	}

	FILE* file = fopen("./output/test.dat","w");
	for(int i = 0; i < binCount; ++i){
		double x = xgrid[i] + dx/2;
		double maxwell = 2*sqrt(x/(pi*cube(kBoltzman*temperature)))*exp(-x/(kBoltzman*temperature));
		double theta = kBoltzman*temperature/(massProton*speed_of_light_sqr);
		double besselK = McDonaldFunction(1/theta, 2.0);
		double gamma = x/(massProton*speed_of_light_sqr);
		double juttner = maxwellJuttnerFunction(gamma, theta, besselK)/(massProton*speed_of_light_sqr);
		fprintf(file, "%g %g %g\n", x/(massProton*speed_of_light_sqr), distribution[i], juttner);
	}
	fclose(file);
}*/

//test solve inverce Juttner functiom
/*int main(){
	double x = 0.5;
	double theta = 1;
	double besselK = McDonaldFunction(1/theta, 2);
	double gamma = solveInverceJuttnerFunction(x, theta, besselK);
	double integral = maxwellJuttnerIntegral(3, theta, besselK);
	printf("x = %lf, theta = %lf, gamma = %lf\n", x, theta, gamma);
	printf("F(3) = %lf\n", integral);
	return 0;
}*/

int main()
{	
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_ALWAYS_DF);
	//omp_set_num_threads(numThreads);
	printf("start\n");
    srand (time(NULL));

	FILE* inputFile = fopen("./input/input.dat","r");
	Simulation simulation = readInput(inputFile);
	fclose(inputFile);

	simulation.simulate();

}

