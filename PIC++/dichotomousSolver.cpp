//
// Created by vadim on 21.10.16.
//
#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include <string>
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "constants.h"
#include "util.h"
#include "dichotomousSolver.h"
#include "paths.h"

BaseDichotomousSolver::BaseDichotomousSolver() {
    maxErrorX = 1E-10;
    maxErrorY = 1E-10;
}

double BaseDichotomousSolver::solve(double minX, double maxX) {
    double minFunction = function(minX);
    double maxFunction = function(maxX);
    if(minFunction == 0){
        return minX;
    }
    if(maxFunction == 0){
        return maxX;
    }
    if(minFunction*maxFunction > 0){
        std::string outputDir = outputDirectory;
        FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "function must have different signes for dichotomous solver\n");
        printf("function must have different signes for dichotomous solver\n");
        fflush(stdout);
        fclose(errorLogFile);
        MPI_Finalize();
        exit(0);
    }

    double middleX = (minX + maxX)/2;

    while(fabs(maxX - minX) > maxErrorX && fabs(maxFunction - minFunction) > maxErrorY){
        middleX = (maxX + minX)/2;
        double middleFunction = function(middleX);
        if(middleFunction == 0){
            return middleX;
        }
        if(middleFunction*minFunction < 0){
            maxX = middleX;
            maxFunction = middleFunction;
        } else {
            minX = middleX;
            minFunction = middleFunction;
        }
    }

    return middleX;
}

double BaseDichotomousSolver::function(double x) {
    std::string outputDir = outputDirectory;
    FILE* errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
    fprintf(errorLogFile, "Using base class for dichotomous solver is restricted\n");
    printf("Using base class for dichotomous solver is restricted\n");
    fflush(stdout);
    fclose(errorLogFile);
    MPI_Finalize();
    exit(0);
}

double TemperatureRelativisticMaxwellSolver::function(double x){
    double McDonald1 = McDonaldFunction(x, 1);
    double McDonald2 = McDonaldFunction(x, 2);

    return 1.0/(x*(1 + ((x/alphaNormal) - 1)*(McDonald1/(x*McDonald2)))) - rightPart;
}

TemperatureRelativisticMaxwellSolver::TemperatureRelativisticMaxwellSolver(double alphaNormalValue, double rightPartValue) {
    alphaNormal = alphaNormalValue;
    rightPart = rightPartValue;
}


