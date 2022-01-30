#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

//new scheme// not new
double evaluateOptimizationFunction5(double* vector, double* time, double** nu, double** observedInu, double** observedError, double* Ee, double**** dFe, int Np, int* Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, double*** psi, int*** thetaIndex, double*** concentrations, double***** nudoppler, double***** Inu, double***** Anu);

//update from simple radiation
void findMinParametersGeneral(double* vector, bool* optPar, const double* grad, double* time, double** nu, double** observedInu, double** observedError, double* Ee, double**** dFe, int Np, int* Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, double*** psi, int*** thetaIndex, double*** concentrations, double***** nudoppler, double***** Inu, double***** Anu, double& currentF);
void optimizeParametersGeneral(double* vector, bool* optPar, double* time,  double** nu, double** observedInu, double** observedError, double* Ee, double**** dFe, int Np, int* Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, double*** psi, int*** thetaIndex, double*** concentrations, double***** nudoppler, double***** Inu, double***** Anu, FILE* logFile);
//void optimizeParametersGeneral(double* vector, bool* optPar, double* time, double** nu, double** observedInu, double* Ee, double**** dFe, int Np, int* Nnu, int Nd, int Nmonth, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double***** Inu, double***** Anu, FILE* logFile);

#endif