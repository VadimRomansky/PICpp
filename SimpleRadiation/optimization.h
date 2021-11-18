#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

double evaluateOptimizationFunction(double* nu, double* observedInu, double* observedError, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, const double& B, const double& R, const double& fraction, const double& epsilonB);

double evaluateOptimizationFunctionBandR(const double* vector, double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length,const double& fraction, const double& epsilonB);
void findMinParametersBandR(double* vector, const double* grad, double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double& currentF, const double& fraction, const double& epsilonB);
void optimizeParametersBandR(double* vector,  double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, const double& fraction, const double& epsilonB, FILE* logFile);

//here vector must be double[4]  B, R, fraction, epsilonB
double evaluateOptimizationFunctionGeneral(const double* vector, double* nu, double* observedInu, double* observedError, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length);
void findMinParametersGeneral(double* vector, int Npar, const double* grad, double* nu, double* observedInu, double* observedError, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, double& currentF);
void optimizeParametersGeneral(double* vector, int Npar,  double* nu, double* observedInu, double* observedError, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, FILE* logFile);

void optimizeParametersGeneral(double* vector, int Npar,  double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, FILE* logFile);

void optimizeParameterB(double& B, double* nu, double* observedInu, double* observedError, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, const double& R, const double& fraction, const double& epsilonB, FILE* logFile);
void optimizeParameterB(double& B, double* nu, double* observedInu, double* Ee, double**** dFe, int Np, int Nnu, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double**** Inu, double**** Anu, double*** area, double*** length, const double& R, const double& fraction, const double& epsilonB, FILE* logFile);

#endif