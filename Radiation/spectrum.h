#ifndef SPECTRUM_H
#define SPECTRUM_H

double evaluateMcDonaldIntegral(const double& nu);
double evaluateMcDonaldFunction5_3(const double& nu);
void evaluateMcDonaldFunctions(const double& nu, double& K1_3, double& K2_3, double& K4_3, double& K5_3);
double criticalNu(const double& E, const double& sinhi, const double& H);
double findEmissivityAt(double* nu, double* Inu, double currentNu, int Nnu);

void evaluateAllEmissivityAndAbsorption(double**** nu, double**** Inu, double**** Anu, int Nnu, double* Ee, double**** Fe, int Np, int Nd, double*** Bn, double*** sintheta, double*** psi, int*** thetaIndex, double*** concentrations, double concentration, double Bfactor, double rfactor, double a, double b, double dopplerBeta);

void evaluateSpectrumSpherical(double* nu, double**** nudoppler, double* I, double**** Inu, double**** Anu, double rmax, int Nnu, double fractionLength, const double& beta);
void evaluateSpectrumSphericalAtNu(double nu, double*** nudoppler, double& I, double*** Inu, double*** Anu, double rmax, double rfactor, double fractionLength, double d, const double& beta);

void evaluateSpectrumFlat(double* nu, double**** nudoppler, double* I, double**** Inu, double**** Anu, double rmax, int Nnu, double fractionLength, const double& beta);
void evaluateSpectrumFlatAtNu(double nu, double*** nudoppler, double& I, double*** Inu, double*** Anu, double rmax, double rfactor, double fractionLength, double d, const double& beta);

double evaluateNextdFe(double* Ee, double* dFe, double dg, int i, int Np);

void evaluateImageSpherical(double*** image, double* nu, double**** nudoppler, double**** Inu, double**** Anu, double rmax, int Nnu, double rfactor, double fractionLength, const double& beta);
void evaluateImageFlat(double*** image, double* nu, double**** nudoppler, double**** Inu, double**** Anu, double rmax, int Nnu, double rfactor, double fractionLength, const double& beta);

void evaluateNuDoppler(double***** NuDoppler, int Nmonth, int* Nnumonth, double** Numonth, const double& beta);
void evaluateNuDoppler(double**** NuDoppler, int Nnu, double* nu, const double& beta);

void updateConcentartions(double*** concentrations3d, double beta);

//for simple flat disk
void evaluateSpectrumFlatSimple(double* nu, double* totalInu, double* Inu, double* Anu, int Nnu, double rmax, double fraction);

void evaluateSpectrumAtNuSimple(double nu, double& totalInu, double Inu, double Anu, double rmax, double fraction);
void evaluateEmissivityAndAbsorptionAtNuSimple(double nu, double& Inu, double& Anu, double* Ee, double* dFe, int Np, double sinhi, double psi, double B, double concentration, double dopplerBeta, double cosbeta, double phi);
void evaluateEmissivityAndAbsorptionFlatSimple(double* nu, double* Inu, double* Anu, double* Ee, double* Fe, int Np, int Nnu, double sinhi, double psi, double B, double concentration, double dopplerBeta);
#endif