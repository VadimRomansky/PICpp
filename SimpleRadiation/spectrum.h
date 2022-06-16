#ifndef SPECTRUM_H
#define SPECTRUM_H

double evaluateMcDonaldIntegral(const double& nu);
double evaluateMcDonaldFunction5_3(const double& nu);
void evaluateMcDonaldFunctions(const double& nu, double& K1_3, double& K2_3, double& K4_3, double& K5_3);
double criticalNu(const double& E, const double& sinhi, const double& H);
double findEmissivityAt(double* nu, double* Inu, double currentNu, int Nnu);
void evaluateVolumeAndLength(double** area, double** length, double rmax, double* Rho, double* Phi, const double& fractionSize);
void evaluateSpectrum(double* nu, double* totalInu, double*** Inu, double*** Anu, double** area, double** length, int Nnu, double rmax, double* Rho, double* Phi);
void evaluateSpectrumFlatSimple(double* nu, double* totalInu, double* Inu, double* Anu, int Nnu, double rmax, double fraction);

void evaluateVolumeAndLength(double*** area, double*** length, double rmax, double fractionSize);
void evaluateAllEmissivityAndAbsorption(double* nu, double**** Inu, double**** Anu, int Nnu, double* Ee, double**** Fe, int Np, int Nd, double*** Bn, double*** sintheta, int*** thetaIndex, double*** concentrations, double concentration, double Bfactor, double rfactor, double dopplerBeta);
void evaluateSpectrum(double* nu, double* I, double**** Inu, double**** Anu, double*** area, double*** length, int Nnu, double rfactor);
void evaluateSpectrumSpherical(double* nu, double* I, double**** Inu, double**** Anu, double rmax, int Nnu, double rfactor, double fractionLength);
void evaluateSpectrumSphericalAtNu(double nu, double& I, double*** Inu, double*** Anu, double rmax, double rfactor, double fractionLength, double d);
void evaluateSpectrumFlat(double* nu, double* I, double**** Inu, double**** Anu, double rmax, int Nnu, double rfactor, double fractionLength);
void evaluateSpectrumFlatAtNu(double nu, double& I, double*** Inu, double*** Anu, double rmax, double rfactor, double fractionLength, double d);
double evaluateNextdFe(double* Ee, double* dFe, double dg, int i, int Np);

void evaluateImageSpherical(double*** image, double* nu, double**** Inu, double**** Anu, double rmax, int Nnu, double rfactor, double fractionLength);
void evaluateImageFlat(double*** image, double* nu, double**** Inu, double**** Anu, double rmax, int Nnu, double rfactor, double fractionLength);

//for simple flat disk
void evaluateSpectrumAtNuSimple(double nu, double& totalInu, double Inu, double Anu, double rmax, double fraction);
void evaluateEmissivityAndAbsorptionAtNuSimple(double nu, double& Inu, double& Anu, double* Ee, double* Fe, int Np, double sinhi, double B, double concentration, double dopplerBeta, double cosbeta);
void evaluateEmissivityAndAbsorptionFlatSimple(double* nu, double* Inu, double* Anu, double* Ee, double* Fe, int Np, int Nnu, double sinhi, double B, double concentration);
#endif