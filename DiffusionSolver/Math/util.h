#ifndef UTIL_H
#define UTIL_H

double*** create3dArray(int x, int y, int z, const double& value = 0.0);
double**** create4dArray(int x, int y, int z, int l, const double& value = 0.0);

void delete3dArray(double*** u, int x, int y, int z);
void delete4dArray(double**** u, int x, int y, int z, int l);

double max(const double& a, const double& b);
double max3(const double& a, const double& b, const double& c);
double min(const double& a, const double& b);
double min3(const double& a, const double& b, const double& c);

#endif
