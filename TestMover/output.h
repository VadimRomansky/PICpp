#ifndef OUTPUT_H
#define OUTPUT_H

void outputDistribution(const char* fileName, double** momentum, int N);
void writeMaxParticles(const char* fileName, double** momentum, int N);
void outputField(const char* fileNameX,const char* fileNameY,const char* fileNameZ, double*** downstreamField, double*** middleField, double*** upstreamField, int downstreamNx, int middleNx, int upstreamNx, int Ny);

#endif