//
// Created by vadim on 19.06.16.
//

#ifndef PIC_MPI_UTIL_H
#define PIC_MPI_UTIL_H

#include <vector>
#include "particle.h"

class Simulation;
class Vector3d;
class Matrix3d;
class Particle;
class MassMatrix;

enum MPI_TAGS{MPI_INPUT_INTEGER_TAG = 1, MPI_INPUT_DOUBLE_TAG = 2, MPI_BFIELD_LEFT = 3, MPI_BFIELD_RIGHT = 4,
    MPI_EFILDZ_RIGHT = 5, MPI_EFIELDX_LEFT = 6, MPI_EFIELDY_LEFT = 7, MPI_EFILDZ_LEFT = 8,
    MPI_TEMPVECTOR_RIGHT = 9, MPI_TEMPVECTOR_LEFT = 10, MPI_SEND_DOUBLE_NUMBER_ALL = 11, MPI_SEND_GENERAL_PARAMETERS = 12,
    MPI_CELL_PARAMETERS_LEFT = 14, MPI_CELL_PARAMETERS_RIGHT = 15, MPI_NODE_PARAMETERS_LEFT = 16, MPI_NODE_PARAMETERS_RIGHT = 17,
    MPI_SEND_INTEGER_NUMBER_LEFT = 18, MPI_SEND_INTEGER_NUMBER_RIGHT = 19, MPI_SEND_DOUBLE_NUMBER_LEFT = 20, MPI_SEND_DOUBLE_NUMBER_RIGHT = 21,
    MPI_SEND_VECTOR_LEFT = 22, MPI_SEND_VECTOR_FIRST_TO_ALL = 23, MPI_SEND_VECTOR_ALL_TO_FIRST = 24, MPI_EFIELD_LEFT = 25, MPI_EFIELD_RIGHT = 26, MPI_SEND_DOUBLE_FIRST_TO_ALL = 27, MPI_SEND_DOUBLE_ALL_TO_FIRST = 30, MPI_SEND_INTEGER_ALL_TO_FIRST = 31, MPI_SEND_INTEGER_FIRST_TO_ALL = 125};

void sendInput(Simulation &simulation, int nprocs);
Simulation recieveInput(MPI_Comm comm);

void exchangeLargeVector(double**** vector, int xnumberAdded, int ynumberAdded, int znumberAdded, int lnumber, int additionalBinNumber, bool periodicX, bool
                         periodicY, bool periodicZ, MPI_Comm& cartComm, int* cartCoord, int* cartDim, double* leftOutGmresBuffer, double* rightOutGmresBuffer, double* leftInGmresBuffer, double* rightInGmresBuffer, double* frontOutGmresBuffer, double* backOutGmresBuffer, double* frontInGmresBuffer, double* backInGmresBuffer, double* bottomOutGmresBuffer, double* topOutGmresBuffer, double* bottomInGmresBuffer, double* topInGmresBuffer);

void sendLargeVectorToRightReceiveFromLeft(double ****tempVector, double* outBuffer, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void sendLargeVectorToLeftReceiveFromRight(double**** tempVector, double* outBuffer, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int lnumber, int additionalBinNumber, MPI_Comm& cartComm);

void sendLargeVectorToRight(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void sendLargeVectorToLeft(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void receiveLargeVectorFromRight(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                                 int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void receiveLargeVectorFromLeft(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                                int lnumber, int additionalBinNumber, MPI_Comm& cartComm);

void sendLargeVectorToBackReceiveFromFront(double ****tempVector, double* outBuffer, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void sendLargeVectorToFrontReceiveFromBack(double**** tempVector, double* outBuffer, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int lnumber, int additionalBinNumber, MPI_Comm& cartComm);

void sendLargeVectorToBack(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                           int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void sendLargeVectorToFront(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void receiveLargeVectorFromBack(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                                int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void receiveLargeVectorFromFront(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                                 int lnumber, int additionalBinNumber, MPI_Comm& cartComm);

void sendLargeVectorToTopReceiveFromBottom(double ****tempVector, double* outBuffer, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void sendLargeVectorToBottomReceiveFromTop(double**** tempVector, double* outBuffer, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                            int lnumber, int additionalBinNumber, MPI_Comm& cartComm);

void sendLargeVectorToTop(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                          int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void sendLargeVectorToBottom(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                             int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void receiveLargeVectorFromTop(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                               int lnumber, int additionalBinNumber, MPI_Comm& cartComm);
void receiveLargeVectorFromBottom(double ****tempVector, double *buffer, int xnumberAdded, int ynumberAdded, int znumberAdded,
                                  int lnumber, int additionalBinNumber, MPI_Comm& cartComm);

void sendCellParametersToLeftReceiveFromRight(double ***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendCellParametersToRightReceiveFromLeft(double ***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendCellParametersLeft(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);
void receiveCellParametersRight(double ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void sendCellParametersRight(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void receiveCellParametersLeft(double ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);

void sendCellParametersToFrontReceiveFromBack(double ***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendCellParametersToBackReceiveFromFront(double ***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendCellParametersFront(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);
void receiveCellParametersBack(double ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void sendCellParametersBack(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void receiveCellParametersFront(double ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);

void sendCellParametersToBottomReceiveFromTop(double ***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendCellParametersToTopReceiveFromBottom(double ***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendCellParametersBottom(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);
void receiveCellParametersTop(double ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void sendCellParametersTop(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void receiveCellParametersBottom(double ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);

void sendCellVectorParametersToLeftReceiveFromRight(Vector3d ***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendCellVectorParametersToRightReceiveFromLeft(Vector3d ***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendCellVectorParametersLeft(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);
void receiveCellVectorParametersRight(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void sendCellVectorParametersRight(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void receiveCellVectorParametersLeft(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);

void sendCellVectorParametersToFrontReceiveFromBack(Vector3d ***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendCellVectorParametersToBackReceiveFromFront(Vector3d ***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendCellVectorParametersFront(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);
void receiveCellVectorParametersBack(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void sendCellVectorParametersBack(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void receiveCellVectorParametersFront(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);

void sendCellVectorParametersToBottomReceiveFromTop(Vector3d ***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendCellVectorParametersToTopReceiveFromBottom(Vector3d ***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendCellVectorParametersBottom(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);
void receiveCellVectorParametersTop(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void sendCellVectorParametersTop(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void receiveCellVectorParametersBottom(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);

void sendCellMatrixParametersToLeftReceiveFromRight(Matrix3d ***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendCellMatrixParametersToRightReceiveFromLeft(Matrix3d ***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendCellMatrixParametersLeft(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);
void receiveCellMatrixParametersRight(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void sendCellMatrixParametersRight(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void receiveCellMatrixParametersLeft(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);

void sendCellMatrixParametersToFrontReceiveFromBack(Matrix3d ***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendCellMatrixParametersToBackReceiveFromFront(Matrix3d ***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendCellMatrixParametersFront(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);
void receiveCellMatrixParametersBack(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void sendCellMatrixParametersBack(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void receiveCellMatrixParametersFront(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);

void sendCellMatrixParametersToBottomReceiveFromTop(Matrix3d ***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendCellMatrixParametersToTopReceiveFromBottom(Matrix3d ***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendCellMatrixParametersBottom(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);
void receiveCellMatrixParametersTop(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void sendCellMatrixParametersTop(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void receiveCellMatrixParametersBottom(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);

void sendNodeParametersToLeftReceiveFromRight(double***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendNodeParametersToRightReceiveFromLeft(double***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendNodeParametersLeft(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);
void receiveNodeParametersRight(double ***array, double * inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void sendNodeParametersRight(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void receiveNodeParametersLeft(double ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);

void sendNodeParametersToFrontReceiveFromBack(double***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendNodeParametersToBackReceiveFromFront(double***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendNodeParametersFront(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);
void receiveNodeParametersBack(double ***array, double * inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void sendNodeParametersBack(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void receiveNodeParametersFront(double ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);

void sendNodeParametersToBottomReceiveFromTop(double***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendNodeParametersToTopReceiveFromBottom(double***array, double * outBuffer, double*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendNodeParametersBottom(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);
void receiveNodeParametersTop(double ***array, double * inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void sendNodeParametersTop(double ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void receiveNodeParametersBottom(double ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);

void sendNodeVectorParametersToLeftReceiveFromRight(Vector3d***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendNodeVectorParametersToRightReceiveFromLeft(Vector3d***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendNodeVectorParametersLeft(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);
void receiveNodeVectorParametersRight(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void sendNodeVectorParametersRight(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void receiveNodeVectorParametersLeft(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);

void sendNodeVectorParametersToFrontReceiveFromBack(Vector3d***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendNodeVectorParametersToBackReceiveFromFront(Vector3d***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendNodeVectorParametersFront(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);
void receiveNodeVectorParametersBack(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void sendNodeVectorParametersBack(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void receiveNodeVectorParametersFront(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);

void sendNodeVectorParametersToBottomReceiveFromTop(Vector3d***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendNodeVectorParametersToTopReceiveFromBottom(Vector3d***array, double * outBuffer, Vector3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendNodeVectorParametersBottom(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);
void receiveNodeVectorParametersTop(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void sendNodeVectorParametersTop(Vector3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void receiveNodeVectorParametersBottom(Vector3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);

void sendNodeMatrixParametersToLeftReceiveFromRight(Matrix3d***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendNodeMatrixParametersToRightReceiveFromLeft(Matrix3d***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendNodeMatrixParametersLeft(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);
void receiveNodeMatrixParametersRight(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void sendNodeMatrixParametersRight(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void receiveNodeMatrixParametersLeft(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);

void sendNodeMatrixParametersToFrontReceiveFromBack(Matrix3d***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendNodeMatrixParametersToBackReceiveFromFront(Matrix3d***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendNodeMatrixParametersFront(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);
void receiveNodeMatrixParametersBack(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void sendNodeMatrixParametersBack(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void receiveNodeMatrixParametersFront(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);

void sendNodeMatrixParametersToBottomReceiveFromTop(Matrix3d***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendNodeMatrixParametersToTopReceiveFromBottom(Matrix3d***array, double * outBuffer, Matrix3d*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendNodeMatrixParametersBottom(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);
void receiveNodeMatrixParametersTop(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void sendNodeMatrixParametersTop(Matrix3d ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void receiveNodeMatrixParametersBottom(Matrix3d ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);

void sendNodeMassMatrixParametersToLeftReceiveFromRight(MassMatrix***array, double * outBuffer, MassMatrix*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendNodeMassMatrixParametersToRightReceiveFromLeft(MassMatrix***array, double * outBuffer, MassMatrix*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);
void sendNodeMassMatrixParametersLeft(MassMatrix ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);
void receiveNodeMassMatrixParametersRight(MassMatrix ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void sendNodeMassMatrixParametersRight(MassMatrix ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int rightRank);
void receiveNodeMassMatrixParametersLeft(MassMatrix ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int leftRank);

void sendNodeMassMatrixParametersToFrontReceiveFromBack(MassMatrix***array, double * outBuffer, MassMatrix*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendNodeMassMatrixParametersToBackReceiveFromFront(MassMatrix***array, double * outBuffer, MassMatrix*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank, int backRank);
void sendNodeMassMatrixParametersFront(MassMatrix ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);
void receiveNodeMassMatrixParametersBack(MassMatrix ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void sendNodeMassMatrixParametersBack(MassMatrix ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int backRank);
void receiveNodeMassMatrixParametersFront(MassMatrix ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int frontRank);

void sendNodeMassMatrixParametersToBottomReceiveFromTop(MassMatrix***array, double * outBuffer, MassMatrix*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendNodeMassMatrixParametersToTopReceiveFromBottom(MassMatrix***array, double * outBuffer, MassMatrix*** tempArray, double* inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);
void sendNodeMassMatrixParametersBottom(MassMatrix ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);
void receiveNodeMassMatrixParametersTop(MassMatrix ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void sendNodeMassMatrixParametersTop(MassMatrix ***array, double * outBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int topRank);
void receiveNodeMassMatrixParametersBottom(MassMatrix ***array, double *inBuffer, int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalNumber, MPI_Comm& cartComm, int rank, int bottomRank);

void sendLeftReceiveRightParticles(std::vector<Particle *> &outParticles, std::vector<Particle *> &inParticles,
                                   ParticleTypeContainer *types, int typesNumber, bool periodic, int verbosity, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);

void sendRightReceiveLeftParticles(std::vector<Particle *> &outParticles, std::vector<Particle *> &inParticles,
                                   ParticleTypeContainer *types, int typesNumber, bool periodic, int verbosity, MPI_Comm& cartComm, int rank, int leftRank, int rightRank);

void sendFrontReceiveBackParticles(std::vector<Particle *> &outParticles, std::vector<Particle *> &inParticles,
                                   ParticleTypeContainer *types, int typesNumber, bool periodic, int verbosity, MPI_Comm& cartComm, int rank, int frontRank, int backRank);

void sendBackReceiveFrontParticles(std::vector<Particle *> &outParticles, std::vector<Particle *> &inParticles,
                                   ParticleTypeContainer *types, int typesNumber, bool periodic, int verbosity, MPI_Comm& cartComm, int rank, int frontRank, int backRank);

void sendBottomReceiveTopParticles(std::vector<Particle *> &outParticles, std::vector<Particle *> &inParticles,
                                   ParticleTypeContainer *types, int typesNumber, bool periodic, int verbosity, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);

void sendTopReceiveBottomParticles(std::vector<Particle *> &outParticles, std::vector<Particle *> &inParticles,
                                   ParticleTypeContainer *types, int typesNumber, bool periodic, int verbosity, MPI_Comm& cartComm, int rank, int bottomRank, int topRank);


#endif //PIC_MPI_UTIL_H
