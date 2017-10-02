//
// Created by vadim on 19.06.16.
//

#ifndef PIC_MPI_UTIL_H
#define PIC_MPI_UTIL_H

#include <vector>

enum MPI_TAGS{MPI_INPUT_INTEGER_TAG = 1, MPI_INPUT_DOUBLE_TAG = 2, MPI_EFIELDX_RIGHT = 3, MPI_EFIELDY_RIGHT = 4,
    MPI_EFILDZ_RIGHT = 5, MPI_EFIELDX_LEFT = 6, MPI_EFIELDY_LEFT = 7, MPI_EFILDZ_LEFT = 8,
    MPI_TEMPVECTOR_RIGHT = 9, MPI_TEMPVECTOR_LEFT = 10, MPI_SEND_DOUBLE_NUMBER_ALL = 11, MPI_SEND_GENERAL_PARAMETERS = 12,
    MPI_CELL_PARAMETERS_LEFT = 14, MPI_CELL_PARAMETERS_RIGHT = 15, MPI_NODE_PARAMETERS_LEFT = 16, MPI_NODE_PARAMETERS_RIGHT = 17,
    MPI_SEND_INTEGER_NUMBER_LEFT = 18, MPI_SEND_INTEGER_NUMBER_RIGHT = 19, MPI_SEND_DOUBLE_NUMBER_LEFT = 20, MPI_SEND_DOUBLE_NUMBER_RIGHT = 21};


void sendGMRESTempVectorToRight(double *tempVector, double *buffer, int number);
void sendGMRESTempVectorToLeft(double *tempVector, double *buffer, int number);
void receiveGMRESTempVectorFromRight(double *tempVector, double *buffer, int number);
void receiveGMRESTempVectorFromLeft(double *tempVector, double *buffer, int number);


#endif //PIC_MPI_UTIL_H
