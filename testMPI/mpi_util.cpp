//
// Created by vadim on 19.06.16.
//
#include <mpi.h>
#include "mpi_util.h"



void sendGMRESTempVectorToRight(double *tempVector, double *buffer, int number) {

    buffer[0] = tempVector[number - 2];
    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rightRank = rank+1;
    if(rightRank >= size){
        rightRank = 0;
    }
    if(rightRank != rank) {
        printf("sent from %d to %d\n", rank, rightRank);
        MPI_Send(buffer, 1, MPI_DOUBLE, rightRank, 0, MPI_COMM_WORLD);
    }

}

void sendGMRESTempVectorToLeft(double *tempVector, double *buffer, int number) {
    buffer[0] = tempVector[1];

    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int leftRank = rank-1;
    if(leftRank < 0){
        leftRank = size-1;
    }

    if(leftRank != rank)
        printf("sent from %d to %d\n", rank, leftRank);
        MPI_Send(buffer, 1, MPI_DOUBLE, leftRank, 0, MPI_COMM_WORLD);


}

void receiveGMRESTempVectorFromRight(double *tempVector, double* buffer, int number) {
    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rightRank = rank+1;
    if(rightRank >= size){
        rightRank = 0;
    }
    if(rightRank != rank) {
        MPI_Status status;
        printf("recieving from %d to %d\n", rightRank, rank);
        MPI_Recv(buffer, 1, MPI_DOUBLE, rightRank, 0, MPI_COMM_WORLD,
                 &status);
        printf("recieved from %d to %d\n", rightRank, rank);

        tempVector[number-1] = buffer[0];
    } else {
        tempVector[number-1] = tempVector[1];
    }
}

void receiveGMRESTempVectorFromLeft(double *tempVector, double *buffer, int number) {
    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int leftRank = rank-1;
    if(leftRank < 0){
        leftRank = size-1;
    }
    if(leftRank != rank) {
        MPI_Status status;
        printf("recieving from %d to %d\n", leftRank, rank);
        MPI_Recv(buffer, 1, MPI_DOUBLE, leftRank, 0, MPI_COMM_WORLD,
                 &status);
        printf("recieved from %d to %d\n", leftRank, rank);

        tempVector[0] = buffer[0];
    } else {
        tempVector[0] = tempVector[number - 2];
    }
}











