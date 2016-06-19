#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>

#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "matrix3d.h"
#include "specialmath.h"
#include "output.h"

void Simulation::evaluateFields() {
    printf("evaluating fields\n");
    printLog("evaluzting fields\n");
    //fopen("./output/outputEverythingFile.dat","a");

    if (solverType == IMPLICIT) {

        evaluateMaxwellEquationMatrix();

        //maxwellMatrixFile = fopen("./output/maxwellMatrixFile.dat", "w");
        //outputMaxwellEquationMatrixFull(maxwellMatrixFile, maxwellEquationMatrix, xnumber, ynumber, znumber, maxwellEquationMatrixSize);
        //fclose(maxwellMatrixFile);
        //outputMaxwellEquationMatrixSimple(maxwellEquationMatrix, xnumber, ynumber, znumber, maxwellEquationMatrixSize);

        double ****gmresOutput = new double ***[xnumber];
        //#pragma omp parallel for
        for (int i = 0; i < xnumber; ++i) {
            gmresOutput[i] = new double **[ynumber];
            for (int j = 0; j < ynumber; ++j) {
                gmresOutput[i][j] = new double *[znumber];
                for (int k = 0; k < znumber; ++k) {
                    gmresOutput[i][j][k] = new double[maxwellEquationMatrixSize];
                }
            }
        }

        //FILE* rightPartFile = fopen("./output/rightPartFile.dat", "w");
        //for(int i = 0; i < xnumber; ++i){
        //fprintf(rightPartFile, "%28.22g %28.22g %28.22g\n", maxwellEquationRightPart[i][0][0][0], maxwellEquationRightPart[i][0][0][1], maxwellEquationRightPart[i][0][0][2]);
        //	}
        //fclose(rightPartFile);

        generalizedMinimalResidualMethod(maxwellEquationMatrix, maxwellEquationRightPart, gmresOutput, xnumber, ynumber,
                                         znumber, maxwellEquationMatrixSize, maxErrorLevel, maxGMRESIterations);
        /*conjugateGradientMethod(maxwellEquationMatrix, maxwellEquationRightPart, gmresOutput, xnumber, ynumber,
                                         znumber, maxwellEquationMatrixSize, maxErrorLevel, maxGMRESIterations);*/
        //#pragma omp parallel for

        //FILE* gmresFile = fopen("./output/gmresFile.dat", "w");
        //for(int i = 0; i < xnumber; ++i){
        //fprintf(gmresFile, "%28.22g %28.22g %28.22g\n", gmresOutput[i][0][0][0], gmresOutput[i][0][0][1], gmresOutput[i][0][0][2]);
        //}
        //fclose(gmresFile);
        for (int i = 0; i < xnumber; ++i) {
            for (int j = 0; j < ynumber; ++j) {
                for (int k = 0; k < znumber; ++k) {
                    for (int l = 0; l < 3; ++l) {
                        tempEfield[i][j][k][l] = gmresOutput[i][j][k][l];
                    }
                    delete[] gmresOutput[i][j][k];
                }
                delete[] gmresOutput[i][j];
            }
            delete[] gmresOutput[i];
        }
        delete[] gmresOutput;

        //updateBoundaries();

        double alfvenV = B0.norm() / sqrt(4 * pi * density);
        double k = 2 * pi / xsize;

        evaluateExplicitDerivative();
        //smoothEderivative();
        for (int i = 0; i < xnumber; ++i) {
            for (int j = 0; j < ynumber; ++j) {
                for (int k = 0; k < znumber; ++k) {
                    explicitEfield[i][j][k] += Ederivative[i][j][k] * deltaT;
                }
            }
        }

        //evaluateMagneticField();
        if (boundaryConditionType == PERIODIC) {
            for (int j = 0; j < ynumber; ++j) {
                for (int k = 0; k < znumber; ++k) {
                    tempEfield[xnumber][j][k] = tempEfield[0][j][k];
                    explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
                }
            }
        } else {
            //tempEfield[xnumber] = E0;

            //tempEfield[xnumber] = Efield[xnumber] + evaluateRotB(xnumber - 1)*speed_of_light_normalized*deltaT*theta;
            for (int j = 0; j < ynumber; ++j) {
                for (int k = 0; k < znumber; ++k) {
                    tempEfield[xnumber][j][k] = tempEfield[xnumber - 1][j][k];

                    explicitEfield[xnumber][j][k] = tempEfield[xnumber][j][k];
                }
            }
        }

        for (int i = 0; i < xnumber + 1; ++i) {
            for (int j = 0; j < ynumber; ++j) {
                tempEfield[i][j][znumber] = tempEfield[i][j][0];
                explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
            }
        }

        for (int i = 0; i < xnumber + 1; ++i) {
            for (int k = 0; k < znumber + 1; ++k) {
                tempEfield[i][ynumber][k] = tempEfield[i][0][k];
                explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
            }
        }

        for (int i = 0; i < xnumber + 1; ++i) {
            for (int j = 0; j < ynumber + 1; ++j) {
                for (int k = 0; k < znumber + 1; ++k) {
                    newEfield[i][j][k] = (tempEfield[i][j][k] - Efield[i][j][k] * (1 - theta)) / theta;
                    //newEfield[i][j][k].x= 0;
                }
            }
        }
    }

    if (solverType == EXPLICIT) {
        evaluateExplicitDerivative();
        //smoothEderivative();
        for (int i = 0; i < xnumber; ++i) {
            for (int j = 0; j < ynumber; ++j) {
                for (int k = 0; k < znumber; ++k) {
                    explicitEfield[i][j][k] += Ederivative[i][j][k] * deltaT;
                }
            }
        }
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
            }
        }

        for (int i = 0; i < xnumber + 1; ++i) {
            for (int k = 0; k < znumber; ++k) {
                explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
            }
        }

        for (int i = 0; i < xnumber + 1; ++i) {
            for (int j = 0; j < ynumber + 1; ++j) {
                explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
            }
        }

        for (int i = 0; i < xnumber + 1; ++i) {
            for (int j = 0; j < ynumber + 1; ++j) {
                for (int k = 0; k < znumber + 1; ++k) {
                    newEfield[i][j][k] = explicitEfield[i][j][k];
                    //newEfield[i].x = 0;
                }
            }
        }
    }
    //fclose(outputEverythingFile);
}

void Simulation::evaluateExplicitDerivative() {
    int maxX = xnumber + 1;
    if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
        maxX = xnumber;
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Ederivative[xnumber][j][k] = Vector3d(0, 0, 0);
            }
        }
    }
    int minX = 0;
    if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
        minX = 1;
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Ederivative[0][j][k] = Vector3d(0, 0, 0);
            }
        }
    }
    for (int i = minX; i < maxX; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                rotB[i][j][k] = evaluateRotB(i, j, k) * speed_of_light_normalized;
                Ederivative[i][j][k] = (evaluateRotB(i, j, k) * speed_of_light_normalized -
                                        (electricFlux[i][j][k] * 4 * pi / fieldScale));
                if (boundaryConditionType == SUPER_CONDUCTOR_LEFT) {
                    if (i == 0) {
                        Ederivative[i][j][k].y = 0;
                        Ederivative[i][j][k].z = 0;
                    }
                    if (i == xnumber) {
                        Ederivative[i][j][k] = Vector3d(0, 0, 0);
                    }
                }
            }
        }
    }
}

void Simulation::updateEfield() {
    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k] = newEfield[i][j][k];
                //Efield[i].y = newEfield[i].y;
                //Efield[i].z = newEfield[i].z;
            }
        }
    }
    if (boundaryConditionType == PERIODIC) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[xnumber][j][k] = Efield[0][j][k];
            }
        }
    } else {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[xnumber][j][k] = E0;
            }
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            Efield[i][j][znumber] = Efield[i][j][0];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int k = 0; k < znumber + 1; ++k) {
            Efield[i][ynumber][k] = Efield[i][0][k];
        }
    }
}

void Simulation::updateBfield() {
    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Bfield[i][j][k] = newBfield[i][j][k];
            }
        }
    }
}

void Simulation::updateFields() {
    updateEfield();
    updateBfield();
}

void Simulation::evaluateMaxwellEquationMatrix() {
    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                for (int l = 0; l < maxwellEquationMatrixSize; ++l) {
                    maxwellEquationMatrix[i][j][k][l].clear();
                }
            }
        }
    }

#pragma omp parallel for
    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                if ((i > 0 && i < xnumber - 1) || boundaryConditionType == PERIODIC) {
                    createInternalEquation(i, j, k);
                } else {
                    if (i == 0) {
                        createSuperConductorLeftEquation(j, k);
                    } else {
                        createFreeRightEquation(j, k);
                    }
                }
            }
        }

    }

    if (debugMode) {
        checkEquationMatrix(maxwellEquationMatrix, maxwellEquationMatrixSize);
    }
}

void Simulation::checkEquationMatrix(std::vector<MatrixElement> ****matrix, int lnumber) {
    //#pragma omp parallel for
    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                for (int l = 0; l < lnumber; ++l) {
                    for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
                        MatrixElement element = matrix[i][j][k][l][m];
                        if (element.i < 0) {
                            printf("element i < 0\n");
                            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
                            fprintf(errorLogFile, "element i = %d < 0\n", element.i);
                            fclose(errorLogFile);
                            exit(0);
                        }
                        if (element.i >= xnumber) {
                            printf("element i >= xnumber");
                            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
                            fprintf(errorLogFile, "element i = %d >= xnumber = %d\n", element.i, xnumber);
                            fclose(errorLogFile);
                            exit(0);
                        }

                        if (element.j < 0) {
                            printf("element j < 0\n");
                            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
                            fprintf(errorLogFile, "element j = %d < 0\n", element.j);
                            fclose(errorLogFile);
                            exit(0);
                        }
                        if (element.j >= ynumber) {
                            printf("element j >= ynumber");
                            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
                            fprintf(errorLogFile, "element j = %d >= ynumber = %d\n", element.j, ynumber);
                            fclose(errorLogFile);
                            exit(0);
                        }

                        if (element.k < 0) {
                            printf("element k < 0\n");
                            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
                            fprintf(errorLogFile, "element k = %d < 0\n", element.k);
                            fclose(errorLogFile);
                            exit(0);
                        }
                        if (element.k >= znumber) {
                            printf("element k >= znumber");
                            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
                            fprintf(errorLogFile, "element k = %d >= xnumber = %d\n", element.k, znumber);
                            fclose(errorLogFile);
                            exit(0);
                        }
                        for (int n = m + 1; n < matrix[i][j][k][l].size(); ++n) {
                            MatrixElement tempElement = matrix[i][j][k][l][n];

                            if (element.equalsIndex(tempElement)) {
                                printf("equals indexes\n");
                                printf("current = %d %d %d %d\n", i, j, k, l);
                                printf("temp = %d %d %d %d\n", tempElement.i, tempElement.j, tempElement.k,
                                       tempElement.l);
                                errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
                                fprintf(errorLogFile, "equal indexes current = %d %d %d %d temp = %d %d %d %d\n", i, j,
                                        k, l, element.i, element.j, element.k, element.l);
                                fclose(errorLogFile);
                                exit(0);
                            }
                        }
                    }
                }
            }
        }
    }
}

void Simulation::createSuperConductorLeftEquation(int j, int k) {
    int i = 0;

    int nextJ = j + 1;
    if (nextJ >= ynumber) {
        nextJ = 0;
    }
    int nextK = k + 1;
    if (nextK >= znumber) {
        nextK = 0;
    }
    //todo!!!

    if ((ynumber) > 1 && (znumber > 1)) {
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaX, i, j, k, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaX, i + 1, j, k, 0));

        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaX, i, nextJ, k, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaX, i + 1, nextJ, k, 0));

        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaX, i, j, nextK, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaX, i + 1, j, nextK, 0));

        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaX, i, nextJ, nextK, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaX, i + 1, nextJ, nextK, 0));

        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaY, i + 1, j, k, 1));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaY, i + 1, nextJ, k, 1));

        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaY, i + 1, j, nextK, 1));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaY, i + 1, nextJ, nextK, 1));

        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaZ, i + 1, j, k, 2));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaZ, i + 1, j, nextK, 2));

        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.25 / deltaZ, i + 1, nextJ, k, 2));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.25 / deltaZ, i + 1, nextJ, nextK, 2));

    } else if (ynumber > 1) {
        /*maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-1.0, i + 1, j, k, 0));*/
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaX, i, j, k, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaX, i + 1, j, k, 0));

        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaX, i, nextJ, k, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaX, i + 1, nextJ, k, 0));

        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaY, i + 1, j, k, 1));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaY, i + 1, nextJ, k, 1));

    } else if (znumber > 1) {
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaX, i, j, k, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaX, i + 1, j, k, 0));

        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaX, i, j, nextK, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaX, i + 1, j, nextK, 0));

        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(0.5 / deltaZ, i + 1, j, k, 2));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-0.5 / deltaZ, i + 1, j, nextK, 2));
    } else {
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0/deltaX, i, j, k, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(-1.0/deltaX, i + 1, j, k, 0));
    }

    maxwellEquationRightPart[i][j][k][0] = -4 * pi * electricDensity[0][j][k] / fieldScale;
    /*if ((ynumber == 1) && (znumber == 1)) {
        maxwellEquationRightPart[i][j][k][0] = 0;
    }*/
    maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(1.0, i, j, k, 1));
    maxwellEquationRightPart[i][j][k][1] = 0;
    maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(1.0, i, j, k, 2));
    maxwellEquationRightPart[i][j][k][2] = 0;
}

void Simulation::createFreeRightEquation(int j, int k) {
    Vector3d rightPart = E0;
    //Vector3d rightPart = Efield[xnumber - 1][j][k];
    Vector3d rightPart2 = Bfield[xnumber - 1][j][k];

    rightPart = rightPart + (evaluateRotB(xnumber - 1, j, k) * speed_of_light_normalized -
                             (electricFlux[xnumber - 1][j][k] * 4 * pi / fieldScale)) * (theta * deltaT) -
                (evaluateGradDensity(xnumber - 1, j, k) * speed_of_light_normalized_sqr * theta * theta * deltaT *
                 deltaT * 4 * pi / fieldScale);

    createFreeRightEquationX(j, k, rightPart);
    createFreeRightEquationY(j, k, rightPart);
    createFreeRightEquationZ(j, k, rightPart);

    /*maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(1.0, xnumber - 1, j, k, 0));
    maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(1.0, xnumber - 1, j, k, 1));
    maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(1.0, xnumber - 1, j, k, 2));*/

    alertNaNOrInfinity(rightPart.x, "right part x = NaN");
    alertNaNOrInfinity(rightPart.y, "right part y = NaN");
    alertNaNOrInfinity(rightPart.z, "right part z = NaN");

    maxwellEquationRightPart[xnumber - 1][j][k][0] = rightPart.x;
    maxwellEquationRightPart[xnumber - 1][j][k][1] = rightPart.y;
    maxwellEquationRightPart[xnumber - 1][j][k][2] = rightPart.z;
}

void Simulation::createFreeRightEquationX(int j, int k, Vector3d &rightPart) {
    Vector3d rightField = E0;
    double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);
    double element = 1.0 - dielectricTensor[xnumber - 1][j][k].matrix[0][0] +
                     c_theta_deltaT2 * ((2.0 - 2 * dielectricTensor[xnumber - 1][j][k].matrix[0][0]) / deltaX2);
    if (ynumber > 1) {
        element += c_theta_deltaT2 * 2.0 / deltaY2;
    }
    if (znumber > 1) {
        element += c_theta_deltaT2 * 2.0 / deltaZ2;
    }

    maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, xnumber - 1, j, k, 0));

    element = -dielectricTensor[xnumber - 1][j][k].matrix[0][1] -
              c_theta_deltaT2 * 2 * dielectricTensor[xnumber - 1][j][k].matrix[0][1] / deltaX2;
    maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, xnumber - 1, j, k, 1));

    element = -dielectricTensor[xnumber - 1][j][k].matrix[0][2] -
              c_theta_deltaT2 * 2 * dielectricTensor[xnumber - 1][j][k].matrix[0][2] / deltaX2;
    maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, xnumber - 1, j, k, 2));

    int prevI = xnumber - 2;

    int nextJ = j + 1;
    if (nextJ >= ynumber) {
        nextJ = 0;
    }
    int prevJ = j - 1;
    if (prevJ < 0) {
        prevJ = ynumber - 1;
    }

    int nextK = k + 1;
    if (nextK >= znumber) {
        nextK = 0;
    }
    int prevK = k - 1;
    if (prevK < 0) {
        prevK = znumber - 1;
    }

    element = c_theta_deltaT2 * (-1.0 + dielectricTensor[xnumber][j][k].matrix[0][0]) / deltaX2;
    rightPart.x -= element * E0.x;
    //maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
    element = c_theta_deltaT2 * dielectricTensor[xnumber][j][k].matrix[0][1] / deltaX2;
    rightPart.x -= element * E0.y;
    //maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 1));
    element = c_theta_deltaT2 * dielectricTensor[xnumber][j][k].matrix[0][2] / deltaX2;
    rightPart.x -= element * E0.z;
    //maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 2));

    element = c_theta_deltaT2 * (-1.0 + dielectricTensor[prevI][j][k].matrix[0][0]) / deltaX2;
    maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));
    element = c_theta_deltaT2 * dielectricTensor[prevI][j][k].matrix[0][1] / deltaX2;
    maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, j, k, 1));
    element = c_theta_deltaT2 * dielectricTensor[prevI][j][k].matrix[0][2] / deltaX2;
    maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, j, k, 2));

    if (ynumber > 1) {
        element = -c_theta_deltaT2 / deltaY2;
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, xnumber - 1, nextJ, k, 0));
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, xnumber - 1, prevJ, k, 0));
    }

    if (znumber > 1) {
        element = -c_theta_deltaT2 / deltaZ2;
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, xnumber - 1, j, nextK, 0));
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, xnumber - 1, j, prevK, 0));
    }

    if (ynumber > 1) {
        element = c_theta_deltaT2 * dielectricTensor[xnumber][nextJ][k].matrix[1][0] / (4 * deltaX * deltaY);
        rightPart.x -= element * E0.x;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[xnumber][nextJ][k].matrix[1][1] / (4 * deltaX * deltaY);
        rightPart.x -= element * E0.y;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber][nextJ][k].matrix[1][2] / (4 * deltaX * deltaY);
        rightPart.x -= element * E0.z;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 2));

        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[1][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[1][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[1][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 2));

        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[1][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 0));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[1][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 1));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[1][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 2));

        element = -c_theta_deltaT2 * dielectricTensor[xnumber][prevJ][k].matrix[1][0] / (4 * deltaX * deltaY);
        rightPart.x -= element * E0.x;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 0));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber][prevJ][k].matrix[1][1] / (4 * deltaX * deltaY);
        rightPart.x -= element * E0.y;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 1));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber][prevJ][k].matrix[1][2] / (4 * deltaX * deltaY);
        rightPart.x -= element * E0.z;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 2));
    }

    if (znumber > 1) {
        element = c_theta_deltaT2 * dielectricTensor[xnumber][j][nextK].matrix[2][0] / (4 * deltaX * deltaZ);
        rightPart.x -= element * E0.x;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 0));
        element = c_theta_deltaT2 * dielectricTensor[xnumber][j][nextK].matrix[2][1] / (4 * deltaX * deltaZ);
        rightPart.x -= element * E0.y;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber][j][nextK].matrix[2][2] / (4 * deltaX * deltaZ);
        rightPart.x -= element * E0.z;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 2));

        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[2][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 0));
        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[2][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 1));
        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[2][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[2][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[2][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[2][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[xnumber][j][prevK].matrix[2][0] / (4 * deltaX * deltaZ);
        rightPart.x -= element * E0.x;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber][j][prevK].matrix[2][1] / (4 * deltaX * deltaZ);
        rightPart.x -= element * E0.y;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber][j][prevK].matrix[2][2] / (4 * deltaX * deltaZ);
        rightPart.x -= element * E0.z;
        //maxwellEquationMatrix[xnumber - 1][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 2));
    }
}

void Simulation::createFreeRightEquationY(int j, int k, Vector3d &rightPart) {
    double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);
    double element = 1.0 - dielectricTensor[xnumber - 1][j][k].matrix[1][1] + c_theta_deltaT2 * (2.0 / deltaX2);
    if (ynumber > 1) {
        element += c_theta_deltaT2 * (2.0 - 2 * dielectricTensor[xnumber - 1][j][k].matrix[1][1]) / deltaY2;
    }
    if (znumber > 1) {
        element += c_theta_deltaT2 * 2.0 / deltaZ2;
    }
    maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, j, k, 1));

    element = -dielectricTensor[xnumber - 1][j][k].matrix[1][0];
    if (ynumber > 1) {
        element -= c_theta_deltaT2 * 2 * dielectricTensor[xnumber - 1][j][k].matrix[1][0] / deltaY2;
    }
    maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, j, k, 0));

    element = -dielectricTensor[xnumber - 1][j][k].matrix[1][2];
    if (ynumber > 1) {
        element -= c_theta_deltaT2 * 2 * dielectricTensor[xnumber - 1][j][k].matrix[1][2] / deltaY2;
    }
    maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, j, k, 2));

    int prevI = xnumber - 2;

    int nextJ = j + 1;
    if (nextJ >= ynumber) {
        nextJ = 0;
    }
    int prevJ = j - 1;
    if (prevJ < 0) {
        prevJ = ynumber - 1;
    }

    int nextK = k + 1;
    if (nextK >= znumber) {
        nextK = 0;
    }
    int prevK = k - 1;
    if (prevK < 0) {
        prevK = znumber - 1;
    }

    if (ynumber > 1) {
        element = c_theta_deltaT2 * (-1.0 + dielectricTensor[xnumber - 1][nextJ][k].matrix[1][1]) / deltaY2;
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, nextJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][k].matrix[1][0] / deltaY2;
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, nextJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][k].matrix[1][2] / deltaY2;
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, nextJ, k, 2));

        element = c_theta_deltaT2 * (-1.0 + dielectricTensor[xnumber - 1][prevJ][k].matrix[1][1]) / deltaY2;
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, prevJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][k].matrix[1][0] / deltaY2;
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, prevJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][k].matrix[1][2] / deltaY2;
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, prevJ, k, 2));
    }

    element = -c_theta_deltaT2 / deltaX2;
    //maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, k, 1));
    rightPart.y -= element * E0.y;
    maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, prevI, j, k, 1));

    if (znumber > 1) {
        element = -c_theta_deltaT2 / deltaZ2;
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, j, nextK, 1));
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, j, prevK, 1));
    }

    if (ynumber > 1) {
        element = c_theta_deltaT2 * dielectricTensor[xnumber][nextJ][k].matrix[0][0] / (4 * deltaX * deltaY);
        rightPart.y -= element * E0.x;
        //maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[xnumber][nextJ][k].matrix[0][1] / (4 * deltaX * deltaY);
        rightPart.y -= element * E0.y;
        //maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber][nextJ][k].matrix[0][2] / (4 * deltaX * deltaY);
        rightPart.y -= element * E0.z;
        //maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 2));

        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[0][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[0][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[0][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 2));

        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[0][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 0));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[0][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 1));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[0][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 2));

        element = -c_theta_deltaT2 * dielectricTensor[xnumber][prevJ][k].matrix[0][0] / (4 * deltaX * deltaY);
        rightPart.y -= element * E0.x;
        //maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 0));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber][prevJ][k].matrix[0][1] / (4 * deltaX * deltaY);
        rightPart.y -= element * E0.y;
        //maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 1));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber][prevJ][k].matrix[0][2] / (4 * deltaX * deltaY);
        rightPart.y -= element * E0.z;
        //maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 2));
    }

    if ((znumber > 1) && (ynumber > 1)) {
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][nextK].matrix[2][0] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, nextJ, nextK, 0));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][nextK].matrix[2][1] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, nextJ, nextK, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][nextK].matrix[2][2] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, nextJ, nextK, 2));

        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][prevK].matrix[2][0] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, prevJ, prevK, 0));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][prevK].matrix[2][1] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, prevJ, prevK, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][prevK].matrix[2][2] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, prevJ, prevK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][nextK].matrix[2][0] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, prevJ, nextK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][nextK].matrix[2][1] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, prevJ, nextK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][nextK].matrix[2][2] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, prevJ, nextK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][prevK].matrix[2][0] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, nextJ, prevK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][prevK].matrix[2][1] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, nextJ, prevK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][prevK].matrix[2][2] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][1].push_back(MatrixElement(element, xnumber - 1, nextJ, prevK, 2));
    }
}

void Simulation::createFreeRightEquationZ(int j, int k, Vector3d &rightPart) {
    double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);
    double element = 1.0 - dielectricTensor[xnumber - 1][j][k].matrix[2][2] + c_theta_deltaT2 * (2.0 / deltaX2);
    if (ynumber > 1) {
        element += c_theta_deltaT2 * 2.0 / deltaY2;
    }
    if (znumber > 1) {
        element += c_theta_deltaT2 * (2.0 - 2 * dielectricTensor[xnumber - 1][j][k].matrix[2][2]) / deltaZ2;
    }
    maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, j, k, 2));

    element = -dielectricTensor[xnumber - 1][j][k].matrix[2][0];
    if (znumber > 1) {
        element -= c_theta_deltaT2 * 2 * dielectricTensor[xnumber - 1][j][k].matrix[2][0] / deltaZ2;
    }
    maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, j, k, 0));

    element = -dielectricTensor[xnumber - 1][j][k].matrix[2][1];
    if (znumber > 1) {
        element -= c_theta_deltaT2 * 2 * dielectricTensor[xnumber - 1][j][k].matrix[2][1] / deltaZ2;
    }
    maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, j, k, 1));

    int prevI = xnumber - 2;

    int nextJ = j + 1;
    if (nextJ >= ynumber) {
        nextJ = 0;
    }
    int prevJ = j - 1;
    if (prevJ < 0) {
        prevJ = ynumber - 1;
    }

    int nextK = k + 1;
    if (nextK >= znumber) {
        nextK = 0;
    }
    int prevK = k - 1;
    if (prevK < 0) {
        prevK = znumber - 1;
    }

    if (znumber > 1) {
        element = c_theta_deltaT2 * (-1.0 + dielectricTensor[xnumber - 1][j][nextK].matrix[2][2]) / deltaZ2;
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, j, nextK, 2));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][j][nextK].matrix[2][1] / deltaZ2;
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, j, nextK, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][j][nextK].matrix[2][0] / deltaZ2;
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, j, nextK, 0));

        element = c_theta_deltaT2 * (-1.0 + dielectricTensor[xnumber - 1][j][prevK].matrix[2][2]) / deltaZ2;
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, j, prevK, 2));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][j][prevK].matrix[2][1] / deltaZ2;
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, j, prevK, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][j][prevK].matrix[2][0] / deltaZ2;
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, j, prevK, 0));
    }

    if (ynumber > 1) {
        element = -c_theta_deltaT2 / deltaY2;
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, nextJ, k, 2));
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, prevJ, k, 2));
    }

    element = -c_theta_deltaT2 / deltaX2;
    rightPart.z -= element * E0.z;
    //maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, k, 2));
    maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, prevI, j, k, 2));

    if ((ynumber > 1) && (znumber > 1)) {
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][nextK].matrix[1][0] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, nextJ, nextK, 0));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][nextK].matrix[1][1] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, nextJ, nextK, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][nextK].matrix[1][2] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, nextJ, nextK, 2));

        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][prevK].matrix[1][0] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, prevJ, prevK, 0));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][prevK].matrix[1][1] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, prevJ, prevK, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][prevK].matrix[1][2] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, prevJ, prevK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][nextK].matrix[1][0] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, prevJ, nextK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][nextK].matrix[1][1] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, prevJ, nextK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][prevJ][nextK].matrix[1][2] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, prevJ, nextK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][prevK].matrix[1][0] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, nextJ, prevK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][prevK].matrix[1][1] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, nextJ, prevK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber - 1][nextJ][prevK].matrix[1][2] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, xnumber - 1, nextJ, prevK, 2));
    }

    if (znumber > 1) {
        element = c_theta_deltaT2 * dielectricTensor[xnumber][j][nextK].matrix[0][0] / (4 * deltaX * deltaZ);
        rightPart.z -= element * E0.x;
        //maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 0));
        element = c_theta_deltaT2 * dielectricTensor[xnumber][j][nextK].matrix[0][1] / (4 * deltaX * deltaZ);
        rightPart.z -= element * E0.y;
        //maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 1));
        element = c_theta_deltaT2 * dielectricTensor[xnumber][j][nextK].matrix[0][2] / (4 * deltaX * deltaZ);
        rightPart.z -= element * E0.z;
        //maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 2));

        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[0][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 0));
        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[0][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 1));
        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[0][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[0][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[0][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[0][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[xnumber - 1][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[xnumber][j][prevK].matrix[0][0] / (4 * deltaX * deltaZ);
        rightPart.z -= element * E0.x;
        //maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber][j][prevK].matrix[0][1] / (4 * deltaX * deltaZ);
        rightPart.z -= element * E0.y;
        //maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[xnumber][j][prevK].matrix[0][2] / (4 * deltaX * deltaZ);
        rightPart.z -= element * E0.z;
        //maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 2));
    }
}

void Simulation::createInternalEquationX(int i, int j, int k) {
    double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);
    double element = 1.0 - dielectricTensor[i][j][k].matrix[0][0];
    if (xnumber > 1) {
        element += c_theta_deltaT2 * ((2.0 - 2 * dielectricTensor[i][j][k].matrix[0][0]) / deltaX2);
    }
    if (ynumber > 1) {
        element += c_theta_deltaT2 * 2.0 / deltaY2;
    }
    if (znumber > 1) {
        element += c_theta_deltaT2 * 2.0 / deltaZ2;
    }

    maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 0));

    element = -dielectricTensor[i][j][k].matrix[0][1];
    if (xnumber > 1) {
        element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[0][1] / deltaX2;
    }
    maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 1));

    element = -dielectricTensor[i][j][k].matrix[0][2];
    if (xnumber > 1) {
        element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[0][2] / deltaX2;
    }
    maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, k, 2));

    int nextI = i + 1;
    if (nextI >= xnumber) {
        nextI = 0;
    }
    int prevI = i - 1;
    if (prevI < 0) {
        prevI = xnumber - 1;
    }

    int nextJ = j + 1;
    if (nextJ >= ynumber) {
        nextJ = 0;
    }
    int prevJ = j - 1;
    if (prevJ < 0) {
        prevJ = ynumber - 1;
    }

    int nextK = k + 1;
    if (nextK >= znumber) {
        nextK = 0;
    }
    int prevK = k - 1;
    if (prevK < 0) {
        prevK = znumber - 1;
    }

    if (xnumber > 1) {
        element = c_theta_deltaT2 * (-1.0 + dielectricTensor[nextI][j][k].matrix[0][0]) / deltaX2;
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[nextI][j][k].matrix[0][1] / deltaX2;
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[nextI][j][k].matrix[0][2] / deltaX2;
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, k, 2));

        element = c_theta_deltaT2 * (-1.0 + dielectricTensor[prevI][j][k].matrix[0][0]) / deltaX2;
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[prevI][j][k].matrix[0][1] / deltaX2;
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[prevI][j][k].matrix[0][2] / deltaX2;
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, k, 2));
    }

    if (ynumber > 1) {
        element = -c_theta_deltaT2 / deltaY2;
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, nextJ, k, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, prevJ, k, 0));
    }

    if (znumber > 1) {
        element = -c_theta_deltaT2 / deltaZ2;
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, nextK, 0));
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i, j, prevK, 0));
    }

    if ((ynumber > 1) && (xnumber > 1)) {
        element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[1][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[1][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[1][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, nextJ, k, 2));

        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[1][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[1][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[1][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, prevJ, k, 2));

        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[1][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 0));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[1][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 1));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[1][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, nextJ, k, 2));

        element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[1][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 0));
        element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[1][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 1));
        element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[1][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, prevJ, k, 2));
    }

    if ((znumber > 1) && (xnumber > 1)) {
        element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[2][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 0));
        element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[2][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 1));
        element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[2][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, nextK, 2));

        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[2][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 0));
        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[2][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 1));
        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[2][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, prevK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[2][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[2][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[2][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, prevI, j, nextK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[2][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[2][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[2][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, nextI, j, prevK, 2));
    }
}

void Simulation::createInternalEquationY(int i, int j, int k) {
    double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);
    double element = 1.0 - dielectricTensor[i][j][k].matrix[1][1];
    if (xnumber > 1) {
        element += c_theta_deltaT2 * (2.0 / deltaX2);
    }
    if (ynumber > 1) {
        element += c_theta_deltaT2 * (2.0 - 2 * dielectricTensor[i][j][k].matrix[1][1]) / deltaY2;
    }
    if (znumber > 1) {
        element += c_theta_deltaT2 * 2.0 / deltaZ2;
    }
    maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 1));

    element = -dielectricTensor[i][j][k].matrix[1][0];
    if (ynumber > 1) {
        element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[1][0] / deltaY2;
    }
    maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 0));

    element = -dielectricTensor[i][j][k].matrix[1][2];
    if (ynumber > 1) {
        element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[1][2] / deltaY2;
    }
    maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, k, 2));

    int nextI = i + 1;
    if (nextI >= xnumber) {
        nextI = 0;
    }
    int prevI = i - 1;
    if (prevI < 0) {
        prevI = xnumber - 1;
    }

    int nextJ = j + 1;
    if (nextJ >= ynumber) {
        nextJ = 0;
    }
    int prevJ = j - 1;
    if (prevJ < 0) {
        prevJ = ynumber - 1;
    }

    int nextK = k + 1;
    if (nextK >= znumber) {
        nextK = 0;
    }
    int prevK = k - 1;
    if (prevK < 0) {
        prevK = znumber - 1;
    }

    if (ynumber > 1) {
        element = c_theta_deltaT2 * (-1.0 + dielectricTensor[i][nextJ][k].matrix[1][1]) / deltaY2;
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[i][nextJ][k].matrix[1][0] / deltaY2;
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[i][nextJ][k].matrix[1][2] / deltaY2;
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, k, 2));

        element = c_theta_deltaT2 * (-1.0 + dielectricTensor[i][prevJ][k].matrix[1][1]) / deltaY2;
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[i][prevJ][k].matrix[1][0] / deltaY2;
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[i][prevJ][k].matrix[1][2] / deltaY2;
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, k, 2));
    }

    if (xnumber > 1) {
        element = -c_theta_deltaT2 / deltaX2;
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, j, k, 1));
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, j, k, 1));
    }

    if (znumber > 1) {
        element = -c_theta_deltaT2 / deltaZ2;
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, nextK, 1));
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, j, prevK, 1));
    }

    if ((ynumber > 1) && (xnumber > 1)) {
        element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[0][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[0][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[nextI][nextJ][k].matrix[0][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, nextJ, k, 2));

        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[0][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 0));
        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[0][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 1));
        element = c_theta_deltaT2 * dielectricTensor[prevI][prevJ][k].matrix[0][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, prevJ, k, 2));

        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[0][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 0));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[0][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 1));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][nextJ][k].matrix[0][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, prevI, nextJ, k, 2));

        element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[0][0] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 0));
        element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[0][1] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 1));
        element = -c_theta_deltaT2 * dielectricTensor[nextI][prevJ][k].matrix[0][2] / (4 * deltaX * deltaY);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, nextI, prevJ, k, 2));
    }

    if ((znumber > 1) && (ynumber > 1)) {
        element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[2][0] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, nextK, 0));
        element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[2][1] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, nextK, 1));
        element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[2][2] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, nextK, 2));

        element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[2][0] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, prevK, 0));
        element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[2][1] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, prevK, 1));
        element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[2][2] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, prevK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[2][0] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, nextK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[2][1] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, nextK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[2][2] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, prevJ, nextK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[2][0] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, prevK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[2][1] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, prevK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[2][2] / (4 * deltaY * deltaZ);
        maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(element, i, nextJ, prevK, 2));
    }
}

void Simulation::createInternalEquationZ(int i, int j, int k) {
    double c_theta_deltaT2 = sqr(speed_of_light_normalized * theta * deltaT);
    double element = 1.0 - dielectricTensor[i][j][k].matrix[2][2];
    if (xnumber > 1) {
        element += c_theta_deltaT2 * (2.0 / deltaX2);
    }
    if (ynumber > 1) {
        element += c_theta_deltaT2 * 2.0 / deltaY2;
    }
    if (znumber > 1) {
        element += c_theta_deltaT2 * (2.0 - 2 * dielectricTensor[i][j][k].matrix[2][2]) / deltaZ2;
    }
    maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 2));

    element = -dielectricTensor[i][j][k].matrix[2][0];
    if (znumber > 1) {
        element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[2][0] / deltaZ2;
    }
    maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 0));

    element = -dielectricTensor[i][j][k].matrix[2][1];
    if (znumber > 1) {
        element -= c_theta_deltaT2 * 2 * dielectricTensor[i][j][k].matrix[2][1] / deltaZ2;
    }
    maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, k, 1));

    int nextI = i + 1;
    if (nextI >= xnumber) {
        nextI = 0;
    }
    int prevI = i - 1;
    if (prevI < 0) {
        prevI = xnumber - 1;
    }

    int nextJ = j + 1;
    if (nextJ >= ynumber) {
        nextJ = 0;
    }
    int prevJ = j - 1;
    if (prevJ < 0) {
        prevJ = ynumber - 1;
    }

    int nextK = k + 1;
    if (nextK >= znumber) {
        nextK = 0;
    }
    int prevK = k - 1;
    if (prevK < 0) {
        prevK = znumber - 1;
    }

    if (znumber > 1) {
        element = c_theta_deltaT2 * (-1.0 + dielectricTensor[i][j][nextK].matrix[2][2]) / deltaZ2;
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 2));
        element = c_theta_deltaT2 * dielectricTensor[i][j][nextK].matrix[2][1] / deltaZ2;
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 1));
        element = c_theta_deltaT2 * dielectricTensor[i][j][nextK].matrix[2][0] / deltaZ2;
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, nextK, 0));

        element = c_theta_deltaT2 * (-1.0 + dielectricTensor[i][j][prevK].matrix[2][2]) / deltaZ2;
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 2));
        element = c_theta_deltaT2 * dielectricTensor[i][j][prevK].matrix[2][1] / deltaZ2;
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 1));
        element = c_theta_deltaT2 * dielectricTensor[i][j][prevK].matrix[2][0] / deltaZ2;
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, j, prevK, 0));
    }

    if (ynumber > 1) {
        element = -c_theta_deltaT2 / deltaY2;
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, k, 2));
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, k, 2));
    }

    if (xnumber > 1) {
        element = -c_theta_deltaT2 / deltaX2;
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, k, 2));
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, k, 2));
    }

    if ((ynumber > 1) && (znumber > 1)) {
        element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[1][0] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, nextK, 0));
        element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[1][1] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, nextK, 1));
        element = c_theta_deltaT2 * dielectricTensor[i][nextJ][nextK].matrix[1][2] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, nextK, 2));

        element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[1][0] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, prevK, 0));
        element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[1][1] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, prevK, 1));
        element = c_theta_deltaT2 * dielectricTensor[i][prevJ][prevK].matrix[1][2] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, prevK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[1][0] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, nextK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[1][1] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, nextK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[i][prevJ][nextK].matrix[1][2] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, prevJ, nextK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[1][0] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, prevK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[1][1] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, prevK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[i][nextJ][prevK].matrix[1][2] / (4 * deltaZ * deltaY);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, i, nextJ, prevK, 2));
    }

    if ((znumber > 1) && (xnumber > 1)) {
        element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[0][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 0));
        element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[0][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 1));
        element = c_theta_deltaT2 * dielectricTensor[nextI][j][nextK].matrix[0][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, nextK, 2));

        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[0][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 0));
        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[0][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 1));
        element = c_theta_deltaT2 * dielectricTensor[prevI][j][prevK].matrix[0][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, prevK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[0][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[0][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[prevI][j][nextK].matrix[0][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, prevI, j, nextK, 2));

        element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[0][0] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 0));
        element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[0][1] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 1));
        element = -c_theta_deltaT2 * dielectricTensor[nextI][j][prevK].matrix[0][2] / (4 * deltaX * deltaZ);
        maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(element, nextI, j, prevK, 2));
    }
}

void Simulation::createInternalEquation(int i, int j, int k) {
    Vector3d rightPart = Efield[i][j][k];
    Vector3d rightPart2 = Bfield[i][j][k];

    //rightPart = rightPart + (evaluateRotB(i)* speed_of_light_normalized - electricFlux[i]*4*pi/fieldScale) * (theta * deltaT);
    rightPart = rightPart +
                (evaluateRotB(i, j, k) * speed_of_light_normalized - (electricFlux[i][j][k] * 4 * pi / fieldScale)) *
                (theta * deltaT) -
                (evaluateGradDensity(i, j, k) * speed_of_light_normalized_sqr * theta * theta * deltaT * deltaT * 4 *
                 pi / fieldScale);
    if (i == debugPoint) {
        //fprintf(outputEverythingFile, "rotB %d = %28.22g %28.22g %28.22g\n", debugPoint, evaluateRotB(i, j, k).x, evaluateRotB(i, j, k).y, evaluateRotB(i, j, k).z);
        //fprintf(outputEverythingFile, "j = %28.22g %28.22g %28.22g\n", electricFlux[i][j][k].x, electricFlux[i][j][k].y, electricFlux[i][j][k].z);
        //fprintf(outputEverythingFile, "grad denasty = %28.22g %28.22g %28.22g\n", evaluateGradDensity(i, j, k).x, evaluateGradDensity(i, j, k).y,evaluateGradDensity(i, j, k).z);
    }
    //rightPart = rightPart + evaluateRotB(i)*speed_of_light_normalized*theta*deltaT - electricFlux[i]*4*pi*theta*deltaT/fieldScale;
    createInternalEquationX(i, j, k);
    createInternalEquationY(i, j, k);
    createInternalEquationZ(i, j, k);

    alertNaNOrInfinity(rightPart.x, "right part x = NaN");
    alertNaNOrInfinity(rightPart.y, "right part y = NaN");
    alertNaNOrInfinity(rightPart.z, "right part z = NaN");

    alertNaNOrInfinity(rightPart2.x, "right part 2 x = NaN");
    alertNaNOrInfinity(rightPart2.y, "right part 2 y = NaN");
    alertNaNOrInfinity(rightPart2.z, "right part 2 z = NaN");


    maxwellEquationRightPart[i][j][k][0] = rightPart.x;
    maxwellEquationRightPart[i][j][k][1] = rightPart.y;
    maxwellEquationRightPart[i][j][k][2] = rightPart.z;

    //printf("right part = %g %g %g\n", rightPart.x, rightPart.y, rightPart.z);

    //maxwellEquationRightPart[i][3] = rightPart2.x;
    //maxwellEquationRightPart[i][4] = rightPart2.y;
    //maxwellEquationRightPart[i][5] = rightPart2.z;
}

void Simulation::evaluateMagneticField() {
    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Vector3d rotEold = evaluateRotE(i, j, k);
                Vector3d rotEnew = evaluateRotNewE(i, j, k);
                newBfield[i][j][k] =
                    Bfield[i][j][k] - (rotEold * (1 - theta) + rotEnew * theta) * (speed_of_light_normalized * deltaT);
                //newBfield[i][j][k] = Bfield[i][j][k] - ((evaluateRotTempE(i, j, k)) * (speed_of_light_normalized * deltaT));
            }
        }
    }
}

void Simulation::updateBoundaries() {
    if (boundaryConditionType == PERIODIC) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                tempEfield[xnumber][j][k] = tempEfield[0][j][k];
                explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
            }
        }
    } else {
        //tempEfield[xnumber] = E0;

        //tempEfield[xnumber] = Efield[xnumber] + evaluateRotB(xnumber - 1)*speed_of_light_normalized*deltaT*theta;
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                tempEfield[xnumber][j][k] = tempEfield[xnumber - 1][j][k];

                explicitEfield[xnumber][j][k] = tempEfield[xnumber][j][k];
            }
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            tempEfield[i][j][znumber] = tempEfield[i][j][0];
            explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int k = 0; k < znumber + 1; ++k) {
            tempEfield[i][ynumber][k] = tempEfield[i][0][k];
            explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
        }
    }
}

void Simulation::updateBoundariesOldField() {
    if (boundaryConditionType == PERIODIC) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Efield[xnumber][j][k] = Efield[0][j][k];
            }
        }
    } else {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Efield[xnumber][j][k] = Efield[xnumber - 1][j][k];
            }
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            Efield[i][j][znumber] = Efield[i][j][0];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int k = 0; k < znumber + 1; ++k) {
            Efield[i][ynumber][k] = Efield[i][0][k];
        }
    }
}

void Simulation::updateBoundariesNewField() {
    if (boundaryConditionType == PERIODIC) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                newEfield[xnumber][j][k] = newEfield[0][j][k];
            }
        }
    } else {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                newEfield[xnumber][j][k] = newEfield[xnumber - 1][j][k];
            }
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            newEfield[i][j][znumber] = newEfield[i][j][0];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int k = 0; k < znumber + 1; ++k) {
            newEfield[i][ynumber][k] = newEfield[i][0][k];
        }
    }
}


Vector3d Simulation::evaluateRotB(int i, int j, int k) {
    if (debugMode) {
        if ((i < 0) || ((i == 0) && (boundaryConditionType == SUPER_CONDUCTOR_LEFT))) {
            printf("i < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "i = %d < 0 in evaluateRotB\n", i);
            fclose(errorLogFile);
            exit(0);
        }

        if ((i > xnumber) || ((i == xnumber) && (boundaryConditionType == SUPER_CONDUCTOR_LEFT))) {
            printf("i >= xnumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "i = %d > xnumber = %d in evaluzteRotB\n", i, xnumber);
            fclose(errorLogFile);
            exit(0);
        }
    }

    Vector3d BrightX;
    Vector3d BleftX;
    Vector3d BrightY;
    Vector3d BleftY;
    Vector3d BrightZ;
    Vector3d BleftZ;

    int curI = i;
    if (curI >= xnumber) {
        if (boundaryConditionType == PERIODIC) {
            curI = 0;
        } else {
            curI = i - 1;
        }
    }

    int prevI = i - 1;
    if (prevI < 0) {
        if (boundaryConditionType == PERIODIC) {
            prevI = xnumber - 1;
        } else {
            prevI = i;
        }
    }

    int prevJ = j - 1;
    if (prevJ < 0) {
        prevJ = ynumber - 1;
    }

    int prevK = k - 1;
    if (prevK < 0) {
        prevK = znumber - 1;
    }
    BrightX = (Bfield[curI][j][k] + Bfield[curI][prevJ][k] + Bfield[curI][j][prevK] + Bfield[curI][prevJ][prevK]) / 4.0;
    BleftX = (Bfield[prevI][j][k] + Bfield[prevI][prevJ][k] + Bfield[prevI][j][prevK] + Bfield[prevI][prevJ][prevK]) / 4.0;

    BrightY = (Bfield[curI][j][k] + Bfield[prevI][j][k] + Bfield[curI][j][prevK] + Bfield[prevI][j][prevK]) / 4.0;
    BleftY = (Bfield[curI][prevJ][k] + Bfield[prevI][prevJ][k] + Bfield[curI][prevJ][prevK] + Bfield[prevI][prevJ][prevK]) / 4.0;

    BrightZ = (Bfield[curI][j][k] + Bfield[prevI][j][k] + Bfield[curI][prevJ][k] + Bfield[prevI][prevJ][k]) / 4.0;
    BleftZ = (Bfield[curI][j][prevK] + Bfield[prevI][j][prevK] + Bfield[curI][prevJ][prevK] + Bfield[prevI][prevJ][prevK]) / 4.0;


    double x = 0;
    double y = 0;
    double z = 0;

    x = (BrightY.z - BleftY.z) / deltaY - (BrightZ.y - BleftZ.y) / deltaZ;
    y = (BrightZ.x - BleftZ.x) / deltaZ - (BrightX.z - BleftX.z) / deltaX;
    z = (BrightX.y - BleftX.y) / deltaX - (BrightY.x - BleftY.x) / deltaY;

    return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotTempE(int i, int j, int k) {
    if (debugMode) {
        if (i < 0) {
            printf("i < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "i = %d < 0 in evaluateRotTempE\n", i);
            fclose(errorLogFile);
            exit(0);
        }

        if (i >= xnumber) {
            printf("x >= xnumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "i = %d >= xnumber = %d in evaluateRotTempE\n", i, xnumber);
            fclose(errorLogFile);
            exit(0);
        }

        if (j < 0) {
            printf("j < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "j = %d < 0 in evaluateRotTempE\n", j);
            fclose(errorLogFile);
            exit(0);
        }

        if (j >= ynumber) {
            printf("y >= ynumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "j = %d >= ynumber = %d in evaluateRotTempE\n", j, ynumber);
            fclose(errorLogFile);
            exit(0);
        }

        if (k < 0) {
            printf("k < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "k = %d < 0 in evaluateRotTempE\n", k);
            fclose(errorLogFile);
            exit(0);
        }

        if (k >= znumber) {
            printf("z >= znumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "k = %d >= znumber = %d in evaluateRotTempE\n", k, znumber);
            fclose(errorLogFile);
            exit(0);
        }
    }

    double x = 0;
    double y = 0;
    double z = 0;

    Vector3d ErightX = (tempEfield[i + 1][j][k] + tempEfield[i + 1][j + 1][k] + tempEfield[i + 1][j][k + 1] +
                        tempEfield[i + 1][j + 1][k + 1]) / 4.0;
    Vector3d EleftX =
        (tempEfield[i][j][k] + tempEfield[i][j + 1][k] + tempEfield[i][j][k + 1] + tempEfield[i][j + 1][k + 1]) / 4.0;

    Vector3d ErightY = (tempEfield[i][j + 1][k] + tempEfield[i + 1][j + 1][k] + tempEfield[i][j + 1][k + 1] +
                        tempEfield[i + 1][j + 1][k + 1]) / 4.0;
    Vector3d EleftY =
        (tempEfield[i][j][k] + tempEfield[i + 1][j][k] + tempEfield[i][j][k + 1] + tempEfield[i + 1][j][k + 1]) / 4.0;
    if (ynumber == 1) {
        ErightY = Vector3d(0, 0, 0);
        EleftY = Vector3d(0, 0, 0);
    }

    Vector3d ErightZ = (tempEfield[i][j][k + 1] + tempEfield[i + 1][j][k + 1] + tempEfield[i][j + 1][k + 1] +
                        tempEfield[i + 1][j + 1][k + 1]) / 4.0;
    Vector3d EleftZ =
        (tempEfield[i][j][k] + tempEfield[i + 1][j][k] + tempEfield[i][j + 1][k] + tempEfield[i + 1][j + 1][k]) / 4.0;
    if (znumber == 1) {
        ErightZ = Vector3d(0, 0, 0);
        EleftZ = Vector3d(0, 0, 0);
    }

    x = (ErightY.z - EleftY.z) / deltaY - (ErightZ.y - EleftZ.y) / deltaZ;
    y = (ErightZ.x - EleftZ.x) / deltaZ - (ErightX.z - EleftX.z) / deltaX;
    z = (ErightX.y - EleftX.y) / deltaX - (ErightY.x - EleftY.x) / deltaY;

    return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotE(int i, int j, int k) {
    if (debugMode) {
        if (i < 0) {
            printf("i < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "i = %d < 0 in evaluateRotTempE\n", i);
            fclose(errorLogFile);
            exit(0);
        }

        if (i >= xnumber) {
            printf("x >= xnumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "i = %d >= xnumber = %d in evaluateRotTempE\n", i, xnumber);
            fclose(errorLogFile);
            exit(0);
        }

        if (j < 0) {
            printf("j < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "j = %d < 0 in evaluateRotTempE\n", j);
            fclose(errorLogFile);
            exit(0);
        }

        if (j >= ynumber) {
            printf("y >= ynumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "j = %d >= ynumber = %d in evaluateRotTempE\n", j, ynumber);
            fclose(errorLogFile);
            exit(0);
        }

        if (k < 0) {
            printf("k < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "k = %d < 0 in evaluateRotTempE\n", k);
            fclose(errorLogFile);
            exit(0);
        }

        if (k >= znumber) {
            printf("z >= znumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "k = %d >= znumber = %d in evaluateRotTempE\n", k, znumber);
            fclose(errorLogFile);
            exit(0);
        }
    }

    double x = 0;
    double y = 0;
    double z = 0;

    Vector3d ErightX =
        (Efield[i + 1][j][k] + Efield[i + 1][j + 1][k] + Efield[i + 1][j][k + 1] + Efield[i + 1][j + 1][k + 1]) / 4.0;
    Vector3d EleftX = (Efield[i][j][k] + Efield[i][j + 1][k] + Efield[i][j][k + 1] + Efield[i][j + 1][k + 1]) / 4.0;

    Vector3d ErightY =
        (Efield[i][j + 1][k] + Efield[i + 1][j + 1][k] + Efield[i][j + 1][k + 1] + Efield[i + 1][j + 1][k + 1]) / 4.0;
    Vector3d EleftY = (Efield[i][j][k] + Efield[i + 1][j][k] + Efield[i][j][k + 1] + Efield[i + 1][j][k + 1]) / 4.0;
    if (ynumber == 1) {
        ErightY = Vector3d(0, 0, 0);
        EleftY = Vector3d(0, 0, 0);
    }

    Vector3d ErightZ =
        (Efield[i][j][k + 1] + Efield[i + 1][j][k + 1] + Efield[i][j + 1][k + 1] + Efield[i + 1][j + 1][k + 1]) / 4.0;
    Vector3d EleftZ = (Efield[i][j][k] + Efield[i + 1][j][k] + Efield[i][j + 1][k] + Efield[i + 1][j + 1][k]) / 4.0;
    if (znumber == 1) {
        ErightZ = Vector3d(0, 0, 0);
        EleftZ = Vector3d(0, 0, 0);
    }

    x = (ErightY.z - EleftY.z) / deltaY - (ErightZ.y - EleftZ.y) / deltaZ;
    y = (ErightZ.x - EleftZ.x) / deltaZ - (ErightX.z - EleftX.z) / deltaX;
    z = (ErightX.y - EleftX.y) / deltaX - (ErightY.x - EleftY.x) / deltaY;

    return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotNewE(int i, int j, int k) {
    if (debugMode) {
        if (i < 0) {
            printf("i < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "i = %d < 0 in evaluateRotTempE\n", i);
            fclose(errorLogFile);
            exit(0);
        }

        if (i >= xnumber) {
            printf("x >= xnumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "i = %d >= xnumber = %d in evaluateRotTempE\n", i, xnumber);
            fclose(errorLogFile);
            exit(0);
        }

        if (j < 0) {
            printf("j < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "j = %d < 0 in evaluateRotTempE\n", j);
            fclose(errorLogFile);
            exit(0);
        }

        if (j >= ynumber) {
            printf("y >= ynumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "j = %d >= ynumber = %d in evaluateRotTempE\n", j, ynumber);
            fclose(errorLogFile);
            exit(0);
        }

        if (k < 0) {
            printf("k < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "k = %d < 0 in evaluateRotTempE\n", k);
            fclose(errorLogFile);
            exit(0);
        }

        if (k >= znumber) {
            printf("z >= znumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "k = %d >= znumber = %d in evaluateRotTempE\n", k, znumber);
            fclose(errorLogFile);
            exit(0);
        }
    }

    double x = 0;
    double y = 0;
    double z = 0;

    Vector3d ErightX = (newEfield[i + 1][j][k] + newEfield[i + 1][j + 1][k] + newEfield[i + 1][j][k + 1] +
                        newEfield[i + 1][j + 1][k + 1]) / 4.0;
    Vector3d EleftX =
        (newEfield[i][j][k] + newEfield[i][j + 1][k] + newEfield[i][j][k + 1] + newEfield[i][j + 1][k + 1]) / 4.0;

    Vector3d ErightY = (newEfield[i][j + 1][k] + newEfield[i + 1][j + 1][k] + newEfield[i][j + 1][k + 1] +
                        newEfield[i + 1][j + 1][k + 1]) / 4.0;
    Vector3d EleftY =
        (newEfield[i][j][k] + newEfield[i + 1][j][k] + newEfield[i][j][k + 1] + newEfield[i + 1][j][k + 1]) / 4.0;
    if (ynumber == 1) {
        ErightY = Vector3d(0, 0, 0);
        EleftY = Vector3d(0, 0, 0);
    }

    Vector3d ErightZ = (newEfield[i][j][k + 1] + newEfield[i + 1][j][k + 1] + newEfield[i][j + 1][k + 1] +
                        newEfield[i + 1][j + 1][k + 1]) / 4.0;
    Vector3d EleftZ =
        (newEfield[i][j][k] + newEfield[i + 1][j][k] + newEfield[i][j + 1][k] + newEfield[i + 1][j + 1][k]) / 4.0;
    if (znumber == 1) {
        ErightZ = Vector3d(0, 0, 0);
        EleftZ = Vector3d(0, 0, 0);
    }

    x = (ErightY.z - EleftY.z) / deltaY - (ErightZ.y - EleftZ.y) / deltaZ;
    y = (ErightZ.x - EleftZ.x) / deltaZ - (ErightX.z - EleftX.z) / deltaX;
    z = (ErightX.y - EleftX.y) / deltaX - (ErightY.x - EleftY.x) / deltaY;

    return Vector3d(x, y, z);
}

double Simulation::evaluateDivE(int i, int j, int k) {
    double ErightX = (Efield[i + 1][j][k].x + Efield[i + 1][j + 1][k].x + Efield[i + 1][j][k + 1].x +
                      Efield[i + 1][j + 1][k + 1].x) / 4.0;
    double EleftX =
        (Efield[i][j][k].x + Efield[i][j + 1][k].x + Efield[i][j][k + 1].x + Efield[i][j + 1][k + 1].x) / 4.0;

    double ErightY = (Efield[i][j + 1][k].y + Efield[i + 1][j + 1][k].y + Efield[i][j + 1][k + 1].y +
                      Efield[i + 1][j + 1][k + 1].y) / 4.0;
    double EleftY =
        (Efield[i][j][k].y + Efield[i + 1][j][k].y + Efield[i][j][k + 1].y + Efield[i + 1][j][k + 1].y) / 4.0;

    double ErightZ = (Efield[i][j][k + 1].z + Efield[i + 1][j][k + 1].z + Efield[i][j + 1][k + 1].z +
                      Efield[i + 1][j + 1][k + 1].z) / 4.0;
    double EleftZ =
        (Efield[i][j][k].z + Efield[i + 1][j][k].z + Efield[i][j + 1][k].z + Efield[i + 1][j + 1][k].z) / 4.0;

    return ((ErightX - EleftX) / deltaX) + ((ErightY - EleftY) / deltaY) + ((ErightZ - EleftZ) / deltaZ);
}

double Simulation::evaluateDivCleaningE(int i, int j, int k) {
    //todo
    double ErightX = (divergenceCleaningField[i + 1][j][k][0] + divergenceCleaningField[i + 1][j + 1][k][0] +
                      divergenceCleaningField[i + 1][j][k + 1][0] + divergenceCleaningField[i + 1][j + 1][k + 1][0]) /
                     4.0;
    double EleftX = (divergenceCleaningField[i][j][k][0] + divergenceCleaningField[i][j + 1][k][0] +
                     divergenceCleaningField[i][j][k + 1][0] + divergenceCleaningField[i][j + 1][k + 1][0]) / 4.0;

    double ErightY = (divergenceCleaningField[i][j + 1][k][1] + divergenceCleaningField[i + 1][j + 1][k][1] +
                      divergenceCleaningField[i][j + 1][k + 1][1] + divergenceCleaningField[i + 1][j + 1][k + 1][1]) /
                     4.0;
    double EleftY = (divergenceCleaningField[i][j][k][1] + divergenceCleaningField[i + 1][j][k][1] +
                     divergenceCleaningField[i][j][k + 1][1] + divergenceCleaningField[i + 1][j][k + 1][1]) / 4.0;

    double ErightZ = (divergenceCleaningField[i][j][k + 1][2] + divergenceCleaningField[i + 1][j][k + 1][2] +
                      divergenceCleaningField[i][j + 1][k + 1][2] + divergenceCleaningField[i + 1][j + 1][k + 1][2]) /
                     4.0;
    double EleftZ = (divergenceCleaningField[i][j][k][2] + divergenceCleaningField[i + 1][j][k][2] +
                     divergenceCleaningField[i][j + 1][k][2] + divergenceCleaningField[i + 1][j + 1][k][2]) / 4.0;

    return (ErightX - EleftX) / deltaX + (ErightY - EleftY) / deltaY + (ErightZ - EleftZ) / deltaZ;
}

double Simulation::evaluateDivTempE(int i, int j, int k) {
    double ErightX = (tempEfield[i + 1][j][k].x + tempEfield[i + 1][j + 1][k].x + tempEfield[i + 1][j][k + 1].x +
                      tempEfield[i + 1][j + 1][k + 1].x) / 4.0;
    double EleftX = (tempEfield[i][j][k].x + tempEfield[i][j + 1][k].x + tempEfield[i][j][k + 1].x +
                     tempEfield[i][j + 1][k + 1].x) / 4.0;

    double ErightY = (tempEfield[i][j + 1][k].y + tempEfield[i + 1][j + 1][k].y + tempEfield[i][j + 1][k + 1].y +
                      tempEfield[i + 1][j + 1][k + 1].y) / 4.0;
    double EleftY = (tempEfield[i][j][k].y + tempEfield[i + 1][j][k].y + tempEfield[i][j][k + 1].y +
                     tempEfield[i + 1][j][k + 1].y) / 4.0;

    double ErightZ = (tempEfield[i][j][k + 1].z + tempEfield[i + 1][j][k + 1].z + tempEfield[i][j + 1][k + 1].z +
                      tempEfield[i + 1][j + 1][k + 1].z) / 4.0;
    double EleftZ = (tempEfield[i][j][k].z + tempEfield[i + 1][j][k].z + tempEfield[i][j + 1][k].z +
                     tempEfield[i + 1][j + 1][k].z) / 4.0;

    return ((ErightX - EleftX) / deltaX) + ((ErightY - EleftY) / deltaY) + ((ErightZ - EleftZ) / deltaZ);
}

double Simulation::evaluateDivNewE(int i, int j, int k) {
    double ErightX = (newEfield[i + 1][j][k].x + newEfield[i + 1][j + 1][k].x + newEfield[i + 1][j][k + 1].x +
                      newEfield[i + 1][j + 1][k + 1].x) / 4.0;
    double EleftX =
        (newEfield[i][j][k].x + newEfield[i][j + 1][k].x + newEfield[i][j][k + 1].x + newEfield[i][j + 1][k + 1].x) /
        4.0;

    double ErightY = (newEfield[i][j + 1][k].y + newEfield[i + 1][j + 1][k].y + newEfield[i][j + 1][k + 1].y +
                      newEfield[i + 1][j + 1][k + 1].y) / 4.0;
    double EleftY =
        (newEfield[i][j][k].y + newEfield[i + 1][j][k].y + newEfield[i][j][k + 1].y + newEfield[i + 1][j][k + 1].y) /
        4.0;

    double ErightZ = (newEfield[i][j][k + 1].z + newEfield[i + 1][j][k + 1].z + newEfield[i][j + 1][k + 1].z +
                      newEfield[i + 1][j + 1][k + 1].z) / 4.0;
    double EleftZ =
        (newEfield[i][j][k].z + newEfield[i + 1][j][k].z + newEfield[i][j + 1][k].z + newEfield[i + 1][j + 1][k].z) /
        4.0;

    return ((ErightX - EleftX) / deltaX) + ((ErightY - EleftY) / deltaY) + ((ErightZ - EleftZ) / deltaZ);
}

double Simulation::evaluateDivFlux(int i, int j, int k) {
    if (debugMode) {
        if (i < 0) {
            printf("i < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "i = %d < 0 in evaluateDivFlux\n", i);
            fclose(errorLogFile);
            exit(0);
        }

        if (i >= xnumber) {
            printf("x >= xnumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "i = %d >= xnumber = %d in evaluateDivFlux\n", i, xnumber);
            fclose(errorLogFile);
            exit(0);
        }

        if (j < 0) {
            printf("i < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "j = %d < 0 in evaluateDivFlux\n", j);
            fclose(errorLogFile);
            exit(0);
        }

        if (j >= ynumber) {
            printf("y >= ynumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "j = %d >= ynumber = %d in evaluateDivFlux\n", j, ynumber);
            fclose(errorLogFile);
            exit(0);
        }

        if (k < 0) {
            printf("k < 0\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "k = %d < 0 in evaluateDivFlux\n", k);
            fclose(errorLogFile);
            exit(0);
        }

        if (k >= znumber) {
            printf("z >= znumber\n");
            errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
            fprintf(errorLogFile, "k = %d >= znumber = %d in evaluateDivFlux\n", k, znumber);
            fclose(errorLogFile);
            exit(0);
        }
    }


    double rightFluxX =
        (electricFlux[i + 1][j][k].x + electricFlux[i + 1][j + 1][k].x + electricFlux[i + 1][j][k + 1].x +
         electricFlux[i + 1][j + 1][k + 1].x) * 0.25;
    double leftFluxX = (electricFlux[i][j][k].x + electricFlux[i][j + 1][k].x + electricFlux[i][j][k + 1].x +
                        electricFlux[i][j + 1][k + 1].x) * 0.25;

    double rightFluxY =
        (electricFlux[i][j + 1][k].x + electricFlux[i + 1][j + 1][k].x + electricFlux[i][j + 1][k + 1].x +
         electricFlux[i + 1][j + 1][k + 1].x) * 0.25;
    double leftFluxY = (electricFlux[i][j][k].x + electricFlux[i + 1][j][k].x + electricFlux[i][j][k + 1].x +
                        electricFlux[i + 1][j][k + 1].x) * 0.25;

    double rightFluxZ =
        (electricFlux[i + 1][j][k + 1].x + electricFlux[i][j][k + 1].x + electricFlux[i][j + 1][k + 1].x +
         electricFlux[i + 1][j + 1][k + 1].x) * 0.25;
    double leftFluxZ = (electricFlux[i + 1][j][k].x + electricFlux[i][j][k].x + electricFlux[i][j + 1][k].x +
                        electricFlux[i + 1][j + 1][k].x) * 0.25;

    return ((rightFluxX - leftFluxX) / deltaX) + ((rightFluxY - leftFluxY) / deltaY) +
           ((rightFluxZ - leftFluxZ) / deltaZ);
}

Vector3d Simulation::evaluateDivPressureTensor(int i, int j, int k) {
    int prevI = i - 1;
    if (prevI < 0) {
        if (boundaryConditionType == PERIODIC) {
            prevI = xnumber - 1;
        } else {
            prevI = 0;
        }
    }

    int prevJ = j - 1;
    if (prevJ < 0) {
        prevJ = ynumber - 1;
    }

    int prevK = k - 1;
    if (prevK < 0) {
        prevK = znumber - 1;
    }

    Matrix3d tensorRightX = (pressureTensor[i][j][k] + pressureTensor[i][prevJ][k] + pressureTensor[i][j][prevK] +
                             pressureTensor[i][prevJ][prevK]) * 0.25;
    Matrix3d tensorLeftX =
        (pressureTensor[prevI][j][k] + pressureTensor[prevI][prevJ][k] + pressureTensor[prevI][j][prevK] +
         pressureTensor[prevI][prevJ][prevK]) * 0.25;

    Matrix3d tensorRightY = (pressureTensor[i][j][k] + pressureTensor[prevI][j][k] + pressureTensor[i][j][prevK] +
                             pressureTensor[prevI][j][prevK]) * 0.25;
    Matrix3d tensorLeftY =
        (pressureTensor[i][prevJ][k] + pressureTensor[prevI][prevJ][k] + pressureTensor[i][prevJ][prevK] +
         pressureTensor[prevI][prevJ][prevK]) * 0.25;

    Matrix3d tensorRightZ = (pressureTensor[i][j][k] + pressureTensor[i][prevJ][k] + pressureTensor[prevI][j][k] +
                             pressureTensor[prevI][prevJ][k]) * 0.25;
    Matrix3d tensorLeftZ =
        (pressureTensor[i][j][prevK] + pressureTensor[i][prevJ][prevK] + pressureTensor[prevI][j][prevK] +
         pressureTensor[prevI][prevJ][prevK]) * 0.25;

    double x = (tensorRightX.matrix[0][0] - tensorLeftX.matrix[0][0]) / deltaX +
               (tensorRightY.matrix[1][0] - tensorRightY.matrix[1][0]) / deltaY +
               (tensorRightZ.matrix[2][0] - tensorLeftZ.matrix[2][0]) / deltaZ;
    double y = (tensorRightX.matrix[0][1] - tensorLeftX.matrix[0][1]) / deltaX +
               (tensorRightY.matrix[1][1] - tensorRightY.matrix[1][1]) / deltaY +
               (tensorRightZ.matrix[2][1] - tensorLeftZ.matrix[2][1]) / deltaZ;
    double z = (tensorRightX.matrix[0][2] - tensorLeftX.matrix[0][2]) / deltaX +
               (tensorRightY.matrix[1][2] - tensorRightY.matrix[1][2]) / deltaY +
               (tensorRightZ.matrix[2][2] - tensorLeftZ.matrix[2][2]) / deltaZ;

    return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateGradDensity(int i, int j, int k) {
    int prevI = i - 1;
    if (prevI < 0) {
        if (boundaryConditionType == PERIODIC) {
            prevI = xnumber - 1;
        } else {
            prevI = 0;
        }
    }

    int prevJ = j - 1;
    if (prevJ < 0) {
        prevJ = ynumber - 1;
    }

    int prevK = k - 1;
    if (prevK < 0) {
        prevK = znumber - 1;
    }
    /*printf("electricDensity[i][j][k] = %g\n", electricDensity[i][j][k]);
    printf("electricDensity[i][j][prevK] = %g\n", electricDensity[i][j][prevK]);
    printf("electricDensity[i][prevJ][k] = %g\n", electricDensity[i][prevJ][k]);
    printf("electricDensity[i][prevJ][prevK] = %g\n", electricDensity[i][prevJ][prevK]);
    printf("electricDensity[prevI][j][k] = %g\n", electricDensity[prevI][j][k]);
    printf("electricDensity[prevI][j][prevK] = %g\n", electricDensity[prevI][j][prevK]);
    printf("electricDensity[prevI][prevJ][k] = %g\n", electricDensity[prevI][prevJ][k]);
    printf("electricDensity[prevI][prevJ][prevK] = %g\n", electricDensity[prevI][prevJ][prevK]);*/
    double densityRightX = (electricDensity[i][j][k] + electricDensity[i][prevJ][k] + electricDensity[i][j][prevK] +
                            electricDensity[i][prevJ][prevK]) * 0.25;
    double densityLeftX =
        (electricDensity[prevI][j][k] + electricDensity[prevI][prevJ][k] + electricDensity[prevI][j][prevK] +
         electricDensity[prevI][prevJ][prevK]) * 0.25;

    double densityRightY = (electricDensity[i][j][k] + electricDensity[prevI][j][k] + electricDensity[i][j][prevK] +
                            electricDensity[prevI][j][prevK]) * 0.25;
    double densityLeftY =
        (electricDensity[i][prevJ][k] + electricDensity[prevI][prevJ][k] + electricDensity[i][prevJ][prevK] +
         electricDensity[prevI][prevJ][prevK]) * 0.25;

    double densityRightZ = (electricDensity[i][j][k] + electricDensity[i][prevJ][k] + electricDensity[prevI][j][k] +
                            electricDensity[prevI][prevJ][k]) * 0.25;
    double densityLeftZ =
        (electricDensity[i][j][prevK] + electricDensity[i][prevJ][prevK] + electricDensity[prevI][j][prevK] +
         electricDensity[prevI][prevJ][prevK]) * 0.25;

    double x = (densityRightX - densityLeftX) / deltaX;
    double y = (densityRightY - densityLeftY) / deltaY;
    double z = (densityRightZ - densityLeftZ) / deltaZ;

    return Vector3d(x, y, z);
}


