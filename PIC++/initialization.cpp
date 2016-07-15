#include "stdlib.h"
#include "stdio.h"
#include <cmath>
#include <omp.h>

#include "simulation.h"
#include "specialmath.h"
#include "util.h"
#include "output.h"
#include "constants.h"
#include "matrix3d.h"
#include "random.h"

Simulation::Simulation() {
    newlyStarted = false;
    preserveChargeGlobal = true;
    outputDir = outputDirectory;

    maxEfield = Vector3d(0, 0, 0);
    maxBfield = Vector3d(0, 0, 0);
    fieldScale = 1.0;

    shockWavePoint = 0;
    chargeBalance = 0;

    Kronecker = Matrix3d(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                LeviCivita[i][j][k] = 0;
            }
        }
    }
    LeviCivita[0][1][2] = 1.0;
    LeviCivita[0][2][1] = -1.0;
    LeviCivita[1][0][2] = -1.0;
    LeviCivita[1][2][0] = 1.0;
    LeviCivita[2][0][1] = 1.0;
    LeviCivita[2][1][0] = -1.0;
    typesNumber = 6;
    types = new ParticleTypeContainer[typesNumber];
    concentrations = new double[typesNumber];
    particlesPerBin = new int[typesNumber];
}

Simulation::Simulation(int xn, int yn, int zn, double xsizev, double ysizev, double zsizev, double temp,
                       double Vx, double Vy, double Vz, double Ex, double Ey, double Ez, double Bx, double By,
                       double Bz, int maxIterations, double maxTimeV, int typesNumberV, int* particlesPerBinV, double* concentrationsV, int inType) {
    outputDir = outputDirectory;
    if (inType == 0) {
        inputType = CGS;
    } else if (inType == 1) {
        inputType = Theoretical;
    } else {
        printf("input type must be 1 or 0\n");
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "input type must be 1 or 0\n");
        fclose(errorLogFile);
        exit(0);
    }
    debugMode = true;
    newlyStarted = true;
    preserveChargeGlobal = true;
    solverType = IMPLICIT; //не явный
    //solverType = EXPLICIT; //явный
    boundaryConditionType = PERIODIC;
    //boundaryConditionType = SUPER_CONDUCTOR_LEFT;
    maxwellEquationMatrixSize = 3;

    typesNumber = typesNumberV;
    if(typesNumber != 6){
        printf("PIC++ support only 6 types of ions, typesNumber must be = 6. if you need less, ypu can initialize them with 0 concentration\n");
        fflush(stdout);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "PIC++ support only 6 types of ions, typesNumber must be = 6. if you need less, ypu can initialize them with 0 concentration\n");
        fclose(errorLogFile);
        exit(0);
    }
    concentrations = concentrationsV;
    particlesPerBin = particlesPerBinV;

    currentIteration = 0;
    time = 0;
    particlesNumber = 0;

    particleEnergy = 0;
    electricFieldEnergy = 0;
    magneticFieldEnergy = 0;
    chargeBalance = 0;

    momentum = Vector3d(0, 0, 0);


    theta = initialTheta;
    eta = theta;

    xnumber = xn;
    ynumber = yn;
    znumber = zn;

    xsize = xsizev;
    ysize = ysizev;
    zsize = zsizev;

    temperature = temp;

    maxIteration = maxIterations;
    maxTime = maxTimeV;

    extJ = 0;

    V0 = Vector3d(Vx, Vy, Vz);

    B0 = Vector3d(Bx, By, Bz);
    E0 = Vector3d(Ex, Ey, Ez);
    if (inputType == CGS) {
        massProton = massProtonReal;
        massElectron = 100 * massElectronReal;
        massAlpha = massAlphaReal;
        massDeuterium = massDeuteriumReal;
        massHelium3 = massHelium3Real;

        double concentration = concentrations[0];

        double gamma = 1 / sqrt(1 - V0.scalarMult(V0) / sqr(speed_of_light));

        double sum = 0;

        sum += particlesPerBin[0]/massElectron;
        sum += particlesPerBin[1]/massProton;
        sum += particlesPerBin[2]/massElectron;
        sum += particlesPerBin[3]/massAlpha;
        sum += particlesPerBin[4]/massDeuterium;
        sum += particlesPerBin[5]/massHelium3;


        double effectiveMass = 1 / (sum/ particlesPerBin[0]);

        //plasma_period = sqrt(effectiveMass / (4 * pi * concentration * sqr(electron_charge))) * (2 * pi) * gamma * sqrt(gamma);
        plasma_period = sqrt(effectiveMass / (4 * pi * concentration * sqr(electron_charge))) * gamma * sqrt(gamma);
        double thermal_momentum;
        if (kBoltzman * temperature > massElectron * speed_of_light * speed_of_light) {
            thermal_momentum = kBoltzman * temperature / speed_of_light;
        } else {
            thermal_momentum = sqrt(2 * massElectron * kBoltzman * temperature);
        }
        thermal_momentum += V0.norm() * massElectron;
        //scaleFactor = thermal_momentum * speed_of_light / (electron_charge * B0.norm());
        //if (B0.norm() <= 0) {
        //scaleFactor = 1.0;
        //}
        scaleFactor = speed_of_light * plasma_period;

        plasma_period = 1.0;
        scaleFactor = 1.0;

        //scaleFactor = xsize;


        E0 = E0 * (plasma_period * sqrt(scaleFactor));
        B0 = B0 * (plasma_period * sqrt(scaleFactor));
        V0 = V0 * plasma_period / scaleFactor;

        fieldScale = max2(B0.norm(), E0.norm());
        if (fieldScale <= 0) {
            fieldScale = 1.0;
        }
        fieldScale = 1.0;

        E0 = E0 / fieldScale;
        B0 = B0 / fieldScale;

        rescaleConstants();

        density = density * cube(scaleFactor);

        xsize /= scaleFactor;
        ysize /= scaleFactor;
        zsize /= scaleFactor;
        printf("xsize/scaleFactor = %lf\n", xsize);
    } else {
        massProton = 1.0;
        massElectron = massProton / 256;
        massAlpha = massProton * massAlphaReal / massProtonReal;
        massDeuterium = massProton * massDeuteriumReal / massProtonReal;
        massHelium3 = massProton * massHelium3Real / massProtonReal;
        double protonScale = massProton / massProtonReal;
        plasma_period = cube(speed_of_light) / (protonScale * electron_charge);
        scaleFactor = speed_of_light * plasma_period;
        rescaleConstantsToTheoretical();
        double densityForUnits = massElectron * (massProton + massElectron) / (4 * pi * electron_charge_normalized * electron_charge_normalized);
        if (fabs(density - densityForUnits) / (densityForUnits + densityForUnits) > 1E-3) {
            printf("density must be changed\n");
            printf("density = %g, must be = %g\n", density, densityForUnits);
        }
        fieldScale = 1.0;
    }

    maxEfield = Vector3d(0, 0, 0);
    maxBfield = Vector3d(0, 0, 0);
    shockWavePoint = 0;

    Kronecker = Matrix3d(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                LeviCivita[i][j][k] = 0;
            }
        }
    }
    LeviCivita[0][1][2] = 1.0;
    LeviCivita[0][2][1] = -1.0;
    LeviCivita[1][0][2] = -1.0;
    LeviCivita[1][2][0] = 1.0;
    LeviCivita[2][0][1] = 1.0;
    LeviCivita[2][1][0] = -1.0;
}

Simulation::~Simulation() {
    delete[] types;

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                delete[] maxwellEquationMatrix[i][j][k];
                delete[] maxwellEquationRightPart[i][j][k];
            }
            delete[] maxwellEquationMatrix[i][j];
            delete[] maxwellEquationRightPart[i][j];
        }
        delete[] maxwellEquationMatrix[i];
        delete[] maxwellEquationRightPart[i];
    }
    delete[] maxwellEquationMatrix;
    delete[] maxwellEquationRightPart;

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                delete[] divergenceCleaningField[i][j][k];
                delete[] divergenceCleaningPotential[i][j][k];
                delete[] divergenceCleanUpMatrix[i][j][k];
                delete[] divergenceCleanUpRightPart[i][j][k];
            }
            delete[] divergenceCleaningField[i][j];
            delete[] divergenceCleaningPotential[i][j];
            delete[] divergenceCleanUpMatrix[i][j];
            delete[] divergenceCleanUpRightPart[i][j];
            delete[] divergenceCleaningPotentialFourier[i][j];
        }
        delete[] divergenceCleaningField[i];
        delete[] divergenceCleaningPotential[i];
        delete[] divergenceCleanUpMatrix[i];
        delete[] divergenceCleanUpRightPart[i];
        delete[] divergenceCleaningPotentialFourier[i];
    }
    delete[] divergenceCleaningField;
    delete[] divergenceCleaningPotential;
    delete[] divergenceCleanUpMatrix;
    delete[] divergenceCleanUpRightPart;
    delete[] divergenceCleaningPotentialFourier;

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            delete[] Efield[i][j];
            delete[] newEfield[i][j];
            delete[] tempEfield[i][j];
            delete[] explicitEfield[i][j];
            delete[] rotB[i][j];
            delete[] Ederivative[i][j];

            delete[] electricFlux[i][j];
            delete[] externalElectricFlux[i][j];
            delete[] divPressureTensor[i][j];
            delete[] dielectricTensor[i][j];
        }
        delete[] Efield[i];
        delete[] newEfield[i];
        delete[] tempEfield[i];
        delete[] explicitEfield[i];
        delete[] rotB[i];
        delete[] Ederivative[i];

        delete[] electricFlux[i];
        delete[] externalElectricFlux[i];
        delete[] divPressureTensor[i];
        delete[] dielectricTensor[i];
    }

    delete[] Efield;
    delete[] newEfield;
    delete[] tempEfield;
    delete[] explicitEfield;
    delete[] rotB;
    delete[] Ederivative;

    delete[] electricFlux;
    delete[] externalElectricFlux;
    delete[] divPressureTensor;
    delete[] dielectricTensor;

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            delete[] Bfield[i][j];
            delete[] newBfield[i][j];

            delete[] electronConcentration[i][j];
            delete[] protonConcentration[i][j];
            delete[] velocityBulkProton[i][j];
            delete[] velocityBulkElectron[i][j];
            delete[] chargeDensity[i][j];
            delete[] pressureTensor[i][j];
            delete[] electricDensity[i][j];
        }
        delete[] Bfield[i];
        delete[] newBfield[i];

        delete[] electronConcentration[i];
        delete[] protonConcentration[i];
        delete[] velocityBulkProton[i];
        delete[] velocityBulkElectron[i];
        delete[] chargeDensity[i];
        delete[] pressureTensor[i];
        delete[] electricDensity[i];
    }

    delete[] Bfield;
    delete[] newBfield;

    delete[] electronConcentration;
    delete[] protonConcentration;
    delete[] velocityBulkProton;
    delete[] velocityBulkElectron;
    delete[] chargeDensity;
    delete[] pressureTensor;
    delete[] electricDensity;

    delete[] xgrid;
    delete[] middleXgrid;
    delete[] ygrid;
    delete[] middleYgrid;
    delete[] zgrid;
    delete[] middleZgrid;
}

void Simulation::rescaleConstants() {
    kBoltzman_normalized = kBoltzman * plasma_period * plasma_period / (scaleFactor * scaleFactor);
    speed_of_light_normalized = speed_of_light * plasma_period / scaleFactor;
    speed_of_light_normalized_sqr = speed_of_light_normalized * speed_of_light_normalized;
    electron_charge_normalized = electron_charge * (plasma_period / sqrt(cube(scaleFactor)));
}

void Simulation::rescaleConstantsToTheoretical() {
    //kBoltzman_normalized = kBoltzman * plasma_period * plasma_period / (scaleFactor * scaleFactor);
    kBoltzman_normalized = kBoltzman * (plasma_period * plasma_period / (scaleFactor * scaleFactor)) * massProton / massProtonReal;
    speed_of_light_normalized = 1.0;
    speed_of_light_normalized_sqr = 1.0;
    electron_charge_normalized = 1.0;
}

void Simulation::initialize() {
    printf("initialization\n");
    printLog("initialization\n");

    deltaX = xsize / (xnumber);
    deltaY = ysize / (ynumber);
    deltaZ = zsize / (znumber);

    deltaX2 = deltaX * deltaX;
    deltaY2 = deltaY * deltaY;
    deltaZ2 = deltaZ * deltaZ;

    for (int i = 0; i <= xnumber; ++i) {
        xgrid[i] = xsize + i * deltaX;
    }
    xgrid[xnumber] = 2 * xsize;

    for (int j = 0; j <= ynumber; ++j) {
        ygrid[j] = ysize + j * deltaY;
    }
    ygrid[ynumber] = 2 * ysize;

    for (int k = 0; k <= znumber; ++k) {
        zgrid[k] = zsize + k * deltaZ;
    }
    zgrid[znumber] = 2 * zsize;

    for (int i = 0; i < xnumber; ++i) {
        middleXgrid[i] = (xgrid[i] + xgrid[i + 1]) / 2;
    }

    for (int j = 0; j < ynumber; ++j) {
        middleYgrid[j] = (ygrid[j] + ygrid[j + 1]) / 2;
    }

    for (int k = 0; k < znumber; ++k) {
        middleZgrid[k] = (zgrid[k] + zgrid[k + 1]) / 2;
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k] = E0;
                newEfield[i][j][k] = Efield[i][j][k];
                tempEfield[i][j][k] = Efield[i][j][k];
                explicitEfield[i][j][k] = Efield[i][j][k];
                rotB[i][j][k] = Vector3d(0, 0, 0);
                Ederivative[i][j][k] = Vector3d(0, 0, 0);
                externalElectricFlux[i][j][k] = Vector3d(0, 0, 0);
            }
        }
    }

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Bfield[i][j][k] = B0;
                //Bfield[i][j][k].x = B0.x;
                newBfield[i][j][k] = Bfield[i][j][k];
            }
        }
    }
    createParticleTypes(concentrations, particlesPerBin);
    density = 0;
    for(int i = 0; i < typesNumber; ++i){
        density += types[i].concentration*types[i].mass;
    }

    checkDebyeParameter();

    double concentration = density / (massProton + massElectron);

    omegaPlasmaProton = sqrt(
        4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massProton);
    omegaPlasmaElectron = sqrt(
        4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massElectron);

    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    if (omegaPlasmaElectron * xsize / speed_of_light_normalized < 5) {
        printf("omegaPlasmaElectron*xsize/speed_of_light_normalized < 5\n");
        fprintf(informationFile, "omegaPlasmaElectron*xsize/speed_of_light_normalized < 5\n");
    }
    printf("omegaPlasmaElectron*xsize/speed_of_light_normalized = %g\n",
           omegaPlasmaElectron * xsize / speed_of_light_normalized);
    fprintf(informationFile, "omegaPlasmaElectron*xsize/speed_of_light_normalized = %g\n",
            omegaPlasmaElectron * xsize / speed_of_light_normalized);

    if (omegaPlasmaElectron * deltaX / speed_of_light_normalized > 1) {
        printf("omegaPlasmaElectron*deltaX/speed_of_light_normalized > 1\n");
        fprintf(informationFile, "omegaPlasmaElectron*deltaX/speed_of_light_normalized > 1\n");
    }
    printf("omegaPlasmaElectron*deltaX/speed_of_light_normalized = %g\n",
           omegaPlasmaElectron * deltaX / speed_of_light_normalized);
    fprintf(informationFile, "omegaPlasmaElectron*deltaX/speed_of_light_normalized = %g\n",
            omegaPlasmaElectron * deltaX / speed_of_light_normalized);
    fclose(informationFile);

    checkDebyeParameter();
    checkGyroRadius();

    omegaGyroProton = electron_charge * B0.norm() * fieldScale / (massProton * speed_of_light);
    omegaGyroProton = electron_charge * B0.norm() * fieldScale / (massElectron * speed_of_light);

}

void Simulation::initializeSimpleElectroMagneticWave() {
    boundaryConditionType = PERIODIC;
    E0 = Vector3d(0, 0, 0);
    B0 = Vector3d(0, 0, 0);
    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Bfield[i][j][k] = Vector3d(0, 0, 0);
                newBfield[i][j][k] = Bfield[i][j][k];
            }
        }
    }
    double kw = 2 * pi / xsize;
    double E = 1E-5;

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Efield[i][j][k].x = 0;
                Efield[i][j][k].y = E * sin(kw * xgrid[i]);
                Efield[i][j][k].z = 0;
                tempEfield[i][j][k] = Efield[i][j][k];
                newEfield[i][j][k] = tempEfield[i][j][k];
                explicitEfield[i][j][k] = Efield[i][j][k];
            }
        }
    }

    for (int k = 0; k < znumber; ++k) {
        for (int j = 0; j < ynumber; ++j) {
            Efield[xnumber][j][k] = Efield[0][j][k];
            tempEfield[xnumber][j][k] = Efield[0][j][k];
            newEfield[xnumber][j][k] = Efield[0][j][k];
            explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            Efield[i][j][znumber] = Efield[i][j][0];
            tempEfield[i][j][znumber] = Efield[i][j][0];
            newEfield[i][j][znumber] = Efield[i][j][0];
            explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
        }
    }

    for (int k = 0; k < znumber + 1; ++k) {
        for (int i = 0; i < xnumber + 1; ++i) {
            Efield[i][ynumber][k] = Efield[i][0][k];
            tempEfield[i][ynumber][k] = Efield[i][0][k];
            newEfield[i][ynumber][k] = Efield[i][0][k];
            explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
        }
    }

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Bfield[i][j][k].z = E * sin(kw * middleXgrid[i]);
                newBfield[i][j][k] = Bfield[i][j][k];
            }
        }
    }

    double t = 2 * pi / (kw * speed_of_light_normalized);
}

void Simulation::initializeRotatedSimpleElectroMagneticWave(int wavesCount) {
    boundaryConditionType = PERIODIC;

    Eyamplitude = 1;
    Ezamplitude = 0;
    Bzamplitude = Eyamplitude;
    Byamplitude = Ezamplitude;

    double kx = wavesCount * 2 * pi / xsize;
    double ky = wavesCount * 2 * pi / ysize;
    double kz = wavesCount * 2 * pi / zsize;
    kz = 0;

    double kw = sqrt(kx * kx + ky * ky + kz * kz);


    double kxy = sqrt(kx * kx + ky * ky);
    double rotatedZortNorm = sqrt(kx * kx + ky * ky + sqr(kx * kx + ky * ky) / (kz * kz));
    double matrixzz = (kx * kx + ky * ky) / (kz * rotatedZortNorm);

    Matrix3d rotationMatrix = Matrix3d(kx / kw, -ky / kxy, -kx / rotatedZortNorm,
                                       ky / kw, kx / kxy, -ky / rotatedZortNorm,
                                       kz / kw, 0, matrixzz);

    if (kz == 0) {
        rotationMatrix = Matrix3d(kx / kw, -ky / kw, 0,
                                  ky / kw, kx / kw, 0,
                                  0, 0, 1);
    }


    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k].x = 0;
                Efield[i][j][k].y = Eyamplitude * cos(kx * xgrid[i] + ky * ygrid[j] + kz * zgrid[k]);
                Efield[i][j][k].z = Ezamplitude * sin(kx * xgrid[i] + ky * ygrid[j] + kz * zgrid[k]);
                Efield[i][j][k] = rotationMatrix * Efield[i][j][k];
                explicitEfield[i][j][k] = Efield[i][j][k];
                tempEfield[i][j][k] = Efield[i][j][k];
                newEfield[i][j][k] = Efield[i][j][k];
            }
        }
    }

    for (int k = 0; k < znumber; ++k) {
        for (int j = 0; j < ynumber; ++j) {
            Efield[xnumber][j][k] = Efield[0][j][k];
            tempEfield[xnumber][j][k] = Efield[0][j][k];
            newEfield[xnumber][j][k] = Efield[0][j][k];
            explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            Efield[i][j][znumber] = Efield[i][j][0];
            tempEfield[i][j][znumber] = Efield[i][j][0];
            newEfield[i][j][znumber] = Efield[i][j][0];
            explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
        }
    }

    for (int k = 0; k < znumber + 1; ++k) {
        for (int i = 0; i < xnumber + 1; ++i) {
            Efield[i][ynumber][k] = Efield[i][0][k];
            tempEfield[i][ynumber][k] = Efield[i][0][k];
            newEfield[i][ynumber][k] = Efield[i][0][k];
            explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
        }
    }

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Bfield[i][j][k].x = 0;
                Bfield[i][j][k].y = Byamplitude * sin(kx * middleXgrid[i] + ky * middleYgrid[j] + kz * middleZgrid[k]);
                Bfield[i][j][k].z = Bzamplitude * cos(kx * middleXgrid[i] + ky * middleYgrid[j] + kz * middleZgrid[k]);
                Bfield[i][j][k] = rotationMatrix * Bfield[i][j][k];
                newBfield[i][j][k] = Bfield[i][j][k];
            }
        }
    }
}

void Simulation::initializeAlfvenWaveX(int wavesCount, double amplitudeRelation) {
    boundaryConditionType = PERIODIC;
    printf("initialization alfven wave\n");
    positronsPerBin = 0;
    alphaPerBin = 0;
    protonsPerBin = electronsPerBin;

    double concentration = density / (massProton + massElectron);
    types[1].particesDeltaX = types[0].particesDeltaX;
    types[1].particlesPerBin = types[0].particlesPerBin;
    types[0].concentration = concentration;
    types[1].concentration = concentration;
    for (int i = 2; i < typesNumber; ++i) {
        types[i].particlesPerBin = 0;
        types[i].concentration = 0;
        types[i].particesDeltaX = xsize;
    }
    createParticles();
    E0 = Vector3d(0, 0, 0);

    informationFile = fopen((outputDir + "information.dat").c_str(), "a");

    double alfvenV = B0.norm() * fieldScale / sqrt(4 * pi * density);
    if (alfvenV > speed_of_light_normalized) {
        printf("alfven velocity > c\n");
        fprintf(informationFile, "alfven velocity > c\n");
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "alfvenV/c = %15.10g > 1\n", alfvenV / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }
    fprintf(informationFile, "alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
    fprintf(informationFile, "alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);
    printf("alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
    printf("alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);

    double kw = wavesCount * 2 * pi / xsize;

    omegaPlasmaProton = sqrt(
        4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massProton);
    omegaPlasmaElectron = sqrt(
        4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
    omegaGyroProton = B0.norm() * fieldScale * electron_charge_normalized / (massProton * speed_of_light_normalized);
    omegaGyroElectron = B0.norm() * fieldScale * electron_charge_normalized / (massElectron * speed_of_light_normalized);

    if (omegaGyroProton < 5 * speed_of_light_normalized * kw) {
        printf("omegaGyroProton < 5*k*c\n");
        fprintf(informationFile, "omegaGyroProton < 5*k*c\n");
        //fclose(informationFile);
        //exit(0);
    }
    printf("omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));
    fprintf(informationFile, "omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));

    if (omegaPlasmaProton < 5 * omegaGyroProton) {
        printf("omegaPlasmaProton < 5*omegaGyroProton\n");
        fprintf(informationFile, "omegaPlasmaProton < 5*omegaGyroProton\n");
        //fclose(informationFile);
        //exit(0);
    }
    printf("omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);
    fprintf(informationFile, "omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);

    //w = q*kw*B/mP * 0.5*(sqrt(d)+-b)/a
    double b = speed_of_light_normalized * kw * (massProton - massElectron) / massProton;
    double discriminant = speed_of_light_normalized_sqr * kw * kw * sqr(
        massProton + massElectron) + 16 * pi * concentration * sqr(
        electron_charge_normalized) * (massProton + massElectron) / sqr(massProton);
    double a = (kw * kw * speed_of_light_normalized_sqr * massProton * massElectron + 4 * pi * concentration * sqr(
        electron_charge_normalized) * (massProton + massElectron)) / sqr(massProton);

    if (discriminant < 0) {
        printf("discriminant < 0\n");
        fprintf(informationFile, "discriminant < 0\n");
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "discriminant = %15.10g\n", discriminant);
        fclose(errorLogFile);
        exit(0);
    }

    double fakeOmega = (kw * electron_charge_normalized * B0.norm() * fieldScale / massProton) * (sqrt(
        discriminant) - b) / (2.0 * a);

    //a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0

    double a4 = sqr(speed_of_light_normalized_sqr * massProton * massElectron);
    a4 = a4 * sqr(sqr(sqr(fakeOmega)));
    double a3 = -2 * cube(speed_of_light_normalized_sqr) * sqr(kw * massElectron * massProton)
                - 8 * pi * concentration * sqr(
        speed_of_light_normalized_sqr * electron_charge_normalized) * massElectron * massProton * (massElectron + massProton)
                - sqr(B0.norm() * fieldScale * speed_of_light_normalized * electron_charge_normalized) * (sqr(
        massProton) + sqr(massElectron));
    a3 = a3 * cube(sqr(fakeOmega));
    double a2 = sqr(sqr(speed_of_light_normalized_sqr * kw) * massProton * massElectron)
                + 8 * pi * cube(speed_of_light_normalized_sqr) * concentration * sqr(
        kw * electron_charge_normalized) * massProton * massElectron * (massProton + massElectron)
                + 16 * sqr(pi * speed_of_light_normalized_sqr * concentration * sqr(
        electron_charge_normalized) * (massProton + massElectron))
                + 2 * sqr(
        B0.norm() * fieldScale * speed_of_light_normalized_sqr * kw * electron_charge_normalized) * (sqr(
        massProton) + sqr(massElectron))
                + 8 * pi * concentration * sqr(B0.norm() * fieldScale * speed_of_light_normalized * sqr(
        electron_charge_normalized)) * (massProton + massElectron)
                + sqr(sqr(B0.norm() * fieldScale * electron_charge_normalized));
    a2 = a2 * sqr(sqr(fakeOmega));
    double a1 = -sqr(
        B0.norm() * fieldScale * cube(speed_of_light_normalized) * kw * kw * electron_charge_normalized) * (sqr(
        massProton) + sqr(massElectron))
                - 8 * pi * concentration * sqr(B0.norm() * fieldScale * speed_of_light_normalized_sqr * kw * sqr(
        electron_charge_normalized)) * (massProton + massElectron)
                - 2 * sqr(speed_of_light_normalized * kw * sqr(B0.norm() * fieldScale * electron_charge_normalized));
    a1 = a1 * sqr(fakeOmega);
    double a0 = sqr(sqr(B0.norm() * fieldScale * speed_of_light_normalized * kw * electron_charge_normalized));

    a4 = a4 / a0;
    a3 = a3 / a0;
    a2 = a2 / a0;
    a1 = a1 / a0;
    a0 = 1.0;

    printf("a4 = %g\n", a4);
    fprintf(informationFile, "a4 = %g\n", a4);
    printf("a3 = %g\n", a3);
    fprintf(informationFile, "a3 = %g\n", a3);
    printf("a2 = %g\n", a2);
    fprintf(informationFile, "a2 = %g\n", a2);
    printf("a1 = %g\n", a1);
    fprintf(informationFile, "a1 = %g\n", a1);
    printf("a0 = %g\n", a0);
    fprintf(informationFile, "a0 = %g\n", a0);

    double fakeOmega1 = kw * alfvenV;
    printf("fakeOmega = %g\n", fakeOmega1 / plasma_period);
    fprintf(informationFile, "fakeOmega = %g\n", fakeOmega1 / plasma_period);
    double realOmega2 = solve4orderEquation(a4, a3, a2, a1, a0, 1.0);
    if (realOmega2 < 0) {
        printf("omega^2 < 0\n");
        fprintf(informationFile, "omega^2 < 0\n");
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "omega^2 = %15.10g > 1\n", realOmega2);
        fclose(errorLogFile);
        exit(0);
    }

    double error = (((a4 * realOmega2 + a3) * realOmega2 + a2) * realOmega2 + a1) * realOmega2 + a0;
    printf("error = %15.10g\n", error);
    fprintf(informationFile, "error = %15.10g\n", error);
    //double
    omega = sqrt(realOmega2) * fakeOmega;
    if (omega < 0) {
        omega = -omega;
    }

    if (omega > speed_of_light_normalized * kw / 5.0) {
        printf("omega > k*c/5\n");
        fprintf(informationFile, "omega > k*c/5\n");
        printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
        fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
        //fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "omega/kc = %15.10g > 0.2\n", omega / (kw * speed_of_light_normalized));
        fclose(errorLogFile);
        //exit(0);
    }
    printf("omega = %g\n", omega / plasma_period);
    fprintf(informationFile, "omega = %g\n", omega / plasma_period);

    printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
    fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));

    if (fabs(omega) > omegaGyroProton / 2) {
        printf("omega > omegaGyroProton/2\n");
        fprintf(informationFile, "omega > omegaGyroProton/2\n");
    }
    printf("omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
    fprintf(informationFile, "omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
    fclose(informationFile);

    checkFrequency(omega);

    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    //checkCollisionTime(omega);
    //checkMagneticReynolds(alfvenV);
    //checkDissipation(kw, alfvenV);

    double epsilonAmplitude = amplitudeRelation;

    double alfvenVReal = omega / kw;

    //double
    Bzamplitude = B0.norm() * epsilonAmplitude;

    double Omegae = omegaGyroElectron;
    double Omegae2 = Omegae * Omegae;
    double Omegap = omegaGyroProton;
    double Omegap2 = Omegap * Omegap;
    double omegae = omegaPlasmaElectron;
    double omegae2 = omegae * omegae;
    double omegap = omegaPlasmaProton;
    double omegap2 = omegap * omegap;

    double kc = kw * speed_of_light_normalized;
    double kc2 = kc * kc;

    double denominator = omega * omega - Omegae2 - (omegae2 * omega * omega / (omega * omega - kc2));

    //double
    VzamplitudeProton = -((1.0 / (4 * pi * concentration * electron_charge_normalized)) * (kc + ((omegae2 + omegap2 - omega * omega) / kc) + (omegae2 * Omegae2 / (kc * denominator))) / ((Omegae * omegae2 * omega / ((kc2 - omega * omega) * denominator)) + (Omegap / omega))) * Bzamplitude * fieldScale;
    //double
    VzamplitudeElectron = (((electron_charge_normalized * omega * Omegae) / (massElectron * kc)) * Bzamplitude * fieldScale + (omegae2 * omega * omega / (kc2 - omega * omega)) * VzamplitudeProton) / denominator;

    //double
    Byamplitude = (4 * pi * concentration * electron_charge_normalized / ((omega * omega / kc) - kc)) * (VzamplitudeElectron - VzamplitudeProton) / fieldScale;

    //double
    VyamplitudeProton = -(Omegap / omega) * VzamplitudeProton - (electron_charge_normalized / (massProton * kc)) * Bzamplitude * fieldScale;
    //double
    VyamplitudeElectron = (Omegae / omega) * VzamplitudeElectron + (electron_charge_normalized / (massElectron * kc)) * Bzamplitude * fieldScale;

    //double
    Eyamplitude = (omega / kc) * Bzamplitude;
    //double
    Ezamplitude = -(omega / kc) * Byamplitude;

    double xshift = 0;
    //double xshift = xsize/4;

    //Eyamplitude = 0.0;
    //VzamplitudeElectron = 0.0;
    //VzamplitudeProton = 0.0;
    //Byamplitude = 0.0;

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k].x = 0;
                Efield[i][j][k].y = Eyamplitude * cos(kw * xgrid[i] - kw * xshift);
                Efield[i][j][k].z = Ezamplitude * sin(kw * xgrid[i] - kw * xshift);
                explicitEfield[i][j][k] = Efield[i][j][k];
                tempEfield[i][j][k] = Efield[i][j][k];
                newEfield[i][j][k] = Efield[i][j][k];
            }
        }
    }

    for (int k = 0; k < znumber; ++k) {
        for (int j = 0; j < ynumber; ++j) {
            Efield[xnumber][j][k] = Efield[0][j][k];
            tempEfield[xnumber][j][k] = Efield[0][j][k];
            newEfield[xnumber][j][k] = Efield[0][j][k];
            explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            Efield[i][j][znumber] = Efield[i][j][0];
            tempEfield[i][j][znumber] = Efield[i][j][0];
            newEfield[i][j][znumber] = Efield[i][j][0];
            explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
        }
    }

    for (int k = 0; k < znumber + 1; ++k) {
        for (int i = 0; i < xnumber + 1; ++i) {
            Efield[i][ynumber][k] = Efield[i][0][k];
            tempEfield[i][ynumber][k] = Efield[i][0][k];
            newEfield[i][ynumber][k] = Efield[i][0][k];
            explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
        }
    }

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Bfield[i][j][k].x = B0.norm();
                Bfield[i][j][k].y = Byamplitude * sin(kw * middleXgrid[i] - kw * xshift);
                Bfield[i][j][k].z = Bzamplitude * cos(kw * middleXgrid[i] - kw * xshift);
                newBfield[i][j][k] = Bfield[i][j][k];
            }
        }
    }

    if (fabs(VzamplitudeProton) > speed_of_light_normalized) {
        printf("VzamplitudeProton > speed_of_light_normalized\n");
        fprintf(informationFile, "VzamplitudeProton > speed_of_light_normalized\n");
        printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
        fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VzamplitudeProton/c = %15.10g > 1\n", VzamplitudeProton / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }
    printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
    fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);

    if (fabs(VzamplitudeElectron) > speed_of_light_normalized) {
        printf("VzamplitudeElectron > speed_of_light_normalized\n");
        fprintf(informationFile, "VzamplitudeElectron > speed_of_light_normalized\n");
        printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
        fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VzamplitudeElectron/c = %15.10g > 1\n", VzamplitudeElectron / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }
    printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
    fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);

    if (fabs(VyamplitudeProton) > speed_of_light_normalized) {
        printf("VyamplitudeProton > speed_of_light_normalized\n");
        fprintf(informationFile, "VyamplitudeProton > speed_of_light_normalized\n");
        printf("VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
        fprintf(informationFile, "VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VyamplitudeProton/c = %15.10g > 1\n", VyamplitudeProton / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }

    if (fabs(VyamplitudeElectron) > speed_of_light_normalized) {
        printf("VyamplitudeElectron > speed_of_light_normalized\n");
        fprintf(informationFile, "VyamplitudeElectron > speed_of_light_normalized\n");
        printf("VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
        fprintf(informationFile, "VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VyamplitudeElectron/c = %15.10g > 1\n", VyamplitudeElectron / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }

    //k > 0, w > 0, kx-wt, Bz > 0, Vz < 0, Vyp > 0, Vye < 0

    for (int pcount = 0; pcount < particles.size(); ++pcount) {
        Particle *particle = particles[pcount];
        Vector3d velocity = particle->velocity(speed_of_light_normalized);
        int xn = (particle->coordinates.x - xgrid[0]) / deltaX;
        double rightWeight = (particle->coordinates.x - xgrid[xn]) / deltaX;
        double leftWeight = (xgrid[xn + 1] - particle->coordinates.x) / deltaX;
        if (particle->type == PROTON) {
            velocity = (Vector3d(0, 0, 1) * (VzamplitudeProton) * cos(kw * particle->coordinates.x - kw * xshift) + Vector3d(0, 1, 0) * VyamplitudeProton * sin(kw * particle->coordinates.x - kw * xshift));
            //velocity = Vector3d(0, 0, 1) * (VzamplitudeProton) * (leftWeight*(cos(kw * (xgrid[xn] - xshift)) + rightWeight*cos(kw*(xgrid[xn+1] - xshift)))) + Vector3d(0, 1, 0) * VyamplitudeProton * (leftWeight*(sin(kw * (xgrid[xn] - xshift)) + rightWeight*sin(kw*(xgrid[xn+1] - xshift))));

        }
        if (particle->type == ELECTRON) {
            velocity = (Vector3d(0, 0, 1) * (VzamplitudeElectron) * cos(kw * particle->coordinates.x - kw * xshift) + Vector3d(0, 1, 0) * VyamplitudeElectron * sin(kw * particle->coordinates.x - kw * xshift));
            //velocity = Vector3d(0, 0, 1) * (VzamplitudeElectron) * (leftWeight*(cos(kw * (xgrid[xn] - xshift)) + rightWeight*cos(kw*(xgrid[xn+1] - xshift)))) + Vector3d(0, 1, 0) * VyamplitudeElectron * (leftWeight*(sin(kw * (xgrid[xn] - xshift)) + rightWeight*sin(kw*(xgrid[xn+1] - xshift))));
        }
        double beta = velocity.norm() / speed_of_light_normalized;
        particle->addVelocity(velocity, speed_of_light_normalized);
        particle->initialMomentum = particle->momentum;
    }

    updateDeltaT();

    printf("dt/Talfven = %g\n", deltaT * omega / (2 * pi));
    printf("dt = %g\n", deltaT * plasma_period);
    fprintf(informationFile, "dt/Talfven = %g\n", deltaT * omega / (2 * pi));
    fprintf(informationFile, "dt = %g\n", deltaT * plasma_period);

    double Vthermal = sqrt(2 * kBoltzman_normalized * temperature / massElectron);
    double thermalFlux = Vthermal * concentration * electron_charge_normalized / sqrt(1.0 * electronsPerBin);
    double alfvenFlux = (VyamplitudeProton - VyamplitudeElectron) * concentration * electron_charge_normalized;
    if (thermalFlux > alfvenFlux / 2) {
        printf("thermalFlux > alfvenFlux/2\n");
        fprintf(informationFile, "thermalFlux > alfvenFlux/2\n");
    }
    printf("alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
    fprintf(informationFile, "alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
    double minDeltaT = deltaX / Vthermal;
    if (minDeltaT > deltaT) {
        printf("deltaT < dx/Vthermal\n");
        fprintf(informationFile, "deltaT < dx/Vthermal\n");

        //printf("deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);
        //fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);

        //fclose(informationFile);
        //exit(0);
    }
    printf("deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
    fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
    fprintf(informationFile, "\n");

    fprintf(informationFile, "Bz amplitude = %g\n", Bzamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "By amplitude = %g\n", Byamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "Vz amplitude p = %g\n", VzamplitudeProton * scaleFactor / plasma_period);
    fprintf(informationFile, "Vz amplitude e = %g\n", VzamplitudeElectron * scaleFactor / plasma_period);
    fprintf(informationFile, "Vy amplitude p = %g\n", VyamplitudeProton * scaleFactor / plasma_period);
    fprintf(informationFile, "Vy amplitude e = %g\n", VyamplitudeElectron * scaleFactor / plasma_period);
    fprintf(informationFile, "Ey amplitude = %g\n", Eyamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "Ez amplitude = %g\n", Ezamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "By/Ez = %g\n", Byamplitude / Ezamplitude);
    fprintf(informationFile, "Bz/Ey = %g\n", Bzamplitude / Eyamplitude);
    fprintf(informationFile, "4*pi*Jy amplitude = %g\n",
            4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "c*rotBy amplitude = %g\n",
            speed_of_light_normalized * kw * Bzamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "4*pi*Jz amplitude = %g\n",
            4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "c*rotBz amplitude = %g\n",
            speed_of_light_normalized * kw * Byamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative By amplitude = %g\n",
            -omega * Byamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "-c*rotEy = %g\n",
            speed_of_light_normalized * kw * Ezamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative Bz amplitude = %g\n",
            omega * Bzamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "-c*rotEz = %g\n",
            speed_of_light_normalized * kw * Eyamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative Ey amplitude = %g\n",
            omega * Eyamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "c*rotBy - 4*pi*Jy = %g\n",
            (speed_of_light_normalized * kw * Bzamplitude * fieldScale - 4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative Ez amplitude = %g\n",
            -omega * Ezamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "c*rotBz - 4*pi*Jz = %g\n",
            (speed_of_light_normalized * kw * Byamplitude * fieldScale - 4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");

    double derivativJy = -electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * omega;
    fprintf(informationFile, "w*Jy amplitude = %g\n",
            derivativJy / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

    double derivativeVelocitiesY = electron_charge_normalized * ((Eyamplitude * fieldScale * ((1.0 / massProton) + (1.0 / massElectron))) + B0.norm() * fieldScale * ((VzamplitudeProton / massProton) + (VzamplitudeElectron / massElectron)) / speed_of_light_normalized);
    fprintf(informationFile, "dJy/dt amplitude = %g\n",
            electron_charge_normalized * concentration * derivativeVelocitiesY / (plasma_period * plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    double derivativJz = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * omega;
    fprintf(informationFile, "w*Jz amplitude = %g\n",
            derivativJz / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

    double derivativeVelocitiesZ = electron_charge_normalized * ((Ezamplitude * fieldScale * ((1.0 / massProton) + (1.0 / massElectron))) - B0.norm() * fieldScale * ((VyamplitudeProton / massProton) + (VyamplitudeElectron / massElectron)) / speed_of_light_normalized);
    fprintf(informationFile, "dJz/dt amplitude = %g\n",
            electron_charge_normalized * concentration * derivativeVelocitiesZ / (plasma_period * plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");

    double derivativVyp = -omega * VyamplitudeProton;
    fprintf(informationFile, "-w*Vyp amplitude = %g\n", derivativVyp * scaleFactor / sqr(plasma_period));

    double derivativeVelocityProtonY = electron_charge_normalized * (Eyamplitude * fieldScale + B0.norm() * fieldScale * VzamplitudeProton / speed_of_light_normalized) / massProton;
    fprintf(informationFile, "dVyp/dt amplitude = %g\n", derivativeVelocityProtonY * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    double derivativVzp = omega * VzamplitudeProton;
    fprintf(informationFile, "w*Vzp amplitude = %g\n", derivativVzp * scaleFactor / sqr(plasma_period));

    double derivativeVelocityProtonZ = electron_charge_normalized * (Ezamplitude * fieldScale - B0.norm() * fieldScale * VyamplitudeProton / speed_of_light_normalized) / massProton;
    fprintf(informationFile, "dVzp/dt amplitude = %g\n", derivativeVelocityProtonZ * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    double derivativVye = -omega * VyamplitudeElectron;
    fprintf(informationFile, "-w*Vye amplitude = %g\n", derivativVye * scaleFactor / sqr(plasma_period));

    double derivativeVelocityElectronY = -electron_charge_normalized * (Eyamplitude * fieldScale + B0.norm() * fieldScale * VzamplitudeElectron / speed_of_light_normalized) / massElectron;
    fprintf(informationFile, "dVye/dt amplitude = %g\n",
            derivativeVelocityElectronY * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    double derivativVze = omega * VzamplitudeElectron;
    fprintf(informationFile, "w*Vze amplitude = %g\n", derivativVze * scaleFactor / sqr(plasma_period));

    double derivativeVelocityElectronZ = -electron_charge_normalized * (Ezamplitude * fieldScale - B0.norm() * fieldScale * VyamplitudeElectron / speed_of_light_normalized) / massElectron;
    fprintf(informationFile, "dVze/dt amplitude = %g\n",
            derivativeVelocityElectronZ * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    fclose(informationFile);
}


void Simulation::initializeAlfvenWaveY(int wavesCount, double amplitudeRelation) {
    boundaryConditionType = PERIODIC;
    printf("initialization alfven wave\n");
    positronsPerBin = 0;
    alphaPerBin = 0;
    protonsPerBin = electronsPerBin;

    double concentration = density / (massProton + massElectron);
    types[1].particesDeltaX = types[0].particesDeltaX;
    types[1].particlesPerBin = types[0].particlesPerBin;
    types[0].concentration = concentration;
    types[1].concentration = concentration;
    for (int i = 2; i < typesNumber; ++i) {
        types[i].particlesPerBin = 0;
        types[i].concentration = 0;
        types[i].particesDeltaX = xsize;
    }
    createParticles();
    E0 = Vector3d(0, 0, 0);

    informationFile = fopen((outputDir + "information.dat").c_str(), "a");

    double alfvenV = B0.norm() * fieldScale / sqrt(4 * pi * density);
    if (alfvenV > speed_of_light_normalized) {
        printf("alfven velocity > c\n");
        fprintf(informationFile, "alfven velocity > c\n");
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "alfvenV/c = %15.10g > 1\n", alfvenV / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }
    fprintf(informationFile, "alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
    fprintf(informationFile, "alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);
    printf("alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
    printf("alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);

    //double kw = wavesCount * 2 * pi / xsize;
    double kw = wavesCount * 2 * pi / ysize;

    omegaPlasmaProton = sqrt(
        4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massProton);
    omegaPlasmaElectron = sqrt(
        4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
    omegaGyroProton = B0.norm() * fieldScale * electron_charge_normalized / (massProton * speed_of_light_normalized);
    omegaGyroElectron = B0.norm() * fieldScale * electron_charge_normalized / (massElectron * speed_of_light_normalized);

    if (omegaGyroProton < 5 * speed_of_light_normalized * kw) {
        printf("omegaGyroProton < 5*k*c\n");
        fprintf(informationFile, "omegaGyroProton < 5*k*c\n");
        //fclose(informationFile);
        //exit(0);
    }
    printf("omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));
    fprintf(informationFile, "omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));

    if (omegaPlasmaProton < 5 * omegaGyroProton) {
        printf("omegaPlasmaProton < 5*omegaGyroProton\n");
        fprintf(informationFile, "omegaPlasmaProton < 5*omegaGyroProton\n");
        //fclose(informationFile);
        //exit(0);
    }
    printf("omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);
    fprintf(informationFile, "omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);

    //w = q*kw*B/mP * 0.5*(sqrt(d)+-b)/a
    double b = speed_of_light_normalized * kw * (massProton - massElectron) / massProton;
    double discriminant = speed_of_light_normalized_sqr * kw * kw * sqr(
        massProton + massElectron) + 16 * pi * concentration * sqr(
        electron_charge_normalized) * (massProton + massElectron) / sqr(massProton);
    double a = (kw * kw * speed_of_light_normalized_sqr * massProton * massElectron + 4 * pi * concentration * sqr(
        electron_charge_normalized) * (massProton + massElectron)) / sqr(massProton);

    if (discriminant < 0) {
        printf("discriminant < 0\n");
        fprintf(informationFile, "discriminant < 0\n");
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "discriminant = %15.10g\n", discriminant);
        fclose(errorLogFile);
        exit(0);
    }

    double fakeOmega = (kw * electron_charge_normalized * B0.norm() * fieldScale / massProton) * (sqrt(
        discriminant) - b) / (2.0 * a);

    //a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0

    double a4 = sqr(speed_of_light_normalized_sqr * massProton * massElectron);
    a4 = a4 * sqr(sqr(sqr(fakeOmega)));
    double a3 = -2 * cube(speed_of_light_normalized_sqr) * sqr(kw * massElectron * massProton)
                - 8 * pi * concentration * sqr(
        speed_of_light_normalized_sqr * electron_charge_normalized) * massElectron * massProton * (massElectron + massProton)
                - sqr(B0.norm() * fieldScale * speed_of_light_normalized * electron_charge_normalized) * (sqr(
        massProton) + sqr(massElectron));
    a3 = a3 * cube(sqr(fakeOmega));
    double a2 = sqr(sqr(speed_of_light_normalized_sqr * kw) * massProton * massElectron)
                + 8 * pi * cube(speed_of_light_normalized_sqr) * concentration * sqr(
        kw * electron_charge_normalized) * massProton * massElectron * (massProton + massElectron)
                + 16 * sqr(pi * speed_of_light_normalized_sqr * concentration * sqr(
        electron_charge_normalized) * (massProton + massElectron))
                + 2 * sqr(
        B0.norm() * fieldScale * speed_of_light_normalized_sqr * kw * electron_charge_normalized) * (sqr(
        massProton) + sqr(massElectron))
                + 8 * pi * concentration * sqr(B0.norm() * fieldScale * speed_of_light_normalized * sqr(
        electron_charge_normalized)) * (massProton + massElectron)
                + sqr(sqr(B0.norm() * fieldScale * electron_charge_normalized));
    a2 = a2 * sqr(sqr(fakeOmega));
    double a1 = -sqr(
        B0.norm() * fieldScale * cube(speed_of_light_normalized) * kw * kw * electron_charge_normalized) * (sqr(
        massProton) + sqr(massElectron))
                - 8 * pi * concentration * sqr(B0.norm() * fieldScale * speed_of_light_normalized_sqr * kw * sqr(
        electron_charge_normalized)) * (massProton + massElectron)
                - 2 * sqr(speed_of_light_normalized * kw * sqr(B0.norm() * fieldScale * electron_charge_normalized));
    a1 = a1 * sqr(fakeOmega);
    double a0 = sqr(sqr(B0.norm() * fieldScale * speed_of_light_normalized * kw * electron_charge_normalized));

    a4 = a4 / a0;
    a3 = a3 / a0;
    a2 = a2 / a0;
    a1 = a1 / a0;
    a0 = 1.0;

    printf("a4 = %g\n", a4);
    fprintf(informationFile, "a4 = %g\n", a4);
    printf("a3 = %g\n", a3);
    fprintf(informationFile, "a3 = %g\n", a3);
    printf("a2 = %g\n", a2);
    fprintf(informationFile, "a2 = %g\n", a2);
    printf("a1 = %g\n", a1);
    fprintf(informationFile, "a1 = %g\n", a1);
    printf("a0 = %g\n", a0);
    fprintf(informationFile, "a0 = %g\n", a0);

    double fakeOmega1 = kw * alfvenV;
    printf("fakeOmega = %g\n", fakeOmega1 / plasma_period);
    fprintf(informationFile, "fakeOmega = %g\n", fakeOmega1 / plasma_period);
    double realOmega2 = solve4orderEquation(a4, a3, a2, a1, a0, 1.0);
    if (realOmega2 < 0) {
        printf("omega^2 < 0\n");
        fprintf(informationFile, "omega^2 < 0\n");
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "omega^2 = %15.10g > 1\n", realOmega2);
        fclose(errorLogFile);
        exit(0);
    }

    double error = (((a4 * realOmega2 + a3) * realOmega2 + a2) * realOmega2 + a1) * realOmega2 + a0;
    printf("error = %15.10g\n", error);
    fprintf(informationFile, "error = %15.10g\n", error);
    //double
    omega = sqrt(realOmega2) * fakeOmega;
    if (omega < 0) {
        omega = -omega;
    }

    if (omega > speed_of_light_normalized * kw / 5.0) {
        printf("omega > k*c/5\n");
        fprintf(informationFile, "omega > k*c/5\n");
        printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
        fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
        //fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "omega/kc = %15.10g > 0.2\n", omega / (kw * speed_of_light_normalized));
        fclose(errorLogFile);
        //exit(0);
    }
    printf("omega = %g\n", omega / plasma_period);
    fprintf(informationFile, "omega = %g\n", omega / plasma_period);

    printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
    fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));

    if (fabs(omega) > omegaGyroProton / 2) {
        printf("omega > omegaGyroProton/2\n");
        fprintf(informationFile, "omega > omegaGyroProton/2\n");
    }
    printf("omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
    fprintf(informationFile, "omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
    fclose(informationFile);

    checkFrequency(omega);

    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    //checkCollisionTime(omega);
    //checkMagneticReynolds(alfvenV);
    //checkDissipation(kw, alfvenV);

    double epsilonAmplitude = amplitudeRelation;

    double alfvenVReal = omega / kw;

    //double
    Bzamplitude = B0.norm() * epsilonAmplitude;

    double Omegae = omegaGyroElectron;
    double Omegae2 = Omegae * Omegae;
    double Omegap = omegaGyroProton;
    double Omegap2 = Omegap * Omegap;
    double omegae = omegaPlasmaElectron;
    double omegae2 = omegae * omegae;
    double omegap = omegaPlasmaProton;
    double omegap2 = omegap * omegap;

    double kc = kw * speed_of_light_normalized;
    double kc2 = kc * kc;

    double denominator = omega * omega - Omegae2 - (omegae2 * omega * omega / (omega * omega - kc2));

    //double
    VzamplitudeProton = -((1.0 / (4 * pi * concentration * electron_charge_normalized)) * (kc + ((omegae2 + omegap2 - omega * omega) / kc) + (omegae2 * Omegae2 / (kc * denominator))) / ((Omegae * omegae2 * omega / ((kc2 - omega * omega) * denominator)) + (Omegap / omega))) * Bzamplitude * fieldScale;
    //double
    VzamplitudeElectron = (((electron_charge_normalized * omega * Omegae) / (massElectron * kc)) * Bzamplitude * fieldScale + (omegae2 * omega * omega / (kc2 - omega * omega)) * VzamplitudeProton) / denominator;

    //double
    Byamplitude = (4 * pi * concentration * electron_charge_normalized / ((omega * omega / kc) - kc)) * (VzamplitudeElectron - VzamplitudeProton) / fieldScale;

    //double
    VyamplitudeProton = -(Omegap / omega) * VzamplitudeProton - (electron_charge_normalized / (massProton * kc)) * Bzamplitude * fieldScale;
    //double
    VyamplitudeElectron = (Omegae / omega) * VzamplitudeElectron + (electron_charge_normalized / (massElectron * kc)) * Bzamplitude * fieldScale;

    //double
    Eyamplitude = (omega / kc) * Bzamplitude;
    //double
    Ezamplitude = -(omega / kc) * Byamplitude;

    double xshift = 0;
    //double xshift = xsize/4;

    //Eyamplitude = 0.0;
    //VzamplitudeElectron = 0.0;
    //VzamplitudeProton = 0.0;
    //Byamplitude = 0.0;

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                /*Efield[i][j][k].x = 0;
                Efield[i][j][k].y = Eyamplitude * cos(kw * xgrid[i] - kw * xshift);
                Efield[i][j][k].z = Ezamplitude * sin(kw * xgrid[i] - kw * xshift);*/
                Efield[i][j][k].x = -Eyamplitude * cos(kw * ygrid[j] - kw * xshift);
                Efield[i][j][k].y = 0;
                Efield[i][j][k].z = Ezamplitude * sin(kw * ygrid[j] - kw * xshift);
                explicitEfield[i][j][k] = Efield[i][j][k];
                tempEfield[i][j][k] = Efield[i][j][k];
                newEfield[i][j][k] = Efield[i][j][k];
            }
        }
    }

    for (int k = 0; k < znumber; ++k) {
        for (int j = 0; j < ynumber; ++j) {
            Efield[xnumber][j][k] = Efield[0][j][k];
            tempEfield[xnumber][j][k] = Efield[0][j][k];
            newEfield[xnumber][j][k] = Efield[0][j][k];
            explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            Efield[i][j][znumber] = Efield[i][j][0];
            tempEfield[i][j][znumber] = Efield[i][j][0];
            newEfield[i][j][znumber] = Efield[i][j][0];
            explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
        }
    }

    for (int k = 0; k < znumber + 1; ++k) {
        for (int i = 0; i < xnumber + 1; ++i) {
            Efield[i][ynumber][k] = Efield[i][0][k];
            tempEfield[i][ynumber][k] = Efield[i][0][k];
            newEfield[i][ynumber][k] = Efield[i][0][k];
            explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
        }
    }

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                /*Bfield[i][j][k].x = B0.norm();
                Bfield[i][j][k].y = Byamplitude * sin(kw * middleXgrid[i] - kw * xshift);
                Bfield[i][j][k].z = Bzamplitude * cos(kw * middleXgrid[i] - kw * xshift);*/
                Bfield[i][j][k].x = -Byamplitude * sin(kw * middleYgrid[j] - kw * xshift);
                Bfield[i][j][k].y = B0.norm();
                Bfield[i][j][k].z = Bzamplitude * cos(kw * middleYgrid[j] - kw * xshift);
                newBfield[i][j][k] = Bfield[i][j][k];
            }
        }
    }

    if (fabs(VzamplitudeProton) > speed_of_light_normalized) {
        printf("VzamplitudeProton > speed_of_light_normalized\n");
        fprintf(informationFile, "VzamplitudeProton > speed_of_light_normalized\n");
        printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
        fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VzamplitudeProton/c = %15.10g > 1\n", VzamplitudeProton / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }
    printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
    fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);

    if (fabs(VzamplitudeElectron) > speed_of_light_normalized) {
        printf("VzamplitudeElectron > speed_of_light_normalized\n");
        fprintf(informationFile, "VzamplitudeElectron > speed_of_light_normalized\n");
        printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
        fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VzamplitudeElectron/c = %15.10g > 1\n", VzamplitudeElectron / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }
    printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
    fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);

    if (fabs(VyamplitudeProton) > speed_of_light_normalized) {
        printf("VyamplitudeProton > speed_of_light_normalized\n");
        fprintf(informationFile, "VyamplitudeProton > speed_of_light_normalized\n");
        printf("VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
        fprintf(informationFile, "VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VyamplitudeProton/c = %15.10g > 1\n", VyamplitudeProton / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }

    if (fabs(VyamplitudeElectron) > speed_of_light_normalized) {
        printf("VyamplitudeElectron > speed_of_light_normalized\n");
        fprintf(informationFile, "VyamplitudeElectron > speed_of_light_normalized\n");
        printf("VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
        fprintf(informationFile, "VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VyamplitudeElectron/c = %15.10g > 1\n", VyamplitudeElectron / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }

    //k > 0, w > 0, kx-wt, Bz > 0, Vz < 0, Vyp > 0, Vye < 0

    for (int pcount = 0; pcount < particles.size(); ++pcount) {
        Particle *particle = particles[pcount];
        Vector3d velocity = particle->velocity(speed_of_light_normalized);
        int yn = (particle->coordinates.y - ygrid[0]) / deltaY;
        double rightWeight = (particle->coordinates.y - ygrid[yn]) / deltaY;
        double leftWeight = (ygrid[yn + 1] - particle->coordinates.y) / deltaY;
        if (particle->type == PROTON) {
            //velocity = Vector3d(0, 0, 1) * (VzamplitudeProton) * (leftWeight*cos(kw * (ygrid[yn] - xshift)) + rightWeight*cos(kw*(ygrid[yn+1] - xshift))) + Vector3d(0, 1, 0) * VyamplitudeProton * (leftWeight*sin(kw * (ygrid[yn] - xshift)) + rightWeight*sin(kw*(ygrid[yn+1] - xshift)));
            velocity = (Vector3d(0, 0, 1) * (VzamplitudeProton) * cos(
                kw * particle->coordinates.y - kw * xshift) - Vector3d(1, 0, 0) * VyamplitudeProton * sin(
                kw * particle->coordinates.y - kw * xshift));

        }
        if (particle->type == ELECTRON) {

            //velocity = Vector3d(0, 0, 1) * (VzamplitudeElectron) * (leftWeight*cos(kw * (ygrid[yn] - xshift)) + rightWeight*cos(kw*(ygrid[yn+1] - xshift))) + Vector3d(0, 1, 0) * VyamplitudeElectron * (leftWeight*sin(kw * (ygrid[yn] - xshift)) + rightWeight*sin(kw*(ygrid[yn+1] - xshift)));
            velocity = (Vector3d(0, 0, 1) * (VzamplitudeElectron) * cos(
                kw * particle->coordinates.y - kw * xshift) - Vector3d(1, 0, 0) * VyamplitudeElectron * sin(
                kw * particle->coordinates.y - kw * xshift));
        }
        double beta = velocity.norm() / speed_of_light_normalized;
        particle->addVelocity(velocity, speed_of_light_normalized);
        particle->initialMomentum = particle->momentum;
    }

    updateDeltaT();

    printf("dt/Talfven = %g\n", deltaT * omega / (2 * pi));
    printf("dt = %g\n", deltaT * plasma_period);
    fprintf(informationFile, "dt/Talfven = %g\n", deltaT * omega / (2 * pi));
    fprintf(informationFile, "dt = %g\n", deltaT * plasma_period);

    double Vthermal = sqrt(2 * kBoltzman_normalized * temperature / massElectron);
    double thermalFlux = Vthermal * concentration * electron_charge_normalized / sqrt(1.0 * electronsPerBin);
    double alfvenFlux = (VyamplitudeProton - VyamplitudeElectron) * concentration * electron_charge_normalized;
    if (thermalFlux > alfvenFlux / 2) {
        printf("thermalFlux > alfvenFlux/2\n");
        fprintf(informationFile, "thermalFlux > alfvenFlux/2\n");
    }
    printf("alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
    fprintf(informationFile, "alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
    double minDeltaT = deltaX / Vthermal;
    if (minDeltaT > deltaT) {
        printf("deltaT < dx/Vthermal\n");
        fprintf(informationFile, "deltaT < dx/Vthermal\n");

        //printf("deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);
        //fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);

        //fclose(informationFile);
        //exit(0);
    }
    printf("deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
    fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
    fprintf(informationFile, "\n");

    fprintf(informationFile, "Bz amplitude = %g\n", Bzamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "By amplitude = %g\n", Byamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "Vz amplitude p = %g\n", VzamplitudeProton * scaleFactor / plasma_period);
    fprintf(informationFile, "Vz amplitude e = %g\n", VzamplitudeElectron * scaleFactor / plasma_period);
    fprintf(informationFile, "Vy amplitude p = %g\n", VyamplitudeProton * scaleFactor / plasma_period);
    fprintf(informationFile, "Vy amplitude e = %g\n", VyamplitudeElectron * scaleFactor / plasma_period);
    fprintf(informationFile, "Ey amplitude = %g\n", Eyamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "Ez amplitude = %g\n", Ezamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "By/Ez = %g\n", Byamplitude / Ezamplitude);
    fprintf(informationFile, "Bz/Ey = %g\n", Bzamplitude / Eyamplitude);
    fprintf(informationFile, "4*pi*Jy amplitude = %g\n",
            4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "c*rotBy amplitude = %g\n",
            speed_of_light_normalized * kw * Bzamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "4*pi*Jz amplitude = %g\n",
            4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "c*rotBz amplitude = %g\n",
            speed_of_light_normalized * kw * Byamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative By amplitude = %g\n",
            -omega * Byamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "-c*rotEy = %g\n",
            speed_of_light_normalized * kw * Ezamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative Bz amplitude = %g\n",
            omega * Bzamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "-c*rotEz = %g\n",
            speed_of_light_normalized * kw * Eyamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative Ey amplitude = %g\n",
            omega * Eyamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "c*rotBy - 4*pi*Jy = %g\n",
            (speed_of_light_normalized * kw * Bzamplitude * fieldScale - 4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative Ez amplitude = %g\n",
            -omega * Ezamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "c*rotBz - 4*pi*Jz = %g\n",
            (speed_of_light_normalized * kw * Byamplitude * fieldScale - 4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");

    double derivativJy = -electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * omega;
    fprintf(informationFile, "w*Jy amplitude = %g\n",
            derivativJy / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

    double derivativeVelocitiesY = electron_charge_normalized * ((Eyamplitude * fieldScale * ((1.0 / massProton) + (1.0 / massElectron))) + B0.norm() * fieldScale * ((VzamplitudeProton / massProton) + (VzamplitudeElectron / massElectron)) / speed_of_light_normalized);
    fprintf(informationFile, "dJy/dt amplitude = %g\n",
            electron_charge_normalized * concentration * derivativeVelocitiesY / (plasma_period * plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    double derivativJz = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * omega;
    fprintf(informationFile, "w*Jz amplitude = %g\n",
            derivativJz / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

    double derivativeVelocitiesZ = electron_charge_normalized * ((Ezamplitude * fieldScale * ((1.0 / massProton) + (1.0 / massElectron))) - B0.norm() * fieldScale * ((VyamplitudeProton / massProton) + (VyamplitudeElectron / massElectron)) / speed_of_light_normalized);
    fprintf(informationFile, "dJz/dt amplitude = %g\n",
            electron_charge_normalized * concentration * derivativeVelocitiesZ / (plasma_period * plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");

    double derivativVyp = -omega * VyamplitudeProton;
    fprintf(informationFile, "-w*Vyp amplitude = %g\n", derivativVyp * scaleFactor / sqr(plasma_period));

    double derivativeVelocityProtonY = electron_charge_normalized * (Eyamplitude * fieldScale + B0.norm() * fieldScale * VzamplitudeProton / speed_of_light_normalized) / massProton;
    fprintf(informationFile, "dVyp/dt amplitude = %g\n", derivativeVelocityProtonY * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    double derivativVzp = omega * VzamplitudeProton;
    fprintf(informationFile, "w*Vzp amplitude = %g\n", derivativVzp * scaleFactor / sqr(plasma_period));

    double derivativeVelocityProtonZ = electron_charge_normalized * (Ezamplitude * fieldScale - B0.norm() * fieldScale * VyamplitudeProton / speed_of_light_normalized) / massProton;
    fprintf(informationFile, "dVzp/dt amplitude = %g\n", derivativeVelocityProtonZ * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    double derivativVye = -omega * VyamplitudeElectron;
    fprintf(informationFile, "-w*Vye amplitude = %g\n", derivativVye * scaleFactor / sqr(plasma_period));

    double derivativeVelocityElectronY = -electron_charge_normalized * (Eyamplitude * fieldScale + B0.norm() * fieldScale * VzamplitudeElectron / speed_of_light_normalized) / massElectron;
    fprintf(informationFile, "dVye/dt amplitude = %g\n",
            derivativeVelocityElectronY * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    double derivativVze = omega * VzamplitudeElectron;
    fprintf(informationFile, "w*Vze amplitude = %g\n", derivativVze * scaleFactor / sqr(plasma_period));

    double derivativeVelocityElectronZ = -electron_charge_normalized * (Ezamplitude * fieldScale - B0.norm() * fieldScale * VyamplitudeElectron / speed_of_light_normalized) / massElectron;
    fprintf(informationFile, "dVze/dt amplitude = %g\n",
            derivativeVelocityElectronZ * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    fclose(informationFile);
}

void Simulation::initializeRotatedAlfvenWave(int wavesCount, double amplitudeRelation) {
    boundaryConditionType = PERIODIC;
    printf("initialization alfven wave\n");
    positronsPerBin = 0;
    alphaPerBin = 0;
    protonsPerBin = electronsPerBin;

    double concentration = density / (massProton + massElectron);
    types[1].particesDeltaX = types[0].particesDeltaX;
    types[1].particlesPerBin = types[0].particlesPerBin;
    types[0].concentration = concentration;
    types[1].concentration = concentration;
    for (int i = 2; i < typesNumber; ++i) {
        types[i].particlesPerBin = 0;
        types[i].concentration = 0;
        types[i].particesDeltaX = xsize;
    }
    createParticles();
    E0 = Vector3d(0, 0, 0);

    informationFile = fopen((outputDir + "information.dat").c_str(), "a");

    double alfvenV = B0.norm() * fieldScale / sqrt(4 * pi * density);
    if (alfvenV > speed_of_light_normalized) {
        printf("alfven velocity > c\n");
        fprintf(informationFile, "alfven velocity > c\n");
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "alfvenV/c = %15.10g > 1\n", alfvenV / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }
    fprintf(informationFile, "alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
    fprintf(informationFile, "alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);
    printf("alfven V = %lf\n", alfvenV * scaleFactor / plasma_period);
    printf("alfven V/c = %lf\n", alfvenV / speed_of_light_normalized);

    double kx = wavesCount * 2 * pi / xsize;
    double ky = wavesCount * 2 * pi / ysize;
    double kz = wavesCount * 2 * pi / zsize;
    kz = 0;

    double kw = sqrt(kx * kx + ky * ky + kz * kz);

    double weight = concentration * volumeB(0, 0, 0) / electronsPerBin;

    omegaPlasmaProton = sqrt(
        4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massProton);
    omegaPlasmaElectron = sqrt(
        4 * pi * concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
    omegaGyroProton = B0.norm() * fieldScale * electron_charge_normalized / (massProton * speed_of_light_normalized);
    omegaGyroElectron = B0.norm() * fieldScale * electron_charge_normalized / (massElectron * speed_of_light_normalized);

    if (omegaGyroProton < 5 * speed_of_light_normalized * kw) {
        printf("omegaGyroProton < 5*k*c\n");
        fprintf(informationFile, "omegaGyroProton < 5*k*c\n");
        //fclose(informationFile);
        //exit(0);
    }
    printf("omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));
    fprintf(informationFile, "omegaGyroProton/kc = %g\n", omegaGyroProton / (kw * speed_of_light_normalized));

    if (omegaPlasmaProton < 5 * omegaGyroProton) {
        printf("omegaPlasmaProton < 5*omegaGyroProton\n");
        fprintf(informationFile, "omegaPlasmaProton < 5*omegaGyroProton\n");
        //fclose(informationFile);
        //exit(0);
    }
    printf("omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);
    fprintf(informationFile, "omegaPlasmaProton/omegaGyroProton = %g\n", omegaPlasmaProton / omegaGyroProton);

    //w = q*kw*B/mP * 0.5*(sqrt(d)+-b)/a
    double b = speed_of_light_normalized * kw * (massProton - massElectron) / massProton;
    double discriminant = speed_of_light_normalized_sqr * kw * kw * sqr(
        massProton + massElectron) + 16 * pi * concentration * sqr(
        electron_charge_normalized) * (massProton + massElectron) / sqr(massProton);
    double a = (kw * kw * speed_of_light_normalized_sqr * massProton * massElectron + 4 * pi * concentration * sqr(
        electron_charge_normalized) * (massProton + massElectron)) / sqr(massProton);

    if (discriminant < 0) {
        printf("discriminant < 0\n");
        fprintf(informationFile, "discriminant < 0\n");
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "discriminant = %15.10g\n", discriminant);
        fclose(errorLogFile);
        exit(0);
    }

    double fakeOmega = (kw * electron_charge_normalized * B0.norm() * fieldScale / massProton) * (sqrt(
        discriminant) - b) / (2.0 * a);

    //a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0

    double a4 = sqr(speed_of_light_normalized_sqr * massProton * massElectron);
    a4 = a4 * sqr(sqr(sqr(fakeOmega)));
    double a3 = -2 * cube(speed_of_light_normalized_sqr) * sqr(kw * massElectron * massProton)
                - 8 * pi * concentration * sqr(
        speed_of_light_normalized_sqr * electron_charge_normalized) * massElectron * massProton * (massElectron + massProton)
                - sqr(B0.norm() * fieldScale * speed_of_light_normalized * electron_charge_normalized) * (sqr(
        massProton) + sqr(massElectron));
    a3 = a3 * cube(sqr(fakeOmega));
    double a2 = sqr(sqr(speed_of_light_normalized_sqr * kw) * massProton * massElectron)
                + 8 * pi * cube(speed_of_light_normalized_sqr) * concentration * sqr(
        kw * electron_charge_normalized) * massProton * massElectron * (massProton + massElectron)
                + 16 * sqr(pi * speed_of_light_normalized_sqr * concentration * sqr(
        electron_charge_normalized) * (massProton + massElectron))
                + 2 * sqr(
        B0.norm() * fieldScale * speed_of_light_normalized_sqr * kw * electron_charge_normalized) * (sqr(
        massProton) + sqr(massElectron))
                + 8 * pi * concentration * sqr(B0.norm() * fieldScale * speed_of_light_normalized * sqr(
        electron_charge_normalized)) * (massProton + massElectron)
                + sqr(sqr(B0.norm() * fieldScale * electron_charge_normalized));
    a2 = a2 * sqr(sqr(fakeOmega));
    double a1 = -sqr(
        B0.norm() * fieldScale * cube(speed_of_light_normalized) * kw * kw * electron_charge_normalized) * (sqr(
        massProton) + sqr(massElectron))
                - 8 * pi * concentration * sqr(B0.norm() * fieldScale * speed_of_light_normalized_sqr * kw * sqr(
        electron_charge_normalized)) * (massProton + massElectron)
                - 2 * sqr(speed_of_light_normalized * kw * sqr(B0.norm() * fieldScale * electron_charge_normalized));
    a1 = a1 * sqr(fakeOmega);
    double a0 = sqr(sqr(B0.norm() * fieldScale * speed_of_light_normalized * kw * electron_charge_normalized));

    a4 = a4 / a0;
    a3 = a3 / a0;
    a2 = a2 / a0;
    a1 = a1 / a0;
    a0 = 1.0;

    printf("a4 = %g\n", a4);
    fprintf(informationFile, "a4 = %g\n", a4);
    printf("a3 = %g\n", a3);
    fprintf(informationFile, "a3 = %g\n", a3);
    printf("a2 = %g\n", a2);
    fprintf(informationFile, "a2 = %g\n", a2);
    printf("a1 = %g\n", a1);
    fprintf(informationFile, "a1 = %g\n", a1);
    printf("a0 = %g\n", a0);
    fprintf(informationFile, "a0 = %g\n", a0);

    double fakeOmega1 = kw * alfvenV;
    printf("fakeOmega = %g\n", fakeOmega1 / plasma_period);
    fprintf(informationFile, "fakeOmega = %g\n", fakeOmega1 / plasma_period);
    double realOmega2 = solve4orderEquation(a4, a3, a2, a1, a0, 1.0);
    if (realOmega2 < 0) {
        printf("omega^2 < 0\n");
        fprintf(informationFile, "omega^2 < 0\n");
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "omega^2 = %15.10g > 1\n", realOmega2);
        fclose(errorLogFile);
        exit(0);
    }

    double error = (((a4 * realOmega2 + a3) * realOmega2 + a2) * realOmega2 + a1) * realOmega2 + a0;
    printf("error = %15.10g\n", error);
    fprintf(informationFile, "error = %15.10g\n", error);
    //double
    omega = sqrt(realOmega2) * fakeOmega;
    if (omega < 0) {
        omega = -omega;
    }

    if (omega > speed_of_light_normalized * kw / 5.0) {
        printf("omega > k*c/5\n");
        fprintf(informationFile, "omega > k*c/5\n");
        printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
        fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
        //fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "omega/kc = %15.10g > 0.2\n", omega / (kw * speed_of_light_normalized));
        fclose(errorLogFile);
        //exit(0);
    }
    printf("omega = %g\n", omega / plasma_period);
    fprintf(informationFile, "omega = %g\n", omega / plasma_period);

    printf("omega/kc = %g\n", omega / (kw * speed_of_light_normalized));
    fprintf(informationFile, "omega/kc = %g\n", omega / (kw * speed_of_light_normalized));

    if (fabs(omega) > omegaGyroProton / 2) {
        printf("omega > omegaGyroProton/2\n");
        fprintf(informationFile, "omega > omegaGyroProton/2\n");
    }
    printf("omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
    fprintf(informationFile, "omega/omegaGyroProton = %g\n", omega / omegaGyroProton);
    fclose(informationFile);

    checkFrequency(omega);

    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    //checkCollisionTime(omega);
    //checkMagneticReynolds(alfvenV);
    //checkDissipation(kw, alfvenV);

    double epsilonAmplitude = amplitudeRelation;

    double alfvenVReal = omega / kw;

    fprintf(informationFile, "alfven V real = %15.10g\n", alfvenVReal * scaleFactor / plasma_period);
    fprintf(informationFile, "alfven V real x = %15.10g\n", alfvenVReal * (kx / kw) * scaleFactor / plasma_period);

    //double
    Bzamplitude = B0.norm() * epsilonAmplitude;

    double Omegae = omegaGyroElectron;
    double Omegae2 = Omegae * Omegae;
    double Omegap = omegaGyroProton;
    double Omegap2 = Omegap * Omegap;
    double omegae = omegaPlasmaElectron;
    double omegae2 = omegae * omegae;
    double omegap = omegaPlasmaProton;
    double omegap2 = omegap * omegap;

    double kc = kw * speed_of_light_normalized;
    double kc2 = kc * kc;

    double denominator = omega * omega - Omegae2 - (omegae2 * omega * omega / (omega * omega - kc2));

    //double
    VzamplitudeProton = -((1.0 / (4 * pi * concentration * electron_charge_normalized)) * (kc + ((omegae2 + omegap2 - omega * omega) / kc) + (omegae2 * Omegae2 / (kc * denominator))) / ((Omegae * omegae2 * omega / ((kc2 - omega * omega) * denominator)) + (Omegap / omega))) * Bzamplitude * fieldScale;
    //double
    VzamplitudeElectron = (((electron_charge_normalized * omega * Omegae) / (massElectron * kc)) * Bzamplitude * fieldScale + (omegae2 * omega * omega / (kc2 - omega * omega)) * VzamplitudeProton) / denominator;

    //double
    Byamplitude = (4 * pi * concentration * electron_charge_normalized / ((omega * omega / kc) - kc)) * (VzamplitudeElectron - VzamplitudeProton) / fieldScale;

    //double
    VyamplitudeProton = -(Omegap / omega) * VzamplitudeProton - (electron_charge_normalized / (massProton * kc)) * Bzamplitude * fieldScale;
    //double
    VyamplitudeElectron = (Omegae / omega) * VzamplitudeElectron + (electron_charge_normalized / (massElectron * kc)) * Bzamplitude * fieldScale;

    //double
    Eyamplitude = (omega / kc) * Bzamplitude;
    //double
    Ezamplitude = -(omega / kc) * Byamplitude;

    double xshift = 0.0;

    //Eyamplitude = 0.0;
    //VzamplitudeElectron = 0.0;
    //VzamplitudeProton = 0.0;
    //Byamplitude = 0.0;


    double kxy = sqrt(kx * kx + ky * ky);
    double rotatedZortNorm = sqrt(kx * kx + ky * ky + sqr(kx * kx + ky * ky) / (kz * kz));
    double matrixzz = (kx * kx + ky * ky) / (kz * rotatedZortNorm);

    Matrix3d rotationMatrix = Matrix3d(kx / kw, -ky / kxy, -kx / rotatedZortNorm,
                                       ky / kw, kx / kxy, -ky / rotatedZortNorm,
                                       kz / kw, 0, matrixzz);

    if (kz == 0) {
        rotationMatrix = Matrix3d(kx / kw, -ky / kw, 0,
                                  ky / kw, kx / kw, 0,
                                  0, 0, 1);
    }
    //Matrix3d inverse = *(rotationMatrix.Inverse());

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k].x = 0;
                Efield[i][j][k].y = Eyamplitude * cos(kx * xgrid[i] + ky * ygrid[j] + kz * zgrid[k]);
                Efield[i][j][k].z = Ezamplitude * sin(kx * xgrid[i] + ky * ygrid[j] + kz * zgrid[k]);
                Efield[i][j][k] = rotationMatrix * Efield[i][j][k];
                //Efield[i][j][k] = inverse * Efield[i][j][k];
                explicitEfield[i][j][k] = Efield[i][j][k];
                tempEfield[i][j][k] = Efield[i][j][k];
                newEfield[i][j][k] = Efield[i][j][k];
            }
        }
    }

    for (int k = 0; k < znumber; ++k) {
        for (int j = 0; j < ynumber; ++j) {
            Efield[xnumber][j][k] = Efield[0][j][k];
            tempEfield[xnumber][j][k] = Efield[0][j][k];
            newEfield[xnumber][j][k] = Efield[0][j][k];
            explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            Efield[i][j][znumber] = Efield[i][j][0];
            tempEfield[i][j][znumber] = Efield[i][j][0];
            newEfield[i][j][znumber] = Efield[i][j][0];
            explicitEfield[i][j][znumber] = explicitEfield[i][j][0];
        }
    }

    for (int k = 0; k < znumber + 1; ++k) {
        for (int i = 0; i < xnumber + 1; ++i) {
            Efield[i][ynumber][k] = Efield[i][0][k];
            tempEfield[i][ynumber][k] = Efield[i][0][k];
            newEfield[i][ynumber][k] = Efield[i][0][k];
            explicitEfield[i][ynumber][k] = explicitEfield[i][0][k];
        }
    }

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Bfield[i][j][k].x = B0.norm();
                Bfield[i][j][k].y = Byamplitude * sin(kx * middleXgrid[i] + ky * middleYgrid[j] + kz * middleZgrid[k]);
                Bfield[i][j][k].z = Bzamplitude * cos(kx * middleXgrid[i] + ky * middleYgrid[j] + kz * middleZgrid[k]);
                Bfield[i][j][k] = rotationMatrix * Bfield[i][j][k];
                //Bfield[i][j][k] = inverse * Bfield[i][j][k];
                newBfield[i][j][k] = Bfield[i][j][k];
            }
        }
    }

    if (fabs(VzamplitudeProton) > speed_of_light_normalized) {
        printf("VzamplitudeProton > speed_of_light_normalized\n");
        fprintf(informationFile, "VzamplitudeProton > speed_of_light_normalized\n");
        printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
        fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VzamplitudeProton/c = %15.10g > 1\n", VzamplitudeProton / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }
    printf("VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);
    fprintf(informationFile, "VzamplitudeProton/c = %g\n", VzamplitudeProton / speed_of_light_normalized);

    if (fabs(VzamplitudeElectron) > speed_of_light_normalized) {
        printf("VzamplitudeElectron > speed_of_light_normalized\n");
        fprintf(informationFile, "VzamplitudeElectron > speed_of_light_normalized\n");
        printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
        fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VzamplitudeElectron/c = %15.10g > 1\n", VzamplitudeElectron / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }
    printf("VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);
    fprintf(informationFile, "VzamplitudeElectron/c = %g\n", VzamplitudeElectron / speed_of_light_normalized);

    if (fabs(VyamplitudeProton) > speed_of_light_normalized) {
        printf("VyamplitudeProton > speed_of_light_normalized\n");
        fprintf(informationFile, "VyamplitudeProton > speed_of_light_normalized\n");
        printf("VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
        fprintf(informationFile, "VyamplitudeProton/c = %g\n", VyamplitudeProton / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VyamplitudeProton/c = %15.10g > 1\n", VyamplitudeProton / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }

    if (fabs(VyamplitudeElectron) > speed_of_light_normalized) {
        printf("VyamplitudeElectron > speed_of_light_normalized\n");
        fprintf(informationFile, "VyamplitudeElectron > speed_of_light_normalized\n");
        printf("VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
        fprintf(informationFile, "VyamplitudeElectron/c = %g\n", VyamplitudeElectron / speed_of_light_normalized);
        fclose(informationFile);
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "VyamplitudeElectron/c = %15.10g > 1\n", VyamplitudeElectron / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }

    //k > 0, w > 0, kx-wt, Bz > 0, Vz < 0, Vyp > 0, Vye < 0

    for (int pcount = 0; pcount < particles.size(); ++pcount) {
        Particle *particle = particles[pcount];
        Vector3d velocity = particle->velocity(speed_of_light_normalized);
        int xn = (particle->coordinates.x - xgrid[0]) / deltaX;
        int yn = (particle->coordinates.y - ygrid[0]) / deltaY;
        double rightWeightX = (particle->coordinates.x - xgrid[xn]) / deltaX;
        double leftWeightX = (xgrid[xn + 1] - particle->coordinates.x) / deltaX;
        double rightWeightY = (particle->coordinates.y - ygrid[yn]) / deltaY;
        double leftWeightY = (ygrid[yn + 1] - particle->coordinates.y) / deltaY;
        if (particle->type == PROTON) {
            velocity = (Vector3d(0, 0, 1) * (VzamplitudeProton) * cos(kx * particle->coordinates.x + ky * particle->coordinates.y + kz * particle->coordinates.z) + Vector3d(0, 1, 0) * VyamplitudeProton * sin(kx * particle->coordinates.x + ky * particle->coordinates.y + kz * particle->coordinates.z));
            /*velocity = Vector3d(0, 0, 1) * (VzamplitudeProton) * (leftWeightX * (leftWeightY * cos(
                kx * (xgrid[xn] - xshift) + ky * ygrid[yn]) + rightWeightY * cos(
                kx * (xgrid[xn] - xshift) + ky * ygrid[yn + 1])) + rightWeightX * (leftWeightY * cos(
                kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn]) + rightWeightY * cos(
                kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn + 1])))
                       + Vector3d(0, 1, 0) * VyamplitudeProton * (leftWeightX * (leftWeightY * sin(
                kx * (xgrid[xn] - xshift) + ky * ygrid[yn]) + rightWeightY * sin(
                kx * (xgrid[xn] - xshift) + ky * ygrid[yn + 1])) + rightWeightX * (leftWeightY * sin(
                kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn]) + rightWeightY * sin(
                kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn + 1])));*/
        }
        if (particle->type == ELECTRON) {
            velocity = (Vector3d(0, 0, 1) * (VzamplitudeElectron) * cos(kx * particle->coordinates.x + ky * particle->coordinates.y + kz * particle->coordinates.z) + Vector3d(0, 1, 0) * VyamplitudeElectron * sin(kx * particle->coordinates.x + ky * particle->coordinates.y + kz * particle->coordinates.z));
            /*velocity = Vector3d(0, 0, 1) * (VzamplitudeElectron) * (leftWeightX * (leftWeightY * cos(
                kx * (xgrid[xn] - xshift) + ky * ygrid[yn]) + rightWeightY * cos(
                kx * (xgrid[xn] - xshift) + ky * ygrid[yn + 1])) + rightWeightX * (leftWeightY * cos(
                kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn]) + rightWeightY * cos(
                kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn + 1])))
                       + Vector3d(0, 1, 0) * VyamplitudeElectron * (leftWeightX * (leftWeightY * sin(
                kx * (xgrid[xn] - xshift) + ky * ygrid[yn]) + rightWeightY * sin(
                kx * (xgrid[xn] - xshift) + ky * ygrid[yn + 1])) + rightWeightX * (leftWeightY * sin(
                kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn]) + rightWeightY * sin(
                kx * (xgrid[xn + 1] - xshift) + ky * ygrid[yn + 1])));*/
        }
        velocity = rotationMatrix * velocity;
        //velocity = inverse * velocity;
        double beta = velocity.norm() / speed_of_light_normalized;
        particle->addVelocity(velocity, speed_of_light_normalized);
    }

    updateDeltaT();

    printf("dt/Talfven = %g\n", deltaT * omega / (2 * pi));
    printf("dt = %g\n", deltaT * plasma_period);
    fprintf(informationFile, "dt/Talfven = %g\n", deltaT * omega / (2 * pi));
    fprintf(informationFile, "dt = %g\n", deltaT * plasma_period);

    double Vthermal = sqrt(2 * kBoltzman_normalized * temperature / massElectron);
    double thermalFlux = Vthermal * concentration * electron_charge_normalized / sqrt(1.0 * electronsPerBin);
    double alfvenFlux = (VyamplitudeProton - VyamplitudeElectron) * concentration * electron_charge_normalized;
    if (thermalFlux > alfvenFlux / 2) {
        printf("thermalFlux > alfvenFlux/2\n");
        fprintf(informationFile, "thermalFlux > alfvenFlux/2\n");
    }
    printf("alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
    fprintf(informationFile, "alfvenFlux/thermalFlux = %g\n", alfvenFlux / thermalFlux);
    double minDeltaT = deltaX / Vthermal;
    if (minDeltaT > deltaT) {
        printf("deltaT < dx/Vthermal\n");
        fprintf(informationFile, "deltaT < dx/Vthermal\n");

        //printf("deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);
        //fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT/minDeltaT);

        //fclose(informationFile);
        //exit(0);
    }
    printf("deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
    fprintf(informationFile, "deltaT/minDeltaT =  %g\n", deltaT / minDeltaT);
    fprintf(informationFile, "\n");

    fprintf(informationFile, "Bz amplitude = %g\n", Bzamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "By amplitude = %g\n", Byamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "Vz amplitude p = %g\n", VzamplitudeProton * scaleFactor / plasma_period);
    fprintf(informationFile, "Vz amplitude e = %g\n", VzamplitudeElectron * scaleFactor / plasma_period);
    fprintf(informationFile, "Vy amplitude p = %g\n", VyamplitudeProton * scaleFactor / plasma_period);
    fprintf(informationFile, "Vy amplitude e = %g\n", VyamplitudeElectron * scaleFactor / plasma_period);
    fprintf(informationFile, "Ey amplitude = %g\n", Eyamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "Ez amplitude = %g\n", Ezamplitude * fieldScale / (plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "By/Ez = %g\n", Byamplitude / Ezamplitude);
    fprintf(informationFile, "Bz/Ey = %g\n", Bzamplitude / Eyamplitude);
    fprintf(informationFile, "4*pi*Jy amplitude = %g\n",
            4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "c*rotBy amplitude = %g\n",
            speed_of_light_normalized * kw * Bzamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "4*pi*Jz amplitude = %g\n",
            4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "c*rotBz amplitude = %g\n",
            speed_of_light_normalized * kw * Byamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative By amplitude = %g\n",
            -omega * Byamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "-c*rotEy = %g\n",
            speed_of_light_normalized * kw * Ezamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative Bz amplitude = %g\n",
            omega * Bzamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "-c*rotEz = %g\n",
            speed_of_light_normalized * kw * Eyamplitude * fieldScale / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative Ey amplitude = %g\n",
            omega * Eyamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "c*rotBy - 4*pi*Jy = %g\n",
            (speed_of_light_normalized * kw * Bzamplitude * fieldScale - 4 * pi * concentration * electron_charge_normalized * (VyamplitudeProton - VyamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    fprintf(informationFile, "derivative Ez amplitude = %g\n",
            -omega * Ezamplitude * fieldScale / (plasma_period * plasma_period * sqrt(scaleFactor)));
    fprintf(informationFile, "c*rotBz - 4*pi*Jz = %g\n",
            (speed_of_light_normalized * kw * Byamplitude * fieldScale - 4 * pi * concentration * electron_charge_normalized * (VzamplitudeProton - VzamplitudeElectron)) / (plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");

    double derivativJy = -electron_charge_normalized * concentration * (VyamplitudeProton - VyamplitudeElectron) * omega;
    fprintf(informationFile, "w*Jy amplitude = %g\n",
            derivativJy / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

    double derivativeVelocitiesY = electron_charge_normalized * ((Eyamplitude * fieldScale * ((1.0 / massProton) + (1.0 / massElectron))) + B0.norm() * fieldScale * ((VzamplitudeProton / massProton) + (VzamplitudeElectron / massElectron)) / speed_of_light_normalized);
    fprintf(informationFile, "dJy/dt amplitude = %g\n",
            electron_charge_normalized * concentration * derivativeVelocitiesY / (plasma_period * plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");
    double derivativJz = electron_charge_normalized * concentration * (VzamplitudeProton - VzamplitudeElectron) * omega;
    fprintf(informationFile, "w*Jz amplitude = %g\n",
            derivativJz / (plasma_period * plasma_period * plasma_period * sqrt(scaleFactor)));

    double derivativeVelocitiesZ = electron_charge_normalized * ((Ezamplitude * fieldScale * ((1.0 / massProton) + (1.0 / massElectron))) - B0.norm() * fieldScale * ((VyamplitudeProton / massProton) + (VyamplitudeElectron / massElectron)) / speed_of_light_normalized);
    fprintf(informationFile, "dJz/dt amplitude = %g\n",
            electron_charge_normalized * concentration * derivativeVelocitiesZ / (plasma_period * plasma_period * plasma_period * sqrt(
                scaleFactor)));
    fprintf(informationFile, "\n");

    double derivativVyp = -omega * VyamplitudeProton;
    fprintf(informationFile, "-w*Vyp amplitude = %g\n", derivativVyp * scaleFactor / sqr(plasma_period));

    double derivativeVelocityProtonY = electron_charge_normalized * (Eyamplitude * fieldScale + B0.norm() * fieldScale * VzamplitudeProton / speed_of_light_normalized) / massProton;
    fprintf(informationFile, "dVyp/dt amplitude = %g\n", derivativeVelocityProtonY * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    double derivativVzp = omega * VzamplitudeProton;
    fprintf(informationFile, "w*Vzp amplitude = %g\n", derivativVzp * scaleFactor / sqr(plasma_period));

    double derivativeVelocityProtonZ = electron_charge_normalized * (Ezamplitude * fieldScale - B0.norm() * fieldScale * VyamplitudeProton / speed_of_light_normalized) / massProton;
    fprintf(informationFile, "dVzp/dt amplitude = %g\n", derivativeVelocityProtonZ * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    double derivativVye = -omega * VyamplitudeElectron;
    fprintf(informationFile, "-w*Vye amplitude = %g\n", derivativVye * scaleFactor / sqr(plasma_period));

    double derivativeVelocityElectronY = -electron_charge_normalized * (Eyamplitude * fieldScale + B0.norm() * fieldScale * VzamplitudeElectron / speed_of_light_normalized) / massElectron;
    fprintf(informationFile, "dVye/dt amplitude = %g\n",
            derivativeVelocityElectronY * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    double derivativVze = omega * VzamplitudeElectron;
    fprintf(informationFile, "w*Vze amplitude = %g\n", derivativVze * scaleFactor / sqr(plasma_period));

    double derivativeVelocityElectronZ = -electron_charge_normalized * (Ezamplitude * fieldScale - B0.norm() * fieldScale * VyamplitudeElectron / speed_of_light_normalized) / massElectron;
    fprintf(informationFile, "dVze/dt amplitude = %g\n",
            derivativeVelocityElectronZ * scaleFactor / sqr(plasma_period));
    fprintf(informationFile, "\n");

    fclose(informationFile);
}

void Simulation::initializeLangmuirWave() {
    boundaryConditionType = PERIODIC;
    positronsPerBin = 0;
    alphaPerBin = 0;
    protonsPerBin = electronsPerBin;
    double concentration = density / (massProton + massElectron);
    types[1].particesDeltaX = types[0].particesDeltaX;
    types[1].particlesPerBin = types[0].particlesPerBin;
    types[0].concentration = concentration;
    types[1].concentration = concentration;
    for (int i = 2; i < typesNumber; ++i) {
        types[i].particlesPerBin = 0;
        types[i].concentration = 0;
        types[i].particesDeltaX = xsize;
    }
    double epsilon = 0.1;
    double kw = 2 * 2 * pi / xsize;
    double omega = 2 * pi;
    double langmuirV = omega / kw;

    checkDebyeParameter();
    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    fprintf(informationFile, "lengmuir V = %lf\n", langmuirV * scaleFactor / plasma_period);
    fclose(informationFile);
    if (langmuirV > speed_of_light_normalized) {
        printf("langmuirV > c\n");
        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
        fprintf(errorLogFile, "langmuireV/c = %15.10g > 1\n", langmuirV / speed_of_light_normalized);
        fclose(errorLogFile);
        exit(0);
    }

    printf("creating particles\n");
    int nproton = 0;
    int nelectron = 0;
    double weight = (concentration / electronsPerBin) * volumeB(0, 0, 0);
    /*for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                double x;
                for (int l = 0; l < particlesPerBin; ++l) {
                    ParticleTypes type;
                    type = PROTON;
                    Particle* particle = createParticle(particlesNumber, i, j, k, weight, type, temperature);
                    nproton++;
                    particles.push_back(particle);
                    particlesNumber++;
                    if (particlesNumber % 1000 == 0) {
                        printf("create particle number %d\n", particlesNumber);
                    }
                }
                for (int l = 0; l < particlesPerBin * (1 + epsilon * cos(kw * middleXgrid[i])); ++l) {
                    ParticleTypes type;
                    type = ELECTRON;
                    Particle* particle = createParticle(particlesNumber, i, j, k, weight, type, temperature);
                    nelectron++;
                    particles.push_back(particle);
                    particlesNumber++;
                    if (particlesNumber % 1000 == 0) {
                        printf("create particle number %d\n", particlesNumber);
                    }
                }
            }
        }
    }
    if (nproton != nelectron) {
        printf("nproton != nelectron\n");
        int n;
        ParticleTypes type;
        if (nproton > nelectron) {
            n = nproton - nelectron;
            type = ELECTRON;
        }
        else {
            n = nelectron - nproton;
            type = PROTON;
        }
        int i = 0;
        while (n > 0) {
            Particle* particle = createParticle(particlesNumber, i, 0, 0, weight, type, temperature);
            particles.push_back(particle);
            particlesNumber++;
            ++i;
            n--;
            if (i >= xnumber) {
                i = 0;
            }
            if (type == PROTON) {
                nproton++;
            }
            else {
                nelectron++;
            }
        }
    }
    if (nproton != nelectron) {
        printf("nproton != nelectron\n");
        errorLogFile = fopen("./output/errorLog.dat", "w");
        fprintf(errorLogFile, "nproton = %d nelectron = %d\n", nproton, nelectron);
        fclose(errorLogFile);
        exit(0);
    }*/
    createParticles();

    double chargeDensityAmplitude = epsilon * concentration * electron_charge_normalized;
    double Eamplitude = -4 * pi * chargeDensityAmplitude / (kw * fieldScale);

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k].x = Eamplitude * sin(kw * xgrid[i]);
                Efield[i][j][k].y = 0;
                Efield[i][j][k].z = 0;

                tempEfield[i] = Efield[i];
                explicitEfield[i] = Efield[i];
            }
        }
    }

    for (int j = 0; j < ynumber + 1; ++j) {
        for (int k = 0; k < znumber + 1; ++k) {
            Efield[xnumber][j][k] = Efield[0][j][k];
            tempEfield[xnumber][j][k] = tempEfield[0][j][k];
            explicitEfield[xnumber][j][k] = explicitEfield[0][j][k];
        }
    }

    double Vamplitude = electron_charge_normalized * Eamplitude * fieldScale / (massElectron * omega);

    for (int pcount = 0; pcount < particles.size(); ++pcount) {
        Particle *particle = particles[pcount];
        if (particle->type == ELECTRON) {
            Vector3d velocity = Vector3d(1, 0, 0) * Vamplitude * cos(kw * particle->coordinates.x);
            particle->addVelocity(velocity, speed_of_light_normalized);
        }
    }
}

void Simulation::initializeFluxFromRight() {
    boundaryConditionType = SUPER_CONDUCTOR_LEFT;
    createParticles();
    E0 = E0 - V0.vectorMult(B0) / (speed_of_light_normalized);
    //initializeAlfvenWaveY(10, 1.0E-4);
    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k] = E0;
                if (i == 0) {
                    Efield[i][j][k].y = 0;
                    Efield[i][j][k].z = 0;
                }
                tempEfield[i][j][k] = Efield[i][j][k];
                newEfield[i][j][k] = Efield[i][j][k];
                explicitEfield[i][j][k] = Efield[i][j][k];
            }
        }
    }

    //fieldsLorentzTransitionX(V0.x);

    for (int i = 0; i < particles.size(); ++i) {
        Particle *particle = particles[i];
        particle->addVelocity(V0, speed_of_light_normalized);
    }

    double magneticEnergy = B0.scalarMult(B0) / (8 * pi);
    double kineticEnergy = density * V0.scalarMult(V0) / 2;

    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    fprintf(informationFile, "magneticEnergy/kineticEnergy = %15.10g\n", magneticEnergy / kineticEnergy);
    printf("magneticEnergy/kinetikEnergy = %15.10g\n", magneticEnergy / kineticEnergy);
    fclose(informationFile);

}

void Simulation::fieldsLorentzTransitionX(const double &v) {
    double gamma = 1.0 / sqrt(1 - v * v / speed_of_light_normalized_sqr);
    for (int i = 0; i < xnumber; ++i) {
        int prevI = i - 1;
        if (prevI < 0) {
            if (boundaryConditionType == PERIODIC) {
                prevI = xnumber - 1;
            } else {
                prevI = 0;
            }
        }
        for (int j = 0; j < ynumber; ++j) {
            int prevJ = j - 1;
            if (prevJ < 0) {
                prevJ = ynumber - 1;
            }
            for (int k = 0; k < znumber; ++k) {
                int prevK = k - 1;
                if (prevK < 0) {
                    prevK = znumber - 1;
                }
                Vector3d middleB = (Bfield[prevI][prevJ][prevK] + Bfield[prevI][j][prevK] + Bfield[prevI][prevJ][k] + Bfield[prevI][j][k]
                                    + Bfield[i][prevJ][prevK] + Bfield[i][j][prevK] + Bfield[i][prevJ][k] + Bfield[i][j][k]) * 0.125;
                newEfield[i][j][k].y = gamma * (Efield[i][j][k].y - v * middleB.z / speed_of_light_normalized);
                newEfield[i][j][k].z = gamma * (Efield[i][j][k].z + v * middleB.y / speed_of_light_normalized);
            }
        }
    }
    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Vector3d middleE = (Efield[i][j][k] + Efield[i][j + 1][k] + Efield[i][j][k + 1] + Efield[i][j + 1][k + 1]
                                    + Efield[i + 1][j][k] + Efield[i + 1][j + 1][k] + Efield[i + 1][j][k + 1] + Efield[i + 1][j + 1][k + 1]) * 0.125;
                newBfield[i][j][k].y = gamma * (Bfield[i][j][k].y + v * middleE.z / speed_of_light_normalized);
                newBfield[i][j][k].z = gamma * (Bfield[i][j][k].z - v * middleE.y / speed_of_light_normalized);
            }
        }
    }

    for (int i = 0; i < xnumber; ++i) {
        for (int k = 0; k < znumber; ++k) {
            newEfield[i][ynumber][k] = newEfield[i][0][k];
        }
    }

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            newEfield[i][j][znumber] = newEfield[i][j][0];
        }
    }

    if (boundaryConditionType == PERIODIC) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                newEfield[xnumber][j][k] = newEfield[0][j][k];
            }
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k] = newEfield[i][j][k];
                tempEfield[i][j][k] = newEfield[i][j][k];
                explicitEfield[i][j][k] = newEfield[i][j][k];
            }
        }
    }

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Bfield[i][j][k] = newBfield[i][j][k];
            }
        }
    }
}

void Simulation::initializeShockWave() {
    boundaryConditionType = FREE_BOTH;

    E0 = Vector3d(0, 0, 0);
    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k] = Vector3d(0, 0, 0);
                tempEfield[i][j][k] = Efield[i][j][k];
                newEfield[i][j][k] = Efield[i][j][k];
                explicitEfield[i][j][k] = Efield[i][j][k];
            }
        }
    }
    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Bfield[i][j][k] = Vector3d(B0.x, 0, 0);
                newBfield[i][j][k] = Bfield[i][j][k];
            }
        }
    }

    printf("creating particles\n");
    double concentration = density / (massProton + massElectron);
    double downstreamTemperature = 1000 * temperature;
    double upstreamTemperature = temperature;
    Vector3d upstreamVelocity = V0;
    Vector3d downstreamVelocity = Vector3d(V0.x / 4, 0, 0);
    double alfvenV = B0.norm() / sqrt(4 * pi * density);
    double soundVelectron = sqrt(5 * kBoltzman_normalized * downstreamTemperature / (3 * massElectron));

    if (alfvenV > V0.norm()) {
        printf("alfvenV > V0\n");
    }

    printf("alfvenV/V0 = %15.10g\n", alfvenV / V0.norm());

    if (soundVelectron > V0.norm()) {
        printf("soundV > V0\n");
    }
    printf("soundV/V0 = %15.10g\n", soundVelectron / V0.norm());
    //Vector3d downstreamVelocity = Vector3d(0, 0, 0);
    shockWavePoint = xnumber / 2;
    int n = 0;
    for (int i = 0; i < xnumber; ++i) {
        for (int typeCounter = 0; typeCounter < typesNumber; ++typeCounter) {
            double x = xgrid[i] + 0.0001 * deltaX;
            int localParticlesPerBin = types[typeCounter].particlesPerBin;
            double localTemperature = upstreamTemperature;
            if (i < shockWavePoint) {
                localParticlesPerBin = localParticlesPerBin * 4;
                localTemperature = upstreamTemperature;
            }
            double deltaXParticles = deltaX / localParticlesPerBin;
            for (int j = 0; j < ynumber; ++j) {
                for (int k = 0; k < znumber; ++k) {
                    double weight = (types[typeCounter].concentration / types[typeCounter].particlesPerBin) * volumeB(i,
                                                                                                                      j,
                                                                                                                      k);
                    for (int l = 0; l < localParticlesPerBin; ++l) {
                        ParticleTypes type = types[typeCounter].type;
                        Particle *particle = createParticle(n, i, j, k, weight, type, types[typeCounter],
                                                            localTemperature,localTemperature,localTemperature);
                        //particle->x = middleXgrid[i];
                        n++;
                        /*if (l % 2 == 0) {
                            x = particle->x;
                        } else {
                            particle->x= x;
                        }*/
                        if (i >= shockWavePoint) {
                            particle->addVelocity(upstreamVelocity, speed_of_light_normalized);
                        } else {
                            particle->addVelocity(downstreamVelocity, speed_of_light_normalized);
                        }
                        particle->coordinates.x = x + deltaXParticles * l;
                        particle->initialMomentum = particle->momentum;
                        particles.push_back(particle);
                        particlesNumber++;
                        if (particlesNumber % 1000 == 0) {
                            printf("create particle number %d\n", particlesNumber);
                        }
                    }
                }
            }
        }
    }

    /*for(int i = 0; i < shockWavePoint; ++i){
        //Bfield[i].y = B0.x*sin(2*20*pi*middleXgrid[i]/xsize);
        double amplitude = 0.1*B0.x;
        Bfield[i].y = amplitude*(uniformDistribution() - 0.5);
        Bfield[i].z = amplitude*(uniformDistribution() - 0.5);
        newBfield[i] = Bfield[i];
    }*/

    initializeKolmogorovSpectrum(0, shockWavePoint);

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                double v = upstreamVelocity.x;
                if (i < shockWavePoint) {
                    //v = V0.x/4;
                    v = downstreamVelocity.x;
                }
                double gamma = 1.0 / sqrt(1 - v * v / speed_of_light_normalized_sqr);
                Vector3d middleE = (Efield[i][j][k] + Efield[i + 1][j][k]) * 0.5;
                newBfield[i][j][k].y = gamma * (Bfield[i][j][k].y + v * middleE.z / speed_of_light_normalized);
                newBfield[i][j][k].z = gamma * (Bfield[i][j][k].z - v * middleE.y / speed_of_light_normalized);
            }
        }
    }

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                double v = upstreamVelocity.x;
                if (i < shockWavePoint) {
                    //v = V0.x/4;
                    v = downstreamVelocity.x;
                }
                Vector3d middleE = (Efield[i][j][k] + Efield[i + 1][j][k]) * 0.5;
                double gamma = 1.0 / sqrt(1 - v * v / speed_of_light_normalized_sqr);
                newBfield[i][j][k].y = gamma * (Bfield[i][j][k].y + v * middleE.z / speed_of_light_normalized);
                newBfield[i][j][k].z = gamma * (Bfield[i][j][k].z - v * middleE.y / speed_of_light_normalized);
            }
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            newEfield[i][j][znumber] = newEfield[i][j][0];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int k = 0; k < znumber + 1; ++k) {
            newEfield[i][ynumber][k] = newEfield[i][0][k];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k] = newEfield[i][j][k];
                tempEfield[i][j][k] = newEfield[i][j][k];
                explicitEfield[i][j][k] = newEfield[i][j][k];
            }
        }
    }

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Bfield[i] = newBfield[i];
            }
        }
    }
}

void Simulation::initializeKolmogorovSpectrum(int start, int end) {
    double turbulenceFraction = 1.0;
    //use if defined shockWavePoint
    double length = xgrid[end] - xgrid[start];

    double minWaveLength = length / 200;
    double maxWaveLength = length / 10;

    int maxHarmonicNumber = length / minWaveLength;
    int minHarmonicNumber = length / maxWaveLength;

    for (int harmCounter = minHarmonicNumber; harmCounter <= maxHarmonicNumber; ++harmCounter) {
        double k = 2 * pi * harmCounter / length;
        double Bamplitude = turbulenceFraction * B0.x * power(k * length / (2 * pi), -5.0 / 6.0);
        double phiY = 2 * pi * uniformDistribution();
        double phiZ = 2 * pi * uniformDistribution();

        for (int i = start; i < end; ++i) {
            for (int j = 0; j < ynumber; ++j) {
                for (int k = 0; k < znumber; ++k) {
                    Bfield[i][j][k].y += Bamplitude * sin(k * middleXgrid[i] + phiY);
                    Bfield[i][j][k].z += Bamplitude * sin(k * middleXgrid[i] + phiZ);
                    newBfield[i][j][k] = Bfield[i][j][k];
                }
            }
        }
    }
}

void Simulation::initializeTwoStream() {
    boundaryConditionType = PERIODIC;
    createParticles();
    collectParticlesIntoBins();
    double u = speed_of_light_normalized / 5;
    Vector3d electronsVelocityPlus = Vector3d(0, u, 0);
    Vector3d electronsVelocityMinus = Vector3d(0, -u, 0);
    B0 = Vector3d(0, 0, 0);

    double Bamplitude = 1E-12 * (plasma_period * sqrt(scaleFactor));
    //Bamplitude = 0;
    double Eamplitude = 0;

    checkDebyeParameter();

    double kw = 2 * pi / xsize;

    informationFile = fopen((outputDir + "information.dat").c_str(), "a");

    if (xsize * omegaPlasmaElectron / speed_of_light_normalized < 5) {
        printf("xsize*omegaPlasmaElectron/speed_of_light_normalized < 5\n");
        fprintf(informationFile, "xsize*omegaPlasmaElectron/speed_of_light_normalized < 5\n");
    }
    printf("xsize*omegaPlasmaElectron/speed_of_light_normalized = %g\n",
           xsize * omegaPlasmaElectron / speed_of_light_normalized);
    fprintf(informationFile, "xsize*omegaPlasmaElectron/speed_of_light_normalized = %g\n",
            xsize * omegaPlasmaElectron / speed_of_light_normalized);

    if (deltaX * omegaPlasmaElectron / speed_of_light_normalized > 0.2) {
        printf("deltaX*omegaPlasmaElectron/speed_of_light_normalized > 0.2\n");
        fprintf(informationFile, "deltaX*omegaPlasmaElectron/speed_of_light_normalized > 0.2\n");
    }
    printf("deltaX*omegaPlasmaElectron/speed_of_light_normalized = %g\n",
           xsize * omegaPlasmaElectron / speed_of_light_normalized);
    fprintf(informationFile, "deltaX*omegaPlasmaElectron/speed_of_light_normalized = %g\n",
            xsize * omegaPlasmaElectron / speed_of_light_normalized);

    if (kw > omegaPlasmaElectron / u) {
        printf("k > omegaPlasmaElectron/u\n");
        fprintf(informationFile, "k > omegaPlasmaElectron/u\n");
    }

    printf("k u/omegaPlasmaElectron = %g\n", kw * u / omegaPlasmaElectron);
    fprintf(informationFile, "k u/omegaPlasmaElectron = %g\n", kw * u / omegaPlasmaElectron);
    fclose(informationFile);

    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                Bfield[i][j][k] = Vector3d(0, 1, 0) * Bamplitude * cos(kw * middleXgrid[i]);
                newBfield[i][j][k] = Bfield[i][j][k];
            }
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        for (int j = 0; j < ynumber + 1; ++j) {
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k] = Vector3d(0, 0, 1) * Eamplitude * cos(kw * xgrid[i]);
                tempEfield[i][j][k] = Efield[i][j][k];
                newEfield[i][j][k] = Efield[i][j][k];
                explicitEfield[i][j][k] = Efield[i][j][k];
            }
        }
    }
    int electronCount = 0;
    for (int pcount = 0; pcount < particles.size(); ++pcount) {
        Particle *particle = particles[pcount];
        if (particle->type == ELECTRON) {
            if (electronCount % 2 == 0) {
                particle->addVelocity(electronsVelocityPlus, speed_of_light_normalized);
            } else {
                particle->addVelocity(electronsVelocityMinus, speed_of_light_normalized);
            }
            electronCount++;
        }
    }
}

void Simulation::initializeExternalFluxInstability() {
    boundaryConditionType = PERIODIC;
    createParticles();
    double alfvenV = B0.norm() * fieldScale / sqrt(4 * pi * density);
    double concentration = density / (massProton + massElectron);
    double phaseV = 2 * alfvenV;
    double kw = 2 * pi / xsize;
    double omega = kw * phaseV;

    extJ = 0.001 * electron_charge_normalized * alfvenV * concentration;
    checkDebyeParameter();

    double Byamplitude = 4 * pi * sqr(alfvenV) * extJ / (speed_of_light_normalized * kw * (sqr(phaseV) - sqr(
        alfvenV))) / (plasma_period * sqrt(scaleFactor));
    double Uyamplitude = B0.norm() * fieldScale * phaseV * extJ / (kw * density * speed_of_light_normalized * (sqr(
        phaseV) - sqr(alfvenV))) * (scaleFactor / plasma_period);
    double cyclothronOmegaElectron = electron_charge_normalized * B0.norm() * fieldScale / (massElectron * speed_of_light_normalized);
    double cyclothronOmegaProton = electron_charge_normalized * B0.norm() * fieldScale / (massProton * speed_of_light_normalized);


    checkGyroRadius();
    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    fprintf(informationFile, "alfven V = %g\n", alfvenV * scaleFactor / plasma_period);
    fprintf(informationFile, "phase V = %g\n", phaseV * scaleFactor / plasma_period);
    fprintf(informationFile, "alfven V/c = %g\n", alfvenV / speed_of_light_normalized);
    fprintf(informationFile, "phase V/c = %g\n", phaseV / speed_of_light_normalized);
    fprintf(informationFile, "external flux = %g\n", extJ * plasma_period * plasma_period * sqrt(scaleFactor));
    fprintf(informationFile, "By max amplitude = %g\n", Byamplitude);
    fprintf(informationFile, "Uy max amplitude = %g\n", Uyamplitude);
    fprintf(informationFile, "omega/omega_plasma = %g\n", omega);
    if (omega > cyclothronOmegaProton) {
        printf("omega > cyclothron Omega Proton\n");
        fprintf(informationFile, "omega > cyclothron Omega Proton\n");
    } else if (omega > cyclothronOmegaProton / 100.0) {
        printf("omega > cyclothrone Omega Proton/100\n");
        fprintf(informationFile, "omega > cyclothron Omega Proton/100\n");
    }
    printf("omega/cyclothronOmega = %g\n", omega / cyclothronOmegaProton);
    fprintf(informationFile, "omega/cyclothronOmega = %g\n", omega / cyclothronOmegaProton);

    fclose(informationFile);
}

void Simulation::initializeAnisotropic() {
    boundaryConditionType = PERIODIC;
    /*double Tx = temperature*10;
    double Ty = temperature*1000;
    double Tz = temperature*1000;*/

    double temperatureIon = 6.9E7;
    double temperatureElectron = 6.9E7;
    double temperatureParallelMinority = 4.63E8;
    double temperatureNormalMinority = 115.9E8;

    for(int i = 0; i < typesNumber; ++i){
        types[i].temperatureX = temperatureIon;
        types[i].temperatureY = temperatureIon;
        types[i].temperatureZ = temperatureIon;
    }

    types[0].temperatureX = temperatureElectron;
    types[0].temperatureY = temperatureElectron;
    types[0].temperatureZ = temperatureElectron;

    types[5].temperatureX = temperatureParallelMinority;
    types[5].temperatureY = temperatureNormalMinority;
    types[5].temperatureZ = temperatureNormalMinority;

    createParticles();

    omegaPlasmaProton = sqrt(
        4 * pi * types[1].concentration * electron_charge_normalized * electron_charge_normalized / massProton);
    omegaPlasmaElectron = sqrt(
        4 * pi * types[0].concentration * electron_charge_normalized * electron_charge_normalized / massElectron);
    double omegaPlasmaAlpha = sqrt(
        4 * pi * types[0].concentration * electron_charge_normalized * electron_charge_normalized / massAlpha);
    omegaGyroProton = B0.norm() * fieldScale * electron_charge_normalized / (massProton * speed_of_light_normalized);
    omegaGyroElectron = B0.norm() * fieldScale * electron_charge_normalized / (massElectron * speed_of_light_normalized);
    double omegaGyroAlpha = B0.norm() * fieldScale * electron_charge_normalized / (massAlpha * speed_of_light_normalized);

    double kmin = 2*pi/xsize;
    double kmax = pi/deltaX;}

void Simulation::createArrays() {
    printf("creating arrays\n");
    printLog("creating arrays\n");
    xgrid = new double[xnumber + 1];
    ygrid = new double[ynumber + 1];
    zgrid = new double[znumber + 1];

    middleXgrid = new double[xnumber];
    middleYgrid = new double[ynumber];
    middleZgrid = new double[znumber];

    Efield = new Vector3d **[xnumber + 1];
    newEfield = new Vector3d **[xnumber + 1];
    tempEfield = new Vector3d **[xnumber + 1];
    explicitEfield = new Vector3d **[xnumber + 1];
    rotB = new Vector3d **[xnumber + 1];
    Ederivative = new Vector3d **[xnumber + 1];
    Bfield = new Vector3d **[xnumber];
    newBfield = new Vector3d **[xnumber];

    for (int i = 0; i < xnumber; ++i) {
        Bfield[i] = new Vector3d *[ynumber];
        newBfield[i] = new Vector3d *[ynumber];
        for (int j = 0; j < ynumber; ++j) {
            Bfield[i][j] = new Vector3d[znumber];
            newBfield[i][j] = new Vector3d[znumber];
            for (int k = 0; k < znumber; ++k) {
                Bfield[i][j][k] = Vector3d(0, 0, 0);
                newBfield[i][j][k] = Vector3d(0, 0, 0);
            }
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        Efield[i] = new Vector3d *[ynumber + 1];
        newEfield[i] = new Vector3d *[ynumber + 1];
        tempEfield[i] = new Vector3d *[ynumber + 1];
        explicitEfield[i] = new Vector3d *[ynumber + 1];
        rotB[i] = new Vector3d *[ynumber + 1];
        Ederivative[i] = new Vector3d *[ynumber + 1];
        for (int j = 0; j < ynumber + 1; ++j) {
            Efield[i][j] = new Vector3d[znumber + 1];
            newEfield[i][j] = new Vector3d[znumber + 1];
            tempEfield[i][j] = new Vector3d[znumber + 1];
            explicitEfield[i][j] = new Vector3d[znumber + 1];
            rotB[i][j] = new Vector3d[znumber + 1];
            Ederivative[i][j] = new Vector3d[znumber + 1];
            for (int k = 0; k < znumber + 1; ++k) {
                Efield[i][j][k] = Vector3d(0, 0, 0);
                newEfield[i][j][k] = Vector3d(0, 0, 0);
                tempEfield[i][j][k] = Vector3d(0, 0, 0);
                explicitEfield[i][j][k] = Vector3d(0, 0, 0);
                rotB[i][j][k] = Vector3d(0, 0, 0);
                Ederivative[i][j][k] = Vector3d(0, 0, 0);
            }
        }
    }

    maxwellEquationMatrix = new std::vector<MatrixElement> ***[xnumber];
    maxwellEquationRightPart = new double ***[xnumber];
    for (int i = 0; i < xnumber; ++i) {
        maxwellEquationMatrix[i] = new std::vector<MatrixElement> **[ynumber];
        maxwellEquationRightPart[i] = new double **[ynumber];
        for (int j = 0; j < ynumber; ++j) {
            maxwellEquationMatrix[i][j] = new std::vector<MatrixElement> *[znumber];
            maxwellEquationRightPart[i][j] = new double *[znumber];
            for (int k = 0; k < znumber; ++k) {
                maxwellEquationMatrix[i][j][k] = new std::vector<MatrixElement>[maxwellEquationMatrixSize];
                maxwellEquationRightPart[i][j][k] = new double[maxwellEquationMatrixSize];
            }
        }
    }

    printf("creating arrays for divergence\n");
    printLog("creating arrays for divergence\n");

    divergenceCleanUpMatrix = new std::vector<MatrixElement> ***[xnumber];
    divergenceCleanUpRightPart = new double ***[xnumber];

    divergenceCleaningField = new double ***[xnumber];
    divergenceCleaningPotential = new double ***[xnumber];
    divergenceCleaningPotentialFourier = new double **[xnumber];

    for (int i = 0; i < xnumber; ++i) {
        divergenceCleanUpMatrix[i] = new std::vector<MatrixElement> **[ynumber];
        divergenceCleanUpRightPart[i] = new double **[ynumber];
        divergenceCleaningField[i] = new double **[ynumber];
        divergenceCleaningPotential[i] = new double **[ynumber];
        divergenceCleaningPotentialFourier[i] = new double *[ynumber];
        for (int j = 0; j < ynumber; ++j) {
            divergenceCleanUpMatrix[i][j] = new std::vector<MatrixElement> *[znumber];
            divergenceCleanUpRightPart[i][j] = new double *[znumber];
            divergenceCleaningField[i][j] = new double *[znumber];
            divergenceCleaningPotential[i][j] = new double *[znumber];
            divergenceCleaningPotentialFourier[i][j] = new double[znumber];
            for (int k = 0; k < znumber; ++k) {
                divergenceCleaningField[i][j][k] = new double[3];
                divergenceCleaningPotential[i][j][k] = new double[1];
                divergenceCleanUpMatrix[i][j][k] = new std::vector<MatrixElement>[3];
                divergenceCleanUpRightPart[i][j][k] = new double[3];
                divergenceCleaningPotentialFourier[i][j][k] = 0;
            }
        }
    }

    particlesInBbin = new std::vector<Particle *> **[xnumber];
    particlesInEbin = new std::vector<Particle *> **[xnumber + 1];

    for (int i = 0; i < xnumber; ++i) {
        particlesInBbin[i] = new std::vector<Particle *> *[ynumber];
        for (int j = 0; j < ynumber; ++j) {
            particlesInBbin[i][j] = new std::vector<Particle *>[znumber];
        }
    }

    for (int i = 0; i < xnumber + 1; ++i) {
        particlesInEbin[i] = new std::vector<Particle *> *[ynumber + 1];
        for (int j = 0; j < ynumber + 1; ++j) {
            particlesInEbin[i][j] = new std::vector<Particle *>[znumber + 1];
        }
    }

    electronConcentration = new double **[xnumber];
    protonConcentration = new double **[xnumber];
    chargeDensity = new double **[xnumber];
    velocityBulkProton = new Vector3d **[xnumber];
    velocityBulkElectron = new Vector3d **[xnumber];
    electricDensity = new double **[xnumber];
    pressureTensor = new Matrix3d **[xnumber];

    for (int i = 0; i < xnumber; ++i) {
        electronConcentration[i] = new double *[ynumber];
        protonConcentration[i] = new double *[ynumber];
        chargeDensity[i] = new double *[ynumber];
        velocityBulkProton[i] = new Vector3d *[ynumber];
        velocityBulkElectron[i] = new Vector3d *[ynumber];
        electricDensity[i] = new double *[ynumber];
        pressureTensor[i] = new Matrix3d *[ynumber];
        for (int j = 0; j < ynumber; ++j) {
            electronConcentration[i][j] = new double[znumber];
            protonConcentration[i][j] = new double[znumber];
            chargeDensity[i][j] = new double[znumber];
            velocityBulkProton[i][j] = new Vector3d[znumber];
            velocityBulkElectron[i][j] = new Vector3d[znumber];
            electricDensity[i][j] = new double[znumber];
            pressureTensor[i][j] = new Matrix3d[znumber];
            for (int k = 0; k < znumber; ++k) {
                electronConcentration[i][j][k] = 0;
                protonConcentration[i][j][k] = 0;
                chargeDensity[i][j][k] = 0;
                velocityBulkProton[i][j][k] = Vector3d(0, 0, 0);
                velocityBulkElectron[i][j][k] = Vector3d(0, 0, 0);

                electricDensity[i][j][k] = 0;
                pressureTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
            }
        }
    }

    electricFlux = new Vector3d **[xnumber + 1];
    dielectricTensor = new Matrix3d **[xnumber + 1];
    externalElectricFlux = new Vector3d **[xnumber + 1];
    divPressureTensor = new Vector3d **[xnumber + 1];

    for (int i = 0; i < xnumber + 1; ++i) {
        electricFlux[i] = new Vector3d *[ynumber + 1];
        dielectricTensor[i] = new Matrix3d *[ynumber + 1];
        externalElectricFlux[i] = new Vector3d *[ynumber + 1];
        divPressureTensor[i] = new Vector3d *[ynumber + 1];
        for (int j = 0; j < ynumber + 1; ++j) {
            electricFlux[i][j] = new Vector3d[znumber + 1];
            dielectricTensor[i][j] = new Matrix3d[znumber + 1];
            externalElectricFlux[i][j] = new Vector3d[znumber + 1];
            divPressureTensor[i][j] = new Vector3d[znumber + 1];
            for (int k = 0; k < znumber + 1; ++k) {
                electricFlux[i][j][k] = Vector3d(0, 0, 0);
                divPressureTensor[i][j][k] = Vector3d(0, 0, 0);
                dielectricTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
            }
        }
    }
}

//add new particle types here
void Simulation::createParticleTypes(double *concentrations, int *particlesPerBin) {
    types = new ParticleTypeContainer[typesNumber];
    ParticleTypeContainer type;

    type.type = ELECTRON;
    type.mass = massElectron;
    type.chargeCount = - 1;
    type.charge = -electron_charge_normalized;
    type.particlesPerBin = particlesPerBin[0];
    type.concentration = concentrations[0];

    types[0] = type;

    type.type = PROTON;
    type.mass = massProton;
    type.chargeCount = 1;
    type.charge = electron_charge_normalized;
    type.particlesPerBin = particlesPerBin[1];
    type.concentration = concentrations[1];

    types[1] = type;

    type.type = POSITRON;
    type.mass = massElectron;
    type.chargeCount = 1;
    type.charge = electron_charge_normalized;
    type.particlesPerBin = particlesPerBin[2];
    type.concentration = concentrations[2];

    types[2] = type;

    type.type = ALPHA;
    type.mass = massAlpha;
    type.chargeCount = 2;
    type.charge = 2 * electron_charge_normalized;
    type.particlesPerBin = particlesPerBin[3];
    type.concentration = concentrations[3];

    types[3] = type;

    type.type = DEUTERIUM;
    type.mass = massDeuterium;
    type.chargeCount = 1;
    type.charge = electron_charge_normalized;
    type.particlesPerBin = particlesPerBin[4];
    type.concentration = concentrations[4];

    types[4] = type;

    type.type = HELIUM3;
    type.mass = massHelium3;
    type.chargeCount = 2;
    type.charge = 2 * electron_charge_normalized;
    type.particlesPerBin = particlesPerBin[5];
    type.concentration = concentrations[5];

    types[5] = type;

    for (int i = 0; i < typesNumber; ++i) {
        if (types[i].particlesPerBin > 0) {
            types[i].particesDeltaX = deltaX / types[i].particlesPerBin;
            types[i].particesDeltaY = deltaY / types[i].particlesPerBin;
            types[i].particesDeltaZ = deltaZ / types[i].particlesPerBin;
        } else {
            types[i].particesDeltaX = xsize;
            types[i].particesDeltaY = ysize;
            types[i].particesDeltaZ = zsize;
        }
        types[i].temperatureX = temperature;
        types[i].temperatureY = temperature;
        types[i].temperatureZ = temperature;
    }
}

void Simulation::createFiles() {
    printf("creating files\n");
    fflush(stdout);
    FILE *logFile = fopen((outputDir + "log.dat").c_str(), "w");
    fflush(logFile);
    fclose(logFile);
    printLog("creatingFiles\n");
    protonTraectoryFile = fopen((outputDir + "trajectory_proton.dat").c_str(), "w");
    fclose(protonTraectoryFile);
    electronTraectoryFile = fopen((outputDir + "trajectory_electron.dat").c_str(), "w");
    fclose(electronTraectoryFile);
    distributionFileProton = fopen((outputDir + "distribution_protons.dat").c_str(), "w");
    fclose(distributionFileProton);
    distributionFileElectron = fopen((outputDir + "distribution_electrons.dat").c_str(), "w");
    fclose(distributionFileElectron);
    anisotropyFileElectron = fopen((outputDir + "anisotropy_electrons.dat").c_str(), "w");
    fclose(anisotropyFileElectron);
    anisotropyFileProton = fopen((outputDir + "anisotropy_protons.dat").c_str(), "w");
    fclose(anisotropyFileProton);
    anisotropyFileAlpha = fopen((outputDir + "anisotropy_alphas.dat").c_str(), "w");
    fclose(anisotropyFileAlpha);
    anisotropyFilePositron = fopen((outputDir + "anisotropy_positrons.dat").c_str(), "w");
    fclose(anisotropyFilePositron);
    anisotropyFileDeuterium = fopen((outputDir + "anisotropy_deuterium.dat").c_str(), "w");
    fclose(anisotropyFileDeuterium);
    anisotropyFileHelium3 = fopen((outputDir + "anisotropy_helium3.dat").c_str(), "w");
    fclose(anisotropyFileHelium3);
    EfieldFile = fopen((outputDir + "Efield.dat").c_str(), "w");
    fclose(EfieldFile);
    BfieldFile = fopen((outputDir + "Bfield.dat").c_str(), "w");
    fclose(BfieldFile);
    velocityFile = fopen((outputDir + "velocity.dat").c_str(), "w");
    fclose(velocityFile);
    velocityElectronFile = fopen((outputDir + "velocity_electron.dat").c_str(), "w");
    fclose(velocityElectronFile);
    Xfile = fopen((outputDir + "Xfile.dat").c_str(), "w");
    fclose(Xfile);
    Yfile = fopen((outputDir + "Yfile.dat").c_str(), "w");
    fclose(Yfile);
    Zfile = fopen((outputDir + "Zfile.dat").c_str(), "w");
    fclose(Zfile);
    generalFile = fopen((outputDir + "general.dat").c_str(), "w");
    fclose(generalFile);
    densityFile = fopen((outputDir + "concentrations.dat").c_str(), "w");
    fclose(densityFile);
    divergenceErrorFile = fopen((outputDir + "divergence_error.dat").c_str(), "w");
    fclose(divergenceErrorFile);
    informationFile = fopen((outputDir + "information.dat").c_str(), "w");
    fclose(informationFile);
    fluxFile = fopen((outputDir + "flux.dat").c_str(), "w");
    fclose(fluxFile);
    rotBFile = fopen((outputDir + "rotBFile.dat").c_str(), "w");
    fclose(rotBFile);
    EderivativeFile = fopen((outputDir + "EderivativeFile.dat").c_str(), "w");
    fclose(EderivativeFile);
    dielectricTensorFile = fopen((outputDir + "dielectricTensorFile.dat").c_str(), "w");
    fclose(dielectricTensorFile);
    maxwellMatrixFile = fopen((outputDir + "maxwellMatrixFile.dat").c_str(), "w");
    fclose(maxwellMatrixFile);
    errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
    fclose(errorLogFile);
    particleProtonsFile = fopen((outputDir + "protons.dat").c_str(), "w");
    fclose(particleProtonsFile);
    particleElectronsFile = fopen((outputDir + "electrons.dat").c_str(), "w");
    fclose(particleElectronsFile);
    particlePositronsFile = fopen((outputDir + "positrons.dat").c_str(), "w");
    fclose(particlePositronsFile);
    particleAlphaFile = fopen((outputDir + "alphas.dat").c_str(), "w");
    fclose(particleAlphaFile);
    particleDeuteriumFile = fopen((outputDir + "deuterium.dat").c_str(), "w");
    fclose(particleDeuteriumFile);
    particleHelium3File = fopen((outputDir + "helium3.dat").c_str(), "w");
    fclose(particleHelium3File);
}

void Simulation::checkFrequency(double omega) {
    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    double cyclothronOmegaElectron = electron_charge_normalized * B0.norm() * fieldScale / (massElectron * speed_of_light_normalized);
    double cyclothronOmegaProton = electron_charge_normalized * B0.norm() * fieldScale / (massProton * speed_of_light_normalized);
    if (omega > cyclothronOmegaProton) {
        printf("omega > cyclothron Omega Proton\n");
        fprintf(informationFile, "omega > cyclothron Omega Proton\n");
    } else if (omega > cyclothronOmegaProton / 100.0) {
        printf("omega > cyclothrone Omega Proton/100\n");
        fprintf(informationFile, "omega > cyclothron Omega Proton/100\n");
    }
    printf("omega/cyclothronOmega = %g\n", omega / cyclothronOmegaProton);
    fprintf(informationFile, "omega/cyclothronOmega = %g\n", omega / cyclothronOmegaProton);

    if (omega > 1.0) {
        printf("omega > omega plasma\n");
        fprintf(informationFile, "omega > omega plasma\n");
    } else if (omega > 0.01) {
        printf("omega > omega plasma/100\n");
        fprintf(informationFile, "omega > omega plasma/100\n");
    }
    printf("omega/omega plasma = %g\n", omega);
    fprintf(informationFile, "omega/omega plasma = %g\n", omega);

    fclose(informationFile);
}

void Simulation::checkDebyeParameter() {
    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    double concentration = density / (massProton + massElectron);
    double weight = concentration * volumeB(0, 0, 0) / electronsPerBin;
    double superParticleCharge = electron_charge_normalized * weight;
    double superParticleConcentration = concentration / weight;
    double superParticleTemperature = temperature * weight;

    double debyeLength = 1 / sqrt(
        4 * pi * electron_charge_normalized * electron_charge_normalized * concentration / (kBoltzman_normalized * temperature));
    double debyeNumber = 4 * pi * cube(debyeLength) * concentration / 3;

    if (debyeLength > deltaX) {
        printf("debye length > deltaX\n");
        fprintf(informationFile, "debye length > deltaX\n");
    }
    printf("debye length/deltaX = %g\n", debyeLength / deltaX);
    fprintf(informationFile, "debye length/deltaX = %g\n", debyeLength / deltaX);

    if (debyeNumber < 1.0) {
        printf("debye number < 1\n");
        fprintf(informationFile, "debye number < 1\n");
    } else if (debyeNumber < 100.0) {
        printf("debye number < 100\n");
        fprintf(informationFile, "debye number < 100\n");
    }
    printf("debye number = %g\n", debyeNumber);
    fprintf(informationFile, "debye number = %g\n", debyeNumber);

    double superParticleDebyeLength = 1 / sqrt(
        4 * pi * superParticleCharge * superParticleCharge * superParticleConcentration / (kBoltzman_normalized * superParticleTemperature));
    double superParticleDebyeNumber = 4 * pi * cube(superParticleDebyeLength) * superParticleConcentration / 3;

    if (superParticleDebyeLength > deltaX) {
        printf("super particle debye length > deltaX\n");
        fprintf(informationFile, "super particle debye length > deltaX\n");
    }
    printf("super particle debye length/deltaX = %g\n", superParticleDebyeLength / deltaX);
    fprintf(informationFile, "super particle debye length/deltaX = %g\n", superParticleDebyeLength / deltaX);

    if (superParticleDebyeNumber < 1.0) {
        printf("superparticle debye number < 1\n");
        fprintf(informationFile, "superparticle debye number < 1\n");
    } else if (superParticleDebyeNumber < 100.0) {
        printf("superparticle debye number < 100\n");
        fprintf(informationFile, "superparticle debye number < 100\n");
    }
    printf("superparticle debye number = %g\n", superParticleDebyeNumber);
    fprintf(informationFile, "superparticle debye number = %g\n", superParticleDebyeNumber);
    fclose(informationFile);
}

void Simulation::checkGyroRadius() {
    informationFile = fopen((outputDir + "information.dat").c_str(), "a");

    if (B0.norm() > 0) {
        double thermalMomentumElectron = sqrt(
            massElectron * kBoltzman_normalized * temperature) + massElectron * V0.norm();
        double gyroRadiusElectron = thermalMomentumElectron * speed_of_light_normalized / (electron_charge_normalized * B0.norm());
        double thermalMomentumProton = sqrt(massProton * kBoltzman_normalized * temperature) + massProton * V0.norm();
        double gyroRadiusProton = thermalMomentumProton * speed_of_light_normalized / (electron_charge_normalized * B0.norm());
        if (deltaX > 0.5 * gyroRadiusElectron) {
            printf("deltaX > 0.5*gyroRadiusElectron\n");
            fprintf(informationFile, "deltaX > 0.5*gyroRadiusElectron\n");
        }

        printf("deltaX/gyroRadiusElectron = %g\n", deltaX / gyroRadiusElectron);
        fprintf(informationFile, "deltaX/gyroRadiusElectron = %g\n", deltaX / gyroRadiusElectron);

        if (xsize < 2 * gyroRadiusProton) {
            printf("xsize < 2*gyroRadiusProton\n");
            fprintf(informationFile, "xsize < 2*gyroRadiusProton\n");
        }

        printf("xsize/gyroRadiusProton= %g\n", xsize / gyroRadiusProton);
        fprintf(informationFile, "xsize/gyroRadiusProton = %g\n", xsize / gyroRadiusProton);
    }

    fclose(informationFile);
}

void Simulation::checkCollisionTime(double omega) {
    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    double concentration = density / (massProton + massElectron);
    double weight = concentration * volumeB(0, 0, 0) / electronsPerBin;
    double superParticleCharge = electron_charge_normalized * weight;
    double superParticleConcentration = concentration / weight;
    double superParticleTemperature = temperature * weight;

    double qLog = 15;
    double nuElectronIon = 4 * sqrt(2 * pi) * qLog * power(electron_charge_normalized, 4) * (concentration) / (3 * sqrt(
        massElectron) * sqrt(cube(kBoltzman_normalized * temperature)));
    double collisionlessParameter = omega / nuElectronIon;

    if (collisionlessParameter < 1.0) {
        printf("collisionlessParameter < 1\n");
        fprintf(informationFile, "collisionlessParameter < 1\n");
    } else if (collisionlessParameter < 100.0) {
        printf("collisionlessParameter < 100\n");
        fprintf(informationFile, "collisionlessParameter < 100\n");
    }
    printf("collisionlessParameter = %g\n", collisionlessParameter);
    fprintf(informationFile, "collisionlessParameter = %g\n", collisionlessParameter);

    double superParticleNuElectronIon = 4 * sqrt(2 * pi) * qLog * power(electron_charge_normalized * weight,
                                                                        4) * (superParticleConcentration) / (3 * sqrt(
        massElectron * weight) * sqrt(cube(kBoltzman_normalized * superParticleTemperature)));
    double superParticleCollisionlessParameter = omega / superParticleNuElectronIon;

    if (superParticleCollisionlessParameter < 1.0) {
        printf("superParticleCollisionlessParameter < 1\n");
        fprintf(informationFile, "superParticleCollisionlessParameter < 1\n");
    } else if (superParticleCollisionlessParameter < 100.0) {
        printf("superParticleCollisionlessParameter < 100\n");
        fprintf(informationFile, "superParticleCollisionlessParameter < 100\n");
    }
    printf("superParticleCollisionlessParameter = %g\n", superParticleCollisionlessParameter);
    fprintf(informationFile, "superParticleCollisionlessParameter = %g\n", superParticleCollisionlessParameter);
    fclose(informationFile);
}

void Simulation::checkMagneticReynolds(double v) {
    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    double concentration = density / (massProton + massElectron);
    double weight = concentration * volumeB(0, 0, 0) / electronsPerBin;
    double superParticleCharge = electron_charge_normalized * weight;
    double superParticleConcentration = concentration / weight;
    double superParticleTemperature = temperature * weight;

    double qLog = 15;
    double conductivity = (3 * massProton * sqrt(massElectron) * sqrt(cube(kBoltzman_normalized * temperature))) / (sqr(
        massElectron) * 4 * sqrt(2 * pi) * qLog * sqr(electron_charge_normalized));

    double magneticReynolds = conductivity * v * 0.1 * xsize;

    if (magneticReynolds < 1.0) {
        printf("magneticReynolds < 1\n");
        fprintf(informationFile, "magneticReynolds < 1\n");
    } else if (magneticReynolds < 100.0) {
        printf("magneticReynolds < 100\n");
        fprintf(informationFile, "magneticReynolds < 100\n");
    }
    printf("magnetic Reynolds = %g\n", magneticReynolds);
    fprintf(informationFile, "magnetic Reynolds = %g\n", magneticReynolds);

    double superParticleConductivity = (3 * massProton * weight * sqrt(massElectron * weight) * sqrt(
        cube(kBoltzman_normalized * superParticleTemperature))) / (sqr(massElectron * weight) * 4 * sqrt(
        2 * pi) * qLog * sqr(electron_charge_normalized * weight));

    double superParticleMagneticReynolds = superParticleConductivity * v * 0.1 * xsize;

    if (superParticleMagneticReynolds < 1.0) {
        printf("superParticleMagneticReynolds < 1\n");
        fprintf(informationFile, "superParticleMagneticReynolds < 1\n");
    } else if (superParticleMagneticReynolds < 100.0) {
        printf("superParticleMagneticReynolds < 100\n");
        fprintf(informationFile, "superParticleMagneticReynolds < 100\n");
    }
    printf("superparticle magnetic Reynolds = %g\n", superParticleMagneticReynolds);
    fprintf(informationFile, "superparticle magnetic Reynolds = %g\n", superParticleMagneticReynolds);
    fclose(informationFile);
}

void Simulation::checkDissipation(double k, double alfvenV) {
    informationFile = fopen((outputDir + "information.dat").c_str(), "a");
    double omega = k * alfvenV;

    double concentration = density / (massProton + massElectron);
    double weight = concentration * volumeB(0, 0, 0) / electronsPerBin;
    double superParticleCharge = electron_charge_normalized * weight;
    double superParticleConcentration = concentration / weight;
    double superParticleTemperature = temperature * weight;

    double qLog = 15;
    double conductivity = (3 * massProton * sqrt(massElectron) * sqrt(cube(kBoltzman_normalized * temperature))) / (sqr(
        massElectron) * 4 * sqrt(2 * pi) * qLog * sqr(electron_charge_normalized));

    double nuMagnetic = speed_of_light_normalized_sqr / (4 * pi * conductivity);

    double kdissipation = omega * omega * nuMagnetic / (2 * cube(alfvenV));

    if (kdissipation > k) {
        printf("kdissipation > k\n");
        fprintf(informationFile, "kdissipation > k\n");
    } else if (kdissipation > 0.1 * k) {
        printf("kdissipation > 0.1*k\n");
        fprintf(informationFile, "kdissipation > 0.1*k\n");
    }
    printf("kdissipation/k = %g\n", kdissipation / k);
    fprintf(informationFile, "kdissipation/k = %g\n", kdissipation / k);

    double superParticleConductivity = (3 * massProton * weight * sqrt(massElectron * weight) * sqrt(
        cube(kBoltzman_normalized * superParticleTemperature))) / (sqr(massElectron * weight) * 4 * sqrt(
        2 * pi) * qLog * sqr(electron_charge_normalized * weight));

    double superParticleNuMagnetic = speed_of_light_normalized_sqr / (4 * pi * superParticleConductivity);

    double superParticleKdissipation = omega * omega * superParticleNuMagnetic / (2 * cube(alfvenV));

    if (superParticleKdissipation > k) {
        printf("super particle kdissipation > k\n");
        fprintf(informationFile, "super particle kdissipation > k\n");
    } else if (superParticleKdissipation > 0.1 * k) {
        printf("super particle kdissipation > 0.1*k\n");
        fprintf(informationFile, "super particle kdissipation > 0.1*k\n");
    }
    printf("super particle kdissipation/k = %g\n", superParticleKdissipation / k);
    fprintf(informationFile, "super particle kdissipation/k = %g\n", superParticleKdissipation / k);
    fclose(informationFile);
}

void Simulation::createParticles() {
    printf("creating particles\n");
    printLog("creating particles\n");
    double concentration = density / (massProton + massElectron);
    int n = 0;
    for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                int maxParticlesPerBin = types[0].particlesPerBin;
                double x = xgrid[i] + 0.0001 * deltaX;
                double y = ygrid[j] + 0.0001 * deltaY;
                double z = zgrid[k] + 0.0001 * deltaZ;
                for (int l = 0; l < maxParticlesPerBin; ++l) {
                    for (int typeCounter = 0; typeCounter < typesNumber; ++typeCounter) {
                        double weight = (types[typeCounter].concentration / types[typeCounter].particlesPerBin) * volumeB(
                            i, j, k);
                        double deltaXParticles = types[typeCounter].particesDeltaX;
                        double deltaYParticles = types[typeCounter].particesDeltaY;
                        double deltaZParticles = types[typeCounter].particesDeltaZ;
                        if (l < types[typeCounter].particlesPerBin) {
                            ParticleTypes type = types[typeCounter].type;
                            Particle *particle = createParticle(n, i, j, k, weight, type, types[typeCounter],
                                                                types[typeCounter].temperatureX, types[typeCounter].temperatureY, types[typeCounter].temperatureZ);
                            n++;
                            particle->coordinates.x = x + deltaXParticles * l;
                            particle->coordinates.y = y + deltaYParticles * l;
                            particle->coordinates.z = z + deltaZParticles * l;
                            //particle->addVelocity(V0, speed_of_light_normalized);
                            particles.push_back(particle);
                            particlesNumber++;
                            if (particlesNumber % 1000 == 0) {
                                printf("create particle number %d\n", particlesNumber);
                            }
                        }
                    }
                }
            }
        }
    }

    if (preserveChargeLocal) {
        moveToPreserveChargeLocal();
    }
}

void Simulation::moveToPreserveChargeLocal() {
    Particle *electron = NULL;
    Particle *notElectron = NULL;
    Particle *particle;
    int electronCount = 0;
    int notElectronCount = 0;
    int particleCount = 0;

    while (particleCount < particles.size()) {
        particle = particles[particleCount];
        if (particle->type != ELECTRON) {
            notElectron = particle;
            notElectronCount = particleCount;

            int necessaryElectrons = 1;

            if (particle->type == ALPHA) {
                necessaryElectrons = 2;
            }

            while (necessaryElectrons > 0) {
                electron = NULL;
                while (electron == NULL) {
                    if (electronCount >= particles.size()) {
                        errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
                        printf("error in preserving charge\n");
                        fprintf(errorLogFile, "error in preserving charge\n");
                        fclose(errorLogFile);
                        exit(0);
                    }
                    if (particles[electronCount]->type == ELECTRON) {
                        electron = particles[electronCount];
                        electron->coordinates = notElectron->coordinates;
                    }
                    electronCount++;
                }
                necessaryElectrons--;
            }
        }
        particleCount++;
    }
}

void Simulation::addToPreserveChargeGlobal() {
    int n = particles.size();
    if (chargeBalance == 0) {
        return;
    }

    if (chargeBalance > 0) {
        double weight = (types[0].concentration / types[0].particlesPerBin) * volumeB(xnumber - 1, 0, 0);
        for (int i = 0; i < chargeBalance; ++i) {
            Particle *particle = createParticle(n, xnumber - 1, 0, 0, weight, types[0].type, types[0], types[0].temperatureX, types[0].temperatureY, types[0].temperatureZ);
            particle->coordinates.x = xgrid[xnumber] - 0.0001 * deltaX;
            particle->coordinates.y = ygrid[0] + ysize * uniformDistribution();
            particle->coordinates.z = zgrid[0] + zsize * uniformDistribution();
            theoreticalEnergy += particle->energy(speed_of_light_normalized) * particle->weight * sqr(
                scaleFactor / plasma_period);
            theoreticalMomentum += particle->momentum * particle->weight * scaleFactor / plasma_period;
            n++;
        }
    } else {
        int typeN = 1;
        for (int i = 1; i < typesNumber; ++i) {
            if ((types[i].particlesPerBin > types[typeN].particlesPerBin) && (types[i].chargeCount == 1)) {
                typeN = i;
            }
        }
        double weight = (types[typeN].concentration / types[typeN].particlesPerBin) * volumeB(xnumber - 1, 0, 0);
        for (int i = 0; i < -chargeBalance; ++i) {
            Particle *particle = createParticle(n, xnumber - 1, 0, 0, weight, types[typeN].type, types[typeN],
                                                types[typeN].temperatureX, types[typeN].temperatureY, types[typeN].temperatureZ);
            particle->coordinates.x = xgrid[xnumber] - 0.0001 * deltaX;
            particle->coordinates.y = ygrid[0] + ysize * uniformDistribution();
            particle->coordinates.z = zgrid[0] + zsize * uniformDistribution();
            theoreticalEnergy += particle->energy(speed_of_light_normalized) * particle->weight * sqr(
                scaleFactor / plasma_period);
            theoreticalMomentum += particle->momentum * particle->weight * scaleFactor / plasma_period;
            n++;
        }
    }

    chargeBalance = 0;
}

Particle *Simulation::getFirstProton() {
    for (int pcount = 0; pcount < particles.size(); ++pcount) {
        Particle *particle = particles[pcount];
        if (particle->type == PROTON) {
            return particle;
        }
    }
    return NULL;
}

Particle *Simulation::getFirstElectron() {
    for (int pcount = 0; pcount < particles.size(); ++pcount) {
        Particle *particle = particles[pcount];
        if (particle->type == ELECTRON) {
            return particle;
        }
    }
    return NULL;
}

Particle *Simulation::getLastProton() {
    for (int pcount = particles.size() - 1; pcount >= 0; --pcount) {
        Particle *particle = particles[pcount];
        if (particle->type == PROTON) {
            return particle;
        }
    }
    return NULL;
}

Particle *Simulation::getLastElectron() {
    for (int pcount = particles.size() - 1; pcount >= 0; --pcount) {
        Particle *particle = particles[pcount];
        if (particle->type == ELECTRON) {
            return particle;
        }
    }
    return NULL;
}

Particle *Simulation::getProton(int n) {
    int count = 0;
    for (int pcount = 0; pcount < particles.size(); ++pcount) {
        Particle *particle = particles[pcount];
        if (particle->type == PROTON) {
            count++;
            if (count == n) {
                return particle;
            }
        }
    }
    errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
    fprintf(errorLogFile, "can not find proton number %d\n", n);
    printf("can not find proton number %d\n", n);
    fclose(errorLogFile);
    exit(0);
    return NULL;
}

Particle *Simulation::getElectron(int n) {
    int count = 0;
    for (int pcount = 0; pcount < particles.size(); ++pcount) {
        Particle *particle = particles[pcount];
        if (particle->type == ELECTRON) {
            count++;
            if (count == n) {
                return particle;
            }
        }
    }
    errorLogFile = fopen((outputDir + "errorLog.dat").c_str(), "w");
    fprintf(errorLogFile, "can not find electron number %d\n", n);
    printf("can not find electron number %d\n", n);
    fclose(errorLogFile);
    exit(0);
    return NULL;
}


//anisotropy only in not relativistic case
Particle *Simulation::createParticle(int n, int i, int j, int k, double weight, ParticleTypes type,
                                     ParticleTypeContainer typeContainer, double localTemparatureX, double localTemparatureY, double localTemparatureZ) {
    int chargeCount = typeContainer.chargeCount;
    double charge = typeContainer.charge;
    double mass = typeContainer.mass;

    double x = xgrid[i] + deltaX * uniformDistribution();
    double y = ygrid[j] + deltaY * uniformDistribution();
    double z = zgrid[k] + deltaZ * uniformDistribution();

    double dx = deltaX / 10;
    double dy = deltaY / 10;
    double dz = deltaZ / 10;

    double energy = mass * speed_of_light_normalized_sqr;
    double p = 0;
    double px, py, pz = 0;

    double thetaParamter = kBoltzman_normalized * localTemparatureX / (mass * speed_of_light_normalized_sqr);

    if (thetaParamter < 0.1) {
        //energy = maxwellDistribution(localTemparature, kBoltzman_normalized);
        px = sqrt(mass*kBoltzman_normalized*localTemparatureX)*normalDistribution();
        py = sqrt(mass*kBoltzman_normalized*localTemparatureY)*normalDistribution();
        pz = sqrt(mass*kBoltzman_normalized*localTemparatureZ)*normalDistribution();
        p = sqrt(px*px + py*py + pz*pz);
    } else {
        energy = maxwellJuttnerDistribution(localTemparatureX, mass, speed_of_light_normalized, kBoltzman_normalized);
        p = sqrt(energy * energy - sqr(mass * speed_of_light_normalized_sqr)) / speed_of_light_normalized;


        //p = 0;

        pz = p * (2 * uniformDistribution() - 1);
        double phi = 2 * pi * uniformDistribution();
        double pnormal = sqrt(p * p - pz * pz);
        px = pnormal * cos(phi);
        py = pnormal * sin(phi);
    }


    Particle *particle = new Particle(n, mass, chargeCount, charge, weight, type, typeContainer, x, y, z, px, py, pz,
                                      dx, dy, dz);

    return particle;
}
