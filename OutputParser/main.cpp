#include <stdio.h>
#include <math.h>
#include <string>

void
outputEfield(int timePointsCount, FILE *XFile, FILE *YFile, FILE *ZFile, int Nx, int Ny, int Nz,
             FILE *EFile, FILE *fluxFile, FILE *gleEFileX, FILE *gleEFileY,
             FILE *gleEFileZ, FILE *gleFluxFileX, FILE *gleFluxFileY, FILE *gleFluxFileZ,
             int Nt, int xpoint, int ypoint, int zpoint, int *timePoints, double * Xgrid, double * middleX,
             double * Ygrid, double * middleY, double * Zgrid, double * middleZ) ;

void
outputBfield(int timePointsCount, int Nx, int Ny, int Nz, int Nt, int xpoint, int ypoint, int zpoint, int *timePoints,
             double *middleX, double *middleY, double *middleZ, FILE *BFile, FILE *divergenceFile, FILE *gleBFileX,
             FILE *gleBFileY, FILE *gleBFileZ, FILE *gleDivergenceFileX, FILE *gleDivergenceFileY,
             FILE *gleDivergenceFileZ);

int min2(int a, int b){
    if(a <= b){
        return a;
    } else {
        return b;
    }
}

int countNumbers(FILE* file){
    int result = 0;
    double a;
    int flag = fscanf(file, "%lf", &a);
    while(flag > 0){
        result++;
        flag = fscanf(file, "%lf", &a);
    }
    return result;
}

int countNumbersLine(FILE* file, int numbersPerLine) {
    int result = 0;
    double a;
    int flag = fscanf(file, "%lf", &a);
    for (int i = 0; i < numbersPerLine - 1; ++i) {
        flag = fscanf(file, "%lf", &a);
    }
    while(flag > 0){
        result++;
        flag = fscanf(file, "%lf", &a);
        for (int i = 0; i < numbersPerLine - 1; ++i) {
            flag = fscanf(file, "%lf", &a);
        }
    }
    return result;
}


int main() {
    const int timePointsCount = 3; //must be >= 2
    const int Np = 1000;
    //for MPI on cluster
    //std::string outputDir = "/home/uikam/romansky/PIC3d5/output/";
    //std::string inputDir = "/home/uikam/romansky/PIC3d5/input/";
    //std::string backupDir = "/home/uikam/romansky/PIC3d5/backup/";

    //for MPI on Ubuntu
    std::string outputDir = "/home/vadim/PICpp/PIC++/output/";
    //std::string inputDir = "/home/vadim/PIC++/trunk/PIC++/input/";
    //std::string backupDir = "/home/vadim/PIC++/trunk/PIC++/backup/";

    //for windowss
    //std::string outputDir = "C:/users/Vadik/Documents/Visual Studio 2010/Projects/v4/PIC++/output/";
    //std::string inputDir = "C:/users/Vadik/Documents/Visual Studio 2010/Projects/v4/PIC++/input/";
    //cstd::string backupDir = "C:/users/Vadik/Documents/Visual Studio 2010/Projects/v4/PIC++/backup/";

    FILE* XFile = fopen((outputDir + "Xfile.dat").c_str(),"r");
    FILE* YFile = fopen((outputDir + "Yfile.dat").c_str(),"r");
    FILE* ZFile = fopen((outputDir + "Zfile.dat").c_str(),"r");

    //FILE* protonsFile = fopen((inputDir + "protons.dat").c_str(), "r");
    //FILE* electronsFile = fopen((inputDir + "electrons.dat").c_str(), "r");
    //FILE* positronsFile = fopen((inputDir + "positrons.dat").c_str(), "r");
    //FILE* alphasFile = fopen((inputDir + "alphas.dat").c_str(), "r");


    //FILE* gleXFile = fopen("../output/Xfile.dat","w");
    //FILE* gleYFile = fopen("../output/Yfile.dat","w");
    //FILE* gleZFile = fopen("../output/Zfile.dat","w");

    //FILE* gleProtonsFile = fopen("../output/protons.dat", "w");
    //FILE* gleElectronsFile = fopen("../output/electrons.dat", "w");
    //FILE* glePositronsFile = fopen("../output/positrons.dat", "w");
    ///FILE* gleAlphasFile = fopen("../output/alphas.dat", "w");

    int Nx = countNumbers(XFile);
    int Ny = countNumbers(YFile);
    int Nz = countNumbers(ZFile);

    FILE* EFile = fopen((outputDir + "Efield.dat").c_str(), "r");
    FILE* fluxFile = fopen((outputDir + "flux.dat").c_str(), "r");

    FILE* gleEFileX = fopen("../output/Efield_x.dat", "w");
    FILE* gleEFileY = fopen("../output/Efield_y.dat", "w");
    FILE* gleEFileZ = fopen("../output/Efield_z.dat", "w");

    FILE* gleFluxFileX = fopen("../output/flux_x.dat", "w");
    FILE* gleFluxFileY = fopen("../output/flux_y.dat", "w");
    FILE* gleFluxFileZ = fopen("../output/flux_z.dat", "w");

    int NE = countNumbersLine(EFile, 3);
    rewind(XFile);
    rewind(YFile);
    rewind(ZFile);
    rewind(EFile);
    int Nt = NE/(Nx*Ny*Nz);
    Nt = 14;
    //int Nt = 10;

    int xpoint = Nx/2 - 1;
    int ypoint = Ny/2 - 1;
    int zpoint = Nz/2 - 1;

    int timePoints[timePointsCount];
    timePoints[0] = 0;
    timePoints[timePointsCount-1] = Nt - 1;
    int deltaNt = Nt/(timePointsCount - 1);
    for(int i = 1; i < timePointsCount - 1; ++i){
        timePoints[i] = timePoints[i-1] + deltaNt;
    }

    /////// arrays with dimensions Nx*Ny*Nz*Nt //////
    double *Xgrid;
    double *middleX;
    double *Ygrid;
    double *middleY;
    double *Zgrid;
    double *middleZ;
    Xgrid= new double[Nx];
    middleX= new double[Nx - 1];
    Ygrid= new double[Ny];
    middleY= new double[Ny - 1];
    Zgrid= new double[Nz];
    middleZ= new double[Nz - 1];
    /*outputEfield(timePointsCount, XFile, YFile, ZFile, Nx, Ny, Nz, EFile, fluxFile,
                 gleEFileX, gleEFileY, gleEFileZ, gleFluxFileX, gleFluxFileY, gleFluxFileZ,
                 Nt, xpoint, ypoint, zpoint, timePoints, Xgrid, middleX, Ygrid, middleY,
                 Zgrid, middleZ);*/


    FILE* BFile = fopen((outputDir + "Bfield.dat").c_str(),"r");
    //FILE* concentrationFile = fopen((inputDir + "concentrations.dat").c_str(), "r");
    //FILE* velocityElectronFile = fopen((inputDir + "velocity_electron.dat").c_str(), "r");
    //FILE* velocityProtonFile = fopen((inputDir + "velocity.dat").c_str(), "r");
    FILE* divergenceFile = fopen((outputDir + "divergence_error.dat").c_str(), "r");

    FILE* gleBFileX = fopen("../output/Bfield_x.dat","w");
    FILE* gleBFileY = fopen("../output/Bfield_y.dat","w");
    FILE* gleBFileZ = fopen("../output/Bfield_z.dat","w");
    //FILE* gleConcentrationFileX = fopen("../output/concentrations_x.dat", "w");
    //FILE* gleConcentrationFileY = fopen("../output/concentrations_y.dat", "w");
    //FILE* gleConcentrationFileZ = fopen("../output/concentrations_z.dat", "w");
    //FILE* gleVelocityFileX = fopen("../output/velocity_x.dat", "w");
    //FILE* gleVelocityFileY = fopen("../output/velocity_y.dat", "w");
    //FILE* gleVelocityFileZ = fopen("../output/velocity_z.dat", "w");
    FILE* gleDivergenceFileX = fopen("../output/divergence_error_x.dat", "w");
    FILE* gleDivergenceFileY = fopen("../output/divergence_error_y.dat", "w");
    FILE* gleDivergenceFileZ = fopen("../output/divergence_error_z.dat", "w");


    /////// arrays with dimensions (Nx - 1)*(Ny - 1)(Nz-1)*Nt //////
    /*outputBfield(timePointsCount, Nx, Ny, Nz, Nt, xpoint, ypoint, zpoint, timePoints, middleX, middleY, middleZ, BFile,
                 divergenceFile, gleBFileX, gleBFileY, gleBFileZ,
                 gleDivergenceFileX, gleDivergenceFileY, gleDivergenceFileZ);*/

    delete[] Xgrid;
    delete[] middleX;
    delete[] Ygrid;
    delete[] middleY;
    delete[] Zgrid;
    delete[] middleZ;


    fclose(BFile);
    //fclose(concentrationFile);
    //fclose(velocityElectronFile);
    //fclose(velocityProtonFile);
    fclose(divergenceFile);

    fclose(gleBFileX);
    fclose(gleBFileY);
    fclose(gleBFileZ);
    //fclose(gleConcentrationFileX);
    //fclose(gleConcentrationFileY);
    //fclose(gleConcentrationFileZ);
    //fclose(gleVelocityFileX);
    //fclose(gleVelocityFileY);
    //fclose(gleVelocityFileZ);
    fclose(gleDivergenceFileX);
    fclose(gleDivergenceFileY);
    fclose(gleDivergenceFileZ);

    /////// arrays with dimensions Nt //////

    FILE* generalFile = fopen((outputDir + "general.dat").c_str(), "r");
    FILE* incrementFile = fopen((outputDir + "increment.dat").c_str(), "r");
    //FILE* trajectoryProtonFile = fopen((inputDir + "trajectory_proton.dat").c_str(), "r");
    //FILE* trajectoryElectronFile = fopen((inputDir + "trajectory_electron.dat").c_str(), "r");

    FILE* gleGeneralFile = fopen("../output/general.dat", "w");
    //FILE* gleTrajectoryProtonFile = fopen("../output/trajectory_proton.dat", "w");
    //FILE* gleTrajectoryElectronFile = fopen("../output/trajectory_electron.dat", "w");

    double** general = new double*[Nt];
    double** protonTraectory = new double*[Nt];
    double** electronTraectory = new double*[Nt];
    for(int i = 0; i < Nt; ++i){
        general[i] = new double[38];
        protonTraectory[i] = new double[8];
        electronTraectory[i] = new double[8];
    }

    for(int i = 0; i < Nt; i++){
        for(int j = 0; j < 38; ++j){
        //for(int j = 0; j < 7; ++j){
            fscanf(generalFile, "%lf", &general[i][j]);
        }
        for(int j = 0; j < 8; ++j){
            //fscanf(trajectoryProtonFile, "%lf", &protonTraectory[i][j]);
            //fscanf(gleTrajectoryElectronFile, "%lf", &electronTraectory[i][j]);
        }
    }

    int NbestIncrement = min2(Nt - 1, 200);
    int Nsaturation = min2(NbestIncrement, 350);

    double gamma;
    fscanf(incrementFile, "%lf", &gamma);

    general[0][18] = general[NbestIncrement][6]/exp(2*gamma*general[NbestIncrement][1]);

    for(int i = 1; i < Nt; ++i){
        if(i <= Nsaturation){
            general[i][18] = general[0][18]*exp(2*gamma*general[i][1]);
        } else {
            general[i][18] = general[Nsaturation][18];
        }
    }

    for(int i = 0; i < Nt; i++){
        for(int j = 1; j < 19; ++j) {
        //for(int j = 1; j < 7; ++j) {
            fprintf(gleGeneralFile, "%20.15g", general[i][j]);
        }
        for(int j = 0; j < 8; ++j){
            //fprintf(gleTrajectoryProtonFile, "%20.15g", protonTraectory[i][j]);
            //fprintf(gleTrajectoryElectronFile, "%20.15g", electronTraectory[i][j]);
        }
        fprintf(gleGeneralFile, "\n");
        //fprintf(gleTrajectoryProtonFile, "\n");
        //fprintf(gleTrajectoryElectronFile, "\n");
    }

    for(int i = 0; i < Nt; ++i){
        delete[] general[i];
        delete[] protonTraectory[i];
        delete[] electronTraectory[i];
    }
    delete[] general;
    delete[] protonTraectory;
    delete[] electronTraectory;

    fclose(generalFile);
    fclose(incrementFile);
    //fclose(trajectoryElectronFile);
    //fclose(trajectoryProtonFile);

    fclose(gleGeneralFile);
    //fclose(gleTrajectoryElectronFile);
    //fclose(gleTrajectoryProtonFile);
    /////// arrays with dimensions Np*Nt //////

    FILE* distributionProtonsFile = fopen((outputDir + "distribution_protons.dat").c_str(), "r");
    FILE* distributionElectronsFile = fopen((outputDir + "distribution_electrons.dat").c_str(), "r");

    FILE* gleDistributionProtonsFile = fopen("../output/distribution_protons.dat", "w");
    FILE* gleDistributionElectronsFile = fopen("../output/distribution_electrons.dat", "w");

    double** distributionProtons = new double*[Np*Nt];
    double** distributionElectrons = new double*[Np*Nt];
    for(int i = 0; i < Np*Nt; ++i){
        distributionProtons[i] = new double[2];
        distributionElectrons[i] = new double[2];
    }

    for(int i = 0; i < Np*Nt; ++i){
        for(int j = 0; j < 2; ++j) {
            fscanf(distributionProtonsFile, "%lf", &distributionProtons[i][j]);
            fscanf(distributionElectronsFile, "%lf", &distributionElectrons[i][j]);
        }
    }

    for(int i = 0; i < Np; ++i){
        double Fp = distributionProtons[i + (Nt - 1) * Np][1] * distributionProtons[i + (Nt - 1) * Np][0] * distributionProtons[i + (Nt - 1) * Np][0];
        if(fabs(Fp) < 1E-100){
            Fp = 1E-100;
        }
        fprintf(gleDistributionProtonsFile, "%20.15g %20.15g\n", distributionProtons[i + (Nt-1)*Np][0], Fp);
        double Fe = distributionElectrons[i + (Nt - 1) * Np][1] * distributionElectrons[i + (Nt - 1) * Np][0] * distributionElectrons[i + (Nt - 1) * Np][0];
        if(fabs(Fe) < 1E-100){
            Fe = 1E-100;
        }
        fprintf(gleDistributionElectronsFile, "%20.15g %20.15g\n", distributionElectrons[i + (Nt - 1) * Np][0], Fe);
    }


    for(int i = 0; i < Np*Nt; ++i){
        delete[] distributionProtons[i];
        delete[] distributionElectrons[i];
    }
    delete[] distributionProtons;
    delete[] distributionElectrons;

    fclose(distributionElectronsFile);
    fclose(distributionProtonsFile);

    fclose(gleDistributionProtonsFile);
    fclose(gleDistributionElectronsFile);



    //fclose(protonsFile);
    //fclose(electronsFile);
    //fclose(positronsFile);
    //fclose(alphasFile);

    //fclose(gleXFile);
    //fclose(gleYFile);
    //fclose(gleZFile);

    //fclose(gleProtonsFile);
    //fclose(gleElectronsFile);
    //fclose(glePositronsFile);
    //fclose(gleAlphasFile);

    return 0;
}

void
outputBfield(int timePointsCount, int Nx, int Ny, int Nz, int Nt, int xpoint, int ypoint, int zpoint, int *timePoints,
             double *middleX, double *middleY, double *middleZ, FILE *BFile, FILE *divergenceFile, FILE *gleBFileX,
             FILE *gleBFileY, FILE *gleBFileZ, FILE *gleDivergenceFileX, FILE *gleDivergenceFileY,
             FILE *gleDivergenceFileZ) {
    double** Bfield = new double*[(Nx - 1) * (Ny - 1) * (Nz - 1) * Nt];
    double** velocity = new double*[(Nx-1)*(Ny-1)*(Nz-1)*Nt];
    double** concentrations = new double*[(Nx-1)*(Ny-1)*(Nz-1)*Nt];
    double** divergence = new double*[(Nx-1)*(Ny-1)*(Nz-1)*Nt];
    for(int i = 0; i < (Nx-1)*(Ny-1)*(Nz-1)*Nt; ++i){
        Bfield[i] =  new double[3];
        velocity[i] =  new double[6];
        concentrations[i] = new double[4];
        divergence[i] = new double[3];
    }

    for(int i = 0; i < (Nx-1)*(Ny-1)*(Nz-1)*Nt; ++i){
        fscanf(BFile, "%lf %lf %lf", &Bfield[i][0], &Bfield[i][1], &Bfield[i][2]);
        //fscanf(velocityProtonFile, "%lf %lf %lf", &velocity[i][0], &velocity[i][1], &velocity[i][2]);
        //fscanf(velocityElectronFile, "%lf %lf %lf", &velocity[i][3], &velocity[i][4], &velocity[i][5]);
        //fscanf(concentrationFile, "%lf %lf %lf %lf", &concentrations[i][0], &concentrations[i][1], &concentrations[i][2], &concentrations[i][3]);
        fscanf(divergenceFile, "%lf %lf %lf", &divergence[i][0], &divergence[i][1], &divergence[i][2]);
    }

    for(int i = 0; i < Nx-1; ++i) {
        fprintf(gleBFileX, "%20.15g", middleX[i]);
        //fprintf(gleVelocityFileX, "%20.15g", Xgrid[i]);
        //fprintf(gleConcentrationFileX, "%20.15g", middleX[i]);
        fprintf(gleDivergenceFileX, "%20.15g", middleX[i]);
        for(int j = 0; j < timePointsCount; ++j){
            int number = timePoints[j]*(Nx-1)*(Ny-1)*(Nz-1) + (Nz-1)*(Ny-1)*i + (Nz-1)*ypoint + zpoint;
            fprintf(gleBFileX, "% 20.15g %20.15g %20.15g", Bfield[number][0], Bfield[number][1], Bfield[number][2]);
            //fprintf(gleVelocityFileX, "% 20.15g %20.15g %20.15g %20.15g %20.15g %20.15g", velocity[number][0], velocity[number][1],velocity[number][2], velocity[number][3], velocity[number][4],velocity[number][5]);
            //fprintf(gleConcentrationFileX, "%20.15g %20.15g %20.15g", concentrations[number][0], concentrations[number][1], concentrations[number][2]);
            fprintf(gleDivergenceFileX, "% 20.15g %20.15g %20.15g", divergence[number][0], divergence[number][1], divergence[number][2]);
        }
        fprintf(gleBFileX, "\n");
        //fprintf(gleVelocityFileX, "\n");
        //fprintf(gleConcentrationFileX, "\r\n");
        fprintf(gleDivergenceFileX, "\n");
    }

    for(int i = 0; i < Ny-1; ++i) {
        fprintf(gleBFileY, "%20.15g ", middleY[i]);
        //fprintf(gleVelocityFileY, "%20.15g ", Ygrid[i]);
        //fprintf(gleConcentrationFileY, "%20.15g ", middleY[i]);
        fprintf(gleDivergenceFileY, "%20.15g ", middleY[i]);
        for(int j = 0; j < timePointsCount; ++j){
            int number = timePoints[j]*(Nx-1)*(Ny-1)*(Nz-1) + (Nz-1)*(Ny-1)*xpoint + (Nz-1)*i + zpoint;
            fprintf(gleBFileY, "%20.15g %20.15g %20.15g", Bfield[number][0], Bfield[number][1], Bfield[number][2]);
            //fprintf(gleVelocityFileY, "%20.15g %20.15g %20.15g %20.15g %20.15g %20.15g", velocity[number][0], velocity[number][1],velocity[number][2], velocity[number][3], velocity[number][4],velocity[number][5]);
            //fprintf(gleConcentrationFileY, "%20.15g %20.15g %20.15g", concentrations[number][0], concentrations[number][1], concentrations[number][2]);
            fprintf(gleDivergenceFileY, "%20.15g %20.15g %20.15g", divergence[number][0], divergence[number][1], divergence[number][2]);
        }
        fprintf(gleBFileY, "\n");
        //fprintf(gleVelocityFileY, "\n");
        //fprintf(gleConcentrationFileY, "\n");
        fprintf(gleDivergenceFileY, "\n");
    }

    for(int i = 0; i < Nz-1; ++i) {
        fprintf(gleBFileZ, "%20.15g ", middleZ[i]);
        //fprintf(gleVelocityFileZ, "%20.15g ", Zgrid[i]);
        //fprintf(gleConcentrationFileZ, "%20.15g ", middleZ[i]);
        fprintf(gleDivergenceFileZ, "%20.15g ", middleZ[i]);
        for(int j = 0; j < timePointsCount; ++j){
            int number = timePoints[j]*(Nx-1)*(Ny-1)*(Nz-1) + (Nz-1)*(Ny-1)*xpoint + (Nz-1)*ypoint + i;
            fprintf(gleBFileZ, "%20.15g %20.15g %20.15g", Bfield[number][0], Bfield[number][1], Bfield[number][2]);
            //fprintf(gleVelocityFileZ, "%20.15g %20.15g %20.15g %20.15g %20.15g %20.15g", velocity[number][0], velocity[number][1],velocity[number][2], velocity[number][3], velocity[number][4],velocity[number][5]);
            //fprintf(gleConcentrationFileZ, "%20.15g %20.15g %20.15g", concentrations[number][0], concentrations[number][1], concentrations[number][2]);
            fprintf(gleDivergenceFileZ, "%20.15g %20.15g %20.15g", divergence[number][0], divergence[number][1], divergence[number][2]);
        }
        fprintf(gleBFileZ, "\n");
        //fprintf(gleVelocityFileZ, "\n");
        //fprintf(gleConcentrationFileZ, "\n");
        fprintf(gleDivergenceFileZ, "\n");
    }

    for(int i = 0; i < (Nx-1)*(Ny-1)*(Nz-1)*Nt; ++i){
        delete[] Bfield[i];
        delete[] velocity[i];
        delete[] concentrations[i];
        delete[] divergence[i];
    }
    delete[] Bfield;
    delete[] velocity;
    delete[] concentrations;
    delete[] divergence;
}

void outputEfield(int timePointsCount, FILE *XFile, FILE *YFile, FILE *ZFile, int Nx, int Ny, int Nz,
             FILE *EFile, FILE *fluxFile, FILE *gleEFileX, FILE *gleEFileY,
             FILE *gleEFileZ, FILE *gleFluxFileX, FILE *gleFluxFileY, FILE *gleFluxFileZ,
             int Nt, int xpoint, int ypoint, int zpoint, int *timePoints, double * Xgrid, double * middleX,
             double * Ygrid, double * middleY, double * Zgrid, double * middleZ) {
    double** Efield = new double*[Nx * Ny * Nz * Nt];
    double** flux = new double*[Nx*Ny*Nz*Nt];
    fscanf(XFile, "%lf", &Xgrid[0]);
    for(int i = 1; i < Nx; ++i){
        double a;
        fscanf(XFile, "%lf", &a);
        Xgrid[i] = a;
        middleX[i-1] = (Xgrid[i-1] + Xgrid[i])/2;
    }

    fscanf(YFile, "%lf", &Ygrid[0]);
    for(int i = 1; i < Ny; ++i){
        double a;
        fscanf(YFile, "%lf", &a);
        Ygrid[i] = a;
        middleY[i-1] = (Ygrid[i-1] + Ygrid[i])/2;
    }

    fscanf(ZFile, "%lf", &Zgrid[0]);
    for(int i = 1; i < Nz; ++i){
        double a;
        fscanf(ZFile, "%lf", &a);
        Zgrid[i] = a;
        middleZ[i-1] = (Zgrid[i-1] + Zgrid[i])/2;
    }

    for(int i = 0; i < Nx*Ny*Nz*Nt; ++i){
        Efield[i] =  new double[3];
        flux[i] = new double[3];
    }

    for(int i = 0; i < Nx*Ny*Nz*Nt; ++i){
        fscanf(EFile, "%lf %lf %lf", &Efield[i][0], &Efield[i][1], &Efield[i][2]);
        fscanf(fluxFile, "%lf %lf %lf", &flux[i][0], &flux[i][1], &flux[i][2]);
    }

    for(int i = 0; i < Nx; ++i) {
        fprintf(gleEFileX, "%20.15g ", Xgrid[i]);
        fprintf(gleFluxFileX, "%20.15g ", Xgrid[i]);
        for(int j = 0; j < timePointsCount; ++j){
            int number = timePoints[j]*Nx*Ny*Nz + Nz*Ny*i + Nz*ypoint + zpoint;
            fprintf(gleEFileX, "%20.15g %20.15g %20.15g", Efield[number][0], Efield[number][1], Efield[number][2]);
            fprintf(gleFluxFileX, "%20.15g %20.15g %20.15g", flux[number][0], flux[number][1], flux[number][2]);
        }
        fprintf(gleEFileX, "\n");
        fprintf(gleFluxFileX, "\n");
    }

    for(int i = 0; i < Ny; ++i) {
        fprintf(gleEFileY, "%20.15g ", Ygrid[i]);
        fprintf(gleFluxFileY, "%20.15g ", Ygrid[i]);
        for(int j = 0; j < timePointsCount; ++j){
            int number = timePoints[j]*Nx*Ny*Nz + Nz*Ny*xpoint + Nz*i + zpoint;
            fprintf(gleEFileY, "%20.15g %20.15g %20.15g", Efield[number][0], Efield[number][1], Efield[number][2]);
            fprintf(gleFluxFileY, "%20.15g %20.15g %20.15g", flux[number][0], flux[number][1], flux[number][2]);
        }
        fprintf(gleEFileY, "\n");
        fprintf(gleFluxFileY, "\n");
    }

    for(int i = 0; i < Nz; ++i) {
        fprintf(gleEFileZ, "%20.15g ", Zgrid[i]);
        fprintf(gleFluxFileZ, "%20.15g ", Zgrid[i]);
        for(int j = 0; j < timePointsCount; ++j){
            int number = timePoints[j]*Nx*Ny*Nz + Nz*Ny*xpoint + Nz*ypoint + i;
            fprintf(gleEFileZ, "%20.15g %20.15g %20.15g", Efield[number][0], Efield[number][1], Efield[number][2]);
            fprintf(gleFluxFileZ, "%20.15g %20.15g %20.15g", flux[number][0], flux[number][1], flux[number][2]);
        }
        fprintf(gleEFileZ, "\n");
        fprintf(gleFluxFileZ, "\n");
    }

    for(int i = 0; i < Nx*Ny*Nz*Nt; ++i){
        delete[] Efield[i];
        delete[] flux[i];
    }
    delete[] Efield;
    delete[] flux;

    fclose(EFile);
    fclose(fluxFile);

    fclose(XFile);
    fclose(YFile);
    fclose(ZFile);

    fclose(gleEFileX);
    fclose(gleEFileY);
    fclose(gleEFileZ);
    fclose(gleFluxFileX);
    fclose(gleFluxFileY);
    fclose(gleFluxFileZ);
}