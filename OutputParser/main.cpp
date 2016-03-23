#include <stdio.h>
#include <math.h>

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
    FILE* XFile = fopen("../../PIC++/output/Xfile.dat","r");
    FILE* YFile = fopen("../../PIC++/output/Yfile.dat","r");
    FILE* ZFile = fopen("../../PIC++/output/Zfile.dat","r");
    FILE* BFile = fopen("../../PIC++/output/Bfield.dat","r");
    FILE* EFile = fopen("../../PIC++/output/Efield.dat", "r");
    FILE* concentrationFile = fopen("../../PIC++/output/concentrations.dat", "r");
    FILE* velocityElectronFile = fopen("../../PIC++/output/velocity_electron.dat", "r");
    FILE* velocityProtonFile = fopen("../../PIC++/output/velocity.dat", "r");
    FILE* traectoryProtonFile = fopen("../../PIC++/output/traectory_proton.dat", "r");
    FILE* traectoryElectronFile = fopen("../../PIC++/output/traectory_electron.dat", "r");
    FILE* protonsFile = fopen("../../PIC++/output/protons.dat", "r");
    FILE* electronsFile = fopen("../../PIC++/output/electrons.dat", "r");
    FILE* positronsFile = fopen("../../PIC++/output/positrons.dat", "r");
    FILE* alphasFile = fopen("../../PIC++/output/alphas.dat", "r");
    FILE* divergenceFile = fopen("../../PIC++/output/divergence_error.dat", "r");
    FILE* distributionProtonsFile = fopen("../../PIC++/output/distribution_protons.dat", "r");
    FILE* distributionElectronsFile = fopen("../../PIC++/output/divergence_electrons.dat", "r");

    FILE* gleXFile = fopen("../output/Xfile.dat","w");
    FILE* gleYFile = fopen("../output/Yfile.dat","w");
    FILE* gleZFile = fopen("../output/Zfile.dat","w");
    FILE* gleBFile = fopen("..output/Bfield.dat","w");
    FILE* gleEFile = fopen("../output/Efield.dat", "w");
    FILE* gleConcentrationFile = fopen("../output/concentrations.dat", "w");
    FILE* gleVelocityFile = fopen("../output/velocity.dat", "w");
    FILE* gleTraectoryProtonFile = fopen("../output/traectory_proton.dat", "w");
    FILE* gleTraectoryElectronFile = fopen("../output/traectory_electron.dat", "w");
    FILE* gleProtonsFile = fopen("../output/protons.dat", "w");
    FILE* gleElectronsFile = fopen("../output/electrons.dat", "w");
    FILE* glePositronsFile = fopen("../output/positrons.dat", "w");
    FILE* gleAlphasFile = fopen("../output/alphas.dat", "w");
    FILE* gleDivergenceFile = fopen("../output/divergence_error.dat", "w");
    FILE* gleDistributionProtonsFile = fopen("../output/distribution_protons.dat", "w");
    FILE* gleDistributionElectronsFile = fopen("../output/divergence_electrons.dat", "w");

    int Nx = countNumbers(XFile);
    int Ny = countNumbers(YFile);
    int Nz = countNumbers(ZFile);

    int NE = countNumbersLine(EFile, 3);
    int Nt = NE/(Nx*Ny*Nz);

    int ypoint = 0;
    int zpoint = 0;

    int timePoints[timePointsCount];
    timePoints[0] = 0;
    timePoints[timePointsCount-1] = Nt - 1;
    int deltaNt = Nt/(timePointsCount - 1);
    for(int i = 1; i < timePointsCount - 1; ++i){
        timePoints[i] = timePoints[i-1] + deltaNt;
    }

    double** Efield = new double*[Nx*Ny*Nz*Nt];
    double* Xgrid = new double[Nx];
    double* middleX = new double[Nx - 1];
    fscanf(XFile, "%lf", Xgrid[0]);
    for(int i = 1; i < Nx; ++i){
        fscanf(XFile, "%lf", Xgrid[i]);
        middleX[i-1] = (Xgrid[i-1] + Xgrid[i])/2;
    }

    for(int i = 0; i < Nx*Ny*Nz*Nt; ++i){
        Efield[i] =  new double[3];
    }

    for(int i = 0; i < Nx*Ny*Nz*Nt; ++i){
        fscanf(EFile, "%lf %lf %lf", Efield[i][0], Efield[i][1], Efield[i][2]);
    }

    for(int i = 0; i < Nx; ++i) {
        for(int j = 0; j < timePointsCount; ++j){
            int number = timePoints[j]*Nx*Ny*Nz + Nz*Ny*i + Nz*ypoint + zpoint;
            fprintf(gleEFile, "%20.15g %20.15g %20.15g", Efield[number][0], Efield[number][1], Efield[number][2]);
        }
        fprintf(gleEFile, "\n");
    }

    for(int i = 0; i < Nx*Ny*Nz*Nt; ++i){
        delete[] Efield[i];
    }
    delete[] Efield;
    delete[] Xgrid;

    double** Bfield = new double*[(Nx-1)*(Ny-1)*(Nz-1)*Nt];
    double** velocity = new double*[(Nx-1)*(Ny-1)*(Nz-1)*Nt];
    double** concentrations = new double*[(Nx-1)*(Ny-1)*(Nz-1)*Nt];
    for(int i = 0; i < (Nx-1)*(Ny-1)*(Nz-1)*Nt; ++i){
        Bfield[i] =  new double[3];
        velocity[i] =  new double[6];
        concentrations[i] = new double[4];
    }

    for(int i = 0; i < (Nx-1)*(Ny-1)*(Nz-1)*Nt; ++i){
        fscanf(BFile, "%lf %lf %lf", Bfield[i][0], Bfield[i][1], Bfield[i][2]);
        fscanf(velocityProtonFile, "%lf %lf %lf", velocity[i][0], velocity[i][1], velocity[i][2]);
        fscanf(velocityElectronFile, "%lf %lf %lf", velocity[i][3], velocity[i][4], velocity[i][5]);
        fscanf(concentrationFile, "%lf %lf %lf %lf", concentrations[i][0], concentrations[i][1], concentrations[i][2], concentrations[i][3]);
    }

    for(int i = 0; i < Nx-1; ++i) {
        for(int j = 0; j < timePointsCount; ++j){
            int number = timePoints[j]*(Nx-1)*(Ny-1)*(Nz-1) + (Nz-1)*(Ny-1)*i + (Nz-1)*ypoint + zpoint;
            fprintf(gleBFile, "%20.15g %20.15g %20.15g", Bfield[number][0], Bfield[number][1], Bfield[number][2]);
            fprintf(gleVelocityFile, "%20.15g %20.15g %20.15g %20.15g %20.15g %20.15g", velocity[number][0], velocity[number][1],velocity[number][2], velocity[number][3], velocity[number][4],velocity[number][5]);
        }
        fprintf(gleBFile, "\n");
        fprintf(gleVelocityFile, "\n");
    }

    fclose(XFile);
    fclose(YFile);
    fclose(ZFile);
    fclose(BFile);
    fclose(EFile);
    fclose(concentrationFile);
    fclose(velocityElectronFile);
    fclose(velocityProtonFile);
    fclose(traectoryElectronFile);
    fclose(traectoryProtonFile);
    fclose(protonsFile);
    fclose(electronsFile);
    fclose(positronsFile);
    fclose(alphasFile);
    fclose(divergenceFile);
    fclose(distributionProtonsFile);
    fclose(distributionElectronsFile);

    fclose(gleXFile);
    fclose(gleYFile);
    fclose(gleZFile);
    fclose(gleBFile);
    fclose(gleEFile);
    fclose(gleConcentrationFile);
    fclose(gleVelocityFile);
    fclose(gleTraectoryElectronFile);
    fclose(gleTraectoryProtonFile);
    fclose(gleProtonsFile);
    fclose(gleElectronsFile);
    fclose(glePositronsFile);
    fclose(gleAlphasFile);
    fclose(gleDivergenceFile);
    fclose(gleDistributionProtonsFile);
    fclose(gleDistributionElectronsFile);
    return 0;
}