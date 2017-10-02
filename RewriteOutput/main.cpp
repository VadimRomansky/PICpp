#include <stdio.h>
#include <math.h>
#include <string>

#include "constants.h"

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
	const int typesNumber = 8;

    //for windowss
    std::string outputDir = std::string(outputDirectory);
	std::string rewriteDir = "C:/users/Vadim/Documents/Visual Studio 2010/Projects/v4/trunk/RewriteOutput/output/";
	const char* timeSuffix = "_50s";
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

    FILE* gleEFileX = fopen((rewriteDir + "Efield.dat").c_str(), "w");


    int NE = countNumbersLine(EFile, 3);
    rewind(XFile);
    rewind(YFile);
    rewind(ZFile);
    rewind(EFile);
    int Nt = NE/(Nx*Ny*Nz);
	Nt = 24;
	int timePoint = Nt - 1;

    int ypoint = 0;
    int zpoint = 0;

    int timePoints[timePointsCount];
    timePoints[0] = 0;
    timePoints[timePointsCount-1] = Nt - 1;
    int deltaNt = Nt/(timePointsCount - 1);
    for(int i = 1; i < timePointsCount - 1; ++i){
        timePoints[i] = timePoints[i-1] + deltaNt;
    }

    /////// arrays with dimensions Nx*Ny*Nz*Nt //////
    double** Efield = new double*[Nx*Ny*Nz*Nt];
    double* Xgrid = new double[Nx];
    double* middleX = new double[Nx - 1];
    double* Ygrid = new double[Ny];
    double* middleY = new double[Ny - 1];
    double* Zgrid = new double[Nz];
    double* middleZ = new double[Nz - 1];
    fscanf(XFile, "%lf", &Xgrid[0]);
    for(int i = 1; i < Nx; ++i){
        double a;
        fscanf(XFile, "%lf", &a);
        Xgrid[i] = a;
        middleX[i-1] = (Xgrid[i-1] + Xgrid[i])/2;
    }
	double shiftX = Xgrid[1];
	for(int i = 0; i < Nx; ++i) {
		Xgrid[i] -= shiftX;
	}
	for(int i = 0; i < Nx-1; ++i) {
		middleX[i] -= shiftX;
	}

    for(int i = 0; i < Nx*Ny*Nz*Nt; ++i){
        Efield[i] =  new double[3];
    }

    for(int i = 0; i < Nx*Ny*Nz*Nt; ++i){
        fscanf(EFile, "%lf %lf %lf", &Efield[i][0], &Efield[i][1], &Efield[i][2]);
    }

	//todo
    for(int i = 1; i < Nx-1; ++i) {
        fprintf(gleEFileX, "%20.15g ", Xgrid[i]);
        int number = timePoint*Nx*Ny*Nz + Nz*Ny*i + Nz*ypoint + zpoint;
        fprintf(gleEFileX, "%20.15g %20.15g %20.15g", Efield[number][0], Efield[number][1], Efield[number][2]);
        
        fprintf(gleEFileX, "\n");
    }

    for(int i = 0; i < Nx*Ny*Nz*Nt; ++i){
        delete[] Efield[i];
    }
    delete[] Efield;

    fclose(EFile);
    fclose(fluxFile);

    fclose(XFile);
    fclose(YFile);
    fclose(ZFile);

    fclose(gleEFileX);


    FILE* BFile = fopen((outputDir + "Bfield.dat").c_str(),"r");
    FILE* concentrationFile = fopen((outputDir + "concentrations.dat").c_str(), "r");
    //FILE* velocityElectronFile = fopen((inputDir + "velocity_electron.dat").c_str(), "r");
    //FILE* velocityProtonFile = fopen((inputDir + "velocity.dat").c_str(), "r");

    FILE* gleBFileX = fopen((rewriteDir + "Bfield.dat").c_str(),"w");
    FILE* gleConcentrationElectronsFile = fopen((rewriteDir + "concentration_electrons.dat").c_str(), "w");
    FILE* gleConcentrationProtonsFile = fopen((rewriteDir + "concentration_protons.dat").c_str(), "w");
    FILE* gleConcentrationAlphasFile = fopen((rewriteDir + "concentration_alphas.dat").c_str(), "w");
    //FILE* gleConcentrationFileY = fopen("../output/concentrations_y.dat", "w");
    //FILE* gleConcentrationFileZ = fopen("../output/concentrations_z.dat", "w");
    //FILE* gleVelocityFileX = fopen("../output/velocity_x.dat", "w");
    //FILE* gleVelocityFileY = fopen("../output/velocity_y.dat", "w");
    //FILE* gleVelocityFileZ = fopen("../output/velocity_z.dat", "w");


    /////// arrays with dimensions (Nx - 1)*(Ny - 1)(Nz-1)*Nt //////
    double** Bfield = new double*[(Nx-1)*(Ny-1)*(Nz-1)*Nt];
    double** velocity = new double*[(Nx-1)*(Ny-1)*(Nz-1)*Nt];
    double** concentrations = new double*[(Nx-1)*(Ny-1)*(Nz-1)*Nt];
    double** divergence = new double*[(Nx-1)*(Ny-1)*(Nz-1)*Nt];
    for(int i = 0; i < (Nx-1)*(Ny-1)*(Nz-1)*Nt; ++i){
        Bfield[i] =  new double[3];
        velocity[i] =  new double[6];
        concentrations[i] = new double[typesNumber + 2];
        divergence[i] = new double[3];
    }

    for(int i = 0; i < (Nx-1)*(Ny-1)*(Nz-1)*Nt; ++i){
        fscanf(BFile, "%lf %lf %lf", &Bfield[i][0], &Bfield[i][1], &Bfield[i][2]);
        //fscanf(velocityProtonFile, "%lf %lf %lf", &velocity[i][0], &velocity[i][1], &velocity[i][2]);
        //fscanf(velocityElectronFile, "%lf %lf %lf", &velocity[i][3], &velocity[i][4], &velocity[i][5]);
		for(int j = 0; j < typesNumber + 2; ++j){
			fscanf(concentrationFile, "%lf", &concentrations[i][j]);
		}
    }

    for(int i = 1; i < Nx-2; ++i) {
        fprintf(gleBFileX, "%20.15g", middleX[i]);
        //fprintf(gleVelocityFileX, "%20.15g", Xgrid[i]);
        fprintf(gleConcentrationElectronsFile, "%20.15g", middleX[i]);
        fprintf(gleConcentrationProtonsFile, "%20.15g", middleX[i]);
        fprintf(gleConcentrationAlphasFile, "%20.15g", middleX[i]);
            int number = timePoint*(Nx-1)*(Ny-1)*(Nz-1) + (Nz-1)*(Ny-1)*i + (Nz-1)*ypoint + zpoint;
            fprintf(gleBFileX, "% 20.15g %20.15g %20.15g", Bfield[number][0], Bfield[number][1], Bfield[number][2]);
            //fprintf(gleVelocityFileX, "% 20.15g %20.15g %20.15g %20.15g %20.15g %20.15g", velocity[number][0], velocity[number][1],velocity[number][2], velocity[number][3], velocity[number][4],velocity[number][5]);
            fprintf(gleConcentrationElectronsFile, "%20.15g %20.15g", concentrations[number][2]);
            fprintf(gleConcentrationProtonsFile, "%20.15g %20.15g", concentrations[number][3]);
            fprintf(gleConcentrationAlphasFile, "%20.15g %20.15g", concentrations[number][5]);
        
        fprintf(gleBFileX, "\n");
        //fprintf(gleVelocityFileX, "\n");
        fprintf(gleConcentrationElectronsFile, "\n");
        fprintf(gleConcentrationProtonsFile, "\n");
        fprintf(gleConcentrationAlphasFile, "\n");
    }

    delete[] Xgrid;
    delete[] middleX;
    delete[] Ygrid;
    delete[] middleY;
    delete[] Zgrid;
    delete[] middleZ;


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

    fclose(BFile);
    fclose(concentrationFile);
    //fclose(velocityElectronFile);
    //fclose(velocityProtonFile);

    fclose(gleBFileX);
    fclose(gleConcentrationElectronsFile);
    fclose(gleConcentrationProtonsFile);
    fclose(gleConcentrationAlphasFile);
    //fclose(gleConcentrationFileY);
    //fclose(gleConcentrationFileZ);
    //fclose(gleVelocityFileX);
    //fclose(gleVelocityFileY);
    //fclose(gleVelocityFileZ);

    /////// arrays with dimensions Nt //////

    FILE* generalFile = fopen((outputDir + "general.dat").c_str(), "r");
    FILE* incrementFile = fopen((outputDir + "increment.dat").c_str(), "r");
    //FILE* trajectoryProtonFile = fopen((inputDir + "trajectory_proton.dat").c_str(), "r");
    //FILE* trajectoryElectronFile = fopen((inputDir + "trajectory_electron.dat").c_str(), "r");

    FILE* gleGeneralFile = fopen((rewriteDir + "general.dat").c_str(), "w");
    //FILE* gleTrajectoryProtonFile = fopen("../output/trajectory_proton.dat", "w");
    //FILE* gleTrajectoryElectronFile = fopen("../output/trajectory_electron.dat", "w");

    double** general = new double*[Nt];
    double** protonTraectory = new double*[Nt];
    double** electronTraectory = new double*[Nt];
    for(int i = 0; i < Nt; ++i){
        general[i] = new double[19];
        protonTraectory[i] = new double[8];
        electronTraectory[i] = new double[8];
    }

    for(int i = 0; i < Nt; i++){
        for(int j = 0; j < 18; ++j){
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
    /////// arrays with dimensions pnumber*Nt //////

    FILE* distributionProtonsFile = fopen((outputDir + "distribution_protons.dat").c_str(), "r");
    FILE* distributionElectronsFile = fopen((outputDir + "distribution_electrons.dat").c_str(), "r");
    FILE* distributionAlphasFile = fopen((outputDir + "distribution_alphas.dat").c_str(), "r");

    FILE* gleDistributionProtonsFile = fopen((rewriteDir + "distribution_protons.dat").c_str(), "w");
    FILE* gleDistributionElectronsFile = fopen((rewriteDir + "distribution_electrons.dat").c_str(), "w");
    FILE* gleDistributionAlphasFile = fopen((rewriteDir + "distribution_alphas.dat").c_str(), "w");

    double** distributionProtons = new double*[pnumber*Nt];
    double** distributionElectrons = new double*[pnumber*Nt];
    double** distributionAlphas = new double*[pnumber*Nt];
    for(int i = 0; i < pnumber*Nt; ++i){
        distributionProtons[i] = new double[2];
        distributionElectrons[i] = new double[2];
        distributionAlphas[i] = new double[2];
    }

    for(int i = 0; i < pnumber*Nt; ++i){
        for(int j = 0; j < 2; ++j) {
            fscanf(distributionProtonsFile, "%lf", &distributionProtons[i][j]);
            fscanf(distributionElectronsFile, "%lf", &distributionElectrons[i][j]);
            fscanf(distributionAlphasFile, "%lf", &distributionAlphas[i][j]);
        }
    }

    for(int i = 0; i < pnumber; ++i){
        double Fp = distributionProtons[i + (Nt - 1) * pnumber][1];    
        fprintf(gleDistributionProtonsFile, "%20.15g %20.15g\n", distributionProtons[i + (Nt-1)*pnumber][0], Fp);

        double Fe = distributionElectrons[i + (Nt - 1) * pnumber][1];
        fprintf(gleDistributionElectronsFile, "%20.15g %20.15g\n", distributionElectrons[i + (Nt - 1) * pnumber][0], Fe);

        double Fa = distributionAlphas[i + (Nt - 1) * pnumber][1];
        fprintf(gleDistributionAlphasFile, "%20.15g %20.15g\n", distributionAlphas[i + (Nt - 1) * pnumber][0], Fa);
    }


    for(int i = 0; i < pnumber*Nt; ++i){
        delete[] distributionProtons[i];
        delete[] distributionElectrons[i];
        delete[] distributionAlphas[i];
    }
    delete[] distributionProtons;
    delete[] distributionElectrons;
    delete[] distributionAlphas;

    fclose(distributionElectronsFile);
    fclose(distributionProtonsFile);
    fclose(distributionAlphasFile);

    fclose(gleDistributionProtonsFile);
    fclose(gleDistributionElectronsFile);
    fclose(gleDistributionAlphasFile);



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