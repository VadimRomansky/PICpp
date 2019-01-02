clear;
load particleTypes.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
load initialParameters.dat;

set(0, 'DefaultLineLineWidth', 2);

concentrations = importdata('concentrationsY_5.dat');

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = Ny;

Nt = size(concentrationsY, 1)/N;
%Nt = 23;
Ntypes = size(particleTypes, 1);

a = 0;
b = fix(Nt/2);
c = fix(Nt)-1;

charge_density(1:Ny, 1:3) = 0;
charge_density_hat(1:Ny, 1:3) = 0;
particle_concentrations(1:Ny, 1:3*Ntypes) = 0;

middleY(1:Ny) = 0;
cv = initialParameters(10);
omega = initialParameters(21);
omegaElectron = initialParameters(20);


for i=1:Ny,
   %middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2))*omegaElectron/cv;
   middleY(i) = (0.5*(Yfile(i) + Yfile(i+1)) - Yfile(2));
   for t = 1:Ntypes,
        particle_concentrations(i, 1 + 3*(t-1)) = concentrationsY(i + a*N, 2 + t);
        particle_concentrations(i, 2 + 3*(t-1)) = concentrationsY(i + b*N, 2 + t);
        particle_concentrations(i, 3 + 3*(t-1)) = concentrationsY(i + c*N, 2 + t);
  
   end;
   charge_density(i, 1) = concentrationsY(i + a*N, 1);
   charge_density(i, 2) = concentrationsY(i + b*N, 1);
   charge_density(i, 3) = concentrationsY(i + c*N, 1);
   
   charge_density_hat(i, 1) = concentrationsY(i + a*N, 2);
   charge_density_hat(i, 2) = concentrationsY(i + b*N, 2);
   charge_density_hat(i, 3) = concentrationsY(i + c*N, 2);
end;

for t = 1:Ntypes,
    if(particleTypes(t) > 0)
    figure(t);
    plot (middleY(1:Ny),particle_concentrations(1:Ny, 1 + 3*(t-1)), 'red', middleY(1:Ny), particle_concentrations(1:Ny, 2 + 3*(t-1)), 'green', middleY(1:Ny), particle_concentrations(1:Ny, 3 + 3*(t-1)), 'blue');
    title ('concentration');
    xlabel ('y\omega_p/c');
    ylabel ('n cm^-3');
    grid ;    
    end;
end;

figure(Ntypes + 1);
plot (middleY(1:Ny),charge_density(1:Ny, 1), 'red', middleY(1:Ny), charge_density(1:Ny, 2), 'green', middleY(1:Ny), charge_density(1:Ny, 3), 'blue');
title ('charge density');
%xlabel ('x\omega_p/c');
xlabel ('y');
ylabel ('rho cgs*cm^-3');
grid ;

figure(Ntypes + 2);
plot (middleY(1:Ny),charge_density_hat(1:Ny, 1), 'red', middleY(1:Ny), charge_density_hat(1:Ny, 2), 'green', middleY(1:Ny), charge_density_hat(1:Ny, 3), 'blue');
title ('charge density hat');
%xlabel ('x\omega_p/c');
xlabel ('y');
ylabel ('rho cgs*cm^-3');
grid ;