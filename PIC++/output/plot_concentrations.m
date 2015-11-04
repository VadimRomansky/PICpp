clear;
load concentrations.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = Nx*Ny*Nz;

Nt = size(concentrations, 1)/N;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

electron_concentration(1:Nx, 1:3) = 0;
proton_concentration(1:Nx, 1:3) = 0;
charge_density(1:Nx, 1:3) = 0;

ynumber = 1;
znumber = 1;


for i=1:Nx,   
   electron_concentration(i, 1) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + a*N, 1);
   electron_concentration(i, 2) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + b*N, 1);
   electron_concentration(i, 3) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + c*N, 1);
   proton_concentration(i, 1) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + a*N, 2);
   proton_concentration(i, 2) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + b*N, 2);
   proton_concentration(i, 3) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + c*N, 2);
   charge_density(i, 1) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + a*N, 3);
   charge_density(i, 2) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + b*N, 3);
   charge_density(i, 3) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + c*N, 3);
end;

figure(1);
plot (Xfile(1:Nx,1),electron_concentration(1:Nx, 1), 'red', Xfile(1:Nx,1), electron_concentration(1:Nx, 2), 'green', Xfile(1:Nx,1), electron_concentration(1:Nx, 3), 'blue');
title ('electrons');
xlabel ('x/r_g');
ylabel ('n cm^-3');
grid ;

figure(2);
plot (Xfile(1:Nx,1),proton_concentration(1:Nx, 1), 'red', Xfile(1:Nx,1), proton_concentration(1:Nx, 2), 'green', Xfile(1:Nx,1), proton_concentration(1:Nx, 3), 'blue');
title ('protons');
xlabel ('x/r_g');
ylabel ('n cm^-3');
grid ;

figure(3);
plot (Xfile(1:Nx,1),charge_density(1:Nx, 1), 'red', Xfile(1:Nx,1), charge_density(1:Nx, 2), 'green', Xfile(1:Nx,1), charge_density(1:Nx, 3), 'blue');
title ('charge density');
xlabel ('x/r_g');
ylabel ('n sgs*cm^-3');
grid ;