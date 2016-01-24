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

electron_concentration(1:Nx, 1:Ny) = 0;
proton_concentration(1:Nx, 1:Ny) = 0;
charge_density(1:Nx, 1:Ny) = 0;

znumber = 1;

middleX(1:Nx) = 0;
middleY(1:Ny) = 0;

for i=1:Nx,  
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
   for j = 1:Ny,
       electron_concentration(i, j) = concentrations(Nz*Ny*(i-1) + Nz*(j-1) + znumber + c*N, 1);
       proton_concentration(i, j) = concentrations(Nz*Ny*(i-1) + Nz*(j-1) + znumber + c*N, 2);
       charge_density(i, j) = concentrations(Nz*Ny*(i-1) + Nz*(j-1) + znumber + c*N, 3);
   end;
end;

for j = 1:Ny,
    middleY(j) = 0.5*(Yfile(j) + Yfile(j+1));
end;

figure(1);
[X, Y] = meshgrid(middleY, middleX);
mesh(X, Y, electron_concentration);
title ('electrons');
xlabel ('y');
ylabel ('x');
zlabel ('n cm^-3');
grid ;

figure(2);
[X, Y] = meshgrid(middleY, middleX);
mesh(X, Y, proton_concentration);
title ('protons');
xlabel ('y');
ylabel ('x');
zlabel ('n cm^-3');
grid ;

figure(3);
[X, Y] = meshgrid(middleY, middleX);
mesh(X, Y, charge_density);
title ('charge density');
xlabel ('y');
ylabel ('x');
zlabel ('n cm^-3');
grid ;