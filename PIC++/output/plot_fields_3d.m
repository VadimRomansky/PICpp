clear;
Efield = importdata('Efield.dat');
Bfield = inportdata('Bfield.dat');
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
load initialParameters.dat;
set(0,'DefaultFigureColormap',feval('jet'));
Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx*Ny*Nz;
NB = (Nx-1)*(Ny-1)*(Nz-1);
Nt = fix(size(Efield, 1)/NE);
NtB = (size(Bfield, 1)/NB);
znumber = 2;

a = 0;
b = fix(Nt/2);
c = fix(Nt)-1;
%c = 0;

Ex(1:Nx, 1:Ny) = 0;
Ey(1:Nx, 1:Ny) = 0;
Ez(1:Nx, 1:Ny) = 0;

Bx(1:Nx-1, 1:Ny-1) = 0;
By(1:Nx-1, 1:Ny-1) = 0;
Bz(1:Nx-1, 1:Ny-1) = 0;

middleX(1:Nx-1) = 0;
middleY(1:Ny-1) = 0;
Xgrid(1:Nx) = 0;
Ygrid(1:Ny) = 0;
cv = initialParameters(10);
omega = initialParameters(21);
omegaElectron = initialParameters(20);



for i=1:Nx,
    Xgrid(i) = (Xfile(i) - Xfile(2))*omegaElectron/cv;
   %Xgrid(i) = (Xfile(i) - Xfile(2));
    for j = 1:Ny,
        Ex(i,j) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(j-1) + znumber + c*NE, 1);
        Ey(i,j) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(j-1) + znumber + c*NE, 2);
        Ez(i,j) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(j-1) + znumber + c*NE, 3);
    end;
end;

for j = 1:Ny,
    Ygrid(j) = (Yfile(j) - Yfile(2))*omegaElectron/cv;
   %Ygrid(i) = (Yfile(i) - Yfile(2));
end;

for i = 1:Nx-1,
   middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2))*omegaElectron/cv;
   %middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2));
   for j = 1:Ny-1,
      Bx(i, j) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(j-1) + znumber) + c*NB, 1);
      By(i, j) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(j-1) + znumber) + c*NB, 2);
      Bz(i, j) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(j-1) + znumber) + c*NB, 3);
   end;
end;

for j = 1:Ny-1,
    middleY(j) = (0.5*(Yfile(j) + Yfile(j+1)) - Yfile(2))*omegaElectron/cv;
   %middleY(i) = (0.5*(Yfile(i) + Yfile(i+1)) - Yfile(2));
end;

set(0, 'DefaultLineLineWidth', 2);

figure(1)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ex);
shading interp;
title ('Ex');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(2)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ey);
shading interp;
title ('Ey');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(3)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ez);
shading interp;
title ('Ez');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(4)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, Bx);
shading interp;
title ('Bx');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');

figure(5)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, By);
shading interp;
title ('By');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');


figure(6)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, Bz);
shading interp;
title ('Bz');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');
