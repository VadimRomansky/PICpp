clear;
%load concentrations.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
load divergence_error.dat;
load initialParameters.dat;
set(0,'DefaultFigureColormap',feval('jet'));
Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = (Nx)*(Ny)*(Nz);
Nt = fix(size(divergence_error, 1)/N)-1;

znumber = 2;

a = 0;
b = fix(Nt/2);
c = Nt;

divergenceError(1:Nx, 1:Ny) = 0;
divergence(1:Nx, 1:Ny) = 0;
density(1:Nx, 1:Ny) = 0;
divergenceB(1:Nx, 1:Ny) = 0;

middleX(1:Nx)=0;
middleY(1:Ny)=0;

cv = initialParameters(10);
omega = initialParameters(21);
omegaElectron = initialParameters(20);

for i=1:Nx,
   middleX(i) = ((Xfile(i)+Xfile(i+1))/2 - Xfile(2))*omegaElectron/cv;
   for j =1:Ny,
      divergenceError(i, j) = divergence_error((Nz*Ny*(i-1) + Nz*(j-1) + znumber) + c*N, 1);
      divergence(i, j) = divergence_error((Nz*Ny*(i-1) + Nz*(j-1) + znumber) + c*N, 2);
      density(i,j) = divergence_error((Nz*Ny*(i-1) + Nz*(j-1) + znumber) + c*N, 3);
      divergenceB(i, j) = divergence_error((Nz*Ny*(i-1) + Nz*(j-1) + znumber) + c*N, 4);
   end;
end;

for j = 1:Ny,
    middleY(j) = (0.5*(Yfile(j) + Yfile(j+1)) - Yfile(2))*omegaElectron/cv;
end;

set(0, 'DefaultLineLineWidth', 2);

figure(1)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, divergenceError);
shading interp;
title ('divergence error');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('divergence error e cgs/cm^3');

figure(2)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, divergence);
shading interp;
title ('divergence');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('divergence e cgs/cm^3');

figure(3)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, density);
shading interp;
title ('4*pi*density');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('4*pi*density e cgs/cm^3');

figure(4);
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, divergenceB);
shading interp;
title ('div B');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('divergence e cgs/cm^3');