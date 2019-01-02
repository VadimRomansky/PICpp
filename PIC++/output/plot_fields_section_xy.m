clear;
EfieldXY = importdata('EfieldXY_5.dat');
BfieldXY = importdata('BfieldXY_5.dat');
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
load initialParameters.dat;
set(0,'DefaultFigureColormap',feval('jet'));
Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx*Ny;
NB = (Nx-1)*(Ny-1);
Nt = fix(size(EfieldXY, 1)/NE);
NtB = (size(BfieldXY, 1)/NB);

a = 0;
b = fix(Nt/2);
c = fix(Nt)-1;
%c = 0;

N1 = 1;
N2 = fix(Nx)-1;

Ex(1:N2-N1+1, 1:Ny) = 0;
Ey(1:N2-N1+1, 1:Ny) = 0;
Ez(1:N2-N1+1, 1:Ny) = 0;

Bx(1:N2-N1, 1:Ny-1) = 0;
By(1:N2-N1, 1:Ny-1) = 0;
Bz(1:N2-N1, 1:Ny-1) = 0;

middleX(1:N2-N1) = 0;
middleY(1:Ny-1) = 0;
Xgrid(1:N2-N1+1) = 0;
Ygrid(1:Ny) = 0;
cv = initialParameters(10);
omega = initialParameters(21);
omegaElectron = initialParameters(20);



for i=1:N2-N1+1,
    %Xgrid(i) = (Xfile(i) - Xfile(2))*omegaElectron/cv;
   Xgrid(i) = (Xfile(i) - Xfile(2));
    for j = 1:Ny,
        Ex(i,j) = EfieldXY((Ny)*(i+N1-1) + j + c*NE, 1);
        Ey(i,j) = EfieldXY((Ny)*(i+N1-1) + j + c*NE, 2);
        Ez(i,j) = EfieldXY((Ny)*(i+N1-1) + j + c*NE, 3);
    end;
end;

for j = 1:Ny,
    %Ygrid(j) = (Yfile(j) - Yfile(2))*omegaElectron/cv;
   Ygrid(j) = (Yfile(j) - Yfile(2));
end;

for i = 1:N2-N1,
   %middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2))*omegaElectron/cv;
   middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2));
   for j = 1:Ny-1,
      Bx(i, j) = BfieldXY(((Ny-1)*(i+N1-1) + j) + c*NB, 1);
      By(i, j) = BfieldXY(((Ny-1)*(i+N1-1) + j) + c*NB, 2);
      Bz(i, j) = BfieldXY(((Ny-1)*(i+N1-1) + j) + c*NB, 3);
   end;
end;

for j = 1:Ny-1,
    %middleY(j) = (0.5*(Yfile(j) + Yfile(j+1)) - Yfile(2))*omegaElectron/cv;
   middleY(j) = (0.5*(Yfile(j) + Yfile(j+1)) - Yfile(2));
end;

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
%title ('B_z G');
%xlabel ('y\omega_e/c');
%ylabel ('x\omega_e/c');

xlabel ('y');
ylabel ('x');
zlabel ('B_z G');
