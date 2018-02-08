clear;
load EfieldXY.dat;
load BfieldXY.dat;
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

Ex1(1:Nx, 1:Ny) = 0;
Ey1(1:Nx, 1:Ny) = 0;
Ez1(1:Nx, 1:Ny) = 0;

Ex2(1:Nx, 1:Ny) = 0;
Ey2(1:Nx, 1:Ny) = 0;
Ez2(1:Nx, 1:Ny) = 0;

Ex3(1:Nx, 1:Ny) = 0;
Ey3(1:Nx, 1:Ny) = 0;
Ez3(1:Nx, 1:Ny) = 0;

Bx1(1:Nx-1, 1:Ny-1) = 0;
By1(1:Nx-1, 1:Ny-1) = 0;
Bz1(1:Nx-1, 1:Ny-1) = 0;

Bx2(1:Nx-1, 1:Ny-1) = 0;
By2(1:Nx-1, 1:Ny-1) = 0;
Bz2(1:Nx-1, 1:Ny-1) = 0;

Bx3(1:Nx-1, 1:Ny-1) = 0;
By3(1:Nx-1, 1:Ny-1) = 0;
Bz3(1:Nx-1, 1:Ny-1) = 0;

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
        Ex1(i,j) = EfieldXY((Ny)*(i-1) + j + a*NE, 1);
        Ey1(i,j) = EfieldXY((Ny)*(i-1) + j + a*NE, 2);
        Ez1(i,j) = EfieldXY((Ny)*(i-1) + j + a*NE, 3);
        Ex2(i,j) = EfieldXY((Ny)*(i-1) + j + b*NE, 1);
        Ey2(i,j) = EfieldXY((Ny)*(i-1) + j + b*NE, 2);
        Ez2(i,j) = EfieldXY((Ny)*(i-1) + j + b*NE, 3);
        Ex3(i,j) = EfieldXY((Ny)*(i-1) + j + c*NE, 1);
        Ey3(i,j) = EfieldXY((Ny)*(i-1) + j + c*NE, 2);
        Ez3(i,j) = EfieldXY((Ny)*(i-1) + j + c*NE, 3);
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
      Bx1(i, j) = BfieldXY(((Ny-1)*(i-1) + j) + a*NB, 1);
      By1(i, j) = BfieldXY(((Ny-1)*(i-1) + j) + a*NB, 2);
      Bz1(i, j) = BfieldXY(((Ny-1)*(i-1) + j) + a*NB, 3);
      Bx2(i, j) = BfieldXY(((Ny-1)*(i-1) + j) + b*NB, 1);
      By2(i, j) = BfieldXY(((Ny-1)*(i-1) + j) + b*NB, 2);
      Bz2(i, j) = BfieldXY(((Ny-1)*(i-1) + j) + b*NB, 3);
      Bx3(i, j) = BfieldXY(((Ny-1)*(i-1) + j) + c*NB, 1);
      By3(i, j) = BfieldXY(((Ny-1)*(i-1) + j) + c*NB, 2);
      Bz3(i, j) = BfieldXY(((Ny-1)*(i-1) + j) + c*NB, 3);
   end;
end;

for j = 1:Ny-1,
    middleY(j) = (0.5*(Yfile(j) + Yfile(j+1)) - Yfile(2))*omegaElectron/cv;
   %middleY(i) = (0.5*(Yfile(i) + Yfile(i+1)) - Yfile(2));
end;

figure(1)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ex1);
shading interp;
title ('Ex1');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(2)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ex2);
shading interp;
title ('Ex2');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(3)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ex3);
shading interp;
title ('Ex3');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(4)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ey1);
shading interp;
title ('Ey1');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(5)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ey2);
shading interp;
title ('Ey2');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(6)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ey3);
shading interp;
title ('Ey3');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(7)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ez1);
shading interp;
title ('Ez1');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(8)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ez2);
shading interp;
title ('Ez2');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(9)
[X, Y] = meshgrid(Ygrid, Xgrid);
surf(X, Y, Ez3);
shading interp;
title ('Ez3');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('E gauss');

figure(10)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, Bx1);
shading interp;
title ('Bx1');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');

figure(11)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, Bx2);
shading interp;
title ('Bx2');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');

figure(12)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, Bx3);
shading interp;
title ('Bx3');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');

figure(13)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, By1);
shading interp;
title ('By1');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');

figure(14)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, By1);
shading interp;
title ('By1');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');

figure(15)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, By1);
shading interp;
title ('By1');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');


figure(16)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, Bz1);
shading interp;
title ('Bz1');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');

figure(17)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, Bz2);
shading interp;
title ('Bz2');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');

figure(18)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, Bz3);
shading interp;
title ('Bz3');
xlabel ('y\omega_e/c');
ylabel ('x\omega_e/c');
zlabel ('B gauss');
