clear;
EfieldXZ = importdata('EfieldXZ_20.dat');
BfieldXZ = importdata('BfieldXZ_20.dat');
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
set(0,'DefaultFigureColormap',feval('jet'));
set(0, 'DefaultLineLineWidth', 2);
Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = (Nx-1)*(Nz-1);
NB = (Nx-1)*(Nz-1);
Nt = fix(size(EfieldXZ, 1)/NE);

a = 0;
b = fix(Nt/2);
c = fix(Nt)-1;
Ex(1:Nx-1, 1:Nz-1) = 0;
Ey(1:Nx-1, 1:Nz-1) = 0;
Ez(1:Nx-1, 1:Nz-1) = 0;

Bx(1:Nx-1, 1:Nz-1) = 0;
By(1:Nx-1, 1:Nz-1) = 0;
Bz(1:Nx-1, 1:Nz-1) = 0;

middleX(1:Nx-1) = 0;
middleZ(1:Nz-1) = 0;
Xgrid(1:Nx-1) = 0;
Zgrid(1:Nz-1) = 0;

for i=1:Nx-1,
    Xgrid(i) = (Xfile(i) - Xfile(2));
    for j = 1:Nz-1,
        Ex(i,j) = EfieldXZ((Nz-1)*(i-1) + j + c*NE, 1);
        Ey(i,j) = EfieldXZ((Nz-1)*(i-1) + j + c*NE, 2);
        Ez(i,j) = EfieldXZ((Nz-1)*(i-1) + j + c*NE, 3);
    end;
end;

for i = 1:Nx-1,
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
   for j = 1:Nz-1,
      Bx(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + c*NB, 1);
      By(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + c*NB, 2);
      Bz(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + c*NB, 3);
   end;
end;

for j = 1:Nz-1,
   Zgrid(j) = (Zfile(j) - Zfile(2));
end;

for j = 1:Nz-1,
    middleZ(j) = 0.5*(Zfile(j) + Zfile(j+1));
end;

figure(1)
[X, Z] = meshgrid(Zgrid, Xgrid);
surf(X, Z, Ex);
shading interp;
title ('Ex');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(2)
surf(X, Z, Ey);
shading interp;
title ('Ey');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(3)
surf(X, Z, Ez);
shading interp;
title ('Ez');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(4)
surf(X, Z, Bx);
shading interp;
title ('Bx');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(5)
surf(X, Z, By);
shading interp;
title ('By');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(6)
surf(X, Z, Bz);
shading interp;
title ('Bz');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;