clear;
load EfieldXZ.dat;
load BfieldXZ.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
set(0,'DefaultFigureColormap',feval('jet'));
Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx*Nz;
NB = (Nx-1)*(Nz-1);
Nt = fix(size(EfieldXZ, 1)/NE);

a = 0;
b = fix(Nt/2);
c = fix(Nt)-1;

Ex1(1:Nx, 1:Nz) = 0;
Ey1(1:Nx, 1:Nz) = 0;
Ez1(1:Nx, 1:Nz) = 0;

Ex2(1:Nx, 1:Nz) = 0;
Ey2(1:Nx, 1:Nz) = 0;
Ez2(1:Nx, 1:Nz) = 0;

Ex3(1:Nx, 1:Nz) = 0;
Ey3(1:Nx, 1:Nz) = 0;
Ez3(1:Nx, 1:Nz) = 0;

Bx1(1:Nx-1, 1:Nz-1) = 0;
By1(1:Nx-1, 1:Nz-1) = 0;
Bz1(1:Nx-1, 1:Nz-1) = 0;

Bx2(1:Nx-1, 1:Nz-1) = 0;
By2(1:Nx-1, 1:Nz-1) = 0;
Bz2(1:Nx-1, 1:Nz-1) = 0;

Bx3(1:Nx-1, 1:Nz-1) = 0;
By3(1:Nx-1, 1:Nz-1) = 0;
Bz3(1:Nx-1, 1:Nz-1) = 0;

middleX(1:Nx-1) = 0;
middleZ(1:Nz-1) = 0;


for i=1:Nx,
    for j = 1:Nz,
        Ex1(i,j) = EfieldXZ((Nz)*(i-1) + j + a*NE, 1);
        Ey1(i,j) = EfieldXZ((Nz)*(i-1) + j + a*NE, 2);
        Ez1(i,j) = EfieldXZ((Nz)*(i-1) + j + a*NE, 3);
        Ex2(i,j) = EfieldXZ((Nz)*(i-1) + j + b*NE, 1);
        Ey2(i,j) = EfieldXZ((Nz)*(i-1) + j + b*NE, 2);
        Ez2(i,j) = EfieldXZ((Nz)*(i-1) + j + b*NE, 3);
        Ex3(i,j) = EfieldXZ((Nz)*(i-1) + j + c*NE, 1);
        Ey3(i,j) = EfieldXZ((Nz)*(i-1) + j + c*NE, 2);
        Ez3(i,j) = EfieldXZ((Nz)*(i-1) + j + c*NE, 3);
    end;
end;

for i = 1:Nx-1,
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
   for j = 1:Nz-1,
      Bx1(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + a*NB, 1);
      By1(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + a*NB, 2);
      Bz1(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + a*NB, 3);
      Bx2(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + b*NB, 1);
      By2(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + b*NB, 2);
      Bz2(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + b*NB, 3);
      Bx3(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + c*NB, 1);
      By3(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + c*NB, 2);
      Bz3(i, j) = BfieldXZ(((Nz-1)*(i-1) + j) + c*NB, 3);
   end;
end;

for j = 1:Nz-1,
    middleZ(j) = 0.5*(Zfile(j) + Zfile(j+1));
end;

figure(1)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ex1);
shading interp;
title ('Ex1');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(2)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ex2);
shading interp;
title ('Ex2');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(3)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ex3);
shading interp;
title ('Ex3');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(4)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ey1);
shading interp;
title ('Ey1');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(5)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ey2);
shading interp;
title ('Ey2');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(6)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ey3);
shading interp;
title ('Ey3');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(7)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ez1);
shading interp;
title ('Ez1');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(8)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ez2);
shading interp;
title ('Ez2');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(9)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ez3);
shading interp;
title ('Ez3');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(10)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, Bx1);
shading interp;
title ('Bx1');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(11)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, Bx2);
shading interp;
title ('Bx2');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(12)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, Bx3);
shading interp;
title ('Bx3');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(13)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, By1);
shading interp;
title ('By1');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(14)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, By2);
shading interp;
title ('By2');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(15)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, By3);
shading interp;
title ('By3');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(16)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, Bz1);
shading interp;
title ('Bz1');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(17)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, Bz2);
shading interp;
title ('Bz2');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(18)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, Bz3);
shading interp;
title ('Bz3');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;