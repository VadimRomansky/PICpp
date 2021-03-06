clear;
EfieldXZ = importdata('EfieldXZ_5.dat');
BfieldXZ = importdata('BfieldXZ_5.dat');
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

Ex(1:Nx, 1:Nz) = 0;
Ey(1:Nx, 1:Nz) = 0;
Ez(1:Nx, 1:Nz) = 0;

Bx(1:Nx-1, 1:Nz-1) = 0;
By(1:Nx-1, 1:Nz-1) = 0;
Bz(1:Nx-1, 1:Nz-1) = 0;

middleX(1:Nx-1) = 0;
middleZ(1:Nz-1) = 0;


for i=1:Nx,
    for j = 1:Nz,
        Ex(i,j) = EfieldXZ((Nz)*(i-1) + j + c*NE, 1);
        Ey(i,j) = EfieldXZ((Nz)*(i-1) + j + c*NE, 2);
        Ez(i,j) = EfieldXZ((Nz)*(i-1) + j + c*NE, 3);
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
    middleZ(j) = 0.5*(Zfile(j) + Zfile(j+1));
end;

figure(1)
%[X, Z] = meshgrid(Zfile, Xfile);
%surf(X, Z, Ex);
quiver (Zfile, Xfile, Ez, Ex);
shading interp;
title ('Exz');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(2)
%[X, Z] = meshgrid(middleZ, middleX);
%surf(X, Z, Bx);
quiver (middleZ, middleX, Bz, Bx);
shading interp;
title ('Bxz');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;