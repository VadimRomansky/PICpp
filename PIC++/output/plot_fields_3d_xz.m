clear;
load Efield.dat;
load Bfield.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx*Ny*Nz;
NB = (Nx-1)*(Ny-1)*(Nz-1);
Nt = fix(size(Efield, 1)/NE);
ynumber = 1;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

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
        Ex(i,j) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + j + c*NE, 1);
        Ey(i,j) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + j + c*NE, 2);
        Ez(i,j) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + j + c*NE, 3);
    end;
end;

for i = 1:Nx-1,
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
   for j = 1:Nz-1,
      Bx(i, j) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + j) + c*NB, 1);
      By(i, j) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + j) + c*NB, 2);
      Bz(i, j) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + j) + c*NB, 3);
   end;
end;

for j = 1:Nz-1,
    middleZ(j) = 0.5*(Zfile(j) + Zfile(j+1));
end;

figure(1)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ex);
title ('Ex');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(2)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ey);
title ('Ey');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(3)
[X, Z] = meshgrid(Zfile, Xfile);
surf(X, Z, Ez);
title ('Ez');
xlabel ('z');
ylabel ('x');
zlabel ('E gauss');
grid ;

figure(4)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, Bx);
title ('Bx');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(5)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, By);
title ('By');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(6)
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, Bz);
title ('Bz');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;