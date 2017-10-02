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
xnumber = 2;

a = 0;
b = fix(Nt/2);
c = fix(Nt)-2;

Ex(1:Ny, 1:Nz) = 0;
Ey(1:Ny, 1:Nz) = 0;
Ez(1:Ny, 1:Nz) = 0;

Bx(1:Ny-1, 1:Nz-1) = 0;
By(1:Ny-1, 1:Nz-1) = 0;
Bz(1:Ny-1, 1:Nz-1) = 0;

middleY(1:Ny-1) = 0;
middleZ(1:Nz-1) = 0;


for i=1:Ny,
    for j = 1:Nz,
        Ex(i,j) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + j + c*NE, 1);
        Ey(i,j) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + j + c*NE, 2);
        Ez(i,j) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + j + c*NE, 3);
    end;
end;

for i = 1:Ny-1,
   middleY(i) = 0.5*(Yfile(i) + Yfile(i+1));
   for j = 1:Nz-1,
      Bx(i, j) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + j) + c*NB, 1);
      By(i, j) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + j) + c*NB, 2);
      Bz(i, j) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + j) + c*NB, 3);
   end;
end;

for j = 1:Nz-1,
    middleZ(j) = 0.5*(Zfile(j) + Zfile(j+1));
end;

figure(1)
[Y, Z] = meshgrid(Zfile, Yfile);
surf(Y, Z, Ex);
shading interp;
title ('Ex');
xlabel ('z');
ylabel ('y');
zlabel ('E gauss');
grid ;

figure(2)
[Y, Z] = meshgrid(Zfile, Yfile);
surf(Y, Z, Ey);
shading interp;
title ('Ey');
xlabel ('z');
ylabel ('y');
zlabel ('E gauss');
grid ;

figure(3)
[Y, Z] = meshgrid(Zfile, Yfile);
surf(Y, Z, Ez);
shading interp;
title ('Ez');
xlabel ('z');
ylabel ('y');
zlabel ('E gauss');
grid ;

figure(4)
[Y, Z] = meshgrid(middleZ, middleY);
surf(Y, Z, Bx);
shading interp;
title ('Bx');
xlabel ('z');
ylabel ('x');
zlabel ('B gauss');
grid ;

figure(5)
[Y, Z] = meshgrid(middleZ, middleY);
surf(Y, Z, By);
shading interp;
title ('By');
xlabel ('z');
ylabel ('y');
zlabel ('B gauss');
grid ;

figure(6)
[Y, Z] = meshgrid(middleZ, middleY);
surf(Y, Z, Bz);
shading interp;
title ('Bz');
xlabel ('z');
ylabel ('y');
zlabel ('B gauss');
grid ;