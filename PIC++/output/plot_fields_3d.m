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
Nt = size(Efield, 1)/NE;
znumber = 1;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

Ex(1:Nx, 1:Ny) = 0;
Ey(1:Nx, 1:Ny) = 0;
Ez(1:Nx, 1:Ny) = 0;

Bx(1:Nx-1, 1:Ny-1) = 0;
By(1:Nx-1, 1:Ny-1) = 0;
Bz(1:Nx-1, 1:Ny-1) = 0;

middleX(1:Nx-1) = 0;
middleY(1:Ny-1) = 0;


for i=1:Nx,
    for j = 1:Ny,
        Ex(i,j) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(j-1) + znumber + c*NE, 1);
        Ey(i,j) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(j-1) + znumber + c*NE, 2);
        Ez(i,j) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(j-1) + znumber + c*NE, 3);
    end;
end;

for i = 1:Nx-1,
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
   for j = 1:Ny-1,
      Bx(i, j) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(j-1) + znumber) + c*NB, 1);
      By(i, j) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(j-1) + znumber) + c*NB, 2);
      Bz(i, j) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(j-1) + znumber) + c*NB, 3);
   end;
end;

for j = 1:Ny-1,
    middleY(j) = 0.5*(Yfile(j) + Yfile(j+1));
end;

figure(1)
[X, Y] = meshgrid(Yfile, Xfile);
surf(X, Y, Ex);
title ('Ex');
xlabel ('x');
ylabel ('y');
zlabel ('E gauss');
grid ;

figure(2)
[X, Y] = meshgrid(Yfile, Xfile);
surf(X, Y, Ey);
title ('Ey');
xlabel ('x');
ylabel ('y');
zlabel ('E gauss');
grid ;

figure(3)
[X, Y] = meshgrid(Yfile, Xfile);
surf(X, Y, Ez);
title ('Ez');
xlabel ('x');
ylabel ('y');
zlabel ('E gauss');
grid ;

figure(4)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, Bx);
title ('Bx');
xlabel ('x');
ylabel ('y');
zlabel ('B gauss');
grid ;

figure(5)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, By);
title ('By');
xlabel ('x');
ylabel ('y');
zlabel ('B gauss');
grid ;

figure(6)
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, Bz);
title ('Bz');
xlabel ('x');
ylabel ('y');
zlabel ('B gauss');
grid ;