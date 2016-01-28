clear;
load fluxFile.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx*Ny*Nz;
Nt = fix(size(fluxFile,1)/NE);

znumber = 1;

a = 0;
b = fix(Nt/2);
c = Nt-1;

Jx(1:Nx, 1:Ny) = 0;
Jy(1:Nx, 1:Ny) = 0;
Jz(1:Nx, 1:Ny) = 0;


for i=1:Nx,
    for j = 1:Ny,
        Jx(i,j) = fluxFile((Nz)*(Ny)*(i-1) + (Nz)*(j-1) + znumber + c*NE, 1);
        Jy(i,j) = fluxFile((Nz)*(Ny)*(i-1) + (Nz)*(j-1) + znumber + c*NE, 2);
        Jz(i,j) = fluxFile((Nz)*(Ny)*(i-1) + (Nz)*(j-1) + znumber + c*NE, 3);
    end;
end;

figure(1)
[X, Y] = meshgrid(Yfile, Xfile);
surf(X, Y, Jx);
title ('Jx');
xlabel ('y');
ylabel ('x');
zlabel ('flux');
grid ;

figure(2)
[X, Y] = meshgrid(Yfile, Xfile);
surf(X, Y, Jy);
title ('Jy');
xlabel ('y');
ylabel ('x');
zlabel ('flux');
grid ;

figure(3)
[X, Y] = meshgrid(Yfile, Xfile);
surf(X, Y, Jz);
title ('Jz');
xlabel ('y');
ylabel ('x');
zlabel ('flux');
grid ;