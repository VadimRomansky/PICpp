clear;

load B.dat;

Nx = 20;
Ny = 20;
Nz = 20;

zpoint = 5;

Bx(1:Nx, 1:Ny) = 0;
By(1:Nx, 1:Ny) = 0;

for i = 1:Nx,
    for j = 1:Ny,
        Bx(i,j) = B((i-1)*Ny*Nz + (j - 1)*Nz + zpoint,1);
        By(i,j) = B((i-1)*Ny*Nz + (j - 1)*Nz + zpoint,2);
    end;
end;


set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
title ('B');
xlabel ('x');
ylabel ('y');

quiver(1:Nx, 1:Ny, By, Bx);

grid ;