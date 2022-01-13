clear;

load Bpoints.dat;
load Npoints.dat;
err = importdata('error.dat');

Nn = size(Npoints,1);
Nb = size(Bpoints,1);

figure(1);
colormap Jet;
[X, Y] = meshgrid(Npoints, Bpoints);
surf(X, Y, err);
shading interp;
title ('dI^2');
xlabel ('n');
ylabel ('B');
zlabel ('Bx');
grid ;
