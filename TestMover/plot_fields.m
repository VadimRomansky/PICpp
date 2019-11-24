clear;
Bx =  importdata('outBx.dat');
By =  importdata('outBy.dat');
Bz =  importdata('outBz.dat');
Ex =  importdata('outEx.dat');
Ey =  importdata('outEy.dat');
Ez =  importdata('outEz.dat');

Nx = size(Bx,1);
Ny = size(Bx,2);

figure(1);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nx));
surf(X, Y, Bx);
shading interp;
title ('Bx');
xlabel ('y');
ylabel ('x');
zlabel ('Bx');
grid ;

figure(2);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nx));
surf(X, Y, By);
shading interp;
title ('By');
xlabel ('y');
ylabel ('x');
zlabel ('By');
grid ;

figure(3);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nx));
surf(X, Y, Bz);
shading interp;
title ('Bz');
xlabel ('y');
ylabel ('x');
zlabel ('Bz');
grid ;

figure(4);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nx));
surf(X, Y, Ex);
shading interp;
title ('Ex');
xlabel ('y');
ylabel ('x');
zlabel ('Ex');
grid ;

figure(5);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nx));
surf(X, Y, Ey);
shading interp;
title ('Ey');
xlabel ('y');
ylabel ('x');
zlabel ('Ey');
grid ;

figure(6);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nx));
surf(X, Y, Ez);
shading interp;
title ('Ez');
xlabel ('y');
ylabel ('x');
zlabel ('Ez');
grid ;

