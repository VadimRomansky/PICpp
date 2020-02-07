clear;

load 'Bx.dat'
load 'By.dat'
load 'Bz.dat'
load 'Theta.dat'

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
surf(X, Y, Theta);
shading interp;
title ('Theta');
xlabel ('y');
ylabel ('x');
zlabel ('Theta');
grid ;