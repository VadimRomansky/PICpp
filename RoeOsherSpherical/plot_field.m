clear;
field = importdata('./output/full_field_1400.dat');
xfile = importdata('./output/xfile.dat');

Nx = size(xfile,1);
Nt = size(field,1)/Nx;

a = Nt - 1;

figure(1);
plot (xfile(1:Nx),field(a*Nx + (1:Nx), 2), 'red');
title ('B');
xlabel ('x');
ylabel ('B');
grid ;

