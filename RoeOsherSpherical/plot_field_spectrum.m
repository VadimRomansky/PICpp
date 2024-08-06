clear;
field = importdata('./output/full_field_1400.dat');
xfile = importdata('./output/xfile.dat');
kfile = importdata('./output/kfile.dat');

Nx = size(xfile,1);
Nk = size(kfile,1);


figure(1);
plot (1:Nk,field(Nx/4, 1:Nk), 'red');
title ('B');
xlabel ('k');
ylabel ('B');
grid ;

