clear;
distribuion = importdata('./output/fullDistribution_19.dat');
pfile = importdata('./output/pfile.dat');

Np = size(pfile,1);
Nt = size(distribuion,1)/Np;

a = Nt - 1;

figure(1);
plot (pfile(1:Np),distribuion(a*Np + (1:Np), 2), 'red');
title ('f');
xlabel ('p');
ylabel ('f');
grid ;
