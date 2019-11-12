clear;
distribuion = importdata('./output/distribution_9.dat');
pfile = importdata('./output/pfile.dat');

Np = size(pfile,1);
Nt = size(distribuion,1)/Np;

a = Nt - 1;

figure(1);
plot (pfile(1:Np),distribuion(a*Np + (1:Np), 2), 'red');
title ('shock');
xlabel ('p');
ylabel ('f');
grid ;

figure(2);
plot (pfile(1:Np),distribuion(a*Np + (1:Np), 4), 'red');
title ('upstream');
xlabel ('p');
ylabel ('f');
grid ;

