clear;
field = importdata('./output/full_field_120.dat');
xfile = importdata('./output/xfile.dat');
kfile = importdata('./output/kfile.dat');

Nx = size(xfile,1);
Nk = size(kfile,1);


figure(1);
plot (kfile(1:Nk),field(Nx/10, 1:Nk), 'red');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('B');
xlabel ('k');
ylabel ('B');
grid ;


figure(2);
plot (kfile(1:Nk),field(2*Nx/10, 1:Nk), 'red');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('B');
xlabel ('k');
ylabel ('B');
grid ;

figure(3);
plot (kfile(1:Nk),field(4*Nx/10, 1:Nk), 'red');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('B');
xlabel ('k');
ylabel ('B');
grid ;

figure(4);
plot (kfile(1:Nk),field(5*Nx/10, 1:Nk), 'red');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('B');
xlabel ('k');
ylabel ('B');
grid ;

