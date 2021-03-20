clear;

radiation = importdata('output.dat');

Nnu = size(radiation,1);


figure(1);
hold on;
title ('{\nu} F_{\nu}');
xlabel ('E eV');
ylabel ('{\nu} F_{\nu}');

plot(radiation(1:Nnu,1),radiation(1:Nnu,2),'red','LineWidth',2);
grid ;
