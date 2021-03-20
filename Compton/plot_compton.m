clear;

radiation = importdata('output.dat');

Nnu = size(radiation,1);


figure(1);
hold on;
title ('I_{\nu}');
xlabel ('E erg');
ylabel ('I');

plot(radiation(1:Nnu,1),radiation(1:Nnu,2),'red','LineWidth',2);
grid ;
