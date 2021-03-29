clear;

radiation1 = importdata('output1.dat');
radiation2 = importdata('output2.dat');

Nnu = size(radiation1,1);


figure(1);
hold on;
title ('{\nu} F_{\nu}');
xlabel ('E eV');
ylabel ('{\nu} F_{\nu}');

plot(radiation1(1:Nnu,1),radiation1(1:Nnu,2),'red',radiation2(1:Nnu,1),radiation2(1:Nnu,2),'blue','LineWidth',2);
grid ;
