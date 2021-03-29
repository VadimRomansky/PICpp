clear;

photons = importdata('photons.dat');

Nnu = size(photons,1);


figure(1);
hold on;
title ('N_{\nu}');
xlabel ('E eV');
ylabel ('N_{\nu} cm^{-3} Hz^{-1}');

plot(photons(1:Nnu,1),photons(1:Nnu,2),'red','LineWidth',2);
grid ;
