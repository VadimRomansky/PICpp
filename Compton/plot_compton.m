clear;


radiationE = importdata('outputE.dat');


radiation = importdata('outputNu.dat');


Nnu = size(radiation,1);


figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\nu}');
xlabel ('E eV');
ylabel ('F_{\nu} mJy');

plot(radiationE(1:Nnu,1),radiation(1:Nnu,2),'red','LineWidth',2);
%plot(radiation1(1:Nnu,1),Fa(1:Nnu),'blue','LineWidth',2);
%legend('Dubus','Uvarov/mc^2','\nu^{-1.5}');
grid ;

