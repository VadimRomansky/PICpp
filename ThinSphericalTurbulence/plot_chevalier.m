clear;

load chevalier.dat;

N = size(chevalier,1);


set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
title ('F');
xlabel ('t days');
ylabel ('F mJy');

loglog(chevalier(1:N,1), chevalier(1:N,2));

grid ;