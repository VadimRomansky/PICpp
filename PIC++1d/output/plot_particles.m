clear;
load protons.dat;
load electrons.dat

N1=size(protons,1);
N2=size(electrons,1);


figure(1);
plot (protons(1:N1,1), protons(1:N1,2), 'red', electrons(1:N2,1), electrons(1:N2,2), 'blue');
title ('particles');
xlabel ('x cm');
ylabel ('P g*cm/s');
legend('protons', 'electrons','Location','northeast');
grid ;

figure(2);
plot (protons(1:N1,1), protons(1:N1,3), 'red', electrons(1:N2,1), electrons(1:N2,3), 'blue');
title ('particles');
xlabel ('x cm');
ylabel ('P g*cm/s');
legend('protons', 'electrons','Location','northeast');
grid ;