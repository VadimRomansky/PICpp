clear;
load distribution_electrons.dat;
load distribution_protons.dat;
load distribution_alphas.dat;

Np = size(distribution_electrons, 1);

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
plot (distribution_electrons(1:Np,1),distribution_electrons(1:Np,2), 'blue');
xlabel ('p g*cm/s');
ylabel ('F_e(p)*p^2');
grid ;

figure(2);
plot (distribution_protons(1:Np,1),distribution_protons(1:Np,2), 'blue');
xlabel ('p g*cm/s');
ylabel ('F_p(p)*p^2');
grid ;

figure(3);
plot (distribution_alphas(1:Np,1),distribution_alphas(1:Np,2), 'blue');
xlabel ('p g*cm/s');
ylabel ('F_{\alpha}(p)*p^2');
grid ;

