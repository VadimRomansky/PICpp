clear;
load concentration_electrons.dat;
load concentration_protons.dat;
load concentration_alphas.dat;

Np = size(concentration_electrons, 1);

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
plot (concentration_electrons(1:Np,1),concentration_electrons(1:Np,2), 'blue');
xlabel ('x cm');
ylabel ('n_e cm^{-3}');
grid ;

figure(2);
plot (concentration_protons(1:Np,1),concentration_protons(1:Np,2), 'blue');
xlabel ('x cm');
ylabel ('n_p cm^{-3}');
grid ;

figure(3);
plot (concentration_alphas(1:Np,1),concentration_alphas(1:Np,2), 'blue');
xlabel ('x cm');
ylabel ('n_{\alpha} cm^{-3}');
grid ;

