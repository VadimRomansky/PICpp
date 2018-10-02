clear;
load distribution_protons_grid.dat;
load distribution_electrons_grid.dat;
load distribution_alphas_grid.dat;
load distribution_positrons_grid.dat;
load initialParameters.dat;
load Xfile.dat;

Np = 500;

Nx = size(Xfile,1) - 1;
Nt = size(distribution_electrons_grid_20, 1)/(Nx + 1);
%Nt = 3;


a = Nt - 1;

Fp(1:Np) = 0;
Fe(1:Np) = 0;
Fa(1:Np) = 0;
Fpos(1:Np) = 0;

Pp(1:Np) = 0;
Pe(1:Np) = 0;
Pa(1:Np) = 0;
Ppos(1:Np) = 0;
Pejuttner(1:Np)=0;

me = initialParameters(36);
mp = 1.67262177*10^-24;
ma = 6.64*10^-24;
v=2.998*10^10;
T = 2*10^15;
kB = 1.3806488*10^-16;
theta = kB*T/(me*v*v);
factor = 1;
minX = 1;
maxX = Nx;
for i=1:Np,   
   Pp(i) = distribution_protons_grid(1 + a*(Nx + 1),i)/(mp*v);
   
   Pe(i) = distribution_electrons_grid(1 + a*(Nx + 1),i)/(me*v);
     
   Pa(i) = distribution_alphas_grid(1 + a*(Nx + 1),i)/(ma*v);
     
   Ppos(i) = distribution_positrons_grid(1 + a*(Nx + 1),i)/(me*v);
   
   for j=minX:maxX,
     factor = Pp(i)*Pp(i);
     Fp(i) = Fp(i) + distribution_protons_grid(1 + j + a*(Nx + 1), i)*factor;

     factor = Pe(i)*Pe(i);
     Fe(i) = Fe(i) + distribution_electrons_grid(1 + j + a*(Nx + 1), i)*factor;

     factor = Pa(i)*Pa(i);
     Fa(i) = Fa(i) + distribution_alphas_grid(1 + j + a*(Nx + 1), i)*factor;

     factor = Ppos(i)*Ppos(i);
     Fpos(i) = Fpos(i) + distribution_positrons_grid(1 + j + a*(Nx + 1), i)*factor;
   end;
   exp1 = exp(-sqrt(1+Pe(i)*Pe(i)/(me*me*v*v))/theta);
   bes = besselk(2, 1/theta);
   p = Pe(i);
   p3 = (p/(me*v))^3;
   Pejuttner(i) = (1.0/(theta*bes))*exp1*p3*Pe(i);
end;
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 2);
figure(1);
plot (Pp(1:Np),Fp(1:Np), 'red');
title ('protons distribution function');
xlabel ('p/mc');
ylabel ('F_p(p) p^4');
%legend('t=0','t=T/2','t=T','Location','southeast');
grid ;

figure(2);
plot (Pe(1:Np),Fe(1:Np), 'red');
%plot (Pe(1:Np,1)/(me*v),Fe(1:Np,1), 'red',Pe(1:Np,2)/(me*v),Fe(1:Np,2), 'green',Pe(1:Np,3)/(me*v),Fe(1:Np,3), 'blue', Pe(1:Np,1)/(me*v), Pejuttner(1:Np), 'black');
title ('electrons distribution function');
xlabel ('p/mc');
ylabel ('F_e(p) p^4');
%legend('t=0','t=T/2','t=T','Location','southeast');
grid ;

figure(3);
plot (Pa(1:Np),Fa(1:Np), 'red');
title ('alphas distribution function');
xlabel ('p/mc');
ylabel ('F_alpha(p) p^4');
%legend('t=0','t=T/2','t=T','Location','southeast');
grid ;

figure(4);
plot (Ppos(1:Np),Fpos(1:Np), 'red');
title ('positrons distribution function');
xlabel ('p/mc');
ylabel ('F_e+(p) p^4');
%legend('t=0','t=T/2','t=T','Location','southeast');
grid ;