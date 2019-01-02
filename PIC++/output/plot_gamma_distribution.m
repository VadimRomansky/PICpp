clear;
distribution_protons = importdata('distribution_protons_5.dat');
distribution_electrons = importdata('distribution_electrons_5.dat');
distribution_alphas = importdata('distribution_alphas_5.dat');
distribution_positrons = importdata('distribution_positrons_5.dat');
load initialParameters.dat;

Np = 500;

Nt = size(distribution_electrons, 1)/Np;
%Nt = 3;


a = 0;
b = fix(Nt/2);
c = Nt - 1;

Fp(1:Np, 1:3) = 0;
Fe(1:Np, 1:3) = 0;
Fa(1:Np, 1:3) = 0;
Fpos(1:Np, 1:3) = 0;

Pp(1:Np, 1:3) = 0;
Pe(1:Np, 1:3) = 0;
Pa(1:Np, 1:3) = 0;
Ppos(1:Np, 1:3) = 0;
Gp(1:Np, 1:3) = 0;
Ge(1:Np, 1:3) = 0;
Ga(1:Np, 1:3) = 0;
Gpos(1:Np, 1:3) = 0;
Pejuttner(1:Np)=0;

me = initialParameters(36);
mp = 1.67262177*10^-24;
ma = 6.64*10^-24;
v=2.998*10^10;
T = 2*10^15;
kB = 1.3806488*10^-16;
theta = kB*T/(me*v*v);
factor = 1;
for i=1:Np,   
   Pp(i,1) = distribution_protons(i + a*Np,1);
   Gp(i,1) = sqrt(Pp(i,1)*Pp(i,1)/(mp*mp*v*v) + 1);
   Pp(i,2) = distribution_protons(i + b*Np,1);
   Gp(i,2) = sqrt(Pp(i,2)*Pp(i,2)/(mp*mp*v*v) + 1);
   Pp(i,3) = distribution_protons(i + c*Np,1);
   Gp(i,3) = sqrt(Pp(i,3)*Pp(i,3)/(mp*mp*v*v) + 1);
   
   factor = mp*v*Gp(i,1)*(Gp(i,1) - 1)/sqrt(Gp(i,1)*Gp(i,1) - 1);
   Fp(i,1) = distribution_protons(i + a*Np, 2)*factor;
   factor = mp*v*Gp(i,2)*(Gp(i,2) - 1)/sqrt(Gp(i,2)*Gp(i,2) - 1);
   Fp(i,2) = distribution_protons(i + b*Np, 2)*factor;
   factor = mp*v*Gp(i,3)*(Gp(i,3) - 1)/sqrt(Gp(i,3)*Gp(i,3) - 1);
   Fp(i,3) = distribution_protons(i + c*Np, 2)*factor;
   
   Pe(i,1) = distribution_electrons(i + a*Np,1);
   Ge(i,1) = sqrt(Pe(i,1)*Pe(i,1)/(me*me*v*v) + 1);
   Pe(i,2) = distribution_electrons(i + b*Np,1);
   Ge(i,2) = sqrt(Pe(i,2)*Pe(i,2)/(me*me*v*v) + 1);
   Pe(i,3) = distribution_electrons(i + c*Np,1);
   Ge(i,3) = sqrt(Pe(i,3)*Pe(i,3)/(me*me*v*v) + 1);
   
   factor = mp*v*Ge(i,1)*(Ge(i,1) - 1)/sqrt(Ge(i,1)*Ge(i,1) - 1);
   Fe(i,1) = distribution_electrons(i + a*Np, 2)*factor;
   factor = mp*v*Ge(i,2)*(Ge(i,2) - 1)/sqrt(Ge(i,2)*Ge(i,2) - 1);
   Fe(i,2) = distribution_electrons(i + b*Np, 2)*factor;
   factor = mp*v*Ge(i,3)*(Ge(i,3) - 1)/sqrt(Ge(i,3)*Ge(i,3) - 1);
   Fe(i,3) = distribution_electrons(i + c*Np, 2)*factor;
   
   Pa(i,1) = distribution_alphas(i + a*Np,1);
   Ga(i,1) = sqrt(Pa(i,1)*Pa(i,1)/(ma*ma*v*v) + 1);
   Pa(i,2) = distribution_alphas(i + b*Np,1);
   Ga(i,2) = sqrt(Pa(i,2)*Pa(i,2)/(ma*ma*v*v) + 1);
   Pa(i,3) = distribution_alphas(i + c*Np,1);
   Ga(i,3) = sqrt(Pa(i,3)*Pa(i,3)/(ma*ma*v*v) + 1);
   
   factor = mp*v*Ga(i,1)*(Ga(i,1) - 1)/sqrt(Ga(i,1)*Ga(i,1) - 1);
   Fa(i,1) = distribution_alphas(i + a*Np, 2)*factor;
   factor = mp*v*Ga(i,2)*(Ga(i,2) - 1)/sqrt(Ga(i,2)*Ga(i,2) - 1);
   Fa(i,2) = distribution_alphas(i + b*Np, 2)*factor;
   factor = mp*v*Ga(i,3)*(Ga(i,3) - 1)/sqrt(Ga(i,3)*Ga(i,3) - 1);
   Fa(i,3) = distribution_alphas(i + c*Np, 2)*factor;
   
   Ppos(i,1) = distribution_positrons(i + a*Np,1);
   Gpos(i,1) = sqrt(Ppos(i,1)*Ppos(i,1)/(me*me*v*v) + 1);
   Ppos(i,2) = distribution_positrons(i + b*Np,1);
   Gpos(i,2) = sqrt(Ppos(i,2)*Ppos(i,2)/(me*me*v*v) + 1);
   Ppos(i,3) = distribution_positrons(i + c*Np,1);
   Gpos(i,3) = sqrt(Ppos(i,3)*Ppos(i,3)/(me*me*v*v) + 1);
   
   factor = mp*v*Gpos(i,1)*(Gpos(i,1) - 1)/sqrt(Gpos(i,1)*Gpos(i,1) - 1);
   Fpos(i,1) = distribution_positrons(i + a*Np, 2)*factor;
   factor = mp*v*Gpos(i,2)*(Gpos(i,2) - 1)/sqrt(Gpos(i,2)*Gpos(i,2) - 1);
   Fpos(i,2) = distribution_positrons(i + b*Np, 2)*factor;
   factor = mp*v*Gpos(i,3)*(Gpos(i,3) - 1)/sqrt(Gpos(i,3)*Gpos(i,3) - 1);
   Fpos(i,3) = distribution_positrons(i + c*Np, 2)*factor;
   
   exp1 = exp(-sqrt(1+Pe(i,1)*Pe(i,1)/(me*me*v*v))/theta);
   bes = besselk(2, 1/theta);
   p = Pe(i,1);
   p3 = (p/(me*v))^3;
   Pejuttner(i) = (1.0/(theta*bes))*exp1*p3*Pe(i,1);
end;
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 2);
figure(1);
plot (Gp(1:Np,1)-1,Fp(1:Np,1), 'red',Gp(1:Np,2)-1,Fp(1:Np,2), 'green',Gp(1:Np,3)-1,Fp(1:Np,3), 'blue');
title ('protons distribution function');
xlabel ('{\gamma}-1');
ylabel ('F_p({\gamma}) ({\gamma}-1)');
legend('t=0','t=T/2','t=T','Location','southeast');
grid ;

figure(2);
plot (Ge(1:Np,1)-1,Fe(1:Np,1), 'red',Ge(1:Np,2)-1,Fe(1:Np,2), 'green',Ge(1:Np,3)-1,Fe(1:Np,3), 'blue');
%plot (Pe(1:Np,1)/(me*v),Fe(1:Np,1), 'red',Pe(1:Np,2)/(me*v),Fe(1:Np,2), 'green',Pe(1:Np,3)/(me*v),Fe(1:Np,3), 'blue', Pe(1:Np,1)/(me*v), Pejuttner(1:Np), 'black');
title ('electrons distribution function');
xlabel ('{\gamma}-1');
ylabel ('F_e({\gamma}) ({\gamma}-1)');
legend('t=0','t=T/2','t=T','Location','southeast');
grid ;

figure(3);
plot (Ga(1:Np,1)-1,Fa(1:Np,1), 'red',Ga(1:Np,2)-1,Fa(1:Np,2), 'green',Ga(1:Np,3)-1,Fa(1:Np,3), 'blue');
title ('alphas distribution function');
xlabel ('{\gamma}-1');
ylabel ('F_alpha({\gamma}) ({\gamma}-1)');
legend('t=0','t=T/2','t=T','Location','southeast');
grid ;

figure(4);
plot (Gpos(1:Np,1)-1,Fpos(1:Np,1), 'red',Gpos(1:Np,2)-1,Fpos(1:Np,2), 'green',Gpos(1:Np,3)-1,Fpos(1:Np,3), 'blue');
title ('positrons distribution function');
xlabel ('{\gamma}-1');
ylabel ('F_{e+}({\gamma}) ({\gamma}-1)');
legend('t=0','t=T/2','t=T','Location','southeast');
grid ;