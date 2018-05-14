clear;
load distribution_protons.dat;
load distribution_electrons.dat;
load distribution_alphas.dat;
load distribution_positrons.dat;
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

me = initialParameters(36);
mp = 1.67262177*10^-24;
v=3*10^10;

for i=1:Np,
   Pp(i,1) = distribution_protons(i + a*Np,1);
   Pp(i,2) = distribution_protons(i + b*Np,1);
   Pp(i,3) = distribution_protons(i + c*Np,1);
   
   Fp(i,1) = distribution_protons(i + a*Np, 2)*Pp(i,1)*Pp(i,1);
   Fp(i,2) = distribution_protons(i + b*Np, 2)*Pp(i,2)*Pp(i,2);
   Fp(i,3) = distribution_protons(i + c*Np, 2)*Pp(i,3)*Pp(i,3);
   
   Pe(i,1) = distribution_electrons(i + a*Np,1);
   Pe(i,2) = distribution_electrons(i + b*Np,1);
   Pe(i,3) = distribution_electrons(i + c*Np,1);
   
   Fe(i,1) = distribution_electrons(i + a*Np, 2)*Pe(i,1)*Pe(i,1);
   Fe(i,2) = distribution_electrons(i + b*Np, 2)*Pe(i,2)*Pe(i,2);
   Fe(i,3) = distribution_electrons(i + c*Np, 2)*Pe(i,3)*Pe(i,3);
   
   Pa(i,1) = distribution_alphas(i + a*Np,1);
   Pa(i,2) = distribution_alphas(i + b*Np,1);
   Pa(i,3) = distribution_alphas(i + c*Np,1);
   
   Fa(i,1) = distribution_alphas(i + a*Np, 2)*Pa(i,1)*Pa(i,1);
   Fa(i,2) = distribution_alphas(i + b*Np, 2)*Pa(i,2)*Pa(i,2);
   Fa(i,3) = distribution_alphas(i + c*Np, 2)*Pa(i,3)*Pa(i,3);
   
   Ppos(i,1) = distribution_positrons(i + a*Np,1);
   Ppos(i,2) = distribution_positrons(i + b*Np,1);
   Ppos(i,3) = distribution_positrons(i + c*Np,1);
   
   Fpos(i,1) = distribution_positrons(i + a*Np, 2)*Ppos(i,1)*Ppos(i,1);
   Fpos(i,2) = distribution_positrons(i + b*Np, 2)*Ppos(i,2)*Ppos(i,2);
   Fpos(i,3) = distribution_positrons(i + c*Np, 2)*Ppos(i,3)*Ppos(i,3);
end;
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
plot (Pp(1:Np,1)/(mp*v),Fp(1:Np,1), 'red',Pp(1:Np,2)/(mp*v),Fp(1:Np,2), 'green',Pp(1:Np,3)/(mp*v),Fp(1:Np,3), 'blue');
%title ('protons distribution function');
xlabel ('p/{m_e c}');
ylabel ('F_p(p) p^4');
legend('t=0','t=T/2','t=T','Location','southeast');
grid ;

figure(2);
plot (Pe(1:Np,1)/(me*v),Fe(1:Np,1), 'red',Pe(1:Np,2)/(me*v),Fe(1:Np,2), 'green',Pe(1:Np,3)/(me*v),Fe(1:Np,3), 'blue');
%title ('electrons distribution function');
xlabel ('p/{m_e c}');
ylabel ('F_e(p) p^4');
legend('t=0','t=T/2','t=T','Location','southeast');
grid ;

figure(3);
plot (Pa(1:Np,1)/(me*v),Fa(1:Np,1), 'red',Pa(1:Np,2)/(me*v),Fa(1:Np,2), 'green',Pa(1:Np,3)/(me*v),Fa(1:Np,3), 'blue');
title ('alphas distribution function');
xlabel ('p/{m_e c}');
ylabel ('F_\alpha(p) p^4');
legend('t=0','t=T/2','t=T','Location','southeast');
grid ;

figure(4);
plot (Ppos(1:Np,1)/(me*v),Fpos(1:Np,1), 'red',Ppos(1:Np,2)/(me*v),Fpos(1:Np,2), 'green',Ppos(1:Np,3)/(me*v),Fpos(1:Np,3), 'blue');
title ('positrons distribution function');
xlabel ('p/{m_e c}');
ylabel ('F_{e+}(p) p^4');
legend('t=0','t=T/2','t=T','Location','southeast');
grid ;