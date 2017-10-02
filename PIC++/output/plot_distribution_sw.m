clear;
load distribution_protons_sw.dat;
load distribution_electrons_sw.dat;
load distribution_alphas_sw.dat;

Np = 500;

Nt = size(distribution_electrons_sw, 1)/Np;
%Nt = 28;


a = 0;
b = fix(Nt/2);
c = Nt - 1;

Fp(1:Np, 1:3) = 0;
Fe(1:Np, 1:3) = 0;
Fa(1:Np, 1:3) = 0;

Pp(1:Np, 1:3) = 0;
Pe(1:Np, 1:3) = 0;
Pa(1:Np, 1:3) = 0;

me = 9.1*10^-26;
v=3*10^10;

for i=1:Np,
   Pp(i,1) = distribution_protons_sw(i + a*Np,1);
   Pp(i,2) = distribution_protons_sw(i + b*Np,1);
   Pp(i,3) = distribution_protons_sw(i + c*Np,1);
   
   Fp(i,1) = distribution_protons_sw(i + a*Np, 2)*Pp(i,1)*Pp(i,1);
   Fp(i,2) = distribution_protons_sw(i + b*Np, 2)*Pp(i,2)*Pp(i,2);
   Fp(i,3) = distribution_protons_sw(i + c*Np, 2)*Pp(i,3)*Pp(i,3);
   
   Pe(i,1) = distribution_electrons_sw(i + a*Np,1);
   Pe(i,2) = distribution_electrons_sw(i + b*Np,1);
   Pe(i,3) = distribution_electrons_sw(i + c*Np,1);
   
   Fe(i,1) = distribution_electrons_sw(i + a*Np, 2)*Pe(i,1)*Pe(i,1);
   Fe(i,2) = distribution_electrons_sw(i + b*Np, 2)*Pe(i,2)*Pe(i,2);
   Fe(i,3) = distribution_electrons_sw(i + c*Np, 2)*Pe(i,3)*Pe(i,3);
   
   Pa(i,1) = distribution_alphas_sw(i + a*Np,1);
   Pa(i,2) = distribution_alphas_sw(i + b*Np,1);
   Pa(i,3) = distribution_alphas_sw(i + c*Np,1);
   
   Fa(i,1) = distribution_alphas_sw(i + a*Np, 2)*Pa(i,1)*Pa(i,1);
   Fa(i,2) = distribution_alphas_sw(i + b*Np, 2)*Pa(i,2)*Pa(i,2);
   Fa(i,3) = distribution_alphas_sw(i + c*Np, 2)*Pa(i,3)*Pa(i,3);
end;
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
plot (Pp(1:Np,1)/(me*v),Fp(1:Np,1), 'red',Pp(1:Np,2)/(me*v),Fp(1:Np,2), 'green',Pp(1:Np,3)/(me*v),Fp(1:Np,3), 'blue');
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