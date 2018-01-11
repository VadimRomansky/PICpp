clear;
%load distribution_protons0.dat;
%load distribution_electrons0.dat;
load distribution_protons1.dat;
load distribution_electrons1.dat;
load distribution_protons2.dat;
load distribution_electrons2.dat;
load distribution_protons3.dat;
load distribution_electrons3.dat;
load distribution_protons4.dat;
load distribution_electrons4.dat;
load distribution_protons5.dat;
load distribution_electrons5.dat;
load distribution_protons6.dat;
load distribution_electrons6.dat;
load distribution_protons7.dat;
load distribution_electrons7.dat;
load distribution_protons8.dat;
load distribution_electrons8.dat;
load distribution_protons9.dat;
load distribution_electrons9.dat;


Np = 500;

Nt(1:11) = 0;

%Nt(1) = fix(size(distribution_protons0,1)/Np) - 1;
Nt(2) = fix(size(distribution_protons1,1)/Np) - 1;
Nt(3) = fix(size(distribution_protons2,1)/Np) - 1;
Nt(4) = fix(size(distribution_protons3,1)/Np) - 1;
Nt(5) = fix(size(distribution_protons4,1)/Np) - 1;
Nt(6) = fix(size(distribution_protons5,1)/Np) - 1;
Nt(7) = fix(size(distribution_protons6,1)/Np) - 1;
Nt(8) = fix(size(distribution_protons7,1)/Np) - 1;
Nt(9) = fix(size(distribution_protons8,1)/Np) - 1;
Nt(10) = fix(size(distribution_protons9,1)/Np) - 1;

Nt(4) = 24;
Nt(5) = 8;
Nt(6) = 9;
Nt(7) = 10;
Nt(8) = 12;
Nt(9) = 14;
Nt(10) = 24;
Fp(1:Np, 1:10) = 0;
Fe(1:Np, 1:10) = 0;

Pp(1:Np, 1:10) = 0;
Pe(1:Np, 1:10) = 0;

upstreamFp(1:Np) = 0;
upstreamFe(1:Np) = 0;
upstreamPp(1:Np) = 0;
upstreamPe(1:Np) = 0;

downstreamFp(1:Np) = 0;
downstreamFe(1:Np) = 0;
downstreamPp(1:Np) = 0;
downstreamPe(1:Np) = 0;

me = 0.910938356*10^-25;
mp = 1.672621*10^-24;
cv = 2.99792458*10^10;
kBoltzman = 1.38064852*10^-16;
T=5*10^8;

gamma = 1.5;

beta = sqrt(1 - 1/(gamma*gamma));
v = beta*cv;

%Color = {'r','green','blue','black','yellow','cyan','magenta','purple','grey'};
Color = {[.7,.3,.3],'red','green','blue','black','cyan','magenta','yellow',[.5,.5,.5],[.3,.7,.3],[1.0,.5,0],[.75,0.0,.7]};

for i=1:Np,
   %Pp(i,1) = distribution_protons0(i + c*Np,1);
   Pp(i,2) = distribution_protons1(i + Nt(2)*Np,1);
   Pp(i,3) = distribution_protons2(i + Nt(3)*Np,1);
   Pp(i,4) = distribution_protons3(i + Nt(4)*Np,1);
   Pp(i,5) = distribution_protons4(i + Nt(5)*Np,1);
   Pp(i,6) = distribution_protons5(i + Nt(6)*Np,1);
   Pp(i,7) = distribution_protons6(i + Nt(7)*Np,1);
   Pp(i,8) = distribution_protons7(i + Nt(8)*Np,1);
   Pp(i,9) = distribution_protons8(i + Nt(9)*Np,1);
   Pp(i,10) = distribution_protons9(i + Nt(10)*Np,1);
   
   %Fp(i,1) = distribution_protons0(i + c*Np, 2)*Pp(i,1)*Pp(i,1);
   Fp(i,2) = distribution_protons1(i + Nt(2)*Np, 2)*Pp(i,2)*Pp(i,2);
   Fp(i,3) = distribution_protons2(i + Nt(3)*Np, 2)*Pp(i,3)*Pp(i,3);
   Fp(i,4) = distribution_protons3(i + Nt(4)*Np, 2)*Pp(i,4)*Pp(i,4);
   Fp(i,5) = distribution_protons4(i + Nt(5)*Np, 2)*Pp(i,5)*Pp(i,5);
   Fp(i,6) = distribution_protons5(i + Nt(6)*Np, 2)*Pp(i,6)*Pp(i,6);
   Fp(i,7) = distribution_protons6(i + Nt(7)*Np, 2)*Pp(i,7)*Pp(i,7);
   Fp(i,8) = distribution_protons7(i + Nt(8)*Np, 2)*Pp(i,8)*Pp(i,8);
   Fp(i,9) = distribution_protons8(i + Nt(9)*Np, 2)*Pp(i,9)*Pp(i,9);
   Fp(i,10) = distribution_protons9(i + Nt(10)*Np, 2)*Pp(i,10)*Pp(i,10);
   
   %Pe(i,1) = distribution_electrons0(i + c*Np,1);
   Pe(i,2) = distribution_electrons1(i + Nt(2)*Np,1);
   Pe(i,3) = distribution_electrons2(i + Nt(3)*Np,1);
   Pe(i,4) = distribution_electrons3(i + Nt(4)*Np,1);
   Pe(i,5) = distribution_electrons4(i + Nt(5)*Np,1);
   Pe(i,6) = distribution_electrons5(i + Nt(6)*Np,1);
   Pe(i,7) = distribution_electrons6(i + Nt(7)*Np,1);
   Pe(i,8) = distribution_electrons7(i + Nt(8)*Np,1);
   Pe(i,9) = distribution_electrons8(i + Nt(9)*Np,1);
   Pe(i,10) = distribution_electrons9(i + Nt(10)*Np,1);
   
   %Fe(i,1) = distribution_electrons0(i + c*Np, 2)*Pe(i,1)*Pe(i,1);
   Fe(i,2) = distribution_electrons1(i + Nt(2)*Np, 2)*Pe(i,2)*Pe(i,2);
   Fe(i,3) = distribution_electrons2(i + Nt(3)*Np, 2)*Pe(i,3)*Pe(i,3);
   Fe(i,4) = distribution_electrons3(i + Nt(4)*Np, 2)*Pe(i,4)*Pe(i,4);
   Fe(i,5) = distribution_electrons4(i + Nt(5)*Np, 2)*Pe(i,5)*Pe(i,5);
   Fe(i,6) = distribution_electrons5(i + Nt(6)*Np, 2)*Pe(i,6)*Pe(i,6);
   Fe(i,7) = distribution_electrons6(i + Nt(7)*Np, 2)*Pe(i,7)*Pe(i,7);
   Fe(i,8) = distribution_electrons7(i + Nt(8)*Np, 2)*Pe(i,8)*Pe(i,8);
   Fe(i,9) = distribution_electrons8(i + Nt(9)*Np, 2)*Pe(i,9)*Pe(i,9);
   Fe(i,10) = distribution_electrons9(i + Nt(10)*Np, 2)*Pe(i,10)*Pe(i,10);
end;

weightP = 0.1;
weightE = 0.1;

upstreamPp(1) = 0.9*beta*gamma*mp*cv;
upstreamPe(1) = 0.9*beta*gamma*me*cv;
for i=1:Np,
    upstreamPp(i) = upstreamPp(1) + (i-1)*0.2*beta*gamma*mp*cv/Np;
    upstreamPe(i) = upstreamPe(1) + (i-1)*0.2*beta*gamma*me*cv/Np;
    downstreamPp(i) = Pp(i,2);
    downstreamPe(i) = Pe(i,2);
end;

denomP = (2*pi*mp*kBoltzman*T)^-1.5;
denomE = (2*pi*me*kBoltzman*T)^-1.5;

for i=1:Np,
    p = upstreamPp(i);
    upstreamFp(i) = weightP*(denomP*pi*2*kBoltzman*T/(beta*gamma*cv))*p*(exp(-(p-beta*gamma*mp*cv)^2/(2*mp*kBoltzman*T))-exp(-(p+beta*gamma*mp*cv)^2/(2*mp*kBoltzman*T)))*p^2;
    p = upstreamPe(i);
    upstreamFe(i) = weightE*(denomE*pi*2*kBoltzman*T/(beta*gamma*cv))*p*(exp(-(p-beta*gamma*me*cv)^2/(2*me*kBoltzman*T))-exp(-(p+beta*gamma*me*cv)^2/(2*me*kBoltzman*T)))*p^2;
end;

%besselKp = 3.669760648*10^-9460;
%besselKe = 3.142420932*10^-517;
Tp = 0.3*10^13;
Te = 0.3*10^12;
besselKp = besselk(2, (mp*cv^2)/(kBoltzman*Tp));
besselKe = besselk(2, (me*cv^2)/(kBoltzman*Te));
denomP2 = 1/(4*pi*(mp^3)*(cv^3)*(kBoltzman*Tp/(mp*cv^2))*besselKp);
denomE2 = 1/(4*pi*(me^3)*(cv^3)*(kBoltzman*Te/(me*cv^2))*besselKe);
%denomP2 = 1/(4*pi*mp^3*cv^3*(kBoltzman*T/(mp*cv^2))*besselKp);
%denomE2 = 1/(4*pi*me^3*cv^3*(kBoltzman*T/(me*cv^2))*besselKe);
for i=1:Np,
    p = downstreamPp(i);
    a = sqrt(1+(p/(mp*cv))^2)*mp*cv^2/(kBoltzman*Tp);
    b = exp(-a);
    downstreamFp(i) = (1-weightP)*denomP2*exp(-sqrt(1+(p/(mp*cv))^2)*mp*cv^2/(kBoltzman*Tp))*p^4;
    p = downstreamPe(i);
    downstreamFe(i) = (1-weightE)*denomE2*exp(-sqrt(1+(p/(me*cv))^2)*me*cv^2/(kBoltzman*Te))*p^4;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
for j=2:10,
    plot (Pp(1:Np,j)/(mp*cv),Fp(1:Np,j),'color',Color{j});
end;
plot(upstreamPp(1:Np)/(mp*cv), upstreamFp(1:Np),'color', Color{11});
plot(downstreamPp(1:Np)/(mp*cv), downstreamFp(1:Np),'color', Color{12});
xlabel ('p/{m_p c}');
ylabel ('F_p(p) p^4');
legend('{\theta} = 10','{\theta} = 20','{\theta} = 30','{\theta} = 40','{\theta} = 50','{\theta} = 60','{\theta} = 70','{\theta} = 80','{\theta} = 90', 'upstream maxwell','downstream maxwell-juttner','Location','southeast');
grid ;

figure(2);
hold on;
for j=2:10,
    plot (Pe(1:Np,j)/(mp*cv),Fe(1:Np,j),'color',Color{j});
end;
plot(upstreamPe(1:Np)/(mp*cv), upstreamFe(1:Np),'color', Color{11});
plot(downstreamPe(1:Np)/(mp*cv), downstreamFe(1:Np),'color', Color{12});
xlabel ('p/{m_p c}');
ylabel ('F_e(p) p^4');
legend('{\theta} = 10','{\theta} = 20','{\theta} = 30','{\theta} = 40','{\theta} = 50','{\theta} = 60','{\theta} = 70','{\theta} = 80','{\theta} = 90','upstream maxwell', 'downstream maxwell-juttner','Location','southeast');
grid ;