clear;
distribution_protons_grid = importdata('distribution_protons_grid_0.dat');
distribution_electrons_grid = importdata('distribution_electrons_grid_0.dat');
%distribution_alphas_grid = importdata('distribution_alphas_grid_5.dat');
%distribution_positrons_grid = importdata('distribution_positrons_grid_5.dat');
load initialParameters.dat;
load Xfile.dat;

Np = 500;

Nx = size(Xfile,1) - 1;
Nt = size(distribution_electrons_grid, 1)/(Nx + 1);
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
Pejuttner(1:Np)=0;

me = initialParameters(36);
mp = 1.67262177*10^-24;
ma = 6.64*10^-24;
v=2.998*10^10;
T = 2*10^15;
kB = 1.3806488*10^-16;
theta = kB*T/(me*v*v);
factor = 1;
minX = 10;
maxX = 80;
for i=1:Np,   
   Pp(i,1) = distribution_protons_grid(1 + a*(Nx + 1),i)/(mp*v);
   Pp(i,2) = distribution_protons_grid(1 + b*(Nx + 1),i)/(mp*v);
   Pp(i,3) = distribution_protons_grid(1 + c*(Nx + 1),i)/(mp*v);
   
   Pe(i,1) = distribution_electrons_grid(1 + a*(Nx + 1),i)/(me*v);
   Pe(i,2) = distribution_electrons_grid(1 + b*(Nx + 1),i)/(me*v);
   Pe(i,3) = distribution_electrons_grid(1 + c*(Nx + 1),i)/(me*v);
     
   %Pa(i,1) = distribution_alphas_grid(1 + a*(Nx + 1),i);
   %Pa(i,2) = distribution_alphas_grid(1 + b*(Nx + 1),i);
   %Pa(i,3) = distribution_alphas_grid(1 + c*(Nx + 1),i);
     
   %Ppos(i,1) = distribution_positrons_grid(1 + a*(Nx + 1),i);
   %Ppos(i,2) = distribution_positrons_grid(1 + b*(Nx + 1),i);
   %Ppos(i,3) = distribution_positrons_grid(1 + c*(Nx + 1),i);
   
   for j=minX:maxX,
     factor = Pp(i,1)*Pp(i,1);
     Fp(i,1) = Fp(i,1) + distribution_protons_grid(1 + j + a*(Nx + 1), i)*factor;
     factor = Pp(i,2)*Pp(i,2);
     Fp(i,2) = Fp(i,2) + distribution_protons_grid(1 + j + b*(Nx + 1), i)*factor;
     factor = Pp(i,3)*Pp(i,3);
     Fp(i,3) = Fp(i,3) + distribution_protons_grid(1 + j + c*(Nx + 1), i)*factor;
   
     factor = Pe(i,1)*Pe(i,1);
     Fe(i,1) = Fe(i,1) + distribution_electrons_grid(1 + j + a*(Nx + 1), i)*factor;
     factor = Pe(i,2)*Pe(i,2);
     Fe(i,2) = Fe(i,2) + distribution_electrons_grid(1 + j + b*(Nx + 1), i)*factor;
     factor = Pe(i,3)*Pe(i,3);
     Fe(i,3) = Fe(i,3) + distribution_electrons_grid(1 + j + c*(Nx + 1), i)*factor;
   
     %factor = Pa(i,1)*Pa(i,1);
     %Fa(i,1) = Fa(i,1) + distribution_alphas_grid(1 + j + a*(Nx + 1), i)*factor;
     %factor = Pa(i,2)*Pa(i,2);
     %Fa(i,2) = Fa(i,1) + distribution_alphas_grid(1 + j + b*(Nx + 1), i)*factor;
     %factor = Pa(i,3)*Pa(i,3);
     %Fa(i,3) = Fa(i,1) + distribution_alphas_grid(1 + j + c*(Nx + 1), i)*factor;
   
     %factor = Ppos(i,1)*Ppos(i,1);
     %Fpos(i,1) = Fpos(i,1) + distribution_positrons_grid(1 + j + a*(Nx + 1), i)*factor;
     %factor = Ppos(i,2)*Ppos(i,2);
     %Fpos(i,2) = Fpos(i,2) + distribution_positrons_grid(1 + j + b*(Nx + 1), i)*factor;
     %factor = Ppos(i,3)*Ppos(i,3);
     %Fpos(i,3) = Fpos(i,3) + distribution_positrons_grid(1 + j + c*(Nx + 1), i)*factor;
   end;
   exp1 = exp(-sqrt(1+Pe(i,1)*Pe(i,1))/theta);
   bes = besselk(2, 1/theta);
   p = Pe(i,1);
   p3 = (p/(me*v))^3;
   Pejuttner(i) = (1.0/(theta*bes))*exp1*p3*Pe(i,1);
end;
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 1);
figure(1);
plot (Pp(1:Np,1),Fp(1:Np,1), 'red',Pp(1:Np,2),Fp(1:Np,2), 'green',Pp(1:Np,3),Fp(1:Np,3), 'blue');
title ('protons distribution function');
xlabel ('{\gamma}-1');
ylabel ('F_p({\gamma}) ({\gamma}-1)');
legend('t=0','t=T/2','t=T','Location','southeast');
grid ;

figure(2);
plot (Pe(1:Np,1),Fe(1:Np,1), 'red',Pe(1:Np,2),Fe(1:Np,2), 'green',Pe(1:Np,3),Fe(1:Np,3), 'blue');
%plot (Pe(1:Np,1)/(me*v),Fe(1:Np,1), 'red',Pe(1:Np,2)/(me*v),Fe(1:Np,2), 'green',Pe(1:Np,3)/(me*v),Fe(1:Np,3), 'blue', Pe(1:Np,1)/(me*v), Pejuttner(1:Np), 'black');
title ('electrons distribution function');
xlabel ('{\gamma}-1');
ylabel ('F_e({\gamma}) ({\gamma}-1)');
legend('t=0','t=T/2','t=T','Location','southeast');
grid ;

%figure(3);
%plot (Ga(1:Np,1)-1,Fa(1:Np,1), 'red',Ga(1:Np,2)-1,Fa(1:Np,2), 'green',Ga(1:Np,3)-1,Fa(1:Np,3), 'blue');
%title ('alphas distribution function');
%xlabel ('{\gamma}-1');
%ylabel ('F_alpha({\gamma}) ({\gamma}-1)');
%legend('t=0','t=T/2','t=T','Location','southeast');
%grid ;

%figure(4);
%plot (Gpos(1:Np,1)-1,Fpos(1:Np,1), 'red',Gpos(1:Np,2)-1,Fpos(1:Np,2), 'green',Gpos(1:Np,3)-1,Fpos(1:Np,3), 'blue');
%title ('positrons distribution function');
%xlabel ('{\gamma}-1');
%ylabel ('F_{e+}({\gamma}) ({\gamma}-1)');
%legend('t=0','t=T/2','t=T','Location','southeast');
%grid ;