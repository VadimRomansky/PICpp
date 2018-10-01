clear;
load distribution_protons_grid_20.dat;
load distribution_electrons_grid_20.dat;
load initialParameters.dat;
load Xfile.dat;

directory_name = '../../../tristan-mp-pitp/output/';
file_name = 'spect';
file_number = '.020';
full_name = strcat(directory_name, file_name, file_number);
Np = 500;
Nx = size(Xfile,1) - 1;

fpB = hdf5read(full_name,'specp');
feB = hdf5read(full_name,'spece');
gB=hdf5read(full_name,'gamma');

NxB = size(fpB,1);
NpB = size(fpB,2);

Fp(1:Np) = 0;
Fe(1:Np) = 0;

Pp(1:Np) = 0;
Pe(1:Np) = 0;
Pejuttner(1:Np)=0;

FpB(1:Np) = 0;
FeB(1:Np) = 0;

PpB(1:Np) = 0;
PeB(1:Np) = 0;

me = initialParameters(36);
mp = 1.67262177*10^-24;
ma = 6.64*10^-24;
v=2.998*10^10;
T = 2*10^15;
kB = 1.3806488*10^-16;
theta = kB*T/(me*v*v);
factor = 1;

minX = 1;
maxX = Nx/2;

minXB = 1;
maxXB = NxB;

a = 0;
for i=1:Np,   
   Pp(i) = distribution_protons_grid_20(1 + a*(Nx + 1),i)/(mp*v);
   
   Pe(i) = distribution_electrons_grid_20(1 + a*(Nx + 1),i)/(me*v);

   
   for j=minX:maxX,
     factor = Pp(i)*Pp(i);
     Fp(i) = Fp(i) + distribution_protons_grid_20(1 + j + a*(Nx + 1), i)*factor;

     factor = Pe(i)*Pe(i);
     Fe(i) = Fe(i) + distribution_electrons_grid_20(1 + j + a*(Nx + 1), i)*factor;
   end;
   exp1 = exp(-sqrt(1+Pe(i)*Pe(i)/(me*me*v*v))/theta);
   bes = besselk(2, 1/theta);
   p = Pe(i);
   p3 = (p/(me*v))^3;
   Pejuttner(i) = (1.0/(theta*bes))*exp1*p3*Pe(i);
end;

for i = 1:NpB,
    %Pp(i) = sqrt((g(i)+1)^2 - 1)*mp*c;
    %Pe(i) = sqrt((g(i)+1)^2 - 1)*me*c;
    PpB(i) = sqrt((gB(i)+1)^2 - 1);
    PeB(i) = sqrt((gB(i)+1)^2 - 1);
    for j = 1:NxB,
        FpB(i) = FpB(i) + fpB(j,i);
        FeB(i) = FeB(i) + feB(j,i);
    end;
    FpB(i)=FpB(i)*(PpB(i)^3)/(1+gB(i));
    FeB(i)=FeB(i)*(PeB(i)^3)/(1+gB(i));
end;

norm = 1.0;

normp = (Fp(1)/(Pp(2)^2))*(Pp(2) - Pp(1));
norme = (Fe(1)/(Pe(2)^2))*(Pe(2) - Pe(1));

for i = 2:Np,
    normp = normp + (Fp(i)/(Pp(i)^2))*(Pp(i) - Pp(i-1));
    norme = norme + (Fe(i)/(Pe(i)^2))*(Pe(i) - Pe(i-1));
end;

for i = 1:Np,
    Fp(i) = Fp(i)*norm/normp;
    Fe(i) = Fe(i)*norm/norme;
end;


normpB = (FpB(1)/(PpB(2)^2))*(PpB(2) - PpB(1));
normeB = (FeB(1)/(PeB(2)^2))*(PeB(2) - PeB(1));

for i = 2:NpB,
    normpB = normpB + (FpB(i)/(PpB(i)^2))*(PpB(i) - PpB(i-1));
    normeB = normeB + (FeB(i)/(PeB(i)^2))*(PeB(i) - PeB(i-1));
end;

for i = 1:NpB,
    FpB(i) = FpB(i)*norm/normpB;
    FeB(i) = FeB(i)*norm/normeB;
end;


set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
plot (Pp(1:Np),Fp(1:Np), 'red',PpB(1:NpB),FpB(1:NpB), 'blue');
title ('protons distribution function');
xlabel ('p/mc');
ylabel ('F_p(p) p^4');
legend('PIDARAC','Tristan','Location','southeast');
grid ;

figure(2);
plot (Pe(1:Np),Fe(1:Np), 'red',PeB(1:NpB),FeB(1:NpB), 'blue');
%plot (Pe(1:Np,1)/(me*v),Fe(1:Np,1), 'red',Pe(1:Np,2)/(me*v),Fe(1:Np,2), 'green',Pe(1:Np,3)/(me*v),Fe(1:Np,3), 'blue', Pe(1:Np,1)/(me*v), Pejuttner(1:Np), 'black');
title ('electrons distribution function');
xlabel ('p/mc');
ylabel ('F_e(p) p^4');
legend('PIDARAC','Tristan','Location','southeast');
grid ;