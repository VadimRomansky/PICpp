clear;
load distribution_protons.dat;
load distribution_electrons.dat;
load distribution_protons_rel.dat;
load distribution_electrons_rel.dat;

Np = 500;

Nt(1:2) = 11;

Nt(1) = fix(size(distribution_protons,1)/Np) - 1;
Nt(2) = fix(size(distribution_protons_rel,1)/Np) - 1;


Fp(1:Np, 1:2) = 0;
Fe(1:Np, 1:2) = 0;


Pp(1:Np, 1:2) = 0;
Pe(1:Np, 1:2) = 0;

me = 0.91*10^-25;
mp = 1.6*10^-24;
v=2.998*10^10;

%Color = {'r','green','blue','black','yellow','cyan','magenta','purple','grey'};
Color = {'red','blue'};

for i=1:Np,
   Pp(i,1) = distribution_protons(i + Nt(1)*Np,1);
   Pp(i,2) = distribution_protons_rel(i + Nt(2)*Np,1);
   
   Fp(i,1) = distribution_protons(i + Nt(1)*Np, 2)*Pp(i,1)*Pp(i,1);
   Fp(i,2) = distribution_protons_rel(i + Nt(2)*Np, 2)*Pp(i,2)*Pp(i,2);
   
   Pe(i,1) = distribution_electrons(i + Nt(1)*Np,1);
   Pe(i,2) = distribution_electrons_rel(i + Nt(2)*Np,1);
   
   Fe(i,1) = distribution_electrons(i + Nt(1)*Np, 2)*Pe(i,1)*Pe(i,1);
   Fe(i,2) = distribution_electrons_rel(i + Nt(2)*Np, 2)*Pe(i,2)*Pe(i,2);
end;
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
for j=1:2,
    plot (Pp(1:Np,j)/(mp*v),Fp(1:Np,j),'color',Color{j});
end;
xlabel ('p/{m_p c}');
ylabel ('F_p(p) p^4');
legend('{\gamma} = 1.5','{\gamma} = 10','Location','southeast');
grid ;

figure(2);
hold on;
for j=1:2,
    plot (Pe(1:Np,j)/(mp*v),Fe(1:Np,j),'color',Color{j});
end;
xlabel ('p/{m_p c}');
ylabel ('F_e(p) p^4');
legend('{\gamma} = 1.5','{\gamma} = 10','Location','southeast');
grid ;