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

Nt = 2;


c = Nt - 1;

Fp(1:Np, 1:10) = 0;
Fe(1:Np, 1:10) = 0;


Pp(1:Np, 1:10) = 0;
Pe(1:Np, 1:10) = 0;

me = 0.91*10^-25;
mp = 1.6*10^-24;
v=2.998*10^10;

%Color = {'r','green','blue','black','yellow','cyan','magenta','purple','grey'};
Color = {[.7,.3,.3],'red','green','blue','black','yellow','cyan','magenta',[.5,.5,.5],[.3,.7,.3]};

for i=1:Np,
   %Pp(i,1) = distribution_protons0(i + c*Np,1);
   Pp(i,2) = distribution_protons1(i + c*Np,1);
   Pp(i,3) = distribution_protons2(i + c*Np,1);
   Pp(i,4) = distribution_protons3(i + c*Np,1);
   Pp(i,5) = distribution_protons4(i + c*Np,1);
   Pp(i,6) = distribution_protons5(i + c*Np,1);
   Pp(i,7) = distribution_protons6(i + c*Np,1);
   Pp(i,8) = distribution_protons7(i + c*Np,1);
   Pp(i,9) = distribution_protons8(i + c*Np,1);
   Pp(i,10) = distribution_protons9(i + c*Np,1);
   
   %Fp(i,1) = distribution_protons0(i + c*Np, 2)*Pp(i,1)*Pp(i,1);
   Fp(i,2) = distribution_protons1(i + c*Np, 2)*Pp(i,2)*Pp(i,2);
   Fp(i,3) = distribution_protons2(i + c*Np, 2)*Pp(i,3)*Pp(i,3);
   Fp(i,4) = distribution_protons3(i + c*Np, 2)*Pp(i,4)*Pp(i,4);
   Fp(i,5) = distribution_protons4(i + c*Np, 2)*Pp(i,5)*Pp(i,5);
   Fp(i,6) = distribution_protons5(i + c*Np, 2)*Pp(i,6)*Pp(i,6);
   Fp(i,7) = distribution_protons6(i + c*Np, 2)*Pp(i,7)*Pp(i,7);
   Fp(i,8) = distribution_protons7(i + c*Np, 2)*Pp(i,8)*Pp(i,8);
   Fp(i,9) = distribution_protons8(i + c*Np, 2)*Pp(i,9)*Pp(i,9);
   Fp(i,10) = distribution_protons9(i + c*Np, 2)*Pp(i,10)*Pp(i,10);
   
   %Pe(i,1) = distribution_electrons0(i + c*Np,1);
   Pe(i,2) = distribution_electrons1(i + c*Np,1);
   Pe(i,3) = distribution_electrons2(i + c*Np,1);
   Pe(i,4) = distribution_electrons3(i + c*Np,1);
   Pe(i,5) = distribution_electrons4(i + c*Np,1);
   Pe(i,6) = distribution_electrons5(i + c*Np,1);
   Pe(i,7) = distribution_electrons6(i + c*Np,1);
   Pe(i,8) = distribution_electrons7(i + c*Np,1);
   Pe(i,9) = distribution_electrons8(i + c*Np,1);
   Pe(i,10) = distribution_electrons9(i + c*Np,1);
   
   %Fe(i,1) = distribution_electrons0(i + c*Np, 2)*Pe(i,1)*Pe(i,1);
   Fe(i,2) = distribution_electrons1(i + c*Np, 2)*Pe(i,2)*Pe(i,2);
   Fe(i,3) = distribution_electrons2(i + c*Np, 2)*Pe(i,3)*Pe(i,3);
   Fe(i,4) = distribution_electrons3(i + c*Np, 2)*Pe(i,4)*Pe(i,4);
   Fe(i,5) = distribution_electrons4(i + c*Np, 2)*Pe(i,5)*Pe(i,5);
   Fe(i,6) = distribution_electrons5(i + c*Np, 2)*Pe(i,6)*Pe(i,6);
   Fe(i,7) = distribution_electrons6(i + c*Np, 2)*Pe(i,7)*Pe(i,7);
   Fe(i,8) = distribution_electrons7(i + c*Np, 2)*Pe(i,8)*Pe(i,8);
   Fe(i,9) = distribution_electrons8(i + c*Np, 2)*Pe(i,9)*Pe(i,9);
   Fe(i,10) = distribution_electrons9(i + c*Np, 2)*Pe(i,10)*Pe(i,10);
end;
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
for j=2:10,
    plot (Pp(1:Np,j)/(mp*v),Fp(1:Np,j),'color',Color{j});
end;
xlabel ('p/{m_p c}');
ylabel ('F_p(p) p^4');
legend('{\theta} = 10','{\theta} = 20','{\theta} = 30','{\theta} = 40','{\theta} = 50','{\theta} = 60','{\theta} = 70','{\theta} = 80','{\theta} = 90','Location','southeast');
grid ;

figure(2);
hold on;
for j=2:10,
    plot (Pe(1:Np,j)/(mp*v),Fe(1:Np,j),'color',Color{j});
end;
xlabel ('p/{m_p c}');
ylabel ('F_e(p) p^4');
legend('{\theta} = 10','{\theta} = 20','{\theta} = 30','{\theta} = 40','{\theta} = 50','{\theta} = 60','{\theta} = 70','{\theta} = 80','{\theta} = 90','Location','southeast');
grid ;