clear;
load 'average_distribution.dat';

Np = size(average_distribution,1);

Pe(1:Np) = 0;
Fe(1:Np) = 0;
for i = 1:Np,
    Pe(i) = average_distribution(i,1);
    Fe(i) = average_distribution(i,2);
end;

figure(1);
loglog(Pe(1:Np),Fe(1:Np),'--','color','red');
%plot (Pp(1:Np),Fp(1:Np), 'red',Pp(1:Np), Fpjuttner(1:Np), 'blue');
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('Fe*p^4');
legend('Fe','Location','southeast');
grid ;

dlmwrite('Pe.dat',Pe,'delimiter','\n');
dlmwrite('Fe.dat',Fe,'delimiter','\n');