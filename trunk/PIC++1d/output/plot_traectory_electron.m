clear;
load traectory_electron.dat;
N1=1;
N2=size(traectory_electron,1);

xteoretic(1:N2) = 0;
Efield = 1E-6;
q = -4.803529695E-10;
m = 0.910938291E-27;

for i = 1:N2;
    xteoretic(i) = traectory_electron(1,2) + (q*Efield*0.5/m)*(traectory_electron(i,1))^2;
end;


figure(1);
%plot (traectory_electron(1:N2,1),traectory_electron(1:N2,2),'red',traectory_electron(1:N2,1),traectory_electron(1:N2,3),'green',traectory_electron(1:N2,1),traectory_electron(1:N2,4),'blue', traectory_electron(1:N2,1), xteoretic(1:N2), 'yellow');
plot (traectory_electron(1:N2,1),traectory_electron(1:N2,2),'red',traectory_electron(1:N2,1),traectory_electron(1:N2,3),'green',traectory_electron(1:N2,1),traectory_electron(1:N2,4),'blue');
title ('r');
xlabel ('t');
ylabel (')');
legend('x', 'y','z','Location','southwest');
grid ;

figure(2);
plot (traectory_electron(1:N2,1),traectory_electron(1:N2,5),'red',traectory_electron(1:N2,1),traectory_electron(1:N2,6),'green',traectory_electron(1:N2,1),traectory_electron(1:N2,7),'blue', traectory_electron(1:N2,1),traectory_electron(1:N2,8),'black');
title ('p');
xlabel ('t');
ylabel ('p');
legend('px', 'py','pz', 'p','Location','southwest');
grid ;

figure(3);
plot (traectory_electron(1:N2,3),traectory_electron(1:N2,4));
title ('r');
xlabel ('y');
ylabel ('z');
grid ;
