clear;
load traectory_electron.dat;
N1=1;
N2=size(traectory_electron,1);

xteoretic(1:N2) = 0;
Efield = 1E-6;
q = -4.803529695E-10;
m = 0.910938291E-27;

for i = 1:N2;
    xteoretic(i) = traectory_electron(1,3) + (q*Efield*0.5/m)*(traectory_electron(i,2))^2;
end;


figure(1);
%plot (traectory_electron(1:N2,1),traectory_electron(1:N2,2),'red',traectory_electron(1:N2,1),traectory_electron(1:N2,3),'green',traectory_electron(1:N2,1),traectory_electron(1:N2,4),'blue', traectory_electron(1:N2,1), xteoretic(1:N2), 'yellow');
plot (traectory_electron(1:N2,1),traectory_electron(1:N2,3),'red',traectory_electron(1:N2,1),traectory_electron(1:N2,4),'green',traectory_electron(1:N2,1),traectory_electron(1:N2,5),'blue');
title ('r');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('r cm');
legend('x', 'y','z','Location','southwest');
grid ;

figure(2);
plot (traectory_electron(1:N2,1),traectory_electron(1:N2,6),'red',traectory_electron(1:N2,1),traectory_electron(1:N2,7),'green',traectory_electron(1:N2,1),traectory_electron(1:N2,8),'blue', traectory_electron(1:N2,1),traectory_electron(1:N2,9),'black');
title ('p');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('p g*cm/s');
legend('px', 'py','pz', 'p','Location','southwest');
grid ;

figure(3);
plot (traectory_electron(1:N2,1),traectory_electron(1:N2,10),'red',traectory_electron(1:N2,1),traectory_electron(1:N2,11),'green',traectory_electron(1:N2,1),traectory_electron(1:N2,12),'blue');
title ('v');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('v cm/s');
legend('vx', 'vy','vz','Location','southwest');
grid ;

figure(4);
plot (traectory_electron(1:N2,4),traectory_electron(1:N2,5));
title ('r');
xlabel ('y cm');
ylabel ('z cm');
grid ;
