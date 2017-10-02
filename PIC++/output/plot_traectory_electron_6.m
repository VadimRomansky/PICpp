clear;
load trajectory_electron_6.dat;
N1=1;
N2=size(trajectory_electron_6,1);

xteoretic(1:N2) = 0;
Efield = 1E-6;
q = -4.803529695E-10;
m = 0.910938291E-27;

for i = 1:N2;
    xteoretic(i) = trajectory_electron_6(1,3) + (q*Efield*0.5/m)*(trajectory_electron_6(i,2))^2;
end;


figure(1);
%plot (traectory_electron(1:N2,1),traectory_electron(1:N2,3),'red',traectory_electron(1:N2,1),traectory_electron(1:N2,4),'green',traectory_electron(1:N2,1),traectory_electron(1:N2,5),'blue', traectory_electron(1:N2,1), xteoretic(1:N2), 'yellow');
plot (trajectory_electron_6(1:N2,1),trajectory_electron_6(1:N2,3),'red',trajectory_electron_6(1:N2,1),trajectory_electron_6(1:N2,4),'green',trajectory_electron_6(1:N2,1),trajectory_electron_6(1:N2,5),'blue');
title ('r');
xlabel ('t');
ylabel ('r');
legend('x', 'y','z','Location','southwest');
grid ;

figure(2);
plot (trajectory_electron_6(1:N2,1),trajectory_electron_6(1:N2,6),'red',trajectory_electron_6(1:N2,1),trajectory_electron_6(1:N2,7),'green',trajectory_electron_6(1:N2,1),trajectory_electron_6(1:N2,8),'blue', trajectory_electron_6(1:N2,2),trajectory_electron_6(1:N2,9),'black');
title ('p');
xlabel ('t');
ylabel ('p');
legend('px', 'py','pz', 'p','Location','southwest');
grid ;

figure(3);
plot (trajectory_electron_6(1:N2,4),trajectory_electron_6(1:N2,5));
title ('r');
xlabel ('y');
ylabel ('z');
grid ;

figure(4);
plot (trajectory_electron_6(1:N2,3),trajectory_electron_6(1:N2,9));
title ('p(x)');
xlabel ('x');
ylabel ('p');
grid ;

figure(5);
plot (trajectory_electron_6(1:N2,3),trajectory_electron_6(1:N2,6));
title ('p_x(x)');
xlabel ('x');
ylabel ('p_x');
grid ;
