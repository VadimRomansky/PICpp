clear;
load trajectory_proton_4.dat;
N1=1;
N2=size(trajectory_proton_4,1);


figure(1);
plot (trajectory_proton_4(1:N2,1),trajectory_proton_4(1:N2,3),'red',trajectory_proton_4(1:N2,1),trajectory_proton_4(1:N2,4),'green',trajectory_proton_4(1:N2,1),trajectory_proton_4(1:N2,5),'blue');
title ('r');
xlabel ('t');
ylabel (')');
legend(4, 'x', 'y','z');
grid ;

figure(2);
plot (trajectory_proton_4(1:N2,1),trajectory_proton_4(1:N2,6),'red',trajectory_proton_4(1:N2,1),trajectory_proton_4(1:N2,7),'green',trajectory_proton_4(1:N2,1),trajectory_proton_4(1:N2,8),'blue', trajectory_proton_4(1:N2,1),trajectory_proton_4(1:N2,9),'black');
title ('p');
xlabel ('t');
ylabel ('p');
legend(4, 'px', 'py','pz', 'p');
grid ;

figure(3);
plot (trajectory_proton_4(1:N2,4),trajectory_proton_4(1:N2,5));
title ('r');
xlabel ('y');
ylabel ('z');
grid ;

figure(4);
plot (trajectory_proton_4(1:N2,3),trajectory_proton_4(1:N2,9));
title ('p(x)');
xlabel ('x');
ylabel ('p');
grid ;

figure(5);
plot (trajectory_proton_4(1:N2,3),trajectory_proton_4(1:N2,6));
title ('p_x(x)');
xlabel ('x');
ylabel ('p_x');
grid ;

