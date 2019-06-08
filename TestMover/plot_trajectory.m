clear;
load trajectory.dat;

Nt = size(trajectory,1);

figure(1);
plot (trajectory(1:Nt,1),1:Nt, 'red');
title ('xt');
xlabel ('x');
ylabel ('t');
grid ;

figure(2);
plot (trajectory(1:Nt,1),trajectory(1:Nt,2), 'red');
title ('xy');
xlabel ('x');
ylabel ('y');
grid ;

figure(3);
plot (1:Nt,trajectory(1:Nt,7), 'red');
title ('tv');
xlabel ('t');
ylabel ('v');
grid ;