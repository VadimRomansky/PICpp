clear;
trajectory = importdata('./output/trajectory_3.dat');

Nt = size(trajectory,1);

figure(1);
plot (trajectory(1:Nt,2),trajectory(1:Nt,1), 'red');
title ('xt');
xlabel ('x');
ylabel ('t');
grid ;

figure(2);
plot (trajectory(1:Nt,2),trajectory(1:Nt,3), 'red');
title ('xy');
xlabel ('x');
ylabel ('y');
grid ;

figure(3);
plot (trajectory(1:Nt,1),trajectory(1:Nt,5), 'red');
title ('pt');
xlabel ('t');
ylabel ('px');
grid ;