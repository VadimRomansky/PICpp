clear;
turbulentB = importdata('turbulentB.dat');
N = size(turbulentB,1);

figure(1);
hold on;
title ('B');
xlabel ('r');
ylabel ('B');
plot(1:N, turbulentB(1:N,1));
grid;

figure(2);
hold on;
title ('{\theta}');
xlabel ('r');
ylabel ('{\theta}');
plot(1:N, turbulentB(1:N,5));
grid;

figure(3);
hold on;
title ('{\lambda}');
xlabel ('r');
ylabel ('{\lambda}');
plot(1:N, turbulentB(1:N,6));
grid;

figure(4);
hold on;
title ('{\delta}');
xlabel ('r');
ylabel ('{\delta}');
plot(1:N, turbulentB(1:N,7));
grid;