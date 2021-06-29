clear;
turbulentB = importdata('turbulentB.dat');
N = size(turbulentB,1);

figure(1);
hold on;
title ('B');
xlabel ('r');
ylabel ('B');
plot(100:N, turbulentB(100:N,2));
grid;

figure(2);
hold on;
title ('{\lambda}');
xlabel ('r');
ylabel ('{\lambda}');
plot(100:N, turbulentB(100:N,3));
grid;

figure(3);
hold on;
title ('{\delta}');
xlabel ('r');
ylabel ('{\delta}');
plot(100:N, turbulentB(100:N,4));
grid;