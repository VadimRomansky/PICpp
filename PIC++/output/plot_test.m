clear;
load test.dat;
N1=1;
N2=size(test,1);


figure(1);
plot (test(1:N2,1),test(1:N2,2),'red',test(1:N2,1),test(1:N2,3),'blue');
title ('f');
xlabel ('x');
ylabel ('f(x)');
legend(3, 'simulation', 'theory');
grid ;
