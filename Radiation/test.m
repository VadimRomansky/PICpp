load Pe8.dat;
load Fe8.dat;
load Pe9.dat;
load Fe9.dat;

N= size(Pe8,2);
figure(1);
hold on;
plot(Pe8(1:N), Fe8(1:N), 'red', Pe9(1:N), Fe9(1:N), 'blue');
grid;