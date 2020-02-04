clear;

load 'B.dat'

N = size(B,1);

figure(1);
hold on;
ylabel('Bz');
xlabel('x');
plot(1:N,B(1:N,3));
grid on;

figure(2);
hold on;
ylabel('Theta');
xlabel('x');
plot(1:N,B(1:N,4));
grid on;