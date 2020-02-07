clear;

B = importdata('B.dat');

N = size(B,1);

Ngood = 0;
for i = 1:N,
    if((B(i,4) < 30) || (B(i,4) > 150))
        Ngood = Ngood + 1;
    end;
end;
fractionGood = Ngood/N;

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