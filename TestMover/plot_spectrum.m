clear;
distribution = importdata('./output/distribution_13.dat');

Np = size(distribution,1);
mp = 1.67262177E-24;
me = mp/100;
c = 2.99792458E10;

F(1:Np) = 0;
for i = 1:Np,
    F(i) = distribution(i,2)*distribution(i,1)*distribution(i,1);
end;

figure(1);
plot (distribution(1:Np,1)/(me*c), F(1:Np), 'red');
title ('F');
xlabel ('p/mc');
ylabel ('F*p^4');
grid ;