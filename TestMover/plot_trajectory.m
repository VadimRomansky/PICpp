clear;
trajectory = importdata('./output/trajectory_1.dat');
Bzfield = importdata('./output/outBz.dat');
Nx = size(Bzfield,1);
Ny = size(Bzfield,2);

mp = 1.67262177E-24;
me = mp/100;
electron_charge = 4.803529695E-10;
c = 2.99792458E10;
n = 1;
ntristan = 2;
sigma = 4.0;
gamma = 1.5;
ctristan = 0.45;
comp = 5;
omp = ctristan/comp;
qtristan = omp*omp*gamma/(ntristan*(1 + me/mp));
metristan = qtristan;
fieldScale = sqrt(4*3.14*(n/ntristan)*(me/metristan)*(c*c/(ctristan*ctristan)));

sampling = 20;
omega_pe = sqrt(4*pi*n*electron_charge*electron_charge/(gamma*me));
dx = 0.2*c/omega_pe;

Nt = size(trajectory,1);
downstreamNx = 300;

E(1:Nt) = 0;
for i = 1:Nt,
    E(i) = sqrt(1.0 + (trajectory(i,5)*trajectory(i,5) + trajectory(i,6)*trajectory(i,6) + trajectory(i,7)*trajectory(i,7))/(me*me*c*c));
end;

figure(1);
plot (trajectory(1:Nt,2),trajectory(1:Nt,1), 'red');
title ('xt');
xlabel ('x');
ylabel ('t');
grid ;

figure(2);
hold on;
%caxis ([0 2.5])
fig = imagesc((1:Ny)*dx*sampling, ((1:Nx) - downstreamNx*3)*dx*sampling, Bzfield);
plot(trajectory(1,3), trajectory(1,2), 'ro', 'MarkerSize', 10,'Color','red');
plot (trajectory(1:Nt,3),trajectory(1:Nt,2), 'red');
title ('xy');
xlabel ('y');
ylabel ('x');
grid ;

figure(3);
hold on;
plot(trajectory(1,3), trajectory(1,4), 'ro', 'MarkerSize', 10,'Color','red');
plot (trajectory(1:Nt,3),trajectory(1:Nt,4), 'red');
title ('yz');
xlabel ('y');
ylabel ('z');
grid ;

figure(4);
plot (trajectory(1:Nt,1),trajectory(1:Nt,5), 'red');
title ('pt');
xlabel ('t');
ylabel ('px');
grid ;

figure(5);
plot (trajectory(1:Nt,1),E(1:Nt), 'red');
title ('gamma');
xlabel ('t');
ylabel ('gamma');
grid ;