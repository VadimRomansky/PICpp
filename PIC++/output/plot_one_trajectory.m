clear;
load particlesTrajectories.dat;
N1=size(particlesTrajectories,1);
N2=(size(particlesTrajectories,2) - 2)/8;
n = 1;
figure(1);
plot(particlesTrajectories(1:N1,1),particlesTrajectories(1:N1,3 + (n-1)*8), 'red');
xlabel ('t s');
ylabel ('x cm');
grid ;

figure(2);
plot(particlesTrajectories(1:N1,1),particlesTrajectories(1:N1,4 + (n-1)*8), 'red');
xlabel ('t s');
ylabel ('y cm');
grid ;

figure(3);
plot(particlesTrajectories(1:N1,1),particlesTrajectories(1:N1,5 + (n-1)*8), 'red');
xlabel ('t s');
ylabel ('z cm');
grid ;

figure(4);
plot(particlesTrajectories(1:N1,1),particlesTrajectories(1:N1,6 + (n-1)*8), 'red');
xlabel ('t s');
ylabel ('px g*cm/s');
grid ;

figure(5);
plot(particlesTrajectories(1:N1,1),particlesTrajectories(1:N1,7 + (n-1)*8), 'red');
xlabel ('t s');
ylabel ('py g*cm/s');
grid ;

figure(6);
plot(particlesTrajectories(1:N1,1),particlesTrajectories(1:N1,8 + (n-1)*8), 'red');
xlabel ('t s');
ylabel ('pz g*cm/s');
grid ;

figure(7);
plot(particlesTrajectories(1:N1,3 + (n-1)*8),particlesTrajectories(1:N1,6 + (n-1)*8), 'red');
xlabel ('x cm');
ylabel ('px g*cm/s');
grid ;

figure(8);
plot(particlesTrajectories(1:N1,3 + (n-1)*8),particlesTrajectories(1:N1,4 + (n-1)*8), 'red');
xlabel ('x cm');
ylabel ('y cm');
grid ;
