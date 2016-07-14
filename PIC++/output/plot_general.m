clear;
load general.dat;
N1=1;
N2=size(general,1);

Nsaturation = N2;
%Nsaturation = 70;
N3 = Nsaturation + 20;
linearField(1:N2) = 0;
linearEnergy(1:N2) = 0;
beta = 0.2;
%gamma = beta*2*3.14*sqrt(sqrt(1 - beta*beta));
%gamma = beta*2*3.14;
gamma = 0.01*1.98*10^8;

linearField(1) = general(Nsaturation, 16)/exp(gamma*general(Nsaturation,2));
linearEnergy(1) = general(Nsaturation, 6)/exp(2*gamma*general(Nsaturation,2));
for i = 2:N2,
    linearField(i) = linearField(1)*exp(gamma*general(i,2));
    linearEnergy(i) = linearEnergy(1)*exp(2*gamma*general(i,2));
end
%for i = Nsaturation+1:N2,
    %linearField(i) = linearField(Nsaturation);
%end


figure(1);
plot (general(1:N2,2), general(1:N2,4), 'green', general(1:N2,2), general(1:N2,5), 'blue', general(1:N2, 2), general(1:N2,6), 'black', general(1:N2, 2), general(1:N2,7), 'red', general(1:N2, 2), general(1:N2,11), 'yellow');
%plot (general(1:N3,2), linearEnergy(1:N3), 'blue', general(1:N2, 2), general(1:N2,6), 'green', general(1:N2, 2), general(1:N2,7), 'red');
title ('energy');
xlabel ('t');
ylabel ('E erg');
legend('particle', 'electric','magnetic', 'full', 'theoretical','Location','southwest');
%legend('linear increment','magnetic', 'full','Location','southwest');
grid ;

figure(2);
plot (general(1:N2,2), general(1:N2,8), 'red', general(1:N2,2), general(1:N2,9), 'green', general(1:N2, 2), general(1:N2,10), 'blue');
title ('momentum');
xlabel ('{{t w_p}/{2\pi}}');
ylabel ('P g*cm/s');
legend('x', 'y','z','Location','northeast');
grid ;

figure(3);
plot (general(1:N2,2), general(1:N2,16), 'red', general(1:N2,2), linearField(1:N2),'blue');
%plot (general(1:N2,2), general(1:N2,15), 'red', general(1:N2,2), general(1:N2,16), 'green');
title ('max field');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('gauss');
%legend('electric field', 'magnetic field','Location','southeast');
legend('magnetic field', 'magnetic field with linear increment {{u}/{c} {\omega_p}}','Location','northwest');
grid ;

figure(4);
plot (general(1:N2,2), general(1:N2,17), 'red');
title ('dt');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('dt');
grid ;