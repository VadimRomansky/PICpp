clear;
load general.dat;
load increment.dat;
load initialParameters.dat;
N1=1;
N2=size(general,1);

Nsaturation = min(N2,120);
Nbestincrement = min(Nsaturation, 50);
linearField(1:N2) = 0;
linearMagneticEnergy(1:N2) = 0;

%omega_plasma = 2*3.14159*general(2,2)/general(2,3);

%omega_gyro_a = 1.982193107*10^8;

gamma = increment(1,1);
shockWaveV(1:N2) = 0;
shockWaveV(1) = 0;
%gamma = 0;
cv = initialParameters(10);
omega_plasma = initialParameters(20);


linearField(1) = general(Nbestincrement, 16)/exp(gamma*general(Nbestincrement,2));
linearMagneticEnergy(1) = (general(Nbestincrement, 6)/exp(2*gamma*general(Nbestincrement,2)));
for i = 2:N2,
    %shockWaveV(i) = (general(i,19) - general(2,19))/general(i,2);
    if(i < Nsaturation + 1)
        linearField(i) = linearField(1)*exp(gamma*general(i,2));
        linearMagneticEnergy(i) = linearMagneticEnergy(1)*exp(2*gamma*general(i,2));
    else
        linearField(i) = linearField(Nsaturation);
        linearMagneticEnergy(i) = linearMagneticEnergy(Nsaturation);
    end
end
shockWave(2) = 0;

set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',16,'DefaultTextFontName','Times New Roman'); 
figure(1);
plot (general(1:N2,3)*omega_plasma, general(1:N2,4), 'green', general(1:N2,3)*omega_plasma, general(1:N2,5), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,6), 'black', general(1:N2, 3)*omega_plasma, general(1:N2,7), 'red', general(1:N2, 3)*omega_plasma, general(1:N2,11), 'yellow');
%plot (2*pi*general(1:N2, 2), general(1:N2,6)/general(1,7), 'red',  2*pi*general(1:Nsaturation, 2), linearMagneticEnergy(1:Nsaturation)/general(1,7),'blue');
%plot (2*pi*general(1:N2, 2), general(1:N2,5)/general(1,7), 'red',  2*pi*general(1:Nsaturation, 2), linearMagneticEnergy(1:Nsaturation)/general(1,7),'blue');
%title ('energy fraction');
xlabel ('{{t {\omega}_p}}');
ylabel ('E/E_0');
%legend('fraction of magnetic energy', 'linear approximation','Location','southeast');
legend('particle', 'electric','magnetic', 'full', 'theoretical','Location','southwest');
grid ;

figure(2);
plot (general(1:N2,3)*omega_plasma, general(1:N2,8), 'red', general(1:N2,3)*omega_plasma, general(1:N2,9), 'green', general(1:N2, 3)*omega_plasma, general(1:N2,10), 'blue');
title ('momentum');
xlabel ('{{t w_p}/{2\pi}}');
ylabel ('P g*cm/s');
legend('x', 'y','z','Location','northeast');
grid ;

figure(3);
plot (general(1:N2,3)*omega_plasma, general(1:N2,16), 'red', general(1:N2,3)*omega_plasma, linearField(1:N2),'blue');
%plot (general(1:N2,2), general(1:N2,15), 'red', general(1:N2,2), general(1:N2,16), 'green');
title ('max field');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('gauss');
%legend('electric field', 'magnetic field','Location','southeast');
legend('magnetic field', 'magnetic field with linear increment {{u}/{c} {\omega_p}}','Location','northwest');
grid ;

figure(4);
plot (general(1:N2,3)*omega_plasma, general(1:N2,17), 'red');
title ('dt');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('dt');
grid ;

figure(5);
plot (general(1:N2,3)*omega_plasma, general(1:N2,19), 'red');
title ('shock wave x');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('dt');
grid ;

figure(6);
plot (general(1:N2,3)*omega_plasma,shockWaveV(1:N2), 'red');
title ('shock wave V');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('dt');
grid ;

figure(7);
plot (general(1:N2,3)*omega_plasma,general(1:N2, 20), 'red');
title ('mean E');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('mean Ex');
grid ;

figure(8);
plot (general(1:N2,3)*omega_plasma,general(1:N2, 21), 'red');
title ('mean E');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('mean Ey');
grid ;

figure(9);
plot (general(1:N2,3)*omega_plasma,general(1:N2, 22), 'red');
title ('mean E');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('mean Ez');
grid ;