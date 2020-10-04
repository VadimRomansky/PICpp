clear;
load general.dat;
load increment.dat;
load initialParameters.dat;
N1=1;
N2=size(general,1);
%N2 = 20;
Nsaturation = min(N2,100);
Nbestincrement = min(Nsaturation, 20);
linearField(1:N2) = 0;
linearMagneticEnergy(1:N2) = 0;

%omega_plasma = 2*3.14159*general(2,2)/general(2,3);

%omega_gyro_a = 1.982193107*10^8;

gamma = increment(1,1);
shockWaveV(1:N2) = 0;
shockWaveV(1) = 0;
%gamma = 0;
cv = initialParameters(10);
lorentzFactor = initialParameters(30);
omega_plasma = initialParameters(21);


linearField(1) = general(Nbestincrement, 16)/exp(gamma*general(Nbestincrement,2));
linearMagneticEnergy(1) = ((general(Nbestincrement, 27) + general(Nbestincrement,28))/exp(2*gamma*general(Nbestincrement,2)));
for i = 2:N2,
    %shockWaveV(i) = (general(i,19) - general(2,19))/general(i,2);
    if(i < Nsaturation + 1)
        linearField(i) = linearField(1)*exp(gamma*general(i,2));
        linearMagneticEnergy(i) = linearMagneticEnergy(1)*exp(2*gamma*general(i,2));
    else
        linearField(i) = linearField(Nsaturation);
        linearMagneticEnergy(i) = linearMagneticEnergy(Nsaturation);
    end
    linearField(i) = 0;
    linearMagneticEnergy(i) = 0;
end
shockWave(2) = 0;

set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',16,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 1);
figure(1);
%plot (general(1:N2,3)*omega_plasma, general(1:N2,4), 'green', general(1:N2,3)*omega_plasma, general(1:N2,5), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,6), 'black', general(1:N2, 3)*omega_plasma, general(1:N2,7), 'red', general(1:N2, 3)*omega_plasma, general(1:N2,11), 'yellow');
plot (general(1:N2,3)*omega_plasma, general(1:N2,4), 'green', general(1:N2,3)*omega_plasma, general(1:N2,5), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,6), 'black', general(1:N2, 3)*omega_plasma, general(1:N2,7), 'red', general(1:N2, 3)*omega_plasma, general(1:N2,11), 'yellow', general(1:Nsaturation, 3)*omega_plasma, linearMagneticEnergy(1:Nsaturation),'cyan');

xlabel ('{{t {\omega}_p}}');
ylabel ('E/E_0');

%legend('particle', 'electric','magnetic', 'full', 'theoretical','Location','southwest');
legend('particle', 'electric','magnetic', 'full', 'theoretical','linear','Location','southwest');
grid ;

figure(2);
%plot (general(1:N2,3)*omega_plasma, general(1:N2,4), 'green', general(1:N2,3)*omega_plasma, general(1:N2,5), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,6), 'black', general(1:N2, 3)*omega_plasma, general(1:N2,7), 'red', general(1:N2, 3)*omega_plasma, general(1:N2,11), 'yellow');
plot (general(1:N2,3)*omega_plasma, general(1:N2,4), 'green', general(1:N2,3)*omega_plasma, general(1:N2,23), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,26), 'black', general(1:N2, 3)*omega_plasma, general(1:N2,7), 'red', general(1:N2, 3)*omega_plasma, general(1:N2,11), 'yellow', general(1:N2, 3)*omega_plasma, general(1:N2, 24) + general(1:N2, 25),'cyan',general(1:N2, 3)*omega_plasma, general(1:N2, 27) + general(1:N2,28),'magenta');

xlabel ('{{t {\omega}_p}}');
ylabel ('E/E_0');

%legend('particle', 'electric','magnetic', 'full', 'theoretical','Location','southwest');
legend('particle', 'electric parallel','magnetic parallel', 'full', 'theoretical','electric normal', 'magnetic normal','Location','southwest');
grid ;

figure(3);
plot (general(1:N2,3)*omega_plasma, general(1:N2,8), 'red', general(1:N2,3)*omega_plasma, general(1:N2,29), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,32), 'green', general(1:N2, 3)*omega_plasma, general(1:N2,12), 'yellow');
title ('momentum x');
xlabel ('{{t w_p}/{2\pi}}');
ylabel ('P g*cm/s');
legend('total', 'electromagnetic','particle' ,'theoretical','Location','northeast');
grid ;

figure(4);
plot (general(1:N2,3)*omega_plasma, general(1:N2,9), 'red', general(1:N2,3)*omega_plasma, general(1:N2,30), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,33), 'green', general(1:N2, 3)*omega_plasma, general(1:N2,13), 'yellow');
title ('momentum y');
xlabel ('{{t w_p}/{2\pi}}');
ylabel ('P g*cm/s');
legend('total', 'electromagnetic','particle' ,'theoretical','Location','northeast');
grid ;

figure(5);
plot (general(1:N2,3)*omega_plasma, general(1:N2,10), 'red', general(1:N2,3)*omega_plasma, general(1:N2,31), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,34), 'green', general(1:N2, 3)*omega_plasma, general(1:N2,14), 'yellow');
title ('momentum z');
xlabel ('{{t w_p}/{2\pi}}');
ylabel ('P g*cm/s');
legend('total', 'electromagnetic','particle' ,'theoretical','Location','northeast');
grid ;

figure(6);
plot (general(1:N2,3)*omega_plasma, general(1:N2,16), 'red', general(1:N2,3)*omega_plasma, linearField(1:N2),'blue');
%plot (general(1:N2,2), general(1:N2,15), 'red', general(1:N2,2), general(1:N2,16), 'green');
title ('max field');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('gauss');
%legend('electric field', 'magnetic field','Location','southeast');
legend('magnetic field', 'magnetic field with linear increment {{u}/{c} {\omega_p}}','Location','northwest');
grid ;

figure(7);
plot (general(1:N2,3)*omega_plasma, general(1:N2,17), 'red');
title ('dt');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('dt');
grid ;

figure(8);
plot (general(1:N2,3)*omega_plasma, general(1:N2,19), 'red');
title ('shock wave x');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('dt');
grid ;

figure(9);
plot (general(1:N2,3)*omega_plasma,shockWaveV(1:N2), 'red');
title ('shock wave V');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('dt');
grid ;

figure(10);
plot (general(1:N2,3)*omega_plasma,general(1:N2, 20), 'red');
title ('mean squared E');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('mean squared Ex');
grid ;

figure(11);
plot (general(1:N2,3)*omega_plasma,general(1:N2, 21), 'red');
title ('mean squared E');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('mean squared Ey');
grid ;

figure(12);
plot (general(1:N2,3)*omega_plasma,general(1:N2, 22), 'red');
title ('mean squared E');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('mean squared Ez');
grid ;

figure(13);
%plot (general(1:N2,3)*omega_plasma, general(1:N2,4), 'green', general(1:N2,3)*omega_plasma, general(1:N2,5), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,6), 'black', general(1:N2, 3)*omega_plasma, general(1:N2,7), 'red', general(1:N2, 3)*omega_plasma, general(1:N2,11), 'yellow');
plot (general(1:N2, 3)*omega_plasma, general(1:N2,27) + general(1:N2,28), 'black', general(1:Nsaturation, 2)*omega_plasma, linearMagneticEnergy(1:Nsaturation),'cyan');

xlabel ('{{t {\omega}_p}}');
ylabel ('E/E_0');

%legend('particle', 'electric','magnetic', 'full', 'theoretical','Location','southwest');
legend('magnetic','linear','Location','southwest');
grid ;

figure(14);
plot (general(1:N2,3)*omega_plasma, general(1:N2,36), 'red');
title ('der Ex point');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('x cm');
grid ;

figure(15);
plot (general(1:N2,3)*omega_plasma, general(1:N2,38), 'red');
title ('const mean level point');
xlabel ('{{t \omega_p}/{2\pi}}');
ylabel ('x cm');
grid ;

%figure(16);
%plot (general(1:N2,3)*omega_plasma,general(1:N2, 39), 'red');
%title ('mean squared B');
%xlabel ('{{t \omega_p}/{2\pi}}');
%ylabel ('mean squared Bx');
%grid ;

%figure(17);
%plot (general(1:N2,3)*omega_plasma,general(1:N2, 40), 'red');
%title ('mean squared B');
%xlabel ('{{t \omega_p}/{2\pi}}');
%ylabel ('mean squared By');
%grid ;

%figure(18);
%plot (general(1:N2,3)*omega_plasma,general(1:N2, 41), 'red');
%title ('mean squared B');
%xlabel ('{{t \omega_p}/{2\pi}}');
%ylabel ('mean squared Bz');
%grid ;