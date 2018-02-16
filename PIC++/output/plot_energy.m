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
omega_plasma = lorentzFactor*initialParameters(21);


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
end
shockWave(2) = 0;

set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',16,'DefaultTextFontName','Times New Roman'); 
figure(1);
%plot (general(1:N2,3)*omega_plasma, general(1:N2,4), 'green', general(1:N2,3)*omega_plasma, general(1:N2,5), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,6), 'black', general(1:N2, 3)*omega_plasma, general(1:N2,7), 'red', general(1:N2, 3)*omega_plasma, general(1:N2,11), 'yellow');
plot (general(1:N2,3)*omega_plasma, general(1:N2,4), 'green', general(1:N2,3)*omega_plasma, general(1:N2,5), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,27) + general(1:N2,28), 'black', general(1:N2, 3)*omega_plasma, general(1:N2,7), 'red', general(1:N2, 3)*omega_plasma, general(1:N2,11), 'yellow', general(1:Nsaturation, 3)*omega_plasma, linearMagneticEnergy(1:Nsaturation),'cyan');

xlabel ('{{t {\omega}_p}}');
ylabel ('E/E_0');

%legend('particle', 'electric','magnetic', 'full', 'theoretical','Location','southwest');
legend('particle', 'electric','magnetic', 'full', 'theoretical','linear','Location','southwest');
grid ;