clear;
load general.dat;
load increment.dat;
load initialParameters.dat;
N1=1;
N2=size(general,1);
%N2 = 20;
Nsaturation = min(N2,400);
Nbestincrement = min(Nsaturation, 400);

linearMagneticEnergy(1:N2) = 0;
normalMagneticEnergy(1:N2) = 0;
parallelMagneticEnergy(1:N2) = 0;
normalElectricEnergy(1:N2) = 0;
parallelElectricEnergy(1:N2) = 0;
deltaParallelMagneticEnergy(1:N2) = 0;
energyDeviation(1:N2)=0;
E0 = general(1,4);
Eb0 = general(1,7);

%omega_plasma = 2*3.14159*general(2,2)/general(2,3);

%omega_gyro_a = 1.982193107*10^8;

gamma = increment(1,1);

gamma = 1;
cv = initialParameters(10);
lorentzFactor = initialParameters(30);
omega_plasma = lorentzFactor*initialParameters(21);

for i = 1:N2,
    normalMagneticEnergy(i) = (general(i,27) + general(i,28))/E0;
    parallelMagneticEnergy(i) = general(i,26)/E0;
    normalElectricEnergy(i) = (general(i,24) + general(i,25))/E0;
    parallelElectricEnergy(i) = general(i,23)/E0;
    deltaParallelMagneticEnergy(i) = (general(i,26) - general(1,26))/E0;
    energyDeviation(i) = (general(i,7)-general(i,11))/E0;
end;


linearMagneticEnergy(1) = (normalMagneticEnergy(Nbestincrement))/exp(2*gamma*general(Nbestincrement,2));
for i = 2:N2,
    if(i < Nsaturation + 1)
        linearMagneticEnergy(i) = linearMagneticEnergy(1)*exp(2*gamma*general(i,2));
    else
        linearMagneticEnergy(i) = linearMagneticEnergy(Nsaturation);
    end
end

set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',16,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 2);

figure(1);
plot (general(2:N2,3)*gamma, normalMagneticEnergy(2:N2), 'green', general(1:N2,3)*gamma, parallelMagneticEnergy(1:N2), 'blue', general(1:N2, 3)*gamma, deltaParallelMagneticEnergy(1:N2), 'black', general(1:N2, 3)*gamma, normalElectricEnergy(1:N2), 'red', general(1:N2, 3)*gamma, parallelElectricEnergy(1:N2), 'yellow', general(1:Nsaturation, 2)*gamma, linearMagneticEnergy(1:Nsaturation),'cyan', general(1:N2,3)*gamma, energyDeviation(1:N2), 'magenta');

xlabel ('{{t {\gamma}_{max}}}');
ylabel ('E/E_{kin}');

legend('normal magnetic', 'parallel magnetic','delta parallel magnetic', 'normal electric', 'parallel electric','linear','full energy deviation','Location','northeast');
grid ;

