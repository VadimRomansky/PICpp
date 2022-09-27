clear;
filenumber = 3;
energy = importdata(strcat('Ee',num2str(filenumber),'.dat'));
distribution = importdata(strcat('Fs',num2str(filenumber),'.dat'));
N1 = size(energy,2);
mp = 100;
me = 1;
for i = 1:N1,
    %energy(i) = energy(i)*me;
    energy(i) = energy(i)*me*sqrt(1800/100);
end;

distribution2 = importdata('GLE_pdf_sf_B1.dat');
N2 = size(distribution2,1);
distribution3(1:N2,2)=0;
for i = 1:N2,
    u = (mp/me)*10^distribution2(i,1);
    gamma = sqrt(u*u + 1);
    distribution3(i,1) = gamma*mp;
    distribution3(i,2) = distribution2(i,2)*gamma/(u*u*sqrt(gamma*gamma - 1));
end;

N3 = N1 + N2-36;
distribution4(1:N3,2) = 0;
for i = 1:N3;
    if (i < 140)
        distribution4(i,1) = energy(i)+1;
        distribution4(i,2) = distribution(i);
    elseif (i < N1)
        distribution4(i,1) = energy(i)+1;
        distribution4(i,2) = distribution(140)*power((energy(i)+1)/(energy(140)+1), -3.5);
    elseif (i < N1 + 10)
        distribution4(i,1) = distribution3(i-N1+36,1);
        distribution4(i,2) = distribution(140)*power(distribution4(i,1)/(energy(140)+1), -3.5);
    else
        distribution4(i,1) = distribution3(i-N1+36,1);
        distribution4(i,2) = distribution3(i-N1+36,2)*(distribution4(N1 + 10-1,2))/distribution3(36+10-1,2);
    end
end;



figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\gamma}');
xlabel ('\gamma');
ylabel ('F_{\gamma}');

plot(energy(1:N1)+1,distribution(1:N1),'blue','LineWidth',2);
plot(distribution3(1:N2,1), distribution3(1:N2,2), 'red', 'LineWidth',2);

figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\gamma}');
xlabel ('\gamma');
ylabel ('F_{\gamma}');

plot(distribution4(1:N3,1), distribution4(1:N3,2), 'red', 'LineWidth',2);

energyOutput(1:N3,1) = 0;
distributionOutput(1:N3) = 0;
dlmwrite(strcat('Em',num2str(filenumber),'.dat'),energyOutput, 'delimiter',' ');
dlmwrite(strcat('Fm',num2str(filenumber),'.dat'),distributionOutput, 'delimiter',' ');
