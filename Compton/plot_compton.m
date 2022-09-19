clear;


radiation1 = importdata('output1.dat');
radiation3 = importdata('output3.dat');

radiation2 = importdata('output2.dat');
radiation4 = importdata('output4.dat');

Nnu = size(radiation1,1);




startPower = 100;
endPower = 110;

Fa(1:Nnu) = 0;

Fa(startPower) = radiation1(startPower,2);
Fa(endPower) = radiation1(endPower,2);

%gammap = log(Fpa(startPower)/Fpa(endPower))/log((me*energy(startPower)+m)/(me*energy(endPower)+m));

polyfitx(1:endPower-startPower + 1) = 0;
polyfity(1:endPower-startPower + 1) = 0;

for i = 1:endPower-startPower + 1,
    polyfitx(i) = log(radiation1(i+startPower - 1,1));
    polyfity(i) = log(radiation1(i+startPower - 1,2));
end;
p = polyfit(polyfitx, polyfity, 1);

%ap = exp(log(Fpa(startPower)) - gammap*log((me*energy(startPower)+m)));

for i = startPower-2:endPower+2,
    %Fpa(i) = ap*((me*energy(i)+m)^gammap);
    Fa(i) = exp(polyval(p, log(radiation1(i))));
end;

figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\nu}');
xlabel ('\nu GHz');
ylabel ('F_{\nu} mJy');

relation = radiation1(startPower,2)/radiation3(startPower,2);

plot(radiation1(1:Nnu,1),radiation1(1:Nnu,2),'red','LineWidth',2);
plot(radiation3(1:Nnu,1),radiation3(1:Nnu,2),'green','LineWidth',2);
plot(radiation1(1:Nnu,1),Fa(1:Nnu),'blue','LineWidth',2);
legend('Dubus','Uvarov/mc^2','\nu^{-1.5}');
grid ;

figure(3);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\nu}');
xlabel ('E eV');
ylabel ('E*F_{E} erg s^{-1}');



plot(radiation2(1:Nnu,1),radiation2(1:Nnu,2),'red','LineWidth',2);
plot(radiation4(1:Nnu,1),radiation4(1:Nnu,2),'green','LineWidth',2);
legend('Dubus','Uvarov/mc^2');
grid ;
