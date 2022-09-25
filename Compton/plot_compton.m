clear;


radiationEu = importdata('outputEu.dat');
radiationEd = importdata('outputEd.dat');

radiationu = importdata('outputNuu.dat');
radiationd = importdata('outputNud.dat');

Nnu = size(radiationu,1);




startPower = 35;
endPower = 50;

Fa(1:Nnu) = 0;

Fa(startPower) = radiationEd(startPower,2);
Fa(endPower) = radiationEd(endPower,2);

%gammap = log(Fpa(startPower)/Fpa(endPower))/log((me*energy(startPower)+m)/(me*energy(endPower)+m));

polyfitx(1:endPower-startPower + 1) = 0;
polyfity(1:endPower-startPower + 1) = 0;

for i = 1:endPower-startPower + 1,
    polyfitx(i) = log(radiationEd(i+startPower - 1,1));
    polyfity(i) = log(radiationEd(i+startPower - 1,2));
end;
p = polyfit(polyfitx, polyfity, 1);

%ap = exp(log(Fpa(startPower)) - gammap*log((me*energy(startPower)+m)));

for i = startPower-2:endPower+2,
    %Fpa(i) = ap*((me*energy(i)+m)^gammap);
    Fa(i) = exp(polyval(p, log(radiationEd(i))));
end;

figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\nu}');
xlabel ('\nu GHz');
ylabel ('F_{\nu} mJy');

relation = radiationEd(startPower,2)/radiationEu(startPower,2);

%plot(radiation1(1:Nnu,1),radiation1(1:Nnu,2),'red','LineWidth',2);
plot(radiationd(1:Nnu,1),radiationd(1:Nnu,2),'red','LineWidth',2);
%plot(radiation1(1:Nnu,1),Fa(1:Nnu),'blue','LineWidth',2);
%legend('Dubus','Uvarov/mc^2','\nu^{-1.5}');
grid ;

figure(3);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\nu}');
xlabel ('E eV');
ylabel ('E F_{E} erg/{s cm^{2}}');


%plot(radiation2(1:Nnu,1),radiation2(1:Nnu,2),'red','LineWidth',2);
plot(radiationEd(1:Nnu,1),radiationEd(1:Nnu,2),'red','LineWidth',2);
plot(radiationEd(1:Nnu,1),Fa(1:Nnu),'blue','LineWidth',2);
%legend('Dubus','Uvarov/mc^2');
%grid ;

dlmwrite('compton.dat',radiationEd,'delimiter',' ');