clear;
radiation = importdata('radiation.dat');
%radiation1 = importdata('radiation1.dat');
%radiation3 = importdata('radiation3.dat');
%radiation4 = importdata('radiation4.dat');
%datar = importdata('at2018radio.dat');
N = size(radiation,1);

%Ndatar = size(datar,1);
%for i = 1:Ndatar,
%    datar(i,1) = datar(i,1)/10^9;
%    datar(i,2) = datar(i,2)/1000;
%end;
%from coppejans css161010
x99(1:4) = 0;
y99(1:4) = 0;
e99(1:4) = 0;
x99(1) = 1.5;
x99(2) = 3;
x99(3) = 6.1;
x99(4) = 9.87;
y99(1) = 1.5;
y99(2) = 4.3;
y99(3) = 6.1;
y99(4) = 4.2;
e99(1) = 0.1;
e99(2) = 0.2;
e99(3) = 0.3;
e99(4) = 0.2;

x162(1:5) = 0;
y162(1:5) = 0;
x162(1) = 1.5;
x162(2) = 2.94;
x162(3) = 6.1;
x162(4) = 9.74;
x162(5) = 22;
y162(1) = 4.7;
y162(2) = 2.9;
y162(3) = 2.3;
y162(4) = 1.74;
y162(5) = 0.56;

x357(1:4) = 0;
y357(1:4) = 0;
x357(1) = 0.33;
x357(2) = 0.61;
x357(3) = 1.5;
x357(4) = 3;
x357(5) = 6.05;
x357(6) = 10;
y357(1) = 0.357;
y357(2) = 0.79;
y357(3) = 0.27;
y357(4) = 0.17;
y357(5) = 0.07;
y357(6) = 0.032;

%at2018 from margutti
x15(1:4) = 0;
y15(1:4) = 0;

x15(1) = 35;
y15(1) = 8;
x15(2) = 225;
y15(2) = 30.8;
x15(3) = 233;
y15(3) = 28.6;
x15(4) = 241;
y15(4) = 27.4;

%at2018cow from nayana

x138(1:3) = 0;
y138(1:3) = 0;

x138(1) = 0.4;
x138(2) = 0.75;
x138(3) = 1.25;

y138(1) = 0.342;
y138(2) = 0.572;
y138(3) = 2.882;

x173(1:3) = 0;
y173(1:3) = 0;

x173(1) = 0.4;
x173(2) = 0.75;
x173(3) = 1.25;

y173(1) = 0.518;
y173(2) = 0.755;
y173(3) = 0.992;

x455(1:3) = 0;
y455(1:3) = 0;

x455(1) = 0.4;
x455(2) = 0.75;
x455(3) = 1.25;

y455(1) = 0.410;
y455(2) = 0.303;
y455(3) = 0.169;

x569(1:3) = 0;
y569(1:3) = 0;

x569(1) = 0.4;
x569(2) = 0.75;
x569(3) = 1.25;

y569(1) = 0.304;
y569(2) = 0.131;
y569(3) = 0.93;

%for ztf2018 z = 0.27
x350(1:4) = 0;
y350(1:4) = 0;
x350(1) = 1.5*(1 + 0.27);
x350(2) = 3*(1 + 0.27);
x350(3) = 6*(1 + 0.27);
x350(4) = 10*(1 + 0.27);

y350(1) = 0.146;
y350(2) = 0.068;
y350(3) = 0.089;
y350(4) = 0.045;

%for ztf20acigmel (AT2020xnd) z = 0.2433 ld = 1261 Mpc angular distance = 816 Mpc\
x40(1:5) = 0;
y40(1:5) = 0;
x40(1) = 10*(1 + 0.2433);
x40(2) = 33*(1 + 0.2433);
x40(3) = 45*(1 + 0.2433);
x40(4) = 79*(1 + 0.2433);
x40(5) = 94*(1 + 0.2433);

y40(1) = 0.079;
y40(2) = 0.497;
y40(3) = 0.675;
y40(4) = 0.912;
y40(5) = 0.825;

x71(1:5) = 0;
y71(1:5) = 0;
x71(1) = 10*(1 + 0.2433);
x71(2) = 15*(1 + 0.2433);
x71(3) = 22*(1 + 0.2433);
x71(4) = 33*(1 + 0.2433);
x71(5) = 45*(1 + 0.2433);

y71(1) = 0.18;
y71(2) = 0.4;
y71(3) = 0.484;
y71(4) = 0.45;
y71(5) = 0.209;

x95(1:4) = 0;
y95(1:4) = 0;
x95(1) = 10*(1 + 0.2433);
x95(2) = 15*(1 + 0.2433);
x95(3) = 22*(1 + 0.2433);
x95(4) = 33*(1 + 0.2433);


y95(1) = 0.168;
y95(2) = 0.278;
y95(3) = 0.301;
y95(4) = 0.213;


figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
plot(radiation(1:N,1), radiation(1:N,2),'Color','red','LineWidth',2);
%plot(radiation1(1:N,1),radiation1(1:N,2),'Color','red','LineWidth',2);
%plot(radiation3(1:N,1),radiation3(1:N,2),'Color','blue','LineWidth',2);
%plot(radiation4(1:N,1),radiation4(1:N,2),'Color','green','LineWidth',2);
sz = 20;
%scatter(x357(1:6), y357(1:6),sz,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','LineWidth',2);
errorbar(x99(1:4), y99(1:4),e99(1:4),'red','LineWidth',2,'LineStyle','--');
%scatter(datar(1:Ndatar,1), datar(1:Ndatar,2),sz,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','LineWidth',2);
%legend('F(E) ~ E^{-3.5}','PIC spectrum','observation','Location','southeast');
title ('F_{\nu}');
xlabel ('{\nu} GHz');
ylabel ('mJy');

nufnu(1:N)=0;
for i = 1:N,
    nufnu(i) = radiation(i,2)*radiation(i,1)/10^17;
end;

% figure(2);
% hold on;
% set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');
% plot(radiation(1:N,1), nufnu(1:N),'Color','red','LineWidth',2);
% title ('{\nu} F_{\nu}');
% xlabel ('{\nu} GHz');
% ylabel ('{\nu} F_{\nu}');

% output(1:N,1:4) = 0;
% for i = 1:N,
%     output(i,1) = log10(radiation1(i,1));
%     output(i,2) = radiation3(i,2);
%     output(i,3) = radiation1(i,2);
%     output(i,4) = radiation4(i,2);
% end;
% output2(1:4, 1:2) = 0;
% output3(1:4, 1:2) = 0;
% for i = 1:4,
%     output2(i,1) = log10(x99(i));
%     output2(i,2) = y99(i);
%     output3(i,1) = log10(x99(i));
%     output3(i,2) = e99(i);
% end;
% 
% dlmwrite('css16.dat',output,'delimiter',' ');
% dlmwrite('css16obs.dat',output2,'delimiter',' ');
% dlmwrite('css16err.dat',output3,'delimiter',' ');