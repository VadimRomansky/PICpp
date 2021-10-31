clear;

radio = importdata('radiation.dat');
compton = importdata('output1.dat');

Nnuc = size(compton,1);
Nnur = size(radio,1);

%AT2018 at 16.5 days
xradio15(1:4) = 0;
yradio15(1:4) = 0;

xradio15(1) = 35;
yradio15(1) = 8;
xradio15(2) = 225;
yradio15(2) = 30.8;
xradio15(3) = 233;
yradio15(3) = 28.6;
xradio15(4) = 241;
yradio15(4) = 27.4;

xcompton15(1:4) = 0;
ycompton15(1:4) = 0;

xcompton15(1) = 1.48*10^8;
ycompton15(1) = 0.00074;
xcompton15(2) = 5.23*10^8;
ycompton15(2) = 0.00055;
xcompton15(3) = 2.04*10^9;
ycompton15(3) = 0.00033;
xcompton15(4) = 7.23*10^9;
ycompton15(4) = 0.00018;

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\nu}');
xlabel ('\nu GHz');
ylabel ('F_{\nu} mJy');


plot(radio(1:Nnur,1),radio(1:Nnur,2),'red','LineWidth',2);
plot(compton(1:Nnuc,1),compton(1:Nnuc,2),'blue','LineWidth',2);
plot(xradio15(1:4),yradio15(1:4),'--o','Color','magenta','LineWidth',2);
plot(xcompton15(1:4),ycompton15(1:4),'--o','Color','green','LineWidth',2);
grid ;

legend('synchrotron', 'inverse compton', 'radio observations', 'X-ray observations');
