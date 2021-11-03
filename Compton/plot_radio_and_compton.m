clear;

radio = importdata('radiation.dat');
compton = importdata('output1.dat');
datar = importdata('at2018radio.dat');
datac = importdata('at2018xray.dat');

q1 = 1.6*10^-9;
h = 6.26*10^-26;

Ndatac = size(datac,1);
Ndatar = size(datar,1);

Nnuc = size(compton,1);
Nnur = size(radio,1);

for i = 1:Ndatac,
    datac(i,1) = (datac(i,1)*q1/h)/10^9;
    datac(i,2) = (datac(i,2)*h*10^26);
end;

for i = 1:Ndatar,
    datar(i,1) = datar(i,1)/10^9;
    datar(i,2) = datar(i,2)/1000;
end;

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

%AT2018 at 7.7 days
xradio7(1:4) = 0;
yradio7(1:4) = 0;

xradio7(1) = 35;
yradio7(1) = 8;
xradio7(2) = 225;
yradio7(2) = 30.8;
xradio7(3) = 233;
yradio7(3) = 28.6;
xradio7(4) = 241;
yradio7(4) = 27.4;

xcompton7(1:4) = 0;
ycompton7(1:4) = 0;

xcompton7(1) = 2.9*10^9;
ycompton7(1) = 0.000144;
xcompton7(2) = 4.3*10^9;
ycompton7(2) = 0.000185;
xcompton7(3) = 6.7*10^9;
ycompton7(3) = 0.000304;
xcompton7(4) = 8.4*10^9;
ycompton7(4) = 0.00039;
xcompton7(5) = 1.12*10^10;
ycompton7(5) = 0.00041;
xcompton7(6) = 1.55*10^10;
ycompton7(6) = 0.0003;

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\nu}');
xlabel ('\nu GHz');
ylabel ('F_{\nu} mJy');


plot(radio(1:Nnur,1),radio(1:Nnur,2),'red','LineWidth',2);
plot(compton(1:Nnuc,1),compton(1:Nnuc,2),'blue','LineWidth',2);
%plot(xradio15(1:4),yradio15(1:4),'--o','Color','magenta','LineWidth',2);
%plot(xcompton15(1:4),ycompton15(1:4),'--o','Color','green','LineWidth',2);
%plot(datac(1:Ndata,1),
%datac(1:Ndata,2),'--o','Color','green','LineWidth',2);
sz = 20;
scatter(datar(1:Ndatar,1), datar(1:Ndatar,2),sz,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','LineWidth',2);
scatter(datac(1:Ndatac,1), datac(1:Ndatac,2),sz,'MarkerEdgeColor','green','MarkerFaceColor','green','LineWidth',2);
grid ;

legend('synchrotron', 'inverse compton', 'radio observations', 'X-ray observations');
