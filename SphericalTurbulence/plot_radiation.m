clear;

load radiation.dat;

N = size(radiation,1);

augx(1:5) = 0;
augx(1) = 0.335;
augx(2) = 0.625;
augx(3) = 1.46;
augx(4) = 4.92;
augx(5) = 8.57;
augy(1:5) = 0;
augy(1) = 3.29;
augy(2) = 7.77;
augy(3) = 8.53;
augy(4) = 2.42;
augy(5) = 1.06;

augmax = 0.886;
augmaxy = 11.2;

junx(1:4) = 0;
junx(1) = 0.628;
junx(2) = 1.45;
junx(3) = 4.89;
junx(4) = 8.53;
juny(1:4) = 0;
juny(1) = 2.98;
juny(2) = 12.3;
juny(3) = 5.79;
juny(4) = 3.15;

junmaxx = 1.65;
junmaxy = 13.2;

mayx(1:3) = 0;
mayx(1) = 1.46;
mayx(2) = 4.94;
mayx(3) = 8.62;
mayy(1:3) = 0;
mayy(1) = 4.91;
mayy(2) = 12.0;
mayy(3) = 6.67;

maymaxx = 2.96;
maymaxy = 15.2;

aprx(1:4) = 0;
aprx(1) = 1.44;
aprx(2) = 4.91;
aprx(3) = 8.58;
aprx(4) = 22.8;
apry(1:4) = 0;
apry(1) = 0.993;
apry(2) = 13.9;
apry(3) = 17.1;
apry(4) = 5.11;

aprmaxx = 6.50;
aprmaxy = 19.3;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
title ('I_{\nu}');
xlabel ('{\nu} GHz');
ylabel ('mJy');

plot(radiation(1:N,1),radiation(1:N,2),'red');
plot(augx(1:5),augy(1:5),'--o');

legend('August theory', 'August observation');

grid ;