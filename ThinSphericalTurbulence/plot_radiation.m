clear;

radiation = importdata('radiation.dat');

N = size(radiation,1);
Nr = size(radiation,2);

jul12x(1:2) = 0;
jul12x(1) = 0.325;
jul12x(2) = 0.61;
jul12y(1:2) = 0;
jul12y(1) = 4.4;
jul12y(2) = 2.0;

jan12x(1:3) = 0;
jan12x(1) = 0.325;
jan12x(2) = 0.61;
jan12x(3) = 1.28;
jan12y(1:3) = 0;
jan12y(1) = 5.9;
jan12y(2) = 1.4;
jan12y(3) = 0.7;

apr11x(1:3) = 0;
apr11x(1) = 0.325;
apr11x(2) = 0.61;
apr11x(3) = 1.28;
apr11y(1:3) = 0;
apr11y(1) = 4.6;
apr11y(2) = 1.4;
ap11y(3) = 0.7;

decx(1:3) = 0;
decx(1) = 0.325;
decx(2) = 0.61;
decx(3) = 1.28;
decy(1:3) = 0;
decy(1) = 6.1;
decy(2) = 1.9;
decy(3) = 0.9;

octx(1:3) = 0;
octx(1) = 0.325;
octx(2) = 0.61;
octx(3) = 1.28;
octy(1:3) = 0;
octy(1) = 12;
octy(2) = 6.3;
octy(3) = 4.4;

augx(1:5) = 0;
augx(1) = 0.332;
augx(2) = 0.617;
augx(3) = 1.43;
augx(4) = 4.86;
augx(5) = 8.46;
augy(1:5) = 0;
augy(1) = 3.3;
augy(2) = 7.9;
augy(3) = 8.68;
augy(4) = 2.47;
augy(5) = 1.084;

augmax = 0.886;
augmaxy = 11.2;

junx(1:4) = 0;
junx(1) = 0.617;
junx(2) = 1.43;
junx(3) = 4.86;
junx(4) = 8.46;
juny(1:4) = 0;
juny(1) = 2.98;
juny(2) = 12.3;
juny(3) = 5.79;
juny(4) = 3.15;

junmaxx = 1.65;
junmaxy = 13.2;

mayx(1:3) = 0;
mayx(1) = 1.43;
mayx(2) = 4.86;
mayx(3) = 8.46;
mayy(1:3) = 0;
mayy(1) = 4.93;
mayy(2) = 12.2;
mayy(3) = 6.82;

maymaxx = 2.96;
maymaxy = 15.2;

aprx(1:4) = 0;
aprx(1) = 1.43;
aprx(2) = 4.86;
aprx(3) = 8.46;
aprx(4) = 22.5;
apry(1:4) = 0;
apry(1) = 1.3;
apry(2) = 12.86;
apry(3) = 17.57;
apry(4) = 5.2;

aprmaxx = 6.50;
aprmaxy = 19.3;

r = 5.17E16;
d = 40*3*10^24;
B = 0.346;
n = 2912;
me = 0.91*10^-27;
c = 3*10^10;
c1 = 6.27*10^18;
c5 = 7.25*10^-24;
c6 = 7.97*10^-41;
g = 3;
E0 = me*c*c;
N0 = n*(g-1)*(E0^(g-1));
f = 0.5;
s = 4*f*r/3;

h = 6.626*10^-27;

%theorRadiation(1:N) = 0;

%nu1 = 2*c1*((s*c6)^(2/(g+4)))*(N0^(2/(g+4)))*B^((g+2)/(g+4));
%for i = 1:N,
    %if (radiation(i,1)*(10^9) < nu1)
%        theorRadiation(i) = (10^26)*(3.14*r*r/d^2)*(c5/c6)*(B^-0.5)*((radiation(i,1)*(10^9)/(2*c1))^2.5)*(1 - exp(-((radiation(i,1)*(10^9)/nu1)^(-(g+4)/2))));
    %else 
    %    theorRadiation(i) = (10^26)*
    %end;
%end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
title ('I_{\nu}');
xlabel ('{\nu} GHz');
ylabel ('mJy');

loglog(radiation(1:N,1),radiation(1:N,2),'red','LineWidth',2);
loglog(radiation(1:N,1),radiation(1:N,3),'green','LineWidth',2);
loglog(radiation(1:N,1),radiation(1:N,4),'blue','LineWidth',2);
loglog(radiation(1:N,1),radiation(1:N,5),'magenta','LineWidth',2);
plot(radiation(1:N,1),radiation(1:N,6),'Color',[1.0,0.6,0],'LineWidth',2);
plot(radiation(1:N,1),radiation(1:N,7),'black','LineWidth',2);
%for i = 3:Nr,
 %   plot(radiation(1:N,1),radiation(1:N,i),'red');
%end;
%plot(radiation(1:N,1),theorRadiation(1:N),'green');
loglog(aprx(1:4),apry(1:4),'--o','Color','red','LineWidth',2);
loglog(mayx(1:3),mayy(1:3),'--o','Color','green','LineWidth',2);
loglog(junx(1:4),juny(1:4),'--o','Color','blue','LineWidth',2);
loglog(augx(1:5),augy(1:5),'--o','Color','magenta','LineWidth',2);
plot(octx(1:3),octy(1:3),'--o','Color',[1.0,0.6,0],'LineWidth',2);
plot(decx(1:3),decy(1:3),'--o','Color','black','LineWidth',2);

%legend('theory', 'shevalier', 'observation');

xlim([0.2 50]);
ylim([0.5 40]);

legend('april','may','june','august');

grid ;

dlmwrite('radiation.dat',radiation,'delimiter',' ');


gamma = - log(apry(3)/apry(4))/log(aprx(3)/aprx(4));

testx(1:3) = 0;
testx(1) = 8.46;
testx(2) = 22.5;
testx(3) = 1000*1.6*10^(-12)*10^(-9)/h;
testy(1:3) = 0;
testy(1) = apry(3);
testy(2) = testy(1)*power(testx(2)/testx(1),-gamma);
testy(3) = testy(1)*power(testx(3)/testx(1),-gamma);

figure(2);
hold on;
title ('I_{\nu}');
xlabel ('{\nu} GHz');
ylabel ('mJy');
loglog(radiation(1:N,1),radiation(1:N,2),'red','LineWidth',2);
loglog(aprx(1:4),apry(1:4),'--o','Color','red','LineWidth',2);