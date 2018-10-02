clear;
load general1.dat;
load general2.dat;
load general3.dat;
load general4.dat;
load general5.dat;

N1=size(general1,1);
N2=size(general2,1);
N3=size(general3,1);
N4=size(general4,1);
N5=size(general5,1);


set(0, 'DefaultLineLineWidth', 2);

figure(1);
plot (2*pi*general1(1:N1, 2), general1(1:N1,6)/general1(1,7), 'red', 2*pi*general2(1:N2, 2), general2(1:N2,6)/general2(1,7), 'blue',2*pi*general3(1:N3, 2), general3(1:N3,6)/general3(1,7), 'black',2*pi*general4(1:N4, 2), general4(1:N4,6)/general4(1,7), 'green',2*pi*general5(1:N5, 2), general5(1:N5,6)/general5(1,7), 'yellow');
title ('energy');
xlabel ('{{t w_p}/{2\pi}}');
ylabel ('E/E0');
grid ;

figure(2);
plot (2*pi*general1(1:N1, 2), general1(1:N1,6), 'red', 2*pi*general2(1:N2, 2), general2(1:N2,6), 'blue',2*pi*general3(1:N3, 2), general3(1:N3,6), 'black',2*pi*general4(1:N4, 2), general4(1:N4,6), 'green',2*pi*general5(1:N5, 2), general5(1:N5,6), 'yellow');
title ('energy');
xlabel ('{{t w_p}/{2\pi}}');
ylabel ('E erg');
grid ;
