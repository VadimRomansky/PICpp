clear;
load general.dat;
N1=1;
N2=size(general,1);


figure(1);
plot (general(1:N2,1), general(1:N2,3), 'green', general(1:N2,1), general(1:N2,4), 'blue', general(1:N2, 1), general(1:N2,5), 'black', general(1:N2, 1), general(1:N2,6), 'red');
title ('energy');
xlabel ('t*w_p');
ylabel ('E erg');
legend(4, 'particle', 'electric','magnetic', 'full');
grid ;

figure(2);
plot (general(1:N2,1), general(1:N2,7), 'red', general(1:N2,1), general(1:N2,8), 'green', general(1:N2, 1), general(1:N2,9), 'blue');
title ('momentum');
xlabel ('t*w_p');
ylabel ('P g*cm/s');
legend(4, 'x', 'y','z');
grid ;