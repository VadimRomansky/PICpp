clear;
load traectory_proton.dat;
N1=1;
N2=size(traectory_proton,1);


figure(1);
plot (traectory_proton(1:N2,1),traectory_proton(1:N2,3),'red',traectory_proton(1:N2,1),traectory_proton(1:N2,4),'green',traectory_proton(1:N2,1),traectory_proton(1:N2,5),'blue');
title ('r');
xlabel ('t');
ylabel (')');
legend(4, 'x', 'y','z');
grid ;

figure(2);
plot (traectory_proton(1:N2,1),traectory_proton(1:N2,6),'red',traectory_proton(1:N2,1),traectory_proton(1:N2,7),'green',traectory_proton(1:N2,1),traectory_proton(1:N2,8),'blue', traectory_proton(1:N2,1),traectory_proton(1:N2,9),'black');
title ('p');
xlabel ('t');
ylabel ('p');
legend(4, 'px', 'py','pz', 'p');
grid ;

figure(3);
plot (traectory_proton(1:N2,4),traectory_proton(1:N2,5));
title ('r');
xlabel ('y');
ylabel ('z');
grid ;
