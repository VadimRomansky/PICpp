general = importdata('output/extra_iterations.dat');
N = size(general,1);

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 1.5);

figure(1);
plot(general(1:N,1),general(1:N,2));
title ('t(N)');
xlabel ('N');
ylabel ('t');
grid ;

figure(2);
plot(general(1:N,2),general(1:N,3));
title ('M');
xlabel ('t');
ylabel ('M');
grid ;

figure(3);
hold on;
plot(general(1:N,2),general(1:N,5),"Color","red");
plot(general(1:N,2),general(1:N,6),"Color","green");
plot(general(1:N,2),general(1:N,7),"Color","blue");
plot(general(1:N,2),general(1:N,8),"Color","black");
plot(general(1:N,2),general(1:N,9),"Color","yellow");
plot(general(1:N,2),general(1:N,12),"Color","magenta");
plot(general(1:N,2),general(1:N,13),"Color","cyan");

legend('total energy', 'kinetic energy', 'thermal energy', 'particle enrgy', 'magnetic energy', 'injected energy', 'u grad p energy exchange');
title ('E');
xlabel ('t');
ylabel ('E');
grid ;

figure(4);
hold on;
plot(general(1:N,2),general(1:N,10),"Color","red");
plot(general(1:N,2),general(1:N,11),"Color","blue");
legend('injected particles','total particles');
title ('N particles');
xlabel ('t');
ylabel ('N particles');
grid ;

