clear;

Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
LegendTitle = {'{\theta} = 0', '{\theta} = 10','{\theta} = 20', '{\theta} = 30', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};

Ee0 = importdata('Ee0.dat');
Fe0 = importdata('Fs0.dat');
Ee1 = importdata('Ee1.dat');
Fe1 = importdata('Fs1.dat');
Ee2 = importdata('Ee2.dat');
Fe2 = importdata('Fs2.dat');
Ee3 = importdata('Ee3.dat');
Fe3 = importdata('Fs3.dat');
Ee4 = importdata('Ee4.dat');
Fe4 = importdata('Fs4.dat');
Ee5 = importdata('Ee5.dat');
Fe5 = importdata('Fs5.dat');
Ee6 = importdata('Ee6.dat');
Fe6 = importdata('Fs6.dat');
Ee7 = importdata('Ee7.dat');
Fe7 = importdata('Fs7.dat');
Ee8 = importdata('Ee8.dat');
Fe8 = importdata('Fs8.dat');
Ee9 = importdata('Ee9.dat');
Fe9 = importdata('Fs9.dat');

N = size(Ee0,2);

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
title ('Fe');
xlabel ('Ee');
ylabel ('Fe');

plot(Ee0(1:N),Fe0(1:N),'Color',Color{1});
plot(Ee1(1:N),Fe1(1:N),'Color',Color{2});
plot(Ee2(1:N),Fe2(1:N),'Color',Color{3});
plot(Ee3(1:N),Fe3(1:N),'Color',Color{4});
plot(Ee4(1:N),Fe4(1:N),'Color',Color{5});
plot(Ee5(1:N),Fe5(1:N),'Color',Color{6});
plot(Ee6(1:N),Fe6(1:N),'Color',Color{7});
plot(Ee7(1:N),Fe7(1:N),'Color',Color{8});
plot(Ee8(1:N),Fe8(1:N),'Color',Color{9});
plot(Ee9(1:N),Fe9(1:N),'Color',Color{10});
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest');

grid ;