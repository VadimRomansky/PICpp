clear;
load Efield.dat;
load Bfield.dat;

NE = size(Efield, 1);
NB = size(Bfield, 1);

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
plot (Efield(1:NE,1),Efield(1:NE,2), 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('E_x gauss');
grid ;

figure(2);
plot (Efield(1:NE,1),Efield(1:NE,3), 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('E_y gauss');
grid ;

figure(3);
plot (Efield(1:NE,1),Efield(1:NE,4), 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('E_z gauss');
grid ;

figure(4);
plot (Bfield(1:NB,1),Bfield(1:NB,2), 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('B_x gauss');
grid ;

figure(5);
plot (Bfield(1:NB,1),Bfield(1:NB,3), 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('B_y gauss');
grid ;

figure(6);
plot (Bfield(1:NB,1),Bfield(1:NB,4), 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('B_z gauss');
grid ;