clear;
load Efield.dat;
load Bfield.dat;
load Xfile.dat;

Nx = size(Xfile, 1);

NE = Nx;
NB = Nx - 1;
Nt = fix(size(Bfield,1)/Nx);

a = 0;
b = fix(Nt/2);
c = Nt;

Ex(1:NE, 1:3) = 0;
Ey(1:NE, 1:3) = 0;
Ez(1:NE, 1:3) = 0;

Bx(1:NB, 1:3) = 0;
By(1:NB, 1:3) = 0;
Bz(1:NB, 1:3) = 0;

middleX(1:NB) = 0;


for i=1:NE,
   Ex(i,1) = Efield(i + a*NE, 1);
   Ex(i,2) = Efield(i + b*NE, 1);
   Ex(i,3) = Efield(i + c*NE, 1);
   Ey(i,1) = Efield(i + a*NE, 2);
   Ey(i,2) = Efield(i + b*NE, 2);
   Ey(i,3) = Efield(i + c*NE, 2);
   Ez(i,1) = Efield(i + a*NE, 3);
   Ez(i,2) = Efield(i + b*NE, 3);
   Ez(i,3) = Efield(i + c*NE, 3);
end;

for i=1:NB,
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
   Bx(i, 1) = Bfield(i + a*NB, 1);
   Bx(i, 2) = Bfield(i + b*NB, 1);
   Bx(i, 3) = Bfield(i + c*NB, 1);
   By(i, 1) = Bfield(i + a*NB, 2);
   By(i, 2) = Bfield(i + b*NB, 2);
   By(i, 3) = Bfield(i + c*NB, 2);
   Bz(i, 1) = Bfield(i + a*NB, 3);
   Bz(i, 2) = Bfield(i + b*NB, 3);
   Bz(i, 3) = Bfield(i + c*NB, 3);
end;
figure(1);
plot (Xfile(1:NE,1),Ex(1:NE,1), 'red',Xfile(1:NE,1),Ex(1:NE,2), 'green',Xfile(1:NE,1),Ex(1:NE,3), 'blue');
title ('Ex');
xlabel ('x/r_g');
ylabel ('E gauss');
grid ;

figure(2);
plot (Xfile(1:NE,1),Ey(1:NE, 1), 'red', Xfile(1:NE,1), Ey(1:NE, 2), 'green',Xfile(1:NE,1),Ey(1:NE, 3), 'blue');
title ('Ey');
xlabel ('x/r_g');
ylabel ('E gauss');
grid ;

figure(3);
plot (Xfile(1:NE,1),Ez(1:NE, 1), 'red', Xfile(1:NE,1), Ez(1:NE, 2), 'green', Xfile(1:NE,1), Ez(1:NE, 3), 'blue');
title ('Ez');
xlabel ('x/r_g');
ylabel ('E gauss');
grid ;

figure(4);
plot (middleX(1:NB),Bx(1:NB, 1), 'red', middleX(1:NB),Bx(1:NB, 2), 'green', middleX(1:NB),Bx(1:NB, 3), 'blue');
title ('Bx');
xlabel ('x cm');
ylabel ('B gauss');
grid ;

figure(5);
plot (middleX(1:NB),By(1:NB, 1), 'red', middleX(1:NB),By(1:NB, 2), 'green', middleX(1:NB),By(1:NB, 3), 'blue');
title ('By');
xlabel ('x cm');
ylabel ('B gauss');
grid ;

figure(6);
plot (middleX(1:NB),Bz(1:NB, 1), 'red', middleX(1:NB),Bz(1:NB, 2), 'green', middleX(1:NB),Bz(1:NB, 3), 'blue');
title ('Bz');
xlabel ('x cm');
ylabel ('B gauss');
grid ;

figure(7);
plot (middleX(1:NB),Bz(1:NB, 1), 'red',Xfile(1:NE,1),Ey(1:NE, 1), 'red', middleX(1:NB),Bz(1:NB, 2), 'green', Xfile(1:NE,1),Ey(1:NE, 2), 'green', middleX(1:NB), Bz(1:NB, 3), 'blue', Xfile(1:NE,1),Ey(1:NE, 3), 'blue');
title ('Bz and Ey');
xlabel ('x cm');
ylabel ('gauss');
grid ;

figure(8);
plot (middleX(1:NB),By(1:NB, 1), 'red',Xfile(1:NE,1),Ez(1:NE, 1), 'red', middleX(1:NB),By(1:NB, 2), 'green', Xfile(1:NE,1),Ez(1:NE, 2), 'green', middleX(1:NB), By(1:NB, 3), 'blue', Xfile(1:NE,1),Ez(1:NE, 3), 'blue');
title ('By and Ez');
xlabel ('x cm');
ylabel ('gauss');
grid ;


