clear;
load Efield.dat;
load Bfield.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx*Ny*Nz;
NB = (Nx-1)*(Ny-1)*(Nz-1);
Nt = (size(Efield, 1)/NE)-1;
NtB = (size(Bfield, 1)/NB)-1;
ynumber = 1;
znumber = 1;
xnumber = 1;

a = 0;
b = fix(Nt/2);
c = Nt-1;

Ex(1:Ny, 1:3) = 0;
Ey(1:Ny, 1:3) = 0;
Ez(1:Ny, 1:3) = 0;

Bx(1:Ny-1, 1:3) = 0;
By(1:Ny-1, 1:3) = 0;
Bz(1:Ny-1, 1:3) = 0;

middleY(1:Ny-1) = 0;


for i=1:Ny,
   Ex(i,1) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + a*NE, 1);
   Ex(i,2) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + b*NE, 1);
   Ex(i,3) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + c*NE, 1);
   Ey(i,1) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + a*NE, 2);
   Ey(i,2) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + b*NE, 2);
   Ey(i,3) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + c*NE, 2);
   Ez(i,1) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + a*NE, 3);
   Ez(i,2) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + b*NE, 3);
   Ez(i,3) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + c*NE, 3);
end;
%ynumber = 1;
%znumber = 1;
for i = 1:Ny-1,
   middleY(i) = 0.5*(Yfile(i) + Yfile(i+1));
   Bx(i, 1) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + a*NB, 1);
   Bx(i, 2) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + b*NB, 1);
   Bx(i, 3) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + c*NB, 1);
   By(i, 1) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + a*NB, 2);
   By(i, 2) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + b*NB, 2);
   By(i, 3) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + c*NB, 2);
   Bz(i, 1) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + a*NB, 3);
   Bz(i, 2) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + b*NB, 3);
   Bz(i, 3) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + c*NB, 3);
end;

figure(1);
plot (Yfile(1:Ny,1),Ex(1:Ny,1), 'red',Yfile(1:Ny,1),Ex(1:Ny,2), 'green',Yfile(1:Ny,1),Ex(1:Ny,3), 'blue');
title ('Ex');
xlabel ('y/r_g');
ylabel ('E gauss');
grid ;

figure(2);
plot (Yfile(1:Ny,1),Ey(1:Ny, 1), 'red', Yfile(1:Ny,1), Ey(1:Ny, 2), 'green',Yfile(1:Ny,1),Ey(1:Ny, 3), 'blue');
title ('Ey');
xlabel ('y/r_g');
ylabel ('E gauss');
grid ;

figure(3);
plot (Yfile(1:Ny,1),Ez(1:Ny, 1), 'red', Yfile(1:Ny,1), Ez(1:Ny, 2), 'green', Yfile(1:Ny,1), Ez(1:Ny, 3), 'blue');
title ('Ez');
xlabel ('y/r_g');
ylabel ('E gauss');
grid ;

figure(4);
plot (middleY(1:Ny-1),Bx(1:Ny-1, 1), 'red', middleY(1:Ny-1),Bx(1:Ny-1, 2), 'green', middleY(1:Ny-1),Bx(1:Ny-1, 3), 'blue');
title ('Bx');
xlabel ('y/r_g');
ylabel ('B gauss');
grid ;

figure(5);
plot (middleY(1:Ny-1),By(1:Ny-1, 1), 'red', middleY(1:Ny-1),By(1:Ny-1, 2), 'green', middleY(1:Ny-1),By(1:Ny-1, 3), 'blue');
title ('By');
xlabel ('y/r_g');
ylabel ('B gauss');
grid ;

figure(6);
plot (middleY(1:Ny-1),Bz(1:Ny-1, 1), 'red', middleY(1:Ny-1),Bz(1:Ny-1, 2), 'green', middleY(1:Ny-1),Bz(1:Ny-1, 3), 'blue');
title ('Bz');
xlabel ('y/r_g');
ylabel ('B gauss');
grid ;

