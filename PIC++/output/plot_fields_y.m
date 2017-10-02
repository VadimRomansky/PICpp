clear;
load Efield.dat;
load Bfield.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
load initialParameters.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx*Ny*Nz;
NB = (Nx-1)*(Ny-1)*(Nz-1);
Nt = (size(Efield, 1)/NE);
NtB = (size(Bfield, 1)/NB);
xnumber = 2;
ynumber = 2;
znumber = 2;
a = 0;
b = fix(Nt/2);
c = fix(Nt)-2;

Ex(1:Ny, 1:3) = 0;
Ey(1:Ny, 1:3) = 0;
Ez(1:Ny, 1:3) = 0;

Bx(1:Ny-1, 1:3) = 0;
By(1:Ny-1, 1:3) = 0;
Bz(1:Ny-1, 1:3) = 0;
Bnorm(1:Ny-1, 1:3) = 0;

B0=initialParameters(19);

middleY(1:Ny-1) = 0;
Ygrid(1:Ny) = 0;
cv = initialParameters(10);
omega = initialParameters(20);


for i=1:Ny,
   %Xgrid(i) = (Xfile(i) - Xfile(2))*omega/cv;
   Ygrid(i) = (Yfile(i) - Yfile(2));
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
   %middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2))*omega/cv;
   middleY(i) = (0.5*(Yfile(i) + Yfile(i+1)) - Yfile(2));
   Bx(i, 1) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + a*NB, 1);
   Bx(i, 2) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + b*NB, 1);
   Bx(i, 3) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + c*NB, 1);
   By(i, 1) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + a*NB, 2);
   By(i, 2) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + b*NB, 2);
   By(i, 3) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + c*NB, 2);
   Bz(i, 1) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + a*NB, 3);
   Bz(i, 2) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + b*NB, 3);
   Bz(i, 3) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(i-1) + znumber) + c*NB, 3);
   Bnorm(i, 1) = sqrt(By(i,1)*By(i,1) + Bz(i,1)*Bz(i,1))/B0;
   Bnorm(i, 2) = sqrt(By(i,2)*By(i,2) + Bz(i,2)*Bz(i,2))/B0;
   Bnorm(i, 3) = sqrt(By(i,3)*By(i,3) + Bz(i,3)*Bz(i,3))/B0;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
plot (Ygrid(1:Ny),Ex(1:Ny,1), 'red',Ygrid(1:Ny),Ex(1:Ny,2), 'green',Ygrid(1:Ny),Ex(1:Ny,3), 'blue');
%title ('E_x');
xlabel ('y');
ylabel ('E_x gauss');
grid ;

figure(2);
plot (Ygrid(1:Ny),Ey(1:Ny, 1), 'red', Ygrid(1:Ny), Ey(1:Ny, 2), 'green',Ygrid(1:Ny),Ey(1:Ny, 3), 'blue');
%title ('E_y');
xlabel ('y');
ylabel ('E_y gauss');
grid ;

figure(3);
plot (Ygrid(1:Ny),Ez(1:Ny, 1), 'red', Ygrid(1:Ny), Ez(1:Ny, 2), 'green', Ygrid(1:Ny), Ez(1:Ny, 3), 'blue');
%title ('E_z');
xlabel ('y');
ylabel ('E_z gauss');
grid ;

figure(4);
plot (middleY(1:Ny-1),Bx(1:Ny-1, 1), 'red', middleY(1:Ny-1),Bx(1:Ny-1, 2), 'green', middleY(1:Ny-1),Bx(1:Ny-1, 3), 'blue');
%title ('B_x');
xlabel ('y');
ylabel ('B_x gauss');
grid ;

figure(5);
plot (middleY(1:Ny-1),By(1:Ny-1, 1), 'red', middleY(1:Ny-1),By(1:Ny-1, 2), 'green', middleY(1:Ny-1),By(1:Ny-1, 3), 'blue');
%title ('B_y');
xlabel ('y');
ylabel ('B_y gauss');
grid ;

figure(6);
plot (middleY(1:Ny-1),Bz(1:Ny-1, 1), 'red', middleY(1:Ny-1),Bz(1:Ny-1, 2), 'green', middleY(1:Ny-1),Bz(1:Ny-1, 3), 'blue');
%title ('B_z');
xlabel ('y');
ylabel ('B_z');
%legend('t=0','t=T/2','t=T','Location','northeast');
grid ;

figure(7);
plot (middleY(1:Ny-1),Bnorm(1:Ny-1, 1), 'red', middleY(1:Ny-1),Bnorm(1:Ny-1, 2), 'green', middleY(1:Ny-1),Bnorm(1:Ny-1, 3), 'blue');
%title ('B_z');
xlabel ('y');
ylabel ('B_{\perp}/B_0');
%legend('t=0','t=T/2','t=T','Location','northeast');
grid ;
