clear;
Efield = importdata('Efield_5.dat');
Bfield = importdata('Bfield_5.dat');
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
c = fix(Nt)-1;

Ex(1:Nz, 1:3) = 0;
Ey(1:Nz, 1:3) = 0;
Ez(1:Nz, 1:3) = 0;

Bx(1:Nz-1, 1:3) = 0;
By(1:Nz-1, 1:3) = 0;
Bz(1:Nz-1, 1:3) = 0;
Bnorm(1:Nz-1, 1:3) = 0;

B0=initialParameters(19);

middleZ(1:Nz-1) = 0;
Zgrid(1:Nz) = 0;
cv = initialParameters(10);
omega = initialParameters(20);


for i=1:Nz,
   %Xgrid(i) = (Xfile(i) - Xfile(2))*omega/cv;
   Zgrid(i) = (Zfile(i) - Zfile(2));
   Ex(i,1) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(ynumber-1) + i + a*NE, 1);
   Ex(i,2) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(ynumber-1) + i + b*NE, 1);
   Ex(i,3) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(ynumber-1) + i + c*NE, 1);
   Ey(i,1) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(ynumber-1) + i + a*NE, 2);
   Ey(i,2) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(ynumber-1) + i + b*NE, 2);
   Ey(i,3) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(ynumber-1) + i + c*NE, 2);
   Ez(i,1) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(ynumber-1) + i + a*NE, 3);
   Ez(i,2) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(ynumber-1) + i + b*NE, 3);
   Ez(i,3) = Efield((Nz)*(Ny)*(xnumber-1) + (Nz)*(ynumber-1) + i + c*NE, 3);
end;
%ynumber = 1;
%znumber = 1;
for i = 1:Nz-1,
   %middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2))*omega/cv;
   middleZ(i) = (0.5*(Zfile(i) + Zfile(i+1)) - Zfile(2));
   Bx(i, 1) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(ynumber-1) + i) + a*NB, 1);
   Bx(i, 2) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(ynumber-1) + i) + b*NB, 1);
   Bx(i, 3) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(ynumber-1) + i) + c*NB, 1);
   By(i, 1) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(ynumber-1) + i) + a*NB, 2);
   By(i, 2) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(ynumber-1) + i) + b*NB, 2);
   By(i, 3) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(ynumber-1) + i) + c*NB, 2);
   Bz(i, 1) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(ynumber-1) + i) + a*NB, 3);
   Bz(i, 2) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(ynumber-1) + i) + b*NB, 3);
   Bz(i, 3) = Bfield(((Nz-1)*(Ny-1)*(xnumber-1) + (Nz-1)*(ynumber-1) + i) + c*NB, 3);
   Bnorm(i, 1) = sqrt(By(i,1)*By(i,1) + Bz(i,1)*Bz(i,1))/B0;
   Bnorm(i, 2) = sqrt(By(i,2)*By(i,2) + Bz(i,2)*Bz(i,2))/B0;
   Bnorm(i, 3) = sqrt(By(i,3)*By(i,3) + Bz(i,3)*Bz(i,3))/B0;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 2);

figure(1);
plot (Zgrid(1:Nz),Ex(1:Nz,1), 'red',Zgrid(1:Nz),Ex(1:Nz,2), 'green',Zgrid(1:Nz),Ex(1:Nz,3), 'blue');
%title ('E_x');
xlabel ('z');
ylabel ('E_x gauss');
grid ;

figure(2);
plot (Zgrid(1:Nz),Ey(1:Nz, 1), 'red', Zgrid(1:Nz), Ey(1:Nz, 2), 'green',Zgrid(1:Nz),Ey(1:Nz, 3), 'blue');
%title ('E_y');
xlabel ('z');
ylabel ('E_y gauss');
grid ;

figure(3);
plot (Zgrid(1:Nz),Ez(1:Nz, 1), 'red', Zgrid(1:Nz), Ez(1:Nz, 2), 'green', Zgrid(1:Nz), Ez(1:Nz, 3), 'blue');
%title ('E_z');
xlabel ('z');
ylabel ('E_z gauss');
grid ;

figure(4);
plot (middleZ(1:Nz-1),Bx(1:Nz-1, 1), 'red', middleZ(1:Nz-1),Bx(1:Nz-1, 2), 'green', middleZ(1:Nz-1),Bx(1:Nz-1, 3), 'blue');
%title ('B_x');
xlabel ('z');
ylabel ('B_x gauss');
grid ;

figure(5);
plot (middleZ(1:Nz-1),By(1:Nz-1, 1), 'red', middleZ(1:Nz-1),By(1:Nz-1, 2), 'green', middleZ(1:Nz-1),By(1:Nz-1, 3), 'blue');
%title ('B_y');
xlabel ('z');
ylabel ('B_y gauss');
grid ;

figure(6);
plot (middleZ(1:Nz-1),Bz(1:Nz-1, 1), 'red', middleZ(1:Nz-1),Bz(1:Nz-1, 2), 'green', middleZ(1:Nz-1),Bz(1:Nz-1, 3), 'blue');
%title ('B_z');
xlabel ('z');
ylabel ('B_z');
%legend('t=0','t=T/2','t=T','Location','northeast');
grid ;

figure(7);
plot (middleZ(1:Nz-1),Bnorm(1:Nz-1, 1), 'red', middleZ(1:Nz-1),Bnorm(1:Nz-1, 2), 'green', middleZ(1:Nz-1),Bnorm(1:Nz-1, 3), 'blue');
%title ('B_z');
xlabel ('z');
ylabel ('B_{\perp}/B_0');
%legend('t=0','t=T/2','t=T','Location','northeast');
grid ;
