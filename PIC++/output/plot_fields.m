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
%Nt=3;
%Nt = 23;
NtB = (size(Bfield, 1)/NB);
ynumber = 2;
znumber = 2;
a = 0;
b = fix(Nt/2);
c = fix(Nt)-1;

Ex(1:Nx, 1:3) = 0;
Ey(1:Nx, 1:3) = 0;
Ez(1:Nx, 1:3) = 0;

Bx(1:Nx-1, 1:3) = 0;
By(1:Nx-1, 1:3) = 0;
Bz(1:Nx-1, 1:3) = 0;
Bnorm(1:Nx-1, 1:3) = 0;

B0=initialParameters(19);

middleX(1:Nx-1) = 0;
Xgrid(1:Nx) = 0;
cv = initialParameters(10);
omega = initialParameters(21);
omegaElectron = initialParameters(20);


for i=1:Nx,
   %Xgrid(i) = (Xfile(i) - Xfile(2))*omegaElectron/cv;
   Xgrid(i) = (Xfile(i) - Xfile(2));
   Ex(i,1) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 1);
   Ex(i,2) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 1);
   Ex(i,3) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 1);
   Ey(i,1) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 2);
   Ey(i,2) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 2);
   Ey(i,3) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 2);
   Ez(i,1) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 3);
   Ez(i,2) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 3);
   Ez(i,3) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 3);
end;
%ynumber = 1;
%znumber = 1;
for i = 1:Nx-1,
    
   %middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2))*omegaElectron/cv;
   middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2));
   Bx(i, 1) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 1);
   Bx(i, 2) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 1);
   Bx(i, 3) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 1);
   By(i, 1) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 2);
   By(i, 2) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 2);
   By(i, 3) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 2);
   Bz(i, 1) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 3);
   Bz(i, 2) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 3);
   Bz(i, 3) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 3);
   Bnorm(i, 1) = sqrt(By(i,1)*By(i,1) + Bz(i,1)*Bz(i,1))/B0;
   Bnorm(i, 2) = sqrt(By(i,2)*By(i,2) + Bz(i,2)*Bz(i,2))/B0;
   Bnorm(i, 3) = sqrt(By(i,3)*By(i,3) + Bz(i,3)*Bz(i,3))/B0;
   %Bz(i, 1) = Bz(i, 1)/B0;
   %Bz(i, 2) = Bz(i, 2)/B0;
   %Bz(i, 3) = Bz(i, 3)/B0;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
plot (Xgrid(1:Nx),Ex(1:Nx,1), 'red',Xgrid(1:Nx),Ex(1:Nx,2), 'green',Xgrid(1:Nx),Ex(1:Nx,3), 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('E_x gauss');
grid ;

figure(2);
plot (Xgrid(1:Nx),Ey(1:Nx, 1), 'red', Xgrid(1:Nx), Ey(1:Nx, 2), 'green',Xgrid(1:Nx),Ey(1:Nx, 3), 'blue');
%title ('E_y');
xlabel ('x');
ylabel ('E_y gauss');
grid ;

figure(3);
plot (Xgrid(1:Nx),Ez(1:Nx, 1), 'red', Xgrid(1:Nx), Ez(1:Nx, 2), 'green', Xgrid(1:Nx), Ez(1:Nx, 3), 'blue');
%title ('E_z');
xlabel ('x');
ylabel ('E_z gauss');
grid ;

figure(4);
plot (middleX(1:Nx-1),Bx(1:Nx-1, 1), 'red', middleX(1:Nx-1),Bx(1:Nx-1, 2), 'green', middleX(1:Nx-1),Bx(1:Nx-1, 3), 'blue');
%title ('B_x');
xlabel ('x');
ylabel ('B_x gauss');
grid ;

figure(5);
plot (middleX(1:Nx-1),By(1:Nx-1, 1), 'red', middleX(1:Nx-1),By(1:Nx-1, 2), 'green', middleX(1:Nx-1),By(1:Nx-1, 3), 'blue');
%title ('B_y');
xlabel ('x');
ylabel ('B_y gauss');
grid ;

figure(6);
plot (middleX(1:Nx-1),Bz(1:Nx-1, 1), 'red', middleX(1:Nx-1),Bz(1:Nx-1, 2), 'green', middleX(1:Nx-1),Bz(1:Nx-1, 3), 'blue');
%title ('B_z');
xlabel ('x\omega_p/c');
ylabel ('B_z');
%legend('t=0','t=T/2','t=T','Location','northeast');
grid ;

figure(7);
plot (middleX(1:Nx-1),Bnorm(1:Nx-1, 1), 'red', middleX(1:Nx-1),Bnorm(1:Nx-1, 2), 'green', middleX(1:Nx-1),Bnorm(1:Nx-1, 3), 'blue');
%title ('B_z');
xlabel ('x\omega_p/c');
ylabel ('B_{\perp}/B_0');
%legend('t=0','t=T/2','t=T','Location','northeast');
grid ;
