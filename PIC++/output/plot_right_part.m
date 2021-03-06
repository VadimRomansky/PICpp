clear;
load rightPart.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
load initialParameters.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NB = (Nx-1)*(Ny-1)*(Nz-1);
Nt = (size(rightPart, 1)/NB);
ynumber = 3;
znumber = 3;
a = 0;
b = fix(Nt/2);
c = fix(Nt)-1;


Bx(1:Nx-1, 1:3) = 0;
By(1:Nx-1, 1:3) = 0;
Bz(1:Nx-1, 1:3) = 0;

B0=initialParameters(19);

middleX(1:Nx-1) = 0;
Xgrid(1:Nx) = 0;
cv = initialParameters(10);
omega = initialParameters(21);
omegaElectron = initialParameters(20);




for i = 1:Nx-1,
   middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2))*omegaElectron/cv;
   %middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2));
   Bx(i, 1) = rightPart(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 1);
   Bx(i, 2) = rightPart(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 1);
   Bx(i, 3) = rightPart(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 1);
   By(i, 1) = rightPart(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 2);
   By(i, 2) = rightPart(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 2);
   By(i, 3) = rightPart(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 2);
   Bz(i, 1) = rightPart(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 3);
   Bz(i, 2) = rightPart(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 3);
   Bz(i, 3) = rightPart(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 3);
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(4);
plot (middleX(1:Nx-1),Bx(1:Nx-1, 1), 'red', middleX(1:Nx-1),Bx(1:Nx-1, 2), 'green', middleX(1:Nx-1),Bx(1:Nx-1, 3), 'blue');
%title ('B_x');
xlabel ('x');
ylabel ('right part x');
grid ;

figure(5);
plot (middleX(1:Nx-1),By(1:Nx-1, 1), 'red', middleX(1:Nx-1),By(1:Nx-1, 2), 'green', middleX(1:Nx-1),By(1:Nx-1, 3), 'blue');
%title ('B_y');
xlabel ('x');
ylabel ('right part y');
grid ;

figure(6);
plot (middleX(1:Nx-1),Bz(1:Nx-1, 1), 'red', middleX(1:Nx-1),Bz(1:Nx-1, 2), 'green', middleX(1:Nx-1),Bz(1:Nx-1, 3), 'blue');
%title ('B_z');
xlabel ('x\omega_p/c');
ylabel ('right part z');
%legend('t=0','t=T/2','t=T','Location','northeast');
grid ;
