clear;
load flux.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx*Ny*Nz;
Nt = fix(size(flux,1)/NE);

ynumber = 2;
znumber = 2;

a = 0;
b = fix(Nt/2);
c = Nt-1;

Jx(1:Nx, 1:3) = 0;
Jy(1:Nx, 1:3) = 0;
Jz(1:Nx, 1:3) = 0;

%extJx(1:Nx, 1:3) = 0;
%extJy(1:Nx, 1:3) = 0;
%extJz(1:Nx, 1:3) = 0;


for i=1:Nx,
   Jx(i,1) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + a*NE, 1);
   Jx(i,2) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + b*NE, 1);
   Jx(i,3) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + c*NE, 1);
   Jy(i,1) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + a*NE, 2);
   Jy(i,2) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + b*NE, 2);
   Jy(i,3) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + c*NE, 2);
   Jz(i,1) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + a*NE, 3);
   Jz(i,2) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + b*NE, 3);
   Jz(i,3) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + c*NE, 3);
 
   %extJx(i,1) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + a*NE, 4);
   %extJx(i,2) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + b*NE, 4);
   %extJx(i,3) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + c*NE, 4);
   %extJy(i,1) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + a*NE, 5);
   %extJy(i,2) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + b*NE, 5);
   %extJy(i,3) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + c*NE, 5);
   %extJz(i,1) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + a*NE, 6);
   %extJz(i,2) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + b*NE, 6);
   %extJz(i,3) = flux(Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber + c*NE, 6);
end;
set(0, 'DefaultLineLineWidth', 2);

figure(1);
plot (Xfile(1:Nx,1),Jx(1:Nx,1), 'red',Xfile(1:Nx,1),Jx(1:Nx,2), 'green',Xfile(1:Nx,1),Jx(1:Nx,3), 'blue');
title ('Jx');
xlabel ('x/r_g');
ylabel ('flux');
grid ;

figure(2);
plot (Xfile(1:Nx,1),Jy(1:Nx, 1), 'red', Xfile(1:Nx,1), Jy(1:Nx, 2), 'green',Xfile(1:Nx,1),Jy(1:Nx, 3), 'blue');
title ('Jy');
xlabel ('x/r_g');
ylabel ('flux');
grid ;

figure(3);
plot (Xfile(1:Nx,1),Jz(1:Nx, 1), 'red', Xfile(1:Nx,1), Jz(1:Nx, 2), 'green', Xfile(1:Nx,1), Jz(1:Nx, 3), 'blue');
title ('Jz');
xlabel ('x/r_g');
ylabel ('flux');
grid ;