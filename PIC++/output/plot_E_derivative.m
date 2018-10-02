clear;
load rotBFile.dat;
load EderivativeFile.dat;
load flux.dat
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx*Ny*Nz;
Nt = (size(flux, 1)/NE)-1;

ynumber = 2;
znumber = 2;

a = 0;
b = fix(Nt/2);
c = fix(Nt);

Jx(1:Nx, 1:3) = 0;
Jy(1:Nx, 1:3) = 0;
Jz(1:Nx, 1:3) = 0;

rotBx(1:Nx, 1:3) = 0;
rotBy(1:Nx, 1:3) = 0;
rotBz(1:Nx, 1:3) = 0;

derEx(1:Nx, 1:3) = 0;
derEy(1:Nx, 1:3) = 0;
derEz(1:Nx, 1:3) = 0;




for i=1:Nx,
   Jx(i,1) = 4*3.14*flux((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 1);
   Jx(i,2) = 4*3.14*flux((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 1);
   Jx(i,3) = 4*3.14*flux((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 1);
   Jy(i,1) = 4*3.14*flux((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 2);
   Jy(i,2) = 4*3.14*flux((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 2);
   Jy(i,3) = 4*3.14*flux((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 2);
   Jz(i,1) = 4*3.14*flux((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 3);
   Jz(i,2) = 4*3.14*flux((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 3);
   Jz(i,3) = 4*3.14*flux((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 3);
   
   rotBx(i,1) = rotBFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 1);
   rotBx(i,2) = rotBFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 1);
   rotBx(i,3) = rotBFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 1);
   rotBy(i,1) = rotBFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 2);
   rotBy(i,2) = rotBFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 2);
   rotBy(i,3) = rotBFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 2);
   rotBz(i,1) = rotBFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 3);
   rotBz(i,2) = rotBFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 3);
   rotBz(i,3) = rotBFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 3);
   
   derEx(i,1) = EderivativeFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 1);
   derEx(i,2) = EderivativeFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 1);
   derEx(i,3) = EderivativeFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 1);
   derEy(i,1) = EderivativeFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 2);
   derEy(i,2) = EderivativeFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 2);
   derEy(i,3) = EderivativeFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 2);
   derEz(i,1) = EderivativeFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 3);
   derEz(i,2) = EderivativeFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 3);
   derEz(i,3) = EderivativeFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 3);
end;

set(0, 'DefaultLineLineWidth', 2);

figure(1);
plot (Xfile(1:Nx,1),Jx(1:Nx,1), 'red',Xfile(1:Nx,1),Jx(1:Nx,2), 'green',Xfile(1:Nx,1),Jx(1:Nx,3), 'blue');
title ('4*pi*x');
xlabel ('x cm');
ylabel ('flux');
grid ;

figure(2);
plot (Xfile(1:Nx,1),Jy(1:Nx, 1), 'red', Xfile(1:Nx,1), Jy(1:Nx, 2), 'green',Xfile(1:Nx,1),Jy(1:Nx, 3), 'blue');
title ('4*pi*Jy');
xlabel ('x cm');
ylabel ('flux');
grid ;

figure(3);
plot (Xfile(1:Nx,1),Jz(1:Nx, 1), 'red', Xfile(1:Nx,1), Jz(1:Nx, 2), 'green', Xfile(1:Nx,1), Jz(1:Nx, 3), 'blue');
title ('4*pi*Jz');
xlabel ('x cm');
ylabel ('flux');
grid ;

figure(4);
plot (Xfile(1:Nx,1),rotBx(1:Nx,1), 'red',Xfile(1:Nx,1),rotBx(1:Nx,2), 'green',Xfile(1:Nx,1),rotBx(1:Nx,3), 'blue');
title ('rotBx');
xlabel ('x cm');
ylabel ('rotB');
grid ;

figure(5);
plot (Xfile(1:Nx,1),rotBy(1:Nx, 1), 'red', Xfile(1:Nx,1), rotBy(1:Nx, 2), 'green',Xfile(1:Nx,1),rotBy(1:Nx, 3), 'blue');
title ('rotBy');
xlabel ('x cm');
ylabel ('rotB');
grid ;

figure(6);
plot (Xfile(1:Nx,1),rotBz(1:Nx, 1), 'red', Xfile(1:Nx,1), rotBz(1:Nx, 2), 'green', Xfile(1:Nx,1), rotBz(1:Nx, 3), 'blue');
title ('rotBz');
xlabel ('x cm');
ylabel ('rotB');
grid ;

figure(7);
plot (Xfile(1:Nx,1),derEx(1:Nx,1), 'red',Xfile(1:Nx,1),derEx(1:Nx,2), 'green',Xfile(1:Nx,1),derEx(1:Nx,3), 'blue');
title ('derEx');
xlabel ('x cm');
ylabel ('derE');
grid ;

figure(8);
plot (Xfile(1:Nx,1),derEy(1:Nx, 1), 'red', Xfile(1:Nx,1), derEy(1:Nx, 2), 'green',Xfile(1:Nx,1),derEy(1:Nx, 3), 'blue');
title ('derEy');
xlabel ('x cm');
ylabel ('derE');
grid ;

figure(9);
plot (Xfile(1:Nx,1),derEz(1:Nx, 1), 'red', Xfile(1:Nx,1), derEz(1:Nx, 2), 'green', Xfile(1:Nx,1), derEz(1:Nx, 3), 'blue');
title ('derEz');
xlabel ('x cm');
ylabel ('derE');
grid ;

