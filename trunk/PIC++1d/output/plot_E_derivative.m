clear;
load rotBFile.dat;
load EderivativeFile.dat;
load fluxFile.dat
load Xfile.dat;

Nx = size(Xfile, 1);

NE = Nx;
Nt = fix(size(fluxFile,1)/NE);

a = 0;
b = fix(Nt/2);
c = Nt - 1;

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
   Jx(i,1) = 4*3.14*fluxFile(i + a*NE, 1);
   Jx(i,2) = 4*3.14*fluxFile(i + b*NE, 1);
   Jx(i,3) = 4*3.14*fluxFile(i + c*NE, 1);
   Jy(i,1) = 4*3.14*fluxFile(i + a*NE, 2);
   Jy(i,2) = 4*3.14*fluxFile(i + b*NE, 2);
   Jy(i,3) = 4*3.14*fluxFile(i + c*NE, 2);
   Jz(i,1) = 4*3.14*fluxFile(i + a*NE, 3);
   Jz(i,2) = 4*3.14*fluxFile(i + b*NE, 3);
   Jz(i,3) = 4*3.14*fluxFile(i + c*NE, 3);
   
   rotBx(i,1) = rotBFile(i + a*NE, 1);
   rotBx(i,2) = rotBFile(i + b*NE, 1);
   rotBx(i,3) = rotBFile(i + c*NE, 1);
   rotBy(i,1) = rotBFile(i + a*NE, 2);
   rotBy(i,2) = rotBFile(i + b*NE, 2);
   rotBy(i,3) = rotBFile(i + c*NE, 2);
   rotBz(i,1) = rotBFile(i + a*NE, 3);
   rotBz(i,2) = rotBFile(i + b*NE, 3);
   rotBz(i,3) = rotBFile(i + c*NE, 3);
   
   derEx(i,1) = EderivativeFile(i + a*NE, 1);
   derEx(i,2) = EderivativeFile(i + b*NE, 1);
   derEx(i,3) = EderivativeFile(i + c*NE, 1);
   derEy(i,1) = EderivativeFile(i + a*NE, 2);
   derEy(i,2) = EderivativeFile(i + b*NE, 2);
   derEy(i,3) = EderivativeFile(i + c*NE, 2);
   derEz(i,1) = EderivativeFile(i + a*NE, 3);
   derEz(i,2) = EderivativeFile(i + b*NE, 3);
   derEz(i,3) = EderivativeFile(i + c*NE, 3);
end;

figure(1);
plot (Xfile(1:Nx,1),Jx(1:Nx,1), 'red',Xfile(1:Nx,1),Jx(1:Nx,2), 'green',Xfile(1:Nx,1),Jx(1:Nx,3), 'blue');
title ('Jx');
xlabel ('x cm');
ylabel ('flux');
grid ;

figure(2);
plot (Xfile(1:Nx,1),Jy(1:Nx, 1), 'red', Xfile(1:Nx,1), Jy(1:Nx, 2), 'green',Xfile(1:Nx,1),Jy(1:Nx, 3), 'blue');
title ('Jy');
xlabel ('x cm');
ylabel ('flux');
grid ;

figure(3);
plot (Xfile(1:Nx,1),Jz(1:Nx, 1), 'red', Xfile(1:Nx,1), Jz(1:Nx, 2), 'green', Xfile(1:Nx,1), Jz(1:Nx, 3), 'blue');
title ('Jz');
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

