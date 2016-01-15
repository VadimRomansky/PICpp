clear;
load fluxFile.dat;
load Xfile.dat;

Nx = size(Xfile, 1);

NE = Nx;
NB = Nx - 1;
Nt = fix(size(fluxFile,1)/NE);

a = 0;
b = fix(Nt/2);
c = Nt - 1;

Jx(1:NE, 1:3) = 0;
Jy(1:NE, 1:3) = 0;
Jz(1:NE, 1:3) = 0;

extJx(1:NE, 1:3) = 0;
extJy(1:NE, 1:3) = 0;
extJz(1:NE, 1:3) = 0;


for i=1:NE,
   Jx(i,1) = fluxFile(i + a*NE, 1);
   Jx(i,2) = fluxFile(i + b*NE, 1);
   Jx(i,3) = fluxFile(i + c*NE, 1);
   Jy(i,1) = fluxFile(i + a*NE, 2);
   Jy(i,2) = fluxFile(i + b*NE, 2);
   Jy(i,3) = fluxFile(i + c*NE, 2);
   Jz(i,1) = fluxFile(i + a*NE, 3);
   Jz(i,2) = fluxFile(i + b*NE, 3);
   Jz(i,3) = fluxFile(i + c*NE, 3);
 
   extJx(i,1) = fluxFile(i + a*NE, 4);
   extJx(i,2) = fluxFile(i + b*NE, 4);
   extJx(i,3) = fluxFile(i + c*NE, 4);
   extJy(i,1) = fluxFile(i + a*NE, 5);
   extJy(i,2) = fluxFile(i + b*NE, 5);
   extJy(i,3) = fluxFile(i + c*NE, 5);
   extJz(i,1) = fluxFile(i + a*NE, 6);
   extJz(i,2) = fluxFile(i + b*NE, 6);
   extJz(i,3) = fluxFile(i + c*NE, 6);
end;
figure(1);
plot (Xfile(1:NE,1),Jx(1:NE,1), 'red',Xfile(1:NE,1),Jx(1:NE,2), 'green',Xfile(1:NE,1),Jx(1:NE,3), 'blue');
title ('Jx');
xlabel ('x cm');
ylabel ('flux');
grid ;

figure(2);
plot (Xfile(1:NE,1),Jy(1:NE, 1), 'red', Xfile(1:NE,1), Jy(1:NE, 2), 'green',Xfile(1:NE,1),Jy(1:NE, 3), 'blue');
title ('Jy');
xlabel ('x cm');
ylabel ('flux');
grid ;

figure(3);
plot (Xfile(1:NE,1),Jz(1:NE, 1), 'red', Xfile(1:NE,1), Jz(1:NE, 2), 'green', Xfile(1:NE,1), Jz(1:NE, 3), 'blue');
title ('Jz');
xlabel ('x cm');
ylabel ('flux');
grid ;

figure(4);
plot (Xfile(1:NE,1),extJx(1:NE,1), 'red',Xfile(1:NE,1),extJx(1:NE,2), 'green',Xfile(1:NE,1),extJx(1:NE,3), 'blue');
title ('external Jx');
xlabel ('x cm');
ylabel ('flux');
grid ;

figure(5);
plot (Xfile(1:NE,1),extJy(1:NE, 1), 'red', Xfile(1:NE,1), extJy(1:NE, 2), 'green',Xfile(1:NE,1),extJy(1:NE, 3), 'blue');
title ('external Jy');
xlabel ('x cm');
ylabel ('flux');
grid ;

figure(6);
plot (Xfile(1:NE,1),extJz(1:NE, 1), 'red', Xfile(1:NE,1), extJz(1:NE, 2), 'green', Xfile(1:NE,1), extJz(1:NE, 3), 'blue');
title ('external Jz');
xlabel ('x cm');
ylabel ('flux');
grid ;


