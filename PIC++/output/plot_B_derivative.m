clear;
load rotEFile.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

NE = Nx*Ny*Nz;
Nt = (size(rotEFile, 1)/NE)-1;

ynumber = 1;
znumber = 1;

a = 0;
b = fix(Nt/2);
c = fix(Nt) - 1;


rotEx(1:Nx, 1:3) = 0;
rotEy(1:Nx, 1:3) = 0;
rotEz(1:Nx, 1:3) = 0;





for i=1:Nx,
   
   rotEx(i,1) = rotEFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 1);
   rotEx(i,2) = rotEFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 1);
   rotEx(i,3) = rotEFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 1);
   rotEy(i,1) = rotEFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 2);
   rotEy(i,2) = rotEFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 2);
   rotEy(i,3) = rotEFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 2);
   rotEz(i,1) = rotEFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 3);
   rotEz(i,2) = rotEFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*NE, 3);
   rotEz(i,3) = rotEFile((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*NE, 3);
   
end;

figure(1);
plot (Xfile(1:Nx,1),rotEx(1:Nx,1), 'red',Xfile(1:Nx,1),rotEx(1:Nx,2), 'green',Xfile(1:Nx,1),rotEx(1:Nx,3), 'blue');
title ('rotEx');
xlabel ('x cm');
ylabel ('rotE');
grid ;

figure(2);
plot (Xfile(1:Nx,1),rotEy(1:Nx, 1), 'red', Xfile(1:Nx,1), rotEy(1:Nx, 2), 'green',Xfile(1:Nx,1),rotEy(1:Nx, 3), 'blue');
title ('rotEy');
xlabel ('x cm');
ylabel ('rotE');
grid ;

figure(3);
plot (Xfile(1:Nx,1),rotEz(1:Nx, 1), 'red', Xfile(1:Nx,1), rotEz(1:Nx, 2), 'green', Xfile(1:Nx,1), rotEz(1:Nx, 3), 'blue');
title ('rotEz');
xlabel ('x cm');
ylabel ('rotE');
grid ;

