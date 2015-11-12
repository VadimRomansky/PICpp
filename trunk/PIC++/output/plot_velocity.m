clear;
load velocity.dat;
load velocity_electron.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = Nx*Ny*Nz;
Nt = size(velocity, 1)/N;
ynumber = 1;
znumber = 1;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

Vx(1:Nx, 1:3) = 0;
Vy(1:Nx, 1:3) = 0;
Vz(1:Nx, 1:3) = 0;

Velectronx(1:Nx, 1:3) = 0;
Velectrony(1:Nx, 1:3) = 0;
Velectronz(1:Nx, 1:3) = 0;

middleX(1:Nx-1) = 0;

for i=1:Nx,
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
   Vx(i,1) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*N, 1);
   Vx(i,2) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*N, 1);
   Vx(i,3) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*N, 1);
   Vy(i,1) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*N, 2);
   Vy(i,2) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*N, 2);
   Vy(i,3) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*N, 2);
   Vz(i,1) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*N, 3);
   Vz(i,2) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*N, 3);
   Vz(i,3) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*N, 3);
   
   Velectronx(i,1) = velocity_electron((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*N, 1);
   Velectronx(i,2) = velocity_electron((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*N, 1);
   Velectronx(i,3) = velocity_electron((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*N, 1);
   Velectrony(i,1) = velocity_electron((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*N, 2);
   Velectrony(i,2) = velocity_electron((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*N, 2);
   Velectrony(i,3) = velocity_electron((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*N, 2);
   Velectronz(i,1) = velocity_electron((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*N, 3);
   Velectronz(i,2) = velocity_electron((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*N, 3);
   Velectronz(i,3) = velocity_electron((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*N, 3);
end;
figure(1);
plot (middleX(1:Nx),Vx(1:Nx,1), 'red',middleX(1:Nx),Vx(1:Nx,2), 'green',middleX(1:Nx),Vx(1:Nx,3), 'blue');
title ('Vx');
xlabel ('x/r_g');
ylabel ('V cm/s');
grid ;

figure(2);
plot (middleX(1:Nx),Vy(1:Nx, 1), 'red', middleX(1:Nx), Vy(1:Nx, 2), 'green',middleX(1:Nx),Vy(1:Nx, 3), 'blue');
title ('Vy');
xlabel ('x/r_g');
ylabel ('V cm/s');
grid ;

figure(3);
plot (middleX(1:Nx),Vz(1:Nx, 1), 'red', middleX(1:Nx), Vz(1:Nx, 2), 'green', middleX(1:Nx), Vz(1:Nx, 3), 'blue');
title ('Vz');
xlabel ('x/r_g');
ylabel ('V cm/s');
grid ;

figure(4);
plot (middleX(1:Nx),Velectronx(1:Nx,1), 'red',middleX(1:Nx),Velectronx(1:Nx,2), 'green',middleX(1:Nx),Velectronx(1:Nx,3), 'blue');
title ('Vx');
xlabel ('x/r_g');
ylabel ('V cm/s');
grid ;

figure(5);
plot (middleX(1:Nx),Velectrony(1:Nx, 1), 'red', middleX(1:Nx), Velectrony(1:Nx, 2), 'green',middleX(1:Nx),Velectrony(1:Nx, 3), 'blue');
title ('Vy');
xlabel ('x/r_g');
ylabel ('V cm/s');
grid ;

figure(6);
plot (middleX(1:Nx),Velectronz(1:Nx, 1), 'red', middleX(1:Nx), Velectronz(1:Nx, 2), 'green', middleX(1:Nx), Velectronz(1:Nx, 3), 'blue');
title ('Vz');
xlabel ('x/r_g');
ylabel ('V cm/s');
grid ;

