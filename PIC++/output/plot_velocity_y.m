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
xnumber = 1;
ynumber = 1;
znumber = 1;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

Vx(1:Ny, 1:3) = 0;
Vy(1:Ny, 1:3) = 0;
Vz(1:Ny, 1:3) = 0;

Velectronx(1:Ny, 1:3) = 0;
Velectrony(1:Ny, 1:3) = 0;
Velectronz(1:Ny, 1:3) = 0;

middleY(1:Ny-1) = 0;

for i=1:Ny,
   middleY(i) = 0.5*(Yfile(i) + Yfile(i+1));
   Vx(i,1) = velocity((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + a*N, 1);
   Vx(i,2) = velocity((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + b*N, 1);
   Vx(i,3) = velocity((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + c*N, 1);
   Vy(i,1) = velocity((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + a*N, 2);
   Vy(i,2) = velocity((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + b*N, 2);
   Vy(i,3) = velocity((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + c*N, 2);
   Vz(i,1) = velocity((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + a*N, 3);
   Vz(i,2) = velocity((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + b*N, 3);
   Vz(i,3) = velocity((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + c*N, 3);
   
   Velectronx(i,1) = velocity_electron((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + a*N, 1);
   Velectronx(i,2) = velocity_electron((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + b*N, 1);
   Velectronx(i,3) = velocity_electron((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + c*N, 1);
   Velectrony(i,1) = velocity_electron((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + a*N, 2);
   Velectrony(i,2) = velocity_electron((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + b*N, 2);
   Velectrony(i,3) = velocity_electron((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + c*N, 2);
   Velectronz(i,1) = velocity_electron((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + a*N, 3);
   Velectronz(i,2) = velocity_electron((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + b*N, 3);
   Velectronz(i,3) = velocity_electron((Nz)*(Ny)*(xnumber-1) + (Nz)*(i-1) + znumber + c*N, 3);
end;
figure(1);
plot (middleY(1:Ny),Vx(1:Ny,1), 'red',middleY(1:Ny),Vx(1:Ny,2), 'green',middleY(1:Ny),Vx(1:Ny,3), 'blue');
title ('Vx');
xlabel ('y/r_g');
ylabel ('V cm/s');
grid ;

figure(2);
plot (middleY(1:Ny),Vy(1:Ny, 1), 'red', middleY(1:Ny), Vy(1:Ny, 2), 'green',middleY(1:Ny),Vy(1:Ny, 3), 'blue');
title ('Vy');
xlabel ('y/r_g');
ylabel ('V cm/s');
grid ;

figure(3);
plot (middleY(1:Ny),Vz(1:Ny, 1), 'red', middleY(1:Ny), Vz(1:Ny, 2), 'green', middleY(1:Ny), Vz(1:Ny, 3), 'blue');
title ('Vz');
xlabel ('y/r_g');
ylabel ('V cm/s');
grid ;

figure(4);
plot (middleY(1:Ny),Velectronx(1:Ny,1), 'red',middleY(1:Ny),Velectronx(1:Ny,2), 'green',middleY(1:Ny),Velectronx(1:Ny,3), 'blue');
title ('Vx');
xlabel ('y/r_g');
ylabel ('V cm/s');
grid ;

figure(5);
plot (middleY(1:Ny),Velectrony(1:Ny, 1), 'red', middleY(1:Ny), Velectrony(1:Ny, 2), 'green',middleY(1:Ny),Velectrony(1:Ny, 3), 'blue');
title ('Vy');
xlabel ('y/r_g');
ylabel ('V cm/s');
grid ;

figure(6);
plot (middleY(1:Ny),Velectronz(1:Ny, 1), 'red', middleY(1:Ny), Velectronz(1:Ny, 2), 'green', middleY(1:Ny), Velectronz(1:Ny, 3), 'blue');
title ('Vz');
xlabel ('y/r_g');
ylabel ('V cm/s');
grid ;

