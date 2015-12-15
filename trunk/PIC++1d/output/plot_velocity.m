clear;
load velocity.dat;
load velocity_electron.dat;
load Xfile.dat;

Nx = size(Xfile, 1) - 1;

N = Nx;
Nt = size(velocity, 1)/N;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

Vx(1:Nx, 1:3) = 0;
Vy(1:Nx, 1:3) = 0;
Vz(1:Nx, 1:3) = 0;

Velectronx(1:Nx, 1:3) = 0;
Velectrony(1:Nx, 1:3) = 0;
Velectronz(1:Nx, 1:3) = 0;

for i=1:Nx,
   Vx(i,1) = velocity(i + a*N, 1);
   Vx(i,2) = velocity(i + b*N, 1);
   Vx(i,3) = velocity(i + c*N, 1);
   Vy(i,1) = velocity(i + a*N, 2);
   Vy(i,2) = velocity(i + b*N, 2);
   Vy(i,3) = velocity(i + c*N, 2);
   Vz(i,1) = velocity(i + a*N, 3);
   Vz(i,2) = velocity(i + b*N, 3);
   Vz(i,3) = velocity(i + c*N, 3);
   
   Velectronx(i,1) = velocity_electron(i + a*N, 1);
   Velectronx(i,2) = velocity_electron(i + b*N, 1);
   Velectronx(i,3) = velocity_electron(i + c*N, 1);
   Velectrony(i,1) = velocity_electron(i + a*N, 2);
   Velectrony(i,2) = velocity_electron(i + b*N, 2);
   Velectrony(i,3) = velocity_electron(i + c*N, 2);
   Velectronz(i,1) = velocity_electron(i + a*N, 3);
   Velectronz(i,2) = velocity_electron(i + b*N, 3);
   Velectronz(i,3) = velocity_electron(i + c*N, 3);
end;
figure(1);
plot (Xfile(1:Nx,1),Vx(1:Nx,1), 'red',Xfile(1:Nx,1),Vx(1:Nx,2), 'green',Xfile(1:Nx,1),Vx(1:Nx,3), 'blue');
title ('Vx');
xlabel ('x cm');
ylabel ('V cm/s');
legend('at t = 0', 'at t = maxT/2','at t = maxT','Location','southwest');
grid ;

figure(2);
plot (Xfile(1:Nx,1),Vy(1:Nx, 1), 'red', Xfile(1:Nx,1), Vy(1:Nx, 2), 'green',Xfile(1:Nx,1),Vy(1:Nx, 3), 'blue');
title ('Vy');
xlabel ('x cm');
ylabel ('V cm/s');
legend('at t = 0', 'at t = maxT/2','at t = maxT','Location','southwest');
grid ;

figure(3);
plot (Xfile(1:Nx,1),Vz(1:Nx, 1), 'red', Xfile(1:Nx,1), Vz(1:Nx, 2), 'green', Xfile(1:Nx,1), Vz(1:Nx, 3), 'blue');
title ('Vz');
xlabel ('x cm');
ylabel ('V cm/s');
legend('at t = 0', 'at t = maxT/2','at t = maxT','Location','southwest');
grid ;

figure(4);
plot (Xfile(1:Nx,1),Velectronx(1:Nx,1), 'red',Xfile(1:Nx,1),Velectronx(1:Nx,2), 'green',Xfile(1:Nx,1),Velectronx(1:Nx,3), 'blue');
title ('Vx');
xlabel ('x cm');
ylabel ('V cm/s');
legend('at t = 0', 'at t = maxT/2','at t = maxT','Location','southwest');
grid ;

figure(5);
plot (Xfile(1:Nx,1),Velectrony(1:Nx, 1), 'red', Xfile(1:Nx,1), Velectrony(1:Nx, 2), 'green',Xfile(1:Nx,1),Velectrony(1:Nx, 3), 'blue');
title ('Vy');
xlabel ('x cm');
ylabel ('V cm/s');
legend('at t = 0', 'at t = maxT/2','at t = maxT','Location','southwest');
grid ;

figure(6);
plot (Xfile(1:Nx,1),Velectronz(1:Nx, 1), 'red', Xfile(1:Nx,1), Velectronz(1:Nx, 2), 'green', Xfile(1:Nx,1), Velectronz(1:Nx, 3), 'blue');
title ('Vz');
xlabel ('x cm');
ylabel ('V cm/s');
legend('at t = 0', 'at t = maxT/2','at t = maxT','Location','southwest');
grid ;

