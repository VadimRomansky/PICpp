clear;
load Xfile.dat;
load divergence_error.dat;

Nx = size(Xfile, 1) - 1;

Nt = size(divergence_error, 1)/Nx;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

divergenceError(1:Nx, 1:3) = 0;
divergence(1:Nx, 1:3) = 0;
density(1:Nx, 1:3) = 0;


for i=1:Nx,   
   divergenceError(i, 1) = divergence_error(i + a*Nx, 1);
   divergenceError(i, 2) = divergence_error(i + b*Nx, 1);
   divergenceError(i, 3) = divergence_error(i + c*Nx, 1);
   divergence(i, 1) = divergence_error(i + a*Nx, 2);
   divergence(i, 2) = divergence_error(i + b*Nx, 2);
   divergence(i, 3) = divergence_error(i + c*Nx, 2);
   density(i,1) = divergence_error(i + a*Nx, 3);
   density(i,2) = divergence_error(i + b*Nx, 3);
   density(i,3) = divergence_error(i + c*Nx, 3);
end;

figure(1);
plot (Xfile(1:Nx,1),divergenceError(1:Nx, 1), 'red', Xfile(1:Nx,1), divergenceError(1:Nx, 2), 'green', Xfile(1:Nx,1), divergenceError(1:Nx, 3), 'blue');
title ('divergence error');
xlabel ('x cm');
ylabel ('rho sgs*cm^-3');
grid ;

figure(2);
plot (Xfile(1:Nx,1),divergence(1:Nx, 1), 'red', Xfile(1:Nx,1), divergence(1:Nx, 2), 'green', Xfile(1:Nx,1), divergence(1:Nx, 3), 'blue');
title ('divergence');
xlabel ('x cm');
ylabel ('rho sgs*cm^-3');
grid ;

figure(3);
plot (Xfile(1:Nx,1),density(1:Nx, 1), 'red', Xfile(1:Nx,1), density(1:Nx, 2), 'green', Xfile(1:Nx,1), density(1:Nx, 3), 'blue');
title ('4*pi*density');
xlabel ('x cm');
ylabel ('rho sgs*cm^-3');
grid ;