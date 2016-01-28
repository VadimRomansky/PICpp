clear;
load concentrations.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
load divergence_error.dat;

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = (Nx)*(Ny)*(Nz);
Nt = fix(size(concentrations, 1)/N)-1;

ynumber = 1;
znumber = 1;

a = 0;
b = fix(Nt/2);
c = Nt;

divergenceError(1:Nx, 1:3) = 0;
divergence(1:Nx, 1:3) = 0;
density(1:Nx, 1:3) = 0;


for i=1:Nx,   
   divergenceError(i, 1) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 1);
   divergenceError(i, 2) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 1);
   divergenceError(i, 3) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 1);
   divergence(i, 1) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 2);
   divergence(i, 2) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 2);
   divergence(i, 3) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 2);
   density(i,1) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 3);
   density(i,2) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 3);
   density(i,3) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 3);
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