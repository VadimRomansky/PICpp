clear;
load concentrations.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
load divergence_error.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

N = Nx*Ny*Nz;
Nt = size(concentrations, 1)/N;

ynumber = 1;
znumber = 1;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

divergence(1:Nx, 1:3) = 0;
divergenceError1(1:Nx, 1:3) = 0;
divergenceError2(1:Nx, 1:3) = 0;


for i=1:Nx,   
   divergence(i, 1) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 1);
   divergence(i, 2) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 1);
   divergence(i, 3) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 1);
   divergenceError1(i, 1) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 2);
   divergenceError1(i, 2) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 2);
   divergenceError1(i, 3) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 2);
   divergenceError2(i, 1) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 3);
   divergenceError2(i, 2) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 3);
   divergenceError2(i, 3) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 3);
end;

figure(1);
plot (Xfile(1:Nx,1),divergence(1:Nx, 1), 'red', Xfile(1:Nx,1), divergence(1:Nx, 2), 'green', Xfile(1:Nx,1), divergence(1:Nx, 3), 'blue');
title ('divergence');
xlabel ('x/r_g');
ylabel ('rho sgs*cm^-3');
grid ;

figure(2);
plot (Xfile(1:Nx,1),divergenceError1(1:Nx, 1), 'red', Xfile(1:Nx,1), divergenceError1(1:Nx, 2), 'green', Xfile(1:Nx,1), divergenceError1(1:Nx, 3), 'blue');
title ('divergence error');
xlabel ('x/r_g');
ylabel ('rho sgs*cm^-3');
grid ;

figure(3);
plot (Xfile(1:Nx,1),divergenceError2(1:Nx, 1), 'red', Xfile(1:Nx,1), divergenceError2(1:Nx, 2), 'green', Xfile(1:Nx,1), divergenceError2(1:Nx, 3), 'blue');
title ('divergence error');
xlabel ('x/r_g');
ylabel ('rho sgs*cm^-3');
grid ;