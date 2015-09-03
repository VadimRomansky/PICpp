clear;
load Xfile.dat;
load divergence_error.dat;

Nx = size(Xfile, 1);

Nt = size(divergence_error, 1)/Nx;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

divergence(1:Nx, 1:3) = 0;
divergence2(1:Nx, 1:3) = 0;


for i=1:Nx,   
   divergence(i, 1) = divergence_error(i + a*Nx, 1);
   divergence(i, 2) = divergence_error(i + b*Nx, 1);
   divergence(i, 3) = divergence_error(i + c*Nx, 1);
   divergence2(i, 1) = divergence_error(i + a*Nx, 2);
   divergence2(i, 2) = divergence_error(i + b*Nx, 2);
   divergence2(i, 3) = divergence_error(i + c*Nx, 2);
end;

figure(1);
plot (Xfile(1:Nx,1),divergence(1:Nx, 1), 'red', Xfile(1:Nx,1), divergence(1:Nx, 2), 'green', Xfile(1:Nx,1), divergence(1:Nx, 3), 'blue');
title ('divergence');
xlabel ('x/r_g');
ylabel ('rho sgs*cm^-3');
grid ;

figure(2);
plot (Xfile(1:Nx,1),divergence2(1:Nx, 1), 'red', Xfile(1:Nx,1), divergence2(1:Nx, 2), 'green', Xfile(1:Nx,1), divergence2(1:Nx, 3), 'blue');
title ('divergence');
xlabel ('x/r_g');
ylabel ('rho sgs*cm^-3');
grid ;