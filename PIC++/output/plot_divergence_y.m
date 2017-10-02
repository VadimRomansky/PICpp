clear;
%load concentrations.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
load divergence_error.dat;

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = (Nx)*(Ny)*(Nz);
Nt = fix(size(divergence_error, 1)/N)-1;

xnumber = 3;
znumber = 2;

a = 0;
b = fix(Nt/2);
c = Nt;

divergenceError(1:Ny, 1:3) = 0;
divergence(1:Ny, 1:3) = 0;
density(1:Ny, 1:3) = 0;
divergenceB(1:Ny, 1:3) = 0;

middleY(1:Ny)=0;


for i=1:Ny,  
   middleY(i) = (Yfile(i)+Yfile(i+1))/2;
   divergenceError(i, 1) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + a*N, 1);
   divergenceError(i, 2) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + b*N, 1);
   divergenceError(i, 3) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + c*N, 1);
   divergence(i, 1) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + a*N, 2);
   divergence(i, 2) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + b*N, 2);
   divergence(i, 3) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + c*N, 2);
   density(i,1) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + a*N, 3);
   density(i,2) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + b*N, 3);
   density(i,3) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + c*N, 3);
   divergenceB(i, 1) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + a*N, 4);
   divergenceB(i, 2) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + b*N, 4);
   divergenceB(i, 3) = divergence_error((Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber) + c*N, 4);
end;

figure(1);
plot (middleY(1:Ny),divergenceError(1:Ny, 1), 'red', middleY(1:Ny), divergenceError(1:Ny, 2), 'green', middleY(1:Ny), divergenceError(1:Ny, 3), 'blue');
title ('divergence error');
xlabel ('y cm');
ylabel ('rho sgs*cm^-3');
grid ;

figure(2);
plot (middleY(1:Ny),divergence(1:Ny, 1), 'red', middleY(1:Ny), divergence(1:Ny, 2), 'green', middleY(1:Ny), divergence(1:Ny, 3), 'blue');
title ('divergence');
xlabel ('y cm');
ylabel ('rho sgs*cm^-3');
grid ;

figure(3);
plot (middleY(1:Ny),density(1:Ny, 1), 'red', middleY(1:Ny), density(1:Ny, 2), 'green', middleY(1:Ny), density(1:Ny, 3), 'blue');
title ('4*pi*density');
xlabel ('y cm');
ylabel ('rho sgs*cm^-3');
grid ;

figure(4);
plot (middleY(1:Ny),divergenceB(1:Ny, 1), 'red', middleY(1:Ny), divergenceB(1:Ny, 2), 'green', middleY(1:Ny), divergenceB(1:Ny, 3), 'blue');
title ('div B');
xlabel ('y cm');
ylabel ('rho sgs*cm^-3');
grid ;