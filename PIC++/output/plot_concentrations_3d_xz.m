clear;
load particleTypes.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

concentratins = importdata('concentrations.dat');

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = Nx*Ny*Nz;

Nt = size(concentrations, 1)/N;
Ntypes = size(particleTypes, 1);

a = 0;
b = fix(Nt/2);
c = Nt - 1;

charge_density(1:Nx, 1:Nz) = 0;
charge_density_hat(1:Nx, 1:Nz) = 0;
particle_concentrations(1:Nx, 1:Nz, 1:Ntypes) = 0;

ynumber = 2;

middleX(1:Nx) = 0;
middleZ(1:Nz) = 0;

for i=1:Nx,  
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
   for j = 1:Nz,
       charge_density(i, j) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + j + c*N, 1);
       charge_density_hat(i, j) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + j + c*N, 2);
       for t = 1:Ntypes,
           particle_concentrations(i, j, t) = concentrations(Nz*Ny*(i-1) + Nz*(ynumber-1) + j + c*N, 2 + t);
       end;
   end;
end;

for j = 1:Nz,
    middleZ(j) = 0.5*(Zfile(j) + Zfile(j+1));
end;

set(0, 'DefaultLineLineWidth', 2);

cons(1:Nx, 1:Nz) = 0;
for t = 1:Ntypes,
    if(particleTypes(t) > 0)
        for i = 1:Nx,
            for j = 1:Nz,
                cons(i,j)=particle_concentrations(i, j, t);
            end;
        end;
        figure(t);
        colormap Jet;
        [X, Z] = meshgrid(middleZ, middleX);
        surf(X, Z, cons);
        shading interp;
        title ('concentration');
        xlabel ('z');
        ylabel ('y');
        zlabel ('n cm^-3');
        grid ;
    end;
end;

figure(Ntypes + 1);
colormap Jet;
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, charge_density);
shading interp;
title ('charge density');
xlabel ('z');
ylabel ('x');
zlabel ('n cm^-3');
grid ;

figure(Ntypes + 2);
colormap Jet;
[X, Z] = meshgrid(middleZ, middleX);
surf(X, Z, charge_density_hat);
shading interp;
title ('charge density hat');
xlabel ('z');
ylabel ('y');
zlabel ('n cm^-3');
grid ;