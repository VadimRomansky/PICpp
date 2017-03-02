clear;
load concentrations.dat;
load particleTypes.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = Nx*Ny*Nz;

Nt = size(concentrations, 1)/N;
Ntypes = size(particleTypes, 1);

a = 0;
b = fix(Nt/2);
c = Nt - 1;

charge_density(1:Nx, 1:Ny) = 0;
charge_density_hat(1:Nx, 1:Ny) = 0;
particle_concentrations(1:Nx, 1:Ny, 1:Ntypes) = 0;

znumber = 1;

middleX(1:Nx) = 0;
middleY(1:Ny) = 0;

for i=1:Nx,  
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
   for j = 1:Ny,
       charge_density(i, j) = concentrations(Nz*Ny*(i-1) + Nz*(j-1) + znumber + c*N, 1);
       charge_density_hat(i, j) = concentrations(Nz*Ny*(i-1) + Nz*(j-1) + znumber + c*N, 2);
       for t = 1:Ntypes,
           particle_concentrations(i, j, t) = concentrations(Nz*Ny*(i-1) + Nz*(j-1) + znumber + c*N, 2 + t);
       end;
   end;
end;

for j = 1:Ny,
    middleY(j) = 0.5*(Yfile(j) + Yfile(j+1));
end;

cons(1:Nx, 1:Ny) = 0;
for t = 1:Ntypes,
    if(particleTypes(t) > 0)
        for i = 1:Nx,
            for j = 1:Ny,
                cons(i,j)=particle_concentrations(i, j, t);
            end;
        end;
        figure(t);
        [X, Y] = meshgrid(middleY, middleX);
        surf(X, Y, cons);
        title ('concentration');
        xlabel ('y');
        ylabel ('x');
        zlabel ('n cm^-3');
        grid ;
    end;
end;

figure(Ntypes + 1);
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, charge_density);
title ('charge density');
xlabel ('y');
ylabel ('x');
zlabel ('n cm^-3');
grid ;

figure(Ntypes + 2);
[X, Y] = meshgrid(middleY, middleX);
surf(X, Y, charge_density_hat);
title ('charge density hat');
xlabel ('y');
ylabel ('x');
zlabel ('n cm^-3');
grid ;