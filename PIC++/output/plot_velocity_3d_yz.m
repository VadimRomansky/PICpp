clear;
load velocity.dat;
load particleTypes.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
set(0,'DefaultFigureColormap',feval('jet'));
Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = Nx*Ny*Nz;

Nt = size(velocity, 1)/N;
Ntypes = size(particleTypes,1);

znumber = 1;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

particleVelocitiesX(1:Ny, 1:Nz, 1:Ntypes) = 0;
particleVelocitiesY(1:Ny, 1:Nz, 1:Ntypes) = 0;
particleVelocitiesZ(1:Ny, 1:Nz, 1:Ntypes) = 0;

xnumber = fix(Nx/2);

middleZ(1:Nz) = 0;
middleY(1:Ny) = 0;

for i=1:Ny,
   middleY(i) = 0.5*(Yfile(i) + Yfile(i+1));
   for j = 1:Nz,
     for t = 1:Ntypes,
       particleVelocitiesX(i, j, t) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + j + c*N, 1 + 3*(t-1));

       particleVelocitiesY(i, j, t) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + j + c*N, 2 + 3*(t-1));
       
       particleVelocitiesZ(i, j, t) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + j + c*N, 3 + 3*(t-1));
     end;
   end;
end;

for j = 1:Nz,
    middleZ(j) = 0.5*(Zfile(j) + Zfile(j+1));
end;

vel(1:Ny,1:Nz) = 0;
for t = 1:Ntypes,
    if (particleTypes(t) > 0)
        for i = 1:Ny,
            for j = 1:Nz,
                vel(i,j)=particleVelocitiesX(i, j, t);
            end;
        end;
        figure(1 + 3*(t-1));
        [Y, Z] = meshgrid(middleZ, middleY);
        surf(Y, Z, vel);
        shading interp;
        title ('velocity x');
        xlabel ('z');
        ylabel ('y');
        zlabel ('v cm/s');
        grid ;  
    
        for i = 1:Ny,
            for j = 1:Nz,
                vel(i,j)=particleVelocitiesY(i, j, t);
            end;
        end;
        figure(2 + 3*(t-1));
        [Y, Z] = meshgrid(middleZ, middleY);
        surf(Y, Z, vel);
        shading interp;
        title ('velocity y');
        xlabel ('z');
        ylabel ('y');
        zlabel ('v cm/s');
        grid ; 
    
        for i = 1:Ny,
            for j = 1:Nz,
                vel(i,j)=particleVelocitiesZ(i, j, t);
            end;
        end;
        figure(3 + 3*(t-1));
        [Y, Z] = meshgrid(middleZ, middleY);
        surf(Y, Z, vel);
        shading interp;
        title ('velocity z');
        xlabel ('z');
        ylabel ('y');
        zlabel ('v cm/s');
        grid ; 
    end;
end;