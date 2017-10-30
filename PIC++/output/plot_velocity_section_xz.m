clear;
load velocityXZ.dat;
load particleTypes.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = Nx*Nz;

Nt = size(velocityXZ, 1)/N;
Ntypes = size(particleTypes,1);

a = 0;
b = fix(Nt/2);
c = Nt - 1;

particleVelocitiesX(1:Nx, 1:Nz, 1:Ntypes) = 0;
particleVelocitiesY(1:Nx, 1:Nz, 1:Ntypes) = 0;
particleVelocitiesZ(1:Nx, 1:Nz, 1:Ntypes) = 0;

middleX(1:Nx) = 0;
middleZ(1:Nz) = 0;

for i=1:Nx,
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
   for j = 1:Nz,
     for t = 1:Ntypes,
       particleVelocitiesX(i, j, t) = velocityXZ(Nz*(i-1) + j + c*N, 1 + 3*(t-1));

       particleVelocitiesY(i, j, t) = velocityXZ(Nz*(i-1) + j + c*N, 2 + 3*(t-1));
       
       particleVelocitiesZ(i, j, t) = velocityXZ(Nz*(i-1) + j + c*N, 3 + 3*(t-1));
     end;
   end;
end;

for j = 1:Nz,
    middleZ(j) = 0.5*(Zfile(j) + Zfile(j+1));
end;

vel(1:Nx,1:Nz) = 0;
for t = 1:Ntypes,
    if (particleTypes(t) > 0)
        for i = 1:Nx,
            for j = 1:Nz,
                vel(i,j)=particleVelocitiesX(i, j, t);
            end;
        end;
        figure(1 + 3*(t-1));
        [X, Z] = meshgrid(middleZ, middleX);
        surf(X, Z, vel);
        shading interp;
        title ('velocity x');
        xlabel ('z');
        ylabel ('x');
        zlabel ('v cm/s');
        grid ;  
    
        for i = 1:Nx,
            for j = 1:Nz,
                vel(i,j)=particleVelocitiesY(i, j, t);
            end;
        end;
        figure(2 + 3*(t-1));
        [X, Z] = meshgrid(middleZ, middleX);
        surf(X, Z, vel);
        shading interp;
        title ('velocity y');
        xlabel ('z');
        ylabel ('x');
        zlabel ('v cm/s');
        grid ; 
    
        for i = 1:Nx,
            for j = 1:Nz,
                vel(i,j)=particleVelocitiesZ(i, j, t);
            end;
        end;
        figure(3 + 3*(t-1));
        [X, Z] = meshgrid(middleZ, middleX);
        surf(X, Z, vel);
        shading interp;
        title ('velocity z');
        xlabel ('z');
        ylabel ('x');
        zlabel ('v cm/s');
        grid ; 
    end;
end;