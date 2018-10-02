clear;
load velocityX.dat;
load particleTypes.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
set(0, 'DefaultLineLineWidth', 2);
Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = Nx;
Nt = size(velocityX, 1)/N;
Ntypes = size(particleTypes,1);

a = 0;
b = fix(Nt/2);
c = Nt - 1;

particleVelocitiesX(1:Nx, 1:(3*Ntypes)) = 0;
particleVelocitiesY(1:Nx, 1:(3*Ntypes)) = 0;
particleVelocitiesZ(1:Nx, 1:(3*Ntypes)) = 0;

middleX(1:Nx) = 0;

for i=1:Nx,
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
   for t = 1:Ntypes,
       particleVelocitiesX(i, 1 + 3*(t-1)) = velocityX(i + a*N, 1 + 3*(t-1));
       particleVelocitiesX(i, 2 + 3*(t-1)) = velocityX(i + b*N, 1 + 3*(t-1));
       particleVelocitiesX(i, 3 + 3*(t-1)) = velocityX(i + c*N, 1 + 3*(t-1));

       particleVelocitiesY(i, 1 + 3*(t-1)) = velocityX(i + a*N, 2 + 3*(t-1));
       particleVelocitiesY(i, 2 + 3*(t-1)) = velocityX(i + b*N, 2 + 3*(t-1));
       particleVelocitiesY(i, 3 + 3*(t-1)) = velocityX(i + c*N, 2 + 3*(t-1));
       
       particleVelocitiesZ(i, 1 + 3*(t-1)) = velocityX(i + a*N, 3 + 3*(t-1));
       particleVelocitiesZ(i, 2 + 3*(t-1)) = velocityX(i + b*N, 3 + 3*(t-1));
       particleVelocitiesZ(i, 3 + 3*(t-1)) = velocityX(i + c*N, 3 + 3*(t-1));
   end;
end;

for t = 1:Ntypes,
    if (particleTypes(t) > 0)
    figure(1 + 3*(t-1));
    plot (Xfile(1:Nx,1),particleVelocitiesX(1:Nx, 1 + 3*(t-1)), 'red', Xfile(1:Nx,1), particleVelocitiesX(1:Nx, 2 + 3*(t-1)), 'green', Xfile(1:Nx,1), particleVelocitiesX(1:Nx, 3 + 3*(t-1)), 'blue');
    title ('velocity x');
    xlabel ('x/r_g');
    ylabel ('V cm/s');
    grid ;    
    
    figure(2 + 3*(t-1));
    plot (Xfile(1:Nx,1),particleVelocitiesY(1:Nx, 1 + 3*(t-1)), 'red', Xfile(1:Nx,1), particleVelocitiesY(1:Nx, 2 + 3*(t-1)), 'green', Xfile(1:Nx,1), particleVelocitiesY(1:Nx, 3 + 3*(t-1)), 'blue');
    title ('velocity y');
    xlabel ('x/r_g');
    ylabel ('V cm/s');
    grid ;
    
    figure(3 + 3*(t-1));
    plot (Xfile(1:Nx,1),particleVelocitiesZ(1:Nx, 1 + 3*(t-1)), 'red', Xfile(1:Nx,1), particleVelocitiesZ(1:Nx, 2 + 3*(t-1)), 'green', Xfile(1:Nx,1), particleVelocitiesZ(1:Nx, 3 + 3*(t-1)), 'blue');
    title ('velocity z');
    xlabel ('x/r_g');
    ylabel ('V cm/s');
    grid ; 
    end;
end;