clear;
load velocityY.dat;
load particleTypes.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = Ny;
Nt = size(velocityY, 1)/N;
Ntypes = size(particleTypes,1);

a = 0;
b = fix(Nt/2);
c = Nt - 1;

particleVelocitiesX(1:Ny, 1:(3*Ntypes)) = 0;
particleVelocitiesY(1:Ny, 1:(3*Ntypes)) = 0;
particleVelocitiesZ(1:Ny, 1:(3*Ntypes)) = 0;

middleY(1:Ny) = 0;

for i=1:Ny,
   middleY(i) = 0.5*(Yfile(i) + Yfile(i+1));
   for t = 1:Ntypes,
       particleVelocitiesX(i, 1 + 3*(t-1)) = velocityY(i + a*N, 1 + 3*(t-1));
       particleVelocitiesX(i, 2 + 3*(t-1)) = velocityY(i + b*N, 1 + 3*(t-1));
       particleVelocitiesX(i, 3 + 3*(t-1)) = velocityY(i + c*N, 1 + 3*(t-1));

       particleVelocitiesY(i, 1 + 3*(t-1)) = velocityY(i + a*N, 2 + 3*(t-1));
       particleVelocitiesY(i, 2 + 3*(t-1)) = velocityY(i + b*N, 2 + 3*(t-1));
       particleVelocitiesY(i, 3 + 3*(t-1)) = velocityY(i + c*N, 2 + 3*(t-1));
       
       particleVelocitiesZ(i, 1 + 3*(t-1)) = velocityY(i + a*N, 3 + 3*(t-1));
       particleVelocitiesZ(i, 2 + 3*(t-1)) = velocityY(i + b*N, 3 + 3*(t-1));
       particleVelocitiesZ(i, 3 + 3*(t-1)) = velocityY(i + c*N, 3 + 3*(t-1));
   end;
end;

for t = 1:Ntypes,
    if (particleTypes(t) > 0)
    figure(1 + 3*(t-1));
    plot (Yfile(1:Ny,1),particleVelocitiesX(1:Ny, 1 + 3*(t-1)), 'red', Yfile(1:Ny,1), particleVelocitiesX(1:Ny, 2 + 3*(t-1)), 'green', Yfile(1:Ny,1), particleVelocitiesX(1:Ny, 3 + 3*(t-1)), 'blue');
    title ('velocity x');
    xlabel ('y/r_g');
    ylabel ('V cm/s');
    grid ;    
    
    figure(2 + 3*(t-1));
    plot (Yfile(1:Ny,1),particleVelocitiesY(1:Ny, 1 + 3*(t-1)), 'red', Yfile(1:Ny,1), particleVelocitiesY(1:Ny, 2 + 3*(t-1)), 'green', Yfile(1:Ny,1), particleVelocitiesY(1:Ny, 3 + 3*(t-1)), 'blue');
    title ('velocity y');
    xlabel ('y/r_g');
    ylabel ('V cm/s');
    grid ;
    
    figure(3 + 3*(t-1));
    plot (Yfile(1:Ny,1),particleVelocitiesZ(1:Ny, 1 + 3*(t-1)), 'red', Yfile(1:Ny,1), particleVelocitiesZ(1:Ny, 2 + 3*(t-1)), 'green', Yfile(1:Ny,1), particleVelocitiesZ(1:Ny, 3 + 3*(t-1)), 'blue');
    title ('velocity z');
    xlabel ('y/r_g');
    ylabel ('V cm/s');
    grid ; 
    end;
end;