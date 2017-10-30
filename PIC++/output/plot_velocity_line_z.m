clear;
load velocityZ.dat;
load particleTypes.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = Nz;
Nt = size(velocityZ, 1)/N;
Ntypes = size(particleTypes,1);

a = 0;
b = fix(Nt/2);
c = Nt - 1;

particleVelocitiesX(1:Nz, 1:(3*Ntypes)) = 0;
particleVelocitiesY(1:Nz, 1:(3*Ntypes)) = 0;
particleVelocitiesZ(1:Nz, 1:(3*Ntypes)) = 0;

middleZ(1:Nz) = 0;

for i=1:Nz,
   middleZ(i) = 0.5*(Zfile(i) + Zfile(i+1));
   for t = 1:Ntypes,
       particleVelocitiesX(i, 1 + 3*(t-1)) = velocityZ(i + a*N, 1 + 3*(t-1));
       particleVelocitiesX(i, 2 + 3*(t-1)) = velocityZ(i + b*N, 1 + 3*(t-1));
       particleVelocitiesX(i, 3 + 3*(t-1)) = velocityZ(i + c*N, 1 + 3*(t-1));

       particleVelocitiesY(i, 1 + 3*(t-1)) = velocityZ(i + a*N, 2 + 3*(t-1));
       particleVelocitiesY(i, 2 + 3*(t-1)) = velocityZ(i + b*N, 2 + 3*(t-1));
       particleVelocitiesY(i, 3 + 3*(t-1)) = velocityZ(i + c*N, 2 + 3*(t-1));
       
       particleVelocitiesZ(i, 1 + 3*(t-1)) = velocityZ(i + a*N, 3 + 3*(t-1));
       particleVelocitiesZ(i, 2 + 3*(t-1)) = velocityZ(i + b*N, 3 + 3*(t-1));
       particleVelocitiesZ(i, 3 + 3*(t-1)) = velocityZ(i + c*N, 3 + 3*(t-1));
   end;
end;

for t = 1:Ntypes,
    if (particleTypes(t) > 0)
    figure(1 + 3*(t-1));
    plot (Zfile(1:Nz,1),particleVelocitiesX(1:Nz, 1 + 3*(t-1)), 'red', Zfile(1:Nz,1), particleVelocitiesX(1:Nz, 2 + 3*(t-1)), 'green', Zfile(1:Nz,1), particleVelocitiesX(1:Nz, 3 + 3*(t-1)), 'blue');
    title ('velocity x');
    xlabel ('z/r_g');
    ylabel ('V cm/s');
    grid ;    
    
    figure(2 + 3*(t-1));
    plot (Zfile(1:Nz,1),particleVelocitiesY(1:Nz, 1 + 3*(t-1)), 'red', Zfile(1:Nz,1), particleVelocitiesY(1:Nz, 2 + 3*(t-1)), 'green', Zfile(1:Nz,1), particleVelocitiesY(1:Nz, 3 + 3*(t-1)), 'blue');
    title ('velocity y');
    xlabel ('z/r_g');
    ylabel ('V cm/s');
    grid ;
    
    figure(3 + 3*(t-1));
    plot (Zfile(1:Nz,1),particleVelocitiesZ(1:Nz, 1 + 3*(t-1)), 'red', Zfile(1:Nz,1), particleVelocitiesZ(1:Nz, 2 + 3*(t-1)), 'green', Zfile(1:Nz,1), particleVelocitiesZ(1:Nz, 3 + 3*(t-1)), 'blue');
    title ('velocity z');
    xlabel ('x/r_g');
    ylabel ('V cm/s');
    grid ; 
    end;
end;