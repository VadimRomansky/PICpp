clear;
load velocity.dat;
load particleTypes.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
set(0, 'DefaultLineLineWidth', 2);

Nx = size(Xfile, 1)-1;
Ny = size(Yfile, 1)-1;
Nz = size(Zfile, 1)-1;

N = Nx*Ny*Nz;
Nt = size(velocity, 1)/N;
Ntypes = size(particleTypes,1);
xnumber = 10;
znumber = 2;

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
       particleVelocitiesX(i, 1 + 3*(t-1)) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber + a*N, 1 + 3*(t-1));
       particleVelocitiesX(i, 2 + 3*(t-1)) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber + b*N, 1 + 3*(t-1));
       particleVelocitiesX(i, 3 + 3*(t-1)) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber + c*N, 1 + 3*(t-1));

       particleVelocitiesY(i, 1 + 3*(t-1)) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber + a*N, 2 + 3*(t-1));
       particleVelocitiesY(i, 2 + 3*(t-1)) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber + b*N, 2 + 3*(t-1));
       particleVelocitiesY(i, 3 + 3*(t-1)) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber + c*N, 2 + 3*(t-1));
       
       particleVelocitiesZ(i, 1 + 3*(t-1)) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber + a*N, 3 + 3*(t-1));
       particleVelocitiesZ(i, 2 + 3*(t-1)) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber + b*N, 3 + 3*(t-1));
       particleVelocitiesZ(i, 3 + 3*(t-1)) = velocity(Nz*Ny*(xnumber-1) + Nz*(i-1) + znumber + c*N, 3 + 3*(t-1));
   end;
end;

for t = 1:Ntypes,
    if (particleTypes(t) > 0)
    figure(1 + 3*(t-1));
    plot (middleY(1:Ny),particleVelocitiesX(1:Ny, 1 + 3*(t-1)), 'red', middleY(1:Ny), particleVelocitiesX(1:Ny, 2 + 3*(t-1)), 'green', middleY(1:Ny), particleVelocitiesX(1:Ny, 3 + 3*(t-1)), 'blue');
    title ('velocity x');
    xlabel ('x/r_g');
    ylabel ('V cm/s');
    grid ;    
    
    figure(2 + 3*(t-1));
    plot (middleY(1:Ny),particleVelocitiesY(1:Ny, 1 + 3*(t-1)), 'red', middleY(1:Ny), particleVelocitiesY(1:Ny, 2 + 3*(t-1)), 'green', middleY(1:Ny), particleVelocitiesY(1:Ny, 3 + 3*(t-1)), 'blue');
    title ('velocity y');
    xlabel ('x/r_g');
    ylabel ('V cm/s');
    grid ;
    
    figure(3 + 3*(t-1));
    plot (middleY(1:Ny),particleVelocitiesZ(1:Ny, 1 + 3*(t-1)), 'red', middleY(1:Ny), particleVelocitiesZ(1:Ny, 2 + 3*(t-1)), 'green', middleY(1:Ny), particleVelocitiesZ(1:Ny, 3 + 3*(t-1)), 'blue');
    title ('velocity z');
    xlabel ('x/r_g');
    ylabel ('V cm/s');
    grid ; 
    end;
end;