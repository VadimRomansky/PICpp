clear;
load anisotropy_protons.dat;
load anisotropy_electrons.dat;
load anisotropy_alphas.dat;
load anisotropy_positrons.dat;
load anisotropy_deuterium.dat;
load anisotropy_helium3.dat;


load protons.dat;
load electrons.dat;
load positrons.dat;
load alphas.dat;
load deuterium.dat;
load helium3.dat;

load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

N2=size(protons,1);
N1=size(electrons,1);
N3=size(positrons,1);
N4=size(alphas,1);
N5 = size(deuterium,1);
N6 = size(helium3, 1);

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NB = (Nx-1)*(Ny-1)*(Nz-1);
Nt = (size(anisotropy_protons, 1)/NB);

a = 0;
b = fix(Nt/2);
c = Nt - 1;

ynumber = 1;
znumber = 1;

electrons(1:Nx, 1:3) = 0;
protons(1:Nx, 1:3) = 0;
positrons(1:Nx, 1:3) = 0;
alphas(1:Nx, 1:3) = 0;
deuterium(1:Nx, 1:3) = 0;
helium3(1:Nx, 1:3) = 0;

middleX(1:Nx-1) = 0;

for i=1:Nx-1,
    middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
    electrons(i, 1) = anisotropy_electrons(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 1);
    electrons(i, 2) = anisotropy_electrons(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 1);
    electrons(i, 3) = anisotropy_electrons(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 1);
    
    protons(i, 1) = anisotropy_protons(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 1);
    protons(i, 2) = anisotropy_protons(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 1);
    protons(i, 3) = anisotropy_protons(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 1);
    
    positrons(i, 1) = anisotropy_positrons(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 1);
    positrons(i, 2) = anisotropy_positrons(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 1);
    positrons(i, 3) = anisotropy_positrons(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 1);
    
    alphas(i, 1) = anisotropy_alphas(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 1);
    alphas(i, 2) = anisotropy_alphas(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 1);
    alphas(i, 3) = anisotropy_alphas(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 1);
    
    deuterium(i, 1) = anisotropy_deuterium(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 1);
    deuterium(i, 2) = anisotropy_deuterium(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 1);
    deuterium(i, 3) = anisotropy_deuterium(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 1);
    
    helium3(i, 1) = anisotropy_helium3(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + a*NB, 1);
    helium3(i, 2) = anisotropy_helium3(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + b*NB, 1);
    helium3(i, 3) = anisotropy_helium3(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber) + c*NB, 1);    
end;

if(N1 > 0)
    figure(1);
    plot (middleX(1:Nx-1),electrons(1:Nx-1, 1), 'red', middleX(1:Nx-1),electrons(1:Nx-1, 2), 'green', middleX(1:Nx-1),electrons(1:Nx-1, 3), 'blue');
    title ('anisotropy electrons');
    xlabel ('x');
    ylabel ('T1/T2 - T2/T1');
    grid ;
end;

if(N2 > 0)
    figure(2);
    plot (middleX(1:Nx-1),protons(1:Nx-1, 1), 'red', middleX(1:Nx-1),protons(1:Nx-1, 2), 'green', middleX(1:Nx-1),protons(1:Nx-1, 3), 'blue');
    title ('anisotropy protons');
    xlabel ('x');
    ylabel ('T1/T2 - T2/T1');
    grid ;
end;

if(N3 > 0)
    figure(3);
    plot (middleX(1:Nx-1),positrons(1:Nx-1, 1), 'red', middleX(1:Nx-1),positrons(1:Nx-1, 2), 'green', middleX(1:Nx-1),positronss(1:Nx-1, 3), 'blue');
    title ('anisotropy positrons');
    xlabel ('x');
    ylabel ('T1/T2 - T2/T1');
    grid ;
end;

if(N4 > 0)
    figure(4);
    plot (middleX(1:Nx-1),alphas(1:Nx-1, 1), 'red', middleX(1:Nx-1),alphas(1:Nx-1, 2), 'green', middleX(1:Nx-1),alphas(1:Nx-1, 3), 'blue');
    title ('anisotropy alphas');
    xlabel ('x');
    ylabel ('T1/T2 - T2/T1');
    grid ;
end;

if(N5 > 0)
    figure(5);
    plot (middleX(1:Nx-1),deuterium(1:Nx-1, 1), 'red', middleX(1:Nx-1),deuterium(1:Nx-1, 2), 'green', middleX(1:Nx-1),deuterium(1:Nx-1, 3), 'blue');
    title ('anisotropy deuterium');
    xlabel ('x');
    ylabel ('T1/T2 - T2/T1');
    grid ;
end;

if(N6 > 0)
    figure(6);
    plot (middleX(1:Nx-1), helium3(1:Nx-1, 1), 'red', middleX(1:Nx-1),helium3(1:Nx-1, 2), 'green', middleX(1:Nx-1),helium3(1:Nx-1, 3), 'blue');
    title ('anisotropy helium3');
    xlabel ('x');
    ylabel ('T1/T2 - T2/T1');
    grid ;
end;
