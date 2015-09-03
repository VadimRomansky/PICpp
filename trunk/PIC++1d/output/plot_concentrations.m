clear;
load Xfile.dat;
load concentrations.dat;

Nx = size(Xfile, 1) - 1;

Nt = size(concentrations, 1)/Nx;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

electron_concentration(1:Nx, 1:3) = 0;
proton_concentration(1:Nx, 1:3) = 0;
charge_density(1:Nx, 1:3) = 0;


for i=1:Nx,   
   electron_concentration(i, 1) = concentrations(i + a*Nx, 1);
   electron_concentration(i, 2) = concentrations(i + b*Nx, 1);
   electron_concentration(i, 3) = concentrations(i + c*Nx, 1);
   proton_concentration(i, 1) = concentrations(i + a*Nx, 2);
   proton_concentration(i, 2) = concentrations(i + b*Nx, 2);
   proton_concentration(i, 3) = concentrations(i + c*Nx, 2);
   charge_density(i, 1) = concentrations(i + a*Nx, 3);
   charge_density(i, 2) = concentrations(i + b*Nx, 3);
   charge_density(i, 3) = concentrations(i + c*Nx, 3);
end;

figure(1);
plot (Xfile(1:Nx,1),electron_concentration(1:Nx, 1), 'red', Xfile(1:Nx,1), electron_concentration(1:Nx, 2), 'green', Xfile(1:Nx,1), electron_concentration(1:Nx, 3), 'blue');
title ('electrons');
xlabel ('x/r_g');
ylabel ('n cm^-3');
grid ;

figure(2);
plot (Xfile(1:Nx,1),proton_concentration(1:Nx, 1), 'red', Xfile(1:Nx,1), proton_concentration(1:Nx, 2), 'green', Xfile(1:Nx,1), proton_concentration(1:Nx, 3), 'blue');
title ('protons');
xlabel ('x/r_g');
ylabel ('n cm^-3');
grid ;

figure(3);
plot (Xfile(1:Nx,1),charge_density(1:Nx, 1), 'red', Xfile(1:Nx,1), charge_density(1:Nx, 2), 'green', Xfile(1:Nx,1), charge_density(1:Nx, 3), 'blue');
title ('charge density');
xlabel ('x/r_g');
ylabel ('n sgs*cm^-3');
grid ;