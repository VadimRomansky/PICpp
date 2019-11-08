clear;
load tamc_radial_profile.dat;
load xfile.dat

Nx = size(xfile,1);
Nt = size(tamc_radial_profile,1)/Nx;

a = Nt - 1;

figure(1);
plot (xfile(1:Nx),tamc_radial_profile(a*Nx + (1:Nx), 1), 'red');
title ('x');
xlabel ('x');
ylabel ('x');
grid ;

figure(2);
plot (xfile(1:Nx),tamc_radial_profile(a*Nx + (1:Nx), 2), 'red');
title (' V');
xlabel ('x');
ylabel ('V');
grid ;

figure(3);
plot (xfile(1:Nx),tamc_radial_profile(a*Nx + (1:Nx), 3), 'red');
title ('{\rho}');
xlabel ('x');
ylabel ('{\rho}');
grid ;

figure(4);
plot (xfile(1:Nx),tamc_radial_profile(a*Nx + (1:Nx), 4), 'red');
title ('P');
xlabel ('x');
ylabel ('P');
grid ;

figure(5);
plot (xfile(1:Nx),tamc_radial_profile(a*Nx + (1:Nx), 5), 'red');
title ('CRp');
xlabel ('x');
ylabel ('CRp');
grid ;

figure(6);
plot (xfile(1:Nx),tamc_radial_profile(a*Nx + (1:Nx), 6), 'red');
title ('T');
xlabel ('x');
ylabel ('T');
grid ;

figure(7);
plot (xfile(1:Nx),tamc_radial_profile(a*Nx + (1:Nx), 7), 'red');
title ('B-B0');
xlabel ('x');
ylabel ('B-B0');
grid ;

figure(8);
plot (xfile(1:Nx),tamc_radial_profile(a*Nx + (1:Nx), 8), 'red');
title ('Eb');
xlabel ('x');
ylabel ('Eb');
grid ;

figure(9);
plot (xfile(1:Nx),tamc_radial_profile(a*Nx + (1:Nx), 9), 'red');
title ('{\tau}');
xlabel ('x');
ylabel ('{\tau}');
grid ;

figure(10);
plot (xfile(1:Nx),tamc_radial_profile(a*Nx + (1:Nx), 10), 'red');
title ('CRn');
xlabel ('x');
ylabel ('CRn');
grid ;

figure(11);
plot (xfile(1:Nx),tamc_radial_profile(a*Nx + (1:Nx), 11), 'red');
title ('Vsc');
xlabel ('x');
ylabel ('Vsc');
grid ;