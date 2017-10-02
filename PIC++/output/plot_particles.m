clear;
load proton.dat;
load electron.dat;
load positron.dat;
load alpha.dat;
load deuterium.dat;
load helium3.dat;
load oxygen_plus_3.dat;
load silicon_plus_1.dat;
load initialParameters.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

cv = initialParameters(10);
omega = initialParameters(21);
omegaElectron = initialParameters(20);

N1=size(electron,1);
N2=size(proton,1);
N3=size(positron,1);
N4=size(alpha,1);
N5=size(deuterium,1);
N6 = size(helium3,1);
N7 = size(oxygen_plus_3,1);
N8 = size(silicon_plus_1,1);
if(N1 > 0)
    figure(1);
    plot ((electron(1:N1,1) - Xfile(2))*omegaElectron/cv, electron(1:N1,7),'blue');
    title ('electrons');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('electrons','Location','northeast');
    grid ; 
    figure(2);
    plot ((electron(1:N1,1) - Xfile(2))*omegaElectron/cv, electron(1:N1,2), 'blue');
    title ('electrons');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('electrons','Location','northeast');
    grid ; 
end;

if(N2 > 0)
    figure(3);
    plot ((proton(1:N2,1) - Xfile(2))*omegaElectron/cv, proton(1:N2,7),'blue');
    title ('protons');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('protons','Location','northeast');
    grid ; 
    figure(4);
    plot ((proton(1:N2,1) - Xfile(2))*omegaElectron/cv, proton(1:N2,2), 'blue');
    title ('protons');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('protons','Location','northeast');
    grid ; 
end;

if(N3 > 0)
    figure(5);
    plot ((positron(1:N3,1) - Xfile(2))*omegaElectron/cv, positron(1:N3,7),'blue');
    title ('positrons');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('positrons','Location','northeast');
    grid ; 
    figure(6);
    plot ((positron(1:N3,1) - Xfile(2))*omegaElectron/cv, positron(1:N3,2), 'blue');
    title ('positrons');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('positrons','Location','northeast');
    grid ; 
end;

if(N4 > 0)
    figure(7);
    plot ((alpha(1:N4,1) - Xfile(2))*omegaElectron/cv, alpha(1:N4,7), 'blue');
    title ('alphas');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('alphas','Location','northeast');
    grid ; 
    figure(8);
    plot ((alpha(1:N4,1) - Xfile(2))*omegaElectron/cv, alpha(1:N4,2), 'blue');
    title ('alphas');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('alphas','Location','northeast');
    grid ; 
end;

if(N5 > 0)
    figure(9);
    plot ((deuterium(1:N5,1) - Xfile(2))*omegaElectron/cv, deuterium(1:N5,7),'blue');
    title ('deuterium');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('deuterium','Location','northeast');
    grid ; 
    figure(10);
    plot ((deuterium(1:N5,1) - Xfile(2))*omegaElectron/cv, deuterium(1:N5,2), 'blue');
    title ('deuterium');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('deuterium','Location','northeast');
    grid ; 
end;

if(N6 > 0)
    figure(11);
    plot ((helium3(1:N6,1) - Xfile(2))*omegaElectron/cv, helium3(1:N6,7),'blue');
    title ('helium3');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('helium3','Location','northeast');
    grid ; 
    figure(12);
    plot ((helium3(1:N6,1) - Xfile(2))*omegaElectron/cv, helium3(1:N6,2), 'blue');
    title ('helium3');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('helium3','Location','northeast');
    grid ; 
end;
if(N7 > 0)
    figure(13);
    plot ((oxygen_plus_3(1:N7,1) - Xfile(2))*omegaElectron/cv, oxygen_plus_3(1:N7,7),'blue');
    title ('oxygen+3');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('helium3','Location','northeast');
    grid ; 
    figure(14);
    plot ((oxygen_plus_3(1:N7,1) - Xfile(2))*omegaElectron/cv, oxygen_plus_3(1:N7,2), 'blue');
    title ('oxygen+3');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('helium3','Location','northeast');
    grid ; 
end;
if(N8 > 0)
    figure(15);
    plot ((silicon_plus_1(1:N8,1) - Xfile(2))*omegaElectron/cv, silicon_plus_1(1:N8,7),'blue');
    title ('silicon+1');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('helium3','Location','northeast');
    grid ; 
    figure(16);
    plot ((silicon_plus_1(1:N8,1) - Xfile(2))*omegaElectron/cv, silicon_plus_1(1:N8,2), 'blue');
    title ('silicon+1');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('helium3','Location','northeast');
    grid ; 
end;