clear;
load protons.dat;
load electrons.dat;
load positrons.dat;
load alphas.dat;
load deuterium.dat;
load helium3.dat

N1=size(electrons,1);
N2=size(protons,1);
N3=size(positrons,1);
N4=size(alphas,1);
N5=size(deuterium,1);
N6 = size(helium3,1);
if(N1 > 0)
    figure(1);
    plot (electrons(1:N1,1), electrons(1:N1,2),'blue');
    title ('electrons');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('electrons','Location','northeast');
    grid ; 
    figure(2);
    plot (electrons(1:N1,1), electrons(1:N1,3), 'blue');
    title ('electrons');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('electrons','Location','northeast');
    grid ; 
end;

if(N2 > 0)
    figure(3);
    plot (protons(1:N2,1), protons(1:N2,2),'blue');
    title ('protons');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('protons','Location','northeast');
    grid ; 
    figure(4);
    plot (protons(1:N2,1), protons(1:N2,3), 'blue');
    title ('protons');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('protons','Location','northeast');
    grid ; 
end;

if(N3 > 0)
    figure(5);
    plot (positrons(1:N3,1), positros(1:N3,2),'blue');
    title ('positrons');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('positrons','Location','northeast');
    grid ; 
    figure(6);
    plot (positrons(1:N3,1), positrons(1:N3,3), 'blue');
    title ('positrons');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('positrons','Location','northeast');
    grid ; 
end;

if(N4 > 0)
    figure(7);
    plot (alphas(1:N4,1), alphas(1:N4,2), 'blue');
    title ('alphas');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('alphas','Location','northeast');
    grid ; 
    figure(8);
    plot (alphas(1:N4,1), alphas(1:N4,3), 'blue');
    title ('alphas');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('alphas','Location','northeast');
    grid ; 
end;

if(N5 > 0)
    figure(9);
    plot (deuterium(1:N5,1), deuterium(1:N5,2),'blue');
    title ('deuterium');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('deuterium','Location','northeast');
    grid ; 
    figure(10);
    plot (deuterium(1:N5,1), deuterium(1:N5,3), 'blue');
    title ('deuterium');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('deuterium','Location','northeast');
    grid ; 
end;

if(N6 > 0)
    figure(11);
    plot (helium3(1:N6,1), helium3(1:N6,2),'blue');
    title ('helium3');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    %legend('helium3','Location','northeast');
    grid ; 
    figure(12);
    plot (helium3(1:N6,1), helium3(1:N6,3), 'blue');
    title ('helium3');
    xlabel ('x cm');
    ylabel ('Px g*cm/s');
    %legend('helium3','Location','northeast');
    grid ; 
end;