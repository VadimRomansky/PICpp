clear;
load protons.dat;
load electrons.dat;
load positrons.dat;
load alphas.dat;

N1=size(protons,1);
N2=size(electrons,1);
N3=size(positrons,1);
N4=size(alphas,1);

if (N3 > 0) && (N4 > 0)
    figure(1);
    plot (protons(1:N1,1), protons(1:N1,2), 'red', electrons(1:N2,1), electrons(1:N2,2), 'blue', positrons(1:N3,1), positrons(1:N3,2), 'green', alphas(1:N4,1), alphas(1:N4,2), 'black');
    title ('particles');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    legend('protons', 'electrons', 'positrons', 'alphas','Location','northeast');
    grid ;

    figure(2);
    plot (protons(1:N1,1), protons(1:N1,3), 'red', electrons(1:N2,1), electrons(1:N2,3), 'blue', positrons(1:N3,1), positrons(1:N3,3), 'green', alphas(1:N4,1), alphas(1:N4,3), 'black');
    title ('particles');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    legend('protons', 'electrons', 'positrons', 'alphas','Location','northeast');
    grid ;
else
    figure(1);
    plot (protons(1:N1,1), protons(1:N1,2), 'red', electrons(1:N2,1), electrons(1:N2,2), 'blue');
    title ('particles');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    legend('protons', 'electrons','Location','northeast');
    grid ;

    figure(2);
    plot (protons(1:N1,1), protons(1:N1,3), 'red', electrons(1:N2,1), electrons(1:N2,3), 'blue');
    title ('particles');
    xlabel ('x cm');
    ylabel ('P g*cm/s');
    legend('protons', 'electrons','Location','northeast');
    grid ;
end;