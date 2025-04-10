clear;
F = importdata('F5.dat');
xgrid = importdata('xgrid.dat');
ygrid = importdata('ygrid.dat');
zgrid = importdata('zgrid.dat');
pgrid = importdata('pgrid.dat');

Nx = length(xgrid);
Ny = length(ygrid);
Nz = length(zgrid);
Np = length(pgrid);

F1(1:Nz, 1:Ny, 1:Nx, 1:Np) = 0;

for k = 1:Nz,
    for j = 1:Ny,
        for i = 1:Nx,
            for l = 1:Np,
                F1(k,j,i,l) = F(Np*Nx*Ny*(k-1) + Np*Nx*(j-1) + Np*(i-1) + l);
            end;
        end;
    end;
end;

Fx(1:Nx) = 0;
for i = 1:Nx,
    Fx(i) = F1(fix(Nz/2)+1, fix(Ny/2)+1, i, 1);
end;

Fp(1:Np) = 0;
for l = 1:Np,
    Fp(l) = F1(fix(Nz/2)+1, fix(Ny/2)+1, fix(Nx/2)+1, l);
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{p}');
xlabel ('p/mc');
ylabel ('F(p)*p^3');
plot(pgrid, Fp);

figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'linear');
title ('F_{x}');
xlabel ('x');
plot(xgrid, Fx);