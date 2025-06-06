clear;
F = importdata('./output/F1.dat');
xgrid = importdata('./output/xgrid.dat');
ygrid = importdata('./output/ygrid.dat');
zgrid = importdata('./output/zgrid.dat');
pgrid = importdata('./output/pgrid.dat');

Nx = length(xgrid);
Ny = length(ygrid);
Nz = length(zgrid);
Np = length(pgrid);

F1(1:Nz, 1:Ny, 1:Nx, 1:Np) = 0;
F2(1:Nz, 1:Ny, 1:Nx, 1:Np) = 0;
F3(1:Nz, 1:Ny, 1:Nx, 1:Np) = 0;

for k = 1:Nz,
    for j = 1:Ny,
        for i = 1:Nx,
            for l = 1:Np,
                F1(k,j,i,l) = F(Np*Nx*Ny*(k-1) + Np*Nx*(j-1) + Np*(i-1) + l, 1);
                F2(k,j,i,l) = F(Np*Nx*Ny*(k-1) + Np*Nx*(j-1) + Np*(i-1) + l, 2);
                F3(k,j,i,l) = F(Np*Nx*Ny*(k-1) + Np*Nx*(j-1) + Np*(i-1) + l, 3);
            end;
        end;
    end;
end;

Fx1(1:Nx) = 0;
Fx2(1:Nx) = 0;
Fx3(1:Nx) = 0;
for i = 1:Nx,
    Fx1(i) = F1(fix(Nz/2)+1, fix(Ny/2)+1, i, 2);
    Fx2(i) = F2(fix(Nz/2)+1, fix(Ny/2)+1, i, 2);
    Fx3(i) = F3(fix(Nz/2)+1, fix(Ny/2)+1, i, 2);
end;

Fxa(1:Nx) = 0;
for i = fix(Nx/2)+1:Nx,
    Fxa(i) = Fx2(fix(Nx/2)+1);
end;
for i = 1:fix(Nx/2),
    Fxa(i) = Fx2(fix(Nx/2)+1)*exp((xgrid(i)-xgrid(fix(Nx/2)+1))/5);
end;

Fp1(1:Np) = 0;
Fp2(1:Np) = 0;
Fp3(1:Np) = 0;
Fpa(1:Np) = 0;
for l = 1:Np,
    Fp1(l) = F1(fix(Nz/2)+1, fix(Ny/2)+1, fix(Nx/2)+1, l);
    Fp2(l) = F2(fix(Nz/2)+1, fix(Ny/2)+1, fix(Nx/2)+1, l);
    Fp3(l) = F3(fix(Nz/2)+1, fix(Ny/2)+1, fix(Nx/2)+1, l);
    Fpa(l) = Fp2(2)*pgrid(2)/pgrid(l);
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
plot(pgrid, Fp1, 'r');
plot(pgrid, Fp2, 'g');
plot(pgrid, Fp3, 'b');
plot(pgrid, Fpa, 'm');

figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'linear');
title ('F_{x}');
xlabel ('x');
plot(xgrid, Fx1, 'r');
plot(xgrid, Fx2, 'g');
plot(xgrid, Fx3, 'b');
plot(xgrid, Fxa, 'm');