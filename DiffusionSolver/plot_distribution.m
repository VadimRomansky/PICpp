clear;
F = importdata('./output/F99.dat');
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

Fx(1:Nx) = 0;
Fx2(1:Nx) = 0;
Fx3(1:Nx) = 0;
for i = 1:Nx,
    Fx(i) = F1(fix(Nz/2)+1, fix(Ny/2)+1, i, 2);
    Fx2(i) = F2(fix(Nz/2) + 1, fix(Ny/2) + 1, i, 2);
    Fx3(i) = F3(fix(Nz/2) + 1, fix(Ny/2) + 1, i, 2);
end;

D = 5;
t = 50*1000*0.01;

Fx1(1:Nx) = 0;
for i = fix(Nx/2)+1:Nx,
    Fx1(i) = Fx2(fix(Nx/2)+1);
end;
for i = 1:fix(Nx/2),
    Fx1(i) = Fx2(fix(Nx/2)+1)*exp((xgrid(i) - xgrid(fix(Nx/2)+1))/D);
end;

%for i = 1:Nx,
%    Fx1(i) = Fx2(fix(Nx/2)+1)*exp(-(xgrid(i) - xgrid(fix(Nx/2)+1))^2/(4*D*t));
%end;

Fp(1:Np) = 0;
Fp1(1:Np) = 0;
Fp2(1:Np) = 0;
Fp3(1:Np) = 0;
for l = 1:Np,
    Fp(l) = F1(fix(Nz/2)+1, fix(Ny/2)+1, fix(Nx/2)+1, l);
    Fp2(l) = F2(fix(Nz/2)+1, fix(Ny/2)+1, fix(Nx/2)+1,l);
    Fp3(l) = F3(fix(Nz/2) + 1, fix(Ny/2) + 1, fix(Nx/2)+1, l);
    Fp1(l) = Fp2(3)*pgrid(3)/pgrid(l);
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
plot(pgrid, Fp, 'b');
plot(pgrid, 1.1*Fp2, 'g');
plot(pgrid, 1.2*Fp3, 'r');
plot(pgrid, 1.3*Fp1, 'm');
legend('explicit','implicit','GMRES', 'analytic');

figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'linear');
title ('F_{x}');
xlabel ('x');
plot(xgrid, Fx, 'b', linewidth = 2);
plot(xgrid, 1.1*Fx2, 'g', linewidth = 2);
plot(xgrid, 1.2*Fx3, 'r', linewidth = 2);
plot(xgrid, 1.3*Fx1, 'm', linewidth = 2);
legend('neumann','integral','dirichlet','analytical');
