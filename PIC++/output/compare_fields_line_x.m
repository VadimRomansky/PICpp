clear;
load EfieldX.dat;
load BfieldX.dat;
load Xfile.dat;

load Yfile.dat;
load Zfile.dat;
load initialParameters.dat;

directory_name = '../../../tristan-mp-pitp/output/';
file_name = 'flds.tot';
file_number = '.008';

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx-1;
NB = (Nx-1);
Nt = (size(EfieldX, 1)/NE);
%Nt=2;
%Nt = 6;
NtB = (size(BfieldX, 1)/NB);
ynumber = 2;
znumber = 2;
%Nt = 10;
c = fix(Nt)-1;

Ex(1:Nx-1) = 0;
Ey(1:Nx-1) = 0;
Ez(1:Nx-1) = 0;

Bx(1:Nx-1) = 0;
By(1:Nx-1) = 0;
Bz(1:Nx-1) = 0;

Bnorm(1:Nx-1) = 0;

B0=initialParameters(19);

middleX(1:Nx-1) = 0;
Xgrid(1:Nx) = 0;
cv = initialParameters(10);
omega = initialParameters(21);
omegaElectron = initialParameters(20);


for i=1:Nx-1,
   Xgrid(i) = (Xfile(i) - Xfile(2))*omega/cv;
   %Xgrid(i) = (Xfile(i) - Xfile(2));
   Ex(i) = EfieldX((i) + c*NE, 1);
   Ey(i) = EfieldX((i) + c*NE, 2);
   Ez(i) = EfieldX((i) + c*NE, 3);
end;
%ynumber = 1;
%znumber = 1;
for i = 1:Nx-1,
    
   middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2))*omega/cv;
   Bx(i) = BfieldX(i + c*NB, 1);
   By(i) = BfieldX(i + c*NB, 2);
   Bz(i) = BfieldX(i + c*NB, 3);
end;

full_name = strcat(directory_name, file_name, file_number);
BxT = hdf5read(full_name,'bx');
ByT = hdf5read(full_name,'by');
BzT = hdf5read(full_name,'bz');
ExT = hdf5read(full_name,'ex');
EyT = hdf5read(full_name,'ey');
EzT = hdf5read(full_name,'ez');

NxT = size(BxT, 1);
NyT = size(ByT, 2);
ypoint = fix(NyT/2) + 1;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

fieldsFactorP = sqrt(Bx(Nx-1)*Bx(Nx-1) + By(Nx-1)*By(Nx-1) + Bz(Nx-1)*Bz(Nx-1));
fieldsFactorT = sqrt(BxT(Nx-1, ypoint)*BxT(Nx-1, ypoint) + ByT(Nx-1, ypoint)*ByT(Nx-1, ypoint) + BzT(Nx-1, ypoint)*BzT(Nx-1, ypoint));

figure(1);
plot (Xgrid(1:Nx-1),Ex(1:Nx-1)/fieldsFactorP, 'red',Xgrid(1:Nx-1),ExT(1:Nx-1,ypoint)/fieldsFactorT, 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('E_x gauss');
legend('PIDARAC','Tristan','Location','northeast');
grid ;

figure(2);
plot (Xgrid(1:Nx-1),Ey(1:Nx-1)/fieldsFactorP, 'red',Xgrid(1:Nx-1),EyT(1:Nx-1,ypoint)/fieldsFactorT, 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('E_y gauss');
legend('PIDARAC','Tristan','Location','northeast');
grid ;

figure(3);
plot (Xgrid(1:Nx-1),Ez(1:Nx-1)/fieldsFactorP, 'red',Xgrid(1:Nx-1),EzT(1:Nx-1,ypoint)/fieldsFactorT, 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('E_z gauss');
legend('PIDARAC','Tristan','Location','northeast');
grid ;

figure(4);
plot (Xgrid(1:Nx-1),Bx(1:Nx-1)/fieldsFactorP, 'red',Xgrid(1:Nx-1),BxT(1:Nx-1,ypoint)/fieldsFactorT, 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('B_x gauss');
legend('PIDARAC','Tristan','Location','northeast');
grid ;

figure(5);
plot (Xgrid(1:Nx-1),By(1:Nx-1)/fieldsFactorP, 'red',Xgrid(1:Nx-1),ByT(1:Nx-1,ypoint)/fieldsFactorT, 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('B_y gauss');
legend('PIDARAC','Tristan','Location','northeast');
grid ;

figure(6);
plot (Xgrid(1:Nx-1),Bz(1:Nx-1)/fieldsFactorP, 'red',Xgrid(1:Nx-1),BzT(1:Nx-1,ypoint)/fieldsFactorT, 'blue');
%title ('E_x');
xlabel ('x');
ylabel ('B_z gauss');
legend('PIDARAC','Tristan','Location','northeast');
grid ;