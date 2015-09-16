clear;
load dielectricTensorFile.dat;
load Xfile.dat;

Nx = size(Xfile, 1);

NE = Nx;
NB = Nx - 1;
Nt = fix(size(dielectricTensorFile,1)/Nx) - 1;

a = 0;
b = fix(Nt/2);
c = Nt;

Mxx(1:NE, 1:3) = 0;
Mxy(1:NE, 1:3) = 0;
Mxz(1:NE, 1:3) = 0;
Myx(1:NE, 1:3) = 0;
Myy(1:NE, 1:3) = 0;
Myz(1:NE, 1:3) = 0;
Mzx(1:NE, 1:3) = 0;
Mzy(1:NE, 1:3) = 0;
Mzz(1:NE, 1:3) = 0;

middleX(1:NB) = 0;


for i=1:NE,
   Mxx(i,1) = dielectricTensorFile(i + a*NE, 1);
   Mxy(i,1) = dielectricTensorFile(i + a*NE, 2);
   Mxz(i,1) = dielectricTensorFile(i + a*NE, 3);

   Myx(i,1) = dielectricTensorFile(i + a*NE, 4);
   Myy(i,1) = dielectricTensorFile(i + a*NE, 5);
   Myz(i,1) = dielectricTensorFile(i + a*NE, 6);
   
   Mzx(i,1) = dielectricTensorFile(i + a*NE, 7);
   Mzy(i,1) = dielectricTensorFile(i + a*NE, 8);
   Mzz(i,1) = dielectricTensorFile(i + a*NE, 9);
   
   Mxx(i,2) = dielectricTensorFile(i + b*NE, 1);
   Mxy(i,2) = dielectricTensorFile(i + b*NE, 2);
   Mxz(i,2) = dielectricTensorFile(i + b*NE, 3);

   Myx(i,2) = dielectricTensorFile(i + b*NE, 4);
   Myy(i,2) = dielectricTensorFile(i + b*NE, 5);
   Myz(i,2) = dielectricTensorFile(i + b*NE, 6);
   
   Mzx(i,2) = dielectricTensorFile(i + b*NE, 7);
   Mzy(i,2) = dielectricTensorFile(i + b*NE, 8);
   Mzz(i,2) = dielectricTensorFile(i + b*NE, 9);
   
   Mxx(i,3) = dielectricTensorFile(i + c*NE, 1);
   Mxy(i,3) = dielectricTensorFile(i + c*NE, 2);
   Mxz(i,3) = dielectricTensorFile(i + c*NE, 3);

   Myx(i,3) = dielectricTensorFile(i + c*NE, 4);
   Myy(i,3) = dielectricTensorFile(i + c*NE, 5);
   Myz(i,3) = dielectricTensorFile(i + c*NE, 6);
   
   Mzx(i,3) = dielectricTensorFile(i + c*NE, 7);
   Mzy(i,3) = dielectricTensorFile(i + c*NE, 8);
   Mzz(i,3) = dielectricTensorFile(i + c*NE, 9);
end;

for i=1:NB,
   middleX(i) = 0.5*(Xfile(i) + Xfile(i+1));
end;
figure(1);
plot(Xfile(1:NE,1),Mxx(1:NE,3), 'Color',[1, 0, 0]);
hold on;
plot(Xfile(1:NE,1),Mxy(1:NE,3), 'Color',[0.7, 0.3, 0]);
hold on;
plot(Xfile(1:NE,1),Mxz(1:NE,3), 'Color',[0.7, 0, 0.3]);
hold on;
plot(Xfile(1:NE,1),Myx(1:NE,3), 'Color',[0.3, 0.7, 0]);
hold on;
plot(Xfile(1:NE,1),Myy(1:NE,3), 'Color',[0, 1, 0]);
hold on;
plot(Xfile(1:NE,1),Myz(1:NE,3), 'Color',[0, 0.7, 0.3]);
hold on;
plot(Xfile(1:NE,1),Mzx(1:NE,3), 'Color',[0.3, 0, 0.7]);
hold on;
plot(Xfile(1:NE,1),Mzy(1:NE,3), 'Color',[0, 0.3, 0.7]);
hold on;
plot(Xfile(1:NE,1),Mzz(1:NE,3), 'Color',[0, 0, 1]);
legend('{{\mu}_{xx}}', '{{\mu}_{xy}}','{{\mu}_{xz}}', '{{\mu}_{yx}}', '{{\mu}_{yy}}','{{\mu}_{yz}}', '{{\mu}_{zx}}', '{{\mu}_{zy}}','{{\mu}_{zz}}','Location','northeast');
title ('dielectric tensor');
xlabel ('x/r_g');
ylabel ('{\mu}');
grid ;

figure(2);
plot(Xfile(1:NE,1),Mxx(1:NE,1), 'r', Xfile(1:NE,1),Mxx(1:NE,2), 'g', Xfile(1:NE,1),Mxx(1:NE,3), 'b');
title ('dielectric tensor');
xlabel ('x/r_g');
ylabel ('{\mu_{xx}}');
grid ;