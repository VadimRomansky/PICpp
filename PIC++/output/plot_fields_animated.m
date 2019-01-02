clear;
Efield = importdata('Efield.dat');
Bfield = inportdata('Bfield.dat');
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
load initialParameters.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx*Ny*Nz;
NB = (Nx-1)*(Ny-1)*(Nz-1);
Nt = (size(Efield, 1)/NE);
NtB = (size(Bfield, 1)/NB);
ynumber = 1;
znumber = 1;

frameTime = 1.0/Nt;

a = 0;

Ex(1:Nx) = 0;
Ey(1:Nx) = 0;
Ez(1:Nx) = 0;
maxEx = 0;
minEx = 0;
maxEy = 0;
minEy = 0;
maxEz = 0;
minEz = 0;

Bx(1:Nx-1) = 0;
By(1:Nx-1) = 0;
Bz(1:Nx-1) = 0;
Bnorm(1:Nx-1) = 0;
maxBx = 0;
minBx = 0;
maxBy = 0;
minBy = 0;
maxBz = 0;
minBz = 0;

B0=initialParameters(19);

middleX(1:Nx-1) = 0;
Xgrid(1:Nx) = 0;
cv = initialParameters(10);
omega = initialParameters(20);


for i=1:Nx,
   Xgrid(i) = (Xfile(i) - Xfile(2))*omega/cv;
   Ex(i) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber, 1);
   Ey(i) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber, 2);
   Ez(i) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber, 3);
   for a = 0:(Nt-1),
       tempEx = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 1);
       if(tempEx > maxEx),
           maxEx = tempEx;
       end;
       if(tempEx < minEx),
           minEx = tempEx;
       end;
       tempEy = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 2);
       if(tempEy > maxEy),
           maxEy = tempEy;
       end;
       if(tempEy < minEy),
           minEy = tempEy;
       end;
       tempEz = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 3);
       if(tempEz > maxEz),
           maxEz = tempEz;
       end;
       if(tempEz < minEz),
           minEz = tempEz;
       end;
   end;
end;

%maxEx = maxEx*1.2;
%minEx = minEx*1.2;
factor = 1.0;
if(maxEx > 0),
    if maxEx > 1,
        while maxEx > 1,
            maxEx = maxEx/10;
            factor = factor*10;
        end;
        maxEx = ceil(maxEx)*factor;
    else
        while maxEx < 1,
            maxEx = maxEx*10;
            factor = factor/10;
        end;
        maxEx = ceil(maxEx)*factor;
    end;
end;

factor = 1.0;
if(minEx < 0),
    if minEx < -1,
        while minEx < -1,
            minEx = minEx/10;
            factor = factor*10;
        end;
        minEx = floor(minEx)*factor;
    else
        while minEx > -1,
            minEx = minEx*10;
            factor = factor/10;
        end;
        minEx = floor(minEx)*factor;
    end;
end;

%maxEy = maxEy*1.2;
%minEy = minEy*1.2;
factor = 1.0;
if(maxEy > 0),
    if maxEy > 1,
        while maxEy > 1,
            maxEy = maxEy/10;
            factor = factor*10;
        end;
        maxEy = ceil(maxEy)*factor;
    else
        while maxEy < 1,
            maxEy = maxEy*10;
            factor = factor/10;
        end;
        maxEy = ceil(maxEy)*factor;
    end;
end;

factor = 1.0;
if(minEy < 0),
    if minEy < -1,
        while minEy < -1,
            minEy = minEy/10;
            factor = factor*10;
        end;
        minEy = floor(minEy)*factor;
    else
        while minEy > -1,
            minEy = minEy*10;
            factor = factor/10;
        end;
        minEy = floor(minEy)*factor;
    end;
end;

%maxEz = maxEz*1.2;
%minEz = minEz*1.2;
factor = 1.0;
if(maxEz > 0),
    if maxEz > 1,
        while maxEz > 1,
            maxEz = maxEz/10;
            factor = factor*10;
        end;
        maxEz = ceil(maxEz)*factor;
    else
        while maxEz < 1,
            maxEz = maxEz*10;
            factor = factor/10;
        end;
        maxEz = ceil(maxEz)*factor;
    end;
end;

factor = 1.0;
if(minEz < 0),
    if minEz < -1,
        while minEz < -1,
            minEz = minEz/10;
            factor = factor*10;
        end;
        minEz = floor(minEz)*factor;
    else
        while minEz > -1,
            minEz = minEz*10;
            factor = factor/10;
        end;
        minEz = floor(minEz)*factor;
    end;
end;

for i = 1:Nx-1,
   middleX(i) = (0.5*(Xfile(i) + Xfile(i+1)) - Xfile(2))*omega/cv;
   Bx(i) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber), 1);
   By(i) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber), 2);
   Bz(i) = Bfield(((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber), 3);
   Bnorm(i) = sqrt(By(i)*By(i) + Bz(i)*Bz(i))/B0;
   for a = 0:(Nt-1),
       tempBx = Bfield((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber + a*NB, 1);
       if(tempBx > maxBx),
           maxBx = tempBx;
       end;
       if(tempBx < minBx),
           minBx = tempBx;
       end;
       tempBy = Bfield((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber + a*NB, 2);
       if(tempBy > maxBy),
           maxBy = tempBy;
       end;
       if(tempBy < minBy),
           minBy = tempBy;
       end;
       tempBz = Bfield((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber + a*NB, 3);
       if(tempBz > maxBz),
           maxBz = tempBz;
       end;
       if(tempBz < minBz),
           minBz = tempBz;
       end;
   end;
end;

%maxBx = maxBx*1.2;
%minBx = minBx*1.2;
factor = 1.0;
if(maxBx > 0),
    if maxBx > 1,
        while maxBx > 1,
            maxBx = maxBx/10;
            factor = factor*10;
        end;
        maxBx = ceil(maxEx)*factor;
    else
        while maxBx < 1,
            maxBx = maxBx*10;
            factor = factor/10;
        end;
        maxBx = ceil(maxBx)*factor;
    end;
end;

factor = 1.0;
if(minBx < 0),
    if minBx < -1,
        while minBx < -1,
            minBx = minBx/10;
            factor = factor*10;
        end;
        minBx = floor(minBx)*factor;
    else
        while minBx > -1,
            minBx = minBx*10;
            factor = factor/10;
        end;
        minBx = floor(minBx)*factor;
    end;
end;

%maxBy = maxBy*1.2;
%minBy = minBy*1.2;
factor = 1.0;
if(maxBy > 0),
    if maxBy > 1,
        while maxBy > 1,
            maxBy = maxBy/10;
            factor = factor*10;
        end;
        maxBy = ceil(maxBy)*factor;
    else
        while maxBy < 1,
            maxBy = maxBy*10;
            factor = factor/10;
        end;
        maxBy = ceil(maxBy)*factor;
    end;
end;

factor = 1.0;
if(minBy < 0),
    if minBy < -1,
        while minBy < -1,
            minBy = minBy/10;
            factor = factor*10;
        end;
        minBy = floor(minBy)*factor;
    else
        while minBy > -1,
            minBy = minBy*10;
            factor = factor/10;
        end;
        minBy = floor(minBy)*factor;
    end;
end;

%maxBz = maxBz*1.2;
%minBz = minBz*1.2;
factor = 1.0;
if(maxBz > 0),
    if maxBz > 1,
        while maxBz > 1,
            maxBz = maxBz/10;
            factor = factor*10;
        end;
        maxBz = ceil(maxBz)*factor;
    else
        while maxBz < 1,
            maxBz = maxBz*10;
            factor = factor/10;
        end;
        maxBz = ceil(maxBz)*factor;
    end;
end;

factor = 1.0;
if(minBz < 0),
    if minBz < -1,
        while minBz < -1,
            minBz = minBz/10;
            factor = factor*10;
        end;
        minBz = floor(minBz)*factor;
    else
        while minBz > -1,
            minBz = minBz*10;
            factor = factor/10;
        end;
        minBz = floor(minBz)*factor;
    end;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 2);

figure(1);
%title ('E_x');
xlabel ('x\omega_p/c');
ylabel ('E_x gauss');
grid on;
hold on;
axis([Xgrid(1) Xgrid(Nx) minEx maxEx]);
fig = plot (Xgrid(1:Nx),Ex(1:Nx), 'red');
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
mov(1:height, 1:width, 1:1, 1:Nt)=0;
f = getframe(gcf);
[mov(:,:,1,1), map]=rgb2ind(f.cdata, colorcube(256));
for a = 1:(Nt-1),
    for i=1:Nx,
        %Xgrid(i) = (Xfile(i) - Xfile(2))*omega/cv;
        Ex(i) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 1);
    end;
    set(fig, 'Ydata', Ex);
    f = getframe(gcf);
    mov(:,:,1,a+1)=rgb2ind(f.cdata, map);
    pause(frameTime);
end;
imwrite(mov, map, 'Ex.gif','DelayTime',0,'LoopCount',1);

figure(2);
%title ('E_y');
xlabel ('x\omega_p/c');
ylabel ('E_y gauss');
grid on;
hold on;
axis([Xgrid(1) Xgrid(Nx) minEy maxEy]);
fig = plot (Xgrid(1:Nx),Ey(1:Nx), 'red');
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
mov(1:height, 1:width, 1:1, 1:Nt)=0;
f = getframe(gcf);
[mov(:,:,1,1), map]=rgb2ind(f.cdata, colorcube(256));
for a = 1:(Nt-1),
    for i=1:Nx,
        %Xgrid(i) = (Xfile(i) - Xfile(2))*omega/cv;
        Ey(i) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 2);
    end;
    set(fig, 'Ydata', Ey);
    f = getframe(gcf);
    mov(:,:,1,a+1)=rgb2ind(f.cdata, map);
    pause(frameTime);
end;
imwrite(mov, map, 'Ey.gif','DelayTime',0,'LoopCount',1);

figure(3);
%title ('E_z');
xlabel ('x\omega_p/c');
ylabel ('E_z gauss');
grid on;
hold on;
axis([Xgrid(1) Xgrid(Nx) minEz maxEz]);
fig = plot (Xgrid(1:Nx),Ez(1:Nx), 'red');
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
mov(1:height, 1:width, 1:1, 1:Nt)=0;
f = getframe(gcf);
[mov(:,:,1,1), map]=rgb2ind(f.cdata, colorcube(256));
for a = 1:(Nt-1),
    for i=1:Nx,
        %Xgrid(i) = (Xfile(i) - Xfile(2))*omega/cv;
        Ez(i) = Efield((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*NE, 3);
    end;
    set(fig, 'Ydata', Ez);
    f = getframe(gcf);
    mov(:,:,1,a+1)=rgb2ind(f.cdata, map);
    pause(frameTime);
end;
imwrite(mov, map, 'Ez.gif','DelayTime',0,'LoopCount',1);

figure(4);
%title ('E_x');
xlabel ('x\omega_p/c');
ylabel ('B_x gauss');
grid on;
hold on;
axis([middleX(1) middleX(Nx-1) minBx maxBx]);
fig = plot (middleX(1:Nx-1),Bx(1:Nx-1), 'red');
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
mov(1:height, 1:width, 1:1, 1:Nt)=0;
f = getframe(gcf);
[mov(:,:,1,1), map]=rgb2ind(f.cdata, colorcube(256));
for a = 1:(Nt-1),
    for i=1:Nx-1,
        %Xgrid(i) = (Xfile(i) - Xfile(2))*omega/cv;
        Bx(i) = Bfield((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber + a*NB, 1);
    end;
    set(fig, 'Ydata', Bx);
    f = getframe(gcf);
    mov(:,:,1,a+1)=rgb2ind(f.cdata, map);
    pause(frameTime);
end;
imwrite(mov, map, 'Bx.gif','DelayTime',0,'LoopCount',1);

figure(5);
%title ('E_x');
xlabel ('x\omega_p/c');
ylabel ('B_y gauss');
grid on;
hold on;
axis([middleX(1) middleX(Nx-1) minBy maxBy]);
fig = plot (middleX(1:Nx-1),By(1:Nx-1), 'red');
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
mov(1:height, 1:width, 1:1, 1:Nt)=0;
f = getframe(gcf);
[mov(:,:,1,1), map]=rgb2ind(f.cdata, colorcube(256));
for a = 1:(Nt-1),
    for i=1:Nx-1,
        %Xgrid(i) = (Xfile(i) - Xfile(2))*omega/cv;
        By(i) = Bfield((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber + a*NB, 2);
    end;
    set(fig, 'Ydata', By);
    f = getframe(gcf);
    mov(:,:,1,a+1)=rgb2ind(f.cdata, map);
    pause(frameTime);
end;
imwrite(mov, map, 'By.gif','DelayTime',0,'LoopCount',1);

figure(6);
%title ('E_x');
xlabel ('x\omega_p/c');
ylabel ('B_z gauss');
grid on;
hold on;
axis([middleX(1) middleX(Nx-1) minBz maxBz]);
fig = plot (middleX(1:Nx-1),Bz(1:Nx-1), 'red');
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
mov(1:height, 1:width, 1:1, 1:Nt)=0;
f = getframe(gcf);
[mov(:,:,1,1), map]=rgb2ind(f.cdata, colorcube(256));
for a = 1:(Nt-1),
    for i=1:Nx-1,
        %Xgrid(i) = (Xfile(i) - Xfile(2))*omega/cv;
        Bz(i) = Bfield((Nz-1)*(Ny-1)*(i-1) + (Nz-1)*(ynumber-1) + znumber + a*NB, 3);
    end;
    set(fig, 'Ydata', Bz);
    f = getframe(gcf);
    mov(:,:,1,a+1)=rgb2ind(f.cdata, map);
    pause(frameTime);
end;
imwrite(mov, map, 'Bz.gif','DelayTime',0,'LoopCount',1);