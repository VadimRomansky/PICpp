clear;
load distribution_protons.dat;
load distribution_electrons.dat;

Np = 1000;


Nt = size(distribution_protons, 1)/Np;

a = 0;
b = fix(Nt/2);
c = Nt - 2;

Fp(1:Np, 1:3) = 0;
Fe(1:Np, 1:3) = 0;
Fp2(1:Np, 1:3) = 0;
Fe2(1:Np, 1:3) = 0;

Pp(1:Np, 1:3) = 0;
Pe(1:Np, 1:3) = 0;

for i=1:Np,
   Pp(i,1) = distribution_protons(i + a*Np,1);
   Pp(i,2) = distribution_protons(i + b*Np,1);
   Pp(i,3) = distribution_protons(i + c*Np,1);
   
   Fp(i,1) = distribution_protons(i + a*Np, 2);
   Fp(i,2) = distribution_protons(i + b*Np, 2);
   Fp(i,3) = distribution_protons(i + c*Np, 2);
   
   Fp2(i,1) = Fp(i,1)*Pp(i,1)*Pp(i,1);
   Fp2(i,2) = Fp(i,2)*Pp(i,2)*Pp(i,2);
   Fp2(i,3) = Fp(i,3)*Pp(i,3)*Pp(i,3);
   
   Pe(i,1) = distribution_electrons(i + a*Np,1);
   Pe(i,2) = distribution_electrons(i + b*Np,1);
   Pe(i,3) = distribution_electrons(i + c*Np,1);
   
   Fe(i,1) = distribution_electrons(i + a*Np, 2);
   Fe(i,2) = distribution_electrons(i + b*Np, 2);
   Fe(i,3) = distribution_electrons(i + c*Np, 2);
   
   Fe2(i,1) = Fe(i,1)*Pe(i,1)*Pe(i,1);
   Fe2(i,2) = Fe(i,2)*Pe(i,2)*Pe(i,2);
   Fe2(i,3) = Fe(i,3)*Pe(i,3)*Pe(i,3);
end;

figure(1);
plot (Pp(1:Np,1),Fp(1:Np,1), 'red',Pp(1:Np,2),Fp(1:Np,2), 'green',Pp(1:Np,3),Fp(1:Np,3), 'blue');
title ('protons distinution function');
xlabel ('p');
ylabel ('Fp');
grid ;

figure(2);
plot (Pe(1:Np,1),Fe(1:Np,1), 'red',Pe(1:Np,2),Fe(1:Np,2), 'green',Pe(1:Np,3),Fe(1:Np,3), 'blue');
title ('electrons distinution function');
xlabel ('p');
ylabel ('Fe');
grid ;

figure(3);
plot (Pp(1:Np,1),Fp2(1:Np,1), 'red',Pp(1:Np,2),Fp2(1:Np,2), 'green',Pp(1:Np,3),Fp2(1:Np,3), 'blue');
title ('protons distinution function');
xlabel ('p');
ylabel ('Fp*p^2');
grid ;

figure(4);
plot (Pe(1:Np,1),Fe2(1:Np,1), 'red',Pe(1:Np,2),Fe2(1:Np,2), 'green',Pe(1:Np,3),Fe2(1:Np,3), 'blue');
title ('electrons distinution function');
xlabel ('p');
ylabel ('Fe*p^2');
grid ;