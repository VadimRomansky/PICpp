clear;
load generalAnisotropy.dat;
load particleTypes.dat;
N1=1;
N2=size(generalAnisotropy,1);
Ntypes = size(particleTypes,1);
linearAnisotropy(1:N2, 1:Ntypes) = 0;

%omega_plasma = 4.432550293*10^10;
omega_plasma = 2*3.14159*generalAnisotropy(2,2)/generalAnisotropy(2,3);
%omega_plasma = 4.37*10^11;
omega_gyro_a = 1.982193107*10^8;
%gamma = 0.01*omega_gyro_a;
gamma = 0.01*2*3.14159*omega_gyro_a/omega_plasma;
set(0, 'DefaultLineLineWidth', 2);
for i = 1:Ntypes,
    if(particleTypes(i) > 0)
        linearAnisotropy(1,i) = generalAnisotropy(1, 4 + 3*(i-1));
        for j = 2:N2,
            linearAnisotropy(j,i) = linearAnisotropy(1, i)*exp(-2*gamma*generalAnisotropy(j,2));
        end
        figure(3*(i-1) +1);
        plot (generalAnisotropy(1:N2,2), generalAnisotropy(1:N2,4 + 3*(i-1)), 'red');
        %if (i == 6)
        %   plot (generalAnisotropy(1:N2,2), generalAnisotropy(1:N2,3 + i), 'red',generalAnisotropy(1:N2,2), linearAnisotropy(1:N2,i), 'blue');
        %end;
        title ('anisotropy');
        xlabel ('{{t w_p}/{2\pi}}');
        ylabel ('T1/T2 - T2/T1');
        grid ;      
        figure(3*(i-1) + 2);
        plot (generalAnisotropy(1:N2,2), generalAnisotropy(1:N2,5 + 3*(i-1)), 'red');
        title ('temperature parallel');
        xlabel ('{{t w_p}/{2\pi}}');
        ylabel ('T2');
        grid ;     
        figure(3*(i-1) + 3);
        plot (generalAnisotropy(1:N2,2), generalAnisotropy(1:N2,6 + 3*(i-1)), 'red');
        title ('temperature normal');
        xlabel ('{{t w_p}/{2\pi}}');
        ylabel ('T1');
        grid ;  
    end;
end;

